# Methodology: Spectral Partitioning and Partition-Pair Matching

> **Historical document.** This file describes an earlier validation workflow
> and is retained for research provenance; it is not current package
> documentation. See the [current documentation index](INDEX.md) and
> [current PCSPM implementation](matching.md).
> Present-tense statements and parameter values below are relative to that
> historical external workflow, including scripts not present in this checkout.

> **Purpose**: Scientific documentation of all processing steps for use in a peer-reviewed article.  
> Generated from code review of the WASP validation pipeline (`valida_wasp`).

---

## Table of Contents

1. [Overview of the Pipeline](#1-overview-of-the-pipeline)
2. [Data Sources and Preprocessing](#2-data-sources-and-preprocessing)
3. [Spectral Partitioning Algorithm (Hanson & Phillips 2001)](#3-spectral-partitioning-algorithm)
   - 3.1 Energy threshold
   - 3.2 Peak identification
   - 3.3 Watershed segmentation
   - 3.4 Peak distance calculation
   - 3.5 Peak spreading parameter
   - 3.6 Partition merging
   - 3.7 Energy integration per partition
   - 3.8 Wave parameter extraction
   - 3.9 Reordering by energy
   - 3.10 Post-processing energy filter
4. [CFOSAT-specific Preprocessing](#4-cfosat-specific-preprocessing)
5. [WW3 Temporal Matching (Steps 02)](#5-ww3-temporal-matching-step-02)
6. [Partition-Pair Matching (Step 07 — Hungarian Algorithm)](#6-partition-pair-matching-step-07)
7. [Sea-State Regime Classification](#7-sea-state-regime-classification)
8. [Verification Metrics](#8-verification-metrics)
9. [Configuration Parameters — Complete Reference](#9-configuration-parameters)
10. [Known Issues and Parameter Notes](#10-known-issues-and-parameter-notes)

---

## 1. Overview of the Pipeline

The validation pipeline consists of five sequential steps:

| Script | Purpose |
|--------|---------|
| `01_partition_sar_new.py` | Partition SAR (Sentinel-1) observed spectra |
| `02_partition_ww3_new.py` | Partition WW3 model spectra (SAR / NDBC / CFOSAT collocations) |
| `03_partition_ndbc_new.py` | Partition NDBC buoy spectra (all available stations) |
| `05_partition_cfosat_new.py` | Partition CFOSAT SWIM spectra |
| `07_partition_system_analysis.py` | Match obs↔WW3 partition pairs and compute statistics |

All partitioning scripts produce individual CSV files per observation, named with the pattern:
```
{sensor}_{ref:03d}_{YYYYmmdd-HHMMSS}.csv
```
Each CSV contains the total-spectrum parameters and up to three partition rows (P1, P2, P3).

---

## 2. Data Sources and Preprocessing

### 2.1 SAR (Sentinel-1, scripts 01 + 02)

- **Input**: NetCDF files, group `obs_params`, variable `E2d[index, freq, dir]`
- **Units**: m²·s·rad⁻¹ (energy spectral density)
- **Directions**: oceanographic convention (0° = N, clockwise positive)
- **Collocations**: matched pairs stored in `auxdata/sar_matches_{case}_track.csv`
- **Quality filter**: `L2_partition_quality_flag == 0`

### 2.2 NDBC Buoys (script 03)

- **Input**: NetCDF per station, reconstructed 2-D spectra via Maximum Entropy Method (MEM)
- **Limitation**: MEM reconstruction introduces directional spreading artifacts; a more aggressive merge_factor is appropriate
- **Temporal sampling**: configurable every N hours (default: every 6 h)

### 2.3 CFOSAT SWIM (script 05)

- **Input**: NetCDF per track, spectra given in wavenumber-direction space  
- **Unit conversion**: from wavenumber $k$ to frequency $f$ using linear deep-water dispersion:  
  $f = \frac{1}{2\pi}\sqrt{g k}$, with Jacobian $\frac{dE}{df} = \frac{dE}{dk}\cdot\left|\frac{dk}{df}\right|$
- **Wavelength filter**: spectral energy at wavelengths < 500 m is set to zero (below SWIM sensitivity)
- **Directional ambiguity**: SWIM is a side-looking instrument and cannot distinguish $\theta$ from $\theta+180°$ (see §4)

### 2.4 WW3 Model (script 02)

- **Input**: NetCDF spectra collocated at observation locations
- **Experiment**: `exp_02-st4-uost-psi-400s-era5-b143-ic5-noref`, year 2020
- **Wind fields**: ERA5 reanalysis forcing

---

## 3. Spectral Partitioning Algorithm

The partitioning algorithm follows **Hanson & Phillips (2001)** as implemented in the `wasp.partition.partition_spectrum()` function. The seven sub-steps are described below.

### 3.1 Energy Threshold (Adaptive Mode)

All production runs use `threshold_mode='adaptive'`:

$$E_\text{thr} = P_p\left(\{E(f,\theta) : E(f,\theta) > 0\}\right)$$

where $P_p$ denotes the $p$-th percentile of the non-zero spectral values.  
In all configurations: **$p = 98$** (98th percentile).

> **Physical meaning**: only the top 2% of spectral energy values qualify as candidate peaks. This makes the threshold relative to each individual spectrum, avoiding arbitrary absolute thresholds.

### 3.2 Peak Identification

For each grid point $(f_i, \theta_j)$ with $E(f_i, \theta_j) \geq E_\text{thr}$, the algorithm examines the 3×3 neighbourhood and computes a direction code:

```
ICOD[i,j] = JY*10 + IX
```

where $(IX, JY)$ encode the index of the largest neighbour. If the point itself is the local maximum ($IX = JY = 2 \Rightarrow \text{ICOD} = 22$), it is classified as a spectral **peak**.

- Peaks are sorted by energy (highest first)
- Only the top `max_partitions = 3` peaks are retained
- The direction dimension is treated as **periodic** (0°/360° wrap)

### 3.3 Watershed Segmentation (ICOD Propagation)

Each grid point is assigned to the nearest peak by following the steepest-descent direction codes recursively:

1. **Phase 1**: iterative forward/backward sweeps along both dimensions until convergence (max 50 iterations). Each unassigned point copies the label of the neighbour indicated by its ICOD code.
2. **Phase 2**: any remaining zero-valued points are filled by majority vote among their 8 neighbours (repeated until no zeros remain). Points with no labelled neighbours default to partition 1.

Result: `MASK[i,j]` ∈ {1, 2, …, `nmask`} — a complete assignment of spectral energy to wave systems.

### 3.4 Peak Distance Calculation

Each peak is mapped to Cartesian spectral space:

$$x_k = f_k\cos(\theta_k), \quad y_k = f_k\sin(\theta_k)$$

The squared Euclidean distance between peaks $i$ and $j$ is:

$$d^2_{ij} = (x_i - x_j)^2 + (y_i - y_j)^2$$

This metric reflects physical separation in the 2-D wavenumber plane.

### 3.5 Peak Spreading Parameter

For each partition $k$, the spectral spreading (variance in Cartesian space) is:

$$E^\text{ip}_k = \frac{\sum_{(i,j)\in k} E_{ij} \Delta f_i \Delta\theta \cdot (f_i^2\cos^2\theta_j + f_i^2\sin^2\theta_j)}{m_0} - \left(\frac{\sum_{(i,j)\in k} E_{ij}\Delta f_i\Delta\theta \cdot f_i\cos\theta_j}{m_0}\right)^2 - \left(\frac{\sum_{(i,j)\in k} E_{ij}\Delta f_i\Delta\theta \cdot f_i\sin\theta_j}{m_0}\right)^2$$

where $m_0 = \sum_{i,j} E_{ij}\Delta f_i\Delta\theta$ is the total spectral energy.  
A small $E^\text{ip}_k$ indicates a narrow, concentrated partition; a large value indicates a broad partition.

### 3.6 Partition Merging

Two partitions $i$ and $j$ are **merged** (j absorbed into i) if their peaks are close relative to their own spreading:

$$d^2_{ij} \leq \alpha \cdot E^\text{ip}_i \quad \text{AND} \quad d^2_{ij} \leq \alpha \cdot E^\text{ip}_j$$

where $\alpha$ is the `merge_factor` parameter. After merging, $j$ takes the label of $i$.

| merge_factor $\alpha$ | Behaviour | Recommended for |
|---|---|---|
| **0.315** | Conservative — splits closely spaced systems | SAR, CFOSAT (high-resolution spectra) |
| **0.5** | Moderate — balanced | WW3 model default |
| **0.7** | Aggressive — combines nearby peaks | NDBC (MEM reconstruction artefacts) |

> **Important**: In all sensitivity experiments, both obs and WW3 use the **same** merge_factor so that the partition count comparison is symmetric. The batch combinations in `07_partition_system_analysis.py` test obs and WW3 with the same or different merge_factors.

### 3.7 Energy Integration per Partition

The zeroth spectral moment for partition $k$ is computed using trapezoidal quadrature in frequency:

$$m_{0,k} = \sum_{i=1}^{N_f}\sum_{j=1}^{N_\theta} E_{ij}\cdot w_i \cdot \Delta\theta \quad \text{where } [i,j]\in k$$

with trapezoidal weights $w_1 = \Delta f_1/2$, $w_i = (\Delta f_{i-1}+\Delta f_i)/2$ for $i>1$, $w_{N_f}=\Delta f_{N_f}/2$.  
Significant wave height: $H_{s,k} = 4\sqrt{m_{0,k}}$.

An energy conservation check is performed: if $|\sum_k m_{0,k} - m_0| > 10^{-4}$ m², a warning is issued.

### 3.8 Wave Parameter Extraction

For each partition $k$, the **peak period** is:

$$T_{p,k} = 1/f_{\text{peak},k}$$

where $f_{\text{peak},k}$ is the frequency of maximum energy within partition $k$.

The **peak direction** is the energy-weighted mean direction at $f_{\text{peak},k}$:

$$D_{p,k} = \arctan\!\left(\frac{\sum_j E(f_{\text{peak},k},\theta_j)\sin\theta_j}{\sum_j E(f_{\text{peak},k},\theta_j)\cos\theta_j}\right) \quad [\text{deg, oceanographic}]$$

### 3.9 Reordering by Energy

Partitions are relabelled so that **P1 always corresponds to the most energetic system**:

$$H_{s,\text{P1}} \geq H_{s,\text{P2}} \geq H_{s,\text{P3}}$$

This provides a consistent ordering convention for statistical comparisons.

### 3.10 Post-Processing Energy Filter

After partitioning, a **minimum energy fraction** filter is applied at save time:

$$\text{keep partition } k \iff m_{0,k} > \epsilon_E \cdot m_0$$

with $\epsilon_E = $ `min_energy_fraction` $= 0.01$ (1% of total energy).

Partitions below this threshold are stored as **zeros** in the CSV (i.e., $H_s = T_p = D_p = 0.0$). This prevents numerical noise from being counted as real wave systems in the subsequent analysis.

> **Critical note for the article**: the 1% threshold is applied at the *partitioning stage* (scripts 01–05), ensuring that CSV files only carry physically meaningful systems. A secondary threshold of $H_s > 0.05$ m is then applied at the *analysis stage* (script 07) to guard against near-zero residuals that survived the first filter.

---

## 4. CFOSAT-specific Preprocessing

CFOSAT SWIM is an along-track, off-nadir radar instrument. Two instrument-specific corrections are applied **before** running the standard partitioning algorithm:

### 4.1 Removal of 180° Directional Ambiguity

SWIM cannot distinguish propagation direction $\theta$ from $\theta + 180°$ because the radar cross-section depends only on the wave slope, not its sign. The spectrum therefore appears mirrored.

**Resolution strategy** (`criterion='ocean_swell'`):  
For each antipodal direction pair $(\theta_i, \theta_i + 180°)$:
1. Compute the integrated energy in each direction: $E_1 = \int E(f, \theta_i)\,df$, $E_2 = \int E(f, \theta_i+180°)\,df$
2. If one direction lies within the preferred oceanic swell quadrant (SW–W–NW: 135°–315°) and the other does not, keep the ocean-quadrant direction and zero out the other
3. If both or neither direction is in the preferred quadrant, keep the one with higher integrated energy

This is applied to all $N_\theta/2$ direction pairs before calling `partition_spectrum`.

### 4.2 Post-Partitioning 180° Correction

After partitioning, a systematic $+180°$ rotation is applied to all peak directions:

$$D_p \leftarrow (D_p + 180°) \mod 360°$$

This corrects for the instrument's measurement convention where the spectrum is expressed in the direction of radar look rather than wave propagation.

---

## 5. WW3 Temporal Matching (Step 02)

For each observation (SAR / NDBC / CFOSAT), the spatially collocated WW3 spectrum is selected as follows:

1. Load the pre-collocated WW3 NetCDF file (one file per observation location)
2. Identify the WW3 timestep $t_\text{WW3}$ closest to the observation time $t_\text{obs}$
3. Accept the match if $|t_\text{obs} - t_\text{WW3}| \leq \Delta t_\text{match}^{(02)}$

The parameter $\Delta t_\text{match}^{(02)}$ = `max_time_diff_hours` **in the config** = **3.0 h** (applies only to NDBC in script 03, which processes all available timesteps and selects within this window).

For SAR and CFOSAT, which have a single target observation time, the closest WW3 time is always used regardless of the offset (the collocation files were pre-screened).

---

## 6. Partition-Pair Matching (Step 07)

After partitioning is complete, script `07_partition_system_analysis.py` matches each observed partition to the most physically compatible WW3 partition.

### 6.1 Temporal Colocation

An obs–WW3 pair is accepted if:

$$|t_\text{obs} - t_\text{WW3}| \leq \Delta t_\text{max} = 1.0\text{ h}$$

This is stricter than the 3.0 h used in step 02/03, because the analysis requires genuine simultaneity for wave-system comparison.

### 6.2 Active Partition Detection

A partition is considered **active** if:

$$H_{s,k} > H_{s,\text{min}} = 0.05 \text{ m}$$

This is a secondary guard against near-zero entries written as zeros by the 1% energy filter.

For SAR observations, an additional spectral period filter excludes short-wave (wind-sea) partitions from matching:

$$T_{p,\text{obs}} \geq T_{p,\text{min}}^{\text{SAR}} = 10.0 \text{ s}$$

This reflects SAR limitations in capturing wind-sea spectra reliably (azimuth cut-off, velocity bunching).

### 6.3 Cost Function

For each candidate obs–WW3 pair $(p_\text{obs}, p_\text{WW3})$, a **matching cost** is computed:

$$C = w_{T_p}\frac{|\Delta T_p|}{T_{p,\text{obs}}} + w_{D_p}\frac{|\Delta D_p|}{180°} + w_{H_s}\frac{|\Delta H_s|}{H_{s,\text{obs}}}$$

with weights $w_{T_p} = w_{D_p} = w_{H_s} = 1.0$ (equal weighting).

Angular difference is computed on the circle:

$$|\Delta D_p| = \left|((D_{p,\text{WW3}} - D_{p,\text{obs}} + 180°) \mod 360°) - 180°\right| \in [0°, 180°]$$

**Hard rejection** (cost set to $\infty$) if either threshold is exceeded:

| Variable | Rejection threshold | Formula |
|---|---|---|
| $T_p$ | Hybrid absolute/relative | $|\Delta T_p| > \max(2.0\text{ s},\ 0.20 \cdot T_{p,\text{obs}})$ |
| $D_p$ | Absolute | $|\Delta D_p| > 60°$ |

The $T_p$ hybrid threshold grows proportionally with swell period, reflecting that longer-period swells have larger absolute uncertainties. Examples: $T_p=10$ s → threshold $= 2.0$ s; $T_p=18$ s → threshold $= 3.6$ s.

### 6.4 Hungarian Algorithm (Optimal Assignment)

The cost matrix $\mathbf{C} \in \mathbb{R}^{N_\text{obs} \times N_\text{WW3}}$ is solved for the minimum-cost assignment using the **Hungarian algorithm** (`scipy.optimize.linear_sum_assignment`).

Only rows and columns with at least one finite cost entry participate in the optimisation. Pairs with $C = \infty$ are excluded.

The result is a set of matched pairs, plus unmatched observations (missed detections) and unmatched WW3 partitions (false alarms).

### 6.5 Match Quality Flag

Each matched pair is annotated with a quality label:

| `match_quality` | Condition |
|---|---|
| `strict` | $|\Delta T_p| \leq \max(2.0\text{ s},\ 0.20 T_p)$ **AND** $|\Delta D_p| \leq 30°$ |
| `relaxed` | Passed hard rejection but outside strict thresholds |

The `relaxed` category captures pairs where model–obs directional differences of 30°–60° exist, which can arise from genuine model direction errors rather than system misidentification (particularly for remote swell systems).

---

## 7. Sea-State Regime Classification

Individual partitions are classified by their peak period:

| Regime | Criterion | Physical interpretation |
|---|---|---|
| **G1** | $T_p < 12$ s | Wind sea / young swell |
| **G2** | $12 \leq T_p < 16$ s | Intermediate swell |
| **G3** | $T_p \geq 16$ s | Long-period / remote swell |

An observation case is classified as **mixed sea** if it contains simultaneously at least one G1 partition **and** at least one G2 or G3 partition. Otherwise the case regime is determined by the dominant partition (P1).

---

## 8. Verification Metrics

### 8.1 Contingency Metrics

Let:
- $N_\text{match}$ = number of matched obs–WW3 partition pairs
- $N_\text{obs}$ = total active obs partitions
- $N_\text{WW3}$ = total active WW3 partitions

$$\text{POD} = \frac{N_\text{match}}{N_\text{obs}}, \quad \text{FAR} = \frac{N_\text{WW3} - N_\text{match}}{N_\text{WW3}}, \quad \text{CSI} = \frac{N_\text{match}}{N_\text{obs} + N_\text{WW3} - N_\text{match}}$$

### 8.2 Scalar Error Metrics

For each matched pair and variable $X \in \{H_s, T_p, D_p\}$:

$$\text{Bias} = \frac{1}{N}\sum_i (X_{\text{WW3},i} - X_{\text{obs},i})$$

$$\text{RMSE} = \sqrt{\frac{1}{N}\sum_i (X_{\text{WW3},i} - X_{\text{obs},i})^2}$$

$$\text{SI} = \frac{\text{RMSE}}{\overline{X}_\text{obs}}$$

Note: For $D_p$, RMSE and Bias use the circular angular difference $|\Delta D_p|$ and only its magnitude (no signed bias for direction since the metric is unsigned).

---

## 9. Configuration Parameters

### 9.1 Partitioning (all sensors, from config YAML)

| Parameter | Config key | Value | Description |
|---|---|---|---|
| Threshold percentile | `threshold_percentile` | **98** | 98th percentile of non-zero spectral energy |
| Merge factor | `merge_factor` | **0.315 / 0.5 / 0.7** | Controls aggressiveness of partition merging |
| Max partitions | `max_partitions` | **3** | Maximum wave systems per spectrum |
| Min energy fraction | `min_energy_fraction` | **0.01** | 1% of total $m_0$ — minimum to save a partition |

### 9.2 Analysis (hardcoded in `07_partition_system_analysis.py`)

| Parameter | Variable | Value | Description |
|---|---|---|---|
| Max time diff (analysis) | `MAX_TIME_DIFF_HOURS` | **1.0 h** | Obs–WW3 temporal tolerance for analysis step |
| Min Hs (active partition) | `HS_MIN` | **0.05 m** | Secondary filter at analysis time |
| $T_p$ min for SAR | `TP_MIN_MATCHING` | **10.0 s** | SAR: exclude wind sea from matching |
| $\Delta T_p$ absolute cap | `TP_MAX_DIFF` | **2.0 s** | Hard rejection threshold |
| $\Delta T_p$ relative cap | `TP_MAX_REL` | **0.20** | 20% of $T_{p,\text{obs}}$ |
| $\Delta D_p$ acceptance | `DP_MAX_DIFF` | **60°** | Hard rejection direction threshold |
| $\Delta D_p$ strict | `DP_STRICT_DIFF` | **30°** | Boundary for match_quality='strict' |

### 9.3 Regime Thresholds

| Regime | Threshold |
|---|---|
| G1 / G2 boundary | $T_p = 12$ s |
| G2 / G3 boundary | $T_p = 16$ s |

---

## 10. Known Issues and Parameter Notes

### 10.1 Two distinct `max_time_diff_hours` values

The config YAML parameter `ndbc.max_time_diff_hours = 3.0 h` is used in script 03 (selecting which WW3 timestep is closest to an NDBC observation). **This is a different parameter** from `MAX_TIME_DIFF_HOURS = 1.0 h` in script 07 (controlling which obs–WW3 file pairs enter the analysis). Both values are intentional and serve distinct purposes:

- **3.0 h** (config, step 02/03): tolerance for WW3 temporal interpolation during spectrum extraction  
- **1.0 h** (script 07): tolerance for including a case in the statistical analysis

If a more relaxed analysis matching is desired, `MAX_TIME_DIFF_HOURS` in script 07 can be increased. The config value does **not** propagate to script 07.

### 10.2 min_energy_fraction two-stage filtering

The 1% energy fraction filter is applied **at partitioning time** (scripts 01–05): partitions below $0.01 \cdot m_0$ are stored as zeros in the CSV. The $H_s > 0.05$ m filter in script 07 is an independent secondary check. Both must be satisfied for a partition to count as "active" in the analysis. Lowering `min_energy_fraction` to 0.01 (from higher values that were removing too many systems) ensures that more marginal systems are retained in the CSV and can be inspected at the analysis stage.

### 10.3 CFOSAT threshold_percentile stored as float

CFOSAT configs use `threshold_percentile: 98.0` (float) vs `98` (int) for other sensors. This causes the partition output directory to be named `partition-cfosat-98.0-0.315` instead of `partition-cfosat-98-0.315`. The `resolve_directories` function in script 07 handles this by fuzzy-matching directory names, but both values must be kept consistent to avoid ambiguity.

### 10.4 NDBC config-95-05.yaml case setting

The file `config-95-05.yaml` has `case: cfosat` / `case_name: cfosat`. This was the setting at the last save. The `case` field must be updated to `ndbc` when running script 03 with this config. The `case` field does not affect partitioning parameters — it only controls which data source (input CSV) is processed by script 02.

---

*Document generated from code review — April 2026*
