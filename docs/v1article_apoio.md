# Results and Discussion

## 1. Convergence Across Observational Datasets and the Reduction of Partition Ambiguity

Our validation framework was applied across three fundamentally different observational datasets: Synthetic Aperture Radar (SAR) measuring surface roughness backscatter, CFOSAT providing nadir-based spectral estimates, and NDBC offering direct in-situ buoy measurements. Despite their distinct sampling characteristics and measurement principles, all three datasets exhibit a consistent convergence toward percentile-preserving configurations—specifically the 98th percentile (P98)—as the optimal solution for wave-system partitioning. 

This convergence addresses a central limitation of conventional partitioning schemes, which frequently generate artificial wave systems, fail to physically merge related systems, or incorrectly merge independent ones. The percentile-preservation strategy successfully mitigates the ambiguity associated with weak secondary peaks. When transitioning from the baseline (no-percentile, merge factor 0.5) to the optimal P98 configurations, contingency-based metrics, which directly quantify partition identification and matching success, demonstrate extraordinary improvements:

*   **SAR**: Critical Success Index (CSI) increased from 0.338 to 0.621 (+84%), with diagonal agreement (Diag%) rising from 39.7% to 60.6% (+20.9 percentage points).
*   **CFOSAT**: CSI increased from 0.222 to 0.264 (+19%), with Diag% rising from 40.3% to 47.5% (+7.2 pp).
*   **NDBC**: CSI increased from 0.468 to 0.580 (+24%), with Diag% rising from 41.1% to 53.0% (+11.9 pp).

Crucially, these dramatic gains occur while traditional bulk-wave statistics (such as the Root Mean Square Error, RMSE, and correlation coefficients, r, of significant wave height, $H_s$) remain almost invariant. For instance, in the SAR dataset, while CSI improved by 84%, RMSE($H_s$) only marginally changed from 0.882 m to 0.910 m, and r($H_s$) remained stable (0.730 to 0.705). This decoupling unequivocally indicates that the P98 methodology acts predominantly as a *spectral-organization filter* rather than an *energy filter*. The technique profoundly enhances the ability to correctly identify and match distinct spectral systems without artificially modifying the bulk-energy content.

## 2. Spectral Organization over Integrated Energy: The Dominant Source of Error

A fundamental hypothesis confirmed by these results is that the principal bottleneck in Wavewatch III (WW3) representation versus observations is associated with spectral organization rather than integrated energy. If the improvements were driven by better energy prediction, one would anticipate proportional enhancements across Taylor-diagram distances and $H_s$-based metrics. Instead, Taylor diagrams display clustered configurations with nearly identical bulk skill statistics, while the contingency matrices undergo radical restructuring.

### 2.1 The Amplification of Error Across Spectral Moments

Direct physical evidence of this organizational deficiency is provided by a comparative analysis of the spectral moments ($m_n = \int_0^\infty \omega^n S(\omega) d\omega$). While $m_0$ characterizes total integrated energy, higher-order moments ($m_1$ and $m_2$) are highly sensitive to spectral width, peak sharpness, and wave-system separation.

Our findings reveal a dataset-specific monotonic amplification of model-to-observation error ratios across spectral orders. For SAR-derived swell systems, the WW3-to-observation ratios scale aggressively: $m_0$ ratios approximate 1.9, $m_1$ scales to 2.7, and $m_2$ amplifies significantly to 4.7. This ~2.5$\times$ amplification from $m_0$ to $m_2$ is the defining signature of a **spectral-broadening bias** inherent to the wave model. 

WW3 does not simply overestimate wave energy; it systematically produces broader, less sharply peaked spectra compared to observations. A broadened spectrum mathematically distributes energy across multiple frequency and directional bins. While it may conserve the correct $m_0$ (total $H_s$), the diminished peak intensity and wider energy spread severely compound matching ambiguity in partition algorithms, leading to artificial multi-system representations. 

### 2.2 The Mitigating Mechanism of the P98 Filter

Conventional spectra contain dominant primary partitions (60–80% energy), secondary partitions (15–35%), and spectral noise (<5%). The matching ambiguity is disproportionately driven by the secondary partitions. Small uncertainties in model spectral discretization (e.g., WW3's ~0.04 Hz grid) versus observational resolution (SAR ~0.1-0.2 Hz) dictate whether a secondary peak is properly matched, erroneously merged, or over-fragmented.

The P98 framework effectively masks this ambiguity. By restricting analysis to environments where $H_s$ exceeds the 98th percentile, the algorithm selects highly energetic conditions. These environments are typically characterized by more coherent, peaked spectral structures where a single dominant swell or organized local sea dominates. Consequently, the filter naturally selects regimes where the spectral-broadening bias of WW3 has less capacity to perturb the structural topology of the wave field, thereby yielding superior partition correspondence.

## 3. The Merge Factor: Statistical Optimization versus Physical Realism

The methodology explicitly utilizes the merge factor ($mf$) to investigate the tension between matching performance (statistical correspondence) and physical realism. As defined by the Hanson & Phillips (2001) algorithm, $mf$ determines the aggressiveness with which adjacent spectral peaks are consolidated based on their spreading variance.

The results reveal an asymmetric sensitivity: the observational merge factor exerts a ~5–7$\times$ stronger influence on successful matching compared to the WW3 model-side merge factor. Applying an aggressive merge factor ($mf=0.7$) systematically maximizes CSI and diagonal agreement by forcibly reducing partition fragmentation. It mitigates discrepancies arising from differing spectral resolutions and provides optimal results for pure predictive validation.

However, case studies indicate that higher merge factors compromise physical individuality. For instance, in mixed-sea regimes captured by CFOSAT and SAR, coexisting primary mature swells ($T_p \approx 14.6$ s) and younger secondary swells ($T_p \approx 11.5$ s) maintain distinct origins and kinematics. An $mf=0.3$ (conservative) correctly preserves these as independent partitions, despite resulting in an apparent statistical "mismatch" against a broadened WW3 output (which blends them into a single system). Conversely, $mf=0.7$ artificially forces the observational data to match the merged WW3 representation, inflating the CSI metric at the direct expense of wave-system individuality.

Therefore, we conclude that the merge factor is not an arbitrary tuning coefficient, but a continuous control parameter regulating the trade-off between statistical validation optics and physical interpretability.

## 4. Dataset-Dependent Behavior and Regime Sensitivities

The heterogeneity of the three datasets illuminates how spectral characteristics modulate model performance and validation responses:

*   **SAR (Swell-Dominated):** Characterized by narrow, long-period systems. SAR exhibited the most dramatic response to the P98 filter (+84% CSI). However, a persistent positive bias ($+0.414$ m to $+0.434$ m) exists across all regimes. While WW3 heavily overestimates mature swell energy, its high $r(T_p)$ values (0.748–0.878) confirm that it correctly identifies the temporal footprint of these systems. SAR's massive CSI boost indicates that swell-dominated environments suffer most from secondary partition ambiguity.
*   **CFOSAT (Broad, Mixed-Sea):** Providing nadir-derived continuous spectra, CFOSAT is the most discriminating benchmark. It captures broader spectral distributions and energetic mixed seas. Uniquely, the P98/mf=0.7 application on CFOSAT simultaneously yielded improved partition metrics (+19% CSI) and enhanced bulk-wave statistics ($r(H_s)$ from 0.839 to 0.876).
*   **NDBC (In-situ Diversity):** Representing the largest sample ($N=51,454$ post-filtering), the NDBC buoy network provides the strongest statistical baseline. Buoy measurements span a diverse range of swell and wind-sea mixtures. Here, the P98 configuration confirmed a +24% CSI enhancement, validating that the improvements are deeply physical rather than artifacts of satellite orbital sampling biases. Interestingly, in extreme high-energy regimes, WW3's characteristic underestimation bias becomes marginally more pronounced, though overall correlations remain exceptional ($r \approx 0.93$).

## 5. Implications for Future Modeling and Machine Learning

The implications of this structural decoupling extend significantly beyond standard model validation. As the wave-modeling community increasingly integrates Machine Learning (ML) for spectral post-processing, data assimilation, and climate emulation, the definition of "truth" in training datasets is paramount.

Configurations optimized purely for statistical correspondence (e.g., $mf=0.7$) will invariably strip training datasets of vital physical features, obfuscating details about independent source regions, crossing swells, and complex energy transfers. Conversely, structure-preserving configurations ($mf=0.3$) generate physically richer descriptors critical for physics-informed neural networks (PINNs) and process-based studies, despite looking statistically "weaker" under traditional L2-norm-based metrics. 

## 6. Conclusions

This investigation demonstrates that the principal bottleneck in contemporary spectral wave partitioning and observation-model matching is rooted in structural ambiguity rather than integrated energy disparities. Through the implementation of a percentile-preservation strategy coupled with an explicit, tunable merge-factor framework, we establish that:

1.  **WW3 exhibits a systematic spectral-broadening bias** that inflates errors at higher-order spectral moments, artificially deteriorating partition-matching capabilities.
2.  **The P98 filter successfully localizes coherent wave systems**, increasing the Critical Success Index by 19% to 84% across all tested sensor suites by suppressing noise-induced secondary partitions.
3.  **A strict divergence exists between statistical optimization and physical realism**: an optimization of partition-matching metrics inherently risks the artificial homogenization of physically distinct wave fields.

Ultimately, instead of prescribing a universally "optimal" validation configuration, this work provides a rigorous, physically-interpretable framework. It equips the oceanographic community to consciously select partitioning parameters based on whether the objective is operational predictive benchmarking or the detailed diagnosis of wave kinematics.