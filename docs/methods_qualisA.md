# Melhorias para Manuscrito Cientifico

> **Documento histórico.** Este texto editorial é mantido como material de
> pesquisa e não define a documentação nem o comportamento atual do pacote.
> Consulte o [índice atual](INDEX.md).

Este documento consolida as melhorias editoriais e as versoes revisadas elaboradas exclusivamente a partir do texto metodologico e da figura fornecidos.

## Prioridades Editoriais

- Eliminar a duplicacao e a redacao informal na secao de preprocessamento CFOSAT.
- Corrigir a numeracao de `2.2.2. NDBC Preprocessing` para `2.2.3. NDBC Preprocessing`.
- Uniformizar a terminologia em todo o manuscrito: *directional spectrum*, *peak period*, *peak direction*, *significant wave height*, *Partition Detection Skill* e *Partition Characterization Skill*.
- Corrigir os indices da equacao de custo: a particao simulada deve ser identificada por `j`, enquanto a observada deve ser identificada por `i`.
- Definir explicitamente que `|Delta Dp|` representa a diferenca direcional circular minima.
- Substituir o marcador bibliografico incompleto `(reference)` associado a parametrizacao IC5 pela citacao bibliografica correspondente antes da submissao.
- Reduzir a quantidade de texto dentro dos blocos da figura e concentrar a explicacao metodologica na legenda.

## Legenda Revisada da Figura

> **Figure 1. Overview of the proposed WASP framework for the analysis of directional wave spectra.** Directional spectra from WAVEWATCH III (WW3), Sentinel-1 SAR, CFOSAT-SWIM, and NDBC buoys are first preprocessed and transformed into a common spectral representation, $E(f,\theta)$, expressed in $m^2\,s\,rad^{-1}$ under the oceanographic propagation-direction convention. The WASP methodology then applies adaptive peak detection, watershed segmentation, partition merging, and bulk-parameter estimation. Observed and simulated partitions are subsequently associated through the Physics-Constrained Spectral Partition Matching (PCSPM) procedure, which combines temporal collocation, physical compatibility filters, a physics-based cost function, and Hungarian one-to-one assignment. The validation framework evaluates both partition detection skill and partition characterization skill.

## Abertura Metodologica Revisada

> A global unstructured grid (Figure 1) was employed to simulate surface ocean wave conditions and to evaluate their representation in terms of two-dimensional directional spectra. Numerical simulations were performed using WAVEWATCH III (WW3), version 7.14 (WW3DG, 2019). The grid resolution varied spatially from approximately 50 km in deep-ocean regions to 2-4 km in coastal areas.
>
> WW3 is a third-generation spectral wave model that solves the spectral action-density balance equation to represent wave generation, propagation, and dissipation across deep- and shallow-water environments. The numerical experiment consisted of a continuous global simulation for 2020, forced with ERA5 atmospheric reanalysis fields, including 10-m wind and sea-ice concentration.
>
> The ST4 source-term package was adopted to represent wind input, whitecapping dissipation, and swell attenuation. Nonlinear quadruplet wave-wave interactions were parameterized using the Discrete Interaction Approximation, whereas depth-induced wave breaking followed the Battjes and Janssen formulation. Additional processes included attenuation by unresolved sub-grid obstacles through the UOST source term and wave-sea-ice interactions using the IC5 parameterization. Wave propagation was computed using the third-order explicit propagation scheme with Garden Sprinkler Effect alleviation. The wind-input growth parameter, $\beta_{\max}$, was maintained at its default value.

## Secao CFOSAT-SWIM Revisada

### 2.2.2 CFOSAT-SWIM Spectral Preprocessing

> The SWIM instrument measures the directional distribution of wave energy using multiple off-nadir beams. Because radar backscatter is sensitive to wave slope rather than to the sign of sea-surface displacement, the instrument cannot distinguish propagation toward $\theta$ from propagation toward $\theta + 180^\circ$. Consequently, the Level-2 product contains 12 azimuthal sectors and presents an intrinsic 180-degree directional ambiguity.
>
> To obtain full directional coverage, the original spectra were expanded from 12 to 24 directional bins spanning 360 degrees. The integrated wave parameters provided in the Level-2 product were used to evaluate the directional convention adopted during processing and the consistency of the reconstructed spectra. Representative cases were subsequently compared with collocated WW3 simulations and NDBC buoy observations to confirm the directional reference and assess the ambiguity-resolution procedure. After establishing the appropriate directional convention, the duplicated component associated with the 180-degree ambiguity was removed, yielding a unique directional spectrum for partitioning analysis.
>
> The measured slope spectrum was converted into an elevation spectrum according to:
>
> $$
> E_\eta(k) = \frac{E_{\mathrm{slope}}(k)}{k^2}.
> $$
>
> The elevation spectrum was then transformed from wavenumber to frequency space using the deep-water Jacobian:
>
> $$
> E(f,\theta) = E_\eta(k)\left|\frac{dk}{df}\right|\frac{\pi}{180},
> $$
>
> resulting in a directional spectrum expressed in $m^2\,s\,rad^{-1}$.

## NDBC: Correcao de Numeracao

Use o seguinte titulo para evitar duplicacao de numeracao:

```text
2.2.3 NDBC Preprocessing
```

## Nucleo do PCSPM Revisado

> The Physics-Constrained Spectral Partition Matching (PCSPM) methodology was developed to establish physically meaningful one-to-one correspondences between observed and simulated wave systems. For each observation, the nearest WW3 spectrum within a maximum temporal separation of one hour was selected. Partitioning was independently applied to the observed and simulated spectra, and only active partitions containing more than 1% of the total spectral energy were retained.
>
> Candidate pairs were first screened using physical compatibility criteria. A pair was rejected when:
>
> $$
> \Delta T_p > \max(3\ \mathrm{s},\,0.30T_{p,\mathrm{obs}}),
> $$
>
> or when the minimum circular difference in peak direction exceeded $60^\circ$. These filters prevent associations between physically distinct systems before optimization.
>
> For compatible pairs, the matching cost was defined as:
>
> $$
> C_{ij} =
> w_{T_p}\frac{|T_{p,j}-T_{p,i}|}{T_{p,i}} +
> w_{D_p}\frac{|\Delta D_{p,ij}|}{180} +
> w_{H_s}\frac{|H_{s,j}-H_{s,i}|}{H_{s,i}},
> $$
>
> where $i$ denotes an observed partition, $j$ denotes a simulated partition, and $|\Delta D_{p,ij}|$ is the minimum circular directional difference. The weights assigned to peak period, peak direction, and significant wave height were 1.5, 1.5, and 0.5, respectively. The resulting cost matrix was optimized with the Hungarian algorithm to obtain the globally optimal one-to-one assignment. Partitions without a physically compatible counterpart were retained as unmatched systems in the subsequent analyses.

## Texto Revisado para Validacao

Substituir a redacao sobre a metrica PCA por:

> The principal metric derived from the contingency matrix was the Partition Count Agreement (PCA), which quantifies the proportion of cases in which observed and simulated spectra contain exactly the same number of active wave systems. PCA therefore provides a direct measure of Partition Detection Skill, whereas the off-diagonal elements identify systematic tendencies toward over-partitioning or under-partitioning.

## Recomendacoes para a Figura

A estrutura da figura comunica adequadamente a sequencia completa:

```text
input data -> preprocessing -> common spectral representation -> WASP -> PCSPM -> validation
```

Para a versao de submissao, recomenda-se:

- Preservar as seis etapas e o fluxo de dados atualmente apresentado.
- Ampliar formulas, siglas e eixos espectrais para garantir legibilidade em pagina de artigo.
- Reduzir frases internas dos blocos a termos metodologicos essenciais.
- Manter a explicacao completa do fluxo na legenda revisada, em vez de concentrar detalhes no interior da figura.
- Garantir consistencia visual entre as cores dos quatro tipos de dados e as setas do fluxo metodologico.
