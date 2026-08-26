# Metodologia WASP: Particionamento Espectral e Matching de Partições

> **Documento histórico.** Este arquivo descreve um fluxo de validação anterior
> e é mantido para proveniência da pesquisa; não define a documentação atual do
> pacote. Consulte o [índice atual](INDEX.md) e a
> [implementação PCSPM atual](matching.md).

> **Propósito**: Documentação científica precisa de todos os passos de processamento.  
> **Baseado em**: Revisão real dos códigos de implementação.  
> **Status**: Corrigido conforme implementação atual (fevereiro 2025)

---

## Sumário

1. [Visão Geral do Pipeline](#1-visão-geral-do-pipeline)
2. [Fontes de Dados e Pré-processamento](#2-fontes-de-dados-e-pré-processamento)
3. [Algoritmo de Particionamento Espectral (Hanson & Phillips 2001)](#3-algoritmo-de-particionamento-espectral)
4. [Pré-processamento Específico CFOSAT](#4-pré-processamento-específico-cfosat)
5. [Matching Temporal WW3 (Passo 02)](#5-matching-temporal-ww3-passo-02)
6. [Matching de Pares de Partições (Análise)](#6-matching-de-pares-de-partições)
7. [Classificação de Regime de Swell](#7-classificação-de-regime-de-swell)
8. [Métricas de Verificação](#8-métricas-de-verificação)

---

## 1. Visão Geral do Pipeline

O pipeline WASP consiste em dois componentes principais:

### 1.1 Biblioteca WASP (Particionamento)

| Módulo | Função |
|--------|--------|
| `wasp.partition.partition_spectrum()` | Particiona espectros 2D usando algoritmo de watershed |
| `wasp.io_sar.load_sar()` | Carrega dados SAR (Sentinel-1) |
| `wasp.io_cfosat.load_cfosat()` | Carrega dados CFOSAT SWIM |
| `wasp.io_ww3.load_ww3()` | Carrega espectros WW3 |
| `wasp.io_ndbc.load_ndbc()` | Carrega espectros NDBC |
| `wasp.wave_params.calculate_wave_parameters()` | Calcula Hs, Tp, Dp por partição |

**Saída do particionamento**: Arquivos CSV individuais por observação com parâmetros de onda para cada partição (P1, P2, P3).

### 1.2 Análise e Matching (Fora da Biblioteca WASP)

O matching entre observações e WW3 é executado por **scripts externos** que:
1. Carregam CSVs particionados de observações e WW3
2. Realizam **matching temporal** entre WW3 e observações
3. Realizam **matching de partições** usando critérios espectrais
4. Geram arquivos de análise com pares matched

---

## 2. Fontes de Dados e Pré-processamento

### 2.1 SAR (Sentinel-1)

- **Input**: Arquivo NetCDF, grupo `obs_params`, variável `E2d[index, freq, dir]`
- **Unidades**: m²·s·rad⁻¹ (densidade espectral de energia)
- **Direções**: Convenção oceanográfica (0° = N, positivo no sentido horário)
- **Filtro de qualidade**: Apenas dados com `L2_partition_quality_flag == 0`

### 2.2 NDBC (Boias)

- **Input**: Arquivo NetCDF por estação, espectros 2D reconstruídos via Maximum Entropy Method (MEM)
- **Limitação**: Artefatos da reconstrução MEM; requer `merge_factor` mais agressivo (0.7)
- **Amostragem temporal**: Configurável, padrão 6 horas

### 2.3 CFOSAT SWIM

- **Input**: Arquivo NetCDF por órbita, espectros em espaço wavenumber-direção
- **Conversão de unidades**: De wavenumber $k$ para frequência $f$ usando dispersão linear em águas profundas
- **Filtro de comprimento de onda**: Energia em comprimentos de onda < 500 m é zerificada
- **Ambiguidade direcional**: Instrumento side-looking; não distingue $\theta$ de $\theta+180°$

### 2.4 WW3 (Modelo)

- **Input**: Arquivo NetCDF com espectro colocalizado em locais de observação
- **Experimento**: `exp_02-st4-uost-psi-400s-era5-b143-ic5-noref`
- **Forçante**: Reanálise ERA5

---

## 3. Algoritmo de Particionamento Espectral

O particionamento segue **Hanson & Phillips (2001)** implementado em `wasp.partition.partition_spectrum()`.

### 3.1 Limiar de Energia (Modo Adaptativo)

$$E_{\text{thr}} = P_{98}\{E(f,\theta) : E(f,\theta) > 0\}$$

Todos os produção utilizam:
- **p = 98º percentil** (apenas os 2% de valores de energia mais altos)
- **Significado físico**: Apenas as cristas espectrais mais energéticas são candidatos a picos
- **Benefício**: Limiar relativo a cada espectro individual

### 3.2 Identificação de Picos

Para cada ponto $(f_i, \theta_j)$ com $E(f_i, \theta_j) \geq E_{\text{thr}}$:

1. Examina vizinhança 3×3
2. Encontra o vizinho com maior energia
3. Se o ponto **é** o máximo local: **ICOD[i,j] = 22** → classificado como pico
4. Senão: ICOD[i,j] = JY*10 + IX codifica a direção do vizinho com maior energia

Depois:
- Ordena picos por energia decrescente
- Retém apenas **max_partitions = 3** picos maiores
- Trata direção como **periódica** (0°/360° wrap)

### 3.3 Segmentação por Watershed (Propagação ICOD)

Cada ponto espectral é atribuído à partição mais próxima via **dois passos**:

**Passo 1 - Propagação ICOD**:
- Iterações forward/backward em frequência e direção
- Máx 50 iterações até convergência
- Cada ponto copia a partição do vizinho indicado por seu ICOD

**Passo 2 - Preenchimento de Zeros**:
- Pontos restantes não-atribuídos (valor = 0)
- Votação entre 8 vizinhos
- Padrão para partição 1 se sem vizinhos rotulados

**Resultado**: `MASK[i,j]` ∈ {1, 2, ..., nmask} — atribuição completa

### 3.4 Cálculo de Distância entre Picos

Converte-se picos para espaço Cartesiano espectral:

$$x_k = f_k \cos(\theta_k), \quad y_k = f_k \sin(\theta_k)$$

Distância euclidiana quadrática:

$$d^2_{ij} = (x_i - x_j)^2 + (y_i - y_j)^2$$

### 3.5 Parâmetro de Espalhamento de Picos

Para cada partição $k$, a variância em espaço Cartesiano:

$$E^{\text{ip}}_k = \text{Var}(x) + \text{Var}(y)$$

Onde:
- $x = f\cos(\theta)$, $y = f\sin(\theta)$
- **Eip pequeno**: Partição estreita e concentrada
- **Eip grande**: Partição ampla e difusa

### 3.6 Merging de Partições Sobrepostas

Duas partições $i$ e $j$ são **merged** se:

$$d^2_{ij} \leq \alpha \cdot E^{\text{ip}}_i \quad \text{E} \quad d^2_{ij} \leq \alpha \cdot E^{\text{ip}}_j$$

Onde $\alpha$ é o `merge_factor`:

| merge_factor | Comportamento | Recomendado para |
|---|---|---|
| **0.315** | Conservador — mantém sistemas distintos | SAR, CFOSAT |
| **0.5** | Moderado — padrão | WW3 |
| **0.7** | Agressivo — combina sistemas próximos | NDBC |

### 3.7 Integração de Energia por Partição

Momento espectral de ordem 0 para partição $k$ usando quadratura trapezoidal:

$$m_{0,k} = \sum_{i,j \in \text{partição } k} E_{ij} \cdot w_i \cdot \Delta\theta$$

Onde $w_i$ são pesos trapezoidais em frequência.

**Altura significativa**: $H_{s,k} = 4\sqrt{m_{0,k}}$

**Verificação**: Se $|\sum_k m_{0,k} - m_0| > 10^{-4}$ m², um aviso é emitido.

### 3.8 Extração de Parâmetros de Onda

Para cada partição $k$:

- **Período de pico**: $T_{p,k} = 1/f_{\text{peak},k}$ (frequência de máxima energia)
- **Direção de pico**: Média ponderada por energia em $f_{\text{peak},k}$:

$$D_{p,k} = \arctan\left(\frac{\sum_j E(f_{\text{peak},k},\theta_j)\sin\theta_j}{\sum_j E(f_{\text{peak},k},\theta_j)\cos\theta_j}\right) \text{ [deg, oceanográfica]}$$

### 3.9 Reordenação por Energia

Partições são renumeradas de forma que:

$$H_{s,\text{P1}} \geq H_{s,\text{P2}} \geq H_{s,\text{P3}}$$

Fornece convenção consistente para comparações.

### 3.10 Filtro de Energia Pós-processamento

Após particionamento, filtro de energia mínima:

$$\text{manter partição } k \iff m_{0,k} > \epsilon_E \cdot m_0$$

Com $\epsilon_E = 0.01$ **(1% de energia total)**.

Partições abaixo desse limiar são **zerificadas** no CSV.

---

## 4. Pré-processamento Específico CFOSAT

CFOSAT SWIM é um instrumento radar side-looking. Duas correções são aplicadas **antes** do particionamento:

### 4.1 Remoção de Ambiguidade Direcional 180°

SWIM não consegue distinguir $\theta$ de $\theta + 180°$ (dependência com inclinação da onda, não sinal).

**Estratégia** (`criterion='ocean_swell'`):

Para cada par antipodal $(\theta_i, \theta_i + 180°)$:

1. Integra energia em cada direção: $E_1 = \int E(f, \theta_i)\,df$, $E_2 = \int E(f, \theta_i+180°)\,df$
2. Se um está no quadrante oceânico (SW–W–NW: 135°–315°) e outro não → mantém oceânico
3. Se ambos ou nenhum está no quadrante → mantém maior energia

Aplicado a todos os pares antes de `partition_spectrum`.

### 4.2 Correção Pós-particionamento 180°

Após particionamento, rotação sistemática:

$$D_p \leftarrow (D_p + 180°) \mod 360°$$

Corrige a convenção do instrumento (look direction → wave propagation direction).

---

## 5. Matching Temporal WW3 (Passo 02)

Para cada observação (SAR/NDBC/CFOSAT), o espectro WW3 colocalizado é selecionado:

1. Carrega arquivo NetCDF WW3 pré-colocalizado (um arquivo por localização de observação)
2. Identifica timestep WW3 $t_{\text{WW3}}$ mais próximo de $t_{\text{obs}}$
3. **Aceita se**: $|t_{\text{obs}} - t_{\text{WW3}}| \leq \Delta t_{\text{match}} = 3.0$ h (NDBC)

Para SAR e CFOSAT (com tempo único), utiliza-se o WW3 mais próximo independentemente do offset.

---

## 6. Matching de Pares de Partições

**IMPORTANTE**: O matching atual **NÃO utiliza algoritmo Hungarian como documentado anteriormente**. 

A implementação real utiliza **greedy matching com critérios espectrais**. Encontra-se em scripts externos (não na biblioteca WASP).

### 6.1 Temporal Colocation

Um par obs–WW3 é aceito se:

$$|t_{\text{obs}} - t_{\text{WW3}}| \leq \Delta t_{\text{max}} = 1.0 \text{ h}$$

(Mais rigoroso que os 3.0 h do passo 02, pois análise requer simultaneidade verdadeira).

### 6.2 Detecção de Partição Ativa

Uma partição é considerada **ativa** se:

$$H_{s,k} > H_{s,\text{min}} = 0.05 \text{ m}$$

Guarda contra entradas próximas de zero que passaram pelo filtro 1%.

**Para SAR**: Filtro adicional de período mínimo para excluir wind sea:

$$T_{p,\text{obs}} \geq T_{p,\text{min}}^{\text{SAR}} = 12.0 \text{ s}$$

Reflete limitações de SAR para capturar espectro wind sea (corte de azimute, velocity bunching).

### 6.3 Estratégia de Matching (Greedy)

Para cada partição SAR com energia $> 1\%$:

1. **Filtra** WW3 com $T_p \geq 12$ s (evita ruído, mantém comparabilidade)
2. **Busca** melhor partição WW3 baseada em proximidade espectral:
   - $|\Delta T_p| \leq 1.0$ s
   - $|\Delta D_p| \leq 15.0°$ (diferença angular mínima)
3. **Rejeita hard** qualquer par fora desses limites
4. **Pontuação** para desempate: $\text{score} = |\Delta T_p| + \frac{|\Delta D_p|}{40}$

**Ordem fixa**: Tenta SAR P1 → WW3 P1, SAR P2 → WW3 P2, SAR P3 → WW3 P3

(NÃO permite reordenação flexível como sugerido em testes anteriores).

### 6.4 Cálculo de Diferença Angular

$$|\Delta D_p| = \left|((D_{\text{WW3}} - D_{\text{obs}} + 180°) \mod 360°) - 180°\right| \in [0°, 180°]$$

Respeita natureza circular das direções.

### 6.5 Flag de Qualidade

Cada par matched é anotado:

| `match_quality` | Condição |
|---|---|
| `strict` | $\|\Delta T_p\| \leq 1.0$ s **E** $\|\Delta D_p\| \leq 15.0°$ |
| `relaxed` | Passou rejeição hard mas fora dos limites strict |

(Diferente da documentação anterior que tinha thresholds adicionais).

---

## 7. Classificação de Regime de Swell

Partições individuais são classificadas por período de pico:

| Regime | Critério | Interpretação |
|---|---|---|
| **G1** | $T_p < 12$ s | Wind sea / swell jovem |
| **G2** | $12 \leq T_p < 16$ s | Swell intermediário |
| **G3** | $T_p \geq 16$ s | Swell de longo período / remoto |

Um **caso mixed sea** contém simultaneamente:
- Pelo menos uma partição G1 **E**
- Pelo menos uma partição G2 ou G3

Senão, o regime do caso é determinado pela partição dominante (P1).

---

## 8. Métricas de Verificação

### 8.1 Métricas a Nível de Caso

Por cada par obs–WW3:

$$\text{Bias} = \overline{(H_{s,\text{WW3}} - H_{s,\text{obs}})} = \frac{1}{N}\sum_i (H_{s,\text{WW3},i} - H_{s,\text{obs},i})$$

$$\text{RMSE} = \sqrt{\frac{1}{N}\sum_i (H_{s,\text{WW3},i} - H_{s,\text{obs},i})^2}$$

$$\text{Normalized BIAS} = \frac{\text{Bias}}{\overline{H_{s,\text{obs}}}}$$

$$\text{Normalized RMSE} = \frac{\text{RMSE}}{\overline{H_{s,\text{obs}}}}$$

$$r = \text{Pearson correlation coefficient}$$

### 8.2 Número de Partições

Para cada caso:

| Métrica | Definição |
|---|---|
| $n_{\text{obs}}$ | Partições com $H_s > 0.05$ m em observação |
| $n_{\text{WW3}}$ | Partições com $H_s > 0.05$ m em WW3 |
| $n_{\text{matched}}$ | Pares obs–WW3 successfully matched |
| $n_{\text{unmat\_obs}}$ | Partições obs sem par WW3 (missed detections) |
| $n_{\text{unmat\_ww3}}$ | Partições WW3 sem par obs (false alarms) |

---

## 9. Resumo de Diferenças da Documentação Anterior

| Aspecto | Documentação Anterior | Implementação Real |
|---|---|---|
| Algoritmo de matching | Hungarian (optimal assignment) | Greedy (first compatible) |
| Matriz de custo | Complexa com pesos | Simples: score = Δtp + Δdp/40 |
| Reordenação | Flexível (P1↔P2↔P3) | Fixa (P1→P1, P2→P2, P3→P3) |
| Tolerância Tp | $\max(2.0\text{ s}, 0.20 \cdot T_p)$ | **Fixa: 1.0 s** |
| Tolerância Dp | 60° (hard rejection) | **15.0° (rejeição total)** |
| Tp mínimo SAR | 10.0 s | **12.0 s** |
| Energia mínima | 0.05 m (Hs) | **1% de energia total** |

---

## 10. Notas Importantes

1. **WASP** é apenas particionamento — matching é externo
2. **Matching atual é determinístico** (não probabilístico)
3. **Greedy matching** pode deixar partições sem pair (accepted by design)
4. **Tolerâncias** são rigorosas (1 s em Tp, 15° em Dp)
5. **SAR wind sea** é explicitamente filtrado por Tp < 12 s
6. **CFOSAT 180°** requer pós-correção após particionamento
7. **Energia < 1%** é zeros fados no CSV (não matching)

---

## Referências

- Hanson, J. L., & O. M. Phillips (2001). Automated analysis of ocean surface directional wave spectra. *J. Atmos. Oceanic Technol.*, 18, 277–293.
- Script original: `04_validate.py` (backup do projeto)
- Implementação WASP: `/src/wasp/partition.py`
