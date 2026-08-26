# Matching de Partições Espectrais: Implementação Atual vs Documentação

> **Documento histórico.** Este arquivo compara fluxos de validação anteriores
> e é mantido para proveniência da pesquisa; não define a implementação atual do
> pacote. Consulte o [índice atual](INDEX.md) e a
> [implementação PCSPM atual](matching.md).

## Resumo Executivo

A documentação anterior (matching.md, methods.md) descreve um algoritmo **Hungarian** (optimal assignment) para matching de partições. **A implementação real utiliza um algoritmo GREEDY muito mais simples**.

Esta documento esclareça:
1. O que a documentação anterior AFIRMA
2. O que REALMENTE está sendo feito
3. As implicações práticas dessa diferença

---

## Seção 1: O Que a Documentação Anterior Diz

### Documentação: Algoritmo Hungarian (6.4 em methods.md)

```
"The cost matrix C ∈ ℝ^(N_obs × N_WW3) is solved for the minimum-cost 
assignment using the Hungarian algorithm (scipy.optimize.linear_sum_assignment)."
```

Implicações da documentação:

1. **Matriz de custo $N_{\text{obs}} \times N_{\text{WW3}}$** criada
2. **Cada célula** contém score de compatibilidade entre par obs-WW3
3. **Algoritmo Hungarian** resolve assinamento ÓTIMO global
4. **Resultado**: Todos os pares possíveis são testados
5. **Ordem flexível**: P1 obs pode match com qualquer WW3 (P1, P2, ou P3)

### Documentação: Função de Custo (6.3 em methods.md)

```
C = w_Tp * |ΔTp|/Tp_obs + w_Dp * |ΔDp|/180° + w_Hs * |ΔHs|/Hs_obs

Weights: w_Tp = w_Dp = w_Hs = 1.0 (equal weighting)

Hard rejection thresholds:
- Tp: |ΔTp| > max(2.0 s, 0.20 * Tp_obs)
- Dp: |ΔDp| > 60°
```

---

## Seção 2: O Que Realmente Está Sendo Feito

### Implementação Real (04_validate.py)

**ALGORITMO**: Greedy sequential matching (busca gulosa)

```python
def find_best_match(sar_partitions, ww3_partitions, sar_pnum):
    """Find best WW3 partition for a single SAR partition"""
    
    # Step 1: Get SAR partition
    sar_data = next((p for p in sar_partitions 
                     if p['partition'] == sar_pnum), None)
    
    # Step 2: Filter SAR by Tp >= 12.0 s (wind sea filter)
    if sar_data['tp'] < 12.0:
        return sar_data, None
    
    # Step 3: Search best match in WW3
    best_ww3 = None
    min_score = float('inf')
    
    for ww3 in ww3_partitions:
        # Step 3a: Reject if Tp outside tolerance
        if abs(ww3['tp'] - sar_data['tp']) > 1.0:  # 1.0 s tolerance
            continue
        
        # Step 3b: Reject if Dp outside tolerance
        if angular_diff(ww3['dp'], sar_data['dp']) > 15.0:  # 15.0° tolerance
            continue
        
        # Step 3c: Accept if both within tolerance
        score = abs(ww3['tp'] - sar_data['tp']) + \
                angular_diff(ww3['dp'], sar_data['dp']) / 40.0
        
        if score < min_score:
            min_score = score
            best_ww3 = ww3
    
    return sar_data, best_ww3

# Main loop: For each SAR partition, find best match
for sar_pnum in [1, 2, 3]:
    sar_match, ww3_match = find_best_match(sar_partitions, 
                                            ww3_partitions, 
                                            sar_pnum)
    if sar_match and ww3_match:
        save_pair(sar_pnum, ww3_match['partition'])
```

### Comparação Detalhada

| Aspecto | Documentação | Implementação Real |
|---|---|---|
| **Estrutura** | Matriz de custo 2D | Loop sequencial P1,P2,P3 |
| **Otimização** | Global (Hungarian) | Local (first compatible) |
| **Ordem de match** | Flexível (qualquer P↔P) | Fixa (P1→P1, P2→P2, P3→P3) |
| **Tolerância Tp** | $\max(2.0, 0.20 \cdot Tp)$ | **Fixa: 1.0 s** |
| **Tolerância Dp** | 60° | **15.0°** |
| **Rejeição Tp < 10s** | Mencionada | **Tp < 12.0 s (SAR)** |
| **Matriz Hs** | Incluída no custo | **NÃO usada no matching** |
| **Algoritmo** | `linear_sum_assignment` | Comparação simples |

---

## Seção 3: Implicações Práticas

### 3.1 Ordem Fixa vs Flexível

**Documentação promete**:
```
Se SAR tem P1 (wind sea, Tp=10s) e P2 (swell, Tp=15s)
E WW3 tem P1 (swell, Tp=16s) e P2 (wind sea, Tp=9s)

Hungarian resolveria como:
- SAR P1 (10s) → WW3 P2 (9s) [match wind sea com wind sea]
- SAR P2 (15s) → WW3 P1 (16s) [match swell com swell]
```

**Implementação real faz**:
```
Loop simples:
  SAR P1 (10s) → Rejeita (Tp < 12s) = NENHUM MATCH
  SAR P2 (15s) → WW3 P1 (16s) [Δtp=1s, dentro do limite] = MATCH
  
Resultado: SAR P1 fica sem pair (missed detection)
           SAR P2 matchado com WW3 P1
```

### 3.2 Tolerâncias Muito Mais Rigorosas

**Exemplo com espectro real**:

Observação tem:
- P1: Hs=2.5m, **Tp=13.0s**, Dp=220°
- P2: Hs=1.2m, Tp=8.5s (filtrado por Tp<12), Dp=195°

WW3 tem:
- P1: Hs=2.3m, **Tp=12.0s**, Dp=210°
- P2: Hs=1.5m, Tp=9.0s, Dp=200°

**Documentação** (`max(2.0, 0.20*13)=2.6s, ΔDp<60°`):
```
|13.0 - 12.0| = 1.0 s < 2.6 s ✓
|220 - 210| = 10° < 60° ✓
→ MATCH VÁLIDO (strict)
```

**Implementação real**:
```
|13.0 - 12.0| = 1.0 s = 1.0 s ✓ EXATO NO LIMITE
|220 - 210| = 10° < 15° ✓
→ MATCH VÁLIDO
```

**Mas se**:
- WW3 P1 tivesse Tp=13.2s:
  - **Documentação**: 13.2-13.0=0.2s < 2.6s → aceita
  - **Implementação**: 13.2-13.0=0.2s > 1.0s → **REJEITA**

### 3.3 Sem Matriz Hs

**Documentação**:
```
Cost = 1.0 * |ΔTp|/Tp + 1.0 * |ΔDp|/180 + 1.0 * |ΔHs|/Hs
```

**Implementação**:
```
Cost = |ΔTp| + |ΔDp|/40
(Hs NÃO é incluído)
```

Implicação: Um match com grande diferença em Hs pode ser aceito se Tp e Dp estão próximos.

Exemplo:
- SAR: Hs=1.0m, Tp=14s, Dp=200°
- WW3: Hs=3.0m, Tp=14.1s, Dp=200.5°
- **Implementação real**: Score = 0.1 + 0.5/40 = 0.1125 → **ACEITA**
- **Documentação**: Com Hs incluído, seria rejeitado com maior custo

---

## Seção 4: Fluxo Completo do Matching Real

### Pseudocódigo Completo

```python
# Carregar dados
df_sar = load_sar_csv(ref_id)
df_ww3_list = load_ww3_csvs(ref_id)  # Múltiplos timesteps

# === STEP 1: Matching Temporal ===
for df_ww3 in df_ww3_list:
    time_diff = abs(df_sar['obs_time'] - df_ww3['obs_time'])
    if time_diff < 1.0 hours:  # MAX_TIME_DIFF
        best_ww3 = df_ww3
        break

if time_diff > 1.0 hours:
    SKIP THIS CASE

# === STEP 2: Extrair Partições ===
sar_partitions = extract_partitions(df_sar)  # P1, P2, P3
ww3_partitions = extract_partitions(df_ww3)  # P1, P2, P3

# === STEP 3: Matching de Partições (GREEDY) ===
matches = []

for sar_pnum in [1, 2, 3]:
    sar_p = sar_partitions[sar_pnum-1]
    
    # Filter 1: Rejeita se Tp < 12s (wind sea)
    if sar_p['Tp'] < 12.0:
        continue
    
    # Filter 2: Rejeita se energia < 1% do total
    if sar_p['energy'] < 0.01 * total_energy:
        continue
    
    best_ww3_p = None
    min_score = inf
    
    # Busca em WW3
    for ww3_pnum, ww3_p in enumerate(ww3_partitions, 1):
        # Hard reject se fora de tolerância
        if abs(sar_p['Tp'] - ww3_p['Tp']) > 1.0:
            continue
        if angular_diff(sar_p['Dp'], ww3_p['Dp']) > 15.0:
            continue
        
        # Score para desempate
        score = abs(sar_p['Tp'] - ww3_p['Tp']) + \
                angular_diff(sar_p['Dp'], ww3_p['Dp']) / 40.0
        
        if score < min_score:
            min_score = score
            best_ww3_p = (ww3_pnum, ww3_p)
    
    if best_ww3_p:
        matches.append({
            'ref_id': ref_id,
            'sar_partition': sar_pnum,
            'ww3_partition': best_ww3_p[0],
            'sar_Hs': sar_p['Hs'],
            'sar_Tp': sar_p['Tp'],
            'sar_Dp': sar_p['Dp'],
            'ww3_Hs': best_ww3_p[1]['Hs'],
            'ww3_Tp': best_ww3_p[1]['Tp'],
            'ww3_Dp': best_ww3_p[1]['Dp'],
            'swap': sar_pnum != best_ww3_p[0],  # Se P2→P1, etc
            'cost': min_score
        })

# === OUTPUT ===
# Salva matches em CSV
# Casos sem nenhum match: 0 linhas no CSV
# Casos com 1 match: 1 linha no CSV
# Casos com 2 matches: 2 linhas no CSV
# Casos com 3 matches: 3 linhas no CSV
```

---

## Seção 5: Estatísticas de Saída

### 5.1 Arquivo `matched_partitions.csv`

Colunas:
- `reference_id`: Identificador único
- `obs_partition`: Número da partição SAR (1, 2, ou 3)
- `ww3_partition`: Número da partição WW3 (1, 2, ou 3)
- `swap`: Boolean indicando se foi reordenação (P1→P2, etc.)
- `match_quality`: "strict" ou "relaxed"
- `obs_Hs`, `ww3_Hs`, `dHs`: Alturas significativas
- `obs_Tp`, `ww3_Tp`, `dTp`: Períodos de pico
- `obs_Dp`, `ww3_Dp`, `dDp`: Direções de pico
- `cost`: Score de matching

### 5.2 Arquivo `case_level_stats.csv`

Colunas:
- `reference_id`: Identificador
- `n_obs`: Número de partições obs com Hs > 0.05m
- `n_ww3`: Número de partições WW3 com Hs > 0.05m
- `n_matched`: Número de pares successfully matched
- `n_unmat_obs`: Partições obs sem pair (missed detections)
- `n_unmat_ww3`: Partições WW3 sem pair (false alarms)

---

## Seção 6: Recomendações para Correção da Documentação

1. **Remover referência ao Hungarian algorithm** em matching.md
2. **Descrever o algoritmo real**: Greedy, sequencial por ordem P1→P2→P3
3. **Atualizar tolerâncias**: Tp=1.0s, Dp=15.0° (não as fórmulas complexas)
4. **Esclarecer o wind sea filter**: Tp >= 12.0s para SAR
5. **Remover matriz Hs do custo** (não é usada)
6. **Documentar que matching é determinístico**, não probabilístico
7. **Exemplos práticos** mostrando casos de rejeição por tolerância

---

## Conclusão

A **documentação anterior é conceitualmente correta mas não descreve a implementação real**. O código atual é mais simples e rigoroso, mas menos flexível. Isso pode resultar em:

- ✅ **Vantagens**: Matching determinístico, fácil de debugar, tolerâncias bem definidas
- ❌ **Desvantagens**: Pode deixar partições sem pair, não otimiza globalmente

Uma **implementação futura** poderia retornar ao Hungarian algorithm se houver necessidade de matching mais sofisticado.
