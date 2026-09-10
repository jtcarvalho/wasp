# Resumo da Revisão de Metodologia WASP

> **Documento histórico.** Este arquivo descreve uma revisão anterior e é
> mantido para proveniência da pesquisa; não define a documentação atual do
> pacote. Consulte o [índice atual](INDEX.md) e a
> [implementação PCSPM atual](matching.md).
> Afirmações no presente e referências a scripts abaixo pertencem ao fluxo
> histórico e não descrevem o pacote atual.

**Data**: Fevereiro 2025  
**Objetivo**: Corrigir e atualizar documentação de matching de partições espectrais  
**Status**: ✅ Completo

---

## Documentos Criados

### 1. `new_methods.md` (Metodologia Completa Corrigida)

**Conteúdo**:
- Descrição completa do pipeline WASP
- Algoritmo de particionamento Hanson & Phillips (correto)
- Processamento específico CFOSAT (correto)
- **Matching temporal WW3** (correto)
- **Matching de partições com algoritmo REAL** (agora correto)
- Classificação de regime de swell
- Métricas de verificação
- Tabela resumida de diferenças

**Principais correções**:
- ❌ Hungarian algorithm → ✅ Greedy sequential matching
- ❌ Matriz de custo complexa → ✅ Score simples: tp_diff + dp_diff/40
- ❌ Tolerância Tp: max(2.0s, 0.20·Tp) → ✅ Fixa: 1.0s
- ❌ Tolerância Dp: 60° → ✅ 15.0°
- ❌ Tp mínimo SAR: 10s → ✅ 12.0s
- ❌ Ordem flexível → ✅ Ordem fixa: P1→P1, P2→P2, P3→P3
- ❌ Hs incluído no custo → ✅ Apenas Tp e Dp

---

### 2. `matching_implementation_vs_documentation.md` (Comparativo Detalhado)

**Conteúdo**:
- Seção 1: O que a documentação anterior afirma
- Seção 2: O que realmente está sendo feito (com código)
- Seção 3: Implicações práticas (exemplos reais)
- Seção 4: Fluxo completo com pseudocódigo
- Seção 5: Estrutura de saída dos dados
- Seção 6: Recomendações para correção

**Propósito**: Deixar cristalino o gap entre documentação e implementação

---

## Investigação Realizada

### Arquivos Revisados

1. **Biblioteca WASP**:
   - ✅ `src/wasp/partition.py` - Particionamento (Hanson & Phillips implementado corretamente)
   - ✅ `src/wasp/matching.py` - Funções SAR-NDBC matching
   - ✅ `src/wasp/wave_params.py` - Cálculo de parâmetros

2. **Scripts de Análise**:
   - ✅ Encontrado script de validação: `/bak_wasp_infos/examples/04_validate.py`
   - ✅ Este arquivo contém a implementação REAL de matching
   - ✅ Analisa múltiplas estratégias de matching em `test_partition_matching.py`

3. **Dados de Saída**:
   - ✅ Estrutura: `data/*/partition_system_analysis_*/`
   - ✅ Arquivos: `matched_partitions.csv` e `case_level_stats.csv`
   - ✅ Validou-se a estrutura contra o código

---

## Conclusões Principais

### O que WASP realmente é

**WASP Library** = Apenas particionamento de espectros 2D
- Não faz matching
- Não faz análise
- Input: Espectro 2D (frequência × direção)
- Output: Partições com Hs, Tp, Dp para cada sistema

### O que está fora de WASP

**Matching e Análise** = Scripts externos
- Carregam CSVs particionados
- Fazem matching temporal (WW3 mais próximo)
- Fazem matching de partições (greedy)
- Geram estatísticas e métricas
- Produzem visualizações

---

## Mudanças Principais Documentadas

### Algoritmo de Matching

| Característica | Documentação Anterior | Agora Correto |
|---|---|---|
| Tipo | Optimal (Hungarian) | Greedy (first compatible) |
| Ordem | Flexível | Fixa (P1→P1, P2→P2, P3→P3) |
| Matriz | 2D com pesos | Sequential linear search |
| Tolerância Tp | Híbrida (relativa+absoluta) | Absoluta: 1.0 s |
| Tolerância Dp | 60° (hard reject) | 15.0° (hard reject) |
| Filtro Tp SAR | 10.0 s | 12.0 s |
| Variáveis no custo | Tp, Dp, Hs | Apenas Tp, Dp |

### Tolerâncias Mais Rigorosas

- **Tp**: De até ±3.6s (18s wave) → **Fixa ±1.0s**
- **Dp**: De até ±60° → **Fixa ±15°**
- **Resultado**: Mais matches rejeitados, menos "missed detections" falsas

### Conceitual

- ✅ WASP é biblioteca pura de particionamento
- ✅ Matching é externo (não parte da "metodologia WASP")
- ✅ Implementação é determinística, não probabilística
- ✅ Resultado pode deixar partições sem pair (by design)

---

## Recomendações

### Curto Prazo (Documentação)
1. ✅ Criar `new_methods.md` - **FEITO**
2. ✅ Criar documento comparativo - **FEITO**
3. ❌ Deprecar `matching.md`, `methods.md` - **PENDENTE**
4. ❌ Atualizar README.md com clarificação - **PENDENTE**

### Médio Prazo (Código)
1. ❌ Mover scripts de análise para repositório principal WASP
2. ❌ Formalizar função `match_partitions()` na biblioteca
3. ❌ Adicionar tests unitários para matching
4. ❌ Documentar em docstrings Python

### Longo Prazo (Futuro)
1. ❌ Considerar Hungarian algorithm se flexibilidade for necessária
2. ❌ Implementar matching probabilístico (e.g., Bayesian)
3. ❌ Adicionar validação cruzada de estratégias de matching

---

## Arquivos de Referência

### Criados
- ✅ `/docs/new_methods.md` - Metodologia completa corrigida
- ✅ `/docs/matching_implementation_vs_documentation.md` - Comparativo

### Antigos (Para Referência)
- 📄 `/docs/matching.md` - Documentação anterior (com erro)
- 📄 `/docs/methods.md` - Documentação anterior (com erro)
- 📄 `/docs/methodology_partitioning_matching.md` - Documentação anterior (com erro)

### De Backup (Encontrados)
- 🔍 `/bak_wasp_infos/examples/04_validate.py` - Script de matching REAL
- 🔍 `/bak_wasp_infos/examples/test_partition_matching.py` - Testes de estratégias

---

## Exemplos de Uso da Nova Documentação

### Para Pesquisadores
- Ler `new_methods.md` seções 1-3: Entender o pipeline
- Ler seção 6: Entender exatamente como matching funciona
- Ler seção 8: Entender as métricas calculadas

### Para Publicação
- Usar `new_methods.md` como base para métodos em artigo
- Seção 9 lista explicitamente o que é diferente de documentação anterior
- Ser preciso sobre: tolerâncias, ordem de matching, filtros aplicados

### Para Reprodução
- Seção 4 em `matching_implementation_vs_documentation.md` tem pseudocódigo
- Pode-se implementar exatamente o mesmo algoritmo em outra linguagem
- Resultados serão idênticos

---

## Validação

✅ **Código revisado**: Todos os scripts Python relevantes  
✅ **Dados validados**: Estrutura de saída corresponde ao código  
✅ **Tolerâncias confirmadas**: Em `04_validate.py` linhas 58-60  
✅ **Algoritmo identificado**: Greedy matching, sequencial, P1→P2→P3  
✅ **Diferenças documentadas**: Tabela completa de comparação  
✅ **Exemplos práticos**: Casos reais de rejeição/aceitação  

---

## Próximos Passos

1. **Revisar** este sumário
2. **Validar** a precisão com código original
3. **Consultar** com autores originais se houver questões
4. **Integrar** `new_methods.md` na documentação oficial
5. **Deprecar** documentação anterior com nota sobre correções
6. **Versionar** mudanças (v1.0 → v1.1 de metodologia)
