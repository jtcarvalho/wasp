# Examples

Esta pasta contém scripts de exemplo que demonstram como usar o pacote WASP para diferentes fontes de dados.

## Pré-requisitos

Instale o pacote WASP em modo editável:

```bash
cd /caminho/para/wasp
pip install -e .
```

## Scripts Disponíveis

### 01_partition_sar.py
Particiona espectros 2D do SAR (Sentinel-1).

**Uso:**
```python
# Edite as configurações no início do script
python 01_partition_sar.py
```

**O que faz:**
- Lê arquivos NetCDF do SAR
- Converte de número de onda para frequência
- Aplica algoritmo de particionamento (watershed)
- Salva resultados em CSV com Hs, Tp, Dp por partição

### 02_partition_ww3.py
Particiona espectros 2D do modelo WaveWatch III.

**Uso:**
```python
# Edite as configurações no início do script
python 02_partition_ww3.py
```

**O que faz:**
- Lê arquivos NetCDF do WW3
- Encontra tempo mais próximo para cada ponto SAR
- Aplica particionamento
- Salva resultados em CSV

### 03_partition_ndbc.py
Particiona dados de bóia NDBC (vazio - para implementar).

**Uso futuro:**
```python
from wasp.partition import partition_spectrum
from wasp.wave_params import calculate_wave_parameters

# Carregue seu espectro de bóia
E_buoy = ...  # matriz 2D [freq x dir]
freq = ...    # array de frequências
dirs = ...    # array de direções

# Aplique particionamento
partitions = partition_spectrum(
    E_buoy, freq, dirs,
    energy_threshold=1e-6,  # ajuste conforme necessário
    max_partitions=3
)

# Calcule parâmetros de onda
for i, partition in enumerate(partitions):
    params = calculate_wave_parameters(partition, freq, dirs)
    print(f"Partição {i+1}:")
    print(f"  Hs = {params['Hs']:.2f} m")
    print(f"  Tp = {params['Tp']:.1f} s")
    print(f"  Dp = {params['Dp']:.1f} deg")
```

### 04_validate.py
Compara resultados SAR vs WW3.

**Uso:**
```python
# Edite as configurações no início do script
python 04_validate.py
```

**O que faz:**
- Carrega CSVs de particionamento
- Faz matching de partições similares
- Gera scatter plots e estatísticas
- Salva figuras e relatórios

## Estrutura de Dados

### Formato de Entrada (SAR/WW3)
- **Espectro 2D**: matriz `E[freq, dir]` com densidade espectral de energia
- **Frequências**: array 1D com bins de frequência (Hz)
- **Direções**: array 1D com bins de direção (graus, convenção oceanográfica)

### Formato de Saída (CSV)
Cada linha representa um ponto/tempo com até N partições:

```csv
lon,lat,time,Hs,Tp,Dp,Hs_p1,Tp_p1,Dp_p1,Hs_p2,Tp_p2,Dp_p2,...
```

## Personalizando o Particionamento

O core do particionamento aceita parâmetros configuráveis:

```python
from wasp.partition import partition_spectrum

partitions = partition_spectrum(
    E, freq, dirs,
    energy_threshold=1e-6,    # Threshold mínimo de energia
    max_partitions=5,         # Número máximo de partições
    min_partition_points=5,   # Pontos mínimos por partição
    min_separation_factor=1.5 # Fator de separação de picos
)
```

Ajuste esses parâmetros conforme sua aplicação:
- **SAR**: geralmente `energy_threshold=1e-6` a `1e-5`
- **WW3**: pode variar com resolução espectral
- **Bóias**: ajuste baseado em energia típica observada

## Aplicação para Bóias NDBC

Para aplicar em dados de bóia NDBC, você precisará:

1. Carregar espectro 2D da bóia (freq x dir)
2. Garantir que as unidades sejam m²/Hz/rad
3. Direções em convenção oceanográfica (para onde vai)
4. Aplicar `partition_spectrum()` com thresholds apropriados
5. Usar `calculate_wave_parameters()` para cada partição

Veja o exemplo comentado em `03_partition_ndbc.py` para um template inicial.
