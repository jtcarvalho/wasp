# **WASP** - **WA**ve **S**pectra **P**artitioning

Watershed Algorithm for partitioning the ocean wave spectra from WW3 and SAR (Sentinel

<!--

**🔗 Companion Repository:** For analysis and validation of partitioned spectra, see [**HIVE** (Hierarchical Integration of Verified wavE partitions)](https://github.com/jtcarvalho/hive)

-->

## 📋 What is WASP?

WASP focuses exclusively on **spectral partitioning** - the process of separating ocean wave spectra into individual wave systems (partitions). Each partition represents a distinct wave system characterized by significant wave height (Hs), peak period (Tp), and direction (Dp).

**WASP handles:**

- ✅ Spectral partitioning using watershed algorithm
- ✅ Processing SAR (Sentinel) and WW3 model spectra
- ✅ Extracting wave parameters (Hs, Tp, Dp) for each partition

👉 **For analysis and validation**, use the repository [**HIVE**](https://github.com/jtcarvalho/hive)

## 🚀 Installation

### Método 1: Instalação Local (Desenvolvimento)

Instale o pacote em modo editável para desenvolvimento ou uso local:

```bash
# Clone o repositório
git clone https://github.com/jtcarvalho/wasp.git
cd wasp

# Instale em modo editável (recomendado)
pip install -e .
```

### Método 2: Ambiente Virtual Tradicional

```bash
# Clone o repositório
git clone https://github.com/jtcarvalho/wasp.git
cd wasp

# Crie ambiente virtual
python -m venv venv

# Ative o ambiente virtual
# No macOS/Linux:
source venv/bin/activate
# No Windows:
# venv\Scripts\activate

# Instale o pacote
pip install -e .
```

### Verificar Instalação

```bash
# Teste a importação
python -c "import wasp; print(f'WASP version: {wasp.__version__}')"

# Teste as funções principais
python -c "from wasp import partition_spectrum, calculate_wave_parameters; print('✓ Instalação bem-sucedida!')"
```

## 📦 Key Dependencies

- **NumPy >= 2.1.0** (required for `np.trapezoid`)
- pandas >= 2.2.0
- xarray >= 2024.11.0
- matplotlib >= 3.8.0
- scipy >= 1.14.0
- scikit-image >= 0.22.0
- netCDF4 >= 1.5.4

> ⚠️ **Importante:** NumPy < 2.1.0 causará erros pois `np.trapezoid` não está disponível.

## 💡 Uso Rápido

### Como Biblioteca Python

```python
import numpy as np
from wasp import partition_spectrum, calculate_wave_parameters

# Seu espectro 2D (freq x dir)
E = np.array(...)  # matriz de energia espectral [m²/Hz/rad]
freq = np.array(...)  # frequências [Hz]
dirs = np.array(...)  # direções [graus, convenção oceanográfica]

# Particionar o espectro
partitions = partition_spectrum(
    E, freq, dirs,
    energy_threshold=1e-6,
    max_partitions=3
)

# Calcular parâmetros de cada partição
for i, partition in enumerate(partitions):
    params = calculate_wave_parameters(partition, freq, dirs)
    print(f"Partição {i+1}:")
    print(f"  Hs = {params['Hs']:.2f} m")
    print(f"  Tp = {params['Tp']:.1f} s")
    print(f"  Dp = {params['Dp']:.1f} deg")
```

### Scripts de Exemplo

Veja a pasta [examples/](examples/) para scripts completos:

- **01_partition_sar.py**: Processar espectros SAR (Sentinel-1)
- **02_partition_ww3.py**: Processar espectros WaveWatch III
- **03_partition_ndbc.py**: Template para processar dados de bóia NDBC
- **04_validate.py**: Comparar e validar resultados SAR vs WW3

```bash
cd examples/
python 01_partition_sar.py
```
