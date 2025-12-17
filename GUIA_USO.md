# Guia de Uso do Pacote WASP

## 📦 Estrutura do Projeto

```
wasp/
├── pyproject.toml          # Configuração do pacote pip
├── LICENSE                 # Licença MIT
├── README.md              # Documentação principal
├── requirements.txt       # [mantido para referência]
├── environment.yml        # [mantido para referência]
│
├── src/                   # ⭐ SRC LAYOUT (melhor prática)
│   └── wasp/              # PACOTE PRINCIPAL (core reutilizável)
│       ├── __init__.py    # Exporta funções principais
│       ├── partition.py   # Algoritmo de particionamento
│       ├── wave_params.py # Cálculo de parâmetros de onda
│       ├── io_sar.py     # I/O para dados SAR
│       ├── io_ww3.py     # I/O para dados WW3
│       ├── plotting.py   # Funções de visualização
│       └── utils.py      # Utilidades gerais
│
├── examples/              # Scripts de exemplo
│   ├── README.md         # Documentação dos exemplos
│   ├── 01_partition_sar.py
│   ├── 02_partition_ww3.py
│   ├── 03_partition_ndbc.py  # [vazio - template para você]
│   └── 04_validate.py
│
├── scripts/               # [mantido - versão original]
│   ├── lib/              # [código original]
│   └── ...
│
├── data/                 # Seus dados
├── auxdata/             # Arquivos auxiliares
├── output/              # Resultados
└── notebooks/           # Jupyter notebooks
```

## 🚀 Como Usar

### 1. Instalação (uma vez)

```bash
cd /Users/jtakeo/googleDrive/myProjects/wasp
pip install -e .
```

O modo `-e` (editável) significa que você pode modificar o código em `wasp/` e as mudanças serão refletidas imediatamente, sem precisar reinstalar.

### 2. Usar como Biblioteca

Agora você pode importar `wasp` de qualquer lugar:

```python
# Em qualquer script Python
from wasp import partition_spectrum, calculate_wave_parameters
from wasp.io_sar import load_sar_spectrum
from wasp.io_ww3 import load_ww3_spectrum
from wasp.utils import spectrum1d_from_2d

# Use as funções normalmente
partitions = partition_spectrum(E, freq, dirs)
```

### 3. Aplicar para Bóia NDBC

Exemplo completo para dados de bóia:

```python
import numpy as np
from wasp import partition_spectrum, calculate_wave_parameters

# 1. Carregue seu espectro de bóia
# (substitua com seu código de leitura)
E_buoy = ...  # matriz [freq x dir] em m²/Hz/rad
freq = ...    # array de frequências [Hz]
dirs = ...    # array de direções [graus, convenção oceanográfica]

# 2. Configure os parâmetros
config = {
    'energy_threshold': 1e-6,    # Ajuste baseado na energia típica
    'max_partitions': 5,         # Máximo de partições
    'min_partition_points': 5,   # Pontos mínimos por partição
}

# 3. Aplique o particionamento
partitions = partition_spectrum(
    E_buoy, freq, dirs,
    **config
)

# 4. Calcule parâmetros para cada partição
results = []
for i, partition in enumerate(partitions):
    params = calculate_wave_parameters(partition, freq, dirs)
    
    result = {
        'partition_id': i + 1,
        'Hs': params['Hs'],      # Altura significativa [m]
        'Tp': params['Tp'],      # Período de pico [s]
        'Dp': params['Dp'],      # Direção de pico [graus]
        'fp': params['fp'],      # Frequência de pico [Hz]
        'E_total': params['E'],  # Energia total [m²]
    }
    results.append(result)
    
    print(f"Partição {i+1}: Hs={params['Hs']:.2f}m, Tp={params['Tp']:.1f}s, Dp={params['Dp']:.0f}°")

# 5. Salve os resultados
import pandas as pd
df = pd.DataFrame(results)
df.to_csv('ndbc_partitions.csv', index=False)
```

### 4. Personalizar Parâmetros

Os parâmetros principais que você pode ajustar:

```python
partitions = partition_spectrum(
    E, freq, dirs,
    
    # Threshold de energia mínima para identificar picos
    # Valores menores = mais sensível a pequenos sistemas
    # SAR: 1e-6 a 1e-5
    # WW3: 1e-6 a 1e-5
    # Bóias: depende da energia típica observada
    energy_threshold=1e-6,
    
    # Número máximo de partições a identificar
    max_partitions=5,
    
    # Número mínimo de pontos espectrais em uma partição
    # Valores maiores = partições mais "robustas"
    min_partition_points=5,
    
    # Fator de separação entre picos
    # Valores maiores = picos precisam estar mais separados
    min_separation_factor=1.5,
)
```

## 🔧 Modificar o Core

Se você precisar modificar a lógica de particionamento:

1. Edite os arquivos em `src/wasp/`
2. As mudanças são refletidas imediatamente (instalação editável)
3. Não precisa reinstalar o pacote

Por exemplo, para mudar `partition.py`:
- Abra [src/wasp/partition.py](src/wasp/partition.py)
- Faça suas modificações
- Teste com `from wasp import partition_spectrum`

> **📝 Nota sobre src/ layout:**  
> Usamos a estrutura `src/wasp/` (ao invés de `wasp/` na raiz) seguindo as melhores práticas modernas de Python. Isso evita imports acidentais do código não-instalado e garante que você sempre está testando o pacote instalado.

## 📝 Convenções Importantes

### Unidades de Energia Espectral
O core espera espectro em **m²/Hz/rad**:

- **SAR**: geralmente vem em m⁴ (k-spectrum), precisa converter
- **WW3**: já está em m²/Hz/rad ✓
- **Bóias**: verificar unidades (geralmente m²/Hz/deg, precisa converter)

### Convenção de Direção
Use **convenção oceanográfica** (para onde vai):

- **SAR**: geralmente vem em convenção meteorológica, precisa converter
- **WW3**: já está em oceanográfica ✓
- **Bóias**: verificar convenção (pode ser meteorológica)

Converter se necessário:
```python
from wasp.utils import convert_meteorological_to_oceanographic
dir_ocean = convert_meteorological_to_oceanographic(dir_met)
```

### Conversão SAR
Para SAR (número de onda → frequência):

```python
from wasp.utils import convert_sar_energy_units

# E_k: espectro em número de onda [m⁴]
# k: número de onda [rad/m]
# phi: direções [graus]
E_freq = convert_sar_energy_units(E_k, k, phi)
```

## 🧪 Testar Suas Modificações

```bash
# Criar um script de teste simples
cat > test_wasp.py << 'EOF'
import numpy as np
from wasp import partition_spectrum, calculate_wave_parameters

# Criar espectro sintético
NF, ND = 32, 36
freq = np.linspace(0.04, 0.5, NF)
dirs = np.linspace(0, 360, ND, endpoint=False)

# Pico Gaussiano simples
E = np.zeros((NF, ND))
E[10, 15] = 10.0  # Um pico em ~0.18 Hz, ~150°
E = np.exp(-((np.arange(NF)[:, None] - 10)**2 + (np.arange(ND)[None, :] - 15)**2) / 10)

# Testar particionamento
partitions = partition_spectrum(E, freq, dirs, energy_threshold=0.01, max_partitions=3)
print(f"Encontradas {len(partitions)} partições")

for i, p in enumerate(partitions):
    params = calculate_wave_parameters(p, freq, dirs)
    print(f"  P{i+1}: Hs={params['Hs']:.2f}m, Tp={params['Tp']:.1f}s, Dp={params['Dp']:.0f}°")
EOF

# Executar teste
python test_wasp.py
```

## 📚 Próximos Passos

1. ✅ Pacote instalado e funcionando
2. ⬜ Implementar `03_partition_ndbc.py` para suas bóias
3. ⬜ Ajustar thresholds conforme sua aplicação
4. ⬜ [Futuro] Adicionar testes automatizados
5. ⬜ [Futuro] Publicar no PyPI quando estiver pronto

## ❓ Dúvidas Comuns

**P: Meus scripts antigos em `scripts/` ainda funcionam?**  
R: Sim! Eles continuam usando `sys.path.insert` e funcionam independentemente.

**P: Preciso reinstalar após modificar o código?**  
R: Não! A instalação editável (`pip install -e .`) faz as mudanças serem refletidas automaticamente.

**P: Posso usar em outro projeto?**  
R: Sim! Com o pacote instalado, você pode importar de qualquer lugar:
```python
from wasp import partition_spectrum
```

**P: Como desinstalar?**  
R: `pip uninstall wasp`

**P: Os notebooks continuam funcionando?**  
R: Sim, mas você pode simplificar os imports usando o pacote ao invés de `sys.path.insert`.
