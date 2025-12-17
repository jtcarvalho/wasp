"""
Template para Particionar Espectros de Bóia NDBC

Este é um template inicial para você adaptar aos seus dados de bóia.
Ajuste conforme necessário para o formato específico dos seus arquivos.
"""

import os
import pandas as pd
import numpy as np

# Importar do pacote wasp
from wasp.wave_params import calculate_wave_parameters
from wasp.partition import partition_spectrum


# ============================================================================
# CONFIGURAÇÃO
# ============================================================================

# Diretórios
INPUT_DIR = '/caminho/para/seus/dados/ndbc'      # ← AJUSTE AQUI
OUTPUT_DIR = '../output/ndbc'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Parâmetros de particionamento
ENERGY_THRESHOLD = 1e-6    # Ajuste baseado na energia típica da sua bóia
MAX_PARTITIONS = 5
MIN_PARTITION_POINTS = 5


# ============================================================================
# FUNÇÕES AUXILIARES
# ============================================================================

def load_ndbc_spectrum(filepath):
    """
    Carrega espectro de bóia NDBC.
    
    VOCÊ PRECISA IMPLEMENTAR ESTA FUNÇÃO baseado no formato dos seus dados.
    
    Returns
    -------
    E : ndarray (NF, ND)
        Espectro 2D em m²/Hz/rad
    freq : ndarray (NF,)
        Frequências em Hz
    dirs : ndarray (ND,)
        Direções em graus (convenção oceanográfica - para onde vai)
    metadata : dict
        Informações adicionais (timestamp, lat, lon, etc.)
    """
    
    # EXEMPLO - AJUSTE PARA SEU FORMATO
    # Opção 1: Se seus dados estão em CSV
    # df = pd.read_csv(filepath)
    # E = df.values.reshape(NF, ND)
    
    # Opção 2: Se seus dados estão em NetCDF
    # import xarray as xr
    # ds = xr.open_dataset(filepath)
    # E = ds['energy'].values  # [freq x dir]
    # freq = ds['frequency'].values
    # dirs = ds['direction'].values
    
    # IMPORTANTE: Verificar unidades e convenções
    # - Energia deve estar em m²/Hz/rad
    # - Se estiver em m²/Hz/deg, multiplique por π/180
    # - Direção deve ser oceanográfica (para onde vai)
    # - Se estiver meteorológica (de onde vem), converta:
    #   from wasp.utils import convert_meteorological_to_oceanographic
    #   dirs = convert_meteorological_to_oceanographic(dirs)
    
    raise NotImplementedError(
        "Você precisa implementar a função load_ndbc_spectrum() "
        "baseado no formato dos seus dados de bóia."
    )
    
    # return E, freq, dirs, metadata


# ============================================================================
# PROCESSAMENTO PRINCIPAL
# ============================================================================

def process_ndbc_station(station_file):
    """
    Processa um arquivo de bóia e retorna partições.
    """
    print(f"\nProcessando: {station_file}")
    
    # 1. Carregar espectro
    E, freq, dirs, metadata = load_ndbc_spectrum(station_file)
    
    print(f"  Espectro: {E.shape[0]} freq x {E.shape[1]} dir")
    print(f"  Freq range: {freq.min():.3f} - {freq.max():.3f} Hz")
    print(f"  Dir range: {dirs.min():.1f} - {dirs.max():.1f} deg")
    print(f"  Energia total: {np.sum(E):.2e} m²")
    
    # 2. Aplicar particionamento
    partitions = partition_spectrum(
        E, freq, dirs,
        energy_threshold=ENERGY_THRESHOLD,
        max_partitions=MAX_PARTITIONS,
        min_partition_points=MIN_PARTITION_POINTS
    )
    
    print(f"  → {len(partitions)} partições identificadas")
    
    # 3. Calcular parâmetros de cada partição
    results = []
    
    for i, partition in enumerate(partitions):
        params = calculate_wave_parameters(partition, freq, dirs)
        
        result = {
            'timestamp': metadata.get('timestamp', ''),
            'station_id': metadata.get('station_id', ''),
            'lat': metadata.get('lat', np.nan),
            'lon': metadata.get('lon', np.nan),
            'partition': i + 1,
            'Hs': params['Hs'],
            'Tp': params['Tp'],
            'Dp': params['Dp'],
            'fp': params['fp'],
            'Tm': params['Tm'],
            'E_total': params['E'],
        }
        results.append(result)
        
        print(f"    P{i+1}: Hs={params['Hs']:.2f}m, Tp={params['Tp']:.1f}s, Dp={params['Dp']:.0f}°")
    
    return results


# ============================================================================
# MAIN
# ============================================================================

def main():
    """
    Processa todos os arquivos de bóia no diretório de entrada.
    """
    print("="*70)
    print("PARTICIONAMENTO DE ESPECTROS NDBC")
    print("="*70)
    
    # Listar arquivos de entrada
    # AJUSTE o padrão de busca conforme seus arquivos
    import glob
    buoy_files = glob.glob(os.path.join(INPUT_DIR, '*.nc'))  # ou *.csv, *.txt, etc.
    
    if not buoy_files:
        print(f"\n⚠️  Nenhum arquivo encontrado em {INPUT_DIR}")
        print("   Ajuste INPUT_DIR e o padrão de busca no código.")
        return
    
    print(f"\nEncontrados {len(buoy_files)} arquivos")
    
    # Processar cada arquivo
    all_results = []
    
    for buoy_file in buoy_files:
        try:
            results = process_ndbc_station(buoy_file)
            all_results.extend(results)
        except Exception as e:
            print(f"  ❌ Erro processando {buoy_file}: {e}")
            continue
    
    # Salvar resultados
    if all_results:
        df = pd.DataFrame(all_results)
        output_file = os.path.join(OUTPUT_DIR, 'ndbc_partitions.csv')
        df.to_csv(output_file, index=False)
        print(f"\n✓ Resultados salvos em: {output_file}")
        print(f"  Total de {len(all_results)} partições processadas")
    else:
        print("\n⚠️  Nenhum resultado para salvar")


if __name__ == '__main__':
    main()
