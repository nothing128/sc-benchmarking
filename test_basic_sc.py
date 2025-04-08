import os
import gc
import sys
import polars as pl
from single_cell import SingleCell

work_dir = 'projects/sc-benchmarking'
data_dir = 'scratch/single-cell/SEAAD'
sys.path.append(work_dir)

from utils_local import TimerCollection, system_info

system_info()
size = '20K'
timers = TimerCollection(silent=False)

# # Note: Load times should not be considered when using $SCRATCH disk 

# with timers('Load data (10X mtx)'):
#     data = SingleCell(f'{data_dir}/SEAAD_raw_{size}/matrix.mtx.gz')
# del data; gc.collect()

# with timers('Load data (h5)'):
#     data = SingleCell(f'{data_dir}/SEAAD_raw_{size}.h5')
# del data; gc.collect()

# # Note: rds, h5ad files contain additional metadata columns

# with timers('Load data (rds)'):
#     data = SingleCell(f'{data_dir}/SEAAD_raw_{size}.rds')
# del data; gc.collect()

with timers('Load data'):
    data = SingleCell(
        f'{data_dir}/SEAAD_raw_{size}.h5ad')

print(f'X num_threads: {data.X._num_threads}')

# Note: QC filters are matched across libraries for timing, then 
# standardized by filtering to single_cell.py QC cells, not timed 

with timers('Quality control'):
    data = data.qc(
        subset=True,
        max_mito_fraction=0.05,
        min_genes=100,
        nonzero_MALAT1=False,
        remove_doublets=False,
        allow_float=True)
    
print(f'cells: {data.shape[0]}, genes: {data.shape[1]}')

with timers('Doublet detection'):
    data = data.find_doublets(batch_column='sample')

data = data.filter_obs(pl.col('passed_QC_tmp'))
print(f'cells: {data.shape[0]}, genes: {data.shape[1]}')

with timers('Feature selection'):
    data = data.hvg(allow_float=True)

with timers('Normalization'):
    data = data.normalize(allow_float=True)

with timers('PCA'):
    data = data.PCA()

with timers('Neighbor graph'):
    data = data.neighbors()
    data = data.shared_neighbors()

with timers('Clustering (3 resolutions)'):
    data = data.cluster(resolution=[1, 0.5, 2])

print(f'cluster_1: {len(data.obs['cluster_1'].unique())}')
print(f'cluster_2: {len(data.obs['cluster_2'].unique())}')
print(f'cluster_3: {len(data.obs['cluster_3'].unique())}')

data = data.cast_obs({'cluster_1': pl.String})

with timers('Embedding'):
    data = data.embed()

with timers('Plot embeddings'):
    data.plot_embedding(
        'cluster_1',
        f'{work_dir}/figures/sc_embedding_cluster_{size}.png')

with timers('Find markers'):
    markers = data.find_markers('cluster_1')

timers.print_summary(sort=False)
timers_df = timers\
    .to_dataframe(sort=False, unit='s')\
    .with_columns(pl.lit('test_basic_sc').alias('test'),
                  pl.lit(size).alias('size'))

print(timers_df)

timers_df.write_csv(f'{work_dir}/output/test_basic_sc_{size}.csv')

'''
--- Timing Summary ---
Load data took 5s 513ms (14.1%)
Quality control took 891ms 254µs (2.3%)
Doublet detection took 15s 870ms (40.6%)
Feature selection took 890ms 739µs (2.3%)
Normalization took 165ms 935µs (0.4%)
PCA took 8s 511ms (21.8%)
Neighbor graph took 307ms 762µs (0.8%)
Clustering (3 resolutions) took 155ms 735µs (0.4%)
Embedding took 1s 121ms (2.9%)
Plot embeddings took 5s 508ms (14.1%)
Find markers took 179ms 444µs (0.5%)

Total time: 39s 116ms
'''

