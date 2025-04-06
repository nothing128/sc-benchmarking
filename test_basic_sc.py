import os
import gc
import sys
import polars as pl
from single_cell import SingleCell

work_dir = 'projects/def-wainberg/karbabi/sc-benchmarking'
data_dir = 'single-cell/SEAAD/subsampled'
sys.path.append(work_dir)

from utils_local import TimerCollection, system_info

system_info()
size = '20K'
timers = TimerCollection(silent=False)

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
    data = SingleCell(f'{data_dir}/SEAAD_raw_{size}.h5ad')

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

print(f'cell_type_1: {len(data.obs['cell_type_1'].unique())}')
print(f'cell_type_2: {len(data.obs['cell_type_2'].unique())}')
print(f'cell_type_3: {len(data.obs['cell_type_3'].unique())}')

data = data.cast_obs({'cell_type_1': pl.String})

with timers('Embedding'):
    data = data.embed()

with timers('Plot embeddings'):
    data.plot_embedding(
        'cell_type_1',
        f'{work_dir}/figures/sc_embedding_cluster_{size}.png')

with timers('Find markers'):
    markers = data.find_markers('cell_type_1')

timers.print_summary(sort=False)
timers_df = timers.to_dataframe()\
    .with_columns(pl.lit('test_basic_sc').alias('test'),
                  pl.lit(size).alias('size'))

'''
--- Timing Summary ---
Load data (h5ad) took 2s 197ms (4.7%)
Quality control took 286ms 565µs (0.6%)
Doublet detection took 9s 903ms (21.3%)
Feature selection took 247ms 793µs (0.5%)
Normalization took 44ms 560µs (0.1%)
PCA took 25s 57ms (53.8%)
Neighbor graph took 814ms 333µs (1.7%)
Clustering took 5s 238ms (11.3%)
Embedding took 738ms 239µs (1.6%)
Plot embeddings took 1s 937ms (4.2%)
Find markers took 79ms 789µs (0.2%)

Total time: 46s 545ms
'''

