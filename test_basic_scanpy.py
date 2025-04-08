import os
import gc
import sys
import scanpy as sc
import polars as pl
import matplotlib.pyplot as plt

work_dir = 'projects/sc-benchmarking'
data_dir = 'scratch/single-cell/SEAAD'
sys.path.append(work_dir)

from utils_local import TimerCollection, system_info

system_info()
size = '400K'
timers = TimerCollection(silent=False)

# # Note: Load times should not be considered when using $SCRATCH disk 

# with timers('Load data (10X mtx)'):
#     data = sc.read_10x_mtx(f'{data_dir}/SEAAD_raw_{size}')
# del data; gc.collect()

# with timers('Load data (h5)'):
#     data = sc.read_10x_h5(f'{data_dir}/SEAAD_raw_{size}.h5')
# del data; gc.collect()

# Note: h5ad file contains additional metadata columns vs above

with timers('Load data'):
    data = sc.read_h5ad(f'{data_dir}/SEAAD_raw_{size}.h5ad')

# Note: QC filters are matched across libraries for timing, then 
# standardized by filtering to single_cell.py QC cells, not timed 

with timers('Quality control'):
    data.var['mt'] = data.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], inplace=True, log1p=True)
    sc.pp.filter_cells(data, min_genes=100)
    data = data[data.obs['pct_counts_mt'] <= 5.0]

print(f'cells: {data.shape[0]}, genes: {data.shape[1]}')

with timers('Doublet detection'):
    sc.pp.scrublet(data, batch_key='sample')

data = data[data.obs['passed_QC_tmp']]
print(f'cells: {data.shape[0]}, genes: {data.shape[1]}')

with timers('Normalization'):
    sc.pp.normalize_total(data)
    sc.pp.log1p(data)

with timers('Feature selection'):
    sc.pp.highly_variable_genes(data, n_top_genes=2000, batch_key='sample')

with timers('PCA'):
    sc.tl.pca(data)

with timers('Neighor graph'):
    sc.pp.neighbors(data)

with timers('Embedding'):
    sc.tl.umap(data)

with timers('Clustering (3 resoltions)'):
    for res in [0.5, 1, 2]:
        sc.tl.leiden(
            data, flavor='igraph', key_added=f'leiden_res_{res:4.2f}',
            resolution=res)

print(f"leiden_res_0.50: {len(data.obs['leiden_res_0.50'].unique())}")
print(f"leiden_res_1.00: {len(data.obs['leiden_res_1.00'].unique())}")
print(f"leiden_res_2.00: {len(data.obs['leiden_res_2.00'].unique())}")

with timers('Plot embeddings'):
    sc.pl.umap(data, color=['leiden_res_1.00'])
    plt.savefig(f'{work_dir}/figures/scanpy_embedding_cluster_{size}.png',
                dpi=300, bbox_inches='tight', pad_inches='layout')
    
with timers('Find markers'):
    sc.tl.rank_genes_groups(data, groupby='leiden_res_1.00', method='wilcoxon')

timers.print_summary(sort=False)
timers_df = timers\
    .to_dataframe(sort=False, unit='s')\
    .with_columns(pl.lit('test_basic_sc').alias('test'),
                  pl.lit(size).alias('size'))

print(timers_df)

timers_df.write_csv(f'{work_dir}/output/test_basic_scanpy_{size}.csv')

'''
--- Timing Summary ---
Load data took 719ms 367µs (0.4%)
Quality control took 3s 386ms (1.9%)
Doublet detection took 1m 4s (36.7%)
Normalization took 645ms 541µs (0.4%)
Feature selection took 2s 847ms (1.6%)
PCA took 2s 412ms (1.4%)
Neighor graph took 28s 563ms (16.2%)
Embedding took 18s 540ms (10.5%)
Clustering (3 resoltions) took 1s 189ms (0.7%)
Plot embeddings took 1s 154ms (0.7%)
Find markers took 51s 765ms (29.4%)

Total time: 2m 55s
'''