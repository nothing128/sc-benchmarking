import gc
import sys
import scanpy as sc
import polars as pl
import matplotlib.pyplot as plt

work_dir = 'projects/sc-benchmarking'
data_dir = 'scratch/single-cell/SEAAD'
sys.path.append(work_dir)

from utils_local import TimerMemoryCollection, system_info

system_info()
size_options = ['20K', '400K', '1M']
size_options = ['20K']
for size in size_options:

    timers = TimerMemoryCollection(silent=False)

    # with timers('Load data (10X mtx)'):
    #     data = sc.read_10x_mtx(f'{data_dir}/SEAAD_raw_{size}')
    # del data; gc.collect()

    # with timers('Load data (h5)'):
    #     data = sc.read_10x_h5(f'{data_dir}/SEAAD_raw_{size}.h5')
    # del data; gc.collect()

    # Note: Loading is much slower from $SCRATCH disk

    with timers('Load data (h5ad/rds)'):
        data = sc.read_h5ad(f'{data_dir}/SEAAD_raw_{size}.h5ad')

    # Note: QC filters are matched across libraries for timing, then 
    # standardized by filtering to single_cell.py QC cells, not timed 

    with timers('Quality control'):
        data.var['mt'] = data.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(
            data, qc_vars=['mt'], inplace=True, log1p=True)
        sc.pp.filter_cells(data, min_genes=100, copy=False)
        sc.pp.filter_genes(data, min_cells=3, copy=True)

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

    with timers('Neighbor graph'):
        sc.pp.neighbors(data)

    with timers('Embedding'):
        sc.tl.umap(data)

    #TODO: The number of clusters needs to match across libraries

    with timers('Clustering (3 resolutions)'):
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
        sc.tl.rank_genes_groups(
            data, groupby='leiden_res_1.00', method='wilcoxon')

    timers.print_summary(sort=False)
    timers_df = timers\
        .to_dataframe(sort=False, unit='s')\
        .with_columns(pl.lit('test_basic_scanpy').alias('test'),
                    pl.lit(size).alias('size'))

    print(timers_df)
    timers_df.write_csv(f'{work_dir}/output/test_basic_scanpy_{size}.csv')

    del timers, timers_df, data; gc.collect()

'''
--- System Information ---
Node: nia0036.scinet.local
CPU: 40 physical cores, 80 logical cores
Memory: 170.6 GB available / 188.6 GB total

--- Timing Summary 20K ---
Load data (h5ad/rds) took 703ms 479µs (0.4%)
Quality control took 5s 288ms (2.7%)
Doublet detection took 1m 22s (42.0%)
Normalization took 653ms 406µs (0.3%)
Feature selection took 2s 879ms (1.5%)
PCA took 2s 435ms (1.2%)
Neighor graph took 29s 321ms (15.0%)
Embedding took 18s 275ms (9.4%)
Clustering (3 resoltions) took 1s 139ms (0.6%)
Plot embeddings took 1s 536ms (0.8%)
Find markers took 51s 102ms (26.2%)

Total time: 3m 15s

--- Timing Summary 400K ---
Load data (h5ad/rds) took 9s 807ms (0.1%)
Quality control took 1m 2s (0.8%)
Doublet detection took 50m 27s (37.5%)
Normalization took 13s 34ms (0.2%)
Feature selection took 35s 798ms (0.4%)
PCA took 44s 335ms (0.5%)
Neighor graph took 56s 655ms (0.7%)
Embedding took 8m 10s (6.1%)
Clustering (3 resoltions) took 45s 215ms (0.6%)
Plot embeddings took 4s 944ms (0.1%)
Find markers took 1h 11m (53.0%)

Total time: 2h 14m
'''