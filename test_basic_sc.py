import gc
import sys
import polars as pl
from single_cell import SingleCell

work_dir = 'projects/sc-benchmarking'
data_dir = 'scratch/single-cell/SEAAD'
sys.path.append(work_dir)

from utils_local import TimerCollection, system_info

system_info()

num_threads = -1

for size in ['20K', '400K', '1M']:

    timers = TimerCollection(silent=False)

    # with timers('Load data (10X mtx)'):
    #     data = SingleCell(f'{data_dir}/SEAAD_raw_{size}/matrix.mtx.gz')
    # del data; gc.collect()

    # with timers('Load data (h5)'):
    #     data = SingleCell(f'{data_dir}/SEAAD_raw_{size}.h5')
    # del data; gc.collect()

    # with timers('Load data (rds)'):
    #     data = SingleCell(f'{data_dir}/SEAAD_raw_{size}.rds')
    # del data; gc.collect()

    # with timers('Load data (h5Seurat)'):
    #     data = SingleCell(f'{data_dir}/SEAAD_raw_{size}.h5Seurat')  
    #     del data; gc.collect()

    # Note: Loading is much slower from $SCRATCH disk
    
    with timers('Load data'):
        data = SingleCell(
            f'{data_dir}/SEAAD_raw_{size}.h5ad',
            num_threads=num_threads)

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
            allow_float=True,
            num_threads=num_threads)
        
    print(f'cells: {data.shape[0]}, genes: {data.shape[1]}')

    with timers('Doublet detection'):
        data = data.find_doublets(
            batch_column='sample',
            num_threads=num_threads)

    data = data.filter_obs(pl.col('passed_QC_tmp'))
    print(f'cells: {data.shape[0]}, genes: {data.shape[1]}')

    with timers('Feature selection'):
        data = data.hvg(
            allow_float=True,
            num_threads=num_threads)

    with timers('Normalization'):
        data = data.normalize(
            allow_float=True,
            num_threads=num_threads)

    with timers('PCA'):
        data = data.PCA(num_threads=num_threads)

    with timers('Neighbor graph'):
        data = data.neighbors(num_threads=num_threads)
        data = data.shared_neighbors(num_threads=num_threads)

    #TODO: The number of clusters needs to match across libraries

    with timers('Clustering (3 resolutions)'):
        data = data.cluster(
            resolution=[1, 0.5, 2],
            num_threads=num_threads)

    print(f'cluster_1: {len(data.obs['cluster_1'].unique())}')
    print(f'cluster_2: {len(data.obs['cluster_2'].unique())}')
    print(f'cluster_3: {len(data.obs['cluster_3'].unique())}')

    data = data.cast_obs({'cluster_1': pl.String})

    with timers('Embedding'):
        data = data.embed(
            num_threads=num_threads)

    with timers('Plot embeddings'):
        data.plot_embedding(
            'cluster_1',
            f'{work_dir}/figures/sc_embedding_cluster_{size}.png')

    with timers('Find markers'):
        markers = data.find_markers(
            'cluster_1',
            num_threads=num_threads)

    timers.print_summary(sort=False)
    timers_df = timers\
        .to_dataframe(sort=False, unit='s')\
        .with_columns(pl.lit('test_basic_sc').alias('test'),
                    pl.lit(size).alias('size'))

    print(timers_df)
    timers_df.write_csv(
        f'{work_dir}/output/test_basic_sc_single_thread_{size}.csv')

    del timers, timers_df, data; gc.collect()

'''
--- System Information ---
Node: nia0036.scinet.local
CPU: 40 physical cores, 80 logical cores
Memory: 170.6 GB available / 188.6 GB total

--- Timing Summary 20K ---
Load data (h5ad/rds) took 10s 441ms (22.0%)
Quality control took 438ms 695µs (0.9%)
Doublet detection took 9s 759ms (20.6%)
Feature selection took 250ms 396µs (0.5%)
Normalization took 98ms 580µs (0.2%)
PCA took 24s 116ms (50.8%)
Neighbor graph took 230ms 84µs (0.5%)
Clustering (3 resolutions) took 108ms 114µs (0.2%)
Embedding took 148ms 844µs (0.3%)
Plot embeddings took 1s 768ms (3.7%)
Find markers took 102ms 201µs (0.2%)

Total time: 47s 462ms

--- Timing Summary 400K ---
Load data (h5ad/rds) took 16s 404ms (19.0%)
Quality control took 5s 490ms (6.3%)
Doublet detection took 9s 988ms (11.5%)
Feature selection took 618ms 243µs (0.7%)
Normalization took 570ms 126µs (0.7%)
PCA took 34s 203ms (39.5%)
Neighbor graph took 5s 714ms (6.6%)
Clustering (3 resolutions) took 4s 474ms (5.2%)
Embedding took 3s 698ms (4.3%)
Plot embeddings took 5s 69ms (5.9%)
Find markers took 300ms 986µs (0.3%)

Total time: 1m 26s

--- Timing Summary 1M ---
Load data (h5ad/rds) took 50s 932ms (29.7%)
Quality control took 12s 550ms (7.3%)
Doublet detection took 19s 255ms (11.2%)
Feature selection took 1s 460ms (0.9%)
Normalization took 1s 464ms (0.9%)
PCA took 30s 932ms (18.0%)
Neighbor graph took 22s 135ms (12.9%)
Clustering (3 resolutions) took 14s 548ms (8.5%)
Embedding took 7s 829ms (4.6%)
Plot embeddings took 10s 48ms (5.9%)
Find markers took 484ms 85µs (0.3%)

Total time: 2m 51s

--- Timing Summary 20K single thread ---
Load data (h5ad/rds) took 2s 758ms (1.9%)
Quality control took 955ms 848µs (0.7%)
Doublet detection took 2m 9s (90.2%)
Feature selection took 815ms 28µs (0.6%)
Normalization took 1s 113ms (0.8%)
PCA took 2s 908ms (2.0%)
Neighbor graph took 779ms 234µs (0.5%)
Clustering (3 resolutions) took 296ms 994µs (0.2%)
Embedding took 1s 685ms (1.2%)
Plot embeddings took 2s 547ms (1.8%)
Find markers took 190ms 362µs (0.1%)


20K
Load data (10X mtx) took 8s 732ms (8.7%)
Load data (h5) took 10s 588ms (10.5%)
Load data (rds) took 18s 250ms (18.1%)

400K
Load data (10X mtx) took 2m 47s (31.8%)
Load data (h5) took 23s 147ms (4.4%)
Load data (rds) took 3m 47s (43.3%)
'''

