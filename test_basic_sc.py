import gc
import sys
import polars as pl  # type: ignore
from single_cell import SingleCell  # type: ignore

work_dir = 'projects/sc-benchmarking'
data_dir = 'single-cell/SEAAD/subsampled'

sys.path.append(work_dir)
from utils_local import TimerMemoryCollection, system_info

num_threads = int(sys.argv[1])
subset = bool(sys.argv[2])
size = str(sys.argv[3])
print('--- Params ---')
print(f'{size=}, {num_threads=}, {subset=}')

system_info()
timers = TimerMemoryCollection(silent=False)

# Note: Loading is much slower from $SCRATCH disk
# TODO: Temporarily setting `num_threads=1` for loading until shared 
# memory benchmarking is possible 

with timers('Load data'):
    data = SingleCell(f'{data_dir}/SEAAD_raw_{size}.h5ad', num_threads=1)
    data = data.set_num_threads(num_threads)

with timers('Quality control'):
    data = data.qc(
        subset=subset,
        remove_doublets=False,
        allow_float=True,
        verbose=False,
        num_threads=num_threads)

with timers('Doublet detection'):
    data = data.find_doublets(batch_column='sample', num_threads=num_threads)
        
# Not timed 
data = data.filter_obs(pl.col('doublet').not_())

with timers('Feature selection'):
    data = data.hvg(num_threads=num_threads)

with timers('Normalization'):
    data = data.normalize(num_threads=num_threads)

with timers('PCA'):
    data = data.PCA(num_threads=num_threads)

# Not timed
if not subset:
    data = data.filter_obs(pl.col('passed_QC'))

with timers('Neighbor graph'):
    data = data.neighbors(num_threads=num_threads)  
    data = data.shared_neighbors(num_threads=num_threads)  

# TODO: The number of clusters needs to match across libraries

with timers('Clustering (3 resolutions)'):
    data = data.cluster(resolution=[1, 0.5, 2], num_threads=num_threads)

print(f'cluster_0: {len(data.obs['cluster_0'].unique())}')
print(f'cluster_1: {len(data.obs['cluster_1'].unique())}')
print(f'cluster_2: {len(data.obs['cluster_2'].unique())}')

with timers('Embedding'):
    data = data.embed(num_threads=num_threads)

with timers('Plot embeddings'):
    data.plot_embedding(
        'cluster_0', f'{work_dir}/figures/sc_embedding_cluster_{size}.png')

with timers('Find markers'):
    markers = data.find_markers('cluster_0', num_threads=num_threads)

timers.print_summary(sort=False)

df = timers.to_dataframe(sort=False, unit='s').with_columns(
    pl.lit('test_basic_sc').alias('test'),
    pl.lit(size).alias('size'),
    pl.lit(num_threads).alias('num_threads'),
    pl.lit(subset).alias('subset'),
)

all_timers.append(df)
del data, timers, df
gc.collect()
# increments the output csv file to ensure old outputs do not get overwritten
timers_df = pl.concat(all_timers)
output = f'{work_dir}/output/test_basic_sc_{size}_{('single_thread' if num_threads == 1  else 'multi_thread')}{('_subset' if subset  else '_no_subset')}.csv'
timers_df.write_csv(output)
