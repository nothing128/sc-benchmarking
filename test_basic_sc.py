import gc
import sys
import polars as pl  # type: ignore
from single_cell import SingleCell  # type: ignore

work_dir = 'projects/sc-benchmarking'
data_dir = '~/single-cell/SEAAD'

sys.path.append(work_dir)
from utils_local import TimerMemoryCollection, system_info

num_threads = int(sys.argv[1])
subset = sys.argv[2].lower() == 'true'
size = sys.argv[3]
output = sys.argv[4]

print('--- Params ---')
print(f'{size=}, {num_threads=}, {subset=}')

system_info()
timers = TimerMemoryCollection(silent=True)

# Note: Loading is much slower from $SCRATCH disk
# TODO: Temporarily setting `num_threads=1` for loading until shared 
# memory benchmarking is possible 

#%% Load data
with timers('Load data'):
    data = SingleCell(f'{data_dir}/SEAAD_raw_{size}.h5ad', num_threads=1)
    data = data.set_num_threads(num_threads)

#%% Quality control
with timers('Quality control'):
    data = data.qc(
        subset=subset,
        remove_doublets=False,
        allow_float=True,
        verbose=False)

#%% Doublet detection
with timers('Doublet detection'):
    data = data.find_doublets(batch_column='sample')
        
# Not timed 
data = data.filter_obs(pl.col('doublet').not_())

#%% Feature selection
with timers('Feature selection'):
    data = data.hvg()

#%% Normalization
with timers('Normalization'):
    data = data.normalize()

#%% PCA
with timers('PCA'):
    data = data.PCA()

# Not timed
if not subset:
    data = data.filter_obs(pl.col('passed_QC'))

#%% Neighbor graph
with timers('Neighbor graph'):
    data = data.neighbors()  
    data = data.shared_neighbors()  

#%% Clustering (3 resolutions)
with timers('Clustering (3 resolutions)'):
    data = data.cluster(resolution=[1, 0.5, 2])

print(f'cluster_0: {len(data.obs['cluster_0'].unique())}')
print(f'cluster_1: {len(data.obs['cluster_1'].unique())}')
print(f'cluster_2: {len(data.obs['cluster_2'].unique())}')

#%% Embedding
with timers('Embedding'):
    data = data.embed()

#%% Plot embeddings
with timers('Plot embeddings'):
    data.plot_embedding(
        'subclass', f'{work_dir}/figures/sc_embedding_subclass_{size}.png')

#%% Find markers
with timers('Find markers'):
    markers = data.find_markers('subclass')

timers.print_summary(sort=False)

df = timers.to_dataframe(sort=False, unit='s').with_columns(
    pl.lit('test_basic_sc').alias('test'),
    pl.lit(size).alias('size'),
    pl.lit(num_threads).alias('num_threads'),
    pl.lit(subset).alias('subset'),
)
df.write_csv(output)

del data, timers, df
gc.collect()

