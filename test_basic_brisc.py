import gc
import sys
import polars as pl  
from single_cell import SingleCell  
sys.path.append('sc-benchmarking')
from utils_local import TimerMemoryCollection, system_info

num_threads = int(sys.argv[1])
subset = sys.argv[2].lower() == 'true'
size = sys.argv[3]
output = sys.argv[4]

print('--- Params ---')
print(f'{size=}, {num_threads=}, {subset=}')

system_info()
timers = TimerMemoryCollection(silent=True)

with timers('Load data'):
    data = SingleCell(
        f'single-cell/SEAAD/SEAAD_raw_{size}.h5ad', 
        num_threads=num_threads)

with timers('Quality control'):
    data = data.qc(
        subset=False,
        remove_doublets=False,
        allow_float=True,
        verbose=False)

with timers('Doublet detection'):
    data = data.find_doublets(batch_column='sample')
        
with timers('Quality control'):
    if subset:
        data = data.filter_obs(
            pl.col('doublet').not_() & pl.col('passed_QC'))

with timers('Feature selection'):
    data = data.hvg()

with timers('Normalization'):
    data = data.normalize()

with timers('PCA'):
    data = data.PCA()

with timers('Nearest neighbors'):
    data = data.neighbors().shared_neighbors()  

with timers('Clustering (3 res.)'):
    data = data.cluster(resolution=[0.5, 1.0, 2.0])

with timers('Embedding'):
    data = data.embed()

with timers('Plot embedding'):
    data.plot_embedding(    
        'subclass', 
        f'sc-benchmarking/figures/sc_embedding_subclass_{size}.png')

with timers('Find markers'):
    markers = data.find_markers('subclass')

timers.print_summary(sort=False)

df = timers.to_dataframe(sort=False, unit='s')\
    .with_columns(
        pl.lit('brisc').alias('library'),
        pl.lit('basic').alias('test'),
        pl.lit(size).alias('size'),
        pl.lit(subset).alias('subset'),
        pl.lit('single-threaded' if num_threads == 1 else 'multi-threaded')
        .alias('num_threads'))
df.write_csv(output)

del data, timers, df
gc.collect()

