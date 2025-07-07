import gc
import sys
import polars as pl  
from single_cell import SingleCell  

work_dir = 'projects/sc-benchmarking'
data_dir = '~/single-cell/SEAAD'

sys.path.append(work_dir)
from utils_local import TimerMemoryCollection, system_info

num_threads = int(sys.argv[1])
size = sys.argv[2]
output = sys.argv[3]

print('--- Params ---')
print(f'{size=}, {num_threads=}')

system_info()
timers = TimerMemoryCollection(silent=True)

# TODO: Temporarily setting `num_threads=1` for loading until shared 
# memory benchmarking is possible 

#%% Load data
with timers('Load data'):
    data = SingleCell(f'{data_dir}/SEAAD_raw_{size}.h5ad', num_threads=1)
    data = data.set_num_threads(num_threads)

#%% Quality control
with timers('Quality control'):
    data = data.qc(
        remove_doublets=False,
        allow_float=True,
        verbose=False)

#%% Doublet detection
with timers('Doublet detection'):
    data = data.find_doublets(batch_column='sample')

#%% Data transformation (pseudobulk / normalization)
with timers('Data transformation (pseudobulk / normalization)'):
    data = data.pseudobulk('sample', 'subclass')
    data = data.qc('ad_dx', verbose=False)

# Not timed, temporary fix for `pmi` column until `prep_data.py` is run again
data = data.with_columns_obs(pl.col('pmi').cast(pl.Float64))

if size == '20K':
    excluded_cell_types = ['L6 IT Car3', 'Sncg']
else:
    excluded_cell_types = ['L6 IT Car3', 'Sncg'] 

#%% Differential expression
with timers('Differential expression'):
    de = data\
        .library_size()\
        .DE(
        '~ ad_dx + apoe4_dosage + sex + age_at_death + '
        'log2(num_cells) + log2(library_size)',
        group='ad_dx',
        excluded_cell_types=excluded_cell_types)
    
timers.print_summary(sort=False)

df = timers.to_dataframe(sort=False, unit='s').with_columns(
    pl.lit('brisc').alias('library'),
    pl.lit('de').alias('test'),
    pl.lit(size).alias('size'),
    pl.when(num_threads == 1).then(pl.lit('single-threaded'))\
        .otherwise(pl.lit('multi-threaded')).alias('num_threads')
)
df.write_csv(output)

del data, de, timers, df
gc.collect()
