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

size = '20K'
num_threads = -1
subset = True

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

#%% Quality control (single cell)
with timers('Quality control (single cell)'):
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

#%% Pseudobulk
with timers('Pseudobulk'):
    data = data.pseudobulk('sample', 'subclass')

#%% Quality control (pseudobulk)
with timers('Quality control (pseudobulk)'):
    data = data.qc('ad_dx', verbose=False)

# Not timed, temporary fix for `pmi` column
data = data.with_columns_obs(pl.col('pmi').cast(pl.Float64))

#%% Differential expression
with timers('Differential expression'):
    de = data\
        .library_size()\
        .DE(
        '~ ad_dx + apoe4_dosage + sex + age_at_death + '
        'log2(num_cells) + log2(library_size)',
        coefficient='ad_dx',
        group='ad_dx',
        excluded_cell_types=['Lamp5 Lhx6', 'Sncg'])
    
timers.print_summary(sort=False)

df = timers.to_dataframe(sort=False, unit='s').with_columns(
    pl.lit('test_basic_sc').alias('test'),
    pl.lit(size).alias('size'),
    pl.lit(num_threads).alias('num_threads'),
    pl.lit(subset).alias('subset'),
)
df.write_csv(output)

'''
--- Timing Summary ---
Load data took 2s 409ms (6.8%) using 1.51 GiB (0.8%)
Quality control (single cell) took 815ms 39µs (2.3%) using 1.58 GiB (0.8%)
Doublet detection took 6s 833ms (19.3%) using 1.61 GiB (0.9%)
Pseudobulk took 447ms 453µs (1.3%) using 1.43 GiB (0.8%)
Quality control (pseudobulk) took 611ms 184µs (1.7%) using 0.71 GiB (0.4%)
Differential expression took 24s 251ms (68.6%) using 0.87 GiB (0.5%)

Total time: 35s 369ms
'''