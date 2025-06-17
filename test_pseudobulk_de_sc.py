import gc
import sys
import polars as pl  # type: ignore
from single_cell import SingleCell  # type: ignore

work_dir = 'projects/sc-benchmarking'
data_dir = 'single-cell/SEAAD'

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

with timers('Load data'):
    data = SingleCell(f'{data_dir}/SEAAD_raw_{size}.h5ad', num_threads=1)
    data = data.set_num_threads(num_threads)

with timers('Quality control (single cell)'):
    data = data.qc(
        subset=subset,
        remove_doublets=False,
        allow_float=True,
        verbose=False)

with timers('Doublet detection'):
    data = data.find_doublets(batch_column='sample')
        
# Not timed 
data = data.filter_obs(pl.col('doublet').not_())

with timers('Pseudobulk'):
    data = data.pseudobulk('sample', 'subclass')

with timers('Quality control (pseudobulk)'):
    data = data.qc('ad_dx', verbose=False)

with timers('Differential expression'):
    de = data\
        .library_size()\
        .DE(
        '~ ad_dx + apoe4_dosage + sex + age_at_death + '
        'log2(num_cells) + log2(library_size)',
        coefficient='ad_dx',
        group='ad_dx')