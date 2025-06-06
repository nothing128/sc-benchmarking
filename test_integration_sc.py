import os
import gc
import sys
import polars as pl
from single_cell import SingleCell

work_dir = 'projects/sc-benchmarking'
data_dir = 'single-cell/SEAAD/subsampled'

sys.path.append(work_dir)
from utils_local import TimerMemoryCollection, system_info

size_options = ["20K", "400K", "1M"]
assert len(sys.argv) == 2

size = sys.argv[1]
assert size in size_options

system_info()
timers = TimerMemoryCollection(silent=False)

with timers("Load data"):
    data = SingleCell(
        f"{data_dir}/SEAAD_raw_{size}.h5ad", num_threads=1)

with timers("Load reference data"):
    data_ref = SingleCell(
        f"{data_dir}/Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad", 
        num_threads=1)

