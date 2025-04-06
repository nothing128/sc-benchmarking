import os
import gc
import sys
import polars as pl
from single_cell import SingleCell

work_dir = 'projects/def-wainberg/karbabi/sc-benchmarking'
sys.path.append(work_dir)

from utils_local import TimerCollection, system_info

system_info()





