from utils import run_slurm
import itertools

num_threads_options = [-1, 1]
subset_options = ["True", "False"]
size_options = ['20K', '400K', '1M']
params = itertools.product(
    size_options, num_threads_options, subset_options
)

for size, num_threads, subset in params:
    run_slurm(f"python test_sc.py {num_threads} {subset} {size}")
