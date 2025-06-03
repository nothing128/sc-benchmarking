from utils import run_slurm
import itertools
import sys

num_threads_options = [-1, 1]
subset_options = ["True", "False"]
size_options = ['20K', '400K', '1M']
params = itertools.product(
    size_options, num_threads_options, subset_options
)

def print_usage():
    print(f"Usage: python {sys.argv[0]} [1-2] ...")
    print("If argv[1]==1, queues an instance for every option combination for sc")
    print("If argv[1]==2, queues argv[2] instances for every pipeline using the 400K dataset (except base seurat)")

if __name__ == "__main__":
    if sys.argv[1]=="1":
        for size, num_threads, subset in params:
            run_slurm(f"python test_sc.py {num_threads} {subset} {size}")
    elif sys.argv[1]=="2" and len(sys.argv)==3:
        sizes = ['400K']
        params = itertools.product(
            sizes, num_threads_options, subset_options
        )
        for i in range(int(sys.argv[2])):
            for size, num_threads, subset in params:
                run_slurm(f"python test_sc.py {num_threads} {subset} {size}")
            run_slurm(f"Rscript test_basic_seurat_bpcells.R {sizes[0]}")
            run_slurm(f"python test_basic_scanpy.R {sizes[0]}")
    else:
        print_usage()
        exit(1)
