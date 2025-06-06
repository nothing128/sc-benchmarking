from utils import run_slurm  # type: ignore
import itertools

size_options = ['20K', '400K', '1M']
pipelines = {
    'sc-python': 'python test_basic_sc.py {threads} {subset} {size}',
    'scanpy': 'python test_basic_scanpy.py {size}',
    'seurat': 'Rscript test_basic_seurat.R {size}',
    'seurat-bpcells': 'Rscript test_basic_seurat_bpcells.R {size}',
}

if __name__ == '__main__':
    for size in size_options:
        for name, command in pipelines.items():
            if name == 'sc-python':
                thread_opts = [-1, 1]
                subset_opts = ['True', 'False']
                params = itertools.product(thread_opts, subset_opts)
                for threads, subset in params:
                    cmd = command.format(
                        threads=threads, subset=subset, size=size
                    )
                    run_slurm(cmd)
            else:
                cmd = command.format(size=size)
                run_slurm(cmd)