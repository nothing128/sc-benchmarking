import itertools
import utils

def main():
    _original_run = utils.run
    def _patched_run(cmd, **kwargs):
        if '.sbatch' in cmd:
            cmd = cmd.replace('.sbatch', 'sbatch')
        return _original_run(cmd, **kwargs)
    utils.run = _patched_run

    sizes = ['20K', '400K', '1M']
    work_dir = 'projects/sc-benchmarking'
    output_dir = f'{work_dir}/output'
    log_dir = f'{work_dir}/logs'
    py_path = 'export PYTHONPATH=$HOME${PYTHONPATH:+:$PYTHONPATH}; '

    scripts = {
        'test_basic_brisc': f'{py_path}python {work_dir}/test_basic_brisc.py',
        'test_basic_scanpy': f'python {work_dir}/test_basic_scanpy.py',
        'test_basic_seurat': f'Rscript {work_dir}/test_basic_seurat.R',
        'test_basic_seurat_bpcells': f'Rscript {work_dir}/test_basic_seurat_bpcells.R',
        'test_transfer_brisc': f'{py_path}python {work_dir}/test_transfer_brisc.py',
        'test_transfer_scanpy': f'python {work_dir}/test_transfer_scanpy.py',
        'test_transfer_seurat': f'Rscript {work_dir}/test_transfer_seurat.R',
        'test_transfer_seurat_bpcells': f'Rscript {work_dir}/test_transfer_seurat_bpcells.R',
        'test_de_brisc': f'{py_path}python {work_dir}/test_de_brisc.py',
        'test_de_scanpy': f'python {work_dir}/test_de_scanpy.py',
        'test_de_seurat': f'Rscript {work_dir}/test_de_seurat.R',
        'test_de_seurat_bpcells': f'Rscript {work_dir}/test_de_seurat_bpcells.R',
    }

    for name, script in scripts.items():
        for size in sizes:
            if (name == 'test_basic_seurat' and size == '1M') or \
               (name.endswith('_bpcells') and size == '20K'):
                continue
            
            if 'brisc' in name and name != 'test_de_brisc':
                for threads, subset in itertools.product([-1, 1], ['True', 'False']):
                    base = f'{name}_{size}_{threads}_{subset}'
                    output = f'{output_dir}/{base}.csv'
                    log = f'{log_dir}/{base}.log'
                    cmd = f'{script} {threads} {subset} {size} {output}'
                    utils.run_slurm(cmd, log_file=log)
            elif name == 'test_de_brisc':
                for threads in [-1, 1]:
                    base = f'{name}_{size}_{threads}'
                    output = f'{output_dir}/{base}.csv'
                    log = f'{log_dir}/{base}.log'
                    cmd = f'{script} {threads} {size} {output}'
                    utils.run_slurm(cmd, log_file=log)
            else:
                base = f'{name}_{size}'
                output = f'{output_dir}/{base}.csv'
                log = f'{log_dir}/{base}.log'
                cmd = f'{script} {size} {output}'
                utils.run_slurm(cmd, log_file=log)

if __name__ == '__main__':
    main() 