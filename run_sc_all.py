import itertools
import utils

# Monkey-patch utils.run to use sbatch instead of .sbatch as a workaround for
# cluster submission issues.
_original_run = utils.run
def _patched_run(cmd, **kwargs):
    if '.sbatch' in cmd:
        cmd = cmd.replace('.sbatch', 'sbatch')
    return _original_run(cmd, **kwargs)
utils.run = _patched_run

size_options = ['20K', '400K', '1M']
work_dir = 'projects/sc-benchmarking'
output_dir = f'{work_dir}/output'
log_dir = f'{work_dir}/logs'

python_path_export = 'export PYTHONPATH=$HOME${{PYTHONPATH:+:$PYTHONPATH}}; '
basic_pipelines = {
    'test_basic_sc': f'{python_path_export}python {work_dir}/test_basic_sc.py {{threads}} {{subset}} {{size}} {{output}}',
    'test_basic_scanpy': f'python {work_dir}/test_basic_scanpy.py {{size}} {{output}}',
    'test_basic_seurat': f'Rscript {work_dir}/test_basic_seurat.R {{size}} {{output}}',
    'test_basic_seurat_bpcells': f'Rscript {work_dir}/test_basic_seurat_bpcells.R {{size}} {{output}}',
}
transfer_pipelines = {
    'test_transfer_sc': f'{python_path_export}python {work_dir}/test_transfer_sc.py {{threads}} {{subset}} {{size}} {{output}}',
    'test_transfer_scanpy': f'python {work_dir}/test_transfer_scanpy.py {{size}} {{output}}',
    'test_transfer_seurat': f'Rscript {work_dir}/test_transfer_seurat.R {{size}} {{output}}',
    'test_transfer_seurat_bpcells': f'Rscript {work_dir}/test_transfer_seurat_bpcells.R {{size}} {{output}}',
}

def run_job(name, command, size, threads=None, subset=None):
    if threads is not None:
        base_name = f'{name}_{size}_{threads}_{subset}'
        cmd_params = {'threads': threads, 'subset': subset}
    else:
        base_name = f'{name}_{size}'
        cmd_params = {}

    output_fname = f'{output_dir}/{base_name}.csv'
    log_fname = f'{log_dir}/{base_name}.log'
    
    cmd = command.format(size=size, output=output_fname, **cmd_params)
    utils.run_slurm(cmd, log_file=log_fname)

if __name__ == '__main__':
    for size in size_options:
        for name, command in basic_pipelines.items():
            if name == 'test_basic_sc':
                for threads, subset in itertools.product([-1, 1], ['True', 'False']):
                    run_job(name, command, size, threads, subset)
            else:
                if name == 'test_basic_seurat' and size == '1M':
                    continue
                if name == 'test_basic_seurat_bpcells' and size == '20K':
                    continue
                run_job(name, command, size)

    for size in size_options:
        for name, command in transfer_pipelines.items():
            if name == 'test_transfer_sc':
                for threads, subset in itertools.product([-1, 1], ['True', 'False']):
                    run_job(name, command, size, threads, subset)
            else:
                if name == 'test_transfer_seurat' and size == '1M':
                    continue
                if name == 'test_transfer_seurat_bpcells' and size == '20K':
                    continue
                run_job(name, command, size)