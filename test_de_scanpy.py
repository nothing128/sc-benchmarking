import gc
import sys
import polars as pl  # type: ignore
import scanpy as sc  # type: ignore
import matplotlib.pyplot as plt  # type: ignore

work_dir = 'projects/sc-benchmarking'
data_dir = '~/single-cell/SEAAD'

sys.path.append(work_dir)
from utils_local import TimerMemoryCollection, system_info

size = sys.argv[1]
output = sys.argv[2]

size = "20K"

system_info()
timers = TimerMemoryCollection(silent=True)

# Note: Loading is much slower from $SCRATCH disk

#%% Load data
with timers('Load data'):
    data = sc.read_h5ad(f'{data_dir}/SEAAD_raw_{size}.h5ad')

#%% Quality control
with timers('Quality control'):
    data.var['mt'] = data.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], inplace=True, log1p=True)
    sc.pp.filter_cells(data, min_genes=100)
    sc.pp.filter_genes(data, min_cells=3)

#%% Doublet detection
with timers('Doublet detection'):
    sc.pp.scrublet(data, batch_key='sample')

# Not timed
data = data[data.obs['predicted_doublet'] == False].copy()

#%% Normalization
with timers('Normalization'):
    sc.pp.normalize_total(data)
    sc.pp.log1p(data)

# Not timed
data.obs['ad_dx'] = data.obs['ad_dx']\
    .astype(str).map({'1': 'AD', '0': 'Control'})

#%% Differential expression
with timers('Differential expression'):
    data.obs['group'] = data.obs['subclass']\
        .astype(str) + '_' + data.obs['ad_dx'].astype(str)

    de = {}
    for subclass in data.obs['subclass'].unique():
        adata_sub = data[data.obs['subclass'] == subclass].copy()
        sc.tl.rank_genes_groups(
            adata_sub, 
            groupby='group', 
            method='wilcoxon',
            key_added=f'de_{subclass}'
        )
        de[subclass] = sc.get.rank_genes_groups_df(
            adata_sub,
            group=subclass + '_AD',
            key=f'de_{subclass}'
        )

timers.print_summary(sort=False)

'''
--- Timing Summary ---
Load data took 588ms 446µs (0.7%) using 0.89 GiB (0.5%)
Quality control took 5s 461ms (6.4%) using 2.58 GiB (1.4%)
Doublet detection took 1m 5s (77.0%) using 2.89 GiB (1.5%)
Normalization took 497ms 736µs (0.6%) using 2.41 GiB (1.3%)
Differential expression took 12s 971ms (15.3%) using 2.08 GiB (1.1%)

Total time: 1m 24s
'''