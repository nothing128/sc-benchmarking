import gc
import sys
import polars as pl  # type: ignore
import scanpy as sc  # type: ignore
import matplotlib.pyplot as plt  # type: ignore

work_dir = 'projects/sc-benchmarking'
data_dir = 'single-cell/SEAAD'

sys.path.append(work_dir)
from utils_local import TimerMemoryCollection, system_info

size = sys.argv[1]
output = sys.argv[2]

system_info()
timers = TimerMemoryCollection(silent=True)

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

timers_df = timers.to_dataframe(sort=False, unit='s').with_columns(
    pl.lit('scanpy').alias('library'),
    pl.lit('de').alias('test'),
    pl.lit(size).alias('size'),
)
timers_df.write_csv(output)

del data, de, adata_sub, timers, timers_df
gc.collect()

