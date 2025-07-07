import gc
import sys
import polars as pl  # type: ignore
import scanpy as sc  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import os

work_dir = 'projects/sc-benchmarking'
data_dir = '~/single-cell/SEAAD'

sys.path.append(work_dir)
from utils_local import TimerMemoryCollection, system_info

size = sys.argv[1]
output = sys.argv[2]

system_info()
timers = TimerMemoryCollection(silent=True)

#%% Load data
with timers('Load data'):
    data = sc.read_h5ad(f'{os.path.expanduser(data_dir)}/SEAAD_raw_{size}.h5ad')

#%% Quality control
with timers('Quality control'):
    data.var['mt'] = data.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], inplace=True, log1p=True)
    sc.pp.filter_cells(data, min_genes=100)
    sc.pp.filter_genes(data, min_cells=3)

#%% Doublet detection
with timers('Doublet detection'):
    sc.pp.scrublet(data, batch_key='sample')

#%% Quality control
with timers('Quality control'):
    data = data[data.obs['predicted_doublet'] == False].copy()

#%% Normalization
with timers('Normalization'):
    sc.pp.normalize_total(data)
    sc.pp.log1p(data)

#%% Feature selection
with timers('Feature selection'):
    sc.pp.highly_variable_genes(data, n_top_genes=2000, batch_key='sample')

#%% PCA
with timers('PCA'):
    sc.tl.pca(data)

#%% Neighbor graph
with timers('Nearest neighbors'):
    sc.pp.neighbors(data)

#%% Embedding
with timers('Embedding'):
    sc.tl.umap(data)

#%% Clustering (3 res.)
with timers('Clustering (3 res.)'):
    for res in [0.5, 1.0, 2.0]:
        sc.tl.leiden(
            data, 
            flavor='igraph',
            key_added=f'leiden_res_{res:4.2f}', 
            resolution=res)

#%% Plot embeddings
with timers('Plot embedding'):
    sc.pl.umap(data, color=['subclass'])
    plt.savefig(
        f'{work_dir}/figures/scanpy_embedding_subclass_{size}.png',
        dpi=300,
        bbox_inches='tight',
        pad_inches='layout',
    )

#%% Find markers
with timers('Find markers'):
    sc.tl.rank_genes_groups(data, groupby='subclass', method='wilcoxon')

timers.print_summary(sort=False)

timers_df = timers.to_dataframe(sort=False, unit='s').with_columns(
    pl.lit('scanpy').alias('library'),
    pl.lit('basic').alias('test'),
    pl.lit(size).alias('size'),
)
timers_df.write_csv(output)

del timers, timers_df, data
gc.collect()

