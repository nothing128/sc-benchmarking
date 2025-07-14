import gc
import sys
import polars as pl  
import scanpy as sc  
import matplotlib.pyplot as plt  
sys.path.append('sc-benchmarking')
from utils_local import TimerMemoryCollection, system_info

size = sys.argv[1]
output = sys.argv[2]

system_info()
timers = TimerMemoryCollection(silent=True)

with timers('Load data'):
    data = sc.read_h5ad(f'single-cell/SEAAD/SEAAD_raw_{size}.h5ad')

with timers('Quality control'):
    data.var['mt'] = data.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], inplace=True, log1p=True)
    sc.pp.filter_cells(data, min_genes=100)
    sc.pp.filter_genes(data, min_cells=3)

with timers('Doublet detection'):
    sc.pp.scrublet(data, batch_key='sample')

with timers('Quality control'):
    data = data[data.obs['predicted_doublet'] == False].copy()

with timers('Normalization'):
    sc.pp.normalize_total(data)
    sc.pp.log1p(data)

with timers('Feature selection'):
    sc.pp.highly_variable_genes(data, n_top_genes=2000, batch_key='sample')

with timers('PCA'):
    sc.tl.pca(data)

with timers('Nearest neighbors'):
    sc.pp.neighbors(data)

with timers('Embedding'):
    sc.tl.umap(data)

with timers('Clustering (3 res.)'):
    for res in [0.5, 1.0, 2.0]:
        sc.tl.leiden(
            data, 
            flavor='igraph',
            key_added=f'leiden_res_{res:4.2f}', 
            resolution=res)

with timers('Plot embedding'):
    sc.pl.umap(data, color=['subclass'])
    plt.savefig(
        f'sc-benchmarking/figures/scanpy_embedding_subclass_{size}.png',
        dpi=300,
        bbox_inches='tight',
        pad_inches='layout',
    )

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

