import sys, gc
import polars as pl
sys.path.append('projects/def-wainberg/karbabi/utils')
from single_cell import SingleCell
from utils import print_df

sc = SingleCell(
    'projects/def-wainberg/single-cell/Green/'
    'p400_qced_shareable.h5ad')
sc.X.sort_indices()
sc = sc\
    .rename_obs({
        '_index': 'cell_id', 'projid': 'sample',
        'state': 'subclass'})\
    .with_columns_obs(
        pl.col('subset').cast(pl.String).replace({
            'CUX2+': 'Excitatory'})
            .alias('class'))\
    .select_obs(
        ['sample', 'class', 'subclass'])\
    .rename_var({'_index': 'gene_symbol'})\
    .filter_obs(
        pl.col('sample').is_not_null() & 
        pl.col('class').is_not_null() &
        pl.col('subclass').is_not_null())

path = 'projects/def-wainberg/single-cell/Subsampled'
sizes = {1000000: '1M', 100000: '100K', 10000: '10K'}
for n, label in sizes.items():
    sc_sub = sc.subsample_obs(n=n)
    sc_sub.save(f'{path}/Green_raw_{label}.h5ad', overwrite=True)
    sc_sub.save(f'{path}/Green_raw_{label}.rds', overwrite=True)
    del sc_sub; gc()


sc = SingleCell(
    'projects/def-wainberg/single-cell/Green/'
    'p400_qced_shareable.h5ad')\
    .filter_obs(pl.col.subset.is_not_null())
print_df(sc.obs['subset'].value_counts())
print(sc.shape[0])
sc_sub1 = sc.subsample_obs(n=10000)
print(sc_sub1.shape[0])
sc_sub2 = sc.subsample_obs(n=10000, by_column='subset')
print(sc_sub2.shape[0])

sc_sub1.to_seurat('tmp')