import polars as pl
import gc
from single_cell import SingleCell

sc = SingleCell(
    'projects/rrg-wainberg/single-cell/Green/'
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

path = 'projects/rrg-wainberg/single-cell/Green/Subsampled'
sizes = { 1000000: '1M', 500000: '500K', 100000: '100K',  20000: '20K'}
for n, label in sizes.items():
    sc_sub = sc.subsample_obs(n=n)
    sc_sub.save(f'{path}/Green_raw_{label}.h5ad', overwrite=True)
    if label != '1M':
        sc_sub.save(f'{path}/Green_raw_{label}.rds', overwrite=True)
    del sc_sub; gc.collect()
