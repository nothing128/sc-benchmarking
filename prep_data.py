import os, gc
import polars as pl
from scipy.io import mmwrite
from single_cell import SingleCell
from utils import run 

def save_to_10x(data, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    mmwrite(os.path.join(output_dir, "matrix.mtx"), data.X)
    df_obs = pl.DataFrame({"barcodes": list(data.obs_names)})
    df_obs.write_csv(os.path.join(output_dir, "barcodes.tsv"),
                     separator="\t", has_header=False)
    df_var = pl.DataFrame({"features": list(data.var_names)})
    df_var.write_csv(os.path.join(output_dir, "features.tsv"),
                     separator="\t", has_header=False)


file = 'single-cell/SEAAD/SEAAD_MTG_RNAseq_all-nuclei.2024-02-13.h5ad'
if not os.path.exists(file):
    run(f'wget https://sea-ad-single-cell-profiling.s3.amazonaws.com/'
        'MTG/RNAseq/SEAAD_MTG_RNAseq_all-nuclei.2024-02-13.h5ad -O {file}')

sc = SingleCell(file)
sc = sc\
    .filter_obs(pl.col('Neurotypical reference').eq('False'))\
    .rename_obs({
        'sample_id': 'cell_id', 'Donor ID': 'sample',
        'Continuous Pseudo-progression Score': 'cp_score',
        'Class': 'class', 'Subclass': 'subclass'})\
    .with_columns_obs(
        pl.when(pl.col('Consensus Clinical Dx (choice=Alzheimers disease)')
            .eq('Checked')).then(1)
            .when(pl.col('Consensus Clinical Dx (choice=Control)')
            .eq('Checked')).then(0)
            .otherwise(None)
            .alias('ad_dx'))\
    .select_obs(
        ['cell_id', 'sample', 'cp_score', 'ad_dx', 'class', 'subclass'])\
    .filter_obs(
        pl.col('sample').is_not_null() &
        pl.col('ad_dx').is_not_null() &
        pl.col('class').is_not_null() &
        pl.col('subclass').is_not_null())

save_to_10x(sc, 'scratch/SEAAD')

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

