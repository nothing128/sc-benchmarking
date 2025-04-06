import os
import gc
import gzip
import h5py
import polars as pl
import numpy as np
from scipy import sparse
from scipy.io import mmwrite
from single_cell import SingleCell
from utils import run 

#region functions

def write_to_mtx(adata, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    with gzip.open(f'{output_dir}/matrix.mtx.gz', 'wb', compresslevel=4) as f:
        mmwrite(f, adata.X.T)
    with gzip.open(f'{output_dir}/barcodes.tsv.gz', 'wb', 
                 compresslevel=4) as f:
        pl.DataFrame({'barcodes': list(adata.obs_names)}).write_csv(
            f, separator='\t', include_header=False)
    with gzip.open(f'{output_dir}/features.tsv.gz', 'wb',
                 compresslevel=4) as f:
        pl.DataFrame({'features': list(adata.var_names)}).write_csv(
            f, separator='\t', include_header=False)
    
def write_to_h5(adata, file):
    if not sparse.isspmatrix_csr(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    if 'feature_types' not in adata.var.columns:
        adata.var['feature_types'] = ['Gene Expression'] * adata.n_vars
    if 'genome' not in adata.var.columns:
        adata.var['genome'] = ['unknown'] * adata.n_vars
    if 'gene_ids' not in adata.var.columns:
        adata.var['gene_ids'] = adata.var.index.values.copy()
    def int_max(x): 
        return int(max(np.floor(len(str(int(max(x)))) / 4), 1) * 4)
    def str_max(x): 
        return max([len(i) for i in x])
    w = h5py.File(file, 'w')
    grp = w.create_group('matrix')
    grp.create_dataset(
        'data',
        data=np.array(adata.X.data, dtype=f'<i{int_max(adata.X.data)}'))
    grp.create_dataset(
        'indices',
        data=np.array(adata.X.indices, dtype=f'<i{int_max(adata.X.indices)}'))
    grp.create_dataset(
        'indptr',
        data=np.array(adata.X.indptr, dtype=f'<i{int_max(adata.X.indptr)}'))
    grp.create_dataset(
        'shape',
        data=np.array([adata.n_vars, adata.n_obs], 
                      dtype=f'<i{int_max([adata.n_vars, adata.n_obs])}'))
    grp.create_dataset(
        'barcodes',
        data=np.array(adata.obs_names, dtype=f'|S{str_max(adata.obs_names)}'))
    ftrs = grp.create_group('features')
    ftrs.create_dataset(
        'feature_type',
        data=np.array(adata.var.feature_types,
                      dtype=f'|S{str_max(adata.var.feature_types)}'))
    ftrs.create_dataset(
        'genome',
        data=np.array(adata.var.genome,
                      dtype=f'|S{str_max(adata.var.genome)}'))
    ftrs.create_dataset(
        'id',
        data=np.array(adata.var.gene_ids,
                      dtype=f'|S{str_max(adata.var.gene_ids)}'))
    ftrs.create_dataset(
        'name',
        data=np.array(adata.var.index,
                      dtype=f'|S{str_max(adata.var.index)}'))
    w.close()


#endregion
    
#region SEAAD

file_data = 'single-cell/SEAAD/SEAAD_MTG_RNAseq_all-nuclei.2024-02-13.h5ad'
file_metadata = 'single-cell/SEAAD/sea-ad_cohort_donor_metadata.xlsx'

if not os.path.exists(file_data):
    run(f'wget https://sea-ad-single-cell-profiling.s3.amazonaws.com/'
        'MTG/RNAseq/SEAAD_MTG_RNAseq_all-nuclei.2024-02-13.h5ad -O {file_data}')
    run(f'wget https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.'
        'divio-media.net/filer_public/b4/c7/b4c727e1-ede1-4c61-b2ee-bf1ae4a3ef68/'
        'sea-ad_cohort_donor_metadata_072524.xlsx -O {file_metadata}')
    
donor_metadata = pl.read_excel(file_metadata)
sc = SingleCell(file_data)

sc = sc\
    .cast_obs({'Donor ID': pl.String})\
    .join_obs(
        donor_metadata.select(['Donor ID'] +
        list(set(donor_metadata.columns).difference(sc.obs.columns))),
        on='Donor ID', validate='m:1')\
    .filter_obs(pl.col('Neurotypical reference').eq('False'))\
    .with_columns_obs(
        pl.when(pl.col('Consensus Clinical Dx (choice=Alzheimers disease)')
            .eq('Checked')).then(1)
            .when(pl.col('Consensus Clinical Dx (choice=Control)')
            .eq('Checked')).then(0)
            .otherwise(None)
            .alias('ad_dx'),
        pl.col('Class').cast(pl.String).replace_strict({
            'Non-neuronal and Non-neural': 'Non-neuronal',
            'Neuronal: Glutamatergic': 'Glutamatergic',
            'Neuronal: GABAergic': 'GABAergic'})
            .alias('Class'))\
    .rename_obs({
        'exp_component_name': '_index',
        'sample_id': 'cell_id', 'Donor ID': 'sample',
        'Continuous Pseudo-progression Score': 'cp_score',
        'PMI': 'pmi', 'Age at Death': 'age_at_death', 'Sex': 'sex',
        'Class': 'class', 'Subclass': 'subclass'})\
    .select_obs([
        'cell_id', 'sample', 'cp_score', 'ad_dx', 
         'pmi', 'age_at_death', 'sex', 'class', 'subclass'])\
    .filter_obs(
        pl.col('sample').is_not_null() &
        pl.col('ad_dx').is_not_null() &
        pl.col('class').is_not_null() &
        pl.col('subclass').is_not_null())\
    .drop_uns(['Great Apes Metadata', 'UW Clinical Metadata',
               'X_normalization', 'batch_condition', 'default_embedding',
               'title', 'normalized', 'QCed'])

path = 'single-cell/SEAAD/subsampled'
sizes = {1000000: '1M', 400000: '400K', 20000: '20K'}

for n, label in sizes.items():
    print(f'Subsampling to {label}')
    sc_sub = sc.subsample_obs(n=n, by_column='subclass')
    print(f'Saving {label} h5ad')
    sc_sub.save(f'{path}/SEAAD_raw_{label}.h5ad', overwrite=True)
    if label != '1M':
        print(f'Saving {label} rds')
        sc_sub.save(f'{path}/SEAAD_raw_{label}.rds', overwrite=True)
    adata = sc_sub.to_anndata()
    print(f'Saving {label} mtx')
    write_to_mtx(adata, f'{path}/SEAAD_raw_{label}')
    print(f'Saving {label} h5')
    write_to_h5(adata, f'{path}/SEAAD_raw_{label}.h5')
    del sc_sub; gc.collect()
    del adata; gc.collect()

#endregion

# sc = SingleCell(
#     'projects/rrg-wainberg/single-cell/Green/'
#     'p400_qced_shareable.h5ad')
# sc.X.sort_indices()
# sc = sc\
#     .rename_obs({
#         '_index': 'cell_id', 'projid': 'sample',
#         'state': 'subclass'})\
#     .with_columns_obs(
#         pl.col('subset').cast(pl.String).replace({
#             'CUX2+': 'Excitatory'})
#             .alias('class'))\
#     .select_obs(
#         ['sample', 'class', 'subclass'])\
#     .rename_var({'_index': 'gene_symbol'})\
#     .filter_obs(
#         pl.col('sample').is_not_null() & 
#         pl.col('class').is_not_null() &
#         pl.col('subclass').is_not_null())

# path = 'projects/rrg-wainberg/single-cell/Green/Subsampled'
# sizes = {1000000: '1M', 500000: '500K', 100000: '100K',  20000: '20K'}
# for n, label in sizes.items():
#     sc_sub = sc.subsample_obs(n=n)
#     sc_sub.save(f'{path}/Green_raw_{label}.h5ad', overwrite=True)
#     if label != '1M':
#         sc_sub.save(f'{path}/Green_raw_{label}.rds', overwrite=True)
#     del sc_sub; gc.collect()

