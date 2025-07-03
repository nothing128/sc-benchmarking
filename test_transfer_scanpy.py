import gc
import sys
import warnings
import polars as pl  
import scanpy as sc 
warnings.filterwarnings('ignore')

work_dir = 'projects/sc-benchmarking'
data_dir = 'single-cell/SEAAD'

sys.path.append(work_dir)
from utils_local import TimerMemoryCollection, system_info

size = sys.argv[1]
output = sys.argv[2]

size_ref = {'1M': '600K', '400K': '200K', '20K': '10K'}

system_info()
timers = TimerMemoryCollection(silent=True)

print('--- Params ---')
print(f'{size=}')

#%% Load data (query)
with timers('Load data (query)'):
    data_query = sc.read_h5ad(f'{data_dir}/SEAAD_raw_{size}.h5ad')

# Not timed 
data_query.obs['subclass_orig'] = data_query.obs['subclass']
data_query.obs = data_query.obs.drop(columns=['subclass'])

#%% Load data (ref)
with timers('Load data (ref)'):
    data_ref = sc.read_h5ad(f'{data_dir}/SEAAD_ref_{size_ref[size]}.h5ad')

#%% Quality control
with timers('Quality control'):
    data_query.var['mt'] = data_query.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(
        data_query, qc_vars=['mt'], log1p=True, inplace=True)
    sc.pp.filter_cells(data_query, min_genes=100, copy=False)
    sc.pp.filter_genes(data_query, min_cells=3, copy=False)
    
#%% Doublet detection
with timers('Doublet detection'):
    sc.pp.scrublet(data_query, batch_key='sample', copy=False)

#%% Quality control
with timers('Quality control'):
    data_query = data_query[
        data_query.obs['predicted_doublet'] == False].copy()

#%% Normalization
with timers('Normalization'):
    sc.pp.normalize_total(data_ref)
    sc.pp.log1p(data_ref)    
    sc.pp.normalize_total(data_query)
    sc.pp.log1p(data_query)

# Note: Highly variable gene selection is not explicitly done in the 
# scanpy vignette, but the exemplar data are pre-filtered.

#%% Feature selection
with timers('Feature selection'):
    var_names = data_ref.var_names.intersection(data_query.var_names)
    data_ref = data_ref[:, var_names].copy()
    data_query = data_query[:, var_names].copy()

    sc.pp.highly_variable_genes(data_ref, n_top_genes=2000, batch_key='sample')
    hvg_genes = data_ref.var[data_ref.var['highly_variable']].index
    data_ref = data_ref[:, hvg_genes].copy()
    data_query = data_query[:, hvg_genes].copy()

# Note: The exemplar data used in the scanpy vignette are scaled
# Not timed 
sc.pp.scale(data_ref)
sc.pp.scale(data_query)

#%% PCA
with timers('PCA'):
    sc.pp.pca(data_ref)

#%% Transfer labels
with timers('Transfer labels'):
    sc.pp.neighbors(data_ref)
    sc.tl.ingest(data_query, data_ref, obs='subclass', embedding_method='pca')

with pl.Config(tbl_rows=-1):
    print('--- Transfer Accuracy ---')
    df = pl.from_pandas(
        data_query.obs[['subclass_orig', 'subclass']].reset_index(drop=True))\
        .with_columns(
            pl.col(["subclass_orig", "subclass"]).cast(pl.String))\
        .group_by('subclass')\
        .agg(
            n_correct=(pl.col("subclass_orig") == pl.col("subclass")).sum(),
            n_total=pl.len())\
        .pipe(lambda df: pl.concat([
            df,
            df.sum().with_columns(subclass=pl.lit("Total"))])\
        .with_columns(
            percent_correct=pl.col("n_correct") / pl.col("n_total") * 100)\
        .sort(  
            pl.when(pl.col("subclass") == "Total").then(1).otherwise(0),
            pl.col("subclass")))
    print(df)
df.write_csv(f'{work_dir}/output/test_transfer_scanpy_{size}_accuracy.csv')

timers.print_summary(sort=False)

df = timers.to_dataframe(sort=False, unit='s').with_columns(
    pl.lit('test_transfer_scanpy').alias('test'),
    pl.lit(size).alias('size'),
)
df.write_csv(output)

del data_query, data_ref, timers, df
gc.collect()
