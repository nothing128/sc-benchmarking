import gc
import sys
import polars as pl  
from single_cell import SingleCell 

work_dir = 'projects/sc-benchmarking'
data_dir = '~/single-cell/SEAAD'

sys.path.append(work_dir)
from utils_local import TimerMemoryCollection, system_info

num_threads = int(sys.argv[1])
subset = sys.argv[2].lower() == "true"
size = sys.argv[3]
output = sys.argv[4]

size_ref = {'1M': '600K', '400K': '200K', '20K': '10K'}

print('--- Params ---')
print(f'{size=}, {num_threads=}, {subset=}')

system_info()
timers = TimerMemoryCollection(silent=False)

#%% Load data (query)
with timers('Load data (query)'):
    data_query = SingleCell(
        f'{data_dir}/SEAAD_raw_{size}.h5ad', num_threads=1)
    data_query = data_query.set_num_threads(num_threads)

#%% Load data (ref)
with timers('Load data (ref)'):
    data_ref = SingleCell(
        f"{data_dir}/SEAAD_ref_{size_ref[size]}.h5ad", num_threads=1)
    data_ref = data_ref.set_num_threads(num_threads)
    
#%% Quality control
with timers('Quality control'):
    data_query = data_query.qc(
        subset=False,
        remove_doublets=False,
        allow_float=True,
        verbose=False)

#%% Doublet detection
with timers('Doublet detection'):
    data_query = data_query.find_doublets(batch_column='sample')
        
#%% Quality control
with timers('Quality control'):
    if subset:
        data_query = data_query.filter_obs(
            pl.col('doublet').not_() & pl.col('passed_QC'))
    
#%% Feature selection
with timers('Feature selection'):
    data_ref, data_query = data_ref.hvg(data_query)

#%% Normalization
with timers('Normalization'):
    data_ref = data_ref.normalize()
    data_query = data_query.normalize()

#%% PCA
with timers('PCA'):
    data_ref, data_query = data_ref.PCA(data_query)

#%% Transfer labels
with timers('Transfer labels'):
    data_ref, data_query = data_ref.harmonize(
        data_query,
        original=True,
        num_threads=1,
        early_stopping=True,
        max_clustering_iterations=20)
    data_query = data_query.label_transfer_from(
        data_ref, 'subclass')

with pl.Config(tbl_rows=-1):
    print('--- Transfer Accuracy ---')
    accuracy_out = data_query.obs\
            .with_columns(
                pl.col(["subclass", "cell_type"]).cast(pl.String))\
            .group_by('subclass')\
            .agg(
                n_correct=(pl.col("cell_type") == pl.col("subclass")).sum(),
                n_total=pl.len())\
            .pipe(lambda df: pl.concat([
                df,
                df.sum().with_columns(subclass=pl.lit("Total"))]))\
            .with_columns(
                percent_correct=pl.col("n_correct") / pl.col("n_total") * 100)\
            .sort(
                pl.when(pl.col("subclass") == "Total").then(1).otherwise(0),
                pl.col("subclass"))
    print(accuracy_out)
    thread_mode = "single_thread" if num_threads == 1 else "multi_thread"
    subset_mode = "subset" if subset else "no_subset"
    output_file = f'{work_dir}/output/test_transfer_sc_{size}_{thread_mode}_{subset_mode}_accuracy.csv'
    accuracy_out.write_csv(output_file)

timers.print_summary(sort=False)

df = timers.to_dataframe(sort=False, unit='s').with_columns(
    pl.lit('brisc').alias('library'),
    pl.lit('transfer').alias('test'),
    pl.lit(size).alias('size'),
    pl.lit(num_threads).alias('num_threads'),
    pl.lit(subset).alias('subset'),
)
df.write_csv(output)

del data_query, data_ref, timers, df
gc.collect()

