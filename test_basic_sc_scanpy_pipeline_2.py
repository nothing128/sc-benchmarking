# Paste pre-Pyrfected Python
import gc
import os
import sys
import polars as pl  # type: ignore
from single_cell import SingleCell # type: ignore
from utils_local import TimerMemoryCollection, system_info
import scanpy as sc  # type: ignore
work_dir = 'projects/sc-benchmarking'
data_dir = '~/single-cell/SEAAD'

os.makedirs(f"{work_dir}/output", exist_ok=True)
os.makedirs(f"{work_dir}/figures", exist_ok=True)
sys.path.append(work_dir)


num_threads_options = [-1, 1]
subset_options = ["True", "False"]
size_options = ["20K", "400K", "1M"]


system_info()
# takes in cmd line args to run one instance of the toolkit
if len(sys.argv) != 4:
    print(
        f"Usage: python {sys.argv[0]} [num_threads_options] [subset_options] [size_options]"
    )
    print("num_threads_options must be -1 or 1")
    print("subset_options must be True or False")
    print("size_options must be 20K or 400K or 1M")
    exit(1)

all_timers = []
delay = 0.15
num_threads = int(sys.argv[1])
subset = True if sys.argv[2] == "True" else False
size = sys.argv[3]

if (
    ("-h" in sys.argv) or
    (num_threads not in num_threads_options) or
    (sys.argv[2] not in subset_options) or
    (size not in size_options)
):
    print(
        f"Usage: python {sys.argv[0]} [num_threads_options] [subset_options] [size_options]"
    )
    print("num_threads_options must be -1 or 1")
    print("subset_options must be True or False")
    print("size_options must be 20K or 400K or 1M")
    exit(1)


timers = TimerMemoryCollection(silent=True)

# Note: Loading is much slower from $scratch disk
# main toolkit steps
with timers('Load data'):
    data = SingleCell(f'{data_dir}/SEAAD_raw_{size}.h5ad', num_threads=1)
    data = data.set_num_threads(num_threads)

with timers('Quality control'):
    data = data.qc(
        subset=subset,
        remove_doublets=False,
        allow_float=True,
        verbose=False)

with timers('Doublet detection'):
    data = data.find_doublets(batch_column='sample')
        
# Not timed 
data = data.filter_obs(pl.col('doublet').not_())

with timers('Feature selection'):
    data = data.hvg()

with timers('Normalization'):
    data = data.normalize()

if not subset:
    print("Filtering cells...")
    data = data.filter_obs(pl.col('passed_QC'))

with timers('PCA'):
    data=data.PCA() 




with timers('Neighbor graph'):
    data = data.neighbors()  
    data = data.shared_neighbors()  

anndata = data.to_scanpy()

print("Moving shared neighbor graph to obsp['connectivities']")
anndata.obsp['connectivities'] = anndata.obsp['shared_neighbors']

# 3. Manually create the metadata dictionary that UMAP was looking for.
print("Manually creating uns['neighbors'] dictionary")
anndata.uns['neighbors'] = {
    'params': {
        'n_neighbors': 20, # Use the value you passed to your function
        'method': 'custom_sctk', # Acknowledge this wasn't Scanpy's method
        'use_rep': 'PCs',
        'metric': 'euclidean'
    },
    'connectivities_key': 'connectivities',
    'distances_key': 'distances'
}

with timers('Embedding'):
    sc.tl.umap(anndata)

# TODO: The number of clusters needs to match across libraries

with timers('Clustering (3 resolutions)'):
    for res in [0.5, 1, 2]:
        sc.tl.leiden(
            anndata, 
            flavor='igraph',
            key_added=f'leiden_res_{res:4.2f}', 
            resolution=res)

print(f'leiden_res_0.50: {len(anndata.obs['leiden_res_0.50'].unique())}')
print(f'leiden_res_1.00: {len(anndata.obs['leiden_res_1.00'].unique())}')
print(f'leiden_res_2.00: {len(anndata.obs['leiden_res_2.00'].unique())}')

#convert to toolkit
data = SingleCell(anndata)
del anndata
with timers("Plot embeddings"):
    data.plot_embedding(
        "subclass", f"{work_dir}/figures/combined_sc_scanpy_embedding_cluster_{size}_pipeline_2.png",
        embedding_key='X_umap'
    )

# Not timed
with timers("Find markers"):
    markers = data.find_markers("leiden_res_0.50", num_threads=num_threads)

# with timers('Save data'):
#    data.save(
#        f'{data_dir}/test_write.h5ad',
#        overwrite=True)

# Not timed
# os.remove(f'{data_dir}/test_write.h5ad')

print("--- Params ---")
print(f"{size=}, {num_threads=}, {subset=}")
timers.print_summary(sort=False)
# store result to polars dataframe
df = timers.to_dataframe(sort=False, unit="s").with_columns(
    pl.lit("test_basic_sc").alias("test"),
    pl.lit(size).alias("size"),
    pl.lit(num_threads).alias("num_threads"),
    pl.lit(subset).alias("subset"),
)

all_timers.append(df)
del data, timers, df
gc.collect()
# increments the output csv file to ensure old outputs do not get overwritten
timers_df = pl.concat(all_timers)
output = f'{work_dir}/output/test_basic_sc_scanpy_2_{size}_{("single_thread" if num_threads == 1  else "multi_thread")}{("_subset" if subset  else "_no_subset")}.csv'
timers_df.write_csv(output)

"""
--- System Information ---
Node: nia0028.scinet.local
CPU: 40 physical cores, 80 logical cores
Memory: 171.3 GB available / 188.6 GB total
leiden_res_0.50: 11
leiden_res_1.00: 14
leiden_res_2.00: 20
--- Params ---
size='20K', num_threads=-1, subset=True

--- Timing Summary ---
Load data took 668ms 80µs (0.2%) using 1.53 GiB (0.8%)
Quality control took 247ms 614µs (0.1%) using 1.67 GiB (0.9%)
Doublet detection took 1s 241ms (0.4%) using 1.72 GiB (0.9%)
Feature selection took 401ms 173µs (0.1%) using 1.38 GiB (0.7%)
Normalization took 158ms 666µs (0.1%) using 1.45 GiB (0.8%)
PCA took 603ms 33µs (0.2%) using 1.49 GiB (0.8%)
Neighbor graph took 4m 14s (89.6%) using 16.71 GiB (8.9%)
Embedding took 21s 951ms (7.7%) using 13.16 GiB (7.0%)
Clustering (3 resolutions) took 1s 837ms (0.6%) using 13.21 GiB (7.0%)
Plot embeddings took 2s 113ms (0.7%) using 13.16 GiB (7.0%)
Find markers took 262ms 281µs (0.1%) using 13.16 GiB (7.0%)


"""
