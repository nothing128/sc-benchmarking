# Paste pre-Pyrfected Python
import gc
import os
import sys
import polars as pl  # type: ignore
from single_cell import SingleCell  # type: ignore
from utils_local import TimerMemoryCollection, system_info
import scanpy as sc  # type: ignore

work_dir = 'projects/sc-benchmarking'
data_dir = 'single-cell/SEAAD'

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
with timers("Load data (h5ad/rds)"):
    data = SingleCell(f"{data_dir}/SEAAD_raw_{size}.h5ad", num_threads=1)
    data = data.set_num_threads(num_threads)
print(f"X num_threads: {data.X._num_threads}")

# Note: QC filters are matched across libraries for timing, then
# standardized by filtering to single_cell.py qc cells
with timers("Quality control"):
    data.qc(
        subset=subset,
        max_mito_fraction=0.05,
        min_genes=100,
        nonzero_MALAT1=False,
        remove_doublets=False,
        allow_float=True,
        verbose=False,
        num_threads=num_threads,
    )
# Not timed
if subset:
    data = data.filter_obs(pl.col("tmp_passed_QC")).with_uns(QCed=True)
else:
    data = data.rename_obs({"tmp_passed_QC": "passed_QC"}).with_uns(QCed=True)
with timers("Doublet detection"):
    data = data.find_doublets(batch_column="sample", num_threads=num_threads)
print(f"cells: {data.shape[0]}, genes: {data.shape[1]}")

with timers("Feature selection"):
    data = data.hvg(num_threads=num_threads)

with timers("Normalization"):
    data = data.normalize(num_threads=num_threads)

with timers("PCA"):
    data = data.PCA(num_threads=num_threads)



# Not timed
if not subset:
    data = data.filter_obs(pl.col("passed_QC"))
    
# convert to scanpy/anndata
data = data.to_scanpy()
with timers("Neighbor graph"):
    sc.pp.neighbors(data)

with timers("Embedding"):
    sc.tl.umap(data)

# TODO: The number of clusters needs to match across libraries

with timers("Clustering (3 resolutions)"):
    for res in [0.5, 1, 2]:
        sc.tl.leiden(
            data, flavor="igraph", key_added=f"leiden_res_{res:4.2f}", resolution=res
        )

print(f"leiden_res_0.50: {len(data.obs['leiden_res_0.50'].unique())}")
print(f"leiden_res_1.00: {len(data.obs['leiden_res_1.00'].unique())}")
print(f"leiden_res_2.00: {len(data.obs['leiden_res_2.00'].unique())}")

#convert to toolkit
data = SingleCell(data)

with timers("Plot embeddings"):
    data.plot_embedding(
        "cluster_0", f"{work_dir}/figures/sc_embedding_cluster_{size}.png"
    )

# Not timed
with timers("Find markers"):
    markers = data.find_markers("cluster_0", num_threads=num_threads)

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
output = f'{work_dir}/output/test_combined_sc_scanpy_{size}_{("single_thread" if num_threads == 1  else "multi_thread")}{("_subset" if subset  else "_no_subset")}.csv'
timers_df.write_csv(output)

"""
--- System Information ---
Node: nl10603.narval.calcul.quebec
CPU: 64 physical cores, 64 logical cores
Memory: 1987.9 GB available / 2015.4 GB total

--- Params ---
size='20K', num_threads=-1, subset=True, drop_X=True

--- Timing Summary ---
Load data took 332ms 514µs (0.8%)
Quality control took 57ms 488µs (0.1%)
Doublet detection took 3s 837ms (8.7%)
Feature selection took 188ms 536µs (0.4%)
Normalization took 11ms 741µs (0.0%)
PCA took 33s 491ms (75.9%)
KNN took 929ms 249µs (2.1%)
SNN took 85ms 344µs (0.2%)
Neighbor graph took 1s 14ms (2.3%)
Clustering (3 resolutions) took 121ms 646µs (0.3%)
Embedding took 1s 112ms (2.5%)
Plot embeddings took 2s 36ms (4.6%)
Find markers took 75ms 189µs (0.2%)
Save data took 810ms 135µs (1.8%)

Total time: 44s 104ms


"""
