# Paste pre-Pyrfected Python
import gc
import os
import sys
import polars as pl  # type: ignore
from single_cell import SingleCell, cython_inline  # type: ignore
from utils_local import TimerMemoryCollection, system_info
import scanpy as sc  # type: ignore
import numpy as np

work_dir = 'projects/sc-benchmarking'
data_dir = '~/single-cell/SEAAD'

os.makedirs(f"{work_dir}/output", exist_ok=True)
os.makedirs(f"{work_dir}/figures", exist_ok=True)
sys.path.append(work_dir)


num_threads_options = [-1, 1]
subset_options = ["True", "False"]
size_options = ["20K", "400K", "1M"]

# Remove self-neighbors from each cell's list of nearest neighbors.
# These are almost always in the 0th column, but occasionally later due
# to the inaccuracy of the nearest-neighbor search. This leaves us with
# `num_neighbors + num_extra_neighbors` nearest neighbors.
remove_self_neighbors = cython_inline(r'''
    from cython.parallel cimport prange

    def remove_self_neighbors(long[:, ::1] neighbors,
                              const unsigned num_threads):
        cdef unsigned i, j, num_cells = neighbors.shape[0], \
            num_neighbors = neighbors.shape[1]
        
        if num_threads == 1:
            for i in range(num_cells):
                # If the cell is its own nearest neighbor (almost always), skip
                
                if <unsigned> neighbors[i, 0] == i:
                    continue
                
                # Find the position where the cell is listed as its own
                # self-neighbor
                
                for j in range(1, num_neighbors):
                    if <unsigned> neighbors[i, j] == i:
                        break
                
                # Shift all neighbors before it to the right, overwriting it
                
                while j > 0:
                    neighbors[i, j] = neighbors[i, j - 1]
                    j = j - 1
        else:
            for i in prange(num_cells, nogil=True,
                            num_threads=num_threads):
                if <unsigned> neighbors[i, 0] == i:
                    continue
                for j in range(1, num_neighbors):
                    if <unsigned> neighbors[i, j] == i:
                        break
                while j > 0:
                    neighbors[i, j] = neighbors[i, j - 1]
                    j = j - 1
        ''')['remove_self_neighbors']


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
    data = data.PCA() 

anndata = data.to_scanpy()

with timers('Neighbor graph'):
    sc.pp.neighbors(anndata, use_rep='PCs', n_neighbors=15)
obs=anndata.n_obs
distance_matrix_sparse = anndata.obsp['distances']
n_neighbors = anndata.uns['neighbors']['params']['n_neighbors'] # A robust way to get 15
neighbor_indices = distance_matrix_sparse.indices.reshape(obs, n_neighbors)
neighbor_indices = neighbor_indices[:, 1:]
anndata.obsm['neighbors'] = neighbor_indices.astype(np.uint32)
data=SingleCell(anndata)
with timers('Clustering (3 resolutions)'):
    data = data.cluster(
        resolution=[1, 0.5, 2],
        shared_neighbors_key='connectivities'
    )

print(f'cluster_0: {len(data.obs['cluster_0'].unique())}')
print(f'cluster_1: {len(data.obs['cluster_1'].unique())}')
print(f'cluster_2: {len(data.obs['cluster_2'].unique())}')
data.obsm['distances']=data.obsp['distances'].toarray().reshape(obs,obs)
with timers('Embedding'):
    data = data.embed(num_threads = 1, num_extra_neighbors=6)

with timers('Plot embeddings'):
    data.plot_embedding(
        'subclass', f'{work_dir}/figures/sc_embedding_cluster_{size}_pipeline_3.png')

with timers('Find markers'):
    markers = data.find_markers('cluster_0')

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
output = f'{work_dir}/output/test_basic_sc_scanpy_{size}_{("single_thread" if num_threads == 1  else "multi_thread")}{("_subset" if subset  else "_no_subset")}.csv'
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
