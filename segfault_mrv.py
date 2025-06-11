# Paste pre-Pyrfected Python
import gc
import os
import sys
import polars as pl  # type: ignore
from single_cell import SingleCell  # type: ignore
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


system_info()
num_threads = -1
subset = True
size = "20K"

# Note: Loading is much slower from $scratch disk
# main toolkit steps
data = SingleCell(f'{data_dir}/SEAAD_raw_{size}.h5ad', num_threads=1)
data = data.set_num_threads(num_threads)

data = data.qc(
        subset=subset,
        remove_doublets=False,
        allow_float=True,
        verbose=False)

data = data.find_doublets(batch_column='sample')

# Not timed
data = data.filter_obs(pl.col('doublet').not_())

data = data.hvg()

data = data.normalize()

if not subset:
    print("Filtering cells...")
    data = data.filter_obs(pl.col('passed_QC'))

data_for_pca = data.copy()

data_with_pcs = data_for_pca.PCA()
pca_result_matrix = data_with_pcs.obsm['PCs']

if pca_result_matrix.dtype != np.float32:
    pca_result_matrix = pca_result_matrix.astype(np.float32)

# Store under both names to satisfy both toolkits
data.obsm['PCs'] = pca_result_matrix
data.obsm['X_pca'] = pca_result_matrix
del data_for_pca, data_with_pcs

anndata = data.to_scanpy()
del data

sc.pp.neighbors(anndata, n_neighbors=15)

dist_matrix = anndata.obsp['distances']
n_obs = anndata.n_obs
neighbor_array = np.zeros((n_obs, 15), dtype=np.uint32)
distance_array = np.zeros((n_obs, 15), dtype=np.float32)

for i in range(n_obs):
    raw_distances = dist_matrix[i].data[:15]
    raw_neighbors = dist_matrix[i].indices[:15]

    # Sort by distance and apply the same sorting to neighbors
    sorting_indices = np.argsort(raw_distances)
    sorted_distances = raw_distances[sorting_indices]
    sorted_neighbors = raw_neighbors[sorting_indices]

    distance_array[i, :] = sorted_distances
    neighbor_array[i, :] = sorted_neighbors

anndata.obsm['neighbors'] = neighbor_array
anndata.obsm['distances'] = distance_array
connectivities_graph = anndata.obsp['connectivities']
del anndata.obsp  # Clean the object for conversion

# --- 5. Convert and Restore ---
data = SingleCell(anndata)
data.obsp['connectivities'] = connectivities_graph

data = data.cluster(
    resolution=[1, 0.5, 2],
    shared_neighbors_key='connectivities'
)

data.obsm['distances']=data.obsp['distances']
data = data.embed(PC_key='X_pca')
print("embed done")
data.plot_embedding(
        'subclass', f'{work_dir}/figures/sc_embedding_cluster_{size}_pipeline_3.png')


markers = data.find_markers('cluster_0')
