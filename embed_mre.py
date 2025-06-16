from single_cell import SingleCell, cython_inline  # type: ignore
import scanpy as sc  # type: ignore
import numpy as np

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


data = SingleCell('~/single-cell/SEAAD/SEAAD_raw_20K.h5ad')\
    .qc(allow_float=True, verbose=False)\
    .hvg()\
    .normalize()\
    .PCA()\
    .to_scanpy()

sc.pp.neighbors(data, use_rep='PCs', n_neighbors=15)
distance_matrix_sparse = data.obsp['distances']
obs= data.n_obs
neighbor_indices = distance_matrix_sparse.indices.reshape(obs, 16)\
    .astype(np.int64)
remove_self_neighbors(neighbor_indices, num_threads=1)
neighbor_indices = neighbor_indices[:, 1:]
data.obsm['neighbors'] = neighbor_indices
data.obsm['neighbors'] = neighbor_indices.astype(np.uint32)
data.obsm['distances']=data.obsp['distances'].toarray().reshape(obs,obs)
data=SingleCell(data)
print(np.max(data.obsm['neighbors']))
data = data.embed()