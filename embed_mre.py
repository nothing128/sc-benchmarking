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





# test_de_brisc mre
import polars as pl  
from single_cell import SingleCell  
data = SingleCell(f'~/single-cell/SEAAD/SEAAD_raw_20K.h5ad', num_threads=32)
data = data.qc(allow_float=True,verbose=False)
data = data.pseudobulk('sample', 'subclass')
data = data.qc('ad_dx', verbose=False)
data = data.with_columns_obs(pl.col('pmi').cast(pl.Float64))
excluded_cell_types = ['L6 IT Car3', 'Sncg']
de = data\
        .library_size()\
        .DE(
        '~ ad_dx + apoe4_dosage + sex + age_at_death + '
        'log2(num_cells) + log2(library_size)',
        group='ad_dx',
        excluded_cell_types=excluded_cell_types)
        
        
# asd

import gc
import sys
import polars as pl  
from single_cell import SingleCell  

work_dir = 'projects/sc-benchmarking'
data_dir = '~/single-cell/SEAAD'

sys.path.append(work_dir)
from utils_local import TimerMemoryCollection, system_info

num_threads=48
size = "20K"
output = "delete_this.csv"

print('--- Params ---')
print(f'{size=}, {num_threads=}')

system_info()
timers = TimerMemoryCollection(silent=True)

# TODO: Temporarily setting `num_threads=1` for loading until shared 
# memory benchmarking is possible 

#%% Load data
with timers('Load data'):
    data = SingleCell(f'{data_dir}/SEAAD_raw_{size}.h5ad', num_threads=num_threads)

#%% Quality control
with timers('Quality control'):
    data = data.qc(
        remove_doublets=False,
        allow_float=True,
        verbose=False)

#%% Doublet detection
with timers('Doublet detection'):
    data = data.find_doublets(batch_column='sample')

#%% Data transformation (pseudobulk / normalization)
with timers('Data transformation (pseudobulk / normalization)'):
    data = data.pseudobulk('sample', 'subclass')
    data = data.qc('ad_dx', verbose=False)

# Not timed, temporary fix for `pmi` column until `prep_data.py` is run again
data = data.with_columns_obs(pl.col('pmi').cast(pl.Float64))

if size == '20K':
    excluded_cell_types = ['L6 IT Car3', 'Sncg']
else:
    excluded_cell_types = ['L6 IT Car3', 'Sncg'] 

#%% Differential expression
with timers('Differential expression'):
    de = data\
        .library_size()\
        .DE(
        '~ ad_dx + apoe4_dosage + sex + age_at_death + '
        'log2(num_cells) + log2(library_size)',
        group='ad_dx',
        excluded_cell_types=excluded_cell_types)
    
timers.print_summary(sort=False)

df = timers.to_dataframe(sort=False, unit='s').with_columns(
    pl.lit('brisc').alias('library'),
    pl.lit('de').alias('test'),
    pl.lit(size).alias('size'),
    pl.when(num_threads == 1).then(pl.lit('single-threaded'))\
        .otherwise(pl.lit('multi-threaded')).alias('num_threads')
)
df.write_csv(output)

del data, de, timers, df
gc.collect()




from single_cell import SingleCell 
data = SingleCell(
        f'single-cell/SEAAD/SEAAD_raw_1.2M.h5ad', 
        num_threads=96)
data = data.qc(
        subset=False,
        remove_doublets=False,
        allow_float=True,
        verbose=False)
data = data.find_doublets(batch_column='sample')
data = data.filter_obs(
            pl.col('doublet').not_() & pl.col('passed_QC'))
data.hvg().normalize().PCA()
data=data.neighbors()
data=data.shared_neighbors()
data = data.cluster(resolution=[0.5, 1.0, 2.0])
data = data.embed()