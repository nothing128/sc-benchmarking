import polars as pl
from single_cell import SingleCell
import scanpy
import matplotlib.pyplot as plt
data=SingleCell(f'single-cell/SEAAD/subsampled/SEAAD_raw_20K.h5ad', num_threads=-1)
data=data.qc(subset=False,remove_doublets=False,allow_float=True,verbose=False)
data=data.find_doublets(batch_column='sample')
data=data.filter_obs(pl.col('doublet' ).not_() & pl.col('passed_QC'))
curr_data=data.hvg().normalize().PCA()

# data=data.to_scanpy()

# for i in [6,8,10,12,14]:
#     for j in [26,28,30,32,34]:
#         data=curr_data
#         data=data.neighbors(num_neighbors=i+j).shared_neighbors(min_shared_neighbors=3)
#         data=data.cluster(min_cluster_size=5,resolution=[0.5,1.0,2.0])
#         data=data.embed(num_neighbors=i,num_extra_neighbors=j,num_mid_near_pairs=int(i*0.5), num_further_pairs=2*i, num_iterations=450)
#         data.plot_embedding('subclass',f'figures2/{i}neighbors_{j}extras.png')
#         print(i,j)
for i in [0.1,0.5,1,1.5,2,2.5]:
    for j in [0.5,1,2,3,4,5]:
        data=curr_data
        data=data.neighbors(num_neighbors=44).shared_neighbors()
        data=data.cluster(resolution=[0.5,1.0,2.0])
        data=data.embed(num_neighbors=10,num_extra_neighbors=34,num_mid_near_pairs=int(10*i),num_further_pairs=int(j*10))
        data.plot_embedding('subclass',f'figures2/{i}mid_near_{j}further.png')
        print(i,j)