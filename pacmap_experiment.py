import polars as pl
from single_cell import SingleCell
# import scanpy
import numpy as np
import pacmap
import matplotlib.pyplot as plt
data=SingleCell(f'single-cell/SEAAD/subsampled/SEAAD_raw_20K.h5ad', num_threads=-1)
data=data.qc(subset=False,remove_doublets=False,allow_float=True,verbose=False)
data=data.find_doublets(batch_column='sample')
data=data.filter_obs(pl.col('doublet' ).not_() & pl.col('passed_QC'))
data=data.hvg().normalize()
data=data.PCA()
pca_data=data.obsm['PCs']
pca_data=np.ascontiguousarray(pca_data,dtype=np.float32)
embedding=pacmap.PaCMAP()
output=embedding.fit_transform(pca_data)
# data=data.neighbors().shared_neighbors()
# output.plot_embedding('subclass',f'figures2/original_pacmap.png')
data.obsm['PaCMAP']=output
data.plot_embedding('subclass','figures2/original_pacmap.png')