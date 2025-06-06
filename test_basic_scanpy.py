import gc
import sys
import polars as pl  # type: ignore
import scanpy as sc  # type: ignore
import matplotlib.pyplot as plt  # type: ignore

work_dir = "projects/sc-benchmarking"
data_dir = "single-cell/SEAAD/subsampled"

sys.path.append(work_dir)
from utils_local import TimerMemoryCollection, system_info

size_options = ["20K", "400K", "1M"]
size_choice = sys.argv[1]
assert len(sys.argv) == 2
assert size_choice in size_options

system_info()
timers = TimerMemoryCollection(silent=False)

# with timers('Load data (10X mtx)'):
#     data = sc.read_10x_mtx(f'{data_dir}/SEAAD_raw_{size}')
# del data; gc.collect()

# with timers('Load data (h5)'):
#     data = sc.read_10x_h5(f'{data_dir}/SEAAD_raw_{size}.h5')
# del data; gc.collect()

# Note: Loading is much slower from $SCRATCH disk

with timers("Load data (h5ad/rds)"):
    data = sc.read_h5ad(f"{data_dir}/SEAAD_raw_{size_choice}.h5ad")

with timers("Quality control"):
    data.var["mt"] = data.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(data, qc_vars=["mt"], inplace=True, log1p=True)
    sc.pp.filter_cells(data, min_genes=100, copy=False)
    sc.pp.filter_genes(data, min_cells=3, copy=True)

with timers("Doublet detection"):
    sc.pp.scrublet(data, batch_key="sample")

with timers("Normalization"):
    sc.pp.normalize_total(data)
    sc.pp.log1p(data)

with timers("Feature selection"):
    sc.pp.highly_variable_genes(data, n_top_genes=2000, batch_key="sample")

with timers("PCA"):
    sc.tl.pca(data)

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

with timers("Plot embeddings"):
    sc.pl.umap(data, color=["leiden_res_1.00"])
    plt.savefig(
        f"{work_dir}/figures/scanpy_embedding_cluster_{size_choice}.png",
        dpi=300,
        bbox_inches="tight",
        pad_inches="layout",
    )

with timers("Find markers"):
    sc.tl.rank_genes_groups(data, groupby="leiden_res_1.00", method="wilcoxon")

timers.print_summary(sort=False)
timers_df = timers.to_dataframe(sort=False, unit="s").with_columns(
    pl.lit("test_basic_scanpy").alias("test"), pl.lit(size_choice).alias("size")
)
# increments the output csv file to ensure old outputs do not get overwritten
print(timers_df)
output = f"{work_dir}/output/test_basic_scanpy_{size_choice}.csv"
timers_df.write_csv(output)

del timers, timers_df, data
gc.collect()

