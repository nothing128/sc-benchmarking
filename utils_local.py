# Paste pre-Pyrfected Python
import gc
import gzip
import io
import os
import socket
import subprocess
import time
from contextlib import contextmanager
from timeit import default_timer

import h5py  # type: ignore
import numpy as np
import polars as pl  # type: ignore
import psutil
from scipy import sparse
from scipy.io import mmwrite

delay = 0.1


def write_to_mtx(adata, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    with gzip.open(f"{output_dir}/matrix.mtx.gz", "wb", compresslevel=4) as f:
        mmwrite(f, adata.X.T)
    with gzip.open(f"{output_dir}/barcodes.tsv.gz", "wb", compresslevel=4) as f:
        pl.DataFrame({"barcodes": list(adata.obs_names)}).write_csv(
            f, separator="\t", include_header=False
        )
    with gzip.open(f"{output_dir}/features.tsv.gz", "wb", compresslevel=4) as f:
        pl.DataFrame({"features": list(adata.var_names)}).write_csv(
            f, separator="\t", include_header=False
        )


def write_to_h5(adata, file):
    if not sparse.isspmatrix_csr(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    if "feature_types" not in adata.var.columns:
        adata.var["feature_types"] = ["Gene Expression"] * adata.n_vars
    if "genome" not in adata.var.columns:
        adata.var["genome"] = ["unknown"] * adata.n_vars
    if "gene_ids" not in adata.var.columns:
        adata.var["gene_ids"] = adata.var.index.values.copy()

    def int_max(x):
        return int(max(np.floor(len(str(int(max(x)))) / 4), 1) * 4)

    def str_max(x):
        return max([len(i) for i in x])

    w = h5py.File(file, "w")
    grp = w.create_group("matrix")
    grp.create_dataset(
        "data", data=np.array(adata.X.data, dtype=f"<i{int_max(adata.X.data)}")
    )
    grp.create_dataset(
        "indices", data=np.array(adata.X.indices, dtype=f"<i{int_max(adata.X.indices)}")
    )
    grp.create_dataset(
        "indptr", data=np.array(adata.X.indptr, dtype=f"<i{int_max(adata.X.indptr)}")
    )
    grp.create_dataset(
        "shape",
        data=np.array(
            [adata.n_vars, adata.n_obs],
            dtype=f"<i{int_max([adata.n_vars, adata.n_obs])}",
        ),
    )
    grp.create_dataset(
        "barcodes",
        data=np.array(adata.obs_names, dtype=f"|S{str_max(adata.obs_names)}"),
    )
    ftrs = grp.create_group("features")
    ftrs.create_dataset(
        "feature_type",
        data=np.array(
            adata.var.feature_types, dtype=f"|S{str_max(adata.var.feature_types)}"
        ),
    )
    ftrs.create_dataset(
        "genome",
        data=np.array(adata.var.genome, dtype=f"|S{str_max(adata.var.genome)}"),
    )
    ftrs.create_dataset(
        "id",
        data=np.array(adata.var.gene_ids, dtype=f"|S{str_max(adata.var.gene_ids)}"),
    )
    ftrs.create_dataset(
        "name", data=np.array(adata.var.index, dtype=f"|S{str_max(adata.var.index)}")
    )
    w.close()


class TimerMemoryCollection:
    def __init__(self, silent=True):
        self.timings = {}
        self.silent = silent

    def __call__(self, message):
        start = default_timer()
        pid = os.getpid()

        @contextmanager
        def timer():

            if not self.silent:
                print(f"{message}...")
            curr_process = subprocess.Popen(
                ["./monitor_mem.sh", "-p", str(pid)],
                shell=False,
                stdout=subprocess.PIPE,
                text=True,
            )

            time.sleep(delay)
            try:
                yield
                aborted = False
            except Exception as e:
                aborted = True
                raise e
            finally:
                duration = default_timer() - start

                if not self.silent:
                    status = "aborted after" if aborted else "took"
                    time_str = self._format_time(duration)
                    print(f"{message} {status} {time_str}\n")
                subprocess.run(["kill", str(curr_process.pid)])

                stdout_output = curr_process.communicate()[0]

                mat = np.loadtxt(io.StringIO(stdout_output), delimiter=",").astype(
                    float
                )
                if mat.ndim == 1:
                    mat = mat[None, None]
                max_mat = np.squeeze(np.max(mat, axis=0))
                curr_process.wait()
                curr_process = None
                self.timings[message] = {
                    "duration": duration,
                    "memory": np.round(max_mat[0] / 1024 / 1024, 2),
                    "%mem": np.round(max_mat[1], 1),
                    "aborted": aborted,
                }
                gc.collect()

        return timer()

    def print_summary(self, sort=True, unit=None):
        print("\n--- Timing Summary ---")
        if sort:
            timings_items = sorted(
                self.timings.items(), key=lambda x: x[1]["duration"], reverse=True
            )
        else:
            timings_items = list(self.timings.items())
        total_time = sum(info["duration"] for _, info in timings_items)

        for message, info in timings_items:
            duration = info["duration"]
            memory = info["memory"]
            percentage = (duration / total_time) * 100 if total_time > 0 else 0
            status = "aborted after" if info["aborted"] else "took"
            time_str = self._format_time(duration, unit)
            print(
                f'{message} {status} {time_str} ({percentage:.1f}%) using {memory} GiB ({info["%mem"]}%)'
            )
        print(f"\nTotal time: {self._format_time(total_time, unit)}")

    def _format_time(self, duration, unit=None):
        if unit is not None:
            converted = duration
            if unit == "s":
                pass
            elif unit == "ms":
                converted = duration * 1000
            elif unit == "us" or unit == "µs":
                converted = duration * 1000000
            elif unit == "ns":
                converted = duration * 1000000000
            elif unit == "m":
                converted = duration / 60
            elif unit == "h":
                converted = duration / 3600
            elif unit == "d":
                converted = duration / 86400
            else:
                raise ValueError(f"Unsupported unit: {unit}")
            return f"{converted}{unit}"

        units = [
            (86400, "d"),
            (3600, "h"),
            (60, "m"),
            (1, "s"),
            (0.001, "ms"),
            (0.000001, "µs"),
            (0.000000001, "ns"),
        ]
        parts = []
        for threshold, suffix in units:
            if duration >= threshold or (not parts and threshold == 0.000000001):
                if threshold >= 1:
                    value = int((duration // threshold))
                    duration %= threshold
                else:
                    value = int((duration / threshold) % 1000)
                if value > 0 or (not parts and threshold == 0.000000001):
                    parts.append(f"{value}{suffix}")
                if len(parts) == 2:
                    break
        return " ".join(parts) if parts else "less than 1ns"

    def to_dataframe(self, sort=True, unit=None):
        if not self.timings:
            return pl.DataFrame(
                {
                    "operation": [],
                    "duration": [],
                    "duration_unit": [],
                    "aborted": [],
                    "percentage": [],
                    "memory": [],
                    "%mem": [],
                }
            )

        ops, durs, aborts, pcts = [], [], [], []
        memory, percent_mem = [], []
        memory_unit = []
        total = sum(info["duration"] for info in self.timings.values())

        items = (
            sorted(self.timings.items(), key=lambda x: x[1]["duration"], reverse=True)
            if sort
            else list(self.timings.items())
        )

        duration_unit = "s"
        if unit is not None:
            duration_unit = unit
            conversion = 1.0
            if unit == "s":
                pass
            elif unit == "ms":
                conversion = 1000
            elif unit == "us" or unit == "µs":
                conversion = 1000000
            elif unit == "ns":
                conversion = 1000000000
            elif unit == "m":
                conversion = 1 / 60
            elif unit == "h":
                conversion = 1 / 3600
            elif unit == "d":
                conversion = 1 / 86400
            else:
                raise ValueError(f"Unsupported unit: {unit}")

            for msg, info in items:
                ops.append(msg)
                durs.append(info["duration"] * conversion)
                aborts.append(info["aborted"])
                pcts.append((info["duration"] / total) * 100 if total > 0 else 0)
                memory.append((info["memory"]))
                memory_unit.append("GiB")
                percent_mem.append((info["%mem"]))
        else:
            for msg, info in items:
                ops.append(msg)
                durs.append(info["duration"])
                aborts.append(info["aborted"])
                pcts.append((info["duration"] / total) * 100 if total > 0 else 0)
                memory.append((info["memory"]))
                memory_unit.append("GiB")
                percent_mem.append((info["%mem"]))

        return pl.DataFrame(
            {
                "operation": ops,
                "duration": durs,
                "duration_unit": [duration_unit] * len(ops),
                "aborted": aborts,
                "percentage": pcts,
                "memory": memory,
                "memory_unit": memory_unit,
                "percent_mem": percent_mem,
            }
        )


def system_info():
    hostname = socket.gethostname()
    cpu_count_physical = psutil.cpu_count(logical=False)
    cpu_count_logical = psutil.cpu_count(logical=True)

    mem = psutil.virtual_memory()
    total_mem_gb = mem.total / (1024**3)
    available_mem_gb = mem.available / (1024**3)

    print("\n--- System Information ---")
    print(f"Node: {hostname}")
    print(
        f"CPU: {cpu_count_physical} physical cores, {cpu_count_logical} "
        f"logical cores"
    )
    print(
        f"Memory: {available_mem_gb:.1f} GB available / {total_mem_gb:.1f} " f"GB total"
    )
