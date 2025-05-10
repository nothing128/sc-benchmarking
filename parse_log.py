import os
import numpy as np
import sys

work_dir = "projects/rrg-wainberg/lamming6/sc-benchmarking"
data_dir = "single-cell/SEAAD/subsampled"


def parse_log(log_file: str, out: str):
    # path = os.path.join(work_dir, log_file)
    path = log_file  # for local testing
    with open(path, "r") as f:
        next(f)
        next(f)
        next(f)  # skip INFO/header
        curr_line = f.readline()
        with open(out, "w") as output:
            output.write("Summary of memory usage per step\n")
            output.write(
                "INFO: Note that if memory used is -1.0, it means the step completed too quickly and no memory check has been run during the step\n"
            )
            output.write(
                "INFO: RSS, VSZ is in KB\n"
            )

            while curr_line != "":  # read until EOF
                if curr_line[0] == "L":
                    params = curr_line.split("(")[1][:-2]
                    output.write(f"Params: {params}\n")
                    curr_line = f.readline()
                    mem = np.array([[-1.0, -1.0, -1.0]])  # RSS, VSZ, %MEM
                    while curr_line[0] != "L":
                        if curr_line[0] != "S":
                            mem = np.append(
                                mem, [curr_line.split(" ")[3:6]], axis=0
                            ).astype(float)
                        else:
                            step_name = curr_line.split(": ")[1].split(" C")[0]
                            max_mem = np.max(mem, axis=0)
                            output.write(f"\t{step_name}:\n\t\tRSS: {int(max_mem[0])}, VSZ: {int(max_mem[1])}, %MEM: {max_mem[2]}\n\n")
                            mem = np.array([[-1.0, -1.0, -1.0]])
                        curr_line = f.readline()

                curr_line = f.readline()


if __name__ == "__main__":
    if len(sys.argv)!=3:
        print("Usage: parse_log.py <log_file> <output file>")
    else:
        parse_log(sys.argv[1],sys.argv[2])
        print(f"output file at {sys.argv[2]}")
