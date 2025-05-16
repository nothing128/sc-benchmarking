import os
import numpy as np
import sys

work_dir = "projects/rrg-wainberg/lamming6/sc-benchmarking"
data_dir = "single-cell/SEAAD/subsampled"


def parse_log(log_file: str, out: str):
    # path = os.path.join(work_dir, log_file)
    path = log_file  # for local testing
    with open(path, "r") as f:
        curr_line = f.readline()
        with open(out, "w") as output:
            output.write("Summary of memory usage per step\n")
            output.write(
                "INFO: Note that if memory used is -1.0, it means the step completed too quickly and no memory check has been run during the step\n"
            )

            while curr_line != "":  # read until EOF
                if curr_line[0]=="L":
                    params=curr_line.split(': (')[1][:-1]
                    output.write(f"\nparams: {params}\n")
                    curr_line = f.readline()
                else:
                    while curr_line[0] != "L":
                        next(f)
                        next(f)
                        curr_line=f.readline()
                        mem = np.array([[-1.0,-1.0]])
                        while curr_line[-2] != ".":
                            print(curr_line)
                            mem=np.append(mem, [[curr_line.split(' ')[3],curr_line.split(' ')[5]]],axis=0).astype(float)
                            curr_line=f.readline()
                        curr_line=f.readline()
                        step = curr_line.split(": ")[1].split(" C")[0]
                        max_mem=np.max(mem, axis=0)
                        output.write(f"{step}: {np.round((max_mem[0]/1024/1024),2)} GB ({max_mem[1]}%) \n")
                        curr_line=f.readline()
                    curr_line=f.readline()

                    
                    
                    
                
                


if __name__ == "__main__":
    if len(sys.argv)!=3:
        print("Usage: parse_log.py <log_file> <output file>")
    else:
        parse_log(sys.argv[1],sys.argv[2])
        print(f"output file at {sys.argv[2]}")
