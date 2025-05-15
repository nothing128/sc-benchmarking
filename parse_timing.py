def parse_timing(timing_path: str, mem_path: str,output: str):
    """
    timing_path should be a csv file
    mem_path should be a txt file generated from parse_log.py
    """
    with open(timing_path, 'r') as timing, open(mem_path, 'r') as mem,open(output, 'w') as f:
        next(mem)
        next(mem)
        next(timing)
        for line in mem:
            if line=="" or line =="\n":
                pass
            elif line[0:6]=="params":
                f.write(f"\n{line}")
            else:
                curr_line=timing.readline()
                if curr_line[0:2]=="Ne":
                    curr_line=timing.readline()
                time_taken=f"{round(float(curr_line.split(',')[1]),2)}{curr_line.split(',')[2]} ({round(float(curr_line.split(',')[4]),1)}%)"
                f.write(f"{line[:-1]}in {time_taken}\n")               
                
if __name__=="__main__":
    parse_timing("timing_output.csv","output.txt","final.txt")