#!/usr/bin/env python3
#Filename: make_pbs_resubmit_script.py

import argparse
from os import listdir
from make_pbs_script import submit

def generate_resubmit_vacc(pimc_output_dir="./OUTPUT/",pbs_script_name="resubmit.pbs",job_name="JOBNAME"):
    # add "/" to end of directory if missing 
    if pimc_output_dir[-1] != "/":
        pimc_output_dir += "/"
    
    # get log files
    output_files = listdir(pimc_output_dir)
    log_files = [f for f in output_files if "-log-" in f]
    
    # get restart commands from log files
    restart_commands = []
    for log_file in log_files:
        with open(pimc_output_dir + log_file,"r") as f:
            f.readline()
            f.readline()
            cmd = f.readline()[2:]
            restart_commands.append(cmd)
    
    # generate resubmit script
    resubmit = submit("resubmit.sh",pbs_script_name)
    resubmit.pbs_script_from_list(restart_commands,jobName=job_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument( "-D",
                         "--pimc_output_dir",
                         type=str,
                         help="output directory containing pimc log files",
                         default="./OUTPUT/")

    parser.add_argument( "-o",
                         "--pbs_script_name",
                         type=str,
                         help="filename for PBS resubmit script",
                         default="resubmit.pbs")

    parser.add_argument( "-N",
                         "--job_name",
                         type=str,
                         help="name of PBS job",
                         default="JOBNAME")
    args = parser.parse_args()

    generate_resubmit_vacc( pimc_output_dir=args.pimc_output_dir,
                            pbs_script_name=args.pbs_script_name,
                            job_name=args.job_name
                            )


