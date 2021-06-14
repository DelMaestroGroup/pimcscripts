#! /usr/bin/env python

# Author:      2020-09-25 Adrian Del Maestro
# =============================================================================

'''rmempty.py

Description:
    Removes all empty files from a directory.
'''

import argparse
import numpy as np
import re
import glob,os
import pimchelp

# Get the number of lines in a file
# -----------------------------------------------------------------------------
def get_num_lines(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

# -----------------------------------------------------------------------------
# Begin Main Program
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description='Generates averages from pimc output data.')
    parser.add_argument('-e', '--estimator', type=str, dest='estimator', default = 'estimator',\
                        help='Which estimator should we check for data? [default: estimator]')
    parser.add_argument("base_dir", help='The base directory\
                        where the data files to be removed are located.',
                        default=os.getcwd(), nargs='?')

    args = parser.parse_args()

    estimator = args.estimator

    # Get the list of all estimator file names
    file_names = glob.glob(f'{args.base_dir}{os.sep}*-{estimator}-*.dat')

    # Iterate and delete all empty files
    for cfile in file_names:
        if get_num_lines(cfile) <= 2:
            pimcid = pimchelp.get_pimcid(cfile)
            file_names_to_delete = glob.glob(f'{args.base_dir}{os.sep}*{pimcid}.dat')
            for dfile in file_names_to_delete:
                print(f'Deleting {os.path.basename(dfile)}')
                os.remove(dfile)
            print('')

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__":
    main()
