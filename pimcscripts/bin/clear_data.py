'''clear_data.py

Description:
Clear all data while leaving headers from a list of estimator files.

Usage: clear_data.py (--estimator=<E>)

Options:
  -h, --help                        Show this help message and exit
  -e <E>, --estimator=<E>           The name of the estimator file, e.g. "estimator"
'''

import glob
from docopt import docopt


# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 
    # parse the command line options
    args = docopt(__doc__)
    estimator_type = args['--estimator']

    # get the lsit of files
    file_list = glob.glob('*-%s-*' % estimator_type)

    if not file_list:
        print("Estimator %s-* not found" % estimator_type)

    # for each estimator file
    for file in file_list:
        
        # get the data
        in_file = open(file,'r')
        in_lines = in_file.readlines()
        in_file.close()
        
        # re-write only the header lines
        with open(file, 'w') as out_file:
            for line in in_lines:
                if line[0] == "#":
                    out_file.write(line)
                else:
                    break

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

if __name__ == "__main__": 
    main()

