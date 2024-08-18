#!/usr/bin/env python

# Filename: convertPIMCuuid.py
# Author: Nathan Nichols
#
# Convert the PIMC datafiles to new UUID.
import uuid
import re
import shlex
import subprocess
import argparse
def ArgParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("DIRECTORY", nargs='?', help="Directory containing PIMC data files to be converted.", type=str, default="./OUTPUT/")
    args = parser.parse_args()
    return args

def convertPIMCuuid(data="./OUTPUT"):
    ls = subprocess.Popen(["ls", data], stdout=subprocess.PIPE)
    process = subprocess.Popen(["grep", "log"], stdin=ls.stdout, stdout=subprocess.PIPE)
    stdout, stderr = process.communicate()
    logFiles = stdout.decode('utf-8').splitlines()

    for f in logFiles:
        fn = f.split(sep='log')
        oldFile = re.split(r'\-',f)
        _id = str(uuid.uuid4())
        _oldid = oldFile[-1].split(sep='.')[0]
        c = "sed -i '3 s/-R {}/-R {}/' {}".format(_oldid,_id,data + f)
        p = subprocess.Popen(shlex.split(c))
        for s in ['log','estimator']:
            c = "sed -i '1 s/^.*$/# PIMCID: {}/' {}".format(_id,data + fn[0] + s + fn[1])
            p = subprocess.Popen(shlex.split(c))
        for s in ['log','estimator','state']:
            oldFile[1] = s
            newFile = oldFile.copy()
            newFile[-1] = _id+'.dat'
            c = "mv {} {}".format(data + '-'.join(oldFile), data + '-'.join(newFile))
            p = subprocess.Popen(shlex.split(c))

if __name__ == '__main__':
    args = ArgParser()
    convertPIMCuuid(args.DIRECTORY)
