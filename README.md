# pimcscripts

This repository includes a number of python scripts that are useful for working
with quantum Monte Carlo data generated via the Del Maestro group PIMC code
which is located at https://code.delmaestro.org. 

## Installation
For now, we have not uploaded to pypi but the scripts can be installed directly
from git:

    pip install --upgrade git+https://github.com/DelMaestroGroup/pimcscripts.git#egg=pimcscripts

This will install the base library `pimcscripts` which includes the modules
`pimcscripts.pimchelp` and `pimcscripts.MCStat` as well as a number of very
useful helper scripts located in [./bin](https://github.com/DelMaestroGroup/pimcscripts/tree/main/pimcscripts/bin). These should be installed to your path and include documentation that can found via:

    script_name.py --help

If you are upgrading after a change, it might be useful to use:

    pip install --upgrade --no-deps --force-reinstall git+https://github.com/DelMaestroGroup/pimcscripts.git#egg=pimcscripts

Below we describe a few useful ones.

## Helper Scripts

### `pimcave.py`

    usage: pimcave.py [-h] [-s SKIP] [-l HEADER_LINES] [-r] file [file ...]

    Generates averages from pimc output data.

    positional arguments:
      file                  File or files to average.

    optional arguments:
      -h, --help            show this help message and exit
      -s SKIP, --skip SKIP  How many input lines should we skip? [default: 0]
      -l HEADER_LINES, --header_lines HEADER_LINES
                            How many header lines to skip?
      -r, --repeated_header
                            deal with duplicate headers

### `pimcplot.py`

    pimcplot

    Description:
      Performs a cumulative average plot of raw Monte Carlo data

    Usage:
        pimcplot.py [options] [--legend=<label>...] --estimator=<name> <file>...

        pimcplot.py -h | --help

    Options:
      -h, --help                    Show this screen.
      --estimator=<name>, -e <name> The estimator to be plotted.
      --skip=<n>, -s <n>            Number of measurements to be skipped [default: 0].
      --period=<m>, -p <m>          The period of the average window [default: 50].
      --truncateid=<t>, -t <t>      Truncate PIMCID to last <t> characters [default: 0].
      --legend=<label>, -l <label>  A legend label
      --period=<m>, -p <m>          The period of the average window [default: 50].
      --error=<units>, -d           Size of the error bars
      --nobin                       Don't use the binned errorbars
      --nolegend                    Turn off the legend
      --ttest                       Perform a ttest
      --hline=<val>                 Include a horizontal line at <val> in the averaged plot
      --hlabel=<hl>                 A legend label for the horizontal line.
      --title=<title>               A title for the plots.
      --savefig=<figure>            A filename for saved plots (extensions supported by active matplotlib backend).
      --quiet                       Suppress output

### `reduce-one.py`

    usage: reduce-one.py [-h] [-T T] [-N N] [-n N] [-t TAU] [-u MU] [-L L] [-V V] -r {T,N,n,u,t,L,V,M} [--canonical] [-R R] [-s SKIP] [-e ESTIMATOR] [-i PIMCID] [base_dir]

    Reduce quantum Monte Carlo output over some parameter.

    positional arguments:
      base_dir              The base directory where the data files to be reduced are located.

    optional arguments:
      -h, --help            show this help message and exit
      -T T, --temperature T
                            simulation temperature in Kelvin
      -N N, --number-particles N
                            number of particles
      -n N, --density N     number density in Angstroms^{-d}
      -t TAU, --imag-time-step TAU
                            imaginary time step
      -u MU, --chemical-potential MU
                            chemical potential in Kelvin
      -L L, --Lz L          Length in Angstroms
      -V V, --volume V      volume in Angstroms^d
      -r {T,N,n,u,t,L,V,M}, --reduce {T,N,n,u,t,L,V,M}
                            variable name for reduction [T,N,n,u,t,L,V,W,M]
      --canonical           are we in the canonical ensemble?
      -R R, --radius R      radius in Angstroms
      -s SKIP, --skip SKIP  number of measurements to skip [0]
      -e ESTIMATOR, --estimator ESTIMATOR
                            specify a single estimator to reduce
      -i PIMCID, --pimcid PIMCID
                            specify a single pimcid
