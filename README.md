# pimcscripts

This repository includes a number of python scripts that are useful for working
with quantum Monte Carlo data generated via the Del Maestro group PIMC code
which is located at https://code.delmaestro.org. 

## Installation
For now, we have not uploaded to pypi but the scripts can be installed directly
from git:

    pip install --upgrade git+git://github.com/DelMaestroGroup/pimcscripts.git#egg=pimcscripts

This will install the base library `pimcscripts` which includes the modules
`pimcscripts.pimchelp` and `pimcscripts.MCStat` as well as a number of very
useful helper scripts located in [./bin](https://github.com/DelMaestroGroup/pimcscripts/tree/main/pimcscripts/bin). These should be installed to your path and include documentation that can found via:

    script_name.py --help

Below we describe a few useful ones.

## Helper Scripts

### `pimcave.py`

    pimcave
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



