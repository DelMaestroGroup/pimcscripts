Script Changes
==============
## 2024-08-18
- Modification to `merge.py` to allow for pre-averaging of seeds using the
  `--seed` cmdline argument.

## 2023-12-19
- Added considerable functionality to use multiple estimators and multiple
  values.  Useful for comparing two different measurements of the same
  estimator.
- Updated pimcave.py to take a single estimator

## 2023-11-15
- Updated merge.py to use less memory for cumulative estimators

## 2023-4-23
- Added parallel to gensubmit.py and rsubmit.py to run parallel tasks on the cluster easier

## 2022-01-10 
- Updated reduce-one.py to deal with "average" estimators

## 2021-12-29
- Added custom output labels to merge.py

## 2021-12-17 
- Adding ssfq files to merge.py

## 2021-12-13 0.3.3
- Added some LJ log parameters to parameter map

## 2021-12-12 0.3.2
- Added a q-resolved `ssfq` to reduce-one

## 2021-12-06 0.3.1
- added static structure factor (`ssf`) to reduce-one

## 2021-08-25
- Bugfix in `pimchelp.get_pimcid()` which only shows up when the chemical
  potential is negative
