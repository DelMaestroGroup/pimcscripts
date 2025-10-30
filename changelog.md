Script Changes
==============

## 2025-10-30 v0.4.9
- Stripped absolute path to pimc.e executable from `rsubmit.py` to avoid issues

## 2025-10-30 v0.4.8
- Fixed a bug in the header of `rsubmit.py`

## 2025-10-27 v0.4.7
- Removed some old code related to ssf/ssfq

## 2024-08-27 v0.4.5
- fixed a bug in pimchelp when reading parameters from a file

## 2024-08-20 v0.4.4
- Updated `merge.py` to deal with possibly empty/corrupt files upon open

## 2024-08-20 v0.4.3
- Modification to `reduce-one.py` to allow for weighted averages and standard
  errors when merging over different seeds.

## 2024-08-18 v0.42
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
