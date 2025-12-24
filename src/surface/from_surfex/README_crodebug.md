
CROCUS debugging option in GEM and SPS
=====================================

For GEM and SPS, activate
```
svs2_crodebug =.true.
```
surface option in namelist and then follow official documentation below.

***Before launching integration export CROCUS debugging environement variables
documented below.***

TIPS:
- Ideally run with `ptopo=1x1` given that the debugging is printed with a simple `print *` and some messages will be printed once per cpu asked
- Ideally run with specific lat/lon - but will still get duplication of some information if `cpu>1`


Documentation for CROCUS debugging
==================================

Original documentation from CNRM/Meteo France: https://www.umr-cnrm.fr/surfex/IMG/pdf/doc_crodebug.pdf

### 1. Print daily informations

To print daily informations, you must define the environment variable `CROCUS_INFO` before running OFFLINE :
```
export CROCUS_INFO=1
```
This will print the number of snow layers for a given point (see section 3),
the snow fraction for the same point, and the number of points with snow on
the ground.

### 2. Print detailed snow profiles at each time step

To print detailed informations necessary for debugging, you must define the
`CROCUS_DEBUG` environment variable before running OFFLINE :
- `export CROCUS_DEBUG=1` to print a vertical profile at the beginning
and at the end of each time step
- `export CROCUS_DEBUG=2` to print a vertical profile after each CROCUS subroutine
- `export CROCUS_DEBUG=3` same behaviour as `CROCUS_DEBUG=1` + print meteorological forcing data
- `export CROCUS_DEBUG=4` same behaviour as `CROCUS_DEBUG=2` + print meteorological forcing data
- `export CROCUS_DEBUG=5` to print a very detailed vertical profile (including variables GRAN1, GRAN2, HIST, AGE) after each CROCUS subroutine

### 3. Choose the ouput point for multi-points simulations

If `CROCUS_INFO` and/or `CROCUS_DEBUG` are defined, you may optionally
define other variables to choose the point which will be printed :

- `export CROCUS_DEBUG_POINT=X` to choose a point by its indice inside CROCUS routine. Be careful, if some areas are snow-free, the point number can vary from a time step to another.
- `export CROCUS_DEBUG_LAT=X`
- `export CROCUS_DEBUG_LON=X` to always print the same geographical
point, defined by latitude and longitude. (recommended)

If there is no spatial information (these variables are not defined), profiles are printed for the first simulation point.

### 4. Choose the dates for printing

-  To activate the debugging mode after a given date : `export CROCUS_DEBUG_DATE=YYYYMMJJ` (recommended if the bug does not occur in the first simulation days)
-  To specify the hour for activating the debugging mode : `export CROCUS_DEBUG_HOUR=HH`
- To stop the debugging mode after a given date : `export CROCUS_DEBUG_DATE_END=YYYYMMJJ`
 
If there is no temporal information, profiles are printed from the first timestep to the end of the run (or the crash). Be careful to not activate `CROCUS_DEBUG` for long simulations without temporal limitations to avoid huge
outputs.