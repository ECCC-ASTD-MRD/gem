This directory contains files to be used in GEM and related modules
needed to compile GEM-MACH (i.e. "GEM WITH_MACH" via `make cmake-mach`)

file commented out are already incorporated and no longer needed to be
over-written.

These files are here only temporarily until they are incorporated as source file
within the respective GEM and its submodules.

---
current file structure follows that of GEM source file directory setup such that
 following steps are used to overwrite the existing file prior to compilation:

```
# cp ${GEM}/src/mach/tools/overridefiles/Makefile ${GEM}/Makefile
# cp ${GEM}/src/mach/tools/overridefiles/src/CMakeLists.txt ${GEM}/src/CMakeLists.txt
# cp ${GEM}/src/mach/tools/overridefiles/src/rpnphy/src/CMakeLists.txt ${GEM}/src/rpnphy/src/CMakeLists.txt
cp ${GEM}/src/mach/tools/overridefiles/src/rpnphy/src/base/*.F90 ${GEM}/src/rpnphy/src/base
```

# GEMMACH setup from gm-integration-test
Add configuration files/directory to the `${GEM}/share/configurations` directory
```
#cp -a ${GEM}/src/mach/tools/overridefiles/GEM-MACH_3p2_cfgs ${GEM}/share/configurations
```

# following steps to aquire working node, compile and test run interactively with the config setup:
```
qsub -I -lselect=1:ncpus=80:mpiprocs=80:ompthreads=1:mem=160gb
cd ${GEM}
. .eccc_setup_intel
make cmake-mach |& tee make.cmake-mach.out
make work -j |& tee make.work.out
source ${GEM}scripts/link-dbase.sh
cd ${GEM}/work-${GEM_ARCH}/
runprep.sh -dircfg configurations/GEM-MACH_3p2_cfgs
runmod.sh -dircfg configurations/GEM-MACH_3p2_cfgs -ptopo 10x8x1 -inorder >& list_mod_chem
```

----

# notes on the run scripts:
Lee,Vivian (ECCC):

`runprep.sh` is a script that "prepares" the GEM run. It creates the directory
called `PREP` (which has directories `bin`, `input`, `output`, `work`).

The `output/cfg_0000` directory in `PREP` contains mainly analysis and dynamic
fields split in time for `GEM` to read at time intervals.


The files/directories of `PREP/output/cfg_0000` will be linked into
`RUNMOD/input/cfg_0000` by the `runmod.sh` but, `runmod.sh` not only uses these
inputs, it adds additional inputs such as geophysical fields, climatology fields
, ozone (all other inputs needed for the GEM run) plus links to the
`gem_settings.nml`, `outcfg.out`.


It is important to remember that `runprep.sh` is a script that calls
`prep_domain.sh` which is the main script that does the main preparation and we
use it mostly for interactive runs.


In your maestro suites, it uses only `prep_domain.sh`. Under `runprep.sh`, it
will source `configexp.cfg` to determine where to find all the input files.
Under maestro, it uses a file called `experiment.cfg` which you probably noticed
in the suites.

In many operational suites, the simple `runprep.sh` is not used because it
is not constructed to handle all the pre-processing needed before you hit the
`prep_domain.sh` part.  For simple tests, `runprep.sh` is ok.


GEMMACH's integration-test script setup (via @ves001) had a test run setup for
this purposes. I created the simple interactive test (which gives the same
result) by taking the "existing" inputs to make it work.

The original integration test script does compiling and running in one shot.
In this setup, we have to compile and build manually before submitting the
runscript.
