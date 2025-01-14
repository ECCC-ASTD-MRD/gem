# How to get, compile and run GEM at the CMC.

For general information, a setup script and information on Maestro, see
README_general and share/gem-maestro/README files

Warning: this repository uses submodules. Make sure you follow the
instructions below.

## Getting gem git repository

### Choose one of the following methods:

1. cloning only the necessary components:
```
git clone git@gitlab.science.gc.ca:MIG/gem.git
cd gem
```

2. or cloning everything, including rpn-si libraries (rmn, vgrid, rpncomm, tdpack) in one step:
```
git clone --recursive git@gitlab.science.gc.ca:MIG/gem.git
cd gem
```

3. or cloning in several steps
```
git clone git@gitlab.science.gc.ca:MIG/gem.git
cd gem
```
Update only cmake_rpn submodules, for example if you want to modify the default compilation flags
```
git submodule update --init cmake_rpn
```
Update everything: rpn-si libraries, utilities and cmake_rpn submodules
```
git submodule update --init --recursive
```

## Choosing a version

```
git branch # what is the current branch
git branch -a # list all branches (look at the list of remote branches to choose from)
git tag # list tags (if you want to select a tagged version)
git checkout <hash|branch|tag> # checkout a branch, a tag, or a specific hash.  Example: git checkout 5.3
```
Before making changes, create your own branch from the current branch
```
git checkout -b mybranch
```

## Linking to GEM database (to be done once)
```
./scripts/link-dbase.sh
```

## Preparing gem compilation for Intel compiler suite
```
. ./.eccc_setup_intel
```

## Or preparing gem compilation for gnu compiler suite

Please note you cannot compile with Intel and then with GNU in the same shell
```
. ./.eccc_setup_gnu
```

Before the first build, or if you made important changes (such as updating
other submodules, or adding or removing source files):
```
. ./.initial_setup
```

### Scripts

Scripts in `scripts/support` and `scripts/rpy` directories are a copy of scripts
already loaded from SSM domains when a `.eccc_setup*` file is called.  By
default, they are not used, but if you want to test or modify them, you can
override SSM scripts by setting `GOAS_SCRIPT_MODE` variable before sourcing
`.eccc_setup_intel` or `.eccc_setup_gnu`:

```
export GOAS_SCRIPT_MODE=true
```

Please also note that if you load maestro, maestro scripts will be used,
either in a maestro suite or when running GEM interactively. Otherwise, goas
task setup files situated in the `scripts` directory will be used instead.

## Building and installing GEM

There is a script called `cado` aimed at replacing the top-level Makefile.
Both coexist, even if we suggest you use the cado script.
See `cado -h` (short help) or `cado help` or the content of the Makefile for
options.
For example: `cado cmake` or `make cmake` generates Makefiles to compile
gem, gemdyn, modelutils and rpnphy.  The cmake command used by cado script
is printed at the end of the process.

Configure for *GEM*:
```
cado cmake
or 
make cmake
```

Or, if you want to compile *GEM-MACH instead of GEM*, use the following instead:
```
cado cmake-mach
```

Compile:
```
cado build -j
or
make -j
```
Install in working directory
```
cado work -j
make -j work
```
cado work -j or make -j work can be used to compile and install in the same step.

In development mode, GEM is compiled using Intel shared libraries: use the
following command to compile with static libraries (for GEM or GEM-MACH):
```
cado cmake-static
or
cado cmake-mach-static
```

See others options with cado -h (short help) or cado help

## Running GEM: example

```
cd $GEM_WORK
runprep.sh -dircfg ./configurations/GEM_cfgs_LU_FISL_H
runmod.sh -dircfg ./configurations/GEM_cfgs_LU_FISL_H
```

If you come back later, and you want to run the executables you compiled
before, you just need to use the following command before going into the
$GEM_WORK directory:
```
. ./.eccc_setup_intel
or, if you compiled with gnu:
. ./.eccc_setup_gnu
and then:
cd $GEM_WORK
```

## Some tips for compilation

When the `cado cmake` command is called, information is printed, among which
the list of compilation flags used, such as (example with Intel on science
side):
```
-- (EC) CMAKE_C_FLAGS=-fp-model precise -traceback -Wtrigraphs -xICELAKE-SERVER -diag-disable=10441 -qmkl 
-- (EC) CMAKE_Fortran_FLAGS=-convert big_endian -align array32byte -assume byterecl -fp-model source -fpe0 -traceback -stand f08 -xICELAKE-SERVER -diag-disable=5268,7025,7373 -qmkl -static-intel
```

If you choose the debug version (`cado cmake-debug`), some flags are added to the previous ones, and, again, printed when `cado cmake-debug` is called:
```
-- (EC) CMAKE_C_FLAGS_DEBUG=-O0 -g -ftrapuv
-- (EC) CMAKE_Fortran_FLAGS_DEBUG=-O0 -g -ftrapuv
```
With `cado cmake-debug-extra`:
```
-- (EC) CMAKE_C_FLAGS=-fp-model precise -traceback -Wtrigraphs -xICELAKE-SERVER -diag-disable=10441 -Wall -qmkl 
-- (EC) CMAKE_Fortran_FLAGS=-convert big_endian -align array32byte -assume byterecl -fp-model source -fpe0 -traceback -stand f08 -xICELAKE-SERVER -diag-disable=5268,7025,7373 -warn all -check all -qmkl -static-intel
```

*Important note*: if you want to change the compilation type, for example, first, you compiled with the debug version (`cado cmake-debug`), and then you want to use the release version (`cado cmake`), you need to remove the contents of the build directory between these two commands. You can use the following command: `. ./.initial_setup` which will empty the build and work directories, and then you can proceed from the start with the `cado cmake` configure command.
 
The compilation flags come from default compiler rules set up by RPN-SI and are applied to all the compilation processes.

If you want to change those flags, you can either:
- update the `cmake_rpn` submodule so that you can edit the files and modify those flags directly:
  - `git submodule update --init cmake_rpn`
  - make the changes in the file corresponding to the platform and compiler used, such as:
    `cmake_rpn/modules/ec_compiler_presets/ECCC/rhel-8-icelake-64/inteloneapi-2022.1.2.cmake`
- or edit the `CMakeLists.txt` file and add the flags at the end of the following lines (for Intel):
```
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qmkl ${STATIC_LINK_INTEL_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qmkl -static-intel -diag-disable 5268 ${STATIC_LINK_INTEL_FLAGS}")
```

If you want to change or add flags for a specific part of GEM, for example rpnphy, you can either:
- change the flags for all sources, by editing the `src/rpnphy/CMakeLists.txt`
  file and add the flags at the end of the following lines (for Intel):
```
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qmkl")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qmkl -static-intel -diag-disable 5268")
```
- or, if you want to change or add flags for a specific source file only,
  edit the `CMakeLists.txt` file situated in the directory where this source
  file is added.
  For example, for the source file `rpnphy/src/utils/sfclayer.F90`, edit the
  `rpnphy/src/CMakeLists.txt`, and modify the following line according to your
  needs (here we are adding the -C flag to the default flags:
```
set_source_files_properties(utils/sfclayer.F90 PROPERTIES COMPILE_OPTIONS "-C")
```

## Structure of the working environment

The structure of the build and work directories is different whether the
$storage_model environment variable exists:

The following environment variables are created (examples):
- gem_DIR = directory where the git clone was created
- GEM_WORK = work directory
- GEM_ARCH = architecture, for example ubuntu-22.04-amd64-64-intel-2022.1.2
- COMPILER_SUITE = compiler suite, for example Intel
- COMPILER_VERSION = compiler version, for example 2022.1.2

- GEM_STORAGE_DIR = where build and work directories are situated
  - Example if $storage_model variable exists:
    - GEM_STORAGE_DIR=/local/storage/gem/ubuntu-22.04-amd64-64-intel-2022.1.2
    - in gem_DIR:
      - build-ubuntu-22.04-amd64-64-intel-2022.1.2 is a link, such as:
        /local/storage/gem/ubuntu-22.04-amd64-64-intel-2022.1.2/build
      - work-ubuntu-22.04-amd64-64-intel-2022.1.2 is a link, such as:
        /local/storage/gem/ubuntu-22.04-amd64-64-intel-2022.1.2/work

  - Example if $storage_model variable doesn't exist:
    - GEM_STORAGE_DIR=$HOME/gem/
    - directories situated in gem_DIR:
      - build-ubuntu-22.04-amd64-64-intel-2022.1.2
      - work-ubuntu-22.04-amd64-64-intel-2022.1.2

## Structure of GEM source code

```
--------------- gemdyn ---------------- TOP
---- rpnphy ---- mach ------- cpl -----  |
------------- modelutils --------------  |
--------------- vgrid -----------------  |
----------rmn---tdpack--rpncomm -------  V
--------- compiler libraries ---------- BOTTOM
```

GEM code is built top-down meaning:

the routines in a library *cannot* call any functions or use any modules above it

ie: any routines in rpnphy *cannot* call a routine or use a module in gemdyn

ie: any routines in gemdyn can call anything below it
