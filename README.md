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

2. or cloning everything, including rpn-si libraries (rmn, vgrid, rpncomm, tdpack) in one step
   (necessary, for now, to compile with gnu)
```
git clone --recursive git@gitlab.science.gc.ca:MIG/gem.git
cd gem
```

3. or cloning in several steps
```
git clone git@gitlab.science.gc.ca:MIG/gem.git
cd gem
```
Update rpn-si libraries and cmake_rpn submodules
```
git submodule update --init --recursive
```

## Choosing a version

```
git branch # what is the current branch
git branch -a # list all branches (look at the list of remote branches to choose from)
git tag # list tags (if you want to select a tagged version)
git checkout <hash|branch|tag> # checkout a branch, a tag, or a specific hash. Example: git checkout 5.3
```
Or, if you want, create your own branch from the current branch
```
git checkout -b mybranch
```

## Preparing gem compilation for Intel compiler
```
./scripts/link-dbase.sh
. ./.eccc_setup_intel
```

### Or for gnu

Please note you cannot compile with Intel and then with GNU in the same shell
```
. ./.eccc_setup_gnu
```

Before the first build, or if you made important changes (such as updating
other submodules, or adding or removing source files):
```
. ./.initial_setup
```

## Building and installing GEM

There is a new script aimed at replacing the top-level Makefile.
For now, both still coexist.
See cado -h (short help) or cado help or the content of the Makefile for options.
For example: cado cmake or make cmake generates Makefiles to compile gem,
gemdyn, modelutils and rpnphy

```
cado cmake
or 
make cmake
```
Compile
```
cado build -j
or
make -j
```
install in working directory
```
cado work -j
make -j work
```
cado work -j or make -j work can be used to compile and install in the same step.

In development mode, gem is compiled using Intel shared libraries: use the
following command to compile with static libraries:
```
cado cmake-static
```

See others options with cado -h (short help) or cado help


## Running GEM: example

```
cd $GEM_WORK
runprep.sh -dircfg ./configurations/GEM_cfgs_LU_FISL_H
runmod.sh -dircfg ./configurations/GEM_cfgs_LU_FISL_H
```

## Structure of the working environment
The structure of the build and work directories is different whether the
$storage_model environment variable exists:

The following environment variables are created (examples):
- gem_DIR = directory where the git clone was created
- GEM_WORK = work directory
- GEM_ARCH = architecture, for example ubuntu-18.04-amd64-64-intel-2022.1.2
- COMPILER_SUITE = compiler suite, for example Intel
- COMPILER_VERSION = compiler version, for example 2022.1.2

- GEM_STORAGE_DIR = where build and work directories are situated
  - Example if $storage_model variable exists:
    - GEM_STORAGE_DIR=/local/storage/gem/ubuntu-18.04-amd64-64-intel-2022.1.2
    - in gem_DIR:
      - build-ubuntu-18.04-amd64-64-intel-2022.1.2 is a link, such as:
        /local/storage/gem/ubuntu-18.04-amd64-64-intel-2022.1.2/build
      - work-ubuntu-18.04-amd64-64-intel-2022.1.2 is a link, such as:
        /local/storage/gem/ubuntu-18.04-amd64-64-intel-19.0.3.199/work

  - Example if $storage_model variable doesn't exist:
    - GEM_STORAGE_DIR=$HOME/gem/
    - directories situated in gem_DIR:
      - build-ubuntu-18.04-amd64-64-intel-2022.1.2
      - work-ubuntu-18.04-amd64-64-intel-19.0.3.199
