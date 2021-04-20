The following instructions explain how to use GEM on a stick at the CMC.

# Instructions in a nutshell (for Ubuntu 18)

See below for extended instructions.

```
git clone git@gitlab.science.gc.ca:cpi001/gem-devel.git
cd gem-devel
. ./.eccc_setup_intel
```
# or, if you want to compile with GNU:
```
. ./.eccc_setup_gnu
```

Before the first build:
```
. ./.eccc_initial_setup
./link-dbase.sh
```

# build and install GEM
```
make cmake
make build
make work
```

# run GEM: examples
```
cd $GEM_WORK
runprep.sh -dircfg ./configurations/GEM_cfgs_LU_FISL_P
runmod.sh -dircfg ./configurations/GEM_cfgs_LU_FISL_P
```

The following environment variables are created (examples):
GEM_DIR=/local/drive1/armn/armncpi/storage/goas/gem-devel_ubuntu-18.04-amd64-64-intel-19.0.3.199
GEM_GIT_DIR=/users/dor/armn/cpi/ords/gem/gem-devel
GEM_WORK=work-ubuntu-18.04-amd64-64-intel-19.0.3.199
GEM_ARCH=ubuntu-18.04-amd64-64-intel-19.0.3.199
COMPILER_suite=intel
COMPILER_version=19.0.3.199

The structure of the build and work directories is different whether the $storage_model environment variable exists:

Example if $storage_model exists:
GEM_DIR=/local/drive1/armn/armncpi/storage/goas/gem-devel_ubuntu-18.04-amd64-64-intel-19.0.3.199
in GEM_GIT_DIR:
build-ubuntu-18.04-amd64-64-intel-19.0.3.199 -> /local/drive1/armn/armncpi/storage/goas/gem-devel_ubuntu-18.04-amd64-64-intel-19.0.3.199/build
work-ubuntu-18.04-amd64-64-intel-19.0.3.199 -> /local/drive1/armn/armncpi/storage/goas/gem-devel_ubuntu-18.04-amd64-64-intel-19.0.3.199/work

Example if $storage_model doesn't exist:
GEM_DIR=/users/dor/armn/cpi/ords/gem/gem-devel
build-ubuntu-18.04-amd64-64-intel-19.0.3.199
work-ubuntu-18.04-amd64-64-intel-19.0.3.199

```
-----------------------------------------------------------------
# Extended instructions:

# Configuring the environment

If you do not know the version of a system on which you are connected,
please run the following command: `lsb_release -a`

## Ubuntu 18.04

The Intel 19 compiler and gnu compilers 7.5 are available on Ubuntu 18.04
systems.  Here is an incomplete list of such systems:

- ppp3.science.gc.ca
- ppp4.science.gc.ca
- daley.science.gc.ca
- banting.science.gc.ca
- rutherford.cmc.ec.gc.ca
- Newer workstations on the cmc.ec.gc.ca domain

To configure your environment for the Intel 19 compiler, from the root of your
GOAS clone, source **.eccc_setup_intel_19**. 

```
   . eccc_setup_intel_19
```

To configure your environment for the gnu compilers, from the root of your
GOAS clone, source **.eccc_setup_gnu_7**.

```
   . eccc_setup_gnu_7
```

It will define environment variables and load the appropriate SSM packages.

## Ubuntu 14.04

The Intel 16 compiler and gfortran 5.1 are available on Ubuntu 14.04 systems.
Older worstations on the cmc.ec.gc.ca domain are on such systems.

To configure your environment for the Intel 16 compiler, from the root of
your GOAS clone, source **.eccc_setup_intel_16**.

```
   . eccc_setup_intel_16
```

To configure your environment for the gnu compilers, from the root of your
GOAS clone, source **.eccc_setup_gnu_5**.  

```
   . eccc_setup_gnu_5
```

It will define environment variables, and load the appropriate SSM packages.

Before the first build, from the root of your GOAS clone, source
**.initial_eccc_setup**, it will create a build directory with the
architecture, compiler and version in its name.
```
. ./eccc_initial_setup
```

Then, link to gem database, and copy some extra files (they will be added
later in the common database).
```
./link-dbase.sh
```

# Compiling and executing the model

Once the environment has been configured, there are a few extra steps to do
depending on how you want to host you clone of the repository.

## Repository in a data file system

The simplest way to work with GOAS is to have it in a data file system such
as ORDS or a local file system on your workstation.  The most important
downside to this approach is that there are no backups or snapshots of the
code you have modified.  Also, if you choose a local file system, it might
not be accessible from all hosts.

To ensure that your modifications are backed-up, you can make sure to regularly
push them to a remote Git repository or to implement you own backup scheme.

When sourced, the `.initial_eccc_setup` will create a build directory following the
pattern below:
`<goas_root>/build-<os-platform>-<compiler>-<compiler-version>`

for example: 
`<goas_root>/build-ubuntu-18.04-amd64-64-intel-19.0.3.199`

## Repository in the $HOME directory

Hosting your GOAS working copy in your home directory has a few advantages,
but it is also slightly more complex.  Being in the home directory, your
repository is automatically backed-up and you also have multiple snapshots.
The downside is that users have very little space in the home directory.  If
you choose to host your repository there, you will have to move your build
and work folders to other file systems. That is why the present system uses
the environment variable storage_model, if it exists, to store your build
and work directories. In that case, you can use the commands listed above in
the Instructions in a nutshell section to compile and install GEM:
```
make cmake
make build
make work
```

If you don't want to use that variable, here is an example of another way:
```
cd $HOME
git clone git@gitlab.science.gc.ca:cpi001/gem-devel.git
cd gem-devel
. ./eccc_setup_intel_19
# We assume that you have a link to your ords folder in your home
mkdir -p $HOME/ords/gem/{build,work}
mv build-ubuntu-18.04-amd64-64-intel-19.0.3.199 $HOME/ords/gem/build
mv work-ubuntu-18.04-amd64-64-intel-19.0.3.199 $HOME/ords/gem/work
rmdir build work
ln -s $HOME/ords/gem/build .
ln -s $HOME/ords/gem/work .
cd build/build-ubuntu-18.04-amd64-64-intel-19.0.3.199
cmake $HOME/gem/project
make -j
make work
```

# Using Visual Studio Code

Please refer to the instructions on the
[Visual Studio Code Wiki page](https://wiki.cmc.ec.gc.ca/wiki/Visual_Studio_Code)
to load the appropriate SSM domains.  It can then be launched by executing the
**vstudio** script.

You can open a terminal within VSCode and go to the working directory to run GEM.

# More about compilers

## PGI 18.4 on einstein

Compiles when debugging is not enabled, but the resulting executables don't work.
In progress on newer versions of PGI (20.11).

## xlf 16 on p9gpu-1 and p9gpu-2

TODO!
