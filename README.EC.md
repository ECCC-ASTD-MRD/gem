The following instructions explain how to use GEM on a stick at the CMC.

# Configuring the environment

If you do not know the version of a system on which you are connected,
please run the following command: `lsb_release -a`

## Ubuntu 14.04

The Intel 16 compiler and gfortran 5.1 are available on Ubuntu 14.04 systems.
Here is an incomplete list of such systems:

- ppp1.science.gc.ca
- ppp2.science.gc.ca
- hare.science.gc.ca
- brooks.science.gc.ca
- Older worstations on the cmc.ec.gc.ca domain

To configure your environment for the Intel 16 compiler, source
**.eccc_setup_intel_16**.  It will define environment variables and load the
appropriate SSM packages and will create a build directory with the compiler
and version in it's name.
```
   . .eccc_setup_intel_16
```


To configure your environment for the gfortran 5.1 compiler, source
**.eccc_setup_gfortran_5.1**.  It will define environment variables and load the
appropriate SSM packages and will create a build directory with the compiler
and version in it's name.
```
   . .eccc_setup_gfortran_5.1
```


## Ubuntu 18.04

The Intel 19 compiler and gfortran 7.4 are available on Ubuntu 18.04 systems.
Here is an incomplete list of such systems:

- ppp3.science.gc.ca
- ppp4.science.gc.ca
- daley.science.gc.ca
- banting.science.gc.ca
- rutherford.cmc.ec.gc.ca
- Newer workstations on the cmc.ec.gc.ca domain

To configure your environment for the Intel 19 compiler, from the root of your
GOAS clone, source **.eccc_setup_intel_19**.  It will define environment
variables and load the appropriate SSM packages and will create a build
directory with the compiler and version in it's name.
```
   . .eccc_setup_intel_19
```

Since gfortran 7.4 is available through the system packages, there should be
no need to load any SSM package to use this compiler.



# Compiling and executing the model

Once the environment has been configured, there are a few extra steps to do
depending on how you where you want to host you clone of the repository.

## Repository in a data file system

The simplest way to work with GOAS is to have it in a data file system such as
ORDS or a local file system on your workstation.  This will lead to a usage
pattern that resembles the one of users outside the CMC.  The most important
downside to this approach is that there are no backups or snapshots of the code
you have modified.  Also, if you choose a local file system, it might not be
accessible from all hosts.

To ensure that your modifications are backed-up, you can make sure to regularly
push them to a remote Git repository or to implement you own backup scheme.

When sourced, the `.eccc_setup_*` will create a build directory following the
pattern below:
`<goas_root>/build/<compiler>-<version>`

You can therefore go to that folder and simply execute the commands bellow to
compile the code:
```
cmake ../../project
make -j
make install
```

## Repository in the $HOME directory

Hosting your GOAS working copy in your home has a few advantages, but it is also
slightly more complex.  Being in the home directory, your repository is
automatically backed-up and you also have multiple snapshots.  The downside is
that users have very little space in the home directory.  If you choose to host
your repository there, you will have to move your build and work folders to
other file systems.

For example:
```
cd $HOME
git clone git@gitlab.com:eccc/gem/gem.git
cd gem
git checkout dev
. .eccc_setup_intel_19
# We assume that you have a link to your ords folder in your home
mkdir -p $HOME/ords/gem/{build,work}
mv build/* $HOME/ords/gem/build
rmdir build
ln -s $HOME/ords/gem/build
ln -s $HOME/ords/gem/work
cd build/intel-19.0.3.199
cmake $HOME/gem/project
make -j
make install
```

Further information about compiling and executing can be found in the main
[README](README.md).


# Using Visual Studio Code

Please refer to the instructions on the
[Visual Studio Code Wiki page](https://wiki.cmc.ec.gc.ca/wiki/Visual_Studio_Code)
to load the appropriate SSM domains.  It can then be launched by executing the
**vstudio** script.

You can open a terminal within VSCode and go to the working directory to run GEM.



# More about compilers

## PGI 17.10 on lorentz

The compiler and libraries can be used by executing the following commands:
```
   PGI=/local/raid/armn/CUDA/PGI/pgi; export PGI
   MANPATH=$MANPATH:$PGI/linux86-64/17.10/man; export MANPATH
   LM_LICENSE_FILE=$PGI/license.dat; export LM_LICENSE_FILE
   PATH=$PGI/linux86-64/17.10/bin:$PATH; export PATH
   export LD_LIBRARY_PATH=$PGI/linux86-64/17.10/lib:$LD_LIBRARY_PATH
   export PATH=$PGI/linux86-64/17.10/mpi/openmpi/bin:$PATH
   export MANPATH=$MANPATH:$PGI/linux86-64/17.10/mpi/openmpi/man
```


## PGI 18.4 on einstein

Compiles when debugging is not enabled, but the resulting executables will not work.
```
   cd gem_on_a_stick_[version]/sources
   rm Linux-x86_64-pgi.cmake; ln -s Linux-x86_64-pgi.cmake_18.4 Linux-x86_64-pgi.cmake
   export MODULEPATH=$MODULEPATH:/opt/pgi/modulefiles
   module load PrgEnv-pgi/18.4
```

## xlf 16 on p9gpu-1 and p9gpu-2

TODO!
