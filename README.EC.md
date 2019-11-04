The following instructions explain how to use GEM on a stick at the CMC.

# Confguring the environment

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
**.eccc_setup_intel_16**.  It will define environement variables and load the
appropriate SSM packages and will create a build directory with the compiler
and version in it's name.
```
   . .eccc_setup_intel_16
```
{: .language-bash}


To configure your environment for the gfortran 5.1 compiler, source
**.eccc_setup_gfortran_5.1**.  It will define environement variables and load the
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
- Newer worstations on the cmc.ec.gc.ca domain

To configure your environment for the Intel 19 compiler, source
**.eccc_setup_intel_19**.  It will define environement variables and load the
appropriate SSM packages and will create a build directory with the compiler
and version in it's name.
```
   . .eccc_setup_intel_19
```

Since gfortran 7.4 is available through the system packages, there should be
no need to load any SSM package to use this compiler.



# Compiling and executing the model

Once the environement has been configured, please refer to the main
[README](README.md) for compiling instructions.


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
