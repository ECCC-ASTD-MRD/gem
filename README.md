# Introduction for benchmark-5.2.0-a13

This document describes the process of building and running 3 GEM cases:

* global benchmark Yin-Yang 15km grid pressure coordinate (GY15_P)
* global benchmark Yin-Yang 5km grid pressure coordinate (GY5_P)
* global test case Yin_Yang (GEM_cfgs)

# Requirements
* Fortran and C compilers
* An MPI implementation such as OpenMPI (with development package),
* OpenMP support (optional)
* (BLAS,LAPACK) or equivalent mathematical/scientific library (ie: MKL)
* fftw3 library (with development package),
* basic Unix utilities such as cmake (version 2.8.7 minimum), bash, sed, etc.

   ## Example on Cray XC50
   
   ```
   export OMP_STACKSIZE=4G
   export FOR_DISABLE_KMP_MALLOC=0
   module unload intel/19.0.3.199
   module load intel/19.0.5.281
   module add craype-hugepages16M
   module add cray-fftw/3.3.8.2
   ```

# Getting the Code
```
git clone https://github.com/ECCC-ASTD-MRD/gem gem
cd gem
git checkout benchmark-5.2.0-a13
git submodule update --init --remote
```

# Getting the input data

* gem_dbase (basic GEM data, required for all cases)

```
./download-dbase.sh .
```
* GY15kmP_dbase (data required for running GY15_P)

```
./download-benchmark-GY15p.sh .
```
* GY5kmP_dbase (data required for running GY5_P)

```
./download-benchmark-GY5p.sh .
```

# Loading the requirements necessary as mentioned above (compiler,libraries,software)

**Important**: in order to export the GEM environment variables, you must
  source the setup script called .common_setup at this level:

```
. ./.common_setup gnu

or

. ./.common_setup intel
```

# Creating a build directory

The build directory has very little value in itself and can be placed
outside the source directory:

```
mkdir build
cd build
```

# Building GEM with cmake
```
    cmake ..
   # or
    cmake -DCOMPILER_SUITE=intel ..
   # or
    cmake -DCOMPILER_SUITE=gnu ..
   # NOTE:on some platforms(like XCs), do not use -DCOMPILER_SUITE, just use:
    cmake ..
   #Compile and create the binaries and the work directory for executing GEM
    make -j work
    cd ..
```

# Running GEM

* You should now have a work directory created in the form of :
work-${OS_NAME}-${COMPILER_NAME}
* runprep.sh prepares the input files in work-${OS_NAME}-${COMPILER_NAME}/PREP
* runmod.sh runs the model in work-${OS_NAME}-${COMPILER_NAME}/RUNMOD

Now test your environment and your executables by running this test case:

```
    cd work-${OS_NAME}-${COMPILER_NAME}
    runprep.sh -dircfg configurations/GEM_cfgs
    runmod.sh -dircfg configurations/GEM_cfgs -ptopo 1x1x1
    runmod.sh -dircfg configurations/GEM_cfgs -ptopo 2x2x2

```

* Note: -ptopo [Npex]x[Npey]x[OMP_threads] for each grid in Yin-Yang.
* And other possible PE topologies are listed in gem/configurations/GEM_cfgs/TOPOLOGIES_possible

* IF the small GEM_cfgs works, then tryout the larger configurations
* Possible PE topologies are listed for each configuration
* A sample PBS job is given to submit to the backend

* For the mini GEM test: look in gem/configurations/GEM_cfgs:
1) file TOPOLOGIES_possible
2) sample pbs file: myjobtest.pbs

* For Global Yin-Yang 15km: look in gem/configurations/GY15km_P:
1) file TOPOs_possible_for_GY15 (possible PE topologies)
2) sample pbs file: myjobbench15.pbs

* For Global Yin-Yang 5km: look in gem/configurations/GY5km_P:
1) file TOPOLOGIES_iterative_solver (possible PE topologies for iterative solver)
2) files TOPOLOGIES_direct_solver (possible PE topologies for direct solver)
3) sample pbs file: myjobbench5.pbs
 
* For each model run on the back end, one must load the compiler,libraries
and activate the GEM environment at the level where GEM was cloned 
```
. ./.common_setup intel
```

* Then go to the working directory:
```
cd work-${OS_NAME}-${COMPILER_NAME}

cp ../configurations/GEM_cfgs/*pbs .
cp ../configurations/GY15km_P/*pbs .
cp ../configurations/GY5km_P/*pbs .
```

* For GEM_cfgs, use the guide in gem/configurations/GEM_cfgs/TOPOLOGIES_possible
to choose a PE topology, modify myjobtest.pbs accordingly and submit job:
```
qsub myjobtest.pbs
```

* For GY15km_P, use the guide in gem/configurations/GY15km_P/TOPOs_possible_for_GY15 to choose a PE topology, modify myjobbench15.pbs accordingly and submit job:
```
qsub myjobbench15.pbs
```

* For GY5km_P, use the guide in gem/configurations/GY5km_P/TOPOLOGIES_iterative_solver to choose a PE topology, modify myjobbench15.pbs accordingly and submit job:
```
qsub myjobbench5.pbs
```

* NOTE: iterative solver is faster for greater number of nodes and allows more choices for the PE topologies.

# Results

* Model outputs will be written into the directory **work-*/RUNMOD/output/cfg_0000/**
* The execution listings will be located in the directory **work-*/**

# Troubleshooting

- If needed, changes to the C and Fortran flags can be done either in intel.cmake or
  gnu.cmake file (according to the compilers you are using) in **cmake_rpn/ec_compiler_presets/default/Linux-x86_64** folder.
  
- If you get error messages (for example, compiler or MPI/OpenMPI not
  found), make sure that the ```PATH``` and ```LD_LIBRARY_PATH``` environment
  variables contain the appropriate paths.  You can also add
  ```-DCMAKE_VERBOSE_MAKEFILE=ON``` to you **cmake** command line to generate
  verbose makefiles which will print the exact compiler command lines issued.

- If the compiler or compile options are not right:
    - Remove the content of the build directory
    - Make appropriate changes to the cmake files corresponding to the
      compilers you are using
    - Re-launch the commands starting at cmake
