# Introduction
This document describes the process of building and running GEM with a global Yin-Yang 15km grid pressure coordinate (GY15)

# Requirements
* Fortran and C compilers
* An MPI implementation such as OpenMPI (with development package),
* Portable Hardware Locality library (hwloc) (with development package),
* OpenMP support (optional)
* BLAS library (with development package),
* LAPACK library (with development package),
* fftw3 library (with development package),
* basic Unix utilities such as cmake (version 2.8.7 minimum), bash, sed, etc.
* on Cray XC, libnuma

   ## Example on Cray XC50
   ```
   export OMP_STACKSIZE=4G
   export FOR_DISABLE_KMP_MALLOC=0
   module unload intel/19.0.3.199
   module load intel/19.0.5.281
   module add craype-hugepages16M
   module add cray-fftw/3.3.8.2
   # on CMC XC50
   . ssmuse-sh -x main/opt/numa/numa-2.0.14
   ```

# Getting the Code
```
git clone https://github.com/ECCC-ASTD-MRD/gem gem
cd gem
git checkout benchmark-5.2.0-a7
git submodule update --init --remote
```

# Getting the data
* GEM input database
```
./download-dbase.sh .
```
* GEM GY15 input database
```
./download-benchmark-GY15.sh .
```

# Building GEM (with intel compiler)
```
. ./.common_setup intel
mkdir build
cd build
cmake ..

# you may have to specify the compiler as follows:
cmake -DCOMPILER_SUITE=intel ..

make -j
make work
cd ..
```

# Running GEM

You can find information in gem/configurations/GEM_cfgs_GY15_P:
1) file TOPOs_possible_for_GY15 (possible PE topologies)
2) sample pbs file: myprep.pbs (job to prepare for model run)
3) sample pbs file: mybatch.pbs (job to run model)
 
```
cd work-*
cp ../configurations/GEM_cfgs_GY15_P/*pbs .
```

* If not already done, load the needed environment got GEM scripts
```
. ./.common_setup intel
```

* Execute preparation job interactively
```
runprep.sh -dircfg configurations/GEM_cfgs_GY15_P
```
* or modify myprep.pbs and submit in batch
```
qsub myprep.pbs
``` 

* Once runprep is done, modify mybatch.pbs (use the guide in gem/configurations/GEM_cfgs_GY15_P/TOPOs_possible_for_GY15 to choose PE topology) and run the model
```
qsub mybatch.pbs
```

# Results

* Results will be written into the directory **work-*/RUNMOD/output/cfg_0000/**
* The execution listing, **list_gy**, will be located in the directory **work-*/**
