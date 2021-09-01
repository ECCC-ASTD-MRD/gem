# Introduction
This document describes the process of building and running GEM with a global Yin-Yang 15km grid pressure coordinate (GY15)

# Requirements
* Fortran and C compiler
* MPI
* fftw3 
* libnuma
* May need to export other variables for MPI to function properly

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
./download_dbase.sh .
```
* GEM GY15 input database
```
./download_benchmark_gy15.sh
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

* Execute preparation job interactively
```
runprep.sh -dircfg configurations/GEM_cfgs_GY15
```
* or modify myprep.pbs and submit in batch
```
qsub myprep.pbs
```

* Once runprep is done, modify mybatch.pbs (use the guide in gem/configurations/GEM_cfgs_GY15_P/TOPOs_possible_for_GY15 to choose PE topology) and run the model
```
qsub mybatch.pbs
```
