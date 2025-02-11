# GEM Benchmarking

Welcome to the GEM Benchmarking System! 

You should have obtained this benchmark from https://github.com/ECCC-ASTD-MRD/gem

# Requirements

* Fortran and C compiler. Theses codes have been tested with compilers from GNU and Intel OneAPI (classic and llvm based)
* An MPI implementation such as OpenMPI, MPICH or Intel MPI (with development package)
* OpenMP support
* BLAS, LAPACK or equivalent mathematical/scientific library (ie: MKL), with development package
* fftw3 library (with development package) (tested with 3.3.10)
* CMake (version >= 3.20)

# Build and run steps

## Compiler specifics

* Compiler specific definitions and flags are defined within the ```cmake_rpn``` submodule of each code repository. If you need to change or add any, 
you can add or modify the rules into `[git source path]/cmake_rpn/modules/ec_compiler_presets/default/[architecture]/`

## Build base library (librmn)

```bash
git clone git@github.com:ECCC-ASTD-MRD/librmn.git
cd librmn
git checkout alpha
git submodule update --init --recursive
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=[rmn install directory] ..
make install
```

## Build verification tool (sverif)

```bash
git clone git@github.com:ECCC-ASTD-MRD/sverif.git
cd sverif
mkdir build
cd build
cmake -Drmn_ROOT=[rmn install directory] -DCMAKE_INSTALL_PREFIX=[sverif install directory] ..
make install
export PATH=[sverif install directory]/bin:$PATH
```

## Build model (GEM)

```bash
git clone git@github.com:ECCC-ASTD-MRD/gem.git
cd gem
git checkout benchmark-5.3
git submodule update --init --recursive
./download-dbase.sh .
./download-dbase-benchmarks.sh .
. ./.common_setup [intel|gnu|nvhpc]
mkdir build
cd build
cmake ..
make -j 5
make work
```

## Run model (GEM)

```bash
cd ../$GEM_WORK
```

* This will give you the possible CPU decomposition for the GEM_cfgs_GY_4km
model configuration:

```
findtopo -npex_low 20 -npex_high 250 -npey_low 20 -npey_high 200 -corespernode 80 -nml $GEM_WORK/configurations/GEM_cfgs_GY_4km/cfg_0000/gem_settings.nml > topo.txt
```

## Run preparation script

```
runprep.sh -dircfg configurations/GEM_cfgs_GY_4km
```

* Run model (or submit to queing system):

The -cpus parameters defines the mpi topology and the openmp number of
threads [X MPI]x[Y MPI]x[OMP threads] (ie: 61x20x8).

This topology is used for 2 simultaneous model run (Yin+Yang), meaning you
will have to use double this number of CPU for the submision/run itself.

In the below example, 9760 cores are needed per sub run, so a total of 19520 cores will be needed.

```
runmod.sh -dircfg configurations/GEM_cfgs_GY_4km -ptopo 61x20x8
```

* If you need a job script so you can submit the job to a queuing system, make sure you define the GOAS_SCRIPT 
variable and load the common setup before running the runmod.sh command:

```bash
cd [gem src path]
. ./.common_setup [intel|gnu|nvhpc]
cd $GEM_WORK
runmod.sh -dircfg configurations/GEM_cfgs_GY_4km -ptopo 61x20x8
```

## Run verification (sverif)

* This script will provide a PASS or FAIL rating

```bash
cd ..
gem_sverif.sh -p $GEM_WORK -f dp2020022915-000-000_006
```

* Expected output

```bash
(INFO) Run configuration found (GEM_cfgs_GY_4km-5.3.0-a11)
sverif_eval.Abs GZ 500 006 GEM/work-rhel-8-icelake-64-intel-2021.5.0//RUNMOD/output/cfg_0000/laststep_0000000180/000-000/dp2020022915-000-000_006 GEM/work-rhel-8-icelake-64-intel-2021.5.0/sverif/GEM_cfgs_GY_4km-5.3.0-a11
PASS T1  (sverif_eval) GZ [ 500mb;    6h; CI=0.01]
PASS NT1 (sverif_eval) GZ [ 500mb;    6h; CI=0.01]
PASS NT5 (sverif_eval) GZ [ 500mb;    6h; CI=0.01]
PASS R   (sverif_eval) GZ [ 500mb;    6h; CI=0.01]
PASS T1  (sverif_eval) TT [ 850mb;    6h; CI=0.01]
PASS NT1 (sverif_eval) TT [ 850mb;    6h; CI=0.01]
PASS NT5 (sverif_eval) TT [ 850mb;    6h; CI=0.01]
PASS R   (sverif_eval) TT [ 850mb;    6h; CI=0.01]
(INFO) Passed (passed 8/8)
```

# Detailed information on Packages

## Building GEM and launching on your platform

* [Building GEM](README-gem.md#outside-cmc-external-users)
* [Running GEM](README-gem.md#running-gem)
 
## Result validation

* [Building sverif](https://github.com/ECCC-ASTD-MRD/sverif/blob/main/README.md)
