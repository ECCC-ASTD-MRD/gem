# Instructions in a nutshell

**The 5.2 branch is a development version.
The 5.1 branch is the stable version, used in production at the Canadian
Meteorological Centre.
Benchmarks must be used with 5.2 version and benchmarks branch.**

See below for extended instructions.  Further details are can be found in
[GEM's manual](doc/GEM-manual.pdf) (PDF).

```
    git clone https://github.com/ECCC-ASTD-MRD/gem.git
    cd gem
    # Optionnaly checkout a specific branch:
    # git checkout <branch_name>
 
    # Download the data files required to run GEM
    ./download-dbase.sh .

    # Important: in order to set environment variables needed to run GEM, use the
    # following setup file, after setting up your compiler environment:
    . ./.common_setup gnu
    or
    . ./.common_setup intel
	
    # Create a build directory
    # The build directory has very little value in itself and can be placed
    # outside the project directory
    mkdir -p build
    cd build

    # If the -DWORK_PREFIX=<path> option isn't given to cmake, the work directory
    # will be created with the name work-${OS_NAME}-${COMPILER_NAME}
    cmake ../project
    # Create an execution environment for GEM
    make -j work

    cd ..
	cd work-${OS_NAME}-${COMPILER_NAME}
    # Configure the model with one of the configurations
    # and execute the model, for example:
    runprep.sh -dircfg configurations/GEM_cfgs_GY_FISL_P
    runmod.sh -dircfg configurations/GEM_cfgs_GY_FISL_P

	# use tools to see list the records in the output files
    voir -iment RUNMOD/output/cfg0000/ ...
    fststat -fst RUNMOD/output/cfg0000/ ...
```

[SPI](http://eer.cmc.ec.gc.ca/software/SPI) can be used to view the output files.
2D fields can also be displayed with [xrec](https://gitlab.com/gem-ec/xoas)

For benchmarks, please see the [benchmarks branch](https://github.com/ECCC-ASTD-MRD/gem/tree/benchmarks)

-----------------------------------------------------------------
# Extended instructions:

## Requirements

To compile and run GEM, you will need:
- Fortran and C compilers
- An MPI implementation such as OpenMPI (with development package),
- Portable Hardware Locality library (hwloc) (with development package),
- OpenMP support (optional)
- BLAS library (with development package),
- LAPACK library (with development package),
- fftw3 library (with development package),
- basic Unix utilities such as cmake (version 2.8.7 minimum), bash, sed, etc.

## Data for examples
After having cloned or downloaded the git tar file of GEM from
[github.com](https://github.com/ECCC-ASTD-MRD/gem), execute the script named
**download-dbase.sh** or download and untar the data archive with the following
link:
[http://collaboration.cmc.ec.gc.ca/science/outgoing/goas/gem_dbase.tar.gz](http://collaboration.cmc.ec.gc.ca/science/outgoing/goas/gem_dbase.tar.gz)

## Compiler specifics

### GNU compiler suite
- By default GEM is configured to use gfortran and gcc compilers
- Changes to the C and Fortran flags can be done in **project/Linux-x86_64-gnu.cmake**
- Make sure the compilers and libraries paths are set in the appropriate
  environment variables (PATH and LD_LIBRARY_PATH).  Here are some examples
  which will need to be adapted for your setup:
    - On Ubuntu:
        ```
          export PATH=/usr/lib/openmpi/bin:${PATH}
          export LD_LIBRARY_PATH=/usr/lib/openmpi/lib:$LD_LIBRARY_PATH
        or
          export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/openmpi/lib:$LD_LIBRARY_PATH

        ```
    - On Fedora:
        ```
          export PATH=/usr/lib64/openmpi/bin:$PATH
          export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH
        ```

### Intel Compilers

- Changes to the C and Fortran flags can be done in **project/Linux-x86_64-intel.cmake**
    - You may need to modify ```-march``` to generate code that can run on
      your system
- Make sure the compilers and libraries are in the appropriate
  environment variables (PATH and LD_LIBRARY_PATH)


## Compiling and installing GEM

You can add extra CMake arguments such as```-DCMAKE_VERBOSE_MAKEFILE=ON```.
You can also add ```-j``` to **make** to launch multiple compile tasks in
parallel.

[OpenMP](https://www.openmp.org/) is enabled by default.  If you wish to build
without OpenMP support, you can add the ```-DWITH_OPENMP=OFF``` argument when
running **cmake**.

The default compiler suite is GNU. If you want to compile with other compilers,
add ```-DCOMPILER=<compiler suite name (gnu|intel)>``` to the CMake
command line.  This release has been tested with GNU and Intel compilers on
Linux x86_64.  Other compilers have also been used in the past but have not been
tested with the current release.  You will likely have to modify the *.cmake
files in the **project** folder.

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

The installation process will create a directory named after the operating system
on which the compilation was executed, and the compiler you used
(work-${OS_NAME}-${COMPILER_NAME}). For example
*work-Fedora-32-x86_64-gnu-10.2.1* could be created in the main directory,
and the following executables are installed in the *bin* sub-folder: 
- cclargs_lite
- checkdmpart
- editfst
- feseri
- flipit
- fstcomp
- fststat
- gem_monitor_end
- gem_monitor_output
- gemgrid
- maingemdm
- pgsm
- prgemnml
- prphynml
- r.fstinfo
- voir
- yy2global


## Running GEM

Go to the working directory, named *work-${OS-NAME}_${COMPILER-NAME}*, for
example *work-Fedora-32-x86_64-gnu-10.2.1*

```
    cd work-${OS-NAME}_${COMPILER-NAME}
    # Configure and execute the model for a specific case, for example:
    runprep.sh -dircfg configurations/GEM_cfgs_GY_FISL_P
    runmod.sh -dircfg configurations/GEM_cfgs_GY_FISL_P
```

*runmod.sh*'s ```-ptopo``` argument can be used to specify the number of CPU to
use.  For example,  ```-ptopo 2x2x1``` will use 4 cpus for a LAM, and
8 cpus for global Yin-Yang - see the [manual](doc/GEM-manual.pdf) for details)

If you get an error message saying runprep.sh or gem_dbase is not found,
make sure to set the environment variables using the setup file situated in
the main directory: 
*./.common_setup gnu* or *./.common_setup intel*

An in-house script (**r.run_in_parallel**) is used to run the model. If you
want to use another command, or if it doesn't work in your environment, edit
the file *scripts/gem_mpirun.sh* to change the script, or move the file
*scripts/r.run_in_parallel* so that the model is run with mpirun directly.


## Working with model outputs

The model stores its outputs in FST files.  The following tools can be used to perform
various tasks on the output files:

- ```voir``` lists the records in FST files:
    ```
        voir -iment RUNMOD.dir/output/cfg0000/laststep0000000024/000-000/dm2009042700-000-000010
    ```

- ```fststat``` produces statistical means of the records in a FST file:
    ```
        fststat -fst RUNMOD.dir/output/cfg0000/laststep0000000024/000-000/dm2009042700-000-000010
    ```

- ```pgsm``` can be used to interpolate records to a different grid
    ```
        pgsm -iment <input FST> -ozsrt <output FST> -i <pgsm.directives>
    ```

- ```editfst``` enables basic record manipulations, such as copying to another
    file.
    ```
        editfst -s <input FST> -d <output FST> -i <editfst.directives>
    ```

[SPI](http://eer.cmc.ec.gc.ca/software/SPI) is a scientific and meteorological
virtual globe offering processing, analysis and visualization capabilities,
with a user interface similar to Google Earth and NASA World Wind, developed by
Environment Canada.

[xrec](https://gitlab.com/gem-ec/xoas) is another visualization program which
can be used to display 2D meteorological fields stored in the FST files,
developed by Research Informatics Services, Meteorological Research Division,
Environment and Climate Change Canada.
