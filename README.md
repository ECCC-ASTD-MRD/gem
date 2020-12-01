# Instructions in a nutshell

See below for extended instructions.  Further details are can be found in
[GEM's manual](doc/GEM-manual.pdf) (PDF).

```
    git clone https://gitlab.com/eccc/gem/gem.git
    cd gem
    # Optionnaly checkout a specific branch:
    # git checkout <branch_name>

    # Create a build directory
    # The build directory has very little value in itself and can be placed
    # outside the project directory
    mkdir -p build
    cd build

    # If the -DWORK_PREFIX=<path> option isn't given to cmake, the work directory
    # will be created under the work directory
    cmake ../project
    # Create an execution environment for GEM
    make -j work

    # Download the data files required to run GEM
    cd ..
    ./download-dbase.sh work

    cd ../work/work-${OS_NAME}-${COMPILER_NAME}
    ./runprep -dircfg configurations/GY_cfgs
    ./runmod -dircfg configurations/GY_cfgs

    ./voir -iment RUNMOD.dir/output/cfg0000/ ...
    ./fststat -fst RUNMOD.dir/output/cfg0000/ ...
```

[SPI](http://eer.cmc.ec.gc.ca/software/SPI) can be used to view the output files.
2D fields can also be displayed with [xrec](https://gitlab.com/gem-ec/xoas)


-----------------------------------------------------------------
# Extended instructions:

## Requirements

To compile and run GEM, you will need:
- Fortran and C compilers, (If using GNU compilers, only versions up to 9.2 will work)
- An MPI implementation such as OpenMPI (with development package),
- OpenMP support (optional)
- BLAS library (with development package),
- LAPACK library (with development package),
- fftw3 library (with development package),
- basic Unix utilities such as cmake (version 2.8.7 minimum), bash, sed, etc.

## Data for examples
After having cloned or download the git tar file of GEM from
[GitLab.com](https://gitlab.com/eccc/gem/gem), execute the script named
**download-dbase.sh** or download and untar the data archive with the following
link:
[http://collaboration.cmc.ec.gc.ca/science/outgoing/goas/gem_dbase.tar.gz](http://collaboration.cmc.ec.gc.ca/science/outgoing/goas/gem_dbase.tar.gz)

## Compiler specifics

### GNU compilers
- By default GEM is configured to use gfortran and gcc compilers
- make sure the compilers and libraries paths are set in the appropriate
  environment variables (PATH and LD_LIBRARY_PATH).  Here are some examples
  which will need to be adapted for your setup:
    - On Ubuntu:
        ```
            export PATH=/usr/lib/openmpi/bin:${PATH}
            export LD_LIBRARY_PATH=/usr/lib/openmpi/lib:$LD_LIBRARY_PATH
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

The default compilers is GNU. If you want to compile with other compilers,
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

The installation process will create directory named after the operating system
on which the compilation was executed, and the compiler you used
(work-${OS_NAME}-${COMPILER_NAME}). For example
*work-Fedora-29-x86_64-gfortran-8.2.1* could be created in the *work*
directory, and the following executables are installed in the *bin* sub-folder:
- maingemdm
- gemgrid
- yy2global
- flipit
- r.fstinfo
- voir
- editfst
- fststat
- cclargs_lite
- pgsm

## Running GEM

Go to the working directory, named *work-${OS-NAME}_${COMPILER-NAME}*, for
example *work-Linux_x86-64-gfortran-8.2.1*

```
    cd <WORK_PREFIX>/work-${OS-NAME}_${COMPILER-NAME}
    # Configure the model for a specific case
    ./runprep -dircfg configurations/${cfgDir}
    # Create a directory named RUNMOD.dir for output files
    # and executes the model
    ./runmod -dircfg configurations/${cfgDir}
```

*runmod*'s ```-ptopo``` argument can be used to specify the number of CPU to
use.  For example,  ```-ptopo 2x2x1``` will use 4 cpus for a LAM, and
8 cpus for global Yin-Yang - see the [manual](doc/GEM-manual.pdf) for details)

**mpirun** is used to run the model.  You can specify another command
by editing *scripts/Um_model.sh either in the original scripts in the
*scripts* folder or in the working directory
*work/work-${OS_NAME}_${COMPILER_NAME}/scripts*. In the latter case, be advised
that any changes done to the scripts will be over-written the next time
```make install``` is executed. Therefore, if you want your changes to the
scripts to be persistent, make your modifications in the *scripts* directory,
and then rerun the ```make install``` in the build directory.

## Working with model outputs

The model stores it's outputs in FST files.  These tools can be used to perform
various tasks on the output files:

- ```voir``` lists the records in FST files:
    ```
        ./voir -iment RUNMOD.dir/output/cfg0000/laststep0000000024/000-000/dm2009042700-000-000010
    ```

- ```fststat``` produces statistical means of the records in a FST file:
    ```
        ./fststat -fst RUNMOD.dir/output/cfg0000/laststep0000000024/000-000/dm2009042700-000-000010
    ```

- ```pgsm``` can be used to interpolate records to a different grid
    ```
        ./pgsm -iment <input FST> -ozsrt <output FST> -i <pgsm.directives>
    ```

- ```editfst``` enables basic record manipulations, such as copying to another
    file.
    ```
        ./editfst -s <input FST> -d <output FST> -i <editfst.directives>
    ```

[SPI](http://eer.cmc.ec.gc.ca/software/SPI) is a scientific and meteorological
virtual globe offering processing, analysis and visualization capabilities,
with a user interface similar to Google Earth and NASA World Wind, developed by
Environment Canada.

[xrec](https://gitlab.com/gem-ec/xoas) is another visualization program which
can be used to display 2D meteorological fields stored in the FST files,
developed by Research Informatics Services, Meteorological Research Division,
Environment and Climate Change Canada.


# Using Visual Studio Code

If [Visual Studio Code](https://code.visualstudio.com/) is available on your
system, the `goas-devel.code-workspace` file can be opened by code to
automatically open relevant files.
