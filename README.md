# Instructions in a nutshell

**The 5.1 branch is the stable version, used in production at the Canadian
Meteorological Centre.**

See below for extended instructions and requirements. 

```
    git clone https://github.com/ECCC-ASTD-MRD/gem.git
    cd gem
    # Checkout 5.1  branch:
    git checkout 5.1-branch

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
    # replace "work" by the WORK_PREFIX path you may have used above
    cd ..
    ./download-dbase.sh work

    cd work/work-${OS_NAME}-${COMPILER_NAME}
    # Configure the model with the default configuration (configurations/GEM_cfgs)
    # and execute the model
    ./runprep
    ./runmod
    # or configure the model for a specific case, for example:
    ./runprep -dircfg configurations/GY_cfgs
    ./runmod -dircfg configurations/GY_cfgs

    # use tools to see list the records in the output files
    ./voir -iment RUNMOD.dir/output/cfg0000/ ...
    ./fststat -fst RUNMOD.dir/output/cfg0000/ ...
```

[SPI](https://github.com/ECCC-ASTD-MRD/SPI) can be used to view the output
files.

2D fields can also be displayed with
[xrec](https://github.com/ECCC-ASTD-MRD/xrec/tree/feat_cmake) (make sure you
use feat_cmake branch)

-----------------------------------------------------------------
# Extended instructions:

## Requirements

To compile and run GEM, you will need:
- Fortran and C compilers
- An MPI implementation such as OpenMPI (with development package),
- OpenMP support (optional)
- (BLAS,LAPACK) or equivalent mathematical/scientific library (ie: MKL), with development package,
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
- By default GEM is configured to use gfortran and gcc compilers, and OpenMPI
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
  environment variables (```PATH``` and ```LD_LIBRARY_PATH```)

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
*work-Fedora-32-x86_64-gnu-10.2.1* could be created in the *work*
directory, and the following executables are installed in the *bin* sub-folder:
- maingemdm
- maingemgrid
- mainyy2global
- flipit
- r.fstinfo
- voir
- editfst
- fststat
- cclargs_lite
- pgsm

## Running GEM

Go to the working directory, named *work-${OS-NAME}_${COMPILER-NAME}*, for
example *work-Fedora-32-x86_64-gnu-10.2.1*

```
    cd <WORK_PREFIX>/work-${OS-NAME}_${COMPILER-NAME}
    # Configure the model with the default configuration (configurations/GEM_cfgs)
    ./runprep
    # Create a directory named RUNMOD.dir for output files
    # and execute the model
    ./runmod
    # or configure and execute the model for a specific case, for example:
    ./runprep -dircfg configurations/GY_cfgs
    ./runmod -dircfg configurations/GY_cfgs
	make sure you give the complete path to the directory in which the
	cfg_000 directory is located, such as:
    ./runprep -dircfg configurations/GEM_theo_cfgs/WILLIAMSON/CASE1_AROUND
    ./runmod -dircfg configurations/GEM_theo_cfgs/WILLIAMSON/CASE1_AROUND
```

*runmod*'s ```-ptopo``` argument can be used to specify the number of CPU to
use.  For example,  ```-ptopo 2x2x1``` will use 4 cpus for a LAM, and
8 cpus for global Yin-Yang.

**mpirun** is used to run the model.  You can specify another command
by editing *scripts/Um_model.sh either in the original scripts in the
*scripts* folder or in the working directory
*work/work-${OS_NAME}_${COMPILER_NAME}/scripts*. In the latter case, be advised
that any changes done to the scripts will be over-written the next time
```make work``` is executed. Therefore, if you want your changes to the
scripts to be persistent, make your modifications in the *scripts* directory,
and then rerun the ```make work``` in the build directory.

See **README.run** in the working directory for other information on the different configurations.

## Working with model outputs

The model stores its outputs in FST files.  The following tools can be used to perform
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

[xrec](https://github.com/ECCC-ASTD-MRD/xrec/tree/feat_cmake) is another visualization program which
can be used to display 2D meteorological fields stored in the FST files,
developed by Research Informatics Services, Meteorological Research Division,
Environment and Climate Change Canada.


## Configurations files

The execution of all three components of GEM is configurable through the use
of three configuration files called:
- gem_settings.nml: file containing some namelists to configure the model execution
- outcfg.out: file use to configure the model output
- configexp.cfg: file use to configure the execution shell environment

Examples of these files can be found in the test cases in the configurations
directory.

A fourth configuration file, named physics_input_table, is used for GEM_cfgs
test cases.

## Running your own configuration

Put the three configurations files (gem_settings.nml, outcfg.out and
configexp.cfg) in a directory structure such as: **experience/cfg_0000** in
the configurations directory.

The master directory name (**experience** in the example above) can be
any valid directory name. However, the second directory must have the name
\textit{cfg\_0000}.

Then run the two scripts with the following commands, to prepare the input,
and then run the model:

```
    cd work-${OS-NAME}_${COMPILER-NAME}
    runprep.sh -dircfg configurations/experience
    runmod.sh -dircfg configurations/experience
```

## Modifying the grid and getting meteorological data}

You can use the script named *grille.sh* in the *scripts* directory to define
your own grid and visualise it with SPI (see above how to get it).

For this, copy one of the *gem_settings.nml* files located in the different
configurations directories, edit it, and then run the command 
```
grille.sh -spi
```

If you want geophysical fields and historical meteorological data for the
region you defined in that new grid, contact us.

## Other Utilities

### pgsm

pgsm is a utility designed to perform horizontal interpolations and basic
arithmetic operations on RPN standard files.

Input files must be RPN standard files. Output files may be RPN standard
files (random or sequential), binary FORTRAN sequential files, or formatted
ASCII files.

PGSM can:
- Interpolate data on various geographical projections, compressed or not.
- Interpolate wind components UU and VV with respect to the scale and orientation of the output grid.
- Perform symmetric or antisymmetric extrapolations from an hemispheric grid.
- Compute thicknesses between 2 levels.
- Compute precipitation amounts between 2 forecast hours.
- Extract low, mid and high clouds.
- Perform mathematical operations on one field or more.
- Compute latitudes and longitudes from the X-Y coordinates of a grid.
- Read navigational records of latitudes-longitudes (grid type Y) or
  grid coordinates (grid type Z) in an input or output file and interpolate
  values from these coordinates.

Example:
````
pgsm -iment <input FST> -ozsrt <output FST> -i <pgsm.directives>
````

### editfst

editfst is a utility used for editing and copying records from RPN standard
files into a new or an existing RPN standard file. It can do a straight
(complete) copy of the input file or it can copy records selectively as
indicated from the standard input or from a file of directives named in the
-i option.

Example:
````
editfst -s <input FST> -d <output FST> -i <editfst.directives>
````
