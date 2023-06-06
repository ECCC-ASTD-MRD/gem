# GEM-MACH Compilation

## Building the binaries

### For the first time

#### Set GEM environment variables

* Set GEM environment
  
Use r.load.dot followed by GEM version from [README\_version.md](README_version.md):

```
 . r.load.dot GEM/x/5.2.0-b1  ## This is an example for version x/5.2.0-b1. Obviously, you'll need to change the version.
```

* Set environmental variable `storage_model`

```
mkdir your_working_directory
export storage_model=your_working_directory
```
**your\_working\_directory** is a directory you need to create somewhere on your /home/ords/ or the supercomputer /space/hall[1,2]/work/ filesystems. It will be used to store the CHEM library and GEM-MACH final binary.

* Unset environment variable `DEBUGMAKE` to enable default optimization

```
unset DEBUGMAKE
```

#### Prepare CHEM library

* Set the RDE build environment for CHEM

```
cd your_checkout_directory/CHEM
ouv_exp_gemmach   # this will use the storage_model
```

* Compile and build the CHEM library

```
cd your_checkout_directory/CHEM
make -j 8 mach
```
#### Compile GEM-MACH

* Set the RDE build environment for GEM

```
cd your_checkout_directory/GEM
ouv_exp_gem      # this will use the storage_model
```

* Create the symlinks to GEM libraries and binary

```
cd your_checkout_directory/GEM
linkit           # this will use the storage_model
```

* Compile and link GEM binaries to create GEM-MACH binary

```
cd your_checkout_directory/GEM
make dep
make -j 8 obj
make gemdm
```
Note: if you want to free some space without removing the binaries, you can use "make clean" in each directory.  This will delete all the temporary files but keep the important one.  The drawback is that you'll have to compile everything again if you make any modification.


### When you've made some modifications


* If you have modified something in the CHEM:

```
cd your_checkout_directory/CHEM
make mach       # this will compile all the changes at once
```
Note: Before you introduce modifications to CHEM code, please get familiar with [GEM-MACH coding standards](https://wiki.cmc.ec.gc.ca/wiki/GEM-MACH/Coding_standards).

* If you have modified something in the GEM part:

```
cd your_checkout_directory/GEM
make -j 8 obj     # if F90 file is modified
```
OR:
```
cd your_checkout_directory/GEM
make dep     # if new dependencies were added in the F90 file
make -j 8 obj
```
OR:
```
cd your_checkout_directory/GEM
make -j 8 obj      # if cdk file is modified, this command has to be repeated
make -j 8 obj      # first time to re-compile all dependencies and second to check that all objects are up to date
```

Note that you can always use "make the_file.o" to compile each file separately, but that command will not build the libraries and necessary dependencies.

* Regardles which part you have modified, at **the end**:

```
cd your_checkout_directory/GEM
make gemdm
```

* Note on RDE 

To retrieve a file from the GEM/PHY directories, use the following command:

```
cd your_checkout_directory/GEM
rde co name_of_file.extension   # use the Tab key for automatic complition!
```

For more information on how to retrieve GEM files and how to compile after your modifications, visit [RDE tutorial](https://wiki.cmc.ec.gc.ca/wiki/RDE/1.0/Tutorial).

<!--
===To speed things up a bit when compiling chemistry part===
Add the "-j" option followed by the number of processors to be used to "make", e.g. in the CHEM, use 'make -j 8 mach'.
-->


## Debugging mode
If you want to debug the code, there are many ways to do so:

### Use CHEM / GEM debug mode

In RDE environment, `DEBUGMAKE` environment variable is used to trigger
compilation with default debug options. 

When `DEBUGMAKE` is set to 1, it will trigger verbose compilation with no optimization. It can be set for both compiling CHEM and GEM code, or it can be set for one and unset for the other code directory. Example below is for both directories.

```
export DEBUGMAKE=1

cd your_checkout_directory/CHEM
make buildclean
make -j8 mach

cd your_checkout_directory/GEM
make buildclean
make -j8 obj
make gemdm
```

One can to further modify, suppress default debug options set by RDE by using
user supplied `COMP_RULES_FILE` in `Makefile.user.mk`.
This file can be modified with the defult rule file, see:
https://wiki.cmc.ec.gc.ca/wiki/RDE/1.0/Ref#User_overrides

Note: Compilation in debug mode will add some flags to the compilation/linkage and put the code in "verbose" mode so expect slower model execution and bigger listing files.

### Using debug variables

It is possible to create up to 18 2D variables and up to 18 3D variables on the volatile bus for debugging purpose.

#### gem_settings.nml

The memory allocation is done automaticaly and is triggered by 2 keys in gem_settings.nml:
```
 chm_debug_2d_i = 0
 chm_debug_3d_i = 0
```

Default value for both is 0.  When chm\_debug\_2d\_i integer value is greater than 0, that number of 2D variables will be created and allocated on the volatile bus.  Same for chm\_debug\_3d\_i.

#### Output

The output name of those variables will be automaticaly generated.<br>
The first 3 characters will be 2DB for 2D variables, and 3DB for 3D variables.<br>
The forth character will depend on the number of variable requested.<br>
Values will be: 1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, C, D, E, F, G, H, I<br>
For a maximum of 18 values.<br>
E.g.:

* 2DB1, 2DB2, ... 2DBA, 2DBB, ...
* 3DB1, 3DB2, ... 3DBA, 3DBB, ...

The next lines are included in the current outcfg.out file:
 sortie_p([2DB1,2DB2,2DB3] , grid, 1, levels, 1, steps, 1)
 sortie_p([3DB1,3DB2,3DB3] , grid, 1, levels, 1, steps, 1)

