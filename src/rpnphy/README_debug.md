
Debugging tips for RPNPhy in GEM/SCM
====================================

Compiling with debug options
----------------------------

Several options are available.

### Using GEM/SCM CMake pre-defined debug options

This will enable debug options on the whole GEM code.
```
rm -rf build-${GEM_ARCH}/*
make cmake-debug && make -j && make -j work
```
Alternatively one can use `cmake-debug-extra`

### Adding user specific for rpnphy
This allows to enable debug options only for the rpnphy code or sub code.
One need to edit `src/rpnphy/CMakeLists.txt` in order to update the `CMAKE_Fortran_FLAGS`.  Then rerun CMake (not `cmake-debug`) as:
```
rm -rf build-${GEM_ARCH}/*
make cmake && make -j && make -j work
```

For example adding some options (`-C -g -fp-speculation=safe -init=snan,arrays -traceback -warn all -warn nointerfaces -check noarg_temp_created -std08`) for the intel compiler for the full physics.
```diff
 if (("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel") AND NOT ("${CMAKE_SYSTEM_NAME}" STREQUAL "CrayLinuxEnvironment"))
   if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 2021)
     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qmkl")
-    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qmkl -static-intel -diag-disable 5268")
+    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qmkl -static-intel -diag-disable 5268 -C -g -fp-speculation=safe -init=snan,arrays -traceback -warn all -warn nointerfaces -check noarg_temp_created -std08")
     set(CMAKE_EXE_LINKER_FLAGS_INIT "${CMAKE_EXE_LINKER_FLAGS_INIT} -qmkl")
   else()
     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mkl")
-    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl -static-intel -diag-disable 5268")
+    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl -static-intel -diag-disable 5268 -C -g -fp-speculation=safe -init=snan,arrays -traceback -warn all -warn nointerfaces -check noarg_temp_created -std08")
     set(CMAKE_EXE_LINKER_FLAGS_INIT "${CMAKE_EXE_LINKER_FLAGS_INIT} -mkl")
   endif()
 elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
```

Example for file or sub-directory specific options

[Pending]

### Seeing the full compiler options list

Compiler options are "saved" in the build directory
```
cat build-${GEM_ARCH}/src/rpnphy/rpnphy/CMakeFiles/rpnphy.dir/flags.make
```

Running in the debugger
------------------------------

There are several debuggers. Here are 2 we have at ECCC/CMC.  
Here is the basics on using them with GEM/SCM

### GDB

```
Running in SCM
runscm -exp expname -debug gdb
```

### DDT
In the following example, DDT is used with its graphical interface and thus need to run on an interactive node with a defined `DISPLAY` (see how we can use this on the back/compute node)
```
# Loading
. ssmuse-sh -x /fs/ssm/main/opt/forge/23.1
# Running in SCM
runscm -exp expname -debug ddt
```

Profiling the code
----------------------

### Basic profiler (`gprof`)

Compile with profiling options:

  Add -p/-pg to compiler/linker options

Run the executable:

  Make sure no cleanup is done afterward to keep the gmon.out files
  May want to add to the env:
    export GMON_OUT_PREFIX=gmon.out-

Sum all gmon.out files:

  gprof -s ${gem_DIR}/${GEM_WORK}/bin/maingemdm* $(find  ${gem_DIR}/${GEM_WORK}/RUNMOD/ -name 'gmon.out*')

Visualize:

  gprof -p ${gem_DIR}/${GEM_WORK}/bin/maingemdm* gmon.sum > gmon.sum-p
  gprof -q ${gem_DIR}/${GEM_WORK}/bin/maingemdm* gmon.sum > gmon.sum-q


Ref:
* https://www.nas.nasa.gov/hecc/support/kb/using-gprof-for-performance-analysis_671.html
* https://web.archive.org/web/20200226032601/http://shwina.github.io/2014/11/profiling-parallel

### Vtune profiler
