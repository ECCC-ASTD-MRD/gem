# GEM-MACH Compilation

## Introduction

Our build process is somewhat complicated by the fact that MACH is a chemistry
library (as that of rpn-physics and gem-dynamics) and part of the GEM (runmod)
compiled binary.  The source code structure of MACH is only to build this
library.  The compilation of workable GEM-MACH binary therefore requires the
full GEM package.

In the current official distribution of mig/GEM.  It contains the default,
official release of MACH code (this is the default branch of the GEM-MACH
repository).  If this is the version of MACH code of interest, you can simply
clone the mig/GEM code and compile utilizing the steps outlined in the GEM's
readme (also see detail steps below).

If you wish to replace the default MACH chemistry library with other MACH code,
either from your own repository or other branches in the GEM-MACH project repo.,
you will need additional steps of replacing the official MACH sub-tree code from
the mig/GEM repository, and linking the mach source code with a new one using
`git-subtree pull` steps prior to compilation (see detail steps below).

Note that, given the MACH library dependency on `rpnphy`, `gemdyn`, there are
instances where code changes are necessary (e.g. for AQ feedback).  These
changes are made directly on their original code of current version and carried
in the `src/gemmod` directory.  These source files are compiled at the same time
as MACH library and object files over-written at the final linking of binary.

Due to the work-around of `src/gemmod`, all source files are required to
be listed in the `src/sourcelist.cmake` and is used during the building process.
The content of the file contains current default release.  Please make relevant
changes if you decide to use a different MACH code branch.

## Quickstart Guide

### Default release version:

Compiling the official mig/GEM release (i.e. current version used in RAQDPS).

Example steps below uses the latest release of GEM 5.2 code (`GEM_BRANCH='5.2'`)
which already includes MACH code from the `3.2` branch of GEM-MACH repository.

```
cd <target directory , e.g., /space/hall5/sitestore/eccc/aq/r1/${USER}/code>
GEM_DIR=gem
GEM_BRANCH='5.2'

# Clone mig/GEM project repo (i.e. super GEM repo)
git clone git@gitlab.science.gc.ca:MIG/gem.git ${GEM_DIR}
cd ${GEM_DIR}
git switch ${GEM_BRANCH}
# Set environment
source ./.eccc_setup_intel
source ./.initial_setup
# Build
cd build-${GEM_ARCH}
cmake -DWITH_MACH=TRUE ../ |& tee ../make.cmake-mach.out
# see "Makefile" wrapper for additional directives
# or "cmake -DCMAKE_BUILD_TYPE=Debug -DWITH_MACH=TRUE" to compile with debug
# Compile and link
make -j work |& tee ../make.work.out
```

The final compiled binary `maingemdm` and other tools exist in `work-${GEM_ARCH}/bin`.

### Other MACH code branch:

To compile MACH library from other sources that is not already as part of
mig/GEM official release, one has to replace the `src/mach` code with your own
version following similar directory structure.

We handle this dependency using `git subtrees` after the full mig/GEM repository
is already available (see earlier steps).

In the example below, the default MACH code is git removed, committed, and
replace with MACH code from the `3.2dev` branch.  One can replace the target MACH
code with other code branch or other remote/local MACH repo.


```
GEM_DIR=gem52mach32dev
GEM_BRANCH='5.2'
MACH_BRANCH='3.2dev'

git clone git@gitlab.science.gc.ca:MIG/gem.git ${GEM_DIR}
cd ${GEM_DIR}
git switch ${GEM_BRANCH}
git checkout -b 5.2mach3.2dev # switch to a new branch name
git rm -rf src/mach
git commit -m "to replace mach subtree with code from git@gitlab.science.gc.ca:ARQI-GEMMACH/GEM-MACH.git ${MACH_BRANCH}"
git subtree add --prefix=src/mach git@gitlab.science.gc.ca:ARQI-GEMMACH/GEM-MACH.git ${MACH_BRANCH}
## optionally push the replaced mach GEM super repo to your own gitlab repo.
#This step assumes you have forked the mig/GEM project to your own space.
git remote add mygem git@gitlab.science.gc.ca:${USER}/gem.git
git push mygemrepo

# now the MACH code is replace, compile following the same steps
source ./.eccc_setup_intel
source ./.initial_setup
cd build-${GEM_ARCH}
cmake -DWITH_MACH=TRUE ../ |& tee ../make.cmake-mach.out
make -j work |& tee ../make.work.out

```
The final compiled binary `maingemdm` and other tools exist in `work-${GEM_ARCH}/bin`.


## Detailed Info on Building the Binaries for the First Time

### Acquiring GEM-MACH Code

The MACH part of GEM-MACH code can be cloned from [ARQI-GEMMACH/GEM-MACH project](https://gitlab.science.gc.ca/ARQI-GEMMACH/GEM-MACH).

```
cd <local_directory_where_one_keeps_different_versions_of_GEM-MACH_code>
git clone git@gitlab.science.gc.ca:ARQI-GEMMACH/GEM-MACH.git <local_path_to_MACH_CMAKE_git_repository>
cd <local_path_to_MACH_CMAKE_git_repository>
```

### Acquiring GEM code

The GEM code, which is necessary for compiling GEM-MACH, can be cloned from the [MIG/gem project](https://gitlab.science.gc.ca/MIG/gem).

If you are working on a feature that involves the GEM code, we recommend to start from the latest version of the code. If you're just working on chemistry, we recommend resetting the repo to the commit that is tagged as the GEM version that corresponds to the chemistry code you are working on.

```
git clone git@gitlab.science.gc.ca:MIG/gem.git <local_path_to_GEM_CMAKE_git_repository>
cd <local_path_to_GEM_CMAKE_git_repository>
git reset --hard $(git rev-list -n 1 ${GEM_version_tag})
```

### Making connection between GEM and GEM-MACH repositories

#### Checking if the MACH code in GEM repository is same as in GEM-MACH repository

There are two ways for this check:
1. Compare `<local_path_to_MACH_CMAKE_git_repository>/src` and `<local_path_to_GEM_CMAKE_git_repository>/src/mach` directories (e.g. by using `xxdiff -BbwirD `)
2. Compare the following two hashes:
 * latest commit in 3.2 branch of `<local_path_to_MACH_CMAKE_git_repository>` or [ARQI-GEMMACH/GEM-MACH](https://gitlab.science.gc.ca/ARQI-GEMMACH/GEM-MACH) repository
 * latest commit that is specified in the latest commit that starts with "Squashed 'src/mach/'" in `<local_path_to_GEM_CMAKE_git_repository>` or [MIG/gem](https://gitlab.science.gc.ca/MIG/gem) repository

#### Connecting MACH with GEM repository

Depending on an interest to keep the git structure of the repositories, there are two ways of connecting MACH repository to GEM repository:

* If the MACH code in GEM repository is not the same as in GEM-MACH repository, or one wants to compile some other branch from `<local_path_to_MACH_CMAKE_git_repository>` repository instead of 3.2, git subtree pull choice of MACH into GEM repository:

```
cd <local_path_to_GEM_CMAKE_git_repository>
git subtree pull --squash -P src/mach <local_path_to_MACH_CMAKE_git_repository> <name_of_branch_chosen_from_local_MACH_CMAKE_git_repository>
```

Note that git subtree pull of the local repository works only if both GEM and MACH repositories are on the same hall (either on hall5 or hall6).

* For development and testing purposes, it is also possible to simply blow away
  the MACH subtree and replace it with a copy of your MACH code directory.
* This is more a hardlink/copy override and assumes that you won't be using git
  tools to version and commit your changes.

```
cd <local_path_to_GEM_CMAKE_git_repository>
rm -r src/mach
cp -a <local_path_to_MACH_CMAKE_git_repository> src/mach
```

### Compile GEM-MACH

There are three options for compiling GEM-MACH:
1. Using `cmake` commands directly, if you are familiar with `cmake` and have used it in other projects
2. Using RPN's many layers of wrapper scripts that obfuscate the compilation process, but under the hood just call the cmake commands
3. Using regression test that compiles, runs and verifies the results, all as batch processes submitted to the machine the test is initiated on

Either of the two first cases requires setting up environment and compilation from GEM repository. The regression test, on the other hand, can also be ran from the MACH repository.

#### Set necessary environment:

```
cd <local_path_to_GEM_CMAKE_git_repository>
source ./.eccc_setup_intel
source ./.initial_setup
```

#### Compile GEM-MACH code on the command line:

* The `cmake` way:

As with a standard `cmake` workflow, run `cmake` command from build subdirectory of GEM repo, using the flag `-DWITH_MACH=TRUE`. The default setup for GEM is to use `build-${GEM_ARCH}` for build directory, and on U2, it translates to `build-rhel-8-icelake-64-intel-2022.1.2` subdirectory:

```
cd <local_path_to_GEM_CMAKE_git_repository>/build-${GEM_ARCH}
cmake -DWITH_MACH=TRUE ../ |& tee ../make.cmake-mach.out
make -j work |& tee ../make.work.out
```

* The RPN way:

Run `Makefile` wrapper for cmake directly from the GEM repository:<br>

```
cd <local_path_to_GEM_CMAKE_git_repository>
make cmake-mach |& tee make.cmake-mach.out
make -j work |& tee make.work.out
```

These save the listings of the compilation into `make.cmake-mach.out` and `make.work.out` files in `<local_path_to_GEM_CMAKE_git_repository>` and the GEM-MACH binary at `<local_path_to_GEM_CMAKE_git_repository>/work-${GEM_ARCH}/bin/maingemdm`

#### Compile GEM-MACH running the regression test:

Regression test can be run from either local MACH repository or from `src/mach` directory of local GEM super repository.

* From local MACH repository:

Use the script that behind the scene, through the set of batch jobs, clones the appropriate tag version from GEM repo, git subtree pulls local MACH repo into the cloned GEM repo, compiles GEM-MACH, runs it and validates the results against the latest version of MACH in the remote 3.2 branch.

```
cd <local_path_to_MACH_CMAKE_git_repository>
tools/regression-test/initialize-regression-test.sh -m -u
```

This will create a directory `/space/<current_hall>/sitestore/eccc/aq/r1/${USER}/maestro/${TRUE_HOST}/<test_version>/gm-test_<MACH_git_directory_basename>_<MACH_git_commit_hash>` and a link to it from `<local_path_to_MACH_CMAKE_git_repository>` named `gm-test-${TRUE_HOST}_<MACH_git_commit_hash>`. The binary will be `gm-test-${TRUE_HOST}_<MACH_git_commit_hash>/build/work-${GEM_ARCH}/bin/maingemdm`.

* From local GEM super repository:

Use the script that behind the scene, through the set of batch jobs, clones local GEM super repository with MACH already pulled in, compiles GEM-MACH, runs it and validates the results against the latest version of MACH in the remote 3.2 branch.

```
cd <local_path_to_GEMM_CMAKE_git_repository>/src/mach
tools/regression-test/initialize-regression-test.sh
```

This will create a directory `/space/<current_hall>/sitestore/eccc/aq/r1/${USER}/maestro/${TRUE_HOST}/<test_version>/gm-test_<GEM_git_directory_basename>_<GEM_git_commit_hash>` and a link to it from `<local_path_to_GEM_CMAKE_git_repository>/src/mach` named `gm-test-${TRUE_HOST}_<GEM_git_commit_hash>`. The binary will be `gm-test-${TRUE_HOST}_<GEM_git_commit_hash>/build/work-${GEM_ARCH}/bin/maingemdm`.

For more information on regression test, please read [readme/README\_regression.md](readme/README_regression.md).


## Building the Binaries After Making Code Modifications

* If code in `<local_path_to_GEM_CMAKE_git_repository>` is edited, repeat compilation:

```
cd <local_path_to_GEM_CMAKE_git_repository>
make -j work |& tee make.work.out
```

Re-running the cmake step should not be necessary unless you've made structural changes to the code. If you notice something weird happening, you can go back to the `make cmake-mach` step.

* If code in `<local_path_to_MACH_CMAKE_git_repository>` is edited, follow one of the two options:

 * Use `cmake`:
  1. commit the changes in `<local_path_to_MACH_CMAKE_git_repository>`
  2. pull in this new version into `<local_path_to_GEM_CMAKE_git_repository>` by repeating `git subtree pull`
  3. repeat the compilation

```
cd <local_path_to_MACH_CMAKE_git_repository>
git add <files_to_be_committed>
git commit
cd <local_path_to_GEM_CMAKE_git_repository>
git subtree pull --squash -P src/mach git@gitlab.science.gc.ca:ARQI-GEMMACH/GEM-MACH.git <name_of_branch_chosen_from_local_MACH_CMAKE_git_repository>
make -j work |& tee make.work.out
```

The same caveats about possibly needing to re-run cmake apply here.

 * Use the regression test:
  1. commit the changes in `<local_path_to_MACH_CMAKE_git_repository>`
  2. run regression test

```
cd <local_path_to_MACH_CMAKE_git_repository>
git add <files_to_be_committed>
git commit
tools/regression-test/initialize-regression-test.sh
```

Note: Before you introduce modifications to MACH code, please get familiar with [GEM-MACH coding standards](https://wiki.cmc.ec.gc.ca/wiki/GEM-MACH/Coding_standards).

## Debugging mode

### Debug Compilation

One can compile in a debug mode by replacing `cmake-mach` with `cmake-mach-debug`
when compiling on command line using `Makefile`, or utilizing `-d` option with
the regression test. The debug binary runs slower.

### Debug Tracing

One can keep track of routines called during runtime regardless of how the code
is compiled by utilizing the following namelist options:
```
&chemistry_cfgs
 chm_debug_trace_l = .true.
```
activates debug tracing only in chemistry and
```
&physics_cfgs
 debug_trace_L = .true.
```
activates debug tracing in physics and chemistry (as a part of physics).

In regression test scripts, the two debug tracing modes are controlled
with options `-c` and `-p`, respectively.

Both of these options slow down model execution and produce large listing files.

### Utilize Debug Variables

Depending on the debugging needs, one can assign values of internal model
variables to the debug variables. To output these debug variables, besides adding
the request in `outcfg.out`, one has to set non-default values for particular
namelist keys in `gem_settings.nml`.

It is possible to create up to 18 2D variables and up to 18 3D variables on the
volatile bus for debugging purpose.

***gem\_settings.nml***

Internal mechanism for automatic allocation of memory for these debug variables
is triggered by setting the value of the following namelist keys to the number
of desired 2D and 3D debug variables:

```
 chm_debug_2d_i
 chm_debug_3d_i
```

Default value for both keys is 0. To request an automatic allocation of memory,
i.e. variable creation and space allocation on the volatile bus, set the desired
key value to be an integer greater than 1:

* For 2D variables, set `chm_debug_2d_i` to the desired number of 2D debug variables
* For 3D variables, set `chm_debug_3d_i` to the desired number of 3D debug variables

***outcfg.out***

Output names for these variables are automatically generated as follows:
* The first 3 characters are `2DB` for 2D variables, and `3DB` for 3D variables.
* The forth character depends on the number of variables requested.
  * For a maximum of 18 variables, the values are `1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, C, D, E, F, G, H, I`
  * For any number less than 18, the values will be a subset from the above list,
    starting with 1 and ending with the character corresponding to the number
    requested in the namelist key

An example of the lines that need to be included in the `outcfg.out` file:

```
 sortie_p([2DB1,2DB2,2DB3] , grid, 1, levels, 1, steps, 1)
 sortie_p([3DB1,3DB2,3DB3] , grid, 1, levels, 1, steps, 1)
```

# See Also
* https://gitlab.science.gc.ca/ARQI-GEMMACH/GEM-MACH/-/wikis/home
* https://wiki.cmc.ec.gc.ca/wiki/GEM-MACH
