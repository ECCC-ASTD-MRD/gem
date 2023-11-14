# GEM-MACH Compilation

GEM-MACH model consists of the meteorological model (GEM) and chemistry library (MACH) called from physics part of GEM.

In CMAKE compiling environment, GEM source code consists of the set of subdirectories and the MACH code is one of these subdirectories. GEM GitLab repository is a superrepository in which dynamics, physics, chemistry, etc are subtrees of individual repositories. In other words, MACH subdirectory with chemistry source code is tracked in [this GEM-MACH](https://gitlab.science.gc.ca/ARQI-GEMMACH/GEM-MACH) repository, and is contained in the GEM repo as a subtree. Depending on the compilation choices GEM can be compiled with and without chemistry. Here, only compilation with chemistry is discussed.

Depending on the level of interest in GEM code, one can choose between different workflows:
1. If you need to edit the GEM code, or if you want to manually build your binary,, the most convenient workflow is to:
 * Clone both GEM and MACH repositories in two independent directories
 * Work on MACH repo
 * Change directory to the GEM repo
 * `git subtree pull` MACH repo into GEM repo
 * Compile GEM-MACH from the command line while in GEM repo
2. Otherwise, you can automate the GEM-related parts of the build process:
 * Clone only MACH repository
 * Use GEM-MACH integration test to compile, run and validate GEM-MACH against the latest GEM-MACH version available in GitLab.

Note:
There are variations of some of the steps in each of the choices and that below-described workflows are in flux. This readme will be updated as required.

## Building the Binaries for the First Time

### Acquiring GEM-MACH Code

The MACH part of GEM-MACH code can be cloned from [ARQI-GEMMACH/GEM-MACH project](https://gitlab.science.gc.ca/ARQI-GEMMACH/GEM-MACH).

At present, the cmake version is only available through the `3.2` branch. To get to this branch, clone the repo and use `git switch`:

```
cd <local_directory_where_one_keeps_different_versions_of_GEM-MACH_code>
git clone git@gitlab.science.gc.ca:ARQI-GEMMACH/GEM-MACH.git <local_path_to_MACH_CMAKE_git_repository>
cd <local_path_to_MACH_CMAKE_git_repository>
git switch 3.2
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

* For development and testing purposes, it is also possible to simply blow away the MACH subtree and replace it with a link to your MACH repo:

```
cd <local_path_to_GEM_CMAKE_git_repository>
rm -r src/mach
ln -s <local_path_to_MACH_CMAKE_git_repository> src/mach
```

### Compile GEM-MACH


There are three options for compiling GEM-MACH:
1. Using `cmake` commands directly, if you are familiar with `cmake` and have used it in other projects
2. Using RPN's many layers of wrapper scripts that obfuscate the compilation process, but under the hood just call the cmake commands
3. Using integration test that compiles, runs and verifies the results, all as batch processes submitted to the machine the test is initiated on

Either of the two first cases requires setting up environment and compilation from GEM repository. The integration test, on the other hand, can also be ran from the MACH repository.

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
cmake -DWITH_MACH=TRUE ../ 2>&1 | tee ../make.cmake-mach.out
make -j 10 work 2>&1 | tee ../make.work.out
```

* The RPN way:

Run RPN wrapper for cmake directly from the GEM repository:<br>

```
cd <local_path_to_GEM_CMAKE_git_repository>
make cmake-mach 2>&1 | tee make.cmake-mach.out
make -j10 work 2>&1 | tee make.work.out
```

These save the listings of the compilation into `make.cmake-mach.out` and `make.work.out` files in `<local_path_to_GEM_CMAKE_git_repository>` and the GEM-MACH binary at `<local_path_to_GEM_CMAKE_git_repository>/work-${GEM_ARCH}/bin/maingemdm`

#### Compile GEM-MACH running the integration test:

Integration test can be ran from either local MACH repository or from src/mach directory of local GEM super repository.

* From local MACH repository:

Use the script that behind the scene, through the set of batch jobs, clones the appropriate tag version from GEM repo, git subtree pulls local MACH repo into the cloned GEM repo, compiles GEM-MACH, runs it and validates the results against the latest version of MACH in the remote 3.2 branch.

```
cd <local_path_to_MACH_CMAKE_git_repository>
tools/gm-integration-test/initialize-gm-integration-test.sh -m
```

This will create a directory `/space/<current_hall>/sitestore/eccc/aq/r1/${USER}/maestro/${TRUE_HOST}/<test_version>/gm-test_<MACH_git_directory_basename>_<MACH_git_commit_hash>` and a link to it from `<local_path_to_MACH_CMAKE_git_repository>` named `gm-test-${TRUE_HOST}_<MACH_git_commit_hash>`. The binary will be `gm-test-${TRUE_HOST}_<MACH_git_commit_hash>/build/work-${GEM_ARCH}/bin/maingemdm`.

* From local GEM super repository:

Use the script that behind the scene, through the set of batch jobs, clones local GEM super repository with MACH already pulled in, compiles GEM-MACH, runs it and validates the results against the latest version of MACH in the remote 3.2 branch.

```
cd <local_path_to_GEMM_CMAKE_git_repository>/src/mach
tools/gm-integration-test/initialize-gm-integration-test.sh
```

This will create a directory `/space/<current_hall>/sitestore/eccc/aq/r1/${USER}/maestro/${TRUE_HOST}/<test_version>/gm-test_<GEM_git_directory_basename>_<GEM_git_commit_hash>` and a link to it from `<local_path_to_GEM_CMAKE_git_repository>/src/mach` named `gm-test-${TRUE_HOST}_<GEM_git_commit_hash>`. The binary will be `gm-test-${TRUE_HOST}_<GEM_git_commit_hash>/build/work-${GEM_ARCH}/bin/maingemdm`.

For detailed information on integration test, please read the comments at the beginning of the scripts in `tools/gm-integration-test`, with special attention to the beginning of `tools/gm-integration-test/initialize-gm-integration-test.sh`. `-h` option of the integration test also provides basic help.


## Building the Binaries After Making Code Modifications

* If code in `<local_path_to_GEM_CMAKE_git_repository>` is edited, repeat compilation:

```
cd <local_path_to_GEM_CMAKE_git_repository>
make -j10 work 2>&1 > make.work.out
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
make -j10 work 2>&1 > make.work.out
```

The same caveats about possibly needing to re-run cmake apply here.

 * Use the integration test:
  1. commit the changes in `<local_path_to_MACH_CMAKE_git_repository>`
  2. run integration test

```
cd <local_path_to_MACH_CMAKE_git_repository>
git add <files_to_be_committed>
git commit
tools/gm-integration-test/initialize-gm-integration-test.sh
```

Note: Before you introduce modifications to MACH code, please get familiar with [GEM-MACH coding standards](https://wiki.cmc.ec.gc.ca/wiki/GEM-MACH/Coding_standards).

## Debugging mode

One can compile in a debug mode by replacing `cmake-mach` with `cmake-mach-debug` when compiling on command line using RPN wrappers, or utilizing `-d` option with the integration test.

Binary produced in debug mode will require updating the model namelist that will be used at the run time by adding `chm_debug_trace_l=.true.` in `physics_cfgs` and `chemistry_cfgs` namelist-keys groups of `gem_settings.nml`. This type of binary will also produce large Runmod listing file.

# Debug variables

Depending on the debugging needs, one can assign values of internal model variables to the debug variables. To output these debug variables, besides adding the request in `outcfg.out`, one has to set non-default values for particular namelist keys in `gem_settings.nml`.

It is possible to create up to 18 2D variables and up to 18 3D variables on the volatile bus for debugging purpose.

### gem\_settings.nml

Internal mechanism for automatic allocation of memory for these debug variables is triggered by setting the value of the following manelist keys to the number of desired 2D and 3D debug variables:

```
 chm_debug_2d_i
 chm_debug_3d_i
```

Default value for both keys is 0. To request an automatic allocation of memory, i.e. variable creation and space allocation on the volatile bus, set the desired key value to be an integer greater than 1:

* For 2D variables, set `chm_debug_2d_i` to the desired number of 2D debug variables
* For 3D variables, set `chm_debug_3d_i` to the desired number of 3D debug variables

### outcfg.out

Output names for these variables are automatically generated as follows:
* The first 3 characters are `2DB` for 2D variables, and `3DB` for 3D variables.
* The forth character depends on the number of variables requested.
 * For a maximum of 18 variables, the values are `1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, C, D, E, F, G, H, I`
 * For any number less than 18, the values will be a subset from the above list, starting with 1 and ending with the character corresponding to the number requested in the namelist key

An example of the lines that need to be included in the `outcfg.out` file:

```
 sortie_p([2DB1,2DB2,2DB3] , grid, 1, levels, 1, steps, 1)
 sortie_p([3DB1,3DB2,3DB3] , grid, 1, levels, 1, steps, 1)
```
