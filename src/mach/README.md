Here, we provide general information about GEM-MACH project. Specific
information on compilation, configuration, etc are provided in
[readme](https://gitlab.science.gc.ca/ARQI-GEMMACH/GEM-MACH/readme)
directory:
- Compilation: [README_compilation.md](https://gitlab.science.gc.ca/ARQI-GEMMACH/GEM-MACH/readme/README_compilation.md)
- Configuration: [README_namelist.md](https://gitlab.science.gc.ca/ARQI-GEMMACH/GEM-MACH/readme/README_namelist.md)
- Available output fields: [README_variables](https://gitlab.science.gc.ca/ARQI-GEMMACH/GEM-MACH/readme/README_variables.md)
- Regression test: [README_regression.md](https://gitlab.science.gc.ca/ARQI-GEMMACH/GEM-MACH/-/blob/3.2/readme/README_regression.md)

# GEM-MACH

GEM-MACH model is composed of the GEM meteorological model and the MACH
chemistry library, which is integrated through the physics component of
GEM. [This project](https://gitlab.science.gc.ca/ARQI-GEMMACH/GEM-MACH/)
congregates development of MACH library.

## Code Structure

In the CMake-based build environment, the GEM source code is structured
in multiple subdirectories representing its components: dynamics,
physics, etc. The MACH chemistry library is an atmospheric gas and
aerosol chemistry module within RPN physics that also carries some
modifications to dynamics and physics components of GEM. As such, the
code of MACH chemistry library exists as one of the GEM subdirectories.

In the git version control system:
- GEM source code is organized as super repository of multiple subtrees
and submodules
- MACH chemistry library is a git subtree of GEM super repository with
its own git repository

In GitLab:
- [GEM-MACH](https://gitlab.science.gc.ca/ARQI-GEMMACH/GEM-MACH)
project contains repository of MACH chemistry library
- [GEM](https://gitlab.science.gc.ca/MIG/gem) project contains
super repository of GEM

Depending on build options, the GEM model can be compiled either with
or without chemistry. The MACH library cannot be compiled independently
from GEM.

## Development

Depending on the level of interest in GEM code, one can choose between
different GEM-MACH-development workflows:

* If you need to edit the GEM code, or if you want to manually build your
binary, the most convenient workflow is to:
  * Clone both GEM and MACH repositories in two independent directories
  * Work on MACH repo in MACH directory
  * Change directory to the GEM repo
  * git subtree pull MACH repo into GEM repo
  * Optionally, work on GEM repo
  * [Compile GEM-MACH](https://gitlab.science.gc.ca/ARQI-GEMMACH/GEM-MACH/-/blob/3.2/readme/README_compilation.md)
    from the command line while in GEM repo

* Otherwise, you can automate the GEM-related parts of the build process:
  * Clone only MACH repository
  * Use GEM-MACH [regression test](https://gitlab.science.gc.ca/ARQI-GEMMACH/GEM-MACH/-/blob/3.2/readme/README_regression.md)
    to compile using available scripts while in MACH repo. This option
    also runs the compiled executable and validates the output against the
    provided control output.

To run the GEM-MACH model, Maestro sequencer is used:
- [gemmach_maestro](https://gitlab.science.gc.ca/gemmach_maestro)
contains projects pertinent to running the GEM-MACH
- [gemmach module](https://gitlab.science.gc.ca/gemmach_maestro/gemmach_module)
is Maestro module that prepares inputs required for running GEM-MACH
- [gem module](https://gitlab.science.gc.ca/MIG/gem-maestro.git)
is Maestro module that runs GEM and is used to run GEM-MACH

# See Also
* For GEM-MACH related information:
  * [MACH repo wiki](https://gitlab.science.gc.ca/ARQI-GEMMACH/GEM-MACH/-/wikis/home)
  * [GEM-MACH wiki](https://wiki.cmc.ec.gc.ca/wiki/GEM-MACH)
  * [GEM-MACH coding standards](https://wiki.cmc.ec.gc.ca/wiki/GEM-MACH/Coding_standards)
* [GEM wiki](https://wiki.cmc.ec.gc.ca/wiki/GEM) for GEM related information
* [HPCS xwiki](https://portal.science.gc.ca/xwiki/bin/view/Main/) for login/environment related information

