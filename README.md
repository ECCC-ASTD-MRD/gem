ModelUtils, RPN Models utility functions, library and tools: RPN, MRD, STB, ECCC, GC, CA
========================================================================================

This document is intended as a notebook and checklist for the ModelUtils
librarian. It is only valid on EC/CMC and GC/Science Networks.

Quick Start: Building, Installing from Version Control
======================================================

This section is just a condensed version of the steps needed to release from
the git repository. It implies the mandatory steps for your SHELL env. setup
and directories structure are already done. It also implies that the code is
fully committed in the git repository and no other changes are needed.
These steps needs to be done in the specified order.

**TODO**: 

General Information
===================

**TODO**: 

Basic profile and directory setup
---------------------------------

**TODO**: 

  * list expected dir structure
  * list expected env var and paths...
  * gmake based
  * supported compilers

Layout
------

This repository is built following some conventions for
RDE (RPN Development Environment) and MIG (Model Infrastructure Group).

It contains the following sub-directories:

  * `bin/       `: used for scripts, will be added to the `PATH`
  * `include/   `: used for included files, will be added to the `INCLUDE_PATH`
  * `lib/       `: will be added to the `LIBRARY_PATH`
  * `lib/python/`: used for python modules, will be added to the `PYTHONPATH`
  * `src/       `: used for source code, will be added to the source path (`VPATH`).
  * `share/     `: used for any other content

It contains the following files:

  **TODO: Complete list with mandatory and optional**
  * `.setenv.dot`:
  * `bin/.env_setup.dot`:
  * `VERSION`:

  **TODO: Complete list of Mandatory Makefile vars and targets**

Closing Issues
--------------

**TODO**: review/close bugzilla issues


Documentation update
--------------------

**TODO**: review doc on wiki


Building, Installing from Version Control
=========================================

Getting the code
----------------

If not already done, you may clone the git repository and checkout
the version you want to run/work on with the following command
(example for version 1.7.0):

        NAME=modelutils
        MYVERSION=1.7.0                                   ## Obviously, you'll need to change this to the desired version
        MYURL=git@gitlab.science.gc.ca:MIG/${NAME}        ## You'll need a GitLab.science account for this URL
        ## MYURL=https://gitlab.science.gc.ca/MIG/${NAME}.git ## You cannot "git push" if you use the http:// URL
        git clone ${MYURL} ${NAME}_${MYVERSION} ${NAME}
        cd ${NAME}
        git checkout -b ${NAME}_${MYVERSION}-${USER}-branch ${NAME}_${MYVERSION}


Update Dependencies
-------------------

**TODO**:

  * Where/How to set change external dependencies (compiler, librmn, vgrid, rpn_comm, ...)
  * Where/How to set compiler options, new compiler rules...


Committing the Code
-------------------

**TODO Pre-commit**:

  * merge in code from other devs (and from other branches if any)
  * update version number
  * update bndl dependencies
  * test (testing section below)

**TODO**:

  * commit changes
  * tag version
  * push code and tags upstream (git push && git push --tags)


Initial Setup
-------------

Setting up the Shell Environment and Working Directories.

The compiling, building and running systems depend on a few Shell
(Bash is the only supported Shell to date) environment variables,
PATHs, directories, links and files.
The following commands perform that setup.

#### Shell Env. SetUp ####

**TODO**:

#### Files, Dir

**TODO**:

Compiling and Building
----------------------

**TODO**:

Testing
-------

**TODO**:

Install
-------

Installation can only be performed by the librarian who has proper permissions.
This would be done on the EC/CMC and GC/Science networks.

**TODO**

#### Post Install ####

**TODO**

  * update doc
  * send emails


Uninstall
---------

Un-installation can only be performed by the librarian who has proper permissions.
This would be done on the EC/CMC and GC/Science networks.

**TODO**

Cleaning up
-----------

To remove all files created by the setup, compile and build process, use the `distclean` target.

        make distclean

Misc
====

Patching the git repo
---------------------

Ref: https://www.devroom.io/2009/10/26/how-to-create-and-apply-a-patch-with-git/

Patch are produced with:

       BASETAG=   #Need to define from what tag (or hash) you want to produce patches
       git format-patch HEAD..${BASETAG}

Before applying the patch, you may check it with:

       git apply --stat PATCH.patch
       git apply --check PATCH.patch

Full apply

       git am --signoff PATCH.patch

Selective application, random list of commands

       git apply --reject PATH/TO/INCLUDE PATCH.patch
       git apply --reject --include PATH/TO/INCLUDE PATCH.patch
       git am    --include PATH/TO/INCLUDE  PATCH.patch
       git apply --exclude PATH/TO/EXCLUDE PATCH.patch
       git am    --exclude PATH/TO/EXCLUDE  PATCH.patch

Fixing apply/am problems

Ref: https://stackoverflow.com/questions/25846189/git-am-error-patch-does-not-apply
Ref: https://www.drupal.org/node/1129120

With --reject: 

  * inspect the reject
  * apply the patch manually (with an editor)
  * add file modified by the patch (git add...)
  * git am --continue
