GEM dynamical core: RPN, MRD, STB, ECCC, GC, CA
===============================================

This document is intended as a notebook and checklist for the GEMDyn
librarian. It is only valid on EC/CMC and GC/Science Networks.


Building, Installing from Version Control
=========================================

Getting the code
----------------

If not already done, you may clone the git repository and checkout
the version you want to run/work on with the following command
(example for version 5.2.0):

NAME=gemdyn
MYVERSION=5.2.0                                   ## Obviously, you'll need to change this to the desired version
MYURL=git@gitlab.science.gc.ca:MIG/${NAME}        ## You'll need a GitLab.science account for this URL
## MYURL=https://gitlab.science.gc.ca/MIG/${NAME}.git ## You cannot "git push" if you use the http:// URL
git clone -b ${NAME}_${MYVERSION} ${MYURL} ${NAME}
cd ${NAME}
git checkout -b ${NAME}_${MYVERSION}-${USER}      ## Create your own branch (example)

Committing the Code
-------------------

**TODO Pre-commit**:

  * merge in code from other devs (and from other branches if any)
  * update version number of gemdyn and other libraries in MANIFEST

**TODO**:

  * commit changes
  * tag version
  * push code and tags upstream (git push && git push --tags)


#### Post Install ####

**TODO**

  * update doc and nml references

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
