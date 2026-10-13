
Updating the rpnphy depot for a GEM release
===========================================

1-Prepare the rpnphy code in the GEM repository
-----------------------------------------------

> _All the following steps (in the section 1) are done
> in your user GEM development clone (created by gemdev)._

... include code from contrubutors & test ...
... see with RPN-SI if there are updates needed to modelutils' CMakeLists.txt ...


### 1.1-Tests

Make sure to test with
* intel
* gnu
* nvhpc
* intel cmake-debug-extra

1st shell
```
. ./.ssmuse_gem intel
. ./.initial_setup
make cmake
make -j work
# ... run tests...
```

2nd shell
```
. ./.ssmuse_gem gnu
. ./.initial_setup
make cmake
make -j work
# ... run tests...
```

3rd shell
```
. ./.ssmuse_gem nvhpc
. ./.initial_setup
make cmake
make -j work
# ... run tests...
```

4th shell
```
. ./.ssmuse_gem intel
. ./.initial_setup
rm -rf build-${GEM_ARCH}/*
make cmake-debug-extra
make -j work
# ... run tests...
```

### 1.2-Finalize

Note that all the steps here assume that we are moving 
from a release version (`$ATM_MODEL_VERSION`). 
If you are doing the release from a tag, then simply 
set `ATM_MODEL_VERSION=[TAG]` and continue with these 
steps.

**1.2.1 Check to see if there are changes to this component.**
```
component=modelutils
git diff --name-only ${ATM_MODEL_VERSION} -- src/${component}
```

**1.2.2 Check if there are merge commits**  
_We need to "squash the merges" before proceeding_
```
FROM=${ATM_MODEL_VERSION:-__scrap__}
git rev-list --min-parents=2 --count ${FROM}..HEAD
```
If >0 then need rebase to linearize history (replace merge commits with "rebase").  
Good luck!

**1.2.3 Update MANIFEST for VERSION**  
_Importante Note: *do not reload env after this**)_
```
emacs src/modelutils/MANIFEST
```

**1.2.4 Commit && Tag**

# Commit
```
component=modelutils
MYVERSION=$(grep ^VERSION src/${component}/MANIFEST | cut -d: -f2)
MYVERSION=${MYVERSION# *}
git commit -a -m "${component}: update VERSION (${component}_${MYVERSION})"
git tag gem_${ATM_MODEL_VERSION}+${component}_${MYVERSION}
```

2-Move patches from GEM dev to modelutils depot
-------------------------------------------

> _All the following steps (in the section 2) are done in your modelutils clone._

**2.1 Create a local clone of the modelutils repo.**
```
git clone git@gitlab.science.gc.ca:MIG/modelutils.git
cd modelutils
```

**2.2 Make sure you are on the right branch & version**
```
git checkout modelutils_1.8-branch
gitk
```
Visually check the your branch (yellow dot) is "ahead"/"on top" of the modelutils version/tag you are 
developping from (e.g. `modelutils_1.8.0-a22`). This is original tag, accessible through `grep modelutils MANIFEST` in the GEM experiment above.

For the first alpha you'll need to create a new branch, for example for the new 1.9.0-a01:
```
git checkout -b modelutils_1.9-branch
```

**2.3 Pull patches from the GEM dev experiment into the local repo.**  
_Set the path and tags in vars bellow_
```
component=modelutils
gemurldir=     # URL or path of the repos where the new modelutils code is (dev experiment above)
FROM_VERSION=  # Original GEM version number on top of which modifications
               # were made (`${ATM_MODEL_VERSION}` in gemdev repos)
gemtag=        # tag of the code to pull from in the repos where the new
               # modelutils code is (above created
               # `gem_${ATM_MODEL_VERSION}+${component}_${MYVERSION}`)
```
_Pull the code from the GEM repos you "finalized" in step 1 above._
```
. r.load.dot eccc/mrd/rpn/MIG/ENV/modeltools/tools
t.gitpull --relative src/${component} -- ${gemurldir} ${gemtag} ${FROM_VERSION}
```

**2.4 Tag**
```
modelutils_version=$(grep ^VERSION MANIFEST | cut -d: -f2)
modelutils_version=${modelutils_version# *}
mubranch=${modelutils_version%.*}
mubranch=modelutils_${mubranch:-6.3}-branch
git tag modelutils_${modelutils_version}
```

**2.5 Push the code and tag to the official repos**  
_Visually inspect the modelutils repos (using `gitk`) to make 
sure everything is in order before pussing to the official repos._
```
git push origin ${mubranch} modelutils_${modelutils_version}
```

3-Documentation
---------------

* Update MIG/modelutils issue and milestones on gitlab by closing them.
* Update modelutils doc
  * Dest: http://wiki.cmc.ec.gc.ca/wiki/modelutils


Appendix: manual patch transfer
==============================

Steps to be done in a GEM dev env.
----------------------------------
### Initial check

> **WARNING: we need to "squash the merges" before proceeding Check if there are merge commits**
```
FROM=${ATM_MODEL_VERSION:-__scrap__}
git rev-list --min-parents=2 --count ${FROM}..HEAD
```
If >0 then need rebase to linearize history.  Good luck.

### Create the patches
```
FROM=${ATM_MODEL_VERSION:-__scrap__}
component=modelutils
# git subtree split --prefix=src/${component}
mkdir -p patches/${component}
git format-patch -o patches/${component} \
    -M -C --find-copies-harder -k --relative=src/${component} \
    ${FROM}..HEAD
``` 

> You might want to do your modelutils update before leaving GEM env
> open `src/modelutils/README_librarian.md` and follow instructions

Steps to be done in the modelutils clone
------------------------------------
### Import the patches
```
mydir=     #<- should be the ${gem_DIR} of the dev directory above
cd ..      #create the clone outside your gem experiment
git clone git@gitlab.science.gc.ca:MIG/modelutils.git
cd modelutils
modelutils_version=$(grep ^VERSION MANIFEST | cut -d: -f2)
modelutils_version=${modelutils_version# *}
mubranch=${modelutils_version%.*}
mubranch=modelutils_${mubranch:-6.3}-branch
git checkout ${mubranch}
component=modelutils
for item in ${mydir}/patches/${component}/*.patch ; do
    if [[ -s ${item} ]] ; then
       git am -3 -k --empty=drop ${item}
       if [[ $? != 0 ]] ; then
          echo "ERROR: problem applying ${item}"
          break
       fi
    else
       echo "Ignoring: ${item}"
    fi
done
```

### Tag & push
```
modelutils_version=$(grep ^VERSION MANIFEST | cut -d: -f2)
modelutils_version=${modelutils_version# *}
mubranch=${modelutils_version%.*}
mubranch=modelutils_${mubranch:-6.3}-branch
git tag modelutils_${modelutils_version}
git push origin ${mubranch} modelutils_${modelutils_version}
```
