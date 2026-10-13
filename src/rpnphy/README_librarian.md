
Updating the rpnphy depot for a GEM release
===========================================

1-Prepare the rpnphy code in the GEM repository
-----------------------------------------------

> _All the following steps (in the section 1) are done
> in your user GEM development clone (created by gemdev)._

... include code from contrubutors & test ...
... see with RPN-SI if there are updates needed to rpnphy's CMakeLists.txt ...


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
component=rpnphy
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
emacs src/rpnphy/MANIFEST
```

**1.2.4 Create nml_ref**
```
. ./.ssmuse_gem
make cmake
make -j work
which prphynml || export PATH=${PATH}:${gem_DIR}/${GEM_WORK}/bin
component=rpnphy
MYVERSION=$(grep VERSION src/${component}/MANIFEST | cut -d: -f2)
MYVERSION=${MYVERSION# *}
reffilename=src/${component}/share/nml_ref/${component}_settings.${MYVERSION}.ref
${component}_nml_mkref ${reffilename}
rpy.nml_get -v -f ${reffilename} > ${reffilename}.kv
cat ${reffilename}.kv | cut -d= -f1 > ${reffilename}.k
```

**1.2.5 Update nml_upd db**
```
component=rpnphy
MYVERSION0=${rpnphy_version}
MYVERSION=$(grep VERSION src/${component}/MANIFEST | cut -d: -f2)
MYVERSION=${MYVERSION# *}
cat >> src/${component}/share/nml_upd/${component}_nml_update_db.txt << EOF
#------
fileVersion: ${MYVERSION0} > ${MYVERSION}
$(diff src/${component}/share/nml_ref/${component}_settings.${MYVERSION0}.ref.k \
       src/${component}/share/nml_ref/${component}_settings.${MYVERSION}.ref.k \
       | egrep '(>|<)' \
       | sed 's/=/ = /g' | sed 's/>/#New: /' | sed 's/</rm: /')
EOF
```

**1.2.6 Update nml doc**
```
cd src
export PATH=$(pwd)/${component}/bin:${PATH}
${component}_ftnnml2wiki --comp ${component} --sort --wiki
${component}_ftnnml2wiki --comp ${component} --sort --md
mv ${component}.namelists.* ${component}/share/doc
cd ..
```

**1.2.7 Commit && Tag**
```
component=rpnphy
MYVERSION=$(grep VERSION src/${component}/MANIFEST | cut -d: -f2)
MYVERSION=${MYVERSION# *}
git add src/${component}/share/nml_ref/${component}_settings.${MYVERSION}.ref*
git commit -a -m "${component}: update VERSION, doc and nml ref (${component}_${MYVERSION})"
git tag gem_${ATM_MODEL_VERSION}+${component}_${MYVERSION}
```


2-Move patches from GEM dev to rpnphy depot
-------------------------------------------

> _All the following steps (in the section 2) are done in your rpnphy clone._

**2.1 Create a local clone of the rpnphy repo.**
```
git clone git@gitlab.science.gc.ca:MIG/rpnphy.git
cd rpnphy
```

**2.2 Make sure you are on the right branch & version**
```
git checkout rpnphy_6.3-branch
gitk
```
Visually check the your branhc (yellow dot) is "ahead"/"on top" of the rpnphy version/tag you are 
developping from (e.g. `rpnphy_6.3.0-b01`). This is original tag, accessible through `rpnphy_${rpnphy_version}` 
in the GEM experiment above.

For the first alpha you'll need to create a new branch, for example for the new 6.4.0-a01`:
```
git checkout -b rpnphy_6.4-branch
```

**2.3 Pull patches from the GEM dev experiment into the local repo.**  
_Set the path and tags in vars bellow_
```
component=rpnphy
gemurldir=     # URL or path of the repos where the new rpnphy code is (dev experiment above)
FROM_VERSION=  # Original GEM version number on top of which modifications
               # were made (`${ATM_MODEL_VERSION}` in gemdev repos)
gemtag=        # tag of the code to pull from in the repos where the new
               # rpnphy code is (above created
					# `gem_${ATM_MODEL_VERSION}+${component}_${MYVERSION}`)
```
_Pull the code from the GEM repos you "finalized" in step 1 above._
```
. r.load.dot eccc/mrd/rpn/MIG/ENV/modeltools/tools
t.gitpull --relative src/${component} -- ${gemurldir} ${gemtag} ${FROM_VERSION}
```

**2.4 Tag**
```
rpnphy_version=$(grep ^VERSION MANIFEST | cut -d: -f2)
rpnphy_version=${rpnphy_version# *}
phybranch=${rpnphy_version%.*}
phybranch=rpnphy_${phybranch:-6.3}-branch
git tag rpnphy_${rpnphy_version}
```

**2.5 Push the code and tag to the official repos**  
_Visually inspect the rpnphy repos (using `gitk`) to make 
sure everything is in order before pussing to the official repos._
```
git push origin ${phybranch} rpnphy_${rpnphy_version}
```

3-Documentation
---------------

* Update MIG/rpnphy issue and milestones on gitlab by closing them.
* Update rpnphy nml doc on the wiki
  * Source: rpnphy repos in file `share/doc/rpnphy.namelists.wiki`
  * Dest: http://wiki.cmc.ec.gc.ca/wiki/RPNPhy/6.2/settings_inc (or 6.3)



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
component=rpnphy
# git subtree split --prefix=src/${component}
mkdir -p patches/${component}
git format-patch -o patches/${component} \
    -M -C --find-copies-harder -k --relative=src/${component} \
    ${FROM}..HEAD
``` 

> You might want to do your modelutils update before leaving GEM env
> open `src/modelutils/README_librarian.md` and follow instructions

Steps to be done in the rpnphy clone
------------------------------------
### Import the patches
```
mydir=     #<- should be the ${gem_DIR} of the dev directory above
cd ..      #create the clone outside your gem experiment
git clone git@gitlab.science.gc.ca:MIG/rpnphy.git
cd rpnphy
rpnphy_version=$(grep ^VERSION MANIFEST | cut -d: -f2)
rpnphy_version=${rpnphy_version# *}
phybranch=${rpnphy_version%.*}
phybranch=rpnphy_${phybranch:-6.3}-branch
git checkout ${phybranch}
component=rpnphy
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
rpnphy_version=$(grep ^VERSION MANIFEST | cut -d: -f2)
rpnphy_version=${rpnphy_version# *}
phybranch=${rpnphy_version%.*}
phybranch=rpnphy_${phybranch:-6.3}-branch
git tag rpnphy_${rpnphy_version}
git push origin ${phybranch} rpnphy_${rpnphy_version}
```
