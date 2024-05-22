Updating the modelutils depot for a GEM release
===========================================

# Steps to be done in a GEM dev env.

... include code from contrubutors & test ...
... see with RPN-SI if there are updates needed to rpnphy's CMakeLists.txt ...

Tests
=====

Make sure to test with GFortran and intel

1st shell
```
. .ssmuse_gem intel
. .intial_setup
make cmake
make -j4
make work
# ... run tests...
```

2nd shell
```
. .ssmuse_gem gnu
. .intial_setup
make cmake
make -j4
make work
# ... run tests...
```

Finalize
========

# Steps to be done in a GEM dev env.

# Update MANIFEST for VERSION
```
emacs src/modelutils/MANIFEST
```

Move patch from GEM dev to modelutils depot
=======================================

## Steps to be done in a GEM dev env.

# WARNING: need to "squash the merges" beforehand
# Check if there are merge commits
```
FROM=${ATM_MODEL_VERSION:-__scrap__}
git rev-list --min-parents=2 --count ${FROM}..HEAD
```

```
FROM=${ATM_MODEL_VERSION:-__scrap__}
component=modelutils
# git subtree split --prefix=src/${component}
mkdir -p patches/${component}
git format-patch -o patches/${component} \
    -M -C --find-copies-harder -k --relative=src/${component} \
    ${FROM}..HEAD
``` 

## Steps to be done in the modelutils clone

```
mydir=
component=modelutils
# patchdir=${gemdir:-no-such-dir}/patches/modelutils
# git am --empty=drop -- ${patchdir}/*.patch
for item in ${mydir:-no-such-dir}/patches/${component}/*.patch ; do
    if [[ -s ${item} ]] ; then
       git am -3 -k ${item}
       if [[ $? != 0 ]] ; then
          echo "ERROR: problem applying ${item}"
          break
       fi
    else
       echo "Ignoring: ${item}"
    fi
done
```

# commit & tag & push


Documentation
=============

# Udpate MIG/modelutils issue and milestones

# Update modelutils doc
# Dest: http://wiki.cmc.ec.gc.ca/wiki/modelutils
