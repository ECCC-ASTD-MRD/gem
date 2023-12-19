
Updating the rpnphy depot for a GEM release
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
emacs src/rpnphy/MANIFEST
```

# Create nml_ref
```
cado cmake
cado build
cado work
component=rpnphy
MYVERSION=$(grep VERSION src/${component}/MANIFEST | cut -d: -f2)
MYVERSION=${MYVERSION# *}
reffilename=src/${component}/share/nml_ref/${component}_settings.${MYVERSION}.ref
${component}_nml_mkref ${reffilename}
rpy.nml_get -v -f ${reffilename} > ${reffilename}.kv
cat ${reffilename}.kv | cut -d= -f1 > ${reffilename}.k
```

# Update nml_upd db
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

# Update nml doc
```
cd src
export PATH=$(pwd)/${component}/bin:${PATH}
${component}_ftnnml2wiki --comp ${component} --sort --wiki
${component}_ftnnml2wiki --comp ${component} --sort --md
mv ${component}.namelists.* ${component}/share/doc
cd ..
```

# Commit
```
component=rpnphy
MYVERSION=$(grep VERSION src/${component}/MANIFEST | cut -d: -f2)
MYVERSION=${MYVERSION# *}
git add src/${component}/share/nml_ref/${component}_settings.${MYVERSION}.ref*
git commit -a -m "${component}: update VERSION, doc and nml ref (${component}_${MYVERSION})"
```


Move patch from GEM dev to rpnphy depot
=======================================

# Steps to be done in a GEM dev env.

# WARNING: need to "squash the merges" beforehand
# Checck if there are merge commits
```
FROM=${ATM_MODEL_VERSION:-__scrap__}
git rev-list --min-parents=2 --count ${FROM}..HEAD
```

```
FROM=${ATM_MODEL_VERSION:-__scrap__}
component=rpnphy
# git subtree split --prefix=src/${component}
mkdir -p patches/${component}
git format-patch -o patches/${component} \
    -M -C --find-copies-harder -k --relative=src/${component} \
    ${FROM}..HEAD
``` 


# Steps to be done in the rpnphy clone

```
mydir=
component=rpnphy
for item in ${mydir}/patches/${component}/*.patch ; do
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

# Udpate MIG/rpnphy issue and milestones

# Update rpnphy nml doc
# Source: `share/doc/rpnphy.namelists.wiki`
# Dest: http://wiki.cmc.ec.gc.ca/wiki/RPNPhy/6.2/settings_inc (or 6.3)
