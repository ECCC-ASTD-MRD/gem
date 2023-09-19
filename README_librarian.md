Updating the modelutils depot for a GEM release
===========================================

# Steps to be done in a GEM dev env.

... include code from contrubutors & test ...


Finalize
========

# Steps to be done in a GEM dev env.

# Update MANIFEST for VERSION
```
emacs src/modelutils/MANIFEST
```

Move patch from GEM dev to modelutils depot
=======================================

# Steps to be done in a GEM dev env.

```
# git subtree split ...
git format-patch FROM..HEAD
```


# Steps to be done in the modelutils clone

```
git am *.patch
```

# commit & tag & push


Documentation
=============

# Udpate MIG/modelutils issue and milestones

# Update modelutils doc
# Dest: http://wiki.cmc.ec.gc.ca/wiki/modelutils
