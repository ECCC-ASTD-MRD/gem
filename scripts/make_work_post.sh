#!/bin/bash -ex
#

if [ -n "${GEM_STORAGE_DIR}" -a -n "${GEM_ARCH}" -a  -n "${GEM_WORK}" ] ; then
   cd ${GEM_WORK}
   /bin/rm -f *.sh *.ksh *.py *.dot grille findtopo r.* vstudio
fi
