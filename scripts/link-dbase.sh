#!/bin/bash

# link to GEM database at CMC
GEM_DBASE=/fs/ssm/eccc/mrd/rpn/MIG/GEM/d/gem-data/gem-data_4.2.0/gem-data_4.2.0_all/share/data/dfiles

if [ ! -d ${GEM_DBASE} ] ; then
    echo "${GEM_DBASE} not found: don't know where database is."
    exit 1
fi

# remove possible existing link or file, and create symbolic link, 
# or display an error message (for example if a directory named gem_dbase already exists)
\rm -f gem_dbase && ln -s ${GEM_DBASE} gem_dbase || ( echo "gem_dbase cannot be removed" && exit 2)
