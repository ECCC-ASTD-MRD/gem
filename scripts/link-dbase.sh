#!/bin/bash

# link to database at CMC

GEM_DBASE=/fs/ssm/eccc/mrd/rpn/MIG/GEM/d/gem-data/gem-data_4.2.0/gem-data_4.2.0_all/share/data/dfiles

if [ -e ${GEM_DBASE} ] &&  [ ! -e gem_dbase ] ; then
    ln -sf ${GEM_DBASE} gem_dbase
else
    echo "Database not found: don't know where database is."
fi
