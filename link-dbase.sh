#!/bin/bash

# link to database at CMC

if [[ "$(hostname -A)" =~ cmc.ec.gc.ca ]]; then
    GEM_DBASE=/home/ordenv/ssm-domains9/release/gem-data_4.2.0/gem-data_4.2.0_all/share/data/dfiles
    [ -e ${GEM_DBASE} ] && [ ! -e gem_dbase ] && ln -sf ${GEM_DBASE} gem_dbase
elif [[ "$(hostname -A)" =~ science.gc.ca ]]; then
    GEM_DBASE=/fs/ssm/eccc/mrd/rpn/MIG/GEM/d/gem-data/gem-data_4.2.0/gem-data_4.2.0_all/share/data/dfiles
    [ -e ${GEM_DBASE} ] && [ ! -e gem_dbase ] && ln -sf ${GEM_DBASE} gem_dbase
else
    echo "hostname not found: don't know where database is."
fi
