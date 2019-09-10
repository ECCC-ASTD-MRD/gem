#!/bin/bash
GEM_DBASE=/fs/ssm/eccc/mrd/rpn/MIG/GEM/d/gem-data/gem-data_4.2.0/gem-data_4.2.0_all/share/data/dfiles;
[ -e ${GEM_DBASE} ] && [ ! -e gem_dbase ] && ln -sf ${GEM_DBASE} gem_dbase;
