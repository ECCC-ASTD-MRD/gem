#!/bin/bash

# script à revoir selon présence storage_model

target=$1

case ${target} in
    "build")
        echo /bin/rm -rf ${GEM_DIR}/build;
        #/bin/rm -rf ${GEM_DIR}/build;
        echo /bin/rm -rf ${GEM_DIR}/build-${GEM_ARCH}
        #/bin/rm -rf ${GEM_DIR}/build-${GEM_ARCH}

        ;;
    "work")
        echo /bin/rm -rf ${GEM_DIR}/work;
        #/bin/rm -rf ${GEM_DIR}/WORK;
        echo /bin/rm -rf ${GEM_DIR}/work-${GEM_ARCH}
        #/bin/rm -rf ${GEM_DIR}/work/work-${GEM_ARCH}
        ;;
    "all")
        #à finaliser!
        echo /bin/rm -rf ${GEM_DIR}/work-${GEM_ARCH}
        echo /bin/rm -rf ${GEM_DIR}/build-${GEM_ARCH}
        ;;
    *)
        echo "please specify a target: build, work or all";
        ;;
esac
