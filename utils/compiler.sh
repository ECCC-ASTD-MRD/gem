#!/bin/sh

COMPILER_VERSION="Unknown_Compiler"

case $1 in
    aocc)
        COMPILER_VERSION=$(clang --version | grep version | sed "s/.*version \([^ ]*\).*$/\1/")
        ;;
    gnu)
        COMPILER_VERSION=$(gfortran --version | head -n 1 | sed "s/.*) \([^ ]*\).*$/\1/")
        ;;
    intel)
        # 1. We have to redirect to &1 as ifx -V and ifort -V output is done on stderr!
        # 2. We don't want an error message if ifx is not found, while ifort is found
        #    So, first, check for ifx and if not found, then check for ifort.
        _VERSIONx=$(ifx -V 2>&1)
        if [ $? == 0 ]; then # on a trouvé ifx
            COMPILER_VERSION=$(echo ${_VERSIONx} | head -n 1 | sed "s/.*Version \([^ ]*\).*$/\1/")
        else
            _VERSIONo=$(ifort -V 2>&1)
            if [ $? == 0 ]; then # on a trouvé ifort
                COMPILER_VERSION=$(echo ${_VERSIONo} | head -n 1 | sed "s/.*Version \([^ ]*\).*$/\1/")
            else  # on n'a ni ifx, ni ifort
                ( echo ${_VERSIONx}; echo ${_VERSIONo} ) >&2
            fi
        fi
        ;;
    llvm)
        COMPILER_VERSION=$(clang --version | grep version | sed "s/.*version \([^ ]*\).*$/\1/")
        ;;
    nvhpc)
        COMPILER_VERSION=$(nvcc --version | grep release | sed "s/.*V\([^ ]*\).*$/\1/")
        ;;
    pgi)
        COMPILER_VERSION=$(pgfortran --version  | sed -n "/fortran/ s/.*pgfortran \([^ ]*\).*$/\1/p")
        ;;
    # Example: "Version: 16.01.0000.0000"
    xlf)
        COMPILER_VERSION=$(xlf -qversion  | sed -n "/Version: / s/.*Version: \([0-9]*.[0-9]*\).*$/\1/p")
        ;;
esac

printf "${COMPILER_VERSION}"
