#!/bin/sh

COMP_VERSION="Unknown_Compiler"

case $1 in
    gfortran)
        COMP_VERSION=`gfortran --version | head -n 1 | sed "s/.*) \([^ ]*\).*$/\1/"`
        ;;
    intel)
        type ifort >/dev/null &&
            COMP_VERSION=`ifort -V 2>&1 | head -n 1 | sed "s/.*Version \([^ ]*\).*$/\1/"`
        ;;
    pgi)
        COMP_VERSION=`pgfortran --version  | sed -n "/fortran/ s/.*pgfortran \([^ ]*\).*$/\1/p"`
        ;;
    # Example: "Version: 16.01.0000.0000"
    xlf)
        COMP_VERSION=`xlf -qversion  | sed -n "/Version: / s/.*Version: \([0-9]*.[0-9]*\).*$/\1/p"`
        ;;
esac

printf "${COMP_VERSION}"
