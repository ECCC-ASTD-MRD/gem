#!/bin/sh

OSNAME=
OSVERSION=unknown
unset OS

# remove white spaces
delSpaces() { echo $* | sed 's/ //g'; }

if [ -n "${ORDENV_PLAT}" ]; then
    OS="${ORDENV_PLAT}";

elif [ -f /etc/os-release ]; then
    . /etc/os-release
    OSNAME=`delSpaces ${NAME}`
    if [ -n "${VERSION_ID}" ]; then
        OSVERSION="${VERSION_ID}"
    elif [ -n "${BUILD_ID}" ]; then
        OSVERSION="${BUILD_ID}"
    fi
    OS="${OSNAME}-${OSVERSION}-`uname -m`";

elif [ -x /usr/bin/lsb_release ]; then
    OSNAME=`/usr/bin/lsb_release -is`;
    OSVERSION=`/usr/bin/lsb_release -rs`;
    OS="${OSNAME}-${OSVERSION}-`uname -m`";

elif [ `uname` = "Darwin" ]; then
    OSNAME="Darwin";
    OSVERSION=`uname -r`;
    OS="${OSNAME}-${OSVERSION}";

else 
    OS="Unknown-OS";
fi

printf "${OS}"
