#!/bin/sh

OSNAME=
OSVERSION=
OS=

# remove white spaces
sansespace() { echo $* | sed 's/ //g'; }

EC_BASE_ARCH=`echo $BASE_ARCH`

# use local variable at EC
if [ -n "${EC_BASE_ARCH}" ];
then
    OS="${EC_BASE_ARCH}";

elif [ -f /etc/os-release ];
then
    . /etc/os-release
    #remove white spaces
    #OSNAME=${NAME// /}
    OSNAME=`sansespace ${NAME}`
    OSVERSION="${VERSION_ID}"
    OS="${OSNAME}-${OSVERSION}-`uname -p`";

elif [ -x /usr/bin/lsb_release ];
    then
        OSNAME=`/usr/bin/lsb_release -is`;
        OSVERSION=`/usr/bin/lsb_release -rs`;
        OS="${OSNAME}-${OSVERSION}-`uname -p`";

elif [ `uname` = "Darwin" ];
    then
        OSNAME="Darwin";
        OSVERSION=`uname -r`;
        OS="${OSNAME}-${OSVERSION}";
else 
    OS="Unknown-OS";
fi

#OSNAME="${OSNAME%\"}"
#OSNAME="${OSNAME#\"}"

printf "${OS}"
