#!/bin/sh

GOASDBASE="gem_dbase.tar.gz"
GOASSITE="http://collaboration.cmc.ec.gc.ca/science/outgoing/goas/gem_dbase.tar.gz"

if [ -x /usr/bin/wget ];
then
    /usr/bin/wget ${GOASSITE};
elif [ -x /usr/bin/curl ];
then
    /usr/bin/curl -O ${GOASSITE}
else
    echo "Error: cannot download database using wget or curl."
    echo "Please download database at: ${GOASSITE}" 
    exit 1
fi

if [ -f ${GOASDBASE} ];
then
    tar -xzvf  ${GOASDBASE}
else
    echo "Error. Cannot find database tar in this directory."
fi
