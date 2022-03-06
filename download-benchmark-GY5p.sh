#!/bin/bash

set -e

BENCHMARK_GY_5="GY5kmP_dbase.tar.gz"
BENCHMARK_GY_5_URL="http://collaboration.cmc.ec.gc.ca/science/outgoing/goas/${BENCHMARK_GY_5}"
BENCHMARK_GY_5_MD5SUM="22c05a81244ce4b46cb9072b8c179701"

printUsage() {
    echo -e "Download a sample database of data files needed to run benchmarks for GEM"
    echo -e "Usage:"
    echo -e "./$(basename $0) <GEM-GIT-DIR>\n"
    echo -e "Usually, GEM-GIT-DIR is the current directory, so use:"
    echo -e "./$(basename $0) ."
}

checkMd5() {
    # $1 File path
    # $2 Expected MD5
    # Return: 0 if matching; 1 otherwise
    md5=$(md5sum "$1" | cut -d' ' -f1)
    [[ "$md5" = "$2" ]]
}

if [[ ! -d "$1" ]]; then
    printUsage
    exit 1
fi

if [[ $# -eq 2 ]]; then
    if [[ ! -r "$2" ]]; then
        echo "Can't read ${2} !"
        exit 1
    else
        tarballPath="$2"
    fi
else
    tarballPath="${1}/${BENCHMARK_GY_5}"
    if [[ -x $(which wget) ]]; then
        wget ${BENCHMARK_GY_5_URL} -O "${tarballPath}";
    elif [[ -x $(which curl) ]]; then
        curl -o "${tarballPath}" ${BENCHMARK_GY_5_URL}
    else
        echo "Error: cannot download using wget or curl."
        echo "Please download database at: ${BENCHMARK_GY_5_URL}" 
        exit 1
    fi
fi

if checkMd5 "$tarballPath" "$BENCHMARK_GY_5_MD5SUM"; then
    echo "MD5 check OK"
else
    echo "The MD5 of $2 does not match what was expected.  The file might be corrupted."
    exit 1
fi

tar -xzvf ${tarballPath} -C "$1"
