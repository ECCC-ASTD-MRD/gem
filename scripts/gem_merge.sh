#!/usr/bin/env bash

PREFIX="dp dm pm pp"
INPATH=`pwd`/RUNMOD
OUTPATH=`pwd`

script=$(basename "${BASH_SOURCE[0]}")
short="i:p:o:h"
long="prefix:,inpath:,outpath:,help"
opts=$(getopt -o $short --long $long --name "$script" -- "$@")
usage="\nMerge per PE output files\n
Usage : ${script} [-p prefix to merge] [-i input model run path] [-o output path]\n
   -p : list of prefix to merge (${PREFIX})
   -i : model run input path (${INPATH})
   -o : merged output path (${OUTPATH})\n
   example : ${script} -p \"pm dm\" -i ../RUNMOD -o \${TMPDIR}\n"

eval set -- "${opts}"

while :; do
    case "${1}" in
        -p | --prefix     ) PREFIX=$2;               shift 2 ;;
        -i | --inpath     ) INPATH=$2;               shift 2 ;;
        -o | --outpath    ) OUTPATH=$2;              shift 2 ;;
        -h | --help       ) echo -e "${usage}" 1>&2; exit ;;
        --                ) shift;                   break ;;
        *                 ) echo "error parsing";    exit 1 ;;
    esac
done

if [[ $(echo $(which editfst >/dev/null 2>&1 ; echo $?)) -gt 0 ]] ; then
   echo "(ERROR) editfst not found in PATH"
   exit 1
fi

# Get list of files
lst=`find ${INPATH} -name dp??????????-*_???`

# Get run date
run=`echo ${lst##*/}`
run=`echo ${run:2:10}`
echo "(INFO) Found run ${run}"

# Get time extensions
exts=`echo $lst | tr ' ' '\n' | awk -F '_' '{print $(NF)}' | sort -u`

for ext in $exts; do
   echo "(INFO) Processing extension ${ext}"
   rm ${OUTPATH}/${run}_${ext}.fst &>/dev/null
   for pre in ${PREFIX}; do
      echo "(INFO)   Processing prefix ${pre}"
      lst=`find ${INPATH} -name ${pre}${run}-*_${ext} | xargs`
      editfst -s $lst -d ${OUTPATH}/${run}_${ext}.fst -i 0 2>/dev/null
   done
done

