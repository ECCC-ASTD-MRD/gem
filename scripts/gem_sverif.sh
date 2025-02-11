#!/usr/bin/env bash

#CONFIG=GEM_cfgs_GY_4km
MD5=85941a37c3fdfb6846c9f97d72222a06
URL=http://collaboration.cmc.ec.gc.ca/science/outgoing/sverif

script=$(basename "${BASH_SOURCE[0]}")
short="t:p:h:f:"
long="type:,runpath:,file:,help"
opts=$(getopt -o $short --long $long --name "$script" -- "$@")
usage="\nValidate GEM benchmark run\n
Usage : ${script} [-c model run configuration] [-p input model run path]\n
   -c : configuration
   -p : model run input path\n
   example : ${script} -c \"GEM_cfgs_GY_4km\" -p ../RUNMOD\n"

eval set -- "${opts}"

while :; do
    case "${1}" in
        -c | --config     ) CONFIG=$2;               shift 2 ;;
        -p | --runpath    ) RUNPATH=$2;              shift 2 ;;
        -f | --file       ) FILE=$2;                 shift 2 ;;
        -h | --help       ) echo -e "${usage}" 1>&2; exit ;;
        --                ) shift;                   break ;;
        *                 ) echo "error parsing";    exit 1 ;;
    esac
done

if [[ -z "${RUNPATH}" ]]; then
    echo "(ERROR) No model run path specified"
    exit 1
fi

# Check for required environment
if [[ $(echo $(which sverif_eval.Abs >/dev/null 2>&1 ; echo $?)) -gt 0 ]] ; then
    echo "(ERROR) sverif_eval.Abs not found in PATH"
    exit 1
fi

# Check for config in work dir
if [[ -z "${CONFIG}" ]]; then
   eval `grep dircfg= ${RUNPATH}/RUNMOD/.setup/task_setup_set.txt`
   eval `grep ATM_MODEL_VERSION= ${RUNPATH}/RUNMOD/.setup/task_setup_set.txt`
   CONFIG=${dircfg##*/}-${ATM_MODEL_VERSION}
   echo "(INFO) Run configuration found (${CONFIG})"
fi

if [[ -z "${CONFIG}" ]]; then
    echo "(ERROR) No run configuration found"
    exit 1
fi

# Verification data if not yet done
statpath=${RUNPATH}/sverif/${CONFIG}

if [[ ! -d ${statpath} ]]; then
   
   echo "(INFO) Downloading validation data to $statpath"
   wget ${URL}/${CONFIG}.tgz -O ${RUNPATH}/${CONFIG}.tgz

   md5=$(md5sum ${RUNPATH}/${CONFIG}.tgz | cut -d' ' -f1)
   if [[ "${md5}" != "${MD5}" ]]; then
      echo "(ERROR) Download failed"
      exit 1
   fi

   tar -zxf ${RUNPATH}/${CONFIG}.tgz -C ${RUNPATH}
fi

# Get forecast hour information
prog=$(find ${statpath}/ -name 'sverif*.dat' | head -n 1 | sed 's/h.dat//' | tail -c 4)
file=$(find ${RUNPATH}/RUNMOD/output -name ${FILE})

if [[ ! -d ${statpath} ]]; then 
    echo "(ERROR) Invalid statistics path (${statpath})"
    exit 1
fi

# Run verification
tmpfile=/tmp/sverif$$
lev=500
# Make sure the stack limit is set to unlimited
ulimit -s unlimited
echo "sverif_eval.Abs GZ ${lev} ${prog} ${file} ${statpath}"
sverif_eval.Abs GZ ${lev} ${prog} ${file} ${statpath} | grep CI=0.01 | grep -v '^\*' >${tmpfile}
hits=$(cat ${tmpfile} | grep 'PASS' | grep -v 'overall' | wc -l)
cat ${tmpfile}

lev=850
sverif_eval.Abs TT ${lev} ${prog} ${file} ${statpath} | grep CI=0.01 | grep -v '^\*' >${tmpfile}
hits=$(( ${hits} + $(cat ${tmpfile} | grep 'PASS' | grep -v 'overall' | wc -l) ))
cat ${tmpfile}; rm -f ${tmpfile}

if [[ ${hits} -ge 6 ]] ; then
    echo "(INFO) Passed (passed ${hits}/8)"
else
    echo "(INFO) Failed (passed ${hits}/8)"
fi
