#!/bin/bash

if [[ $# -lt 1 ]]
then
  echo "${0} Error: Please provide a full path to GEM-MACH II directory"
  exit 1
fi

cd ${1}
chk_c=`ls -1 | grep CHEM` ; chk_g=`ls -1 | grep GEM`

if [[ -z ${chk_c} || -z ${chk_g}  ]]
then
  echo "${0} Error: One of CHEM or GEM directories is not present " 
  echo "In ${1}"
  exit 1
fi

if  [[ "${TRUE_HOST}" == "hadar" ||  "${TRUE_HOST}" == "spica" ]]
then
   echo "Host platform is: ${TRUE_HOST} ... continuing ... "
else
   echo "${0} Error: Must execute this script on either hadar or spica (IBM P7 clusters)" 
   exit 1
fi

if [[ ! -z ${GEM_VER} ]]
then
   echo "GEM_VER define as: ${GEM_VER}"
   if [[ -z ${gemdyn}  ]]
   then
      . s.ssmuse.dot GEM/x/${GEM_VER}
   fi
else
   echo '${0} ERROR: GEM_VER environment variable is empty, please define it  (eg. 4.8.rc8)'
   exit 1
fi


cd ${1}/GEM
ouv_exp_gem
linkit

cd ${1}/CHEM
make

cd ${1}/GEM
make dep
make -j obj
make gemdm
