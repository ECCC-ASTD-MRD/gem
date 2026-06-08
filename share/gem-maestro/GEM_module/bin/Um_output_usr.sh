#!/bin/bash
#
arguments=$*
. r.entry.dot

#====> Obtaining the arguments:
eval `cclargs_lite $0 \
     -src      ""       ""       "[Source directory                   ]"\
     -type     ""       ""       "[type of output must be usr*        ]"\
     -dst      "output" "output" "[Destination directory for output   ]"\
     -_nf2t    "0"      "0"      "[# of files treated                 ]"\
     -_nerr    "0"      "0"      "[# of errors detected               ]"\
     -nthreads "1"       "1"     "[Number of bemol process to run in parallel]"\
  ++ $arguments`

set -ex
if [ -z "$src" -o -z "$type" ] ; then
  printf "\n ##### ABORT in Um_output_usr.sh: -src UNDEFINED #####\n\n"
  . r.return.dot
  exit 1
fi

printf "\n    =====> `basename $0` $arguments\n"

# Obtain list of files to process
here=${PWD}
cd ${src}
flist=${here}/file_list_${type}
find -L ./${type} -type f >${flist}
cd ${here}

ls -1 $src/${type} | grep '[0-9][0-9][0-9]*' > dir_list$$
rep_search=$(head -n 1 dir_list$$)
nrep=$(cat dir_list$$ | wc -l)
rm dir_list$$

find -L ${src}/${type}/${rep_search} -type f -name "[0-9]*" > files_found_${type}
nfiles=0
if [ -s files_found_${type} ] ; then
   nfiles=$(cat files_found_${type} | wc -l)
fi
nfiles=$((${nfiles}*nrep))
echo ${nfiles} > file_count_${type}
_nerr=0

if [ ${nfiles} -gt 0 ] ; then

   echo "Building list of TIMEFRAME for ${type} ..."
   dliste=""
   for i in $(cat files_found_${type}) ; do
       fname=${i##*/}
       step=${fname#*_}
       if [[ $(echo ${dliste} | grep -- ^${step} | wc -l) -lt 1 ]] ; then
           dliste="${dliste} ${step}"
       fi
   done
   cnt=0

   for i in ${dliste} ; do
       
       cnt=$(( cnt + 1 ))
       lis=assemble_${type}_${i}
   
       ${TASK_BIN}/Um_reassemble.sh -src ${src} -dst ${dst}\
      	          -progh =$i  -type ${type}\
                  -assemble ${assemble} -flist ${flist} > ${lis}.lis 2> ${lis}.err &
       if [[ $cnt -eq $nthreads ]] ; then
           date ; wait ; date
           cnt=0
       fi
   done
   date ; wait ; date          
fi

cnt=$(head -n 1 file_count_${type})
_nf2t=$((${_nf2t} + ${cnt}))
_nerr=0

. r.return.dot

exit 0

