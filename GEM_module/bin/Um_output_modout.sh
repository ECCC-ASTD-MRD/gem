#!/bin/bash
#
arguments=$*
. r.entry.dot

#====> Obtaining the arguments:
eval `cclargs_lite $0 \
     -src      ""       ""       "[Source directory                   ]"\
     -dst      "output" "output" "[Destination directory for output   ]"\
     -liste    ""       ""       "[list of file types to process      ]"\
     -assemble "0"      "1"      "[Reassemble or not                  ]"\
     -dplusp   "0"      "1"      "[Combine dynamics and physics output]"\
     -_nf2t    "0"      "0"      "[# of files treated                 ]"\
     -_nerr    "0"      "0"      "[# of errors detected               ]"\
     -nthreads "1"       "1"     "[Number of bemol process to run in parallel]"\
  ++ $arguments`

set -ex
if [ -z "$src" ] ; then
  printf "\n ##### ABORT in Um_output_modout.sh: -src UNDEFINED #####\n\n"
  . r.return.dot
  exit 1
fi

printf "\n    =====> `basename $0` $arguments\n"

process_modout(){
  set -ex
  type=${1}
  ls -1 $src/ | grep [0-9]*-[0-9]* > dir_list$$
  rep_search=$(head -n 1 dir_list$$)
  nrep=$(cat dir_list$$ | wc -l)
  rm dir_list$$

  find -L ${src}/${rep_search}/ -name "${type}[0-9]*" > files_found_${type}
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
    echo "TIMEFRAME (${type}): ${dliste}"

    cnt=0
    for i in ${dliste} ; do
      cnt=$(( cnt + 1 ))

      lis=assemble_${type}_${i}
      ${TASK_BIN}/Um_reassemble.sh -src ${src} -dst ${dst}     \
	                                 -progh =$i  -type ${type}   \
                          -assemble ${assemble} -flist ${flist} \
	                       -dplusp ${dplusp} 1> ${lis}.lis 2> ${lis}.err &
      if [[ $cnt -eq $nthreads ]] ; then
         date ; wait ; date
         cnt=0
      fi
    done
    date ; wait ; date
  fi
}

# Obtain list of files to process
here=${PWD}
flist=${here}/file_list
cd ${src}
find -L ./ >${flist}
cd ${here}

# Loop over file types to process
touch .thread_init
for ftype in ${liste} ; do
  printf "\n  ==> PROCESSING ${ftype} files \n\n"
  process_modout ${ftype}
done
date ; wait ; date

for ftype in ${liste} ; do
  cnt=$(head -n 1 file_count_${ftype})
  _nf2t=$((${_nf2t} + ${cnt}))
done
_nerr=0

. r.return.dot

exit 0

