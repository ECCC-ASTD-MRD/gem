#!/bin/bash

set -x
# set -o noclobber is used to avoid writing in existing "content" files
#    This was needed for parallel execution of this script by the r1 prep_lbcs_loop task
set -o noclobber

rep_in=${1}
wild=${2}
rep_ou=${3}
nthreads=${4}

model_inputs ()
{
set -x
date
file=$1
rep_ou=$2
lnk_dest=$3
bname=$(basename ${file})
if [ "${lnk_dest}" != "${bname}" ] ; then lnk_dest=${lnk_dest}/${bname} ; fi
list_A=$(r.fstliste -izfst $file -typvar "A" |cut -d ":" -f 11 | sort -u )
list_P=$(r.fstliste -izfst $file -typvar "P" |cut -d ":" -f 11 | sort -u )
list_I=$(r.fstliste -izfst $file -typvar "I" |cut -d ":" -f 11 | sort -u )
for i in $(echo "${list_A} ${list_P} ${list_I}" | tr ' ' '\n' | sort -u | tr '\n' ' ') ; do
  valid=$(echo $i | cut -c1-8).$(echo $i | cut -c9-14)
  dir=${rep_ou}/VALID_${valid}
  mkdir -p ${dir}
  ln -sf ../../${lnk_dest} ${dir}/GEM_input_file-${bname}
done
}

if [ "${wild}" == "@NIL@" ] ; then
  wild=
fi

if [ -d $rep_in ] ; then

  cd $rep_in
  count=0
  for file in $(ls -1 ${wild}) ; do
    if [ ! -d $file ] ; then
      count=$(( count + 1 ))
      model_inputs $file ${rep_ou} MODEL_inrep &
    fi
    if [ $count -eq ${nthreads} ] ; then
      wait
      count=0
    fi
  done
  wait

else
  model_inputs $rep_in ${rep_ou} ANALYSIS
fi

for i in ${rep_ou}/VALID_*  ; do
  if [ -d $i ] ; then
    cd $i
    cnt=$(ls -1 GEM_input_file* | wc -l)
    set +e
    cat > content <<EOF
$(echo $cnt)
$(ls -1 GEM_input_file*)
EOF
    set -e
  fi
done
