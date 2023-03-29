#!/bin/bash
#
arguments=$*

eval `cclargs_lite $0 \
  -inrep    ""             ""             "[input directory     ]"\
  -nthreads "1"            "1"            "[# of simultaneous threads]"\
  -rsync_cmd  "rsync"      ""             "[command used for rsync]" \
  -abortf   "upload_abort" "upload_abort" "[abort file          ]"\
  ++ $arguments`

printf "\n=====>  Um_upload.sh starts: `date` ###########\n\n"

abort_file=${TASK_WORK}/${abortf}_$$
touch ${abort_file}
set -ex
ici=$PWD

if [ -n "${inrep}" ] ; then

   bname=$(basename ${inrep})

   this_task=Um_upload_${bname}_$$
   workdir=${TASK_WORK}/${this_task}
   cd $inrep
   upload_list=$(find . -type l | sed 's-^\./--')
   mkdir -p $workdir ; cd $workdir

   cnt=0 ; abort_prefix=remote_copy_${this_task}_abort_

#   set +ex
   for file in ${upload_list} ; do
      printf "\n ===> Treating $file\n\n"
      target_base=${TASK_OUTPUT}/$(basename $inrep)
	   out_mach=${TRUE_HOST}
	   out_dst=${target_base}/${file}
	   mkdir -p ${target_base}/$(dirname ${file})
      cnt=$(( cnt + 1 ))
      lis=remote_copy_${this_task}_$(basename ${file}).lis
      ${TASK_BIN}/Um_upload_data.sh \
         -src $(readlink ${inrep}/${file}) \
	      -dst ${out_mach}:${out_dst} \
         -rsync_cmd ${rsync_cmd} \
         -abort ${abort_prefix}$(basename ${file}) \
         1> $lis 2>&1 &
      if [[ $cnt -eq $nthreads ]] ; then
        date ; wait ; date
        cnt=0
      fi
   done

   date ; wait ; date

# Check for aborted functions
   if [[ $(find . -name "${abort_prefix}*" | wc -l) -gt 0 ]] ; then
     echo "ERROR: One or more ${TASK_BIN}/Um_upload_data.sh function calls aborted ... see listings ${PWD}/remote_copy_${this_task}_*.lis for details"
     cat remote_copy_${this_task}_*.lis
     exit 1
   fi
fi

cd $ici ; /bin/rm -f ${abort_file}

printf "\n=====>  Upload.sh ends: `date` ###########\n\n"
