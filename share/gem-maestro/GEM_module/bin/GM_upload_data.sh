#!/bin/bash
#
arguments=$*
eval `cclargs_lite $0 \
   -src   ""       ""    "[ ]" \
   -dst   ""       ""    "[ ]" \
   -fcp   "0"      "1"   "[ ]" \
   -rsync_cmd  "rsync"  ""  "[command used for rsync]" \
   -abort ""       ""    "[Abort signal file prefix]"\
 ++ $arguments`

set -ex
this_process=$$
if [[ -n "${abort}" ]] ; then 
   abort_file=${abort}-${this_process}
   touch ${abort_file}
fi

src_mach=$(echo $src | awk '{print $1}')
src_file=$(echo $src | awk '{print $2}')
dst_mach=$(echo $dst | awk '{print $1}')
dst_file=$(echo $dst | awk '{print $2}')

if [ -z "$src" -o -z "$dst_mach" -o -z "$dst_file" ] ; then
  exit 0
fi

printf "\n$(basename $0) $arguments\n"

if [ -z "$src_file" ] ; then
  src_file=$src_mach
  src_mach=$dst_mach
fi 

cnt=$(echo $src_file | grep @@ | wc -l)
if [ $cnt == 1 ] ; then
  opt=${src_file##*@@}
  src_file=${src_file%%@@*}
  fcp=1
  if [ "${opt}" == "lnk" ] ; then fcp=0 ; fi
fi

if [ ${src_mach} == ${dst_mach} -a $fcp == 0 ] ; then

  if [ ${src_file} == ${dst_file} ] ; then
    echo "Attempt to create a circular link for ${src_file} in $0 ... exiting"
    exit 1
  fi
  if [ ${dst_mach} == ${TRUE_HOST} ] ; then
     ln -sf $src_file ${dst_file}
  else
     ssh ${dst_mach} "ln -sf $src_file ${dst_file}"
  fi

else

   # We here take advantage of the fact that ppp5-6:sitestore are mounted on robert and underhill
   machinelist=":robert:underhill:ppp5:ppp6:"
   is_visible() {
      if [[ -n "$(echo ${machinelist} | grep :${1}:)" ]] ; then
         echo 1
      else
         echo 0
      fi
   }
   HOST_src="${src_mach}:"
   if [[ $(is_visible ${src_mach}) == 1 ]] ; then
      unset HOST_src
   fi

  # Add a trailing '/' extension to copy directories and check for target type match for directory source

  if [ -z "${HOST_src}" ] ; then
     ext=$(if [[ -d $src_file ]] ; then echo '/' ; fi)
  else
     ext=$(ssh $src_mach "if [[ -d $src_file ]] ; then echo '/' ; fi")
  fi
  if [ ${ext}x == '/x' ] ; then
     ssh ${dst_mach} "[[ -d ${dst_file} ]] || rm -f ${dst_file}"
  fi

  local_lis=$(basename ${src_file})_${this_process}
  rsync_status=rsync_${local_lis}.lis

  ${rsync_cmd} --copy-unsafe-links --delete-before --chmod=u+w -ruvHl ${HOST_src}${src_file}${ext} ${dst_file} ; echo STATUS=\$? 1> ${rsync_status} 2>&1

  printf "\n DESTINATION: ${dst_mach}:${dst_file}\n"
  cat ${rsync_status}
  eval $(grep "STATUS=" ${rsync_status})
  printf "Final rsync exit status: $STATUS\n"
  if [ $STATUS != 0 ] ; then exit $STATUS ; fi

fi

if [[ -n "${abort_file}" ]] ; then rm -f ${abort_file} ; fi

