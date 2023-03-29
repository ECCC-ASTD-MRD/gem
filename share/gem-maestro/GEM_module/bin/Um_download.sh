#!/bin/bash
#
arguments=$*
eval `cclargs_lite $0 \
  -src       ""	  ""         "[]"\
  -dst_mach   ""	  ""         "[]"\
  -dst_dir   ""	  ""         "[]"\
  -abortf    "Um_output" "Um_output" "[abort file     ]"\
  -scp_cmd   "sscp"  ""        "[remote copy command ]"\
  -scp_opt   "=-r"   ""        "[remote copy option  ]"\
  ++ $arguments`
#  -domain    ""	  ""         "[]"\

printf "\n=====>  Um_download.sh starts: `date` ###########\n\n"

abort_file=${TASK_WORK}/${abortf}_$$
touch ${abort_file}
set -ex

#  CHECKSUM=$(which Um_checksum.sh)
#  src_md5=source_md5_${DOMAIN}_$(basename ${src_file}).lis
#  dst_md5=destination_md5_${DOMAIN}_$(basename ${src_file}).lis
#  ssh $src_mach "${CHECKSUM} ${src_file}" > ${src_md5}

#SRC=$(r.read_link ${src})
SRC=$(readlink ${src})
dst=${dst_dir}/$(basename ${src})

if [ "${dst_mach}" == "${TRUE_HOST}" ] ; then
  # printf "download: dst_mach is TRUE_HOST\n"
  # ls -l ${SRC} 
  # ls -l ${dst}
  if [[ -d ${dst} && x${SRC} != x${dst} ]] ; then
     /bin/rm -rf ${dst}
  fi
  if [[ -L ${dst} ]] ; then
     /bin/rm -f ${dst}
  fi
  ln -s ${SRC} ${dst}
elif [[ -d ${src} && ${MOD_GEM_ln_optimiz:-0} == 1 ]] ; then
  # printf "download: dst_mach can see src\n"
  # ls -l ${SRC} 
  # ls -l ${dst}
  if [[ -d ${dst} && x${SRC} != x${dst} ]] ; then
     /bin/rm -rf ${dst}
  fi
  if [[ -L ${dst} ]] ; then
     /bin/rm -f ${dst}
  fi
  ln -s ${SRC} ${dst}
else
  #printf "download: dst_mach cannot see src\n"
  # We here take advantage of the fact that ppp5-6:sitestore are mounted on robert and underhill
  machinelist=":robert:underhill:ppp5:ppp6:"
  is_visible() {
     if [[ -n "$(echo ${machinelist} | grep :${1}:)" ]] ; then
        echo 1
     else
        echo 0
     fi
  }
  if [[ $(is_visible ${dst_mach}) == 1 ]] ; then
    /bin/rm -rf ${dst}
    ${scp_cmd} ${scp_opt} ${SRC} ${dst}
  else
    #TODO: do we need sscp on U2 systems?
    ${scp_cmd} ${scp_opt} ${SRC} ${dst_mach}:${dst}
  fi
fi

#  Um_checksum.sh ${TASK_WORK}/${DOMAIN}/last* > ${dst_md5}
#  cnt=$(diff ${src_md5} ${dst_md5} | wc -l)
#  if [ ${cnt} -gt 0 ] ; then
#     printf "\n Problem with data transfer from ${src_mach} - ABORT\n\n"
#     exit 1
#  fi
#  rep_in=${TASK_WORK}/${DOMAIN}/last*

/bin/rm -f ${abort_file}

printf "\n=====>  Um_download.sh ends: `date` ###########\n\n"



