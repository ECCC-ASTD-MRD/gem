#!/bin/bash

arguments=$*
echo $0 $arguments
date
eval `cclargs_lite -D " " $0 \
  -bindir       ""  ""   "[path to namelist        ]"\
  -destination  ""  ""  "[to perform check or not ]"\
  -scp_cmd      "sscp"  ""  "[remote copy command ]"\
  -scp_opt      "=-r"   ""  "[remote copy option ]"\
  ++ $arguments`

set -ex

if [ -z "${bindir}"      ] ; then exit 0 ; fi
if [ -z "${destination}" ] ; then exit 0 ; fi

if [ "${bindir}" == "release" -o \
     "${bindir}" == "RELEASE" ] ; then exit 0 ; fi

dest_dir=${destination}
dest_mach=${TRUE_HOST}

ici=${PWD}
bin_mach=${bindir%%:*}
bin_dir=${bindir##*:}

if [ "${bin_mach}" == "${bin_dir}" ] ; then
  cd $bin_dir
  mkdir -p ${dest_dir}  
  find . -type f -exec ln -s ${PWD}/{} ${dest_dir} \;
  for elem in $(find . -type l) ; do
    ls -l ${elem}
    src_bin=$(readlink ${elem})
	 bin_mach=${src_bin%%:*}
	 bin_dir=${src_bin##*:}
	 if [ "${bin_mach}" == "${bin_dir}" ] ; then
      ln -s ${src_bin} ${dest_dir}/${elem}
    else
      if [ "${bin_mach}" == "${TRUE_HOST}" ] ; then
        rep=$(true_path ${bin_dir})
        find ${rep} -type f -exec ln -s {} ${dest_dir} \;
      else
        ${scp_cmd} ${scp_opt} ${src_bin}/* ${dest_dir}
      fi
    fi
  done
else
  ssh ${dest_mach} mkdir ${dest_dir}
  ${scp_cmd} ${scp_opt} ${bin_mach}:${bin_dir}/*/*.Abs ${dest_mach}:${dest_dir} 2> /dev/null || true
  ${scp_cmd} ${scp_opt} ${bin_mach}:${bin_dir}/*.Abs   ${dest_mach}:${dest_dir} 2> /dev/null || true
fi
cd ${ici}

