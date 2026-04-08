#!/bin/bash
#
#====> Obtaining the arguments:
arguments=$*
eval `cclargs_lite $0 \
     -progh    ""            ""     "[Source directory                    ]"\
     -type     ""            ""     "[dm pm dh ph dp pp usr[1-9] file type]"\
     -assemble "0"           "1"    "[Reassemble or not                   ]"\
     -dplusp   "0"           "1"    "[Combine dynamics and physics output ]"\
     -src      "input"   "input"    "[Source directory                    ]"\
     -dst      "output"  "output"   "[Destination directory for output    ]"\
     -flist    ""            ""     "[List of all model output files      ]"\
  ++ $arguments`

# Preliminary setup
set -ex

is_usr=0
if echo ${type} | grep -q usr[1-9];then
    # User output
    is_usr=1
    bliste=$(grep "/.*_${progh}$" ${flist} | xargs)
else
    # Regular model output
    bliste=$(grep "/${type}.*_${progh}$" ${flist} | xargs)
fi

# Determine number of filest to process
ienati=$(echo ${bliste} | wc -w)
echo 'exclure (-1,[">>","^^","^>","!!"])' > e1.dir
echo 'desire (-1,[">>","^^","^>","!!"])' > e2.dir

if [ ${ienati} -gt 0 ] ; then

  # Common output header
  printf "    reassembleur $progh: START: $(date)\n"
  printf "    source               : ${bliste}\n"

  # Reassemble files
  if [ ${assemble} -gt 0 ] ; then

     abort_file=assemble.abort_$$
     touch ${abort_file}
     for ii in ${bliste} ; do
        destination=${ii##*/}
        if [ ${dplusp} -gt 0 ] ; then
           destination=$(echo $destination | sed 's/\(^.\)\(.*\)/\2/')
        fi
        if [ ${is_usr} = 1 ];then
	        destination=${type: -1}${destination}
        fi
        destination=$(echo ${destination} | sed 's/-[0-9]*-[0-9]*//')
        break
     done
     printf "    destination          : ${destination}\n"

     fplis=" "
     for j in ${bliste} ; do
        fplis=${fplis}" "${src}/$j
     done
     bliste=${fplis}
     
     nfiles=$(echo ${bliste} | wc -w)
     del=80 ; upv=1 ; cur=1
     for k in $(seq 1 $del $nfiles) ; do
	     cur=$k
        upv=$((cur+del-1))
        upv=$((upv < nfiles ? upv : nfiles))
        subl=$(echo ${bliste} | cut -d " " -f ${cur}-${upv})
        editfst -s ${subl} -d ${dst}/${destination} -i e1.dir
        editfst -s ${subl} -d ${dst}/${destination} -e -i e2.dir
     done

     /bin/rm -f ${abort_file}
  else

    # No file reassembly
    printf "    destination          : ${dst}\n"

    mkdir -p ${dst}

    for i in ${bliste} ; do
      fn=$(basename $i)
      dyn=$(echo $fn | grep ^d | wc -l)
      phy=$(echo $fn | grep ^p | wc -l)
      if [ $dplusp -gt 0 ] ; then
        destination=$(echo $fn | sed 's/\(^.\)\(.*\)/\2/')
	     if [ $dyn -gt 0 ] ; then
	       cp ${src}/$i ${dst}/${destination}
	     else
	       editfst -e -s ${src}/${i} -d ${dst}/${destination} -i /dev/null
	     fi
      else
	     ln -sf $(readlink -e ${src})/$i ${dst}
      fi
    done

  fi

  printf "    Um_reassemble.sh: END: $(date)\n"

fi

