#!/bin/bash
#
#====> Obtaining the arguments:
arguments=$*
eval `cclargs_lite $0 \
     -progh    ""            ""     "[Source directory                   ]"\
     -type     ""            ""     "[dm pm dh ph dp pp file type        ]"\
     -assemble "0"           "1"    "[Reassemble or not                  ]"\
     -dplusp   "0"           "1"    "[Combine dynamics and physics output]"\
     -src      "input"   "input"    "[Source directory                   ]"\
     -dst      "output"  "output"   "[Destination directory for output   ]"\
     -flist    ""            ""     "[List of all model output files     ]"\
  ++ $arguments`

# Preliminary setup
set -ex
bliste=$(grep "/${type}.*_${progh}$" ${flist} | xargs)

# Determine number of filest to process
ienati=$(echo ${bliste} | wc -w)

if [ ${ienati} -gt 0 ] ; then

  # Common output header
  printf "    reassembleur $progh: START: $(date)\n"
  printf "    source               : ${bliste}\n"

  # Reassemble files
  if [ ${assemble} -gt 0 ] ; then

    for ii in ${bliste} ; do
      destination=${ii##*/}
      if [ ${dplusp} -gt 0 ] ; then
        destination=$(echo $destination | sed 's/\(^.\)\(.*\)/\2/')
      fi
      destination=${destination%%-*}_${destination#*_}
      break
    done
    printf "    destination          : ${destination}\n"

    fplis=" "
    for j in ${bliste} ; do
      fplis=${fplis}" "${src}/$j
    done
    bliste=${fplis}

    nfiles=$(echo ${bliste} | wc -w)
    echo "editfst -s ${bliste} -d ${dst}/${destination} #input files=$nfiles"
    editfst -s ${bliste} -d ${dst}/${destination} -i <<EOF
exclure (-1,[">>","^^","^>","!!"])
EOF
    editfst -s ${bliste} -d ${dst}/${destination} -e -i <<EOF
desire (-1,[">>","^^","^>","!!"])
EOF

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

