#!/bin/bash
#
arguments=$*
eval `cclargs_lite $0 \
  -domain    ""	""         "[input/output domain                          ]"\
  -xcasc     ""   ""         "[backend directory to link cascade files      ]"\
  -d2z       "0"  "1"        "[assemble all # files into one E file         ]"\
  -cleanup   "0"  "1"        "[remove input files once done and _status=OK  ]"\
  -dplusp    "0"  "0"        "[move physics output into dynamics output     ]"\
  -prefix    ""	""         "[add prefix to name of all files              ]"\
  -yyoutgrid "U"  "U"        "[output grid for yin-yang: U, or GLB            ]"\
  -xferl     ""   ""         "[transfer model listings along with model output]"\
  -nthreads  "1"  "1"        "[Number of bemol process to run in parallel     ]"\
  -abortf    "Um_output" "Um_output" "[abort file                             ]"\
  -serdate   "0"  "1"        "[Keep date suffix for time series files         ]"\
  ++ $arguments`

printf "\n=====>  Um_output.sh starts: `date` ###########\n\n"
set -ex

abort_file=${TASK_WORK}/${abortf}_$$
touch ${abort_file} ; unset _ERR_from

DOMAIN=${domain}
mkdir ${TASK_WORK}/${DOMAIN}
cd ${TASK_WORK}/${DOMAIN}

rep_in=${TASK_INPUT}/${DOMAIN}
rep_ou=${TASK_WORK}/${DOMAIN}/files_2_xfer

. ${TASK_INPUT}/status_${DOMAIN##*_}

GEM_YINYANG=${GEM_YINYANG:-0}
if [ $GEM_YINYANG -gt 0 ] ; then
  if [ "$yyoutgrid" == "U"    ] ; then yy2e=1; fi
  if [ "$yyoutgrid" == "GLB"  ] ; then yy2e=2; fi
  yy2e=${yy2e:-1}
  if [ ${yy2e} -eq 2 ] ; then d2z=1 ; fi
fi

if [ -n "${xcasc}" ] ; then
  find -L ${rep_in} -type f -name "casc_*" 1> findcasc.lis 2> /dev/null
  cnt=`cat findcasc.lis | wc -l`
  if [ $cnt -gt 0 ] ; then
    mkdir -p ${xcasc} 
    printf "\n##### Linking casc files to ${xcasc} #####\n"
    for i in `cat findcasc.lis` ; do
      ln -sf $i ${xcasc}
    done
  fi
else
  REPCASC=casc_${RUNSTART:-0000000000}
fi

/bin/rm -rf $rep_ou prep ; mkdir -p $rep_ou prep ; cd prep

liste_of_prefix_to_treat="dm dp dh pm pp ph"

date
. r.call.dot ${TASK_BIN}/Um_output_prep.sh                  \
                     -input $rep_in -output $rep_ou          \
                     -assemble ${d2z:-0} -dplusp ${dplusp:-0}\
                     -liste ${liste_of_prefix_to_treat}      \
                     -listm ${xferl} -nthreads ${nthreads}   \
                     -repcasc ${REPCASC} -serdate ${serdate}
date

STATUS_prep=$_status
printf "STATUS_prep=$STATUS_prep\n"

if [ ${dplusp} -gt 0 ] ; then liste_of_prefix_to_treat="m p h" ; fi

if [ "$STATUS_prep" == "OK" ] ; then

  if [ ${GEM_YINYANG} -gt 0 ] ; then

    STATUS_yyoutg='ABORT'

    if [ ${yy2e} -eq 1 ] ; then
      STATUS_yyoutg='OK'
    else
      date
      . r.call.dot ${TASK_BIN}/Um_output_yyoutg.sh          \
                          -s ${rep_ou} -nthreads ${nthreads} \
                          -liste ${liste_of_prefix_to_treat}
      date
      STATUS_yyoutg=$_status	
    fi

    printf "STATUS_yyoutg=$STATUS_yyoutg\n"
    if [ "$STATUS_yyoutg" == "ABORT"  ] ; then
      _ERR_from='Um_output_yyoutg.sh'
    fi

  fi

else
  _ERR_from='Um_output_prep.sh'
fi

if [ -z "$_ERR_from" ] ; then
  /bin/rm -f ${abort_file}
fi

printf "\n=====>  Um_output.sh ends: `date` ###########\n\n"

