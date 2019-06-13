#!/bin/bash

arguments=$*

eval `cclargs_lite0 $0 \
  -npex     "1"            "1"         "[Domain partitioning along-x]" \
  -npey     "1"            "1"         "[Domain partitioning along-y]" \
  -nomp     "1"            "1"         "[Number of OMP threads      ]" \
  -ndom     "1"            "1"         "[Number of domains          ]" \
  -debug    "0"            "1"         "[Debug session if debug=1   ]"\
  -_status  "ABORT"        "ABORT"     "[return status              ]" \
  -_endstep ""             ""          "[last time step performed   ]" \
  ++ $arguments`


set ${SETMEX:-+ex}

export OMP_NUM_THREADS=$nomp

npe_total=$(( npex * npey * ndom ))

printf "\n Running  ${TASK_BIN}/ATM_MOD on $npe_total ($npex x $npey) PEs:\n"
printf " OMP_NUM_THREADS=$OMP_NUM_THREADS\n\n"
if [ $nomp -gt 1 ] ; then 
export OMP_STACKSIZE=4G
printf " OMP_STACKSIZE=$OMP_STACKSIZE\n\n"
fi

printf " ##### UM_TIMING: Um_model.sh STARTING AT: `date`\n"

#for CRAY
#CMD="aprun -j 1 -n $((npex*npey)) -d $nomp ${TASK_BIN}/ATM_MOD"
#for Linux
CMD="mpirun -np $((npex*npey)) ${TASK_BIN}/ATM_MOD"

printf "MPIRUN CMD is `echo $CMD` \n"
if [[ x$debug != x0 ]] ; then
    CMD="${CMD} -gdb"
fi
$CMD

set +ex
printf " ##### UM_TIMING: Um_model.sh ENDING AT: `date`\n"

set ${SETMEX:-+ex}

status_file=./status_MOD.dot
nb_abort=0
nb_restart=0
nb_end=0

printf " ##### UM_TIMING: POST Um_model.sh STARTING AT: `date`\n"
for i in cfg_* ; do
  fn=${TASK_OUTPUT}/${i}/${status_file}
  if [ -s ${fn} ] ; then
    . ${fn}
    printf "STATUS_FROM_DOMAIN: ${i} $_status\n"
    if [ "$_status" = "ABORT" ] ; then ((nb_abort=nb_abort+1))    ; fi
    if [ "$_status" = "RS"    ] ; then ((nb_restart=nb_restart+1)); fi
    if [ "$_status" = "ED"    ] ; then ((nb_end=nb_end+1))        ; fi
  fi

# Deal with special files: time_series.bin, zonaux_* and *.hpm*
# Files will be transfered from ${TASK_WORK}/$i to ${TASK_OUTPUT}/$i

  cd ${i}  
  /bin/rm -rf busper
  if [ "$_status" = "ED" ] ; then
    REP=${TASK_OUTPUT}/${i}/`cat ${TASK_OUTPUT}/${i}/output_ready_MASTER | grep "\^last" | cut -d " " -f3 | sed 's/\^last//g'`/endstep_misc_files
    mkdir -p ${REP}
    if [ -d YIN ] ; then
      mkdir ${REP}/YIN       ${REP}/YAN   2> /dev/null || true
      mv YIN/time_series.bin ${REP}/YIN   2> /dev/null || true
      mv YAN/time_series.bin ${REP}/YAN   2> /dev/null || true
      mv YIN/[0-9]*/*.hpm    ${REP}/YIN   2> /dev/null || true
      mv YAN/[0-9]*/*.hpm    ${REP}/YAN   2> /dev/null || true
    else
      mv time_series.bin ${REP} 2> /dev/null || true
      mv [0-9]*/*.hpm    ${REP} 2> /dev/null || true
    fi
    liste_busper=`find ./ -type f -name "BUSPER4spinphy*"`
    fn=`echo $liste_busper| awk '{print $1}'`
    if [ -n "${fn}" ] ; then
      mkdir -p ${REP}
      fn=`basename $fn`
      tar cvf ${REP}/${fn}.tar $liste_busper
      /bin/rm -f $liste_busper
    fi
  fi
  cd ../
  
done
printf " ##### UM_TIMING: POST Um_model.sh ENDING AT: `date`\n"

if [ $nb_abort -gt 0 ] ; then 
  _status="ABORT"
else
  if [ $nb_restart -eq $ndom ] ; then
    _status="RS"
  else
    if [ $nb_end -eq $ndom ] ; then
      _status="ED"
    fi
  fi
fi
set +ex

# End of task
