#!/bin/bash
printf "\n=====>  rungem.sh starts: $(date) ###########\n"

# Store command line arguments
arguments=$*
printf "$0 ${arguments}\n\n"

# Process command line arguments
. r.entry.dot
eval `cclargs_lite -D "" $0 \
  -npex      "1"     "1"     "[Block partitioning along-x     ]"\
  -npey      "1"     "1"     "[Block partitioning along-y     ]"\
  -nomp      "1"     "1"     "[Number of OMP threads          ]"\
  -mimd      ""      ""      "[mimd configuration file        ]"\
  -dom_start "1"     "1"     "[Starting domain number         ]"\
  -dom_end   "1"     "1"     "[Ending domain number           ]"\
  -inorder   "0"     "5"     "[Ordered listing                ]"\
  -barrier   "0"     "0"     "[DO NOT run binary              ]"\
  -debug     "0"     "gdb"   "[Debug option: gdb, ddt         ]"\
  -_status   "ABORT" "ABORT" "[return status                  ]"\
  -_endstep  ""      ""      "[return last time step performed]"\
  ++ $arguments`

_status="ABORT"

ndomains=$((dom_end - dom_start + 1))
npe_total=$(( npex * npey * ndomains ))
unset cfglist
idom=${dom_start}
while [ ${idom} -le ${dom_end} ] ; do
   cfglist="${cfglist} cfg_$(printf "%4.4d" ${idom})"
   idom=$((idom+1))
done

export DOMAINS_this_instance=${cfglist}

set -x
#export OMP_STACKSIZE=4G
export OMP_NUM_THREADS=$nomp
export MKL_NUM_THREADS=1
set ${SETMEX:-+x}

printf "\n Running `readlink ${TASK_BIN}/ATM_MOD.Abs` on $npe_total ($npex x $npey) PEs:\n"
printf " OMP_STACKSIZE=$OMP_STACKSIZE\n"
printf " OMP_NUM_THREADS=$OMP_NUM_THREADS\n\n"
printf " ##### UM_TIMING: Um_model.sh STARTING AT: `date`\n"

if [ $barrier -gt 0 ] ; then

  printf "GEM_Mtask 1\n"
  r.barrier
  printf "GEM_Mtask 2\n"
  r.barrier
  printf "\n =====> Um_model.sh CONTINUING after last r.barrier\n\n"

else

  unset INORDER
  if [ ${inorder} -gt 0 ] ; then INORDER="-inorder -tag"; fi
  if [ -n "${mimd}" ] ; then
     unset all_apps
     total_cpus=0
     while read line ; do
        app=$(echo $line | awk '{print $1}')
        pes=$(echo $line | awk '{print $2}')
        pes=${pes:-1}
        if [ -n "${app}" ] ; then
        if [ $(echo $app | cut -c1) != "#" ] ;then
        if [ "${app}" == "gemdm" -o "${app}" == "gem" ] ; then
           app=${TASK_BIN}/ATM_MOD.Abs
           pes=$((npex*npey))
        fi
        total_cpus=$((total_cpus+pes))
      #  all_apps="${all_apps}"" -np ${pes} ${app} :"
        all_apps="${all_apps}"" ${app} +${pes}"
        fi
        fi
     done < ${mimd}
     #all_apps=$(echo $all_apps | sed -E 's/(.*):/\1/')
     #CMD="mpirun $all_apps"
      CMD="${TASK_BIN}/r.mpirun -npex ${total_cpus} -npey $ndomains -pgm ${all_apps} $INORDER -minstdout ${inorder} -nocleanup"
  else
     CMD="${TASK_BIN}/r.mpirun -pgm ${TASK_BIN}/ATM_MOD.Abs -npex $((npex*npey)) -npey $ndomains $INORDER -minstdout ${inorder} -nocleanup"
  fi
  if [[ "x${debug}" != "x0" ]] ; then
     [[ "x${debug}" == "xgdb" || "x${debug}" == "x1"  ]] && export debug=gdb || true
     if [[ "x$(which ${debug} 2>/dev/null)" == "x" ]] ; then
        printf "ERROR: cannot find requested debug tool '${debug}'\n"
        if [[ "x${debug}" == "xddt" ]] ; then
           printf "    Maybe you forgot to load forge? Try:\n    . ssmuse-sh -x main/opt/forge/20.0.3\n"
        fi
        exit 1
     fi
     CMD="${CMD} -debug ${debug}" ;
  fi
  
  printf "\n EXECUTING: $CMD\n\n"
  $CMD

fi

set +x
printf " ##### UM_TIMING: Um_model.sh ENDING AT: `date`\n"

status_file=./status_MOD.dot
nb_abort=0
nb_restart=0
nb_end=0

printf " ##### UM_TIMING: POST Um_model.sh STARTING AT: `date`\n"
set -x
for i in ${cfglist} ; do
  fn=${TASK_OUTPUT}/${i}/${status_file}
  if [ -s ${fn} ] ; then
    . ${fn}
    printf "STATUS_FROM_DOMAIN: ${i} $_status\n"
    if [ "$_status" = "ABORT" ] ; then ((nb_abort=nb_abort+1))    ; fi
    if [ "$_status" = "RS"    ] ; then ((nb_restart=nb_restart+1)); fi
    if [ "$_status" = "ED"    ] ; then ((nb_end=nb_end+1))        ; fi
  else
    _status="ABORT"
    ((nb_abort=nb_abort+1))
  fi

# Deal with special files: time_series*.bin, zonaux_* and *.hpm*
# Files will be transfered from ${TASK_WORK}/$i to ${TASK_OUTPUT}/$i

  cd ${i}  
  /bin/rm -rf busper
  if [ "$_status" = "ED" ] ; then
    REP=${TASK_OUTPUT}/${i}/`cat ${TASK_OUTPUT}/${i}/output_ready_MASTER | grep "\^last" | cut -d " " -f3 | sed 's/\^last//g'`/endstep_misc_files
    mkdir -p ${REP}
    if [ -d YIN ] ; then
      mkdir ${REP}/YIN         ${REP}/YAN   2> /dev/null || true
      mv YIN/time_series*.bin* ${REP}/YIN   2> /dev/null || true
      mv YAN/time_series*.bin* ${REP}/YAN   2> /dev/null || true
      mv YIN/[0-9]*/*.hpm      ${REP}/YIN   2> /dev/null || true
      mv YAN/[0-9]*/*.hpm      ${REP}/YAN   2> /dev/null || true
    else
      mv time_series*.bin* ${REP} 2> /dev/null || true
      mv [0-9]*/*.hpm      ${REP} 2> /dev/null || true
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
  if [ $nb_restart -eq $ndomains ] ; then
    _status="RS"
  else
    if [ $nb_end -eq $ndomains ] ; then
      _status="ED"
    fi
  fi
fi
set +ex

# End of task
. r.return.dot
