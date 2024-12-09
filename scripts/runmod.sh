#!/bin/bash
printf "\n=====>  runmod.sh starts: $(date) ###########\n"

# Store command line arguments
arguments=$*
printf "$0 ${arguments}\n\n"

# Process command line arguments
. r.entry.dot
eval `cclargs_lite -D " " $0 \
   -cfg           "0:0"        "0:0"       "[Configurations number/range (START:END)]"\
   -dircfg        "configurations/GEM_cfgs"   "configurations/GEM_cfgs"  "[Location of config files]"\
   -barrier       "0"          "0"         "[DO NOT run binary]"\
   -timing        "0"          "0"         "[Report performance timers]"\
   -ptopo         "1x1x1"      "1x1x1"     "[MPI & OMP PEs topology (NPEXxNPEYx NOMP)]"\
   -smt           ""           ""          "[SMT controler (AIX) (smtdyn x smtphy)]"\
   -along_Y       "1"          "0"         "[Distribute PEs alog Y axis first]"\
   -nodespec      "NoNe"       "NoNe"      "[Node distribution specification]"\
   -inorder       "0"          "5"         "[Order listing]"\
   -debug         "0"          "gdb"       "[Debug session: gdb, ddt]"\
   -task_basedir  "RUNMOD"     "RUNMOD"    "[Task dir name]"\
   -no_setup      "0"          "1"         "[Do not run setup]"\
   -_status       "ABORT"      "ABORT"     "[Return status]"\
   -_endstep      ""           ""          "[Last time step performed]"\
   -_npe          "1"          "1"         "[Number of subdomains]"\
  ++ ${arguments}`

export CMCCONST=${CMCCONST:-${ATM_MODEL_DFILES}/datafiles/constants}
unset MIMD_cfg

restart=0

if [ ${no_setup} = 0 ] ; then
   set ${SETMEX:-+ex}
   if [ -d ${task_basedir}/work ] ; then
      restart=$(find -L ${task_basedir}/work -type f -name "gem_restart" | wc -l)
   fi
   if [ ${restart} -gt 0 ] ; then
      export TASK_BASEDIR=$(readlink -e ${task_basedir})
      export TASK_INPUT=${TASK_BASEDIR}/input
      export TASK_BIN=${TASK_BASEDIR}/bin
      export TASK_WORK=${TASK_BASEDIR}/work
      export TASK_OUTPUT=${TASK_BASEDIR}/output
      printf "Running in NORMAL RESTART mode\n"
   else
      datadir=${TMPDIR}/modeldata.$$
      /bin/rm -fr ${task_basedir:-yenapas}/* ${task_basedir:-yenapas}/.setup ${datadir}
      mkdir -p ${datadir}
      TASK_CFGFILE=${TMPDIR}/mod$$.cfg
      setmod.sh \
         -cfg $(echo ${cfg} | cut -d : -f 1):$(echo ${cfg} | cut -d : -f 2) \
         -dircfg ${dircfg} -tsk_cfgfile ${TASK_CFGFILE} -dirdata ${datadir}
      if [ -z "${TASK_SETUP}" ] ; then
         cat <<EOF
  WARNING in $(basename $0)
  Env variable TASK_SETUP not defined, using task_setup.dot
EOF
      fi
      export TASK_SETUP=${TASK_SETUP:-task_setup.dot}
      printf "\n##### EXECUTING TASK_SETUP ##### ${TASK_SETUP}\n"
      set -e
      . ${TASK_SETUP} -f $TASK_CFGFILE --base $(pwd)/${task_basedir} \
         --verbose --clean
      set ${SETMEX:-+ex}
      printf "##### EXECUTING TASK_SETUP DONE...#####\n"
      printf "\n##### RESULT OF TASK_SETUP #####\n"
      ls -l ${TASK_BIN} ${TASK_INPUT}/cfg_*
      /bin/rm -fr ${TASK_CFGFILE} ${datadir}
   fi
   if [ -s ${dircfg}/MIMD.cfg ] ; then
      export MIMD_cfg=$(true_path ${dircfg}/MIMD.cfg)
   fi
fi

set ${SETMEX:-+ex}

ptopo=$(echo ${ptopo} | tr "X" "x")
smt=$(  echo ${smt}   | tr "X" "x")
npex=$(echo ${ptopo} | cut -d "x" -f 1)
npey=$(echo ${ptopo} | cut -d "x" -f 2)
nomp=$(echo ${ptopo} | cut -d "x" -f 3)
npex=${npex:-1}
npey=${npey:-1}
nomp=${nomp:-1}
smtdyn=$(echo ${smt} | cut -d "x" -f 1)
smtphy=$(echo ${smt} | cut -d "x" -f 2)
smtdyn=${smtdyn:-0}
smtphy=${smtphy:-0}

_npe=$((npex*npey))

ngrids=1
for i in ${TASK_INPUT}/cfg_* ; do
   GRDTYP=$(rpy.nml_get -u -f ${i}/model_settings.nml -- grid/grd_typ_s 2>/dev/null)
   OPSCFG=$(rpy.nml_get -u -f ${i}/model_settings.nml -- ops_cfgs/Ops_configuration_S 2>/dev/null)
#  GRDTYP=$(fetchnml.sh grd_typ_s grid ${i}/model_settings.nml)
#  OPSCFG=$(fetchnml.sh Ops_configuration_S ops_cfgs ${i}/model_settings.nml)
   if [ -n "${GRDTYP}" ] ; then
      if [ "$GRDTYP" == "GY" ] ; then ngrids=2 ; fi
   else
      if [ -n "${OPSCFG}" ] ; then ngrids=${OPSCFG##*:} ; fi
   fi
   unset GRDTYP OPSCFG
   break
done

for i in ${TASK_INPUT}/cfg_* ; do
   dname=$(basename $i)
   mkdir -p ${TASK_OUTPUT}/${dname} ${TASK_WORK}/${dname}
   if [ -e ${TASK_INPUT}/${dname}/configexp.cfg ] ; then
      cp ${TASK_INPUT}/${dname}/configexp.cfg ${TASK_OUTPUT}/${dname}
   fi
   if [[ ! -e ${TASK_WORK}/${dname}/model_settings.nml ]] ; then
      cp ${TASK_INPUT}/${dname}/model_settings.nml ${TASK_WORK}/${dname}
   fi
   if [[ -e ${TASK_INPUT}/${dname}/output_settings ]] ; then
      cp ${TASK_INPUT}/${dname}/output_settings ${TASK_WORK}/${dname}
   fi
   if [[ -e ${TASK_INPUT}/${dname}/coupleur_settings.nml ]] ; then
      cp ${TASK_INPUT}/${dname}/coupleur_settings.nml ${TASK_WORK}/${dname}
   fi
   if [[ -e ${TASK_INPUT}/${dname}/SERVER_IO ]] ; then
      export MIMD_cfg=${TASK_INPUT}/${dname}/SERVER_IO
   fi
   if [ -s ${TASK_INPUT}/${dname}/BUSPER.tar ] ; then
      (mkdir -p ${TASK_WORK}/$dname/busper ; \
         cd ${TASK_WORK}/$dname/busper ; \
         tar xvf ${TASK_INPUT}/${dname}/BUSPER.tar)
   fi
done
if [ -n "${MIMD_cfg}" ] ; then
   printf "\n  Content of file ${MIMD_cfg}:\n"
   cat ${MIMD_cfg}
fi

export DOMAIN_start=$(echo ${cfg} | cut -d: -f1)
export DOMAIN_end=$(  echo ${cfg} | cut -d: -f2)
export DOMAIN_total=$((DOMAIN_end - DOMAIN_start + 1))
DOMAIN_wide=$(echo ${cfg} | cut -d : -f3)
DOMAIN_wide=${DOMAIN_wide:-${DOMAIN_total}}
if [ ${DOMAIN_wide} -lt 1 ] ; then
  DOMAIN_wide=1
fi
export DOMAIN_wide=${DOMAIN_wide}
alongYfirst=.false.
if [ $along_Y -gt 0  ] ; then alongYfirst=.true. ; fi
# Use performance timers on request
if [ ${timing} -gt 0 ] ; then export TMG_ON=YES      ; fi

if [[ "x${debug}" != "x0" ]] ; then
   export RPN_COMM_DIAG=2
   [[ "x${debug}" == "xgdb" || "x${debug}" == "x1"  ]] && export debug=gdb || true
   if [[ "x$(which ${debug} 2>/dev/null)" == "x" ]] ; then
      printf "ERROR: cannot find requested debug tool '${debug}'\n"
      if [[ "x${debug}" == "xddt" ]] ; then
         printf "    Maybe you forgot to load forge? Try:\n    . ssmuse-sh -x main/opt/forge/20.0.3\n"
      fi
      exit 1
   fi
fi

cd $TASK_WORK

if [ -f ${TASK_INPUT}/cfg_0000/restart.tar ] ; then
  (cd cfg_0000 ; tar xvf ${TASK_INPUT}/cfg_0000/restart.tar)
  printf "\nRunning in FORCED RESTART mode\n\n"
fi

export CCARD_OPT='ABORT'
nb_abort=0
DOM=${DOMAIN_start}
while [ ${DOM} -le ${DOMAIN_end} ] ; do

   last_domain=$((DOM+DOMAIN_wide-1))
   if [ ${last_domain} -gt ${DOMAIN_end} ] ;then
      last_domain=${DOMAIN_end}
   fi
   export CCARD_ARGS="-dom_start ${DOM} -dom_end ${last_domain} -dom_last ${DOMAIN_end} -npex ${npex} -npey ${npey} -ngrids ${ngrids} -smtdyn $smtdyn -smtphy $smtphy -along_Y ${alongYfirst} -input ${TASK_INPUT} -output ${TASK_OUTPUT}"
   export GEM_NDOMAINS=${DOM}:${last_domain} # for launch_sortie.sh

   domain_number=$(printf "%04d" $DOM)
   file2watch=${TASK_BASEDIR}/output/cfg_${domain_number}/output_ready
   MONLIS=gem_monitor_output_cfg_${domain_number}.lis
   mkdir -p ${TASK_WORK}/post_process_output_cfg_${domain_number}
   gem_monitor_output ${file2watch} ${TASK_BIN}/launch_sortie.sh 1>> ${MONLIS} 2>&1

   # Run main program wrapper 'rungem.sh'
   set +ex
   printf "\n LAUNCHING rungem.sh for domain: cfg_${domain_number} $(date)\n\n"
   . r.call.dot ${TASK_BIN}/rungem.sh \
      -npex $((npex*ngrids)) -npey $npey -nomp $nomp \
      -mimd ${MIMD_cfg} \
      -dom_start ${DOM} -dom_end ${last_domain} -debug $debug \
      -barrier ${barrier} -inorder ${inorder}

   cnt=$(ls -1 ${TASK_WORK}/post_process_output_*/output_ready_work*.active 2> /dev/null | wc -l)
   while [ ${cnt} -gt 0 ] ; do
      printf "Waiting for: $cnt launch_sortie.sh process to end\n"
      sleep 10
      cnt=$(ls -1 ${TASK_WORK}/post_process_output_*/output_ready_work*.active 2> /dev/null| wc -l)
   done

   set ${SETMEX:-+ex}
   if [ "$_status" = "ED" -o "$_status" = "RS" ] ; then
      if [ "$_status" = "ED" ] ; then echo ${Runmod} > ${TASK_OUTPUT}/last_npass ; fi
      ${TASK_BIN}/launch_sortie.sh ${file2watch}_MASTER 1>> ${MONLIS} 2>&1
   fi
   if [ "$_status" = "ABORT" ] ; then ((nb_abort=nb_abort+1)) ; fi

   DOM=$((DOM+DOMAIN_wide))
done

if [ $nb_abort -gt 0 ] ; then 
  _status="ABORT"
fi
set +ex

printf "\n DONE LAUNCHING all domains $(date)\n\n"

# Config file cleanup
/bin/rm -rf ${TASK_WORK}/busper

printf "\n=====>  rungem.sh ends: $(date) ###########\n\n"

# End of task
. r.return.dot

