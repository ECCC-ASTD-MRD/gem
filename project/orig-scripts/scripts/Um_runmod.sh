#!/bin/bash

arguments=$*
echo $0 $arguments

eval `cclargs_lite0 -D " " $0 \
  -cfg           "0:0"          "0:0"       "[configurations to consider                ]"\
  -dircfg        "GEM_cfgs"     "GEM_cfgs"  "[location of config files                  ]"\
  -theoc         "${theoc:-0}"  "1"         "[theoretical case flag                     ]"\
  -timing        "0"            "0"         "[report performance timers                 ]"\
  -ptopo         ""             ""          "[MPI&OMP PEs topology (npex x npey x nomp) ]"\
  -instances     ""             ""          "[(yinyang x ndomains)                      ]"\
  -bind          "0"            "1"         "[Processor binding logical controler       ]"\
  -debug         "0"            "1"         "[Debug session if debug=1                  ]"\
  -task_basedir  "RUNMOD"       "RUNMOD"    "[name of task dir                          ]"\
  -no_setup      "0"            "1"         "[do not run setup                          ]"\
  -_status       "ABORT"        "ABORT"     "[return status                             ]"\
  -_endstep      ""             ""          "[last time step performed                  ]"\
  -_npe          "1"            "1"         "[number of subdomains                      ]"\
  ++ $arguments`

printf "\n=====>  Um_runmod.sh starts: cfg_$cfg `date` ###########\n\n"

buf=`echo ${ptopo:-1x1x1} | sed 's/x/ /g'`
npex=`echo   $buf | awk '{print $1}'`
npey=`echo   $buf | awk '{print $2}'`
nomp=`echo   $buf | awk '{print $3}'`
npex=${npex:-1}
npey=${npey:-1}
nomp=${nomp:-1}
smtdyn=0
smtphy=0
if [ $bind -gt 0 ] ; then
  bind=true
else
  bind=false
fi
_npe=$((npex*npey))

ngrids=1
for i in ${TASK_INPUT}/cfg_* ; do
  GRDTYP=`Um_fetchnml2.sh Grd_typ_s grid ${i}/model_settings.nml`
  if [ "$GRDTYP" == "GY" ] ; then
    export GEM_YINYANG=YES
    ngrids=2
  fi
  break
done

for i in ${TASK_INPUT}/cfg_* ; do
  dname=`basename $i`
  mkdir -p ${TASK_OUTPUT}/$dname ${TASK_WORK}/$dname
  if [ "$GRDTYP" == "GY" ] ; then
    mkdir -p ${TASK_WORK}/$dname/YIN ${TASK_WORK}/$dname/YAN
  fi
  if [ -e ${TASK_INPUT}/${dname}/configexp.cfg ] ; then
    cp ${TASK_INPUT}/${dname}/configexp.cfg ${TASK_OUTPUT}/$dname
  fi
  if [ -e ${dircfg}/cfg_all/BATCH_config.cfg ] ; then
    cat ${dircfg}/cfg_all/BATCH_config.cfg >> ${TASK_OUTPUT}/$dname/configexp.cfg
  fi
  /bin/rm -f ${TASK_WORK}/$dname/theoc
  if [ ${theoc} -gt 0 ] ; then
    touch ${TASK_WORK}/$dname/theoc
  fi

  if [[ ! -e ${TASK_WORK}/$dname/model_settings.nml ]] ; then echo TRUE for model_settings.nml in ${TASK_WORK}/$dname ; fi
  # Do not overwrite an existing model_settings.nml since a parent script (i.e. Maestro task Runmod) may have modified it
  if [[ ! -e ${TASK_WORK}/$dname/model_settings.nml ]] ; then cp ${TASK_INPUT}/${dname}/model_settings.nml ${TASK_WORK}/$dname ; fi
  if [[ -e ${TASK_INPUT}/${dname}/output_settings ]] ; then
    cp ${TASK_INPUT}/${dname}/output_settings ${TASK_WORK}/$dname
  fi
  if [[ -e ${TASK_INPUT}/${dname}/coupleur_settings.nml ]] ; then 
    cp ${TASK_INPUT}/${dname}/coupleur_settings.nml ${TASK_WORK}/$dname
  fi


  chmod u+w ${TASK_WORK}/${dname}/model_settings.nml
  cat >> ${TASK_WORK}/$dname/model_settings.nml <<EOF

 &cpus
  Cpus_npex    = $npex    ,  Cpus_npey    = $npey
  Cpus_nthreads_dyn= $smtdyn,  Cpus_nthreads_phy= $smtphy
/

EOF

  if [ -s ${TASK_INPUT}/${dname}/BUSPER.tar ] ; then
    (mkdir -p ${TASK_WORK}/$dname/busper ; cd ${TASK_WORK}/$dname/busper ; tar xvf ${TASK_INPUT}/${dname}/BUSPER.tar)
  fi

done

set ${SETMEX:-+ex}
#typeset -Z4 domain_number

export DOMAIN_start=$(echo $cfg | cut -d : -f1)
export DOMAIN_end=$(  echo $cfg | cut -d : -f2)
export DOMAIN_total=$((DOMAIN_end - DOMAIN_start + 1))
DOMAIN_wide=$( echo $cfg | cut -d : -f3)
DOMAIN_wide=${DOMAIN_wide:-${DOMAIN_total}}
if [ $DOMAIN_wide -lt 1 ] ; then
  DOMAIN_wide=1
fi
export DOMAIN_wide=${DOMAIN_wide}

export BATCH_launch_execdir=${BATCH_launch_execdir:-${PWD}}
export GEMBNDL_VERSION=${ATM_MODEL_VERSION}

# Use performance timers on request
if [ ${timing} -gt 0 ] ; then export TMG_ON=YES; fi

cd $TASK_WORK

if [ -f ${TASK_INPUT}/cfg_0000/restart.tar ] ; then
  (cd cfg_0000 ; tar xvf ${TASK_INPUT}/cfg_0000/restart.tar)
  printf "\nRunning in FORCED RESTART mode\n\n"
fi

DOM=$DOMAIN_start
echo "THE DOMAIN_START IS",$DOM, "AND DOMAIN_END IS", $DOMAIN_end
while [ $DOM -le $DOMAIN_end ] ; do

  last_domain=$((DOM+DOMAIN_wide-1))
  if [ $last_domain -gt $DOMAIN_end ] ;then
    last_domain=$DOMAIN_end
  fi
  loop_cfg=${DOM}:${last_domain}
  ndomains=$((last_domain - DOM + 1))
  if [ -n "$instances" ] ; then 
    Um_check_instances.sh $instances $((ngrids*ndomains))
    if [ $? -ne 0 ] ; then 
       echo "Abort in Um_check_instances.sh"
      exit 1
    fi
  fi
  export GEM_NDOMAINS=$loop_cfg

  domain_number=$DOM
  file2watch=${TASK_BASEDIR}/output/cfg_${domain_number}/output_ready
  mkdir -p ${TASK_WORK}/post_process_output_cfg_${domain_number}

# Run main program wrapper
  set +ex
  
  printf "\n LAUNCHING Um_model.sh for domain: cfg_${domain_number} $(date)\n\n"
  ${TASK_BIN}/UM_MODEL -npex $((npex*ngrids)) -npey $npey -nomp $nomp \
                                    -ndom $ndomains -debug $debug 

  set ${SETMEX:-+ex}

  DOM=$((DOM+DOMAIN_wide))
done
set +ex

printf "\n DONE LAUNCHING all domains $(date)\n\n"

# Config file cleanup
/bin/rm -rf ${TASK_WORK}/busper

printf "\n=====>  Um_runmod.sh ends: `date` ###########\n\n"

# End of task

