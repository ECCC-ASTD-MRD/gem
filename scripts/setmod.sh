#!/bin/bash

arguments=$*
printf "\n$0 ${arguments}\n\n"

eval `cclargs_lite -D " " $0 \
  -cfg          "0:0"        "0:0"       "[configurations to consider]"\
  -dircfg       "GEM_cfgs"   "GEM_cfgs"  "[location of config files  ]"\
  -dirdata      ""           ""          "[directory to assemble inputs]"\
  -tsk_cfgfile   ""          ""          "[task_setup config file    ]"\
  ++ $arguments`
DOMAIN_start=`echo $cfg | sed 's/:/ /' | awk '{print $1}'`
DOMAIN_end=`  echo $cfg | sed 's/:/ /' | awk '{print $2}'`

TASK_CFGFILE=$tsk_cfgfile

############ FUNCTION tskcfg ###############
tskcfg ()
{
arguments=$*
eval `cclargs_lite -D " " $0 \
  -gem_cfgfile   ""  ""    "[gem config file       ]"\
  -tsk_cfgfile   ""  ""    "[task_setup config file]"\
  -nml           ""  ""    "[gem namelist file     ]"\
  -phytbl        ""  ""    "[physics input table file]"\
  -out_cfg       ""  ""    "[gem outcfg.out file   ]"\
  -rep2cfg       ""  ""    "[configuration to set  ]"\
  -data_cfg      ""  ""    "[input data for domain ]"\
  ++ $arguments`
set -x
default_DIR_ATM=${ATM_MODEL_DFILES}
default_PATH_ATMDATA=${ATM_MODEL_DFILES}/datafiles/modeldata
default_DIR_SIO=${CMCCONST}

geophy=${default_DIR_ATM}/bcmk/geophy
climato=${default_DIR_ATM}/bcmk/climato

PREP_dir=${PWD}/PREP/output/${rep2cfg}
IRTAB=${default_DIR_SIO}/irtab5_std
OZONE=${default_DIR_SIO}/ozoclim_phy45
CONST=${default_DIR_SIO}/thermoconsts

outcfg=${out_cfg:-'<no value>'}
phytbl=${phytbl:-'<no value>'}

CACHEDIR='<no value>'
iaurep='<no value>'
restart='<no value>'
busper='<no value>'
analysis='<no value>'
mod_analysis='<no value>'
mod_inrep='<no value>'
if [ -e ${PREP_dir}/ANALYSIS ] ; then
  analysis=${PREP_dir}/ANALYSIS
fi
if [ -e ${PREP_dir}/analysis ] ; then
  mod_analysis=${PREP_dir}/analysis
fi
if [ -e ${PREP_dir}/model_inrep ] ; then
  mod_inrep=${PREP_dir}/model_inrep
fi
if [ -e ${PREP_dir}/IAUREP ] ; then
  iaurep=${PREP_dir}/IAUREP
fi

ATMMOD=$(which maingemdm)
BINMOD=$(dirname ${ATMMOD})

if [ -n "$gem_cfgfile" ] ; then
  . $gem_cfgfile
  config=$gem_cfgfile
  if [ -n "${GEM_model_input}" ] ; then
    input=${GEM_model_input}"/*"
  fi
  geophy=${GEM_geophy:-$geophy}
  climato=${GEM_climato:-$climato}
  iaurep=${GEM_iaurep:-$iaurep}
  restart=${GEM_restart_tarbal:-$restart}
  busper=${GEM_busper:-$busper}
  CONST=${GEM_const:-${CONST}}
  IRTAB=${GEM_radtab:-${IRTAB}}
  OZONE=${GEM_ozone:-$OZONE}
  if [ ! -e ${OZONE} ] ; then
    OZONE=${default_DIR_SIO}/${GEM_ozone}
  fi
  PHYTB=${GEM_phy_intable:-${phytbl}}
  CACHEDIR=${GEM_cache:-$CACHEDIR}
  BINMOD=${GEM_ovbin:-${BINMOD}}
  ATMMOD=${BINMOD}/maingemdm
else
  config='<no value>'
fi

cat >> $tsk_cfgfile <<EOF
# ${rep2cfg}/configexp.cfg $config
# ${rep2cfg}/model_settings.nml ${nml}
# ${rep2cfg}/output_settings $outcfg
# ${rep2cfg}/ANALYSIS $analysis
# ${rep2cfg}/GEOPHY $geophy
# ${rep2cfg}/CLIMATO $climato
# ${rep2cfg}/MODEL_ANALYSIS $mod_analysis
# ${rep2cfg}/MODEL_INREP $mod_inrep
# ${rep2cfg}/MODEL_INPUT ${input} ${default_PATH_ATMDATA}/*
# $rep2cfg/RESTART.tar $restart
# $rep2cfg/BUSPER.tar $busper
# $rep2cfg/CACHE $CACHEDIR
# ${rep2cfg}/IAUREP ${iaurep}
# ${rep2cfg}/constantes $CONST
# ${rep2cfg}/ozone_clim.fst $OZONE
# ${rep2cfg}/rad_table.fst $IRTAB
# ${rep2cfg}/physics_input_table $PHYTB
EOF
}
############ END FUNCTION tskcfg ###############

#=====> Main scripts

cat > $TASK_CFGFILE <<EOF
#############################################
# <input>
EOF

domain_number=$(printf "%04d" $DOMAIN_start)
while [ $domain_number -le $DOMAIN_end ] ; do
  unset CFGFILE OUTCFG CPLNML PHYTABLE
  dname=cfg_$domain_number
  if [ -e ${dircfg}/${dname}/configexp.cfg ] ; then
    CFGFILE=${dircfg}/${dname}/configexp.cfg
  else
    if [ -e ${PWD}/configexp.cfg ] ; then
      CFGFILE=${PWD}/configexp.cfg
    fi
  fi
  if [ -e ${dircfg}/${dname}/outcfg.out ] ; then
    OUTCFG=${dircfg}/${dname}/outcfg.out
  fi
  if [ -e ${dircfg}/${dname}/physics_input_table ] ; then
    PHYTABLE=${dircfg}/${dname}/physics_input_table
  else
#   default_phytab=${gem_DIR:+${gem_DIR}/src/rpnphy/include}
    default_phytab=${GEM_STORAGE_DIR:+${gem_DIR}/src/rpnphy/include}
    default_phytab=${default_phytab:-${gem_DIR}/share/rpnphy}
    PHYTABLE=${default_phytab}/physics_input_table
  fi
  tskcfg -rep2cfg ${dname} -gem_cfgfile $CFGFILE  \
         -nml ${dircfg}/${dname}/gem_settings.nml \
         -phytbl ${PHYTABLE} -out_cfg $OUTCFG     \
         -tsk_cfgfile $TASK_CFGFILE -data_cfg ${dirdata}/${dname}
  domain_number=$(printf "%04d" $(( domain_number+1 )))
done
cat $TASK_CFGFILE

cat >> $TASK_CFGFILE <<EOF
# </input>
# <executables>
# ATM_MOD.Abs     ${ATMMOD}
# rungem.sh       rungem.sh 
# r.mpirun        gem_mpirun.sh
# launch_sortie.sh Um_process_output.sh
# </executables>
# <output>
${XCHGDIR_cfg}
# </output>
#############################################
EOF

printf "\n### Content of config file to TASK_SETUP ####\n\n"
cat $TASK_CFGFILE 2>/dev/null

