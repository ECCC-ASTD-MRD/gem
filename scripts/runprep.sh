#!/bin/bash
# Prepare input data for an interactive GEM integration and perform
# basic checks for configuration errors.
printf "\n=====>  runprep.sh starts: $(date) ###########\n"

# Store command line arguments
arguments=$*
printf "$0 ${arguments}\n\n"

# Process command line arguments
. r.entry.dot
eval `cclargs_lite -D " " $0 \
  -cfg           "0:0"          "0:0"       "[multi domains to run      ]"\
  -dircfg        "configurations/GEM_cfgs"  "configurations/GEM_cfgs"  "[location of config files  ]"\
  -npe           "1"            "1"         "[# of simultaneous threads ]"\
  -checkpart     "1x1"         "1x1"        "[MPI topology to check     ]"\
  -verbose       "0"            "1"         "[verbose mode              ]"\
  -_status       "ABORT"      "ABORT"       "[Return status             ]"\
  ++ $arguments`

export CMCCONST=${CMCCONST:-${ATM_MODEL_DFILES}/datafiles/constants}

DOMAIN_start=$(echo ${cfg} | cut -d":" -f1)
DOMAIN_end=$(  echo ${cfg} | cut -d":" -f2)

ici=${PWD}
cd PREP ; cd ${ici}
rm -fr PREP/*
work=${PWD}/PREP/work
mkdir -p ${work}

ptopo=$(echo ${checkpart} | tr "X" "x")
npex=$(echo ${ptopo} | cut -d "x" -f 1)
npey=$(echo ${ptopo} | cut -d "x" -f 2)

abort_prefix=err_
domain_number=$(printf "%04d" $DOMAIN_start)
while [ ${domain_number} -le ${DOMAIN_end} ] ; do
   dname=cfg_${domain_number}
   DEST=${PWD}/PREP/output/${dname}
   /bin/rm -rf ${DEST}
   mkdir -p ${DEST}
   . ${dircfg}/${dname}/configexp.cfg
   if [[ -n ${GEM_headscript_E} ]] ; then
      headscript_E=$(which ${GEM_headscript_E})
      if [[ ! -x ${headscript_E:-__No_Such_File_-} ]] ; then
         echo "ERROR: cannot find/execute headscript ${GEM_headscript_E}" 1>&2
         _status='ABORT'
         . r.return.dot
      fi
   else
      headscript_E=''
   fi
   GEM_check_settings=${GEM_check_settings:-0}
   . r.call.dot prep_domain.sh -anal ${GEM_anal} -input \'${GEM_inrep}\' \
      -o ${DEST} -work ${work} -headscript ${headscript_E} -bin \
      -check_namelist ${GEM_check_settings} \
      -npex ${npex} -npey ${npey} -cache $GEM_cache \
      -nmlfile ${dircfg}/${dname}/gem_settings.nml -nthreads $npe\
      -verbose ${verbose} -abort ${abort_prefix}${domain_number}
   if [[ -n ${GEM_inrep} ]] ; then
      ln -s ${GEM_inrep} ${DEST}/MODEL_inrep || true
   else
      mkdir -p ${DEST}/MODEL_inrep || true
   fi
   domain_number=$(printf "%04d" $(( domain_number+1 )))
done

if [[ $(find PREP/work -name "${abort_prefix}*" | wc -l) -gt 0 ]] ; then
   printf "One or more prep_domain.sh function calls aborted ...\n\n"
   _status='ABORT'
else
   _status='ED'
fi

. r.return.dot
