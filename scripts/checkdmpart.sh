#!/bin/bash
#
printf "\n=====>  checkdmpart.sh starts: $(date) ###########\n"

# Store command line arguments
arguments=$*
printf "$0 ${arguments}\n\n"

# Process command line arguments
. r.entry.dot
eval `cclargs_lite -D " " $0 \
  -cfg         "cfg_0000"          "cfg_0000"         "[Analysis file or archive    ]"\
  -gemnml      "gem_settings.nml"  "gem_settings.nml" "[Model namelist settings file]"\
  -npex           "1"  "1"  "[# of Pes along-x        ]"\
  -npey           "1"  "1"  "[# of Pes along-y        ]"\
  -cache          ""   ""   "[GEM_cache               ]"\
  -verbose        "0"  "1"  "[verbose mode            ]"\
  -_status        ""   ""   "[output status           ]"\
  ++ $arguments`

set -ex

BIN=$(which checkdmpart_${BASE_ARCH}.Abs)
if [[ -z "${BIN}" ]] ; then
   BIN=$(which checkdmpart)
fi

ici=${PWD}
ROOT_WORK=${PWD}/checkdmpart$$
WORKDIR=${ROOT_WORK}/${cfg}
/bin/rm -rf ${ROOT_WORK}
mkdir -p ${WORKDIR}
cp ${gemnml}  ${WORKDIR}/model_settings.nml
domain=$(echo ${cfg##*_} | sed 's/^0*//')
if [[ -z "${domain}" ]] ; then domain=0 ;fi

cd ${ROOT_WORK}
cat > checkdm.nml <<EOF
&cdm_cfgs
cdm_npex = ${npex}
cdm_npey = ${npey}
cdm_eigen_S='${cache}'
/
EOF

GRDTYP=$(fetchnml.sh grd_typ_s grid ${WORKDIR}/model_settings.nml)
#GRDTYP=$(rpy.nml_get -u -f ${WORKDIR}/model_settings.nml grid/grd_typ_s 2> /dev/null)
if [[ "$GRDTYP" == "GY" ]] ; then 
   mkdir -p ${WORKDIR}/YIN/000-000 ${WORKDIR}/YAN/000-000
   ln -s ${ATM_MODEL_DFILES:-${AFSISIO:-/home/binops/afsi/sio}}/datafiles/constants/thermoconsts ${WORKDIR}/YIN/000-000/constantes
   ln -s ${ATM_MODEL_DFILES:-${AFSISIO:-/home/binops/afsi/sio}}/datafiles/constants/thermoconsts ${WORKDIR}/YAN/000-000/constantes
   cp checkdm.nml ${WORKDIR}/YIN/000-000
   mv checkdm.nml ${WORKDIR}/YAN/000-000
   ngrids=2
else
   mkdir -p ${WORKDIR}/000-000
   ln -s ${ATM_MODEL_DFILES:-${AFSISIO:-/home/binops/afsi/sio}}/datafiles/constants/thermoconsts ${WORKDIR}/000-000/constantes
   mv checkdm.nml ${WORKDIR}/000-000
   ngrids=1
fi

export DOMAIN_start=${domain}
export DOMAIN_end=${domain}
export DOMAIN_total=1
export TASK_INPUT=${PWD}
export TASK_WORK=${PWD}
export CCARD_OPT='ABORT'
export CCARD_ARGS="-dom_start ${domain} -dom_end ${domain} -npex 1 -npey 1 -ngrids ${ngrids} -input ${TASK_INPUT} -output ${TASK_OUTPUT}"
if [[ -n "${cache}" ]] ; then
   domain_number=$(printf "%04d" $domain)
   ln -s $cache ${TASK_INPUT}/cfg_${domain_number}/CACHE
fi

set +ex
printf "\n RUNNING ${BIN} \n"
lis=checkdmpartlis$$
echo checkdmpart_status='ABORT' > checkdmpart_status.dot
gem_mpirun.sh -pgm ${BIN} -npex ${ngrids} -inorder 1> $lis 2>&1

. checkdmpart_status.dot
grep topo_allowed checkdmpart_status.dot > $TMPDIR/listopoallowed$$
cnt=$(cat $TMPDIR/listopoallowed$$ | wc -l)

if [ "${checkdmpart_status}" != 'OK' -o ${verbose} -gt 0 ] ; then cat $lis ; fi

cd ${ici}
/bin/rm -rf ${ROOT_WORK} || true

if [ "${checkdmpart_status}" != 'OK' ] ; then
   bin=$(basename ${BIN})
   printf "\n  Error: Problem with ${bin}\n\n"
   _status="ABORT_${bin}"
else
     printf "\n  MPI topology allowed\n"
     cat $TMPDIR/listopoallowed$$
     if [ -n "${MAX_PES_IO}" ] ; then
        printf "\n  MAXIMUM number of I/O PES for this configuration is: $(echo ${MAX_PES_IO} | sed 's/^0*//')\n\n"
     fi
     _status='OK'
  if [ "${SOLVER}" != 'OK' ] ; then
   printf "\n  Error: VERTICAL LAYERING IS INCOMPATIBLE WITH THE TIMESTEP"
   printf "\n         THE SOLVER WILL NOT WORK\n\n"
   _status='ABORT_solver'
  fi
fi
/bin/rm -f $TMPDIR/listopoallowed$$ 
. r.return.dot

