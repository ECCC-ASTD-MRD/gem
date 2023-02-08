#!/bin/bash
printf "\n=====>  prep_domain.sh starts: $(date) ###########\n"

# Store command line arguments
arguments=$*
printf "$0 ${arguments}\n\n"

# Process command line arguments
eval `cclargs_lite -D " " $0 \
  -anal        ""          ""       "[Analysis file or archive     ]"\
  -input       ""          ""       "[Model input file path        ]"\
  -o           ""          ""       "[Domain output path           ]"\
  -work        "${TASK_WORK}" ""    "[Work directory path          ]"\
  -bin         "${TASK_BIN}"  ""    "[Binaries directory path      ]"\
  -headscript  ""          ""       "[Headscript to run            ]"\
  -check_namelist "0"      "0"      "[Check namelist validity      ]"\
  -nmlfile     ""          ""       "[Model namelist settings file ]"\
  -npex        "1"         "1"      "[# of processors along x      ]"\
  -npey        "1"         "1"      "[# of processors along y      ]"\
  -cache       ""          ""       "[GEM_cache                    ]"\
  -nthreads    "1"         "1"      "[# of simultaneous threads    ]"\
  -verbose     "0"         "1"      "[verbose mode                 ]"\
  -abort       ""          ""       "[Abort signal file prefix     ]"\
  ++ $arguments`

# Check for required inputs
set ${SETMEX:-+ex}

if [[ -z "${o}" ]] ; then
   echo "Error: output path (-o) must be defined for $0" >&2
   exit 1
fi

mydomain=$(basename ${o})

# Normalize path names
mkdir -p ${work} ${o}
for fpath in anal input o work headscript nmlfile ; do
   target=$(eval echo \$${fpath})
   if [[ -e ${target} ]] ; then 
      eval ${fpath}=$(readlink -e ${target})
   fi
done
# Prepare abort file in case of early return
if [[ -n "${abort}" ]] ; then 
   abort_file=${work}/${abort}-$$
   touch ${abort_file}
fi

if [ -n "${bin}" ] ; then
  bin=${bin}/
fi

set -ex
work=${work}/${mydomain}
mkdir -p ${work} ; cd ${work}

date

if [ -n "${nmlfile}" ] ; then
if [ -e "${nmlfile}" ] ; then
   # Verify namelist entries on request
   if [ ${check_namelist} -gt 0 ] ; then
      nml_to_check="convection_cfgs dyn_fisl dyn_kernel gem_cfgs grid hvdif init out physics_cfgs series step surface_cfgs vert_layers ensembles"
      #${bin}checknml --nml="${nml_to_check}" -r -- ${nmlfile}
      ${bin}checknml --nml="${nml_to_check}" -- ${nmlfile}
   fi
   if [ ${npex} -gt 1 -o ${npey} -gt 1 ] ; then
      . r.call.dot ${bin}checkdmpart.sh -gemnml ${nmlfile} -cfg ${mydomain} -cache "${cache}" -npex ${npex} -npey ${npey} -verbose $verbose
      if [ "${_status}" != 'OK' ] ; then
         echo $+status ; exit 1
      fi
   fi
fi
fi
date

# Prepare output and work space
target_dir=${o}
tmp_analysis_path=tmp_analysis ; mkdir -p ${tmp_analysis_path}
tmp_analysis_path=$(readlink -e ${tmp_analysis_path})
cd ${tmp_analysis_path}

# Extract cmc archive file input or link in surface analysis
local_anal_file=${tmp_analysis_path}/ANALYSIS
if [ -e "${anal}" ] ; then
   if [[ -L ${anal} ]] ; then
      analysis=$(readlink ${anal})
   else
      analysis=${anal}
   fi

   if [ -d ${analysis} ] ; then
     set -x
     cd ${analysis}
     for item in $(ls -1 .) ; do
        if [[ -f ${item} ]] ; then
           editfst -s ${item} -d ${local_anal_file} -i 0
        elif [[ -d ${item} ]] ; then
           ln -s $(readlink -e $item) ${o}/IAUREP
        fi
     done
     set +x
   else
      is_cmcarc=${work}/.is_cmcarc ; rm -f ${is_cmcarc}
      if r.filetype ${analysis} -t 36 37 ; then touch ${is_cmcarc} ; fi
      if [[ -e ${is_cmcarc} ]] ; then
        mkdir -p ${work}/unpack_archive
        cd ${work}/unpack_archive
        set -x
        cmcarc -x -f ${analysis}
        for item in $(ls -1 .) ; do
          if [[ -f ${item} ]] ; then
            editfst -s ${item} -d ${local_anal_file} -i 0
            /bin/rm -f ${item}
          elif [[ -d ${item} ]] ; then
            mv $item ${o}/IAUREP
          fi
        done
        set +x
      else
        set -x
        cp ${analysis} ${local_anal_file}
        set +x
      fi
   fi
   cd ${work}

   if r.filetype ${local_anal_file} -t 1 33 ; then

   # Run the user-defined headscript
   final_file=${o}/ANALYSIS
   if [[ -x ${headscript} ]] ; then
      local_touchup_file=${local_anal_file}_touchup
      ${headscript} ${local_anal_file} ${local_touchup_file}
      mv ${local_touchup_file} ${final_file}
   else 
      mv ${local_anal_file} ${final_file}
   fi

   else
      printf "\n   INCORRECT FORMAT FOR FILE ${analysis} - ABORT  \n\n"
      exit 1
   fi
else
      printf "\n   FILE ${anal} - NOT FOUND - ABORT  \n\n"
      exit 1
fi
date
set -ex
gtype=$(rpy.nml_get -u -f ${nmlfile} -- grid/grd_typ_s 2>/dev/null)
#gtype=$(fetchnml.sh grd_typ_s grid ${nmlfile})
if [ "${gtype}" == "'GU'" -o "${gtype}" == "'GY'" ] ; then
  input=''
fi

# Temporary check
#HEIGHT=$(fetchnml.sh dynamics_Kernel_S dyn_kernel ${nmlfile})
HEIGHT=$(rpy.nml_get -u -f ${nmlfile} -- dyn_kernel/dynamics_Kernel_S 2>/dev/null)
#HYDRO=$(fetchnml.sh dynamics_hydro_l dyn_kernel ${nmlfile})
HYDRO=$(rpy.nml_get  -u -f ${nmlfile} -- dyn_kernel/dynamics_hydro_l 2>/dev/null)
if [ -n "${HEIGHT}" -a -n "${HYDRO}" ] ; then
   HEIGHT=$(echo $HEIGHT | cut -d"_" -f3)
   HYDRO=$(echo $HYDRO | cut -c2)
   if [ "${HEIGHT}" == "H" -o "${HEIGHT}" == "h" ] ; then
      if [ "${HYDRO}" == "T" -o "${HYDRO}" == "t" ] ; then
         printf "\n ====> Dynamics_hydro_L= .TRUE. Not allowed in GEM-H - ABORT\n\n"
         exit 1
      fi
   fi
fi

# Preparation splitting of input files
splitdir=${o}/analysis
inrepdir=${o}/model_inrep
mkdir -p ${splitdir} ${inrepdir}
wild='@NIL@'
input=$(echo $input | sed "s/'//g")
if [ -n "${input}" ] ; then
  if [ ! -d ${input%%\**} ] ; then
    wild=${input##*/}
    input=$(dirname ${input%%\**})
  fi
fi

if [[ -e ${final_file} ]] ; then
   ${bin}GEM_trim_input.sh ${final_file} '@NIL@' ${splitdir} ${nthreads}
   dir=${splitdir}/$(ls -1 ${splitdir} | head -1)
   for i in $(ls -1 ${dir}/GEM_input_file*) ; do
      varname='TT'
      if [ $(r.fstliste -izfst $i -nomvar ${varname} | wc -l) -le 0 ] ; then
         varname='TT1'
      fi
      if [ $(r.fstliste -izfst $i -nomvar ${varname} | wc -l) -le 0 ] ; then
         varname='VT'
      fi
      valid=$(r.fstliste -izfst $i -nomvar ${varname} | head -1 | cut -d ":" -f 11)
      if [ -n "${valid}" ] ; then
         echo $(echo $valid | cut -c1-8).$(echo $valid | cut -c9-14) > ${splitdir}/analysis_validity_date
         break
      fi
   done
fi
date
if [[ -d ${input} ]] ; then
   ${bin}GEM_trim_input.sh ${input} ${wild} ${inrepdir} ${nthreads}
fi
date
/bin/rm -rf ${work}/tmp_analysis ${work}/unpack_archive

# Final cleanup for return
if [[ -n "${abort_file}" ]] ; then rm -f ${abort_file} ; fi

