#!/bin/bash

# Prepare a "domain" (model config instance) for execution
arguments=$*
echo $0 $arguments
date
eval `cclargs_lite0 -D " " $0 \
  -anal        ""          ""       "[Analysis file or archive     ]"\
  -input       ""          ""       "[Model input file path        ]"\
  -o           ""          ""       "[Domain output path           ]"\
  -work        "${TASK_WORK}" ""    "[Work directory path          ]"\
  -bin         "${TASK_BIN}"  ""    "[Binaries directory path      ]"\
  -headscript  ""          ""       "[Headscript to run            ]"\
  -nmlfile     ""          ""       "[Model namelist settings file ]"\
  -npex        "1"         "1"      "[# of processors along x      ]"\
  -npey        "1"         "1"      "[# of processors along y      ]"\
  -nthreads    "1"         "1"      "[# of simultaneous threads    ]"\
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
       echo "file for ${fpath} is ok:" ${target}
   elif [[ -d ${target} ]] ; then
       echo "directory for ${fpath} is ok:" ${target}
   else
       echo "target for ${fpath} is NOT FOUND:" ${target}
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

work=${work}/${mydomain}
mkdir -p ${work} ; cd ${work}

date

# Prepare output and work space
target_dir=${o}
tmp_analysis_path=$PWD/tmp_analysis ; mkdir -p ${tmp_analysis_path}
#echo "Whole path of tmp_analysis_path is",${tmp_analysis_path}
cd ${tmp_analysis_path}

# Extract cmc archive file input or link in surface analysis
local_anal_file=${tmp_analysis_path}/ANALYSIS
if [ -e "${anal}" ] ; then
      analysis=${anal}

   if [ -f ${analysis} ] ; then
        set -x
        touch ${local_anal_file}
        cp ${analysis} ${local_anal_file}
        set +x
   else
        echo ${analysis} " file not found"
        exit 1
   fi
   cd ${work}

   # Run the user-defined headscript
   final_file=${o}/ANALYSIS
   if [[ -x ${headscript} ]] ; then
      local_touchup_file=${local_anal_file}_touchup
      ${headscript} ${local_anal_file} ${local_touchup_file}
      mv ${local_touchup_file} ${final_file}
   else 
      mv ${local_anal_file} ${final_file}
   fi
fi
date
set -ex
gtype=$(${bin}Um_fetchnml2.sh grd_typ_s grid ${nmlfile})
if [ "${gtype}" == "'GU'" -o "${gtype}" == "'GY'" ] ; then
  input=''
fi

# Preparation splitting of input files
splitdir=${o}/analysis.d
inrepdir=${o}/model_inrep.d
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
      if [ $(r.fstliste.sh -izfst $i -nomvar ${varname} | wc -l) -le 0 ] ; then
         varname='TT1'
      fi
      if [ $(r.fstliste.sh -izfst $i -nomvar ${varname} | wc -l) -le 0 ] ; then
         varname='VT'
      fi
      valid=$(r.fstliste.sh -izfst $i -nomvar ${varname} | head -1 | cut -d ":" -f 11)
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

