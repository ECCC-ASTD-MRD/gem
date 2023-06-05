#/bin/bash

set -e
echo "---------------> ${0} BEGIN <---------------"
###
# load RPN utils environment (r.entry.dot)
###
. r.load.dot /fs/ssm/eccc/mrd/rpn/utils/20220509
###
# read in arguments:
###
arguments=$*
echo ${0} ${arguments}
. r.entry.dot
eval `cclargs_lite -D " " ${0}                                                                     \
  -s          "$PWD"    "nil"    "[Directory FROM where GEMMACHv2/CHEM code shall be imported ]"   \
  -v          "0"       "1"      "[verbose switch; triggered if -v is typed in ]"                  \
  ++ ${arguments}`
###
# check whether the provided CHEM directory exists
###
if [ ! -d ${s} ]; then
   echo "${0} ERROR --> Import Directory ${s} DOES NOT EXIST : ABORT"
   exit 1
else
   sc=`true_path ${s}`
   echo "True path for CHEM code is ${sc}"
   cd ${sc}
fi
###
# prepare working directories for updating headers modules
###
mkdir -p imported_headers regenerated_headers
cp ${sc}/*header*.F90 imported_headers
if [ $? -eq 0 ]; then
   if [ $(ls -1 imported_headers/*header*.F90 | wc -l) -gt 0 ]
   then
      head_list=$(ls -1 imported_headers/*header*.F90)
   else
      echo "Error: no header files found in  ${s} "
      echo "ABORT ${0}"
      exit 1
   fi
else
   echo "Error: no header files found in  ${s} "
   echo "ABORT ${0}"
   exit 1
fi
###
# loop through header files to update them
#  in each file, loop through the existing subroutines and functions and update them
###
for hf in ${head_list}; do
   if [ ${v} -gt 0 ] ; then echo "Processing ${hf}" ; fi
   sublist=$(grep "end subroutine" ${hf} | sed 's/^\s*end\s*subroutine\s*//g')
   funlist=$(grep "end function"   ${hf} | sed 's/^\s*end\s*function\s*//g')
   if [ ${v} -gt 0 ] ; then echo ${sublist} ; fi
   if [ ${v} -gt 0 ] ; then echo ${funlist} ; fi
   dest_hf=$(basename ${hf})
   dest_hf_nosuffix=$(basename ${hf} .F90)
   sed -n '/\!begin trap head/,/\!end trap head/p' ${hf}               >  regenerated_headers/${dest_hf}
   echo -en "\n"                                                       >> regenerated_headers/${dest_hf}
   for sub in ${sublist}; do
      if [[ ${sub} == @("chm_exe2"|"chm_load_emissions2") ]] ; then
         subfn=$(basename ${sub} 2)
      elif [[ ${sub} == "mach_kpp_integrator" ]] ; then
         subfn=mach_kpp_Integrator
      else
         subfn=${sub}
      fi
      sed -n '/\!\!if_on/{:a;n;/\!\!if_off/b;p;ba}' ${sc}/${subfn}.F90 >> regenerated_headers/${dest_hf}
      echo "end subroutine ${sub}"                                     >> regenerated_headers/${dest_hf}
      echo -en "\n"                                                    >> regenerated_headers/${dest_hf}
   done
   for fun in ${funlist}; do
      sed -n '/\!\!if_on/{:a;n;/\!\!if_off/b;p;ba}' ${sc}/${fun}.F90 >> regenerated_headers/${dest_hf}
      echo "end function ${fun}"                                       >> regenerated_headers/${dest_hf}
      echo -en "\n"                                                    >> regenerated_headers/${dest_hf}
   done
   echo "end interface"                                                >> regenerated_headers/${dest_hf}
   echo "end module ${dest_hf_nosuffix}"                               >> regenerated_headers/${dest_hf}
done

echo "---------------> ${0} NORMAL END <---------------"
exit 0
