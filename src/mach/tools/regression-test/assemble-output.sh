#!/bin/bash

###
#
#  Script for assembling output of GEM-MACH run and submitting the 
#      script for validation of the regression test results.  
#      Generated combined output files can be used as the new control.
#
#  USAGE: 
#      ord_soumet ${TASK_BIN}/assemble-output.sh -args "${TASK_BASEDIR}" -mach ${Test_JobsMach} -cpus ${AssmblJobNcpu} -cm ${AssmblJobMemory} -t ${AssmblJobTime} -queue ${Test_JobsQueue} -jn ${AssmblJobName} -listing ${BASE_LIST}
#      or
#      assemble-output.sh ${TASK_BASEDIR} 2>&1 > ${TASK_LIST}/gm-test-${TRUE_HOST}-assemble-listings.txt
#
#      Note that the variables provided under quotes after "-args"
#      will be read by this script to fill in the values of apropriate
#      variables.
#
#  The script relies on:
#      1. Existence of directory structure needed for Runmod task with
#         root directory being placed at ${TASK_BASEDIR}
#      2. Existence of output files from the successful GEM-MACH run 
#         saved in ${TASK_OUTPUT} 
#      3. Existence of configuration file named 
#         ${TASK_BASEDIR}/assemble_env_arg
#  If any of these are missing, or not up to date, the script will fail.
#  The error messages are provided throughout the script to indicate any
#  obvious issues and failure of the script to behave as expected.
#
#  Script:
#     - Following the succesfful model run, combines the output in the 
#       fst files named as YYYYMMDDHH_hhh, which can be used as new 
#       control output.
#     - Based on value of the ${cntrl_fl_opt} variable submits the script 
#       for validation of the model results against the existing control 
#       output.
#
#  Author: Balbir Pabla
#  Date:   August 2025
#
#         Jack Chen and Verica Savic-Jovcic contributed to development 
#         of the scripts from which this one was created.
#
###

scriptstartdate=$(date '+%C%y%m%d%H%M%S')

# Read in arguments
export TASK_BASEDIR=$1

# Set up environment
source ${TASK_BASEDIR}/assemble_env_arg
cd ${TASK_OUTPUT}

# Update log
echo -e "\n== Starting assemble-output script: $scriptstartdate == \n" | tee -a ${gmtestinfo}

# Create a directory where to save the assembled output
[[ ! -d ${TASK_OUTPUT}/model ]] && mkdir -p ${TASK_OUTPUT}/model

# Specify forecast hours for which the output is assembled
if [ "${cmpl_opt}" == "dbg" ] ; then 
 fhrs=(-003 000 009)
else
 fhrs=(-003 000 009 012 021 024)
fi

# Loop through the forecast hours for which the output is assembled by combining physics and dynamics output files
for fhr in ${fhrs[*]} ; do
 run_step=$(basename $(dirname $(dirname ${TASK_OUTPUT}/cfg_0000/laststep_0000000???/000-000/pm??????????-000-000_${fhr})))
 rundate=$(basename ${TASK_OUTPUT}/cfg_0000/${run_step}/000-000/pm??????????-000-000_${fhr} | cut -c3-12)
 dmoutfile=${TASK_OUTPUT}/cfg_0000/${run_step}/000-000/dm${rundate}-000-000_${fhr}
 pmoutfile=${TASK_OUTPUT}/cfg_0000/${run_step}/000-000/pm${rundate}-000-000_${fhr}
 outfile=${TASK_OUTPUT}/model/${rundate}_${fhr}
 if [ ! -f ${dmoutfile} ] ; then
  echo -e "\n\n ERROR: dynamics output, ${dmoutfile}, unavailable. \n\n" | tee -a ${gmtestinfo}
  exit 1
 else
  echo -e "* dynamics output:\n ${dmoutfile} \n" | tee -a ${gmtestinfo}
 fi
 if [ ! -f ${pmoutfile} ] ; then
  echo -e "\n\n ERROR: physics output, ${pmoutfile}, unavailable. \n\n" | tee -a ${gmtestinfo}
  exit 1
 else
  echo -e "* physics output:\n ${pmoutfile} \n" | tee -a ${gmtestinfo}
 fi
 # Combine physics (pm) and dynamics (dm) outputs on hybrid levels to one output file
 for f in ${TASK_OUTPUT}/cfg_0000/${run_step}/00*/?m${rundate}-*_${fhr} ; do
  ${TASK_BIN}/editfst -s $f -d ${outfile} -i -m errors <<EOD
 exclure(-1,['>>','^^','!!','>^'],-1)
EOD
 done
 # Add coordinate information to the assembled output file
 ${TASK_BIN}/editfst -s ${dmoutfile} -d ${outfile} -i -m errors <<EOD
 desire(-1,['>>','^^','!!','>^'],-1)
EOD
 ${TASK_BIN}/editfst -s ${pmoutfile} -d ${outfile} -i -m errors <<EOD
 desire(-1,['>>','^^','!!','>^'],-1)
EOD
 echo -e "*** combined output: ${outfile} \n"  | tee -a ${gmtestinfo}
done

# Submit job to verify results against control files
if [[ "${sequence_opt}" == "true" ]] ; then
 if [[ "${cntrl_fl_opt}" != "new" ]] ; then
  echo -e "\n Validate results against a control run \n" | tee -a ${gmtestinfo}
  ord_soumet ${TASK_BIN}/validate-results.sh -args "${TASK_BASEDIR}" \
             -mach ${Test_JobsMach} -cpus ${PostProcJobProcTopo} -cm ${VldtJobMemory} \
             -t ${PostProcJobTime} -queue ${Test_JobsQueue} -jn ${PostProcJobName} \
             -listing ${TASK_LIST}
 else
  # Create expected control-directory structure
  ln -s ${TASK_BASEDIR}/input ${TASK_BASEDIR}/gm-input
  ln -s ${TASK_BASEDIR}/output ${TASK_BASEDIR}/gm-output
  echo -e "\n *** Location of new control output: $(dirname ${outfile}) *** \n " | tee -a ${gmtestinfo}
  echo -e "\n The following path can be used for 'DIR' in '-l DIR' option of regression test:\n${TASK_BASEDIR} \n " | tee -a ${gmtestinfo}
 fi
 # Tell the world how long it took to run GEM-MACH
 echo -e "\n== It took $(r.date -n -MM -L $(date '+%C%y%m%d%H%M%S') ${scriptstartdate}) seconds to run the script that assembles the output == \n" | tee -a ${gmtestinfo}
fi

