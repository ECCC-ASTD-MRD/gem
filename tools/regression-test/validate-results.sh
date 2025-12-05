#!/bin/bash

###
#
#  Script for validating GEM-MACH regression test results
#
#  USAGE: 
#      ord_soumet ${TASK_BIN}/validate-results.sh -args "${TASK_BASEDIR}" -mach ${Test_JobsMach} -cpus ${VldtJobNcpu} -cm ${Test_JobsMemory} -t ${VldtJobTime} -queue ${Test_JobsQueue} -jn ${VldtJobName} -listing ${TASK_LIST}
#      or
#      validate-results.sh ${TASK_BASEDIR} 2>&1 > ${TASK_LIST}/gm-test-${TRUE_HOST}-validation-listings.txt
#
#      Note that the variables provided under quotes after "-args"
#      will be read by this script to fill in the values of apropriate
#      variables.
#
#      If the test is successful, the script will print:
#      ***
#      New binary reproduces reference output.
#      ***
#      on the screen or in the listings file.
#
#  The script relies on:
#      1. Existence of directory structure needed for Runmod task with 
#         root directory being placed at ${TASK_BASEDIR}
#      2. Existence of configuration file named 
#         ${TASK_BASEDIR}/validate_env_arg
#      3. Existence of reference output file named 
#         ${control_dir}/gm-output/model/${rundate}_${fhr}
#  If any of these are missing, or not up to date, the script will fail.
#  Error messages are provided throughout the script to indicate any
#  obvious issues and failure of the script to behave as expected.
#
#  Script:
#      - Compares results of the run with the previously saved results 
#        of the control run performed with the same input and a binary
#        built with the last released GEM-MACH code.
#
# Author: Verica Savic-Jovcic
# Date:   March 2022
#
#         Jack Chen, Balbir Pabla and Verica contributed to updates afterwards.
#
###

scriptstartdate=$(date '+%C%y%m%d%H%M%S')

# Read in arguments 
TASK_BASEDIR=$1

# Set up environment
source ${TASK_BASEDIR}/validate_env_arg
cd ${TASK_OUTPUT}

# Update log
echo -e "\n== Starting validate-results script: $scriptstartdate == \n" | tee -a ${gmtestinfo}

# Prepare necessary tools
cp -a $(which fstcomp) ${TASK_BIN}/.

# Define forecast hour for which the output is validated
fhr=024
[[ "${cmpl_opt}" == "dbg" ]] && fhr=009

# Define location of control output
ctrl_output=${control_dir}/gm-output/model
[[ "${cmpl_opt}" == "dbg" ]] && ctrl_output=${control_dir}/gm-output_${cmpl_opt}/model

# Define forecast date
rundate=$(basename ${ctrl_output}/??????????_${fhr} |cut -c1-10)

# Define the contorl and model output files to compare
ctrl_outfile=${ctrl_output}/${rundate}_${fhr}
test_outfile=${TASK_OUTPUT}/model/${rundate}_${fhr}

# Inform the world about the location of the files that are compared
[[ ! -f ${ctrl_outfile} ]] && \
  echo -e "\n ERROR: control file not found: ${ctrl_outfile}" |tee -a ${gmtestinfo} && exit 1
[[ ! -f ${test_outfile} ]] && \
  echo -e "\n ERROR: run output not found: ${test_outfile}" |tee -a ${gmtestinfo} && exit 1
echo -e "\n ref. control: ${ctrl_outfile}" | tee -a ${gmtestinfo}
echo -e "\n test outfile: ${test_outfile}" | tee -a ${gmtestinfo}

# Compare the model output with the control files
complistfile=${TASK_BASEDIR}/fstcomp_listing
echo -e "\n fstcomp listing: ${complistfile}" | tee -a ${gmtestinfo}
${TASK_BIN}/fstcomp -a ${ctrl_outfile} -b ${test_outfile} > ${complistfile}

# Inform the world about the results of comparison
if [ -f ${complistfile} ] && [[ $(grep -i  "error" ${complistfile}) == "" ]] ; then
 if [[ $(grep '<' ${complistfile}) != "" ]] ; then
  echo -e "\n *** ERROR: New binary does not reproduce control test. ***\n" | tee -a ${gmtestinfo}
  exit 1
 elif [[ $(grep 'PAS TROUVE' ${complistfile}) != "" ]] ; then
  echo -e "\n *** ERROR: New binary does not reproduce control test. ***\n" | tee -a ${gmtestinfo}
  exit 1
 elif [[ $(grep "0.0000E+00  0.0000E+00  0.0000E+00" ${complistfile}) != "" ]] || \
      [[ $(grep "0.0000E+00 0.0000E+00 0.0000E+00" ${complistfile}) != "" ]] ; then
  echo -e "\n *** New binary reproduces reference output. ***\n" | tee -a ${gmtestinfo}
 else
  echo -e "\n *** ERROR: Comparison seem not to produce any results. ***\n" | tee -a ${gmtestinfo}
  exit 1
 fi
else
 echo -e "\n *** ERROR: Comparison of the output failed. ***\n" | tee -a ${gmtestinfo}
 exit 1
fi

if [[ "${ci_opt}" == "false" ]] ; then
 # Tell the world how long it took to validate GEM-MACH regression test results
 echo -e "\n== It took $(r.date -n -MM -L $(date '+%C%y%m%d%H%M%S') ${scriptstartdate}) seconds to run the validation script. == \n" | tee -a ${gmtestinfo}
fi

