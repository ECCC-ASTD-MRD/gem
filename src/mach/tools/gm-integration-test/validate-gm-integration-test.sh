#!/bin/bash

###
#
#   Script for for validating GEM-MACH integration test results.
#
#   USAGE: ord_soumet ${TASK_BIN}/validate-gm-integration-test.sh -args "${IntegrationTest_rpn_utils} ${IntegrationTest_version} ${control_dir} ${TASK_BASEDIR} ${gmtestinfo}" -mach ${GMJobMach} -cpus ${VldtJobProcTopo} -cm ${GMJobMemory} -t ${VldtJobTime} -queue ${GMJobQueue} -jn ${VldtJobName} -listing ${TASK_BASEDIR}/listing
#          or
#          validate-gm-integration-test.sh ${IntegrationTest_rpn_utils} ${IntegrationTest_version} ${control_dir} ${TASK_BASEDIR} ${gmtestinfo} 2>&1 > gm-test-${TRUE_HOST}-validation-listings.txt
#
#   If the test is successful, the script will print:
#   ***
#   Congratulations: New binary reproduces control test.
#   ***
#   on the screen or in the listings file.
#
#   The script relies on:
#       1. Existence of directory structure needed for Runmod task with root directory being placed at ${TASK_BASEDIR}
#       2. Existence of reference output file named ${control_dir}/gm-output/${IntegrationTest_version}_${TRUE_HOST}/model/${rundate}_${fhr}
#   If any of these are missing, or not up to date, the script will fail. The error messages are provided throughout the
#   script to indicate the obvious issues and failure of the script to behave as expected.
#
#   The script will:
#       - Compare results of the run with the previously saved results of the control run performed with the same input and a
#         binary built with the last released GEM-MACH code.
#       - Stop and provide the error message when there are issues with any of the described steps.
#       - Time the execution.
#
# Author: Verica Savic-Jovcic
# Date:   March 2022
# Update: January 2023
#
# 2023-Apr. Jack C - source GEM environment, minor tweak, simplify,
#         remove HARDCODD $rundate etc.
# 2023-May. Jack C - simplify to just fstcomp of output with ref./cntrl output
#
###

scriptstartdate=$(date '+%C%y%m%d%H%M%S')

# Read in arguments
GEM_version=$1
IntegrationTest_version=$2
control_dir=$3
TASK_BASEDIR=$4
gmtestinfo=$5
cmpl_opt=$6

# update log
echo -e "\n== Starting script: validate-gm-integration-test: $scriptstartdate == \n" | tee -a ${gmtestinfo}

# Check and load GEM environment
gembndl=/fs/ssm/eccc/mrd/rpn/models/gem/bundles/${GEM_version}.bndl
[[ ! -f $gembndl ]] && echo -d " ERROR: GEM environment file not available: $gembndl \n" && exit 1
source r.load.dot /fs/ssm/eccc/mrd/rpn/models/gem/bundles/${GEM_version}

# Define task directory structure
export TASK_BIN=${TASK_BASEDIR}/bin
export TASK_INPUT=${TASK_BASEDIR}/input
export TASK_WORK=${TASK_BASEDIR}/work
export TASK_OUTPUT=${TASK_BASEDIR}/output

cp -a $(which fstcomp) ${TASK_BIN}

# output fst files to compare
ctrl_output=${control_dir}/gm-output_${TRUE_HOST}/model
[[ "${cmpl_opt}" == "dbg" ]] && ctrl_output=${control_dir}/gm-output_${TRUE_HOST}_${cmpl_opt}/model
fhr=024
rundate=$(basename ${ctrl_output}/??????????_${fhr} |cut -c1-10)
ctrl_outfile=${ctrl_output}/${rundate}_${fhr}
test_outfile=${TASK_OUTPUT}/${rundate}_${fhr}

[[ ! -f ${ctrl_outfile} ]] && \
  echo -e "\n ERROR: control file not found: ${ctrl_outfile}" |tee -a ${gmtestinfo} && exit 1
[[ ! -f ${test_outfile} ]] && \
  echo -e "\n ERROR: run output not found: ${test_outfile}" |tee -a ${gmtestinfo} && exit 1
echo -e "\n ref. control: ${ctrl_outfile}" | tee -a ${gmtestinfo}
echo -e "\n test outfile: ${test_outfile}" | tee -a ${gmtestinfo}

# Compare the output with the existing resutls
complistfile=${TASK_BASEDIR}/fstcomp_listing
${TASK_BIN}/fstcomp -a ${ctrl_outfile} -b ${test_outfile} > ${complistfile}
if [ -f ${complistfile} ] && [[ $(grep -i  "error" ${complistfile}) == "" ]] ; then
 if [[ $(grep '<' ${complistfile}) != "" ]] ; then
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

# Tell the world how long it took to validate GEM-MACH integration test results
echo -e "\n== DONE with validation script == \n" | tee -a ${gmtestinfo}

