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
#       2. Existence of reference output file named ${control_dir}/gm-output/${IntegrationTest_version}_${TRUE_HOST}/model/${test_date}_${test_fhr}
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
#         remove HARDCODD $test_date etc.
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

ln -sf $(which editfst) ${TASK_BIN}

# Define the output hour to compare
test_fhr=024

# Define run output
TASK_OUTPUT=${TASK_BASEDIR}/output
case ${test_fhr} in
 000) test_step=072 ;;
 009) test_step=144 ;;
 012) test_step=144 ;;
 021) test_step=288 ;;
 024) test_step=288 ;;
esac
test_date=$(basename ${TASK_OUTPUT}/cfg_0000/laststep_0000000${test_step}/000-000/pm*-000-000_${test_fhr} |cut -c3-12)
dmoutfile=${TASK_OUTPUT}/cfg_0000/laststep_0000000${test_step}/000-000/dm${test_date}-000-000_${test_fhr}
pmoutfile=${TASK_OUTPUT}/cfg_0000/laststep_0000000${test_step}/000-000/pm${test_date}-000-000_${test_fhr}
outfile=${TASK_OUTPUT}/${test_date}_${test_fhr}
if [ ! -f ${dmoutfile} ] ; then
 echo -e "\n\n ERROR: GEM dynamics output file, ${dmoutfile}, is unavailable. \n\n" | tee -a ${gmtestinfo}
 exit 1
else
 echo -e "\n Location of GEM dynamics output:\n ${dmoutfile}" | tee -a ${gmtestinfo}
fi
if [ ! -f ${pmoutfile} ] ; then
 echo -e "\n\n ERROR: GEM physics output file, ${pmoutfile}, is unavailable. \n\n" | tee -a ${gmtestinfo}
 exit 1
else
 echo -e "\n Location of GEM physics output:\n ${pmoutfile}" | tee -a ${gmtestinfo}
fi

# combine pm,dm level outputs in one
for f in ${TASK_OUTPUT}/cfg_0000/laststep_0000000${test_step}/00*/?m${test_date}-*_${test_fhr} ; do
 ${TASK_BASEDIR}/bin/editfst -s $f -d ${outfile} -i <<EOD
 exclure(-1,['>>','^^','!!','>^'],-1)
EOD
done
${TASK_BIN}/editfst -s ${dmoutfile} -d ${outfile} -i <<EOD
 desire(-1,['>>','^^','!!','>^'],-1)
EOD
if [ ! -f ${outfile} ] ; then
 echo -e "\n\n ERROR: The output is not successfully combined \n\n" | tee -a ${gmtestinfo}
 exit 1
fi

# input reference control output
test_outfile=$outfile
ctrl_outfile=${control_dir}/gm-output_${TRUE_HOST}/model/${test_date}_${test_fhr}
[[ "${cmpl_opt}" == "dbg" ]] && \
 ctrl_outfile=${control_dir}/gm-output_${TRUE_HOST}_${cmpl_opt}/model/${test_date}_${test_fhr}
if [ ! -f ${ctrl_outfile} ] ; then
 echo -e "\n\n ERROR: control output file, ${ctrl_outfile}, is unavailable. \n\n" | tee -a ${gmtestinfo}
 exit 1
else
 echo -e "\n Location of ref. control: ${ctrl_outfile}" | tee -a ${gmtestinfo}
 echo -e "\n Location of test outfile: ${test_outfile}" | tee -a ${gmtestinfo}
fi

# Compare the output with the existing resutls
complistfile=${TASK_BASEDIR}/fstcomp_listing
fstcomp -a ${ctrl_outfile} -b ${test_outfile} > ${complistfile}
if [ -f ${complistfile} ] && [[ $(grep -i  "error" ${complistfile}) == "" ]] ; then
 if [[ $(grep '<' ${complistfile}) != "" ]] ; then
  echo -e "\n\n ***\n ERROR: New binary does not reproduce control test.\n ***\n\n" | tee -a ${gmtestinfo}
  exit 1
 elif [[ $(grep "0.0000E+00  0.0000E+00  0.0000E+00" ${complistfile}) != "" ]] || \
      [[ $(grep "0.0000E+00 0.0000E+00 0.0000E+00" ${complistfile}) != "" ]] ; then
  echo -e "\n\n ***\n Congratulations: New binary reproduces control test.\n ***\n\n" | tee -a ${gmtestinfo}
 else
  echo -e "\n\n ***\n ERROR: Comparison seem not to produce any results.\n ***\n\n" | tee -a ${gmtestinfo}
  exit 1
 fi
else
 echo -e "\n\n ***\n ERROR: Comparison of the output failed.\n ***\n\n" | tee -a ${gmtestinfo}
 exit 1
fi

# Tell the world how long it took to validate GEM-MACH integration test results
echo -e "\n== It took $(r.date -n -MM -L $(date '+%C%y%m%d%H%M%S') ${scriptstartdate}) seconds for the validation script to run. == \n" | tee -a ${gmtestinfo}

