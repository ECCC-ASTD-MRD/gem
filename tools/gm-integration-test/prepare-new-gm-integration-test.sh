#!/bin/bash

###
#
#   Script for for preparing new control output for the GEM-MACH integration test.
#
#   USAGE from the root of GEM-MACH git repository:
#          ord_soumet TOOLS/gm-integration-test/prepare-new-gm-integration-test.sh -args "${IntegrationTest_rpn_utils} ${TASK_BASEDIR} ${gmtestinfo}" -mach ${GMJobMach} -cpus ${NUTOJobProcTopo} -cm ${GMJobMemory} -t ${NUTOJobTime} -queue ${GMJobQueue} -jn ${NUTOJobName} -listing ${TASK_BASEDIR}/listing
#          or
#          TOOLS/gm-integration-test/prepare-new-gm-integration-test.sh ${IntegrationTest_rpn_utils} ${TASK_BASEDIR} ${gmtestinfo} 2>&1 > gm-test-${TRUE_HOST}-validation-listings.txt
#
#   If the test is successful, the script will produce the output for 000, 009, 012, 021 and 024 forecast hours.
#
#   The script relies on:
#       1. Existence of directory structure needed for Runmod task with root directory being placed at ${TASK_BASEDIR}
#       2. Existence of the output files in ${TASK_BASEDIR}/output directory
#   If any of these are missing the script will fail. The error messages are provided throughout the
#   script to indicate the obvious issues and failure of the script to behave as expected.
#
#   The script will:
#       - Combine the output of run-gm-integration-test.sh script (i.e. pm${test_date}-000-000_${test_fhr} and dm${test_date}-000-000_${test_fhr} files)
#         into {test_date}_${test_fhr} files in ${TASK_BASEDIR}/new_control/model directory.
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
TASK_BASEDIR=$2
gmtestinfo=$3

# update log
echo -e "\n== Starting script: prepare-new-gm-integration-test: $scriptstartdate == \n" | tee -a ${gmtestinfo}

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

# Create new control output
NewCntrlDir=${TASK_BASEDIR}/new_control/model
[[ ! -d ${NewCntrlDir} ]] && mkdir -p ${NewCntrlDir}

for test_fhr in 000 009 012 021 024 ; do
 case ${test_fhr} in
  000) test_step=072 ;;
  009) test_step=144 ;;
  012) test_step=144 ;;
  021) test_step=288 ;;
  024) test_step=288 ;;
 esac
 test_date=$(basename ${TASK_OUTPUT}/cfg_0000/laststep_0000000${test_step}/000-000/pm*-000-000_${test_fhr} |cut -c3-12)
 ctrl_outfile=${NewCntrlDir}/${test_date}_${test_fhr}

 # Define new control output configuration
 dmoutfile=${TASK_OUTPUT}/cfg_0000/laststep_0000000${test_step}/000-000/dm${test_date}-000-000_${test_fhr}
 pmoutfile=${TASK_OUTPUT}/cfg_0000/laststep_0000000${test_step}/000-000/pm${test_date}-000-000_${test_fhr}
 outfile=${TASK_OUTPUT}/${test_date}_${test_fhr}
 if [ ! -f ${dmoutfile} ] ; then
  echo -e "\n\n ERROR: GEM dynamics output file, ${dmoutfile}, is unavailable. \n\n" | tee -a ${gmtestinfo}
  exit 1
 else
  echo -e "\n Location of GEM dynamics output: ${dmoutfile} \n" | tee -a ${gmtestinfo}
 fi
 if [ ! -f ${pmoutfile} ] ; then
  echo -e "\n\n ERROR: GEM physics output file, ${pmoutfile}, is unavailable. \n\n" | tee -a ${gmtestinfo}
  exit 1
 else
  echo -e "\n Location of GEM physics output: ${pmoutfile} \n" | tee -a ${gmtestinfo}
 fi

 # combine pm,dm level outputs in one
 for f in ${TASK_OUTPUT}/cfg_0000/laststep_0000000${test_step}/00*/?m${test_date}-*_${test_fhr} ; do
  ${TASK_BIN}/editfst -s $f -d ${outfile} -i <<EOD
 exclure(-1,['>>','^^','!!','>^'],-1)
EOD
 done
 ${TASK_BIN}/editfst -s ${dmoutfile} -d ${outfile} -i <<EOD
 desire(-1,['>>','^^','!!','>^'],-1)
EOD
 if [ ! -f ${outfile} ] ; then
  echo -e "\n\n ERROR: The output is not successfully combined. \n\n" | tee -a ${gmtestinfo}
  exit 1
 else
  mv ${outfile} ${ctrl_outfile}
  echo "***" | tee -a ${gmtestinfo}
  echo "Location of new control output: ${ctrl_outfile}" | tee -a ${gmtestinfo}
 fi
done

# Tell the world how long it took to prepare the control output for GEM-MACH integration test
echo -e "\n== It took $(r.date -n -MM -L $(date '+%C%y%m%d%H%M%S') ${scriptstartdate}) seconds for the script to run. == \n"

