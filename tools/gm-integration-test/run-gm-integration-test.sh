#!/bin/bash

###
#
#   Script for running the GEM-MACH integration test and submitting one of the post-processing scripts:
#        for validation of the integration test results, or for creation of the new control output.
#
#   USAGE: ord_soumet ${TASK_BIN}/run-gm-integration-test.sh -args "${GEM_version} ${IntegrationTest_version} ${control_dir} ${TASK_BASEDIR} ${GMJobMach} ${GMJobMemory} ${GMJobQueue} ${GMJobProcTopo} ${gmtestinfo} ${cmpl_opt} ${cntrl_fl_opt}" -mach ${GMJobMach} -cpus ${GMJobProcTopo} -cm ${GMJobMemory} -t ${GMJobTime} -mpi 1 -queue ${GMJobQueue} -jn ${GMJobName} -listing ${TASK_BASEDIR}/listing
#
#        Note that the variables provided under quotes after "-args" will be
#        read by this script as $1, $2, etc to fill in the values of apropriate variables
#
#   The script relies on:
#       1. Existence of directory structure needed for Runmod task with root directory being placed at ${TASK_BASEDIR}
#       2. Existence of GEM-MACH binary, and run-gm-integration-test.sh and validate-gm-integration-test.sh scripts to be in ${TASK_BIN} directory
#       3. Existence of input files to be in ${TASK_INPUT} directory
#   If any of these are missing, or not up to date, the script will fail. The error messages are provided throughout the
#   script to indicate the obvious issues and failure of the script to behave as expected.
#
#   Script:
#     - Submits the model run
#     - Following the succesfful model run and based on value of the ${cntrl_fl_opt} variable submits
#       either the script for validation of the model results against the existing control output,
#       or the script for preparing new control output.
#
#   Note:
#     Script passes ${cmpl_opt} variable to the validation script to ensure suitability of validation.
#     I.e., that the validation of the debug mode run is done against the debug control output.
#
# Author: Verica Savic-Jovcic
# Date:   March 2022
# Update: January 2023
#
# 2023-Apr. Jack C - update to run with cmake compilation.  Note that there is
#   no version check, and this is not compatible with RDE
#   (only work with cmake, after >GEM5.2_b2)
###

scriptstartdate=$(date '+%C%y%m%d%H%M%S')

# Read in arguments
export GEM_version=$1
IntegrationTest_version=$2
control_dir=$3
export TASK_BASEDIR=$4
GMJobMach=$5
GMJobMemory=$6
GMJobQueue=$7
GMJobProcTopo=$8
gmtestinfo=$9
cmpl_opt=${10}
cntrl_fl_opt=${11}

# update log
echo -e "\n== Starting script: run-gm-integration-test: $scriptstartdate == \n" | tee -a ${gmtestinfo}

# Check and load GEM environment
gembndl=/fs/ssm/eccc/mrd/rpn/models/gem/bundles/${GEM_version}.bndl
[[ ! -f $gembndl ]] && echo -d " ERROR: GEM environment file not available: $gembndl \n" && exit 1
source r.load.dot /fs/ssm/eccc/mrd/rpn/models/gem/bundles/${GEM_version}

# Define task directory structure
export TASK_BIN=${TASK_BASEDIR}/bin
export TASK_INPUT=${TASK_BASEDIR}/input
export TASK_WORK=${TASK_BASEDIR}/work
export TASK_OUTPUT=${TASK_BASEDIR}/output

### Copy the scripts necessary for running GEM-MACH into the bin directory
cp $(which r.run_in_parallel) ${TASK_BIN} ; ln -s ${TASK_BIN}/r.run_in_parallel ${TASK_BIN}/r.mpirun
cp $(which rungem.sh) ${TASK_BIN}
cp $(which runmod.sh) ${TASK_BIN}

# Submit the model run
cd ${TASK_BASEDIR}
. r.call.dot ${TASK_BIN}/runmod.sh -task_basedir ${TASK_BASEDIR} -ptopo ${GMJobProcTopo} -smt 0x0 -inorder 1 -barrier -timing 0 -cfg 0:0 -no_setup -debug 0
if [ "$_status" == "ABORT" ]; then
 echo -e "\n\n ERROR: There is an issue with the execution of runmod.sh. \n\n" | tee -a ${gmtestinfo}
 exit 1
fi

# Set the job resources
export PostProcJobTime=1200
export PostProcJobProcTopo=40

# Submit post processing job for model verification against the control run, or to prepare a new control-run results

if [[ "${cntrl_fl_opt}" == "new" ]]; then
   export PostProcJobName=prepgm
   echo -e "\n Submit post-processing job to ${GMJobMach} to prepare new control output \n"
   ord_soumet ${TASK_BIN}/prepare-new-gm-integration-test.sh \
              -args "${GEM_version} ${TASK_BASEDIR} ${gmtestinfo}" \
              -mach ${GMJobMach} -cpus ${PostProcJobProcTopo} -cm ${GMJobMemory} \
              -t ${PostProcJobTime} -queue ${GMJobQueue} -jn ${PostProcJobName} \
              -listing ${TASK_BASEDIR}/listing
else
   export PostProcJobName=vldtgm
   echo -e "\n Submit post-processing job to ${GMJobMach} to validate the model against a control run \n"
   ord_soumet ${TASK_BIN}/validate-gm-integration-test.sh \
              -args "${GEM_version} ${IntegrationTest_version} ${control_dir} ${TASK_BASEDIR} ${gmtestinfo} ${cmpl_opt}" \
              -mach ${GMJobMach} -cpus ${PostProcJobProcTopo} -cm ${GMJobMemory} \
              -t ${PostProcJobTime} -queue ${GMJobQueue} -jn ${PostProcJobName} \
              -listing ${TASK_BASEDIR}/listing
fi

# Tell the world how long it took to run GEM-MACH
echo -e "\n== It took $(r.date -n -MM -L $(date '+%C%y%m%d%H%M%S') ${scriptstartdate}) seconds for the GEM-MACH integration test to run on ${GMJobProcTopo} processors. == \n" | tee -a ${gmtestinfo}

