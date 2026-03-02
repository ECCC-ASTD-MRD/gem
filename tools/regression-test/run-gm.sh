#!/bin/bash

###
#
#  Script for running GEM-MACH and submitting the script
#      for assembling of the regression test results.
#
#  USAGE: 
#      ord_soumet ${TASK_BIN}/run-gm.sh -args "${TASK_BASEDIR}" -mach ${Test_JobsMach} -cpus ${GMRunJobProcTopo} -cm ${GMRunJobMemory} -t ${GMRunJobTime} -mpi 1 -queue ${Test_JobsQueue} -jn ${GMRunJobName} -listing ${TASK_LIST}
#
#      Note that the variables provided in quotes to ord_soumet after
#      "-args" will be read by this script, not the ord_soumet.
#
#  The script relies on:
#      1. Existence of directory structure needed for Runmod task with
#         root directory being placed at ${TASK_BASEDIR}
#      2. Existence of GEM super repository in ${TASK_BASEDIR}/GEM-MACH 
#         directory
#      3. Existence of GEM-MACH binary and ${TASK_BIN}/assemble-output.sh 
#         ${TASK_BIN}/assemble-output.sh script.
#      4. Existence of input files to be in ${TASK_INPUT} directory
#      5. Existence of configuration file named 
#         ${TASK_BASEDIR}/run_env_arg
#  If any of these are missing, or not up to date, the script will fail. 
#  The error messages are provided throughout the script to indicate any 
#  obvious issues and failure of the script to behave as expected.
#
#  Script:
#      - Submits the model run
#      - Following the succesfful model run, submits the script that 
#        assembles the output in one fst file per requested forecast hour.
#
# Author: Verica Savic-Jovcic
# Date:   March 2022
#
#         Jack Chen and Balbir Pabla contributed to the updates of this script. 
#
###

scriptstartdate=$(date '+%C%y%m%d%H%M%S')

# Read in arguments
export TASK_BASEDIR=$1

# Set up environment
source ${TASK_BASEDIR}/run_env_arg
cd ${TASK_WORK}

# Update log
echo -e "\n== Starting run-gm script: $scriptstartdate == \n" | tee -a ${gmtestinfo}

# Submit the model run
. r.call.dot ${TASK_BIN}/runmod.sh -task_basedir ${TASK_BASEDIR} -ptopo ${GMRunJobProcTopo} \
            -smt 0x0 -inorder 1 -barrier -timing 0 -cfg 0:0 -no_setup -debug 0
if [ "$_status" == "ABORT" ]; then
 echo -e "\n\n ERROR: There is an issue with the execution of runmod.sh. \n\n" | tee -a ${gmtestinfo}
 exit 1
fi

if [[ ${sequence_opt} == "true" ]]; then
 # Submit job to assemble outputs to one file per forecast hour
 echo -e "\n Assemble output to one file per requested forecast hour \n" | tee -a ${gmtestinfo}
 ord_soumet ${TASK_BIN}/assemble-output.sh -args "${TASK_BASEDIR}" \
            -mach ${Test_JobsMach} -cpus ${AssmblJobNcpu} -cm ${AssmblJobMemory} \
            -t ${AssmblJobTime} -queue ${Test_JobsQueue} -jn ${AssmblJobName} \
            -listing ${TASK_LIST}

 # Tell the world how long it took to run GEM-MACH
 echo -e "\n== It took $(r.date -n -MM -L $(date '+%C%y%m%d%H%M%S') ${scriptstartdate}) seconds to run the script that runs GEM-MACH. == \n" | tee -a ${gmtestinfo}
fi
