#!/bin/bash

###
#
#   Script for running the GEM-MACH, combining the output and submitting  the post-processing script
#        for validation of the integration test results. The combined output can be used as the new control output.
#
#   USAGE: 
#        ord_soumet ${TASK_BIN}/run-gm-integration-test.sh -args "${IntegrationTest_version} ${control_dir} ${TASK_BASEDIR} ${GMJobMach} ${GMJobMemory} ${GMJobQueue} ${GMJobProcTopo} ${gmtestinfo} ${cmpl_opt} ${cntrl_fl_opt}" -mach ${GMJobMach} -cpus ${GMJobProcTopo} -cm ${GMJobMemory} -t ${GMJobTime} -mpi 1 -queue ${GMJobQueue} -jn ${GMJobName} -listing ${TASK_BASEDIR}/listings
#
#        Note that the variables provided under quotes after "-args" will be
#        read by this script as $1, $2, etc to fill in the values of apropriate variables
#
#   The script relies on:
#       1. Existence of directory structure needed for Runmod task with root directory being placed at ${TASK_BASEDIR}
#       2. Existence of GEM super repository in ${TASK_BASEDIR}/GEM-MACH directory
#       3. Existence of GEM-MACH binary and validate-gm-integration-test.sh scripts to be in ${TASK_BIN} directory
#       4. Existence of input files to be in ${TASK_INPUT} directory
#   If any of these are missing, or not up to date, the script will fail. The error messages are provided throughout the
#   script to indicate the obvious issues and failure of the script to behave as expected.
#
#   Script:
#     - Submits the model run
#     - Following the succesfful model run, combines the output in the fst files named as YYYYMMDDHH_hhh, 
#       which can be used as new control output.
#     - Based on value of the ${cntrl_fl_opt} variable submits the script for validation of the model results 
#       against the existing control output
#
#   Note:
#     Script passes ${cmpl_opt} variable to the validation script to ensure suitability of validation.
#     I.e., that the validation of the debug mode run is done against the debug control output.
#
# Author: Verica Savic-Jovcic
# Date:   March 2022
# Update: January 2023
#
#         April 2023, Jack C - update to run with cmake compilation.  Note that there is
#          no version check, and this is not compatible with RDE
#          (only work with cmake, after >GEM5.2_b2)
#         May 2023, Jack C - combine dm/pm output files after runmod finishes, and calls
#          the validate script when "${cntrl_fl_opt}" != "new"
#         August 2023, Verica S-J - update documentation
#
###

scriptstartdate=$(date '+%C%y%m%d%H%M%S')

# Read in arguments
IntegrationTest_version=$1
control_dir=$2
export TASK_BASEDIR=$3
GMJobMach=$4
GMJobMemory=$5
GMJobQueue=$6
GMJobProcTopo=$7
gmtestinfo=$8
cmpl_opt=$9
cntrl_fl_opt=${10}

# Update log
echo -e "\n== Starting script: run-gm-integration-test: $scriptstartdate == \n" | tee -a ${gmtestinfo}

# Load needed environment (instead of GEM environment)
source ${TASK_BASEDIR}/GEM-MACH/.eccc_setup_intel

# Define task directory structure
export TASK_BIN=${TASK_BASEDIR}/bin
export TASK_INPUT=${TASK_BASEDIR}/input
export TASK_WORK=${TASK_BASEDIR}/work
export TASK_OUTPUT=${TASK_BASEDIR}/output

# Submit the model run
cd ${TASK_BASEDIR}
. r.call.dot ${TASK_BIN}/runmod.sh -task_basedir ${TASK_BASEDIR} -ptopo ${GMJobProcTopo} \
            -smt 0x0 -inorder 1 -barrier -timing 0 -cfg 0:0 -no_setup -debug 0
if [ "$_status" == "ABORT" ]; then
 echo -e "\n\n ERROR: There is an issue with the execution of runmod.sh. \n\n" | tee -a ${gmtestinfo}
 exit 1
fi

# Combine physics and dynamics outputs
if [ "${cmpl_opt}" == "dbg" ] ; then 
 fhrs=(-003 000 009)
else
 fhrs=(-003 000 009 012 021 024)
fi
for fhr in ${fhrs[*]} ; do
 run_step=$(basename $(dirname $(dirname ${TASK_OUTPUT}/cfg_0000/laststep_0000000???/000-000/pm??????????-000-000_${fhr})))
 rundate=$(basename ${TASK_OUTPUT}/cfg_0000/${run_step}/000-000/pm??????????-000-000_${fhr} | cut -c3-12)
 dmoutfile=${TASK_OUTPUT}/cfg_0000/${run_step}/000-000/dm${rundate}-000-000_${fhr}
 pmoutfile=${TASK_OUTPUT}/cfg_0000/${run_step}/000-000/pm${rundate}-000-000_${fhr}
 outfile=${TASK_OUTPUT}/${rundate}_${fhr}
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
 # Combine pm and dm outputs in one
 for f in ${TASK_OUTPUT}/cfg_0000/${run_step}/00*/?m${rundate}-*_${fhr} ; do
  ${TASK_BIN}/editfst -s $f -d ${outfile} -i <<EOD
 exclure(-1,['>>','^^','!!','>^'],-1)
EOD
 done
 ${TASK_BIN}/editfst -s ${dmoutfile} -d ${outfile} -i <<EOD
 desire(-1,['>>','^^','!!','>^'],-1)
EOD
 echo -e "*** combined output: ${outfile} \n"  | tee -a ${gmtestinfo}
done

# Set the validation job resources
export PostProcJobTime=1200
export PostProcJobProcTopo=40

# Submit job to verify results against control files
if [[ "${cntrl_fl_opt}" != "new" ]]; then
   export PostProcJobName=vldtgm
   echo -e "\n Validate results against a control run \n" | tee -a ${gmtestinfo}
   ord_soumet ${TASK_BIN}/validate-gm-integration-test.sh \
              -args "${IntegrationTest_version} ${control_dir} ${TASK_BASEDIR} ${gmtestinfo} ${cmpl_opt}" \
              -mach ${GMJobMach} -cpus ${PostProcJobProcTopo} -cm ${GMJobMemory} \
              -t ${PostProcJobTime} -queue ${GMJobQueue} -jn ${PostProcJobName} \
              -listing ${TASK_BASEDIR}/listings
else
  echo -e "\n *** Location of new control output: $(dirname ${outfile}) *** \n " | tee -a ${gmtestinfo}
fi

# Tell the world how long it took to run GEM-MACH
echo -e "\n== It took $(r.date -n -MM -L $(date '+%C%y%m%d%H%M%S') ${scriptstartdate}) seconds for the script that runs GEM-MACH and combines the output to run == \n" | tee -a ${gmtestinfo}
