#!/bin/bash

###
#
# Script for submitting the runmod.sh script with wall-clock time predefined
#        in the ord_soumet call for this script.
#
# To call the script, set
#        - the locationss of rungm.sh script (${TASK_BIN}/rungm.sh), script's working (${TASK_BASEDIR}) 
#          and listings (${TASK_BASEDIR}/listing) directories,
#        - GEM version (${GEM_version}),
#        - processor topoplogy ($GMJobProcTopo}),
#        - computer name that will execute the script ($GMJobMach}),
#        - amount of memory each processor will use, including the units (e.g. 2G), ($GMJob}),
#        - wall-clock time ($GMJobTime}),
#        - queue ($GMJobQueue}),
#        - name of the job that will be used to make the listing and to show in the queue ($GMJob}),
#        and then type in the command line, or use in the script the following:
#
#        ord_soumet ${TASK_BIN}/rungm.sh -args "${GEM_version} ${GMJobProcTopo} ${TASK_BASEDIR}" -mach ${GMJobMach} -cpus ${GMJobProcTopo} -cm ${GMJobMemory} -t ${GMJobTime} -mpi 1 -queue ${GMJobQueue} -jn ${GMJobName} -listing ${TASK_BASEDIR}/listing
#
#        Note that the variables provided under quotes after "-args" will be 
#        read by this script as $1, $2 and $3 to fill in the values of apropriate variables
#
# Author: Verica Savic-Jovcic
# Date:   July 2021
#
###

# Set GEM environment
export GEM_version=$1
. r.load.dot eccc/mrd/rpn/MIG/GEM/${GEM_version}

# Set the task's processor topology
GMJobProcTopo=$2

# Set task directory structure
export TASK_BASEDIR=$3
export TASK_BIN=${TASK_BASEDIR}/bin
export TASK_INPUT=${TASK_BASEDIR}/input
export TASK_WORK=${TASK_BASEDIR}/work
export TASK_OUTPUT=${TASK_BASEDIR}/output

# Initiate the task
cd ${TASK_BASEDIR}
. r.call.dot ${TASK_BIN}/runmod.sh -task_basedir ${TASK_BASEDIR} -ptopo ${GMJobProcTopo} -smt 0x0 -inorder 1 -barrier -timing 0 -cfg 0:0 -no_setup -debug 0
if [ "$_status" == "ABORT" ]; then 
 echo "ERROR: There is an issue with the execution of runmod.sh."
 exit 1
fi

