#!/bin/bash
set -e

###
#
#  Script for compiling GEM-MACH as a part of GEM-MACH integration test
#
#  Usage: ord_soumet ${TASK_BIN}/compile-gm-for-integration-test.sh -args "${GEM_version} ${IntegrationTest_version} ${control_dir} ${TASK_BASEDIR} ${GMJobMach} ${GMJobMemory} ${GMJobQueue} ${GMJobTopo} ${cmpl_opt} ${cntrl_fl_opt}" -mach ${GMJobMach} -cpus ${CmplJobProcTopo} -cm ${GMJobMemory} -t ${CmplJobTime} -mpi 1 -queue ${GMJobQueue} -jn ${CmplJobName} -listing ${TASK_BASEDIR}/listing
#
#   The script relies on:
#       1. Existence of directory structure needed for Runmod task with root directory being placed at ${TASK_BASEDIR}
#       2. Availability of requested GEM version
#   If any of these are missing, or not up to date, the script will fail. The error messages are provided throughout the
#   script to indicate the obvious issues and failure of the script to behave as expected.
#
#  Script:
#    1. Loads GEM environment for GEM version provided as an argument
#    2. Based on the value of ${cmpl_opt}:
#       - Compiles GEM-MACH code available in ${TASK_BASEDIR}/GEM-MACH directory in debug mode or in usually utilized optimized mode
#       - Sets &physics_cfgs/debug_trace_L and &chemistry_cfgs/chm_debug_trace_l namelist keys to .true. for the test that is set to run in the debug mode
#    3. Copies the binary in ${TASK_BASEDIR}/bin directory
#    4. Copies GEM scripts needed for running GEM-MACH integration test into ${TASK_BASEDIR}/bin directory
#    5. Submits the GEM-MACH integration test
#    6. Times the execution and saves it in the listings
#
#  Note:
#    When compiling in the debug mode, the script modifies the gem_settings.nml to set physics and chemistry debug keys.
#
# Author: Verica Savic-Jovcic
# Date:   March 2022
# Update: January 2023
# Update: by Jack Chen and Verica Savic-Jovcic, March 2023
#
# 2023-Apr/May Jack C - update to compile using cmake with mach library as part
#   of MIG's GEM git subtree.
#   Note there is no version check, thus will not work with RDE.
#
###
scriptstartdate=$(date '+%C%y%m%d%H%M%S')

# Read in arguments
GEM_version=$1
IntegrationTest_version=$2
control_dir=$3
TASK_BASEDIR=$4
GMJobMach=$5
GMJobMemory=$6
GMJobQueue=$7
GMJobTopo=$8
gmtestinfo=$9
cmpl_opt=${10}
cntrl_fl_opt=${11}

# update log
echo -e "\n== Strating script: compile-gm-for-integration-test. ${scriptstartdate}  == \n" | tee -a ${gmtestinfo}

# Compile GEM-MACH with cmake (with system RPN library)
echo -e "\n GEM-MACH compilation location: ${TASK_BASEDIR}/build \n"
cmake_dir=${TASK_BASEDIR}/build
cd ${cmake_dir}
source ${cmake_dir}/.eccc_setup_intel
source ${cmake_dir}/.initial_setup

# Tell the world about GEM environment and compiler
cat << EOF | tee -a ${gmtestinfo}
GEM-MACH compilation location: ${cmake_dir}
ORDENV_DIST: ${ORDENV_DIST}
COMPILER_SUITE: ${COMPILER_SUITE}
COMPILER_VERSION: ${COMPILER_VERSION}
GEM_ARCH: ${GEM_ARCH}
ATM_MODEL_VERSION: ${ATM_MODEL_VERSION}
gemdyn_version: ${gemdyn_version}
rpnphy_version: ${rpnphy_version}
modelutils_version: ${modelutils_version}
mach_version: $(cat ${cmake_dir}/src/mach/MANIFEST | sed -n '/VERSION/s/.*: //p')

compile build directory: ${cmake_dir}/build-${GEM_ARCH}
link binary: ${cmake_dir}/work-${GEM_ARCH}/bin/maingemdm

EOF

# build Makefile, and compile/link binary
if [[ "${cmpl_opt}" == "dbg" ]] ; then
   (time make VERBOSE=1 cmake-mach-debug) |& tee ${cmake_dir}/make.cmake-mach-debug.out
   (time make VERBOSE=1 -j work) |& tee make.work.out
else
   (time make cmake-mach) |& tee ${cmake_dir}/make.cmake-mach.out
   (time make -j work) |& tee make.work.out
fi

### link compiled binary to TASK_BASEDIR/bin
TASK_BIN=${TASK_BASEDIR}/bin
gemmach_abs=${cmake_dir}/work-${GEM_ARCH}/bin/maingemdm
[[ ! -f ${gemmach_abs} ]] && echo -e "\nERROR: GEM-MACH binary not available $gemmach_abs \n" && exit 1
echo -e "\n == GEM-MACH compiled at $(date) == \n"
ln -s ${gemmach_abs} ${TASK_BIN}/ATM_MOD.Abs
echo -e "\nGEM-MACH binary is copied to ${TASK_BIN}" | tee -a ${gmtestinfo}
echo -e " GEM_ovbin=${TASK_BIN} \n\n"

# Set the job resources
export GMJobTime=1200
export GMJobProcTopo=$GMJobTopo
export GMJobName=rungm
[[ "${cmpl_opt}" == "dbg" ]] && export GMJobTime=2400

# Submit GEM-MACH integration test to run
echo -e "\n == Submit the GEM-MACH integration-test run at $(date) == \n"
ord_soumet ${TASK_BIN}/run-gm-integration-test.sh \
           -args "${GEM_version} ${IntegrationTest_version} ${control_dir} ${TASK_BASEDIR} ${GMJobMach} ${GMJobMemory} ${GMJobQueue} ${GMJobProcTopo} ${gmtestinfo} ${cmpl_opt} ${cntrl_fl_opt}" \
           -mach ${GMJobMach} -cpus ${GMJobProcTopo} -cm ${GMJobMemory} -t ${GMJobTime} \
           -mpi 1 -queue ${GMJobQueue} -jn ${GMJobName} -listing ${TASK_BASEDIR}/listing

# Update status
echo -e "\n== It took $(r.date -n -MM -L $(date '+%C%y%m%d%H%M%S') ${scriptstartdate}) seconds for the compilation script to run. == \n" | tee -a ${gmtestinfo}

