#!/bin/bash

###
#
#  Script for compiling GEM-MACH as a part of GEM-MACH regression test
#      and for preparing all necessary scripts for running GEM-MACH.
#
#  Usage:
#      ord_soumet ${TASK_BIN}/compile-gm.sh -args "${TASK_BASEDIR}" -mach ${Test_JobsMach} -cpus ${CmplJobNcpu} -cm ${Test_JobsMemory} -t ${CmplJobTime} -mpi 1 -queue ${Test_JobsQueue} -jn ${CmplJobName} -listing ${TASK_LIST}
#
#  The script relies on:
#      1. Existence of directory structure needed for Runmod task with 
#         root directory being placed at ${TASK_BASEDIR}
#      2. Availability of GEM super repository in the 
#         ${TASK_BASEDIR}/GEM-MACH directory
#      3. Existence of configuration file named 
#         ${TASK_BASEDIR}/compile_env_arg
#  If any of these are missing, or not up to date, the script will fail. 
#  The error messages are provided throughout the script to indicate any
#  obvious issues and failure of the script to behave as expected.
#
#  Script:
#      1. Loads environment provided in the GEM super repository
#      2. Based on the value of ${cmpl_opt}:
#         - Compiles GEM-MACH code available in ${TASK_BASEDIR}/GEM-MACH 
#           directory in debug mode or in usually utilized optimized mode
#         - Sets &physics_cfgs/debug_trace_L and 
#           &chemistry_cfgs/chm_debug_trace_l namelist keys to .true. for
#           the test that is set to run in the debug mode
#      3. Links the binary to ${TASK_BASEDIR}/bin directory
#      4. Copies GEM scripts needed for running GEM-MACH executable into 
#         ${TASK_BASEDIR}/bin directory
#      5. Submits the GEM-MACH run
#      6. Times the execution and saves it in the listings
#
#  Note:
#    When compiling in the debug mode, the script modifies 
#    gem_settings.nml to set physics and chemistry debug keys.
#
# Author: Verica Savic-Jovcic
# Date:   March 2022
#
#         Jack Chen and Balbir Pabla contributed to updating this script.
#
###
scriptstartdate=$(date '+%C%y%m%d%H%M%S')

# Read in arguments
TASK_BASEDIR=$1
source ${TASK_BASEDIR}/compile_env_arg

# Update log
echo -e "\n== Strating compile-gm script: ${scriptstartdate}  == \n" | tee -a ${gmtestinfo}

# Compile GEM-MACH with cmake (with system RPN library)
echo -e "\n GEM-MACH compilation location: ${TASK_BASEDIR}/GEM-MACH \n"
cmake_dir=${TASK_BASEDIR}/GEM-MACH
cd ${cmake_dir}
#source ${cmake_dir}/.eccc_setup_intel_2025.1.0 -> Commented out because it is done through common environment
source ${cmake_dir}/.initial_setup
${cmake_dir}/scripts/link-dbase.sh

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

# Build Makefile, and compile/link binary
if [[ "${cmpl_opt}" == "dbg" ]] ; then
   (time make VERBOSE=1 cmake-mach-debug) |& tee ${cmake_dir}/make.cmake-mach-debug.out
   (time make VERBOSE=1 -j work) |& tee ${cmake_dir}/make.work.out
elif [[ "${cmpl_opt}" == "dbe" ]] ; then
   (time make VERBOSE=1 cmake-mach-debug-extra) |& tee ${cmake_dir}/make.cmake-mach-debug-extra.out
   (time make VERBOSE=1 -j work) |& tee ${cmake_dir}/make.work.out
elif [[ "${cmpl_opt}" == "str" ]] ; then
   (time make VERBOSE=1 cmake-mach-strict) |& tee ${cmake_dir}/make.cmake-mach-strict.out
   (time make VERBOSE=1 -j work) |& tee ${cmake_dir}/make.work.out
elif [[ "${cmpl_opt}" == "pkg" ]] ; then
   (time make cmake-mach-static) |& tee ${cmake_dir}/make.cmake-mach-static.out
   (time make -j work) |& tee ${cmake_dir}/make.work.out
   (time make package) |& tee ${cmake_dir}/make.package.out
else
   (time make cmake-mach) |& tee ${cmake_dir}/make.cmake-mach.out
   (time make -j work) |& tee ${cmake_dir}/make.work.out
fi

### Link compiled binary and copy all supporting programs to the bin directory
gemmach_abs=${cmake_dir}/work-${GEM_ARCH}/bin/maingemdm
[[ ! -f ${gemmach_abs} ]] && echo -e "\nERROR: GEM-MACH binary not available $gemmach_abs \n" && exit 1
echo -e "\n == GEM-MACH compiled at $(date) == \n"
ln -s ${gemmach_abs} ${TASK_BIN}/ATM_MOD.Abs
echo -e "\nGEM-MACH binary is copied to ${TASK_BIN}" | tee -a ${gmtestinfo}
echo -e " GEM_ovbin=${TASK_BIN} \n\n"

### Copy the scripts necessary for running GEM-MACH into the bin directory
gemmach_bin=${cmake_dir}/work-${GEM_ARCH}/bin
cp -r ${gemmach_bin}/* ${TASK_BIN}/
cp $(which r.run_in_parallel) ${TASK_BIN} ; ln -s ${TASK_BIN}/r.run_in_parallel ${TASK_BIN}/r.mpirun
cp $(which rungem.sh) ${TASK_BIN}
cp $(which runmod.sh) ${TASK_BIN}
cp $(which editfst) ${TASK_BIN}

if [[ ${sequence_opt} == true ]] ;then
 # Submit GEM-MACH run
 echo -e "\n == Submit the GEM-MACH run at $(date) == \n"
 ord_soumet ${TASK_BIN}/run-gm.sh \
           -args "${TASK_BASEDIR}" \
           -mach ${Test_JobsMach} -cpus ${GMRunJobProcTopo} -cm ${Test_JobsMemory} -t ${GMRunJobTime} \
           -mpi 1 -queue ${Test_JobsQueue} -jn ${GMRunJobName} -listing ${TASK_LIST}

 # Update status
 echo -e "\n== It took $(r.date -n -MM -L $(date '+%C%y%m%d%H%M%S') ${scriptstartdate}) seconds to run the compilation script. == \n" | tee -a ${gmtestinfo}
fi
