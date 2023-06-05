#!/bin/bash
set -e

###
#
#   Script for initiating GEM-MACH integration test from the root directory of the GEM-MACH repository.
#
#   USAGE: While in the root directory of GEM-MACH repository, type:
#          TOOLS/gm-integration-test/initialize-gm-integration-test.sh
#      or, if you want to save the listings that the script produces into gm-test-${TRUE_HOST}-listings.txt file, type:
#          TOOLS/gm-integration-test/initialize-gm-integration-test.sh 2>&1 > gm-test-${TRUE_HOST}-listings.txt
#      or, if you want to compile the model in debug mode, type:
#          TOOLS/gm-integration-test/initialize-gm-integration-test.sh -d
#      or, if you want to prepare new control files, type:
#          TOOLS/gm-integration-test/initialize-gm-integration-test.sh -n
#      or, if you want to use different control directory, type:
#          TOOLS/gm-integration-test/initialize-gm-integration-test.sh -l LOCAL_DIRECTORY_PATH
#      where LOCAL_DIRECTORY_PATH is the path to the control directory of your choice,
#      or, if you want all of the above, type:
#          TOOLS/gm-integration-test/initialize-gm-integration-test.sh -d -n -l LOCAL_DIRECTORY_PATH 2>&1 > gm-test-${TRUE_HOST}-listings.txt
#
#   If the integration test is successful in reproducing the control output, the validation script will print:
#   ***
#   Congratulations: New binary reproduces control test.
#   ***
#   into the listing file named gm-test-${TRUE_HOST}_${latest_commit}/listing/vldtgm*out file
#
#   Optionally, if the integration test is set to prepare new control output, the preparation script will print:
#   ***
#   Congratulations: New control output is prepared.
#   ***
#   into the listing file named gm-test-${TRUE_HOST}_${latest_commit}/listing/vldtgm*out file
#
#   The script relies on:
#       1. Being started from the root of GEM-MACH git repository
#       2. Existence of README_version.md file in the root directory of GEM-MACH git repository
#       3. Correct version of GEM being stored in the mentioned README_version.md file
#       4. Existence of TOOLS/gm-integration-test/compile-gm-for-integration-test.sh script in the GEM-MACH git repository
#   If any of these are missing, or not up to date, the script will fail. The error messages are provided throughout the
#   script to indicate the obvious issues and failure of the script to behave as expected.
#
#   The script will:
#       - Build a directory ${TASK_BASEDIR} to be a working directory for the script and link it to gm-test-${TRUE_HOST}_${latest_commit} in
#         the root directory of the GEM-MACH repository.
#       - Copy current status of GEM-MACH local repository's CHEM and GEM directories into a new directory
#         ${TASK_BASEDIR}/GEM-MACH, where GEM-MACH will be compiled using ${TASK_BASEDIR}/storage_model as a storage model.
#       - Clean GEM-MACH code directory within the integration test work directory (i.e. ${TASK_BASEDIR}/GEM-MACH) from the
#         residuals of any previous compilation in the GEM-MACH git repository.
#       - Build a directory structure necessary for GEM's runmod.sh script to run GEM-MACH:
#         - input, where the script will copy the input files for running predefined configuration of GEM-MACH, and where
#           GEM's runmod.sh script expects the input to be in
#         - bin, where the script will copy GEM scripts necessary for compiling, running and validating GEM-MACH integration test,
#           and where the scripts for running and validating will copy necessary scripts
#         - work, which GEM's runmod.sh script uses as a work directory, and where it produces tmpdir${TRUE_HOST}* directory
#           with the listings from each processor running GEM-MACH
#         - output, where GEM's runmod.sh script saves the output of the GEM-MACH run in the format of Maestro's Runmod
#           task's output
#         - listing, where GEM's compile-gm-for-integration-test.sh, run-gm-integration-test.sh and validate-gm-integration-test.sh scripts save
#           the listings of compiling, running and validating GEM-MACH integration test
#       - Submit compilation script that compiles the code in the current state of GEM-MACH local repository.
#       - Stop and provide the error message when there are issues with any of the described steps.
#       - Save the information about the model versions, repository, computer and location of the input and control output
#         files into gm-test-${TRUE_HOST}-info.txt_${latest_commit} file.
#
#   Optionally, the script can:
#       - Compile the model in debug mode if "-d" option is provided. In that case, the script will link integration test working directory
#         to gm-test-${TRUE_HOST}_${latest_commit}_${dbg} in the root of the GEM-MACH git repository. Note that the validation script will compare
#         the test results with the control-output files saved in gm-output/${IntegrationTest_version}_${TRUE_HOST}_dbg/model subdirectory.
#       - Prepare new control files if "-n" option is provided. In that case, the script will pass "new" value for "cntrl_fl_opt" variable to
#         the compilation script, which will pass it further to the script for running the model. Script for running GEM-MACH will use this
#         information to start the script for preparing new control output in ${TASK_BASEDIR}/new_control/model. Note that this option replaces
#         validation of test results against the control output with the creation of the new output that can be utilized to prepare new control output.
#       - Utilize the control input and output of your choice if "-l" option followed by "LOCAL_DIRECTORY_PATH" is provided. In that case, the script
#         will use the provided path instead of the default one on the common disk space. Note that the script requires that the control directory
#         contains two subdirectories: gm-input for the input files and gm-output/${Test_version}_${TRUE_HOST}/model for control output used for
#         verification.
#
# Author: Verica Savic-Jovcic
# Date:   March 2022
# Update: January 2023
#
# 2023-Apr/May Jack C - update to use git clone from local git repo and to
#   compile using cmake, there's no version check, thus no compatible with RDE
#   (only work with cmake compile version), also some simplificaiton of script.
#
###

scriptstartdate=$(date '+%C%y%m%d%H%M%S')

usage="\n
USAGE: TOOLS/gm-integration-test/initialize-gm-integration-test.sh [-h | <-d> <-n> <-l DIR> 2>&1 > gm-test-${TRUE_HOST}-listings.txt] \n
OPTIONS: \n
  -h      help\n
  -d      compilation in debug mode\n
  -n      prepare new control file\n
  -l DIR  sets the control input and output directory to be 'DIR'\n
'< >' brackets surrounding options indicate that you can choose any combination of these options.\n\n
For description of the integration test, please read the notes at the beginning of the script.\n"

# Inform the world how to get information about the integration test
cmpl_opt=optm
cntrl_fl_opt=vldt
local_control_opt=false
while getopts "hdnl:" opt; do
  case $opt in
    h)
      echo -e ${usage}
      exit 0
      ;;
    d)
      echo -e "Compile in debug mode and compare with debug base case"
      cmpl_opt=dbg
      ;;
    n)
      echo -e "Prepare new control file"
      cntrl_fl_opt=new
      ;;
    l)
      local_control_opt=true
      echo -e "Use local control directory"
      local_control_dir=${OPTARG}
      ;;
    *)
      echo "Invalid option '-${OPTARG}'"
      echo -e ${usage}
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))

# Define integration test version:
Test_version=gm320rc1

# Define runmod cpu topo ('-ptopo 1x1x1' MPIxOMP )
export GMJobPtopo=10x8x1  # default GMJobTopo=10x8x1

# Determine hall of working directory
case ${TRUE_HOST} in
 ppp5 | underhill )
 current_hall=hall5 ;;
 ppp6 | robert )
 current_hall=hall6 ;;
esac

# Get GEM git repo info
[[  ! -f ./MANIFEST ]]  && echo "Error: GEM MANIFEST file not found!" && exit
[[ $(cat MANIFEST | sed -n '/NAME/s/.*: //p') != "gem" ]] && echo "Error: wrong project MANIFEST" && exit
GEM_version=$(cat MANIFEST | sed -n '/VERSION/s/.*: //p')
gmgit_dir=$(git rev-parse --show-toplevel)  ## to do, error out if not found
gmgit_bn=$(basename ${gmgit_dir})
gmgit_br=$(git branch --show-current)
gmgit_rm=$(git config --get remote.origin.url)
gmgit_id=$(git rev-parse --verify --short HEAD)
mach_dir=${gmgit_dir}/src/mach

## Version tag
vertag=${gmgit_id} # option tag ${GMJobPtopo}
[[ "${cmpl_opt}" == "dbg" ]] && vertag=${vertag}_dbg

# Prepare the file for saving the information about the test
gmtestinfo=${gmgit_dir}/gm-test-${TRUE_HOST}-info.txt_${vertag}

# Check and load GEM environment
gembndl=/fs/ssm/eccc/mrd/rpn/models/gem/bundles/${GEM_version}.bndl
[[ ! -f $gembndl ]] && echo -d " ERROR: GEM environment file not available: $gembndl \n" && exit 1
source r.load.dot /fs/ssm/eccc/mrd/rpn/models/gem/bundles/${GEM_version}

# Set test work directory structure
export TASK_BASEDIR=/space/${current_hall}/sitestore/eccc/aq/r1/${USER}/maestro/${TRUE_HOST}/${Test_version}/gm-test_${gmgit_bn}_${vertag}
[[ -d ${TASK_BASEDIR} ]] && rm -rf ${TASK_BASEDIR} && mkdir -p ${TASK_BASEDIR}
ln -sf ${TASK_BASEDIR} $gmgit_dir
export TASK_BIN=${TASK_BASEDIR}/bin
export TASK_INPUT=${TASK_BASEDIR}/input
export TASK_WORK=${TASK_BASEDIR}/work
export TASK_OUTPUT=${TASK_BASEDIR}/output
mkdir -p ${TASK_BIN} ${TASK_INPUT} ${TASK_WORK} ${TASKOUTPUT} ${TASK_BASEDIR}/listing

# Tell the world about the repository, compiler and computer
cat << EOF > ${gmtestinfo}

GEM local git location: ${gmgit_dir}
GEM version: ${GEM_version}
GEM remote repo: ${gmgit_rm}
GEM branch: ${gmgit_br}
GEM commit hash: ${gmgit_id}

Integration test version: ${Test_version}

GEM environment: ${gembndl}
Current machine: ${TRUE_HOST}
Operating system and architecture: ${ORDENV_PLAT}

Integration test location: ${TASK_BASEDIR}
EOF

# git clone local and checkout current commit to working directory
gemmach_dir=${TASK_BASEDIR}/build
[[ -d ${gemmach_dir} ]]  && rm -rf ${gemmach_dir}
git clone -l --single-branch $gmgit_dir ${gemmach_dir}
cd ${gemmach_dir}; git reset --hard ${gmgit_id}

# Set input directory and copy input files
control_dir=/space/${current_hall}/sitestore/eccc/aq/r1/sarq000/gmtest/${Test_version}
[[ "${local_control_opt}" == "true" ]] && control_dir=${local_control_dir}
control_input_dir=${control_dir}/gm-input
[[ ! -d $control_input_dir ]] && echo "Error: control_input_dir not exist: $control_input_dir" && exit 1
cp -r ${control_input_dir}/* ${TASK_INPUT}/

## override config files e.g:
##  'physics_input_table->gm_phy_intable'; 'model_settings.nml->gem_settings.nml'
##  'output_settings->outcfg.out'
for file in ${mach_dir}/tools/gemmach_cfg/*; do
   [ -e "$file" ] || continue
   echo -e "override file: ${file} \n"
   cp -fv ${file} ${TASK_INPUT}
   filename=$(basename $file)
   [[ $filename == 'gem_settings.nml' ]] && ln -sfv ${TASK_INPUT}/${filename} ${TASK_INPUT}/cfg_0000/model_settings.nml
   [[ $filename == 'outcfg.out' ]] && ln -sfv ${TASK_INPUT}/${filename} ${TASK_INPUT}/cfg_0000/output_settings
   [[ $filename == 'gm_phy_intable' ]] && ln -sfv ${TASK_INPUT}/${filename} ${TASK_INPUT}/cfg_0000/physics_input_table
done

## if dbg, turn on debug_trace- warning: produce large runmod listing
if [[ "${cmpl_opt}" == "dbg" ]] ; then
 [[ -L ${TASK_INPUT}/cfg_0000/model_settings.nml ]] && \
   /bin/cp --remove-destination `readlink ${TASK_INPUT}/cfg_0000/model_settings.nml` ${TASK_INPUT}/cfg_0000/model_settings.nml
 echo -e "\n DBG: turn on debug_trace_L and chm_debug_trace_l in the namelist \n"
 ln -s $(which setnml) ${TASK_BIN}/setnml
 ${TASK_BIN}/setnml -f ${TASK_INPUT}/cfg_0000/model_settings.nml physics_cfgs/debug_trace_L=\.true\.
 ${TASK_BIN}/setnml -f ${TASK_INPUT}/cfg_0000/model_settings.nml chemistry_cfgs/chm_debug_trace_l=\.true\.
fi

# Copy the integration-test related scripts into the task's bin directory
cp ${mach_dir}/tools/gm-integration-test/initialize-gm-integration-test.sh ${TASK_BIN}
cp ${mach_dir}/tools/gm-integration-test/compile-gm-for-integration-test.sh ${TASK_BIN}
cp ${mach_dir}/tools/gm-integration-test/run-gm-integration-test.sh ${TASK_BIN}
cp ${mach_dir}/tools/gm-integration-test/validate-gm-integration-test.sh ${TASK_BIN}
if [[ "${cntrl_fl_opt}" == "new" ]] ; then cp ${mach_dir}/tools/gm-integration-test/prepare-new-gm-integration-test.sh ${TASK_BIN} ; fi

# Set the task resources
export GMJobMach=${TRUE_HOST}
export GMJobMemory=2G
export GMJobQueue=development
export GMJobTopo=$GMJobPtopo
export CmplJobTime=600
export CmplJobNcpu=80
export CmplJobName=cmplgm

# Submit the GEM-MACH compilation script
echo -e "\n == Submit the GEM-MACH integration-test compilation script at $(date) == \n"
ord_soumet ${TASK_BIN}/compile-gm-for-integration-test.sh \
           -args "${GEM_version} ${Test_version} ${control_dir} ${TASK_BASEDIR} ${GMJobMach} ${GMJobMemory} ${GMJobQueue} ${GMJobTopo} ${gmtestinfo} ${cmpl_opt} ${cntrl_fl_opt}" \
           -mach ${GMJobMach} -cpus ${CmplJobNcpu} -cm ${GMJobMemory} -t ${CmplJobTime} \
           -mpi 1 -queue ${GMJobQueue} -jn ${CmplJobName} -listing ${TASK_BASEDIR}/listing

# Inform the world about the integration test
cat << EOF | tee -a ${gmtestinfo}

== Location of listings: ==

GEM-MACH compilation submitted to ${GMJobMach} at $(date).
Listings for the compilation job are in ${TASK_BASEDIR}/listing/cmplgm*out file.

When compilation job finishes, automatically initiated GEM-MACH run by submitting the job named "rungm" to ${TRUE_HOST}.
Listings for the GEM-MACH run are in ${TASK_BASEDIR}/listing/rungm*out file.

EOF

if [[ "${cntrl_fl_opt}" != "new" ]] ; then
 cat << EOF | tee -a ${gmtestinfo}
When GEM-MACH run job finishes automoatically initiated validation of GEM-MACH run results by submitting the job named "vldtgm" to ${TRUE_HOST}.
Listings for the validation job are in ${TASK_BASEDIR}/listing/vldtgm*out file.

When validation job finishes, look for word "Congratulations" in the listings of validation job.

EOF
else
 cat << EOF | tee -a ${gmtestinfo}
When GEM-MACH run job finishes automoatically initiated preparation of new control output submitting the job named "prepgm" to ${TRUE_HOST}.
Listings for the new contol-output preparation job are in ${TASK_BASEDIR}/listing/prepgm*out file.

When preparation job finishes, look for "Congratulations" and/or "Location of new control output:" phrase in the listings of preparation job.

EOF
fi
