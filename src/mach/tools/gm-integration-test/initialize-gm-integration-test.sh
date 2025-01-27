#!/bin/bash
#set -e

###
#
#   Script for initiating GEM-MACH integration test from the root directory of the MACH repository.
#
#   USAGE: While in the root directory of MACH repository, or src/mach directory of GEM super repository, type:
#          tools/gm-integration-test/initialize-gm-integration-test.sh
#      or, if you want to save the listings that the script produces into gm-test-${TRUE_HOST}-listings.txt file, type:
#          tools/gm-integration-test/initialize-gm-integration-test.sh 2>&1 > gm-test-${TRUE_HOST}-listings.txt
#      or, if you want to compile the model in debug mode, type:
#          tools/gm-integration-test/initialize-gm-integration-test.sh -d
#      or, if you want to run with debug tracing only in chemistry turned on, type:
#          tools/gm-integration-test/initialize-gm-integration-test.sh -c
#      or, if you want to run with debug tracing only in physics turned on, type:
#          tools/gm-integration-test/initialize-gm-integration-test.sh -p
#      or, if you want to prepare new control files, type:
#          tools/gm-integration-test/initialize-gm-integration-test.sh -n
#      or, if you want to use different control directory, type:
#          tools/gm-integration-test/initialize-gm-integration-test.sh -l LOCAL_DIRECTORY_PATH
#      where LOCAL_DIRECTORY_PATH is the path to the control directory of your choice,
#      or, if you want to run the test from the local MACH repository instead of default local GEM super repository:
#          tools/gm-integration-test/initialize-gm-integration-test.sh -m
#      or, if you want to run the test from local MACH repository with different history then the one in remote GEM super repository:
#          tools/gm-integration-test/initialize-gm-integration-test.sh -m -u
#      or, if you want to run the test with GEM super repository and branch different from MIG/gem and GEMs_version defined in MANIFEST:
#          tools/gm-integration-test/initialize-gm-integration-test.sh -g GEM_remote,GEM_version
#      where GEM_remote is location of GEM super repository of your choice and GEM_version is your choice of a branch, or a tag in your 
#      chosen super repository; Note that GEM_remote and GEM_version have to be separated with comma (,).
#      or, if you want to compile and create an ssm package:
#          tools/gm-integration-test/initialize-gm-integration-test.sh -s
#      or, if you want all of the above, type:
#          tools/gm-integration-test/initialize-gm-integration-test.sh -d -c -p -n -l LOCAL_DIRECTORY_PATH -m -u -s 2>&1 > gm-test-${TRUE_HOST}-listings.txt
#
#   If the integration test is successful in reproducing the control output, the validation script will print:
#   ***
#   New binary reproduces reference output.
#   ***
#   into the listing files named gm-test-${TRUE_HOST}_${latest_commit}/listings/vldtgm*out and gm-test-${TRUE_HOST}-info.txt_${latest_commit}
#
#   Optionally, if the integration test is set to prepare new control output, the script for running GEM-MACH and combining the output will print:
#   ***
#   Location of new control output: <path to directory with the model output>
#   ***
#   into the listing files named gm-test-${TRUE_HOST}_${latest_commit}/listings/rungm*out and gm-test-${TRUE_HOST}-info.txt_${latest_commit}
#
#   The script relies on:
#       1. Being started from the root of MACH git repository, or from src/mach directory of GEM git super repository
#       2. Existence of MANIFEST file in the root directory of MACH repository (i.e. src/mach directory of GEM super repository)
#       3. Existence of tools/gm-integration-test/compile-gm-for-integration-test.sh script in the MACH repo (i.e. src/mach directory of GEM super repo)
#   If any of these are missing, or not up to date, the script will fail. The error messages are provided throughout the
#   script to indicate the obvious issues and failure of the script to behave as expected.
#
#   The script will:
#       - Create a directory ${TASK_BASEDIR} to be a work directory for the script and link it to gm-test-${TRUE_HOST}_${latest_commit} in
#         the root directory of MACH repository (i.e. src/mach directory of GEM super repo).
#       - If initiated from MACH repository (with -m option):
#          * Read MACH version from the MANIFEST and use it to git clone correct version of GEM to the GEM-MACH subdirectory of work directory
#          * git subtree pull current branch of the local MACH repository into the GEM super repo within work directory
#       - If initiated from GEM super repository:
#          * git clone local GEM super repo that contains MACH code into the GEM-MACH subdirectory of work directory
#       - Build a directory structure necessary for GEM's runmod.sh script to run GEM-MACH:
#         - input, where the script will copy the input files for running predefined configuration of GEM-MACH, and where
#           GEM's runmod.sh script expects the input to be in
#         - bin, where the scripts of the test will copy GEM scripts necessary for compiling, running and validating GEM-MACH integration test,
#         - work, which GEM's runmod.sh script uses as a work directory, and where it produces tmpdir${TRUE_HOST}* directory
#           with the listings from each processor running GEM-MACH
#         - output, where GEM's runmod.sh script saves the output of the GEM-MACH run in the format of Maestro's Runmod
#           task's output
#         - listings, where GEM-MACH's compile-gm-for-integration-test.sh, run-gm-integration-test.sh and validate-gm-integration-test.sh scripts save
#           the listings of compiling, running and validating GEM-MACH integration test
#       - Submit compilation script that compiles the code prepared in the GEM-MACH directory of work directory.
#       - Stop and provide the error message when there are issues with any of the described steps.
#       - Save the information about the model versions, repository, computer and location of the input and control output
#         files into gm-test-${TRUE_HOST}-info.txt_${latest_commit} file.
#
#   Optionally, the script can set integration test to:
#       - Compile the model in debug mode if "-d" option is provided. In that case, the script will link integration test working directory
#         to gm-test-${TRUE_HOST}_${latest_commit}_${dbg} in the root of the GEM-MACH git repository. Note that the validation script will compare
#         the test results with the control-output files saved in gm-output/${IntegrationTest_version}_${TRUE_HOST}_dbg/model subdirectory.
#       - Run model with debug tracing turned on in chemistry if "-c" option is provided. In that case, namelist will be updated to set 
#         "chm_debug_trace_l" to true.
#       - Run model with debug tracing turned on in physics (including parts of chemistry code) if "-p" option is provided. In that case, 
#         namelist will be updated to set "debug_trace_l" to true.
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
#         Apr/May 2023, Jack C - update to use git clone from local git repo and to
#          compile using cmake, there's no version check, thus no compatible with RDE
#          (only work with cmake compile version), also some simplification of script.
#         May 2023, Jack C - run outputs is combined in `run-gm-integration-test` and
#          submit `validate-gm-integration-test` when '$cntrl_fl_opt!=new'
#          Header notes of have not been revised.
#         August 2023, Verica S-J - update to run the test from MACH repository instead
#          of from GEM repository; update documentation; rename build directory in work
#          directory to GEM-MACH directory
#         October/November 2023, Verica S-J - update to allow for running the test either
#          from MACH repo or GEM super repo and to load environments from GEM super repo
#          instead of released package
#         June 2024, Verica S-J - update to assure that test works in debug mode fully and
#          to allow for tracing regardless of compilation mode
###

scriptstartdate=$(date '+%C%y%m%d%H%M%S')

usage="\n
USAGE:\n\n
cd <mach>\n
./tools/gm-integration-test/initialize-gm-integration-test.sh [-h -d -c -p -n -l DIR -m -u -g GEM_remote,GEM_version -s]\n\n
<mach> is the location of src/mach subdirectory of local GEM super repository, or when -m option is used of local MACH repository.\n\n
OPTIONS:\n
  -h      -> help\n
  -d      -> compilation in debug mode\n
  -c      -> run the model with debug tracing turned on only in chemistry\n
  -p      -> run the model with debug tracing turned on in physics\n
  -n      -> create new control output, i.e. skip fstcomp with control output in validate-gm-integration-test.sh\n 
  -l DIR  -> sets the control input and output directory to be 'DIR'\n
  -m      -> run the test from local MACH repository\n
  -u      -> different MACH history in MACH repository and MIG super repository\n
  -g      -> provide path and branch, or tag of GEM super repository\n
  -s      -> compile and create an ssm package \n\n
You can choose any combination of these options.\n\n
NOTES:\n
- Integration test is done on the code that is committed to the local git repository last.\n
- For description of the integration test, please read the notes at the beginning of the script.\n"

# Inform the world how to get information about the integration test
cmpl_opt=optm
chm_trcng_opt=off
phy_trcng_opt=off
cntrl_fl_opt=vldt
local_control_opt=false
mach_repo=false
common_mach_history=true
gem_repo_opt=false
gem_branch_opt=false
while getopts "hdcpnl:mug:s" opt; do
  case $opt in
    h)
      echo -e ${usage}
      exit 0
      ;;
    d)
      echo -e "Compile in debug mode and run with tracing turned off.\n"
      cmpl_opt=dbg
      ;;
    c)
      echo -e "Run with tracing turned on only in chemistry code. \n  Warning: This produces very long listings.\n"
      chm_trcng_opt=on
      ;;
    p)
      echo -e "Run with tracing turned on in physics code. \n  Warning: This produces very long listings.\n"
      phy_trcng_opt=on
      ;;
    n)
      echo -e "Prepare new control files.\n"
      cntrl_fl_opt=new
      ;;
    l)
      local_control_opt=true
      echo -e "Use local control directory.\n"
      local_control_dir=${OPTARG}
      ;;
    m)
      echo -e "Run the test from local MACH repository using remote GEM super repository.\n"
      mach_repo=true
      ;;
    u)
      echo -e "Different MACH hisotries in local MACH repository and remote GEM super repository.\n"
      common_mach_history=false
      ;;
    g)
      echo -e "Provide GEM super repository and branch, or tag to compile GEM-MACH with.\n"
      gem_repo_opt=true
      gem_branch_opt=true
      GEM_remote=$(echo ${OPTARG} | cut -d "," -f1)
      GEM_version=$(echo ${OPTARG} | cut -d "," -f2)
      ;;
    s)
      echo -e "Compile and create an ssm package"
      cmpl_opt=pkg
      ;;
    *)
      echo "Invalid option '-${OPTARG}'"
      echo -e ${usage}
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))

# Show reminder
echo -e "\n
REMINDER:\n
Test is done on the last-committed code in your local git repository.\n"


# Define remote GEM super repository
if [[ "${gem_repo_opt}" == "false" ]] ; then
 GEM_remote='git@gitlab.science.gc.ca:MIG/gem.git'
fi

# Get information about local git repository
[[  ! -f ./MANIFEST ]]  && echo -e "Error: GEM MANIFEST file not found.\n" && exit 1
[[ $(cat MANIFEST | sed -n '/NAME/s/.*: //p') != "mach" ]] && echo -e "Error: wrong project MANIFEST.\n" && exit 1
repo_dir=$(git rev-parse --show-toplevel)  ## to do, error out if not found
repo_bn=$(basename ${repo_dir})
repo_br=$(git branch --show-current)
repo_id=$(git rev-parse --verify --short HEAD)

# Determine hall of working directory
if [[ -n "$(echo ${repo_dir} | grep /homeu2/ )" ]] ; then
 current_hall=hall5
else
 current_site=$(echo ${repo_dir} | cut -d "/" -f3)
 current_hall=$(echo "${current_site//site/hall}")
fi

if [[ "${mach_repo}" != "true" ]] ; then
   gmgit_dir=${repo_dir}/src/mach
   [[ ! -d  ${gmgit_dir} ]] && echo -e "Error: When running test from MACH repository, have to use -m option.\n" && exit 1 
   repo_nm="GEM-MACH"
elif [[ "${mach_repo}" == "true" ]] ; then
   gmgit_dir=${repo_dir}
   repo_nm="MACH"
fi

# Get GEM-MACH, GEM and integration test versions
MACH_version=$(cat ${gmgit_dir}/MANIFEST | sed -n '/^VERSION/s/.*: //p')
if [[ "${gem_branch_opt}" == "false" ]] ; then
 GEM_version=$(cat ${gmgit_dir}/MANIFEST | sed -n '/^GEM_VERSION/s/.*: //p')
fi
Test_version=gm$(sed 's/[.-]//g' <<< ${MACH_version})

## Version tag (i.e. commit hash)
vertag=${repo_id} # option tag ${GMJobPtopo}
[[ "${cmpl_opt}" == "dbg" ]] && vertag=${vertag}_dbg

# Name the integration-test listings file
gmtestinfo=${gmgit_dir}/gm-test-${TRUE_HOST}-info.txt_${vertag}

# Set work directory structure
export TASK_BASEDIR=$(true_path -n ${HOME}/data_maestro/${TRUE_HOST})/maestro/${TRUE_HOST}/${Test_version}/gm-test_${repo_bn}_${vertag}
[[ -d ${TASK_BASEDIR} ]] && rm -rf ${TASK_BASEDIR} && mkdir -p ${TASK_BASEDIR}
ln -sf ${TASK_BASEDIR} ${gmgit_dir}/gm-test-${TRUE_HOST}_${vertag}
export TASK_BIN=${TASK_BASEDIR}/bin
export TASK_INPUT=${TASK_BASEDIR}/input
export TASK_WORK=${TASK_BASEDIR}/work
export TASK_OUTPUT=${TASK_BASEDIR}/output
mkdir -p ${TASK_BIN} ${TASK_INPUT} ${TASK_WORK} ${TASKOUTPUT} ${TASK_BASEDIR}/listings

# Tell the world about the integration test, model repositories, compiler and computer
cat << EOF > ${gmtestinfo}

Integration-test version: ${Test_version}
Integration-test work directory: ${TASK_BASEDIR}

MACH version: ${MACH_version}
${repo_nm} local git repository: ${gmgit_dir}
${repo_nm} branch: ${repo_br}
${repo_nm} commit hash: ${repo_id}

GEM version: ${GEM_version}
GEM remote git repository: ${GEM_remote}

Current machine: ${TRUE_HOST}
Operating system and architecture: ${ORDENV_PLAT}

EOF

# Using git commands, copy GEM-MACH code to the directory where the model will be compiled
gemmach_dir=${TASK_BASEDIR}/GEM-MACH
[[ -d ${gemmach_dir} ]]  && rm -rf ${gemmach_dir}

if [[ "${mach_repo}" != "true" ]] ; then
   # git clone local GEM super repo and checkout current commit to working directory
   git clone -l --single-branch --no-hardlinks ${repo_dir} ${gemmach_dir}
   cd ${gemmach_dir}; git reset --hard ${repo_id}
elif [[ "${mach_repo}" == "true" ]] ; then
   # git clone remote GEM super repo to working directory, reset it to the tagged GEM version, and
   # git subtree pull local MACH repo into the working-directory GEM super repo
   mach_dir=${gemmach_dir}/src/mach
   git clone -b ${GEM_version} ${GEM_remote} ${gemmach_dir}
   cd ${gemmach_dir}
   git switch -c ${Test_version}
   export GIT_MERGE_AUTOEDIT=no
   if [[ "${common_mach_history}" == "true" ]] ; then
      git subtree pull --squash -P src/mach ${gmgit_dir} ${repo_br}
   else
      git rm -rf src/mach
      git commit -m "Remove src/mach subtree"
      git subtree add --squash -P src/mach ${gmgit_dir} ${repo_br}
   fi
fi

# Load environments needed for compiling and running GEM (instead of loading GEM environment)
source ${gemmach_dir}/.eccc_setup_intel

# Set location of control directory and copy input files to work directory
control_dir=/space/${current_hall}/sitestore/eccc/aq/r1/sarq000/gmtest/${Test_version}
[[ "${local_control_opt}" == "true" ]] && control_dir=${local_control_dir}
control_input_dir=${control_dir}/gm-input
[[ ! -d $control_input_dir ]] && echo "Error: control_input_dir not exist: $control_input_dir" && exit 1
cp -r ${control_input_dir}/* ${TASK_INPUT}/
echo -e "Integration-test input file: ${control_input_dir}" | tee -a ${gmtestinfo}

## If local config files are provided, override config files e.g:
##  'physics_input_table->gm_phy_intable'; 'model_settings.nml->gem_settings.nml'
##  'output_settings->outcfg.out'
for file in ${gmgit_dir}/tools/gemmach_cfg/*; do
   [ -e "$file" ] || continue
   echo -e "override file: ${file} \n"
   cp -fv ${file} ${TASK_INPUT}
   filename=$(basename $file)
   [[ $filename == 'gem_settings.nml' ]] && ln -sfv ${TASK_INPUT}/${filename} ${TASK_INPUT}/cfg_0000/model_settings.nml
   [[ $filename == 'outcfg.out' ]] && ln -sfv ${TASK_INPUT}/${filename} ${TASK_INPUT}/cfg_0000/output_settings
   [[ $filename == 'gm_phy_intable' ]] && ln -sfv ${TASK_INPUT}/${filename} ${TASK_INPUT}/cfg_0000/physics_input_table
done

## For either debug compilation or running with tracing, namelist has to be modified locally
if [[ "${cmpl_opt}" == "dbg" ]] || [[ "${chm_trcng_opt}" == "on" ]] || [[ "${phy_trcng_opt}" == "on" ]] ; then
 [[ -L ${TASK_INPUT}/cfg_0000/model_settings.nml ]] && \
   /bin/cp --remove-destination `readlink ${TASK_INPUT}/cfg_0000/model_settings.nml` ${TASK_INPUT}/cfg_0000/model_settings.nml
 ln -s $(which setnml) ${TASK_BIN}/setnml
fi

if [[ "${cmpl_opt}" == "dbg" ]] ; then
 ${TASK_BIN}/setnml -f ${TASK_INPUT}/cfg_0000/model_settings.nml step/Fcst_end_S='9h'
fi

if [[ "${chm_trcng_opt}" == "on" ]] ; then
 ${TASK_BIN}/setnml -f ${TASK_INPUT}/cfg_0000/model_settings.nml chemistry_cfgs/chm_debug_trace_l=\.true\.
fi

if [[ "${phy_trcng_opt}" == "on" ]] ; then
 ${TASK_BIN}/setnml -f ${TASK_INPUT}/cfg_0000/model_settings.nml physics_cfgs/debug_trace_L=\.true\.
fi

# Copy the integration-test related scripts into the task's bin directory
cp ${gmgit_dir}/tools/gm-integration-test/compile-gm-for-integration-test.sh ${TASK_BIN}
cp ${gmgit_dir}/tools/gm-integration-test/run-gm-integration-test.sh ${TASK_BIN}
cp ${gmgit_dir}/tools/gm-integration-test/validate-gm-integration-test.sh ${TASK_BIN}

# Set the resources for running GEM-MACH
export GMJobMach=${TRUE_HOST}
export GMJobMemory=2G
export GMJobQueue=development
export GMJobTopo=10x8x1
export CmplJobTime=600
export CmplJobNcpu=80
export CmplJobName=cmplgm

# Submit the GEM-MACH compilation script
echo -e "\n == Submit the GEM-MACH integration-test compilation script at $(date) == \n"
ord_soumet ${TASK_BIN}/compile-gm-for-integration-test.sh \
           -args "${Test_version} ${control_dir} ${TASK_BASEDIR} ${GMJobMach} ${GMJobMemory} ${GMJobQueue} ${GMJobTopo} ${gmtestinfo} ${cmpl_opt} ${cntrl_fl_opt}" \
           -mach ${GMJobMach} -cpus ${CmplJobNcpu} -cm ${GMJobMemory} -t ${CmplJobTime} \
           -mpi 1 -queue ${GMJobQueue} -jn ${CmplJobName} -listing ${TASK_BASEDIR}/listings

# Inform the world about the integration test
cat << EOF | tee -a ${gmtestinfo}

== Location of listings: ==

GEM-MACH compilation submitted to ${GMJobMach} at $(date).
Listings for the compilation job are in ${TASK_BASEDIR}/listings/cmplgm*out file.

At the end of compilation job, it initiates GEM-MACH run by submitting the job named "rungm" to ${TRUE_HOST}.
Listings for the GEM-MACH run are in ${TASK_BASEDIR}/listings/rungm*out file.

EOF

if [[ "${cntrl_fl_opt}" != "new" ]] ; then
 cat << EOF | tee -a ${gmtestinfo}
At the end of GEM-MACH run job, it initiates validation of GEM-MACH run results by submitting the job named "vldtgm" to ${TRUE_HOST}.
Listings for the validation job are in ${TASK_BASEDIR}/listings/vldtgm*out file.

When validation job finishes, look for "New binary reproduces reference output" phrase in the listings of vldtgm job. It signifies a successsful execution of the integration test.

EOF
else
 cat << EOF | tee -a ${gmtestinfo}
The GEM-MACH run job also combines the model output to expected YYYYMMDDHH_hhh files, which can be used as a new control output.

When preparation job finishes, look for "Location of new control output" phrase in the listings of rungm job. It signifies a successsful execution of the integration test.

EOF
fi
