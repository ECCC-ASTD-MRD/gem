#!/bin/bash

###
#
#  Script for initiating GEM-MACH regression test 
#
#  USAGE: 
#
#    tools/regression-test/initialize-regression-test.sh [-h -d -c -p -n -l DIR -m -u -g GEM_remote,GEM_version -s -o -i]
#
#    Run this command either from the root of the MACH subtree
#    repository, or from the `src/mach` directory in the GEM super 
#    repository:
#      
#    Description of options is provided when script is called with 
#    "-h" option.
#
#  The script relies on:
#      1. Being started from the root of MACH git repo, or from 
#         src/mach directory of GEM git super repo
#      2. Existence of MANIFEST file and tools/regression-test/*sh
#         scripts in the root directory of MACH repo (i.e. 
#         src/mach directory of GEM super repo)
#      3. Existence of ${HOME}/data_maestro directory with links to
#         sitestore disk space
#  If any of these are missing, or not up to date, the script will fail.
#  Error messages are provided throughout the script to indicate any
#  obvious issues and failure of the script to behave as expected.
#
#  Script:
#      - Creates a directory ${TASK_BASEDIR} to be a work directory for
#        the regression test and links it to 
#        gm-test-${TRUE_HOST}_${latest_commit} in the root directory 
#        of MACH repo (i.e. src/mach directory of GEM super repo).
#      - Builds GEM-MACH source code in ${TASK_BASEDIR} either by 
#        clonning local GEM super repo, or by cloning a GEM super repo
#        and replacing its mach directory with the local MACH repo.
#      - Checks if the header modules in MACH library have been updated.
#      - Builds a directory structure necessary for GEM's runmod.sh 
#        script to run GEM-MACH (bin, input, work, output, listings).
#      - Links and copies relevant inputs, configuration files and
#        scripts to respective directories. 
#      - Saves the information about the model versions, repo, computer
#        and location of the input and control output files into 
#        gm-test-${TRUE_HOST}-info.txt_${latest_commit} file.
#      - Saves locations of directories and chosen options into the
#        ${TASK_BASEDIR}/*_env_arg configuration files to be used by 
#        other scripts of regression test.
#      - Submits compilation script.
#
# Author: Verica Savic-Jovcic
# Date:   March 2022
#
#         Jack Chen and Balbir Pabla contributed to updating this script.
#
###

scriptstartdate=$(date '+%C%y%m%d%H%M%S')

# Inform the world how to get information about the regression test
usage="\n
USAGE:\n\n
cd <mach>\n
./tools/regression-test/initialize-regression-test.sh [-h -d -c -p -n -l DIR -m -u -g GEM_remote,GEM_version -s -o -i]\n\n
<mach> is the location of src/mach subdirectory of local GEM super repository, or when -m option is used of local MACH repository.\n\n
OPTIONS:\n
  -h      -> help\n
  -d      -> configures compilation script to compile GEM-MACH in debug mode\n
  -c      -> configures GEM-MACH to run with debug tracing turned on only in chemistry\n
  -p      -> configures GEM-MACH to run with debug tracing turned on in physics\n
  -n      -> configures test to create new control output and skip validation of results\n
  -l DIR  -> sets the control input and output directory to be 'DIR'\n
  -m      -> starts the test from local MACH repository, i.e. builds GEM-MACH source code incorporating local MACH repository\n
  -u      -> builds GEM-MACH source code incorporating MACH repository with history divergent from one in adopted GEM super repository\n
  -g GEM_remote,GEM_version  -> builds GEM-MACH source code adopting GEM super repository located at GEM_remote, using GEM_version branch, or tag \n
  -s      -> configures compilation script to additionally create an ssm package\n
  -o      -> configures regression-test scripts to disable automatic submission of subsequent scripts\n
  -i      -> configures regression-test scripts to be utilized in GitLab continuous integration \n\n
You can choose any combination of these options.\n\n
NOTES:\n
- Regression test is done on the code that is committed to the local git repository last.\n
- For description of the regression test, please read readme/README_regression.md.\n\n"

# Default options
cmpl_opt=optm
chm_trcng_opt=off
phy_trcng_opt=off
cntrl_fl_opt=vldt
mach_repo=false
common_mach_history=true
sequence_opt=true
ci_opt=false

# Chosen options
while getopts "hdcpnl:mug:soi" opt; do
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
   control_dir=${OPTARG}
   echo -e "Use ${control_dir} for control directory.\n"
   ;;
  m)
   echo -e "Run the test from local MACH repository using remote GEM super repository.\n"
   mach_repo=true
   ;;
  u)
   echo -e "Different MACH histories in local MACH repository and remote GEM super repository.\n"
   common_mach_history=false
   ;;
  g)
   GEM_remote=$(echo ${OPTARG} | cut -d "," -f1)
   GEM_version=$(echo ${OPTARG} | cut -d "," -f2)
   echo -e "Compile GEM-MACH with GEM from super repository ${GEM_remote} using branch/tag ${GEM_version}.\n"
   ;;
  s)
   echo -e "Compile and create an ssm package"
   cmpl_opt=pkg
   ;;
  o)
   echo -e "Run only initialization script"
   sequench_opt=false
   ;;
  i)
   echo -e "Applying Continuous Integration"
   echo -e "Combine work of -m, -u and -o options"
   ci_opt=true
   mach_repo=true
   common_mach_history=false
   sequence_opt=false
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
GEM_remote=${GEM_remote:-'git@gitlab.science.gc.ca:MIG/gem.git'}

# Get information about local git repository
[[  ! -f ./MANIFEST ]]  && echo -e "ERROR: GEM MANIFEST file not found.\n" && exit 1
[[ $(cat MANIFEST | sed -n '/NAME/s/.*: //p') != "mach" ]] && echo -e "ERROR: wrong project MANIFEST.\n" && exit 1
repo_dir=$(git rev-parse --show-toplevel)  ## to do, error out if not found
repo_bn=$(basename ${repo_dir})
repo_id=$(git rev-parse --verify --short=8 HEAD)
if [[ "${ci_opt}" == "true" ]] ; then
 git switch -c test_${repo_id}
fi
repo_br=$(git branch --show-current)

# Check if the choice of options works from the local git repository
if [[ "${mach_repo}" != "true" ]] ; then
 gmgit_dir=${repo_dir}/src/mach
 [[ ! -d  ${gmgit_dir} ]] && echo -e "ERROR: When running test from MACH repository, have to use -m option.\n" && exit 1 
 repo_nm="GEM-MACH"
elif [[ "${mach_repo}" == "true" ]] ; then
 gmgit_dir=${repo_dir}
 repo_nm="MACH"
fi

# Get GEM-MACH, GEM and regression test versions
MACH_version=$(cat ${gmgit_dir}/MANIFEST | sed -n '/^VERSION/s/.*: //p')
GEM_version=${GEM_version:-$(cat ${gmgit_dir}/MANIFEST | sed -n '/^GEM_VERSION/s/.*: //p')}
Test_version=gm$(sed 's/[.-]//g' <<< ${MACH_version})

## Version tag (i.e. commit hash)
vertag=${repo_id} # option tag ${GMJobPtopo}
[[ "${cmpl_opt}" == "dbg" ]] && vertag=${vertag}_dbg

# Set base directory for the test (location and structure)
TASK_BASEDIR=$(true_path -n ${HOME}/data_maestro/${TRUE_HOST})/maestro/gm-test_${repo_bn}_${vertag}
[[ -d ${TASK_BASEDIR} ]] && rm -rf ${TASK_BASEDIR} && mkdir -p ${TASK_BASEDIR}
TASK_BIN=${TASK_BASEDIR}/bin
TASK_INPUT=${TASK_BASEDIR}/input
TASK_WORK=${TASK_BASEDIR}/work
TASK_OUTPUT=${TASK_BASEDIR}/output
TASK_LIST=${TASK_BASEDIR}/listings
mkdir -p ${TASK_BIN} ${TASK_INPUT} ${TASK_WORK} ${TASKOUTPUT} ${TASK_LIST}

# Make it easier for the world to access info about the test
if [[ ${ci_opt} == false ]] ; then
 # Link the test directory to the local git repository and name the file where to store information about the test
 ln -sf ${TASK_BASEDIR} ${gmgit_dir}/gm-test-${TRUE_HOST}_${vertag}
 gmtestinfo=${gmgit_dir}/gm-test-${TRUE_HOST}-info.txt_${vertag}
else
 # When the script is used by CI, runner controls local git repository and the file with the test info has to be in the test directory
 gmtestinfo=${TASK_BASEDIR}/gm-test-${TRUE_HOST}-info.txt_${vertag}
fi

# Tell the world about the regression test, model repositories, compiler and computer
cat << EOF > ${gmtestinfo}

Regression-test version: ${Test_version}
Regression-test directory: ${TASK_BASEDIR}

MACH version: ${MACH_version}
${repo_nm} local git repository: ${gmgit_dir}
${repo_nm} branch: ${repo_br}
${repo_nm} commit hash: ${repo_id}

GEM version: ${GEM_version}
GEM remote git repository: ${GEM_remote}

Current machine: ${TRUE_HOST}
Operating system and architecture: ${ORDENV_PLAT}

EOF

# Prepare directory for building GEM-MACH
gemmach_dir=${TASK_BASEDIR}/GEM-MACH
[[ -d ${gemmach_dir} ]]  && rm -rf ${gemmach_dir}

# Using git commands, build GEM-MACH source code in the new test directory
if [[ "${mach_repo}" != "true" ]] ; then
 # git clone local GEM super repo to the test directory and checkout current commit
 git clone -l --single-branch --no-hardlinks ${repo_dir} ${gemmach_dir}
 cd ${gemmach_dir}; git reset --hard ${repo_id}
elif [[ "${mach_repo}" == "true" ]] ; then
 # git clone remote GEM super repo to the test directory, reset it to the tagged GEM version, and
 # git subtree pull local MACH repo into the GEM super repo in the test directory
 git clone -b ${GEM_version} ${GEM_remote} ${gemmach_dir}
 cd ${gemmach_dir}
 git switch -c ${Test_version}
 sed -i "s|MACH_VERSION : .*|MACH_VERSION : ${MACH_version}|g" MANIFEST
 sed -i "s|mach = .*|mach = ${MACH_version}|g" MANIFEST 
 git add MANIFEST
 git commit -m "Update mach version in MANIFEST"
 mach_subtree_path="src/mach"
 export GIT_MERGE_AUTOEDIT=no
 if [[ "${common_mach_history}" == "true" ]] ; then
  git subtree pull --squash -P ${mach_subtree_path} ${gmgit_dir} ${repo_br}
  [[ $(grep ${MACH_version} ${mach_subtree_path}/MANIFEST) == '' ]] && echo -e "ERROR: src/mach is not updated to ${MACH_version} version.\n" | tee -a ${gmtestinfo} && exit 1
 else
  git rm -rf src/mach
  git commit -m "Remove ${mach_subtree_path} subtree"
  git subtree add --squash -P ${mach_subtree_path} ${gmgit_dir} ${repo_br}
  [[ ! -d ./src/mach ]] && echo -e "ERROR: ${mach_subtree_path} is missing.\n" | tee -a ${gmtestinfo} && exit 1
 fi
 sed -i "s|^--prefix=${mach_subtree_path}.*|--prefix=${mach_subtree_path} ${gmgit_dir} ${repo_br}|" .git-subtree
 git add .git-subtree ; git commit -m "Update info on ${mach_subtree_path} source"
fi

# Check if the headers in MACH library have been updated
cd ${gemmach_dir}
mach_dir=${gemmach_dir}/src/mach
${mach_dir}/tools/regen_headers.sh -s ${mach_dir}/src/base
[[ -n $(diff -q ${mach_dir}/src/base/imported_headers/ ${mach_dir}/src/base/regenerated_headers/) ]] && (echo -e "ERROR: Heather modules in MACH library are not updated.\n" | tee -a ${gmtestinfo} ; echo -e "Review headers in ${mach_dir}/src/base/regenerated_headers directory.\n" | tee -a ${gmtestinfo}) && exit 1

# Load environment needed for compiling and running GEM 
#common_env_script=.eccc_setup_intel
common_env_script=.eccc_setup_intel_2025.1.0
source ${gemmach_dir}/${common_env_script}

# Determine hall of test work directory
if [[ -n "$(echo ${repo_dir} | grep /homeu3/ )" ]] ; then
 current_hall=hall7
else
 current_site=$(echo ${repo_dir} | cut -d "/" -f3)
 current_hall=$(echo "${current_site//site/hall}")
fi
# Set location of control directory 
control_dir=${control_dir:-/space/${current_hall}/sitestore/eccc/aq/r1/sarq000/gmtest/${Test_version}}
control_input_dir=${control_dir}/gm-input
[[ ! -d $control_input_dir ]] && echo "ERROR: control_input_dir does not exist: $control_input_dir" | tee -a ${gmtestinfo} && exit 1
echo -e "Regression-test input dir: ${control_input_dir}" | tee -a ${gmtestinfo}
# Link control input files to the test directory
mkdir -p ${TASK_INPUT}/cfg_0000
for fl in ${control_input_dir}/cfg_0000/* ; do 
 ln -s ${fl} ${TASK_INPUT}/cfg_0000/$(basename $fl)
done

## If local config files are provided, override config files e.g:
##  'physics_input_table->gm_phy_intable'; 'model_settings.nml->gem_settings.nml'
##  'output_settings->outcfg.out'
if [[ -d ${gmgit_dir}/tools/gemmach_cfg ]] ; then 
 if [[ -n $(ls ${gmgit_dir}/tools/gemmach_cfg) ]] ; then
  for file in $(ls ${gmgit_dir}/tools/gemmach_cfg/* | grep -v README.md) ; do
   [ -e "$file" ] || continue
   echo -e "override file: ${file} \n"
   cp -fv ${file} ${TASK_INPUT}
   filename=$(basename $file)
   [[ $filename == 'gem_settings.nml' ]] && ln -sfv ${TASK_INPUT}/${filename} ${TASK_INPUT}/cfg_0000/model_settings.nml
   [[ $filename == 'outcfg.out' ]] && ln -sfv ${TASK_INPUT}/${filename} ${TASK_INPUT}/cfg_0000/output_settings
   [[ $filename == 'gm_phy_intable' ]] && ln -sfv ${TASK_INPUT}/${filename} ${TASK_INPUT}/cfg_0000/physics_input_table
  done
 fi
fi

# For either debug compilation or running with tracing, namelist has to be modified locally
if [[ "${cmpl_opt}" == "dbg" ]] || [[ "${chm_trcng_opt}" == "on" ]] || [[ "${phy_trcng_opt}" == "on" ]] ; then
 [[ -L ${TASK_INPUT}/cfg_0000/model_settings.nml ]] && \
   cp --remove-destination `readlink ${TASK_INPUT}/cfg_0000/model_settings.nml` ${TASK_INPUT}/cfg_0000/model_settings.nml
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

# Copy regression-test related scripts into the test bin directory
if [[ "${ci_opt}" == "false" ]] ; then
 cp ${gmgit_dir}/tools/regression-test/compile-gm.sh ${TASK_BIN}
 cp ${gmgit_dir}/tools/regression-test/run-gm.sh ${TASK_BIN}
 cp ${gmgit_dir}/tools/regression-test/assemble-output.sh ${TASK_BIN}
 cp ${gmgit_dir}/tools/regression-test/validate-results.sh ${TASK_BIN}
fi

# Save configuration and resources to be used by all regression test scripts
cat << EOF >${TASK_BASEDIR}/common_env_arg
# Set information about the regression test
export Test_version=${Test_version}
export control_dir=${control_dir}
export gmtestinfo=${gmtestinfo}
export cmpl_opt=${cmpl_opt}
export cntrl_fl_opt=${cntrl_fl_opt}
export sequence_opt=${sequence_opt}
export ci_opt=${ci_opt}
# Set information about directory structure needed for regression test
export TASK_BIN=${TASK_BIN}
export TASK_INPUT=${TASK_INPUT}
export TASK_WORK=${TASK_WORK}
export TASK_OUTPUT=${TASK_OUTPUT}
export TASK_LIST=${TASK_LIST}
# Set the commonly use job resources
export Test_JobsMach=${TRUE_HOST}
export Test_JobsQueue=development
# Set environment for compiling and running GEM, and for working with fst files
cd ${gemmach_dir}
source ./${common_env_script}
EOF

# Save configuration and resources specific for compiling GEM-MACH 
 cat << EOF >${TASK_BASEDIR}/compile_env_arg
source ${TASK_BASEDIR}/common_env_arg
export GMRunJobProcTopo=10x8x1
export GMRunJobMemory=2G
[[ "${cmpl_opt}" == "dbg" ]] && export GMRunJobTime=3600 || export GMRunJobTime=1200
export GMRunJobName=rungm
EOF

# Save configuration and resources specific for running GEM-MACH
 cat << EOF >${TASK_BASEDIR}/run_env_arg
source ${TASK_BASEDIR}/common_env_arg
export GMRunJobProcTopo=10x8x1
export AssmblJobNcpu=1
export AssmblJobMemory=2G
export AssmblJobTime=300
export AssmblJobName=assmblgmo
EOF

# Save configuration and resources specific for assembling the output
 cat << EOF >${TASK_BASEDIR}/assemble_env_arg
source ${TASK_BASEDIR}/common_env_arg
export VldtJobNcpu=1
export VldtJobMemory=2G
export VldtJobTime=300
export VldtJobName=vldtgm
EOF

# Save configuration and resources specific for validating the output
 cat << EOF >${TASK_BASEDIR}/validate_env_arg
source ${TASK_BASEDIR}/common_env_arg
EOF

if [[ "${sequence_opt}" == "true" ]] ; then
 # Set the resources specific for compilation of GEM-MACH
 export Test_JobsMach=${TRUE_HOST}
 export Test_JobsQueue=development
 export CmplJobNcpu=1
 export CmplJobMem=160G
 export CmplJobTime=300
 export CmplJobName=cmplgm

 # Submit the GEM-MACH compilation script
 echo -e "\n == Submit the GEM-MACH compilation script at $(date) == \n"
 set -x
 ord_soumet ${TASK_BIN}/compile-gm.sh -args "${TASK_BASEDIR}" \
            -mach ${Test_JobsMach} -cpus ${CmplJobNcpu} -cm ${CmplJobMem} \
            -t ${CmplJobTime} -queue ${Test_JobsQueue} -jn ${CmplJobName} \
            -listing ${TASK_LIST}
 set +x

 # Inform the world about the regression test
 cat << EOF | tee -a ${gmtestinfo}

== Location of listings: ==

GEM-MACH compilation submitted to ${Test_JobsMach} at $(date).
Listings for the compilation job are in ${TASK_BASEDIR}/listings/cmplgm*out file.

Compilation job at the end initiates GEM-MACH run by submitting the job named "rungm" to ${TRUE_HOST}.
Listings for the GEM-MACH run are in ${TASK_BASEDIR}/listings/rungm*out file.

GEM-MACH run job at the end initiates assembly of GEM-MACH output by submitting the job named "assmblgmo" to ${TRUE_HOST}.
Listings for the assemble output job are in ${TASK_BASEDIR}/listings/assmblgmo*out file.

EOF

 if [[ "${cntrl_fl_opt}" != "new" ]] ; then
  cat << EOF | tee -a ${gmtestinfo}
At the end of assemble GEM-MACH output job, it initiates validation of GEM-MACH run results by submitting job named "vldt" to ${TRUE_HOST}.
Listings for the validation job are in ${TASK_BASEDIR}/listings/vldtgm*out file.

If the listings of vldtgm job contain "New binary reproduces reference output" phrase, regression test has successfully finished all stages.

EOF
 else
  cat << EOF | tee -a ${gmtestinfo}
Assemble GEM-MACH output job combines the model output to expected YYYYMMDDHH_hhh files, which can be used as a new control output.

If the listings of assmblgmo job contain "Location of new control output" phrase, regression test has successfully finished all requested stages.

EOF
 fi
fi

