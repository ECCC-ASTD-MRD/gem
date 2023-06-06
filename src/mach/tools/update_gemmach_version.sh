#!/bin/bash

# update_gemmach_version.sh
#
# Script to facilitate update to newer GEM-MACH version by
#        updating GEM-MACH version numbers in all relevant files
#
# Author: Verica Savic-Jovcic (25/11/2019)
#
# Usage: TOOLS/update_gemmach_version.sh <current GEM-MACH version> <new GEM-MACH version>
#        where:
#               <current GEM-MACH version> and <new GEM-MACH version> must contain
#               "x\/" or "test\/" if GEM-MACH version is not in GEM/ directory;
#               "\" is necessary in front of "/" !


# Make sure that all necesary command-line arguments are provided:

shift $((OPTIND -1))
if [ "$#" -ne 2 ]; then
 echo "Restart the script with correct GEM-MACH versions:"
 echo "TOOLS/update_gemmach_version.sh <current GEM-MACH version> <new GEM-MACH version>"
 exit 1
fi

# Read in GEM-MACH versions
gemmach_v=$1
new_gemmach_v=$2

# Replace GEM-MACH versions in all necessary files
echo "Update GEM-MACH version from ${gemmach_v} to ${new_gemmach_v}? [y/n] "
read gvopt
if [ "${gvopt}" = "y" ] ; then
 vi README_version.md
 sed -i "s|${gemmach_v}|${new_gemmach_v}|g" .gitlab-ci.yml
 sed -i "s|${gemmach_v}|${new_gemmach_v}|g" TOOLS/gm-integration-test/initialize-gm-integration-test.sh
 grep ${new_gemmach_v} README_version.md .gitlab-ci.yml TOOLS/gm-integration-test/initialize-gm-integration-test.sh
else
 echo "Restart the script with correct GEM-MACH versions:"
 echo "TOOLS/update_gemmach_version.sh <current GEM-MACH version> <new GEM-MACH version>"
 exit 1
fi
