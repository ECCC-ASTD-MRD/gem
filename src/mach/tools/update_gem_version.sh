#!/bin/bash

# update_gem_version.sh
#
# Script to facilitate update to newer GEM version by
#        updating GEM version numbers in all relevant files
#
# Author: Verica Savic-Jovcic (25/11/2019)
#
# Usage: TOOLS/update_gem_version.sh <current gem version> <new gem version>
#        where:
#               <current gem version> and <new gem version> must contain
#               "x\/" or "test\/" if GEM version is not in GEM/ directory;
#               "\" is necessary in front of "/" !


# Make sure that all necesary command-line arguments are provided:

shift $((OPTIND -1))
if [ "$#" -ne 2 ]; then
 echo "Restart the script with correct GEM versions:"
 echo "TOOLS/update_gem_version.sh <current gem version> <new gem version>"
 exit 1
fi

# Read in GEM versions
gem_v=$1
new_gem_v=$2

# Replace GEM versions in all necessary files
echo "Update GEM version from ${gem_v} to ${new_gem_v}? [y/n] "
read gvopt
if [ "${gvopt}" = "y" ] ; then
 sed -i "s|${gem_v}|${new_gem_v}|g" README_compilation.md
 sed -i "s|${gem_v}|${new_gem_v}|g" README_version.md
 sed -i "s|${gem_v}|${new_gem_v}|g" README_special.md
 sed -i "s|${gem_v}|${new_gem_v}|g" README.md
 vi README_version.md
else
 echo "Restart the script with correct GEM versions:"
 echo "TOOLS/update_gem_version.sh <current gem version> <new gem version>"
 exit 1
fi
