#!/bin/bash

# Script to clean GEM-MACH code directories from the files and directories built in CHEM ang GEM directories during compilation.
# The script will not clean the ${storage_model} and /home/ords/.../storage_model files and directories.
# The script is to be called from the GEM-MACH's code root directory (not CHEM and not GEM).
# Author: Verica Savic-Jovcic
# Date:   Nov. 2019

# Set the colors for display:
fcred=`tput setaf 1`
fcyellow=`tput setaf 3`
fccyan=`tput setaf 6`
clrreset=`tput sgr0`

# Set directory to be cleaned:
cgmcd=`pwd`
echo "${fcred} Do you want to clean ${cgmcd} directory? [y/n] ${clrreset}"
read cgmcdopt

if [ "${cgmcdopt}" = "y" ] ; then

# Clean CHEM directory:
 cd ${cgmcd}/CHEM
 echo "${fcyellow} Clean `pwd` directory ${clrreset}"
 rm -r build-* Makefile.user.root.mk Makefile.rules.mk Makefile.dep.Linux_x86-64.mk Makefile.build.mk Makefile .rde.* share/ lib/ include/.restricted bin/

# Clean GEM directory:
 cd ${cgmcd}/GEM
 echo "${fccyan} Clean `pwd` directory ${clrreset}"
 rm -r .linkit.log .rde.* BINMOD Makefile Makefile.[bdr]* Makefile.user.root.mk PREP RUNMOD bin/ build-* include/ lib/ share/

else

 exit 1

fi

