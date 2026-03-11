#!/bin/bash

set -ex

echo "We are in $0 -> $(true_path ${0})"

which touchup_analysis.sh

editfst -s $1 -d ${1}_tmp -i 0
touchup_analysis.sh -i ${1}_tmp -o $2 -zapdate "TM LG" -sd -pgsm pgsm
#touchup_analysis.sh -i $1 -o $2 -var2chk P0 -zapdate "TM LG" -pgsm pgsm
