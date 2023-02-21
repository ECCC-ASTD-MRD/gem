#!/bin/bash
#
rep2check=${1:-./}
fs_c=FSFS`df $rep2check 2> /dev/null | awk '{print $NF}'`
fs_h=FSFS`df ${HOME}    2> /dev/null | awk '{print $NF}'`

if [ -d $rep2check ] ; then
  if [[ $fs_c == $fs_h ]] ; then
     printf "\n Filesystem can not be $rep2check -- ABORT --\n\n"
     echo rep2check=${rep2check} : $(true_path ${rep2check})
     df ${rep2check}
     df ${HOME}
    exit 1
  fi
else
  printf "\n Directory $rep2check is unavailable -- ABORT --\n\n"
  exit 1
fi
printf "\n Directory $rep2check is OK (not in your \$HOME) \n\n"
df $rep2check
exit 0
