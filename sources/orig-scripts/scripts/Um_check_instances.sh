#!/bin/bash
#
if [ $1 -ne $2 ] ; then
  printf "\n Mismatch between required and provided # of instances:\n"
  printf "      provided: $1\n"
  printf "      required: $2\n"
  printf " ------ ABORT ------\n\n"
  exit 1
fi

exit 0
