#!/bin/sh
#
# Store command line arguments
arguments=$*
printf "$0 ${arguments}\n\n"

# Process command line arguments
eval `cclargs_lite -D " " $0 \
   -pgm     ""      ""   "[Program name]"\
   -npex    "1"     "1"  "[]"\
   -npey    "1"     "1"  "[]"\
   -tag       ""     ""  "[]"\
   -gdb       ""     ""  "[]"\
   -inorder   ""     ""  "[]"\
   -minstdout ""     ""  "[]"\
   -nocleanup ""     ""  "[]"\
  ++ ${arguments}`

which r.run_in_parallel 2> /dev/null
if [ $? == 0 ] ; then
   r.run_in_parallel ${arguments}
else
   mpirun -np $((npex*npey)) ${pgm}
fi
