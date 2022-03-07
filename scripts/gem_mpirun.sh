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
   echo "exporting OMP_STACKSIZE=4G"
   export OMP_STACKSIZE=4G
   echo "setting ulimit stack to unlimited in lance.sh for all PEs"
   echo "which mpirun"
   which mpirun
   mpirun -np $((npex*npey)) lance.sh ${pgm}
fi
