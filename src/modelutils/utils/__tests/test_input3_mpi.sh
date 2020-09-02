#!/bin/bash

if [[ x"$1" == x-h ]] ; then
   echo "USAGE: $0 [-v [LEVEL]]"
   exit
fi

export TEST_VERBOSITY=critical
if [[ x"$1" == x-v ]] ; then
   export TEST_VERBOSITY=info
   [[ x"$2" != x ]]&& export TEST_VERBOSITY=$2
fi

bindir=$(rdevar build/bin)
testname=${0%.*}
testname=${testname##*/}
mybin=$(pwd)/${bindir}/${testname}.Abs

if [[ ! -x $mybin ]] ; then
   echo "ERROR: ${testname}.mpiAbs Not Found"
   echo "==== Abort ===="
   exit 1
fi


runtest() {
   export MPI_NPEX=$1
   export MPI_NPEY=$2
   export MPI_NBLOCX=${3-$MPI_NPEX}
   export MPI_NBLOCY=${4-$MPI_NPEY}

   rmpirun="$(which r.run_in_parallel || true)"
   [[ x$rmpirun == x ]] && rmpirun="r.mpirun"

   echo "P=${MPI_NPEX}x${MPI_NPEY} B=${MPI_NBLOCX}x${MPI_NBLOCY}"
   $rmpirun -pgm $mybin -npex $MPI_NPEX -npey $MPI_NPEY
}

#runtest 1 2
runtest 1 1
#runtest 2 4
#runtest 4 2

# runtest 1 1
# runtest 2 1
# runtest 1 3
# runtest 3 2

# runtest 3 1 1 1
# runtest 1 2 1 1
# runtest 2 3 1 1
