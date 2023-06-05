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

testname=${0%.*}
testname=${testname##*/}
mybin=$(pwd)/malib${EC_ARCH}/${testname}.mpiAbs
[[ ! -x $mybin ]] && \
   mybin=$(pwd)/${EC_ARCH}/build/${testname}.mpiAbs

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

   echo "P=${MPI_NPEX}x${MPI_NPEY} B=${MPI_NBLOCX}x${MPI_NBLOCY}"
   r.mpirun -pgm $mybin -npex $MPI_NPEX -npey $MPI_NPEY
}

#runtest 1 1
#runtest 2 1
#runtest 1 3
runtest 3 2

