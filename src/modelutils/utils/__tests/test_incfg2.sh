#!/bin/bash

if [[ x"$1" == x-h ]] ; then
   echo "USAGE: $0 [-v [LEVEL]]"
   exit
fi

export TEST_INORDER=
export TEST_VERBOSITY=critical
export TEST_DEBUG=
if [[ x"$1" == x-inorder ]] ; then
   export TEST_INORDER=-inorder
   shift
fi

if [[ x"$1" == x-d ]] ; then
   export TEST_VERBOSITY=debug
   export TEST_DEBUG=-gdb
   shift
fi

if [[ x"$1" == x-v ]] ; then
   export TEST_VERBOSITY=info
   [[ x"$2" != x ]]&& export TEST_VERBOSITY=$2
fi

testname=${0%.*}
testname=${testname##*/}
mybin=$(pwd)/malib${EC_ARCH}/${testname}.mpiAbs
[[ ! -x $mybin ]] && \
   mybin=$(pwd)/${EC_ARCH}/build/${testname}.mpiAbs
[[ ! -x $mybin ]] && \
   mybin=$(pwd)/$(rdevar build/bin)/${testname}.Abs

if [[ ! -x $mybin ]] ; then
   echo "ERROR: ${testname}.mpiAbs Not Found"
   echo "==== Abort ===="
   exit 1
fi


runtest() {
   export MPI_NPEX=${1:-1}
   export MPI_NPEY=${2:-1}
   export MPI_NPEIO=${3:-1}
   export MPI_NBLOCX=1
   export MPI_NBLOCY=1
   export MPI_NDOM=1
   export MPI_IDOM=0

   echo "P=${MPI_NPEX}x${MPI_NPEY} NIO=${MPI_NPEIO}"
   r.run_in_parallel ${TEST_DEBUG} ${TEST_INORDER} -pgm ${mybin} -npex ${MPI_NPEX} -npey ${MPI_NPEY}
}

runtest 1 1 1
