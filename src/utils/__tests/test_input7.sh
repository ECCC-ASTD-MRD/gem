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
   set +x
   export MPI_NPEX=${1:-1}
   export MPI_NPEY=${2:-1}
   export MPI_NBLOCX=${3-$MPI_NPEX}
   export MPI_NBLOCY=${4-$MPI_NPEY}
   export MPI_NBLOCX=1
   export MPI_NBLOCY=1
   export MPI_NDOM=1
   export MPI_IDOM=0

   echo "P=${MPI_NPEX}x${MPI_NPEY} B=${MPI_NBLOCX}x${MPI_NBLOCY}"
   r.run_in_parallel ${TEST_DEBUG} ${TEST_INORDER} -pgm ${mybin} -npex ${MPI_NPEX} -npey ${MPI_NPEY}
}


set -x
runtest 1 1
# runtest 2 1
# runtest 1 2
# runtest 2 2
# runtest 3 1 1 1
# runtest 1 2 1 1
# runtest 2 3 1 1
