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
testdir=$(pwd)/__${testname}-to-delete__
mybin=$(pwd)/malib${EC_ARCH}/${testname}.mpiAbs
[[ ! -x $mybin ]] && \
   mybin=$(pwd)/${EC_ARCH}/build/${testname}.mpiAbs

if [[ ! -x $mybin ]] ; then
   echo "ERROR: ${testname}.mpiAbs Not Found"
   echo "==== Abort ===="
   exit 1
fi

setfilebase=$(pwd)/${testname}

runtest() {
   export MPI_NPEX=$1
   export MPI_NPEY=$2
   export MPI_NBLOCX=${3-$MPI_NPEX}
   export MPI_NBLOCY=${4-$MPI_NPEY}

   mysubdir=${MPI_NPEX}x${MPI_NPEY}x${MPI_NBLOCX}x${MPI_NBLOCY}
   mkdir $mysubdir
   cd $mysubdir
   ln -s ${setfilebase}.out .
   rm -rf __test_outout_mpi-to-delete__
   echo "P=${MPI_NPEX}x${MPI_NPEY} B=${MPI_NBLOCX}x${MPI_NBLOCY}"
   r.mpirun -pgm $mybin -npex $MPI_NPEX -npey $MPI_NPEY
   cd ..
}

here=$(pwd)
testdir=$(pwd)/__test_outout_mpi-to-delete__
mkdir $testdir
cd $testdir

#runtest 1 1
#runtest 2 1
#runtest 1 3
#runtest 3 2
runtest 3 2 1 1

#runtest 3 1 1 1
#runtest 1 2 1 1
#runtest 2 3 1 1

#cd $here
#rm -rf $testdir
