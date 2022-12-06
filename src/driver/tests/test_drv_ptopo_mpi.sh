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
   export MPI_NDOMS=$1
   export MPI_NGRIDS=$2
   export MPI_NPEX=$3
   export MPI_NPEY=$4
   export MPI_NBLOCX=${5-$MPI_NPEX}
   export MPI_NBLOCY=${6-$MPI_NPEY}
   cat > ${setfilebase}.cfg <<EOF
version=100
@ptopo_cfgs
npx = $MPI_NPEX
npy = $MPI_NPEY
@
EOF

   rm -rf $testdir
   mkdir $testdir
   cd $testdir

   export UM_EXEC_NDOMAINS=1
   [[ x$MPI_NDOMS == x2 ]] && export UM_EXEC_NDOMAINS=2:3
   export UM_EXEC_NGRIDS=$MPI_NGRIDS
   export UM_EXEC_CONFIG_BASENAME=${setfilebase##*/}
   export UM_EXEC_CONFIG_DIR=${setfilebase%/*}
   export TASK_WORK=$testdir/work
   export TASK_INPUT=$testdir/input
   export TASK_OUTPUT=$testdir/output

   echo "D=${UM_EXEC_NDOMAINS} G=${UM_EXEC_NGRIDS} P=${MPI_NPEX}x${MPI_NPEY} B=${MPI_NBLOCX}x${MPI_NBLOCY}"
   ((MPI_NPEX =MPI_NDOMS * MPI_NGRIDS * MPI_NPEX))
   r.mpirun -pgm $mybin -npex $MPI_NPEX -npey $MPI_NPEY

   rm -f ${setfilebase}.cfg
}

#runtest 1 1 1 1
#runtest 1 1 2 1
#runtest 1 1 1 3
#runtest 1 2 1 1
#runtest 1 2 2 1
#runtest 2 1 1 1
#runtest 2 1 2 1
#runtest 2 2 1 1
runtest 2 2 2 1

cd /
rm -rf $testdir
