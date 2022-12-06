#!/bin/bash

if [[ skip_this == x ]] ; then
   if [[ ! -d log ]] ; then
      mkdir -p build-${ORDENV_PLAT}/log/
      ln -s build-${ORDENV_PLAT}/log/ log
   fi
   mymemfn() {
      _file=$1
      _proc=$2
      _step=$3
      grep "${_step} stat" ${_file} | grep oe-${_proc} | cut -d';' -f2 | cut -d" " -f3
   }
   mymemloop() {
      _file=$1
      _proc=$2
      # Read every 3h, interpolate in between
      # 1h (00012), 3h (00036), 12h (00144)
      for _step in 00001 00035 00036 00037 00143 00144 ; do
         #mymemfn ${_file} ${_proc} ${_step}
         grep "${_step} stat" ${_file} | grep oe-${_proc} | cut -d';' -f2 | cut -d" " -f3 | tr '\n' ' '
      done
   }
   typelist="dist bloc gem48"
   # typelist="dist"
   topolist="1x1x1:1x1:1 6x6x1:1x1:1 6x6x1:3x2:6"
   # topolist="4x4x1:1x1:1"
   # topolist="2x1x1:1x1:1"
   # inorder="--inorder -v debug"
   inorder="--inorder"
   for item in ${typelist} ; do
      for np in ${topolist} ; do
         npl=($(echo ${np} | tr ":" " "))
         cmd="test_inputio8.sh -t ${item} --npe ${npl[0]} --nbloc ${npl[1]} --npio ${npl[2]} ${inorder}"
         logfile=log/test_inputio8.${item}.${np}.log
         echo + ${cmd}
         # ${cmd} --nsteps 1
         ${cmd} > ${logfile} 2>&1
         echo "${item}.${np} Time: $(grep real ${logfile})"
         echo "${item}.${np} Mem0: $(mymemloop ${logfile} 00000)"
         echo "${item}.${np} Mem1: $(mymemloop ${logfile} 00001)"
         echo emacs ${logfile} \&
         unset cmd logfile
      done
   done
   unset typelist topolist inorder mymemfn mymemloop

   # test_inputio8.sh --nsteps 1 --nhours 1 -t dist --npe 2x1x1 --nbloc 1x1 --npio 1 --inorder -v debug >log.log 2>&1
fi

eval $(rpy.cclargparse \
          -D ":" \
          ${0##*/} \
          "Run input tests" \
          "Timings will also be displayed" \
          "--nsteps"      'type=int'  '144'    '[Number of steps]'\
          "--nhours"      'type=int'  '12'     '[Number of hours (max 12)]'\
          "--testtype,-t" 'nargs=1'   'dist'   '[test type (dist, blod, gem48)]' \
          "--npe"         'nargs=1'   '1x1x1'  '[NPXxNPYxNOMP]'\
          "--nbloc"       'nargs=1'   '1x1'    '[NBLOCXxNBLOCY]'\
          "--npio"        'type=int'  '1'      '[NPE for Dist IO]'\
          "-v,--verbose"  'nargs=1'               'critical'  '[verbosity (debug, plus, info, warning, error, system, critical)]' \
          "+v"            'action=count, dest=v2' '0'         '[increase verbosity]' \
          "--debug,-d"    'action=store_true'  'false'  '[debug mode]' \
          "--inorder,-i"  'action=store_true'  'false'  '[output in order]' \
          "--keep,-k"     'action=store_true'  'false'  '[Keep input and test files]' \
             ++++ "$@")

          # "--noprep"      'action=store_true'  'false'  '[use previously prepared input if any]' \

TESTLIST=('dist' 'bloc' 'gem48')
VERBLVLS=('critical' 'error' 'warning' 'info' 'plus' 'debug')

if [[ ! "${TESTLIST[@]}" =~ "${testtype:-x}" ]] ; then
   cat 1>&2 <<EOF
ERROR: testtype must be one of: ${TESTLIST[@]}
EOF
   exit 1
fi
if [[ ! "${VERBLVLS[@]}" =~ "${verbose:-x}" ]] ; then
   cat 1>&2 <<EOF
ERROR: verbose must be one of: ${VERBLVLS[@]}
EOF
   exit 1
fi
if [[ ${nhours} -gt 12 ]] ; then
   cat 1>&2 <<EOF
ERROR: nhours must be <= 12
EOF
   exit 1
fi

export TEST_SUBNAME=${testtype}
export TEST_NSTEP=${nsteps}
export TEST_NHR=${nhours}
if [[ x"${verbose}" != x"" ]] ; then
   export TEST_VERBOSITY=${verbose}
fi
if [[ ${v2} -gt 0 ]] ; then
   [[ ${v2} -ge ${#VERBLVLS[@]} ]] && v2=$((${#VERBLVLS[@]} -1)) || true
   export TEST_VERBOSITY=${VERBLVLS[${v2}]}
fi
[[ x"${inorder}" == x"True" ]] && export TEST_INORDER=-inorder || true
if [[ x"${debug}" == x"True" ]] ; then
   export TEST_VERBOSITY=debug
   export TEST_DEBUG=-gdb
   export TEST_INORDER=
fi

#TODO: check ${npe} formating
npe123=($(echo ${npe} | tr 'x' ' '))
if [[ ${#npe123[@]} != 3 ]] ; then
      cat 1>&2 <<EOF
ERROR: npe must be of format: 9x9x9
EOF
   exit 1
fi
npex=${npe123[0]}
npey=${npe123[1]}
nomp=${npe123[2]}

#TODO: check ${nbloc} formating
npe12=($(echo ${nbloc} | tr 'x' ' '))
if [[ ${#npe12[@]} != 2 ]] ; then
      cat 1>&2 <<EOF
ERROR: nbloc must be of format: 9x9
EOF
   exit 1
fi
nblocx=${npe12[0]}
nblocy=${npe12[1]}

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


prepdata() {
   DFILES2=$1
   mkdir -p ${DFILES2}/geophy
   cp test_inputio8_table ${DFILES2}
   ln -s ${ATM_MODEL_DFILES}/bcmk/climato ${DFILES2}/
   ln -s ${ATM_MODEL_DFILES}/bcmk/geophy/Gem_geophy.fst ${DFILES2}/geophy
   ln -s ${ATM_MODEL_DFILES}/bcmk_toctoc/2009042700_000 ${DFILES2}/
   ln -s ${ATM_MODEL_DFILES}/bcmk_toctoc/2009042700_012 ${DFILES2}/
   for hr in 3 6 9 ; do
      myfile=${DFILES2}/2009042700_00${hr}
      #echo ==== $myfile
      cp ${ATM_MODEL_DFILES}/bcmk_toctoc/2009042700_000 ${myfile}
      rpy.fstzap -i ${myfile} ++nomvar P0 TT UU VV EN --zap datev=$(r.date -S -n 2009042700 +${hr}) >/dev/null
      #rpy.fstlist2 -i ${myfile} ++nomvar P0
   done
}

runtest() {
   set +x
   export MPI_NPEX=${1:-1}
   export MPI_NPEY=${2:-1}
   export MPI_NPEIO=${3:-1}
   export MPI_NBLOCX=${4:-1}
   export MPI_NBLOCY=${5:-1}
   export MPI_NDOM=1
   export MPI_IDOM=0
   redirect=${6}
   
   mydir=__to-delete_${testname}
   rm -rf ${mydir}
   mkdir ${mydir}
   cd ${mydir}
   
   echo "P=${MPI_NPEX}x${MPI_NPEY} NIO=${MPI_NPEIO}"
   set -x
   if [[ x"${redirect}" != x"" ]] ; then
      r.run_in_parallel ${TEST_DEBUG} ${TEST_INORDER} -pgm ${mybin} -npex $((${MPI_NPEX}*${MPI_NPEY})) -npey ${MPI_NDOM} > ${redirect} 2>&1
   else
      r.run_in_parallel ${TEST_DEBUG} ${TEST_INORDER} -pgm ${mybin} -npex $((${MPI_NPEX}*${MPI_NPEY})) -npey ${MPI_NDOM} 2>&1 | grep -v "label 1,idate2"
   fi
   set +x
   cd ..
   [[ x"${keep}" == x"False" ]] && rm -rf ${mydir} || true
}

DFILES2=${TMPDIR}/${testname}-$$
prepdata ${DFILES2}
export ATM_MODEL_DFILES=${DFILES2}

# runtest 1 1 1
# runtest 2 1
# runtest 1 3
# runtest 3 2

# runtest 2 2 2
# runtest 4 1 2
# runtest 1 4 2
# runtest 4 4 2
# runtest 4 4 4

#devnull="/dev/null"
set -x
time runtest ${npex:-1} ${npey:-1}  ${npio:-1}  ${nblocx:-1} ${nblocy:-1} ${devnull}
# time runtest 1 1  1  1 1 ${devnull}
# time runtest 2 1  1  1 1 ${devnull}
# time runtest 1 2  1  1 1 ${devnull}
# time runtest 2 2  1  1 1 ${devnull}
# time runtest 2 2  2  1 1 ${devnull}
# time runtest 4 1  1  1 1 ${devnull}
# time runtest 1 4  1  1 1 ${devnull}
# time runtest 4 2  1  1 1 ${devnull}
# time runtest 2 4  1  1 1 ${devnull}
# time runtest 4 2  2  1 1 ${devnull}
# time runtest 2 4  2  1 1 ${devnull}
# time runtest 4 4  1  1 1 ${devnull}
# time runtest 4 4  2  1 1 ${devnull}
# time runtest 4 4  4  1 1 ${devnull}

## runtest 3 1 1
## runtest 1 3 1
## runtest 3 3 1
## runtest 3 3 2
## runtest 3 3 3

[[ x"${keep}" == x"False" ]] && rm -rf ${DFILES2} || true
