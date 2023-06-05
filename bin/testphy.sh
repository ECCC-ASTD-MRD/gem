#!/usr/bin/env bash

# Establish temporary directory name
basedir=${TMPDIR}/testphy

# Retrieve command line arguments
args=$*
eval `cclargs_lite $0 \
  -debug    ""           "valgrind"                        "[run with a debugger            ]"\
  -supp     ""           "${basedir}/input/valgrind.supp"  "[debug error suppression file   ]"\
  -cfg      "${rpnphy}/share/testphy.cfg"    ""            "[path to test configuration     ]"\
  -o        "${basedir}/output/listing.txt"  "STDOUT"      "[print output to the screen     ]"\
  ++ $args`

# Preliminary setup
set -e
if [[ -n "${debug}" || ${o} == "STDOUT" ]] ; then
  . task_setup.dot -v --clean --base=${basedir} --file=${cfg}
else
  build_tmp=${TMPDIR}/testphy.cfgout.$$
  . task_setup.dot -v --clean --base=${basedir} --file=${cfg} >${build_tmp} 2>&1
  mv ${build_tmp} ${o}
fi

# Run test
if [[ "${debug}" == "valgrind" ]] ; then
  vsupp=""
  if [[ -n "${supp}" ]] ; then vsupp="--suppressions=${supp}" ; fi
  large_stack=20000000
  ${debug} --main-stacksize=${large_stack} --max-stackframe=${large_stack} ${vsupp} ${TASK_BIN}/testphy.Abs
else
  if [[ -n "${debug}" || ${o} == "STDOUT" ]] ; then
    ${debug} ${TASK_BIN}/testphy.Abs
  else
    ${debug} ${TASK_BIN}/testphy.Abs >${o} 2>&1
    grep 'CRITICAL ERROR' ${o} 2>/dev/null || true
    grep 'ABORTING' ${o} 2>/dev/null || true
    grep '*FAIL' ${o} 2>/dev/null || true
    passed=$(grep '*PASS' ${o} 2>/dev/null || true)
    if [[ -z "${passed}" ]] ; then
      printf "Test Failed.  Rerun with the following options to see output:\n  $0 ${args} -o\n"
      exit 1
    fi
    printf "${passed}\n"
  fi
fi

