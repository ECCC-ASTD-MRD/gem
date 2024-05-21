#!/bin/bash
#
arguments=$*
. r.entry.dot
#
eval `cclargs_lite $0 \
     -s        ""      ""       "[source                                      ]"\
     -liste    ""      ""       "[Liste of model output type to treat         ]"\
     -nthreads "1"     "1"      "[Number of yyoutg process to run in parallel ]"\
     -_status  "ABORT" "ABORT"  "[return status ]"\
  ++ $arguments`
#
printf "\n=====>  Um_output_yyoutg.sh starts: `date` ###########\n\n"

ici=`pwd`

laliste=""
cd $s
for i in $liste ; do
  laliste=${laliste}" "`find ./ -name "${i}[0-9]*" | sed 's/^..//' | xargs`
done

count=0
for i in ${laliste} ; do
  (( count = count + 1 ))
  echo Processing $i in background $count $nthreads
  lis=${TASK_WORK}/yy2e_${count}

  tmpfile=ufile_${count}
  mv ./${i} $tmpfile
  ${TASK_BIN}/yy2global -iyy $tmpfile -oglob ${i} 1> ${lis}.lis 2> ${lis}.err &

  if [ $count = $nthreads ]; then
    echo WAITING for end of current threads $(date)
    wait
    echo WAITING DONE $(date)
    count=0
  fi
done

echo WAITING for last threads $(date)
wait
echo WAITING DONE $(date)

/bin/rm -f ufile_*

cd ${ici}

_status='OK'

printf "\nUm_output_yyoutg.sh ends: `date` \n\n"

# End of task
. r.return.dot
