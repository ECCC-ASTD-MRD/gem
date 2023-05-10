#!/bin/bash -ex

launch_output_job ()
{
printf "\n ===> launch_output_job: STARTS $(date)\n"
set -ex

REP_OUTPUT_SRC=$(readlink -e ${2}/$(echo $1 | sed 's/\^last//g'))
stepnumber=$(basename $REP_OUTPUT_SRC)
ID_LOJ=${3}

DIRname=$(dirname $REP_OUTPUT_SRC)
DOMAIN=$(basename $DIRname)
config_file=configexp_total_${ID_LOJ}

touch ${config_file}

STATUS_file=${2}/status_MOD.dot
if [ -s ${STATUS_file} ] ; then
  . ${STATUS_file}
  cat ${STATUS_file} > ${config_file}
fi

for i in ${DOMAINS_this_instance} ; do
   dir=${DIRname}/../$i
   cp ${config_file} ${dir}/${stepnumber}/configexp.cfg
done

start_DOM=$(echo $GEM_NDOMAINS | cut -d : -f1)
cnt=$((start_DOM-DOMAIN_start))
cnt=$((cnt/DOMAIN_wide))
domain_cnt=$(printf "%02d" $((cnt+1)) )

NODE=$(basename ${SEQ_NODE:-'undefined'})
MODULE=$(basename $(dirname ${SEQ_NODE:-'undefined'}))
if [ $(tty -s;echo $?) -ne 0 -a -n "${SEQ_EXP_HOME}" -a "${MODULE}" == "${SEQ_MODULE:-'NoNe'}" -a "${NODE}" == "Runmod" ] ; then
  name=S${domain_cnt}_$(echo ${1} | cut -f 2 -d "_" )
  if [[ -n "$SEQ_CONTAINER_LOOP_ARGS" ]]; then
       $SEQ_BIN/maestro -n $(dirname $SEQ_NODE)/Sortie -s submit ${SEQ_CONTAINER_LOOP_ARGS},Sortie=${name}
  else 
       $SEQ_BIN/maestro -n $(dirname $SEQ_NODE)/Sortie -s submit -l Sortie=${name}
  fi
fi
printf "\n ===> launch_output_job: TERMINATE $(date)\n"
}

process_output ()
{
printf "\n ===> Process_output: $1 STARTS $(date)\n"
set -ex
file2process=${1}
OUTREP=${2}
ID_po=${3}
OUT2PRC=output2process_${ID_po}
cat ${file2process} | grep NORMAL_OUTPUT | cut -d " " -f3 > ${OUT2PRC}
cat ${OUT2PRC}
file_done=${OUTREP}/output_done

#printf "\n WAITING FOR A LOCK on file ${file_done} $(date)\n"
#lockfile ${file_done}.lock
#date

unset laliste
if [ -s ${OUT2PRC} ] ; then
  if [ -s $file_done ] ; then
    cat $file_done
    for i in $(cat ${OUT2PRC}) ; do
      if [ $(grep $i $file_done | wc -l) -eq 0 ] ; then
        laliste=${laliste}" "$i
      fi
    done
  else
    laliste=$(cat ${OUT2PRC})
  fi
fi

echo LALISTE: $laliste

for REP in $laliste ; do
  printf "\n ===> Process_output: treating $REP $(date)\n"
  time launch_output_job ${REP} ${OUTREP} ${ID_po}
  echo $REP >> $file_done
done
echo wait STARTS $(date)
wait
echo wait DONE $(date)

mv ${file2process} $(echo ${file2process} | cut -d"." -f 1).done
printf "\n ===> Process_output: TERMINATED $(date)\n"
}

# MAIN SCRIPTS

if [ ! -s ${1:-pasdefichier} ]; then exit; fi

printf "\n ===> ########## Um_process_output STARTS $(date) ##########\n"

rep=${TASK_WORK}/post_process_output_$(basename $(dirname ${1}))
cd $rep ; id=$$
fn=$(basename ${1})_work_${id}
cp ${1} ${fn}.active
lis=${fn}.lis
date
process_output ${fn}.active $(dirname ${1}) ${id} 1> ${lis} 2>&1
# &

printf "\n ===> ########## Um_process_output DONE $(date) ##########\n"
