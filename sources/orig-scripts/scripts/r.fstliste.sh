#!/bin/bash
# set -x
eval `cclargs_lite0 $0 "[ Utilitaire pour extraires des parametres des fichiers standard]"\
 -izfst bidon "[rep et nom du fichier standard]"\
 -datev =-1 "[date de validite (stamp)]"\
 -vdatev =-1 "[date de validite (true-date)]"\
 -etiket ' ' "[Etiket]"\
 -ip1 =-1 "[Niveaux]"\
 -ip2 =-1 "[Heures]"\
 -ip3 =-1 "[IP3]"\
 -typvar ' ' "[Type de variable]"\
 -nomvar ' ' "[Nom de la variable]"\
 -col 0 "1-20" "[colomne(s) desiree(s)]"\
 -liste 1 0 "[sortie sous forme de liste]"\
 -unique 1 0 "[enleve les selections identiques]"\
 -champs 0 1 "[Liste les colonnes et variables]"\
 -nosort 0 1 "[Ne pas trier les sorties lorsque unique = 1]"\
 ++ $*`
#
#set -ex
if test ${champs} -eq 1
then
  r.fstinfo0 -champs
  exit
fi
#
list=""

#r.fstinfo0 -izfst ${izfst}
#echo $? jljl
#ls -l
#exit

r.fstinfo0 -izfst ${izfst} -datev =${datev} -vdatev =${vdatev} -etiket "${etiket}" -ip1 =${ip1} -ip2 =${ip2} -ip3 =${ip3} -typvar "${typvar}" -nomvar "${nomvar}" \
 -otxt ${TMPDIR}/fstlist.$$ > /dev/null 2>&1
cat ${TMPDIR}/fstlist.$$
exit
integer infoexit=$?
echo $? jljl
if [ $infoexit -ne 0 ]
then
    echo ${list}
    exit $infoexit
fi
#
if test "${col}" != "0"
then
  if test ${unique} -ne 0
  then
    if [ ${nosort} -eq 1 ] ; then
      lastarg=""
      for arg in `cut -f${col} -d':' ${TMPDIR}/fstlist.$$ | sed 's/:/ /g'` ; do
        if [ "${arg}" != "${lastarg}" ] ; then
          if [ ${liste} -eq 0 ] ; then
            echo ${arg}
          else
            list="${list} ${arg}"
          fi
          lastarg=${arg}
        fi
      done
      if [ ${liste} -eq 1 ] ; then echo ${list} ; fi
    else
      if test ${liste} -eq 0
      then
        cut -f${col} -d':' ${TMPDIR}/fstlist.$$ | sed 's/:/ /g' | sort -u
      else
        for i in `cut -f${col} -d':' ${TMPDIR}/fstlist.$$ | sed 's/:/ /g' | sort -u`
        do
          list="${list} $i"
        done
        echo ${list}
      fi
    fi
  else
    if test ${liste} -eq 0
    then
      cut -f${col} -d':' ${TMPDIR}/fstlist.$$ | sed 's/:/ /g'
    else
      for i in `cut -f${col} -d':' ${TMPDIR}/fstlist.$$ | sed 's/:/ /g'`
      do
        list="${list} $i"
      done
      echo ${list}
    fi
  fi
else
  cat ${TMPDIR}/fstlist.$$
fi

rm ${TMPDIR}/fstlist.$$
exit 0

