#!/bin/bash
#
arguments=$*
. r.entry.dot

#====> Obtaining the arguments:
eval `cclargs_lite $0 \
     -input    ""      ""       "[Input directory                    ]"\
     -output   ""      ""       "[Output directory                   ]"\
     -assemble "0"     "1"      "[Reassemble or not                  ]"\
     -dplusp   "0"     "1"      "[Combine dynamics and physics output]"\
     -liste    ""      ""       "[Liste of model output type to treat]"\
     -listm    ""      ""       "[Xfer model listings along          ]"\
     -repcasc  ""      ""       "[Xfer cascade files                 ]"\
     -serdate  "0"     "1"      "[Keep date suffix for time series files]"\
     -_status  "ABORT" "ABORT"  "[return status                      ]"\
     -nthreads "1"     "1"      "[Number of bemol process to run in parallel]"\
  ++ $arguments`

set ${SETMEX:-+ex}

printf "\n    =====> `basename $0` $arguments: `date`\n"

ici=`pwd`

nbfiles=0

##### We first deal with subdirectories ${input}/[0-9]*-[0-9]* #####

# casc_* cascade files
laliste="casc_"
if [ -n "${repcasc}" ] ; then
   find -L ${input} -type f -name "casc_*" > casc_files_found
   cnt=`cat casc_files_found | wc -l`
   if [ $cnt -gt 0 ] ; then
     REP=${output}/${repcasc}
     mkdir ${REP}
     for i in `cat casc_files_found` ; do
       ln -s $i ${REP}
     done
   fi
   nbfiles=$((nbfiles+cnt))
fi

# Process all types of model output
printf "\n  =====> MODOUT starts ${nbfiles} `date`\n"
. r.call.dot ${TASK_BIN}/Um_output_modout.sh -src ${input} -dst ${output} \
    -liste ${liste} -assemble ${assemble} -dplusp ${dplusp} -nthreads ${nthreads}
nbfiles=$((nbfiles+_nf2t))
assemble_error=${_nerr}
flag_err=${assemble_error}
printf "  =====> MODOUT ends: $nbfiles files treated: $assemble_error errors detected `date`\n"
laliste=${laliste}" dm pm dp pp dh ph"

##### Then we deal with other special files
prefix=''

# BUSPER4spinphy file to recycle BUSPER

if [ -s ${input}/endstep_misc_files/BUSPER4spinphy* ] ; then
  mv ${input}/endstep_misc_files/BUSPER4spinphy* ${output}
fi
laliste=${laliste}" BUSPER4spinphy"

# Restart files
find -L ${input} -type d -name "restart*" > restart_files_found
cnt=`cat restart_files_found | wc -l`
if [ $cnt -gt 0 ] ; then
  for i in `cat ${ici}/restart_files_found` ; do
     restart_name=$(basename ${i})
     cd ${i}
     tar cvf ${output}/${restart_name}.tar ${prefix}*
  done
  nbfiles=$((nbfiles+cnt))
  cd ${ici}
fi
laliste=${laliste}" restart"

# Time series file
find -L ${input} -type f -name 'time_series*.bin*' > series_found
cnt=`cat series_found | wc -l`
series_ok=1
if [ $cnt -gt 0 ] ; then
  in_serie=`cat series_found`
  #ln -sf ${in_serie} .
  printf "\n  =====>  feseri -iserial time_series.bin -omsorti time_series.fst \n"
  for item in ${in_serie} ; do
     if [[ ! -s ${item} ]] ; then
        echo "Skipping empty file: ${item}"
        ls -l ${item}
        rm -f ${item}
        continue
     fi
     itemnum=${item##*_}_YIN
     if [[ x$(echo ${item} | sed 's|/YAN/|/|') != x${item} ]] ; then
        itemnum=${item##*_}_YAN
     fi
     if [[ ${serdate} == 0 ]] ; then
        itemout=time_series.fst_${itemnum}
     else
        itembase=${item%_*}
        itemout=${itembase%.bin}.fst_${itemnum}
     fi
     itemout=${itemout##*/}
     itemin=time_series.bin_${itemnum}
     ln -s ${item} ${itemin}
     # itemin=${item}
     # printf "\n  =====>  feseri -iserial ${in_serie} -omsorti ${itemout} \n"
     
     ${TASK_BIN}/feseri -iserial ${itemin} -omsorti ${itemout} 1> ${itemin}.lis
     if [ $? -ne 0 ]; then
        printf "\n feseri aborted \n\n"
        /bin/rm -f ${itemout} ${itemin}
        flag_err=$((flag_err+1))
        series_ok=0
        break
     else
        /bin/rm -f ${item} ${itemin} 2> /dev/null || true
     fi
  done
  if [[ $series_ok == 1 ]] ; then
     
     if [[ ${serdate} == 0 ]] ; then

        for yinyan in YIN YAN; do
           if [[ x"$(ls time_series.fst_*_${yinyan} 2>/dev/null)" != x"" ]] ; then
              editfst -i 0 \
                 -s time_series.fst_*_${yinyan} \
                 -d time_series.fst_${yinyan}
              mv time_series.fst_${yinyan} ${output}
              nbfiles=$((nbfiles+1))
           fi
        done
        mv ${output}/time_series.fst_YIN \
           ${output}/time_series.fst 2>/dev/null || true

     elif [[ x"$(ls time_series_*.fst_*_* 2>/dev/null)" != x"" ]] ; then

        datelist="$(ls time_series_*.fst_*_* 2>/dev/null | sed 's/time_series_//' | cut -d'.' -f1-2 | tr ' ' '\n' | sort -u | tr '\n' ' ')"
        for serdate in ${datelist} ; do
           for yinyan in YIN YAN; do
            if [[ x"$(ls time_series_${serdate}.fst_*_${yinyan} 2>/dev/null)" != x"" ]] ; then 
              editfst -i 0 \
                 -s time_series_${serdate}.fst_*_${yinyan} \
                 -d time_series_${serdate}_${yinyan}.fst
              mv time_series_${serdate}_${yinyan}.fst ${output}
              nbfiles=$((nbfiles+1))
            fi    
           done
           mv ${output}/time_series_${serdate}_YIN.fst \
              ${output}/time_series_${serdate}.fst 2>/dev/null || true
        done
     fi

  fi
fi
laliste=${laliste}" time_series"

# Listings files
if [ "${listm}" ] ; then
  lalistetomv=`ls ${listm}* 2> /dev/null | xargs`
  if [ -n "$lalistetomv" ] ; then
    cp ${lalistetomv} ${output} 2> /dev/null
    nbfiles=$((nbfiles+1))
  fi
fi

# other untreated files
for i in ${laliste} config ; do
   echo /'\'/$i/d >> commands.txt
done
ici=$PWD ; cd ${input}
find -L ./ -type f > ${ici}/listoffiles$$
cd $ici
sed -f commands.txt listoffiles$$ > list_of_untreated_files$$
while read line ; do
   mv ${input}/$line ${output}
done < list_of_untreated_files$$

if [ $flag_err -eq 0 ] ; then
  _status='OK'
  cnt=`ls -l ${output} | wc -l`
  cnt=$((cnt-1))
  printf "\n=====> Um_output_prep: ${nbfiles} elements treated\n"
  printf "=====> Um_output_prep: ${cnt} elements produced for XFER\n"
else
  printf "\n $flag_err ERRORS found in Um_output_prep.sh - ABORT\n"
fi

printf "\n=====>  Um_output_prep.sh ends: `date` ###########\n\n"

# End of task
. r.return.dot



