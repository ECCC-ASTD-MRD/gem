#!/bin/bash

lis=${PWD}/md5sum$$.lis

if [ -d $1 ] ; then

  cd $1
  find -L ./ -type f -exec md5sum {} \; > ${lis}.d

else

  md5sum $1 > ${lis}.f

fi

if [ -e ${lis}.d ] ; then
  while read line ; do
    f=$(echo $line | awk '{print $2}')
    s=$(echo $line | awk '{print $1}')
    echo $f $s
  done < ${lis}.d >> ${lis}_final
  rm ${lis}.d
fi

if [ -e ${lis}.f ] ; then
  while read line ; do
    f=$(echo $line | awk '{print $2}')
    s=$(echo $line | awk '{print $1}')
    if [ -n "${2}" ] ; then
      echo $2 $s
    else
      echo $(basename $f) $s
    fi
  done < ${lis}.f >> ${lis}_final
  rm ${lis}.f
fi

cat ${lis}_final | sort

rm ${lis}_final
