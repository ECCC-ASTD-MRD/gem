#!/bin/sh

cd $1
mylist=$(ls -1 .)
cd $2

for item in $mylist ; do
    [ -d ../configurations/$item  -a ! -e $item ] && ln -sf ../configurations/$item $item
done
