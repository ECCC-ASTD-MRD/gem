#!/bin/bash
mkdir -p local1 local2
echo "'model'" >local1/id.txt
echo "'assim'" >local2/id.txt
echo "'wrong'" >id.txt
export RPN_COMM_DOM=2,0,1,1,2,1,5 
export RPN_COMM_DIRS="' ','./local1','./local2'"  
mpirun -n 2 ../Build/${EC_ARCH}/TEST_004.Abs : -n 3 ../Build/${EC_ARCH}/TEST_004.Abs
rm -f id.txt local1/id.txt local2/id.txt
