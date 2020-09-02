#!/bin/bash

export MPI_NPEX=2
export MPI_NPEY=1

export MPI_NBLOCX=${MPI_NPEX}
export MPI_NBLOCY=${MPI_NPEY}
echo "P=${MPI_NPEX}x${MPI_NPEY} B=${MPI_NBLOCX}x${MPI_NBLOCY}"
r.mpirun -pgm malib${EC_ARCH}/${0%.*}.Abs -npex $MPI_NPEX -npey $MPI_NPEY

export MPI_NBLOCX=1
export MPI_NBLOCY=1
echo "P=${MPI_NPEX}x${MPI_NPEY} B=${MPI_NBLOCX}x${MPI_NBLOCY}"
r.mpirun -pgm malib${EC_ARCH}/${0%.*}.Abs -npex $MPI_NPEX -npey $MPI_NPEY

export MPI_NPEX=1
export MPI_NPEY=2

export MPI_NBLOCX=${MPI_NPEX}
export MPI_NBLOCY=${MPI_NPEY}
echo "P=${MPI_NPEX}x${MPI_NPEY} B=${MPI_NBLOCX}x${MPI_NBLOCY}"
r.mpirun -pgm malib${EC_ARCH}/${0%.*}.Abs -npex $MPI_NPEX -npey $MPI_NPEY

export MPI_NBLOCX=1
export MPI_NBLOCY=1
echo "P=${MPI_NPEX}x${MPI_NPEY} B=${MPI_NBLOCX}x${MPI_NBLOCY}"
r.mpirun -pgm malib${EC_ARCH}/${0%.*}.Abs -npex $MPI_NPEX -npey $MPI_NPEY
