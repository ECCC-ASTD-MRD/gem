#!/bin/bash

export MPI_NPEX=1
export MPI_NPEY=1
export MPI_NBLOCX=1
export MPI_NBLOCY=1

r.mpirun -pgm malib${EC_ARCH}/${0%.*}.Abs -npex $MPI_NPEX -npey $MPI_NPEY

export MPI_NPEX=1
export MPI_NPEY=3
export MPI_NBLOCX=1
export MPI_NBLOCY=1

r.mpirun -pgm malib${EC_ARCH}/${0%.*}.Abs -npex $MPI_NPEX -npey $MPI_NPEY

export MPI_NPEX=1
export MPI_NPEY=3
export MPI_NBLOCX=1
export MPI_NBLOCY=3

r.mpirun -pgm malib${EC_ARCH}/${0%.*}.Abs -npex $MPI_NPEX -npey $MPI_NPEY
