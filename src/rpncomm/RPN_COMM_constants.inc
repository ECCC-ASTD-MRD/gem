!/* RPN_COMM - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2012  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
! error status
        integer, parameter :: RPN_COMM_OK    =  0
        integer, parameter :: RPN_COMM_ERROR = -1
!
        integer, parameter :: RPN_COMM_MAGIC = int( Z'13579BDF', kind=kind(RPN_COMM_MAGIC))
! command op codes
        integer, parameter :: RPN_COMM_MASTER = 0
        integer, parameter :: RPN_COMM_GET=12345
        integer, parameter :: RPN_COMM_SET=12346
        integer, parameter :: RPN_COMM_DEL=12347
! constants for the transpose operations along X/Y
        integer, parameter :: RPN_COMM_FORWARD_X=+1
        integer, parameter :: RPN_COMM_BACKWARD_X=-1
        integer, parameter :: RPN_COMM_FORWARD_y=+2
        integer, parameter :: RPN_COMM_BACKWARD_Y=-2
!
        logical,parameter :: RPN_COMM_TOPO_X      = .true.
        logical,parameter :: RPN_COMM_TOPO_Y      = .false.
        logical,parameter :: RPN_COMM_TOPO_FILL   = .true.
        logical,parameter :: RPN_COMM_TOPO_NOFILL = .false.
! communicators
        character(len=*),parameter :: RPN_COMM_WORLD      = 'WORLD'      ! PLS do not use
        character(len=*),parameter :: RPN_COMM_ALLGRIDS   = 'ALLGRIDS'   ! use rather than RPN_COMM_WORLD
        character(len=*),parameter :: RPN_COMM_ALL        = 'ALL'
        character(len=*),parameter :: RPN_COMM_MULTIGRID  = 'MULTIGRID'
        character(len=*),parameter :: RPN_COMM_GRIDPEERS  = 'GRIDPEERS'
        character(len=*),parameter :: RPN_COMM_ALLDOMAINS = 'ALLDOMAINS'
        character(len=*),parameter :: RPN_COMM_GRID       = 'GRID'       ! use rather than RPN_COMM_DOMM
        character(len=*),parameter :: RPN_COMM_DOMM       = 'DOMM'
        character(len=*),parameter :: RPN_COMM_DEFAULT    = 'DEFO'
        character(len=*),parameter :: RPN_COMM_BLOCMASTER = 'BLOCMASTER'
        character(len=*),parameter :: RPN_COMM_BLOCK      = 'BLOC'
        character(len=*),parameter :: RPN_COMM_BLOC_COMM  =  RPN_COMM_BLOCK 
        character(len=*),parameter :: RPN_COMM_EW         = 'EW'
        character(len=*),parameter :: RPN_COMM_ROW        = 'EW'
        character(len=*),parameter :: RPN_COMM_NS         = 'NS'
        character(len=*),parameter :: RPN_COMM_COLUMN     = 'NS'
        character(len=*),parameter :: RPN_COMM_UNIVERSE   = 'UNIVERSE'
        character(len=*),parameter :: RPN_COMM_NULL       = 'RPN_COMM_NULL'
!        character(len=*),parameter :: RPN_COMM_PE0        = 'PE_00'
! miscellaneous NULL/ANY constants
        character (len=32), parameter :: RPN_COMM_ANY_SOURCE =           &                                   
     &                                  'MPI_ANY_SOURCE'
        character (len=32), parameter :: RPN_COMM_ANY_TAG =              &                                   
     &                                  'MPI_ANY_TAG'
        character (len=32), parameter :: RPN_COMM_DATATYPE_NULL =        &                                   
     &                                  'MPI_DATATYPE_NULL'
        character (len=32), parameter :: RPN_COMM_COMM_NULL =            &                                   
     &                                  'MPI_COMM_NULL'
        character (len=32), parameter :: RPN_COMM_GROUP_NULL =           &                                   
     &                                  'MPI_GROUP_NULL'
        character (len=32), parameter :: RPN_COMM_COMM_WORLD =           &                                   
     &                                  'MPI_COMM_WORLD'
        character (len=32), parameter :: RPN_COMM_SUCCESS =              &                                   
     &                                  'MPI_SUCCESS'
! data types
        character (len=32), parameter :: RPN_COMM_BYTE =                 &                                   
     &                                  'MPI_BYTE'
        character (len=32), parameter :: RPN_COMM_PACKED =               &                                   
     &                                  'MPI_PACKED'
        character (len=32), parameter :: RPN_COMM_UB =                   &                                   
     &                                  'MPI_UB'
        character (len=32), parameter :: RPN_COMM_LB =                   &                                   
     &                                  'MPI_LB'
        character (len=32), parameter :: RPN_COMM_CHARACTER =            &                                   
     &                                  'MPI_CHARACTER'
        character (len=32), parameter :: RPN_COMM_LOGICAL =              &                                   
     &                                  'MPI_LOGICAL'
        character (len=32), parameter :: RPN_COMM_INTEGER =              &                                   
     &                                  'MPI_INTEGER'
        character (len=32), parameter :: RPN_COMM_INTEGER1 =             &                                   
     &                                  'MPI_INTEGER1'
        character (len=32), parameter :: RPN_COMM_INTEGER2 =             &                                   
     &                                  'MPI_INTEGER2'
        character (len=32), parameter :: RPN_COMM_INTEGER4 =             &                                   
     &                                  'MPI_INTEGER4'
        character (len=32), parameter :: RPN_COMM_INTEGER8 =             &                                   
     &                                  'MPI_INTEGER8'
        character (len=32), parameter :: RPN_COMM_INTEGER16 =            &                                   
     &                                  'MPI_INTEGER16'
        character (len=32), parameter :: RPN_COMM_REAL =                 &                                   
     &                                  'MPI_REAL'
        character (len=32), parameter :: RPN_COMM_REAL4 =                &                                   
     &                                  'MPI_REAL4'
        character (len=32), parameter :: RPN_COMM_REAL8 =                &                                   
     &                                  'MPI_REAL8'
        character (len=32), parameter :: RPN_COMM_REAL16 =               &                                   
     &                                  'MPI_REAL16'
        character (len=32), parameter :: RPN_COMM_DOUBLE_PRECISION =     &                                   
     &                                  'MPI_DOUBLE_PRECISION'
        character (len=32), parameter :: RPN_COMM_COMPLEX =              &                                   
     &                                  'MPI_COMPLEX'
        character (len=32), parameter :: RPN_COMM_COMPLEX8 =             &                                   
     &                                  'MPI_COMPLEX8'
        character (len=32), parameter :: RPN_COMM_COMPLEX16 =            &                                   
     &                                  'MPI_COMPLEX16'
        character (len=32), parameter :: RPN_COMM_COMPLEX32 =            &                                   
     &                                  'MPI_COMPLEX32'
        character (len=32), parameter :: RPN_COMM_DOUBLE_COMPLEX =       &                                   
     &                                  'MPI_DOUBLE_COMPLEX'
        character (len=32), parameter :: RPN_COMM_2REAL =                &                                   
     &                                  'MPI_2REAL'
        character (len=32), parameter :: RPN_COMM_2DOUBLE_PRECISION =    &                                   
     &                                  'MPI_2DOUBLE_PRECISION'
        character (len=32), parameter :: RPN_COMM_2INTEGER =             &                                   
     &                                  'MPI_2INTEGER'
        character (len=32), parameter :: RPN_COMM_2COMPLEX =             &                                   
     &                                  'MPI_2COMPLEX'
        character (len=32), parameter :: RPN_COMM_2DOUBLE_COMPLEX =      &                                   
     &                                  'MPI_2DOUBLE_COMPLEX'
! operations
        character (len=32), parameter :: RPN_COMM_MAX =                  &                                   
     &                                  'MPI_MAX'
        character (len=32), parameter :: RPN_COMM_MIN =                  &                                   
     &                                  'MPI_MIN'
        character (len=32), parameter :: RPN_COMM_SUM =                  &                                   
     &                                  'MPI_SUM'
        character (len=32), parameter :: RPN_COMM_PROD =                 &                                   
     &                                  'MPI_PROD'
        character (len=32), parameter :: RPN_COMM_LAND =                 &                                   
     &                                  'MPI_LAND'
        character (len=32), parameter :: RPN_COMM_BAND =                 &                                   
     &                                  'MPI_BAND'
        character (len=32), parameter :: RPN_COMM_LOR =                  &                                   
     &                                  'MPI_LOR'
        character (len=32), parameter :: RPN_COMM_BOR =                  &                                   
     &                                  'MPI_BOR'
        character (len=32), parameter :: RPN_COMM_LXOR =                 &                                   
     &                                  'MPI_LXOR'
        character (len=32), parameter :: RPN_COMM_BXOR =                 &                                   
     &                                  'MPI_BXOR'
        character (len=32), parameter :: RPN_COMM_MAXLOC =               &                                   
     &                                  'MPI_MAXLOC'
        character (len=32), parameter :: RPN_COMM_MINLOC =               &                                   
     &                                  'MPI_MINLOC'
        character (len=32), parameter :: RPN_COMM_REPLACE =              &                                   
     &                                  'MPI_REPLACE'
        character (len=32), parameter :: RPN_COMM_OP_NULL =              &                                   
     &                                  'MPI_OP_NULL'
