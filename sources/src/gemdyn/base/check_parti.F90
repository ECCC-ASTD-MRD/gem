!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

!**s/r check_parti - Checks domain partitioning
!
      logical function check_parti (F_si,F_n,F_npe)
!
      implicit none
#include <arch_specific.hf>
!
      integer F_si,F_n,F_npe
!
!author
!     Michel Desgagne
!
!revision
! v2_10 - Desgagne M.       - initial version
! v2_21 - Desgagne M.       - rpn_comm stooge for MPI
!
!object
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
!  F_si         I         Starting global index of local sub-domain
!  F_n          I         Global dimension
!  F_npe        I         Number of PEs
!----------------------------------------------------------------

!
      integer part_ok,part_ok_,err
!
!-------------------------------------------------------------------
!
      check_parti = .true.
!
      part_ok = 0
      if (F_si >= F_n) part_ok = 1
!
      if (part_ok == 1) write (6,900) F_n,F_npe,F_si
!
      call rpn_comm_ALLREDUCE (part_ok,part_ok_,1,"MPI_INTEGER", &
                                          "MPI_BOR","grid",err)
!
      if (part_ok_ > 0) check_parti = .false.
!

 900  format (/'  =====> ILLEGAL MPI PARTITIONNING: ',I5,' /',I5, &
               '   Starting index= ',i5)
!
!-------------------------------------------------------------------
!
      return
      end
