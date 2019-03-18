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

!**s/r set_ops8 - prepares tri-diagonal basic operators
!

!
      subroutine set_ops8 (F_oper_8, F_delt_8, F_wt_8, F_period_L,  &
                                                   NN, MXDM, F_case )
      implicit none
#include <arch_specific.hf>
!
      integer NN, MXDM, F_case
      real*8  F_oper_8(MXDM,3), F_delt_8(NN), F_wt_8
      logical F_period_L
!
!author
!     jean cote - 1990
!
!revision
! v2_00 - Desgagne/Lee      - initial MPI version (from setops v1_03)
!
!object
!     See above id.
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_oper_8     O    - operators
! F_delt_8     I    - distances between grid points
! F_wt_8       I    - weight (0.0,one,2.0,5.0)=>(explicit,pseudo,f.e,fourth)
! F_period_L   I    - .true. if periodic
! NN           I    - number of grid points
! F_case       I    - 1: second derivative
!                     2: identity projector
!                     3: first derivative
!                     4: first derivative with boundary condition modified
!                        by integration by part
!
!N.B.: periodic case => F_delt_8(NN) different from zero
!
!*
      real*8 zero, half, one, two
      parameter( zero = 0.0 )
      parameter( half = 0.5 )
      parameter( one  = 1.0 )
      parameter( two  = 2.0 )
!
      integer j
      real*8  pds, pdt
!
      if      ( F_case == 1 ) then
!               Second Derivative
         if ( F_period_L ) then
            F_oper_8(NN,3) = one/F_delt_8(NN)
            F_oper_8(1,1)   = F_oper_8(NN,3)
         else
            F_oper_8(NN,3) = zero
            F_oper_8(1,1)   = zero
         end if
         do j=1,NN-1
            F_oper_8(j,3)   = one/F_delt_8(j)
            F_oper_8(j+1,1) = F_oper_8(j,3)
            F_oper_8(j,2)   = - ( F_oper_8(j,1) + F_oper_8(j,3) )
         end do
            F_oper_8(NN,2) = - ( F_oper_8(NN,1) + F_oper_8(NN,3) )

      else if ( F_case == 2 ) then
!               Identity projector
         pds = one/( two + two * F_wt_8 )
         pdt = F_wt_8 * pds
         if ( F_wt_8 == zero ) pds = zero
         if ( F_wt_8 == zero ) pdt = half
         if ( F_period_L ) then
            F_oper_8(NN,3) = pds * F_delt_8(NN)
            F_oper_8(1,1)   = F_oper_8(NN,3)
            F_oper_8(1,2)   = pdt * ( F_delt_8(1) + F_delt_8(NN) )
         else
            F_oper_8(NN,3) = zero
            F_oper_8(1,1)   = zero
            F_oper_8(1,2)   = pdt * F_delt_8(1)
         end if
         do j=2,NN-1
            F_oper_8(j-1,3) = pds * F_delt_8(j-1)
            F_oper_8(j,2)   = pdt * ( F_delt_8(j-1) + F_delt_8(j) )
            F_oper_8(j,1)   = F_oper_8(j-1,3)
         end do
         j = NN
         if ( F_period_L ) then
            F_oper_8(j-1,3) = pds * F_delt_8(j-1)
            F_oper_8(j,2)   = pdt * ( F_delt_8(j-1) + F_delt_8(j) )
            F_oper_8(j,1)   = F_oper_8(j-1,3)
         else
            F_oper_8(j-1,3) = pds * F_delt_8(j-1)
            F_oper_8(j,2)   = pdt * F_delt_8(j-1)
            F_oper_8(j,1)   = F_oper_8(j-1,3)
         end if

      else if ( F_case == 3 ) then
!               First Derivative
         do j=1,NN-1
            F_oper_8(j,3)   =   half
            F_oper_8(j,2)   =   zero
            F_oper_8(j+1,1) = - half
         end do
         if ( F_period_L ) then
            F_oper_8(NN,3) =   half
            F_oper_8(1,1)   = - half
            F_oper_8(1,2)   =   zero
            F_oper_8(NN,2) =   zero
         else
            F_oper_8(NN,3) =   zero
            F_oper_8(1,1)   =   zero
            F_oper_8(1,2)   = - half
            F_oper_8(NN,2) =   half
         end if

      else if ( F_case == 4 ) then
!               First Derivative with B.C. modified by integration by part
         do j=1,NN-1
            F_oper_8(j,3)   = - half
            F_oper_8(j,2)   =   zero
            F_oper_8(j+1,1) =   half
         end do
         if ( F_period_L ) then
            F_oper_8(NN,3) = - half
            F_oper_8(1,1)   =   half
            F_oper_8(1,2)   =   zero
            F_oper_8(NN,2) =   zero
         else
            F_oper_8(NN,3) =   zero
            F_oper_8(1,1)   =   zero
            F_oper_8(1,2)   = - half
            F_oper_8(NN,2) =   half
         end if
      else if ( F_case == 5 ) then
!               Average operator
         F_oper_8(NN,3) = zero
         F_oper_8(1,1)   = zero
         do j=1,NN-1
            F_oper_8(j,3)   = one/F_delt_8(j)
            F_oper_8(j+1,1) = F_oper_8(j,3)
            F_oper_8(j,2)   = F_oper_8(j,1) + F_oper_8(j,3)
         end do
          F_oper_8(NN,2) = F_oper_8(NN,1) + F_oper_8(NN,3)

      end if
      return
      end
