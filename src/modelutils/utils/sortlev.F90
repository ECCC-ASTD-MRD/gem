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

!opyright (C) 2001  MSC-RPN COMM  %%%MC2%%%
!**s/r e_sortlev - to eliminate redundant levels and sort in increasing order
!
      subroutine sortlev (F_rna,F_ip1,LV,LVO)

      implicit none
!!!#include <arch_specific.hf>
!
      integer LV,F_ip1(LV),LVO
      real F_rna(LV)
!
!REVISION
! v4_03 - Lee V.            - Adapted for GEMDM
!
!ARGUMENTS
!    NAMES       I/O  TYPE  DESCRIPTION
!    F_rna       level value 
!    F_ip1       ip1code
!    LV          number of incoming levels
!    LVO         number of outgoing levels
!
!
      integer k,i,m,j,n
      real x1
!*
!----------------------------------------------------------------------
!
      n = lv
      do i = 1, n-1
         k = i
         do j = i+1, n
            if (F_rna(k) .gt. F_rna(j))  k=j
         enddo
         if (k .ne. i) then
            x1     = F_rna(k)
            m      = F_ip1(k)
            F_rna(k) = F_rna(i)
            F_ip1(k)  = F_ip1(i)
            F_rna(i) = x1
            F_ip1(i)  = m
         endif
      enddo

!     eliminate levels that are redundant in LISTE
      i = 1
      do j=2,n
         if (F_rna(i) .ne. F_rna(j)) then
             i = i+1
             if (i .ne. j) then
                 F_rna(i) = F_rna(j)
                 F_ip1(i) =  F_ip1(j)
             endif
         endif
      enddo
      lvo = i

!
!----------------------------------------------------------------------
      return
      end

