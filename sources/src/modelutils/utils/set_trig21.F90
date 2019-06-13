!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------


subroutine set_trig21( F_ai_8, F_bi_8, F_ci_8, F_di_8, &
     F_a_8, F_b_8, F_c_8, &
     F_idel, F_jdel, fni, fnj, F_per_L)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
   include "rmnlib_basics.inc"
   !@object Preparation pour un eventuel solveur tridiagonal
   !        Le solveur est soit parallel ou non
   !        Periodicite est envisagee aussi.
   !@arguments
   logical F_per_L
   integer F_idel, F_jdel, fni, fnj
   real(REAL64) :: F_ai_8(F_idel,fni),F_bi_8(F_idel,fni),F_ci_8(F_idel,fni), &
        F_di_8(F_idel,fni), F_a_8(F_idel,fni), F_b_8(F_idel,fni), &
        F_c_8 (F_idel,fni),tempo
   !@author Abdessamad Qaddouri
   !@revision
   ! v2_00 - Qaddouri A.       - initial MPI version
   ! v3_11 - Gravel S.         - modify for theoretical cases

   integer i,l,m
   real(REAL64), parameter :: one = 1.0
   !     ---------------------------------------------------------------
   do l= 1,fnj*F_jdel,F_jdel
      F_bi_8(l,1) = one/F_b_8(l,1)
      F_ci_8(l,1) = F_c_8(l,1) * F_bi_8(l,1)
   enddo

   if (F_per_L) then
      m =fni-1
   else
      m=fni
   endif

   if ( m .gt. 1) then    ! to avoid calculation with theoretical cases
      do i=2, m-1
         do l= 1,fnj*F_jdel,F_jdel
            F_bi_8(l,i) = one/( F_b_8(l,i) - F_a_8(l,i) * F_ci_8(l,i-1) )
            F_ci_8(l,i) = F_bi_8(l,i) * F_c_8(l,i)
            F_ai_8(l,i) = F_bi_8(l,i) * F_a_8(l,i)
         enddo
      enddo

      do l= 1,fnj*F_jdel,F_jdel
         tempo=F_b_8(l,m)-F_a_8(l,m) *F_ci_8(l,m-1)
         if(F_b_8(l,m)-F_a_8(l,m) *F_ci_8(l,m-1).eq.0.)then
            tempo=epsilon(tempo)
         endif
         F_bi_8(l,m)=one/tempo
         F_ai_8(l,m)=F_bi_8(l,m)*F_a_8(l,m)
      enddo
   endif                  ! end of exception for theoretical cases

   if (F_per_L) then

      do l= 1,fnj*F_jdel,F_jdel
         F_di_8(l,1) = -F_bi_8(l,1) * F_a_8(l,1)
      enddo

      do i=2,m
         do l= 1,fnj*F_jdel,F_jdel
            F_di_8(l,i) = -F_ai_8(l,i) * F_di_8(l,i-1)
         enddo
      enddo

      do l= 1,fnj*F_jdel,F_jdel
         F_di_8(l,m)= F_di_8(l,m) - F_bi_8(l,m)*F_c_8(l,m)
      enddo

      do i=m-1,1,-1
         do l= 1,fnj*F_jdel,F_jdel
            F_di_8(l,i)= F_di_8(l,i) - F_ci_8(l,i)*F_di_8(l,i+1)
         enddo
      enddo

      do l= 1,fnj*F_jdel,F_jdel

         F_bi_8(l,fni)= one/(F_b_8(l,fni)+F_c_8(l,fni)*F_di_8(l,1) &
              + F_a_8(l,fni)* F_di_8(l,m) )
         F_ai_8(l,1)  = - F_a_8(l,fni)*F_bi_8(l,fni)
         F_ci_8(l,fni)= - F_c_8(l,fni)* F_bi_8(l,fni)
      enddo

   endif
   !     ---------------------------------------------------------------
   return
end subroutine set_trig21

