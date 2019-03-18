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

!*s/r prever - prepares projection matrix for the vertical
!       exclusively for staggered model symmetric and  nonsymmetric_ln(Z) version.
!       no special treatement for singularity. Matrices are non-singulars by construction.
!
      subroutine preverln2 ( F_eval_8, F_levec_8, F_evec_8, &
                             F_nk, KDIM, F_errcode )
      use cstv
      use glb_ld
      use lun
      use opr
      implicit none
#include <arch_specific.hf>

      integer F_nk, KDIM, F_errcode
      real*8 F_eval_8(KDIM), F_levec_8(KDIM,KDIM), F_evec_8(KDIM,KDIM)
!
!author  Abdessamad Qaddouri - 2007
!
!revision
! v4_00 - Plante & Girard   - Log-hydro-pressure coord on Charney-Phillips grid
! v4.60 - Qaddouri A.       - prevents F_evec_8 deviation from "to be negative definite"
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_eval_8     O    - eigenvalues
! F_evec_8     O    - Right eigenvector matrix
!F_levec_8     O    - left eigenvector matrix
! F_wk1_8           - work field
! F_nk         I    - number of vertical levels
!


      integer i, j, err
      real*8 zero, one, xxx, yyy
      parameter(zero=.0d0, one=1.d0)
      real*8 wk1(KDIM,KDIM), B1(KDIM,KDIM)
!
! --------------------------------------------------------------------
!
      F_evec_8= 0. ; F_levec_8= 0. ; wk1= 0. ; B1= 0. ; err= 0

      xxx = - Cstv_hco2_8
      yyy = - Cstv_hco1_8

      do j=1,F_nk

         i = j - 1
         if ( i > 0 ) then
!           A wing
            F_evec_8(i,j) =     Opr_opszp2_8(2*G_nk+i) &
                          +     Opr_opszpl_8(2*G_nk+i) &
                          + xxx*Opr_opszpm_8(2*G_nk+i)
            if (F_evec_8(i,j)< 0.0) err= err-1
         end if

         i = j
!        B: positive definit
         wk1(i,j) = Opr_opszp0_8(G_nk+i)
!        A diag
         F_evec_8(i,j) =     Opr_opszp2_8(G_nk+i) &
                       +     Opr_opszpl_8(G_nk+i) &
                       + xxx*Opr_opszpm_8(G_nk+i) &
                       + yyy*Opr_opszp0_8(G_nk+i)
         i=j+1
         if(i < (F_nk+1)) then
            F_evec_8(i,j) =     Opr_opszp2_8(i) &
                          +     Opr_opszpl_8(i) &
                          + xxx*Opr_opszpm_8(i)
            if (F_evec_8(i,j)< 0.0) err= err-1
         end if

      end do

! note:  B1 is modified by nsyeigl
      do j=1,F_nk
      do i=1,F_nk
         F_levec_8(i,j)= F_evec_8(i,j)
         B1(i,j)= wk1(i,j)
      end do
      end do

      F_errcode= 0
      if ((err < 0).and.(lun_out > 0)) then
         write(lun_out,9000)
         F_errcode = -1
      else
         call nsyeigl (F_eval_8,F_levec_8,F_evec_8,B1,F_nk,KDIM,8*F_nk)
      end if

 9000 format (/x,42('*')/2x,'VERTICAL LAYERING (vertical resolution)',&
              /2x,'IS INCOMPATIBLE WITH THE TIMESTEP.',&
              /2x,'THE SOLVER WILL NOT WORK --- ABORTING',/x,42('*')/)
!
! --------------------------------------------------------------------
!
      return
      end
