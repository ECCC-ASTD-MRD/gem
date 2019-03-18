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

!**s/r set_sol - Computes matrices a,b,c for the elliptic solver

      subroutine set_sol
      use lam_options
      use dyn_fisl_options
      use glb_ld
      use cstv
      use lun
      use sol
      use opr
      use prec
      use trp
      implicit none
#include <arch_specific.hf>

      integer k,k1,k0
      real*8  wk(G_nk)
      real*8 yg_8(G_nj)
!     __________________________________________________________________
!
      do k =1, G_nk-Lam_gbpil_T
         do k1=1, G_nk-Lam_gbpil_T
            wk(k) = Opr_zevec_8 ((k-1)*G_nk+k1)
         end do
         wk(k)= Cstv_hco0_8*Opr_zeval_8(k)
      end do

      sol_nk = trp_12sn-east*Lam_gbpil_T

      if ( Sol_type_S(1:9) == 'ITERATIVE' ) then

         if (Sol_type_S(11:12) == '2D') then
            if (Lun_out > 0) write (Lun_out,1001) trim(sol2D_precond_S)
         else
            if (Lun_out > 0) write (Lun_out,1002) trim(Sol3D_krylov_S), trim(sol3D_precond_S)
         end if

         sol_pil_w= pil_w ; sol_pil_e= pil_e
         sol_pil_s= pil_s ; sol_pil_n= pil_n
         sol_niloc= (l_ni-pil_e)-(1+pil_w)+1
         sol_njloc= (l_nj-pil_n)-(1+pil_s)+1
         sol_nloc = sol_niloc*sol_njloc*Schm_nith

         allocate (Prec_xevec_8(sol_niloc*sol_niloc)          ,&
                   Prec_xeval_8(sol_niloc),Prec_ai_8(sol_nloc),&
                   Prec_bi_8(sol_nloc),Prec_ci_8(sol_nloc))

         do k =1,G_nk
            do k0=1,G_nk
               wk(k) = Opr_zevec_8 ((k-1)*G_nk+k0)
            end do
            if ( k <= Schm_nith ) then
               wk(k) = (Cstv_hco1_8+Cstv_hco0_8*Opr_zeval_8(k))
            end if
         end do

         call eigenabc_local (Prec_xeval_8,Prec_xevec_8,Prec_ai_8,&
                              Prec_bi_8,Prec_ci_8,l_ni,l_nj      ,&
                        sol_niloc,sol_njloc,Schm_nith,l_i0,l_j0,wk)
      else
         yg_8(1:G_nj)=G_yg_8(1:G_nj)
         call sol_abc ( wk,yg_8,Opr_opsyp0_8, &
                        Opr_opsyp2_8,Opr_xeval_8 , &
             trp_12smin, trp_12smax,  sol_nk, trp_12sn0, &
             trp_22min , trp_22max , trp_22n, trp_22n0 , &
             G_ni,G_nj,G_nk, Sol_ai_8, Sol_bi_8, Sol_ci_8 )
      end if

 1001 format(/,'WILL USE FGMRES 2D ITERATIVE SOLVER WITH ',a, &
               ' PRECONDITIONNER' &
             /,'=============================================')
 1002 format(/,'WILL USE ',a,' 3D ITERATIVE SOLVER WITH ',a, &
               ' PRECONDITIONNER' &
             /,'=============================================')
!
!     ---------------------------------------------------------------
!
      return
      end
