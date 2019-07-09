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

!**s/r Boundary- 3D_elliptic boundary RHS computation for GEM_H
!

      subroutine  boundary ( F_rhs,F_rt,F_rf,F_nt,Minx,Maxx,Miny,Maxy, &
                             Nk,ni,nj,i0,j0,in,jn )
      use geomh
      use gem_options
      use HORgrid_options
      use tdpack
      use glb_ld
      use cstv
      use ver
      use sol
      use opr
      use metric
      use ptopo

      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
!
      integer, intent(in) :: Minx,Maxx,Miny,Maxy,Nk,ni,nj
      integer, intent(in) :: i0,j0,in,jn
      real(kind=REAL64), dimension(ni,nj,Nk),             intent(out) :: F_rhs
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(in)  :: F_rt, F_nt, F_rf
!author
!       Abdessamad Qaddouri -  2017
!
      integer :: i,j
      real(kind=REAL64)  :: add_v8,bdd_v8,cdd_v8
      real(kind=REAL64), parameter :: one=1.d0, half=0.5d0, zero=0.d0
      real(kind=REAL64), dimension(l_minx:l_maxx,l_miny:l_maxy) :: Afdg1,Bfdg1,Cfdg1
      real, dimension(l_minx:l_maxx,l_miny:l_maxy) :: fdg2
      integer sol_pil_w_ext,sol_pil_e_ext,sol_pil_s_ext,sol_pil_n_ext
!
!     ---------------------------------------------------------------
!
      fdg2  = 0.0
      Afdg1 = 0.0d0
      Bfdg1 = 0.0d0
      Cfdg1 = 0.0d0
      sol_pil_s_ext=sol_pil_s-1
      sol_pil_n_ext= sol_pil_n
      sol_pil_w_ext=sol_pil_w-1
      sol_pil_e_ext=sol_pil_e

      if(.not.Grd_yinyang_L) then
         if (l_west) sol_pil_w_ext=sol_pil_w
         if (l_east) sol_pil_e_ext=sol_pil_e+1
         if (l_south) sol_pil_s_ext=sol_pil_s
         if (l_north) sol_pil_n_ext= sol_pil_n+1
      endif


      do j=j0,jn
         do i=i0,in
            fdg2(i,j) = mc_css_H_8(i,j)   * (F_rt(i,j,l_nk)-Ver_wmstar_8(G_nk)*F_rt(i,j,l_nk-1) &
                     + Cstv_invT_m_8*(F_rf(i,j,l_nk)-Ver_wmstar_8(G_nk)*F_rf(i,j,l_nk-1))) &
                      - mc_css_H_8(i,j)   * (F_nt(i,j,l_nk ) - Ver_wmstar_8(G_nk)*F_nt(i,j,l_nk-1))
         end do
      end do

      if ((Grd_yinyang_L)) then
         call yyg_xchng (fdg2, l_minx,l_maxx,l_miny,l_maxy, &
                               l_ni,l_nj, 1, .false., 'CUBIC', .true.)
      else
         call rpn_comm_xch_halo(fdg2,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,1, &
                             G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      endif

      do j=j0,jn
         do i=i0,in
            Cfdg1(i,j) =  gama_8*(fdg2(i,j)*mc_iJz(i  ,j,Nk )&
                   -mu_8*half*fdg2(i,j))
         end do
      end do

      do j=1+sol_pil_s_ext, l_nj-sol_pil_n
         do i=1+sol_pil_w_ext, l_ni-sol_pil_e_ext
            Afdg1(i,j) = - mc_Jx(i,j,Nk) * &
                          Ver_wp_8%m(Nk)*half*( fdg2(i+1,j)*mc_iJz(i+1,j,Nk) &
                                              + fdg2(i  ,j)*mc_iJz(i  ,j,Nk) )
         end do
      end do

      do j=1+sol_pil_s_ext, l_nj-sol_pil_n_ext
         do i=1+sol_pil_w_ext, l_ni-sol_pil_e
            Bfdg1(i,j) = - mc_Jy(i,j,Nk) * &
                        Ver_wp_8%m(Nk)*half*( fdg2(i,j+1)*mc_iJz(i,j+1,Nk) &
                                            + fdg2(i,j  )*mc_iJz(i,j  ,Nk) )
         end do
      end do

      do j=j0,jn
         do i=i0,in
            add_v8 = (Afdg1 (i,j)-Afdg1 (i-1,j))*geomh_invDXM_8(j) &
                     + half * (mc_Ix(i,j,Nk)*(Afdg1(i,j)+Afdg1(i-1,j)))
            bdd_v8 = (Bfdg1 (i,j)*geomh_cyM_8(j)-Bfdg1 (i,j-1)*&
                     geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                     + half * (mc_Iy(i,j,Nk)*(Bfdg1(i,j)+Bfdg1(i,j-1)))
            cdd_v8=Cfdg1(i,j)*Ver_idz_8%m(Nk)&
                     +(mc_Iz(i,j,Nk)-epsi_8)*(Ver_wp_8%m(Nk)*Cfdg1(i,j))
            F_rhs(i,j,Nk)=F_rhs(i,j,Nk)-Cstv_hco0_8*(add_v8+bdd_v8+cdd_v8)

         end do
      end do

      return
      end

