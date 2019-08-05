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
!
!**s/r pre_H -  Combine some rhs obtaining Rt' and Rc", the linear
!               contributions to the rhs of Helmholtz equation
!            -  Height-type vertical coordinate

!
      subroutine fislh_pre ( F_ru ,F_rv ,F_rt ,F_rw ,F_rc, F_rf, F_rb, F_nest_t, F_fis, &
                            Minx, Maxx, Miny, Maxy, i0, j0, k0, in, jn, Nk )
      use HORgrid_options
      use gem_options
      use geomh
      use tdpack
      use glb_ld
      use lun
      use cstv
      use ver
      use metric
      use fislh_sol
      use dyn_fisl_options

      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Minx,Maxx,Miny,Maxy, i0, j0, k0, in, jn, Nk
      real, dimension(Minx:Maxx,Miny:Maxy),     intent(out)   :: F_rb
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(in)    :: F_nest_t
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(inout) :: F_ru,F_rv,F_rt
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(inout) :: F_rw,F_rc,F_rf
      real, dimension(Minx:Maxx,Miny:Maxy),     intent(in)    :: F_fis
!
!Author: Claude Girard, July 2017 (initial version)
!        Syed Husain, June 2019 (major revision)
!        Abdessamad Qaddouri, July 2019 (open-top)
!
      integer :: i, j, k, km,k0t
      real    :: w0
      real(kind=REAL64)  :: div, c1, w1, w2
      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0, half=0.5d0
!
!     ---------------------------------------------------------------
!
      if (Lun_debug_L) write (Lun_out,1000)
      k0t = k0
      if(Schm_opentop_L) k0t= k0-1


      c1 = grav_8 * Cstv_tau_8

      call rpn_comm_xch_halo( F_ru, l_minx,l_maxx,l_miny,l_maxy, l_niu, l_nj,G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      call rpn_comm_xch_halo( F_rv, l_minx,l_maxx,l_miny,l_maxy, l_ni, l_njv,G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )

!
      if (Schm_opentop_L) then
!$omp do
         do j= j0, jn
            do i= i0, in
               F_rb(i,j) = F_rt(i,j,k0t)-mu_8*Cstv_tau_nh_8* F_rw(i,j,k0t)

            end do
         end do
!$omp enddo
      end if

!*************************************
! Combination of governing equations *
!*************************************

!$omp parallel private(km,div,w1,w2)

!$omp do
      do k=k0, l_nk

         do j= j0, jn
            do i= i0, in

               F_rc(i,j,k) = F_rc(i,j,k) - Cstv_invT_8 * mc_logJz(i,j,k)

!              Compute Rc' combining Ru, Rv and Rc
!              ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               div  = (F_ru(i,j,k)-F_ru(i-1,j,k))*geomh_invDXM_8(j) &
                    + (F_rv(i,j,k)*geomh_cyM_8(j)-F_rv(i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                    + isol_i*half*(mc_Ix(i,j,k)*(F_ru(i,j,k)+F_ru(i-1,j,k)) &
                               +mc_Iy(i,j,k)*(F_rv(i,j,k)+F_rv(i,j-1,k)) )

               F_rc(i,j,k) = div - Cstv_invT_m_8 * F_rc(i,j,k) - Cstv_bar0_8 * F_fis(i,j)

            end do
         end do

      end do
!$omp end do

!$omp do
      do k=k0t, l_nk
         do j= j0, jn
            do i= i0, in
!              Compute Rf'
!              ~~~~~~~~~~~
               w0 = Cstv_invT_nh_8*(ztht(i,j,k)-Ver_z_8%t(k))*Cstv_bar1_8
               F_rf(i,j,k) = F_rf(i,j,k) - w0

!              Compute Rt'
!              ~~~~~~~~~~~
               F_rt(i,j,k) = gama_8 * ( c1 * F_rt(i,j,k) + F_rw(i,j,k) )

            end do
         end do

      end do
!$omp end do


!$omp do
      do k=k0,l_nk
         km=max(k-1,1)
         do j= j0, jn
            do i= i0, in
               w1= (Ver_idz_8%m(k) + (isol_i*mc_Iz(i,j,k) - epsi_8)*Ver_wp_8%m(k))
               w2= (Ver_idz_8%m(k) - (isol_i*mc_Iz(i,j,k) - epsi_8)*Ver_wm_8%m(k))*Ver_onezero(k)

!              Compute Rc"
!              ~~~~~~~~~~~
               F_rc(i,j,k) = F_rc(i,j,k) + Cstv_invT_m_8* &
                            (F_rf(i,j,k)-Ver_onezero(k)*F_rf(i,j,km))*Ver_idz_8%m(k)

!              Compute Rc"', combining Rc" and Rt"
!              ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               F_rc(i,j,k) = F_rc(i,j,k) + w1 * F_rt(i,j,k) - w2 * F_rt(i,j,km)

            end do
         end do
      end do
!$omp end do
      if (Schm_opentop_L) then
!        Apply opentop boundary conditions

         w1=Cstv_invT_8/Ver_Tstar_8%t(k0t)
!$omp do
         do j= j0, jn
            do i= i0, in
               F_rb(i,j) = F_rb(i,j) - (one/Cstv_tau_8+mu_8*Cstv_tau_nh_8*grav_8)*(F_nest_t(i,j,k0t)&
                                                  -Ver_Tstar_8%t(k0t))/Ver_Tstar_8%t(k0t)
               F_rc(i,j,k0) = F_rc(i,j,k0  ) + Ver_cstp_8 * F_rb(i,j)
            end do
         end do
!$omp enddo

      end if

!     Apply lower boundary conditions
!
!$omp do
      do j= j0, jn
         do i= i0, in
            F_rc(i,j,l_nk) = F_rc(i,j,l_nk) - isol_d*mc_cssp_H_8(i,j) * &
                     (F_rt(i,j,l_nk)-Ver_wmstar_8(G_nk)*F_rt(i,j,l_nk-1) &
                    + Cstv_invT_m_8*(F_rf(i,j,l_nk)-Ver_wmstar_8(G_nk)*F_rf(i,j,l_nk-1)))
         end do
      end do
!$omp end do

!$omp end parallel

1000  format(3X,'UPDATE  THE RIGHT-HAND-SIDES: (S/R PRE_H)')
!
!     ---------------------------------------------------------------
!
      return
      end
