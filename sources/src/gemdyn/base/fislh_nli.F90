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
!**s/r nli_H  - compute non-linear terms:  Nu, Nv, Nt, Nc, Nw,
!             - compute full right-hand side of Helmholtz eqn: Rp=Rc-Nc
!             - Height-type vertical coordinate
!
!**********************************************************************
!
      subroutine fislh_nli ( F_nu , F_nv , F_nt , F_nw , F_nc ,        &
                             F_u  , F_v  , F_t  , F_zd , F_q  , &
                             F_rc , F_rt , F_rhs, &
                             Minx, Maxx, Miny, Maxy, Nk, ni, nj, i0, j0, in, jn, icln)
      use HORgrid_options
      use gem_options
      use dyn_fisl_options
      use coriolis
      use geomh
      use tdpack
      use glb_ld
      use cstv
      use dcst
      use ptopo
      use ver
      use lun
      use metric
      use fislh_sol
      implicit none

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk, ni, nj, i0, j0, in, jn, icln
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(out)   :: F_nu,F_nv,F_nw
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(out)   :: F_nt,F_nc
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(inout) :: F_u, F_v, F_t
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(in)    :: F_zd
      real, dimension(Minx:Maxx,Miny:Maxy,Nk+1),intent(inout) :: F_q
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(inout) :: F_rc,F_rt
      real*8, dimension(ni,nj,Nk),              intent(out)   :: F_rhs

!     Author: Claude Girard, July 2017

#include <arch_specific.hf>

      integer :: i, j, k, km, kp, i0u, inu, j0v, jnv, nij, onept
      real*8  :: c0,c1,div,w1,w2,barz,barzp, t_interp, u_interp, v_interp
      real*8, dimension(i0:in,j0:jn) :: xtmp_8, ytmp_8
      real*8, parameter :: one=1.d0, zero=0.d0, half=0.5d0
!     __________________________________________________________________
!
      if (Lun_debug_L)  write(Lun_out,1000)

      if(icln > 1) then
         call rpn_comm_xch_halo( F_u, l_minx, l_maxx, l_miny, l_maxy, l_niu, l_nj ,G_nk, &
                                 G_halox, G_haloy, G_periodx, G_periody, l_ni, 0 )
         call rpn_comm_xch_halo( F_v, l_minx, l_maxx, l_miny, l_maxy, l_ni, l_njv, G_nk, &
                                 G_halox, G_haloy, G_periodx, G_periody, l_ni, 0 )
         call rpn_comm_xch_halo( F_t, l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj, G_nk, &
                                 G_halox, G_haloy, G_periodx, G_periody, l_ni, 0 )
         call rpn_comm_xch_halo( F_q, l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj, G_nk+1, &
                                 G_halox, G_haloy, G_periodx, G_periody, l_ni, 0 )
      end if

      c0 = Dcst_rayt_8**2

      c1 = grav_8 * Cstv_tau_8

      nij = (in - i0 +1)*(jn - j0 +1)

      onept= 0
      if (Grd_yinyang_L) onept = 1

!
!***********************************************************
! The nonlinear deviation of horizontal momentum equations *
!***********************************************************
!
!     Indices

      i0u = i0-1
      j0v = j0-1
      inu = l_niu-pil_e
      jnv = l_njv-pil_n

      if (l_west ) i0u = i0 - onept
      if (l_south) j0v = j0 - onept
      if (l_east ) inu = inu + onept
      if (l_north) jnv = jnv + onept

!$omp parallel private(km,kp,barz,barzp,div, &
!$omp w1,w2,t_interp,u_interp,v_interp,xtmp_8,ytmp_8)

!$omp do
      do k=1,l_nk
         km=max(k-1,1)

!        Compute Nu
!        ~~~~~~~~~~

         do j= j0, jn
            do i= i0u, inu

               barz  = Ver_wp_8%m(k)*F_t(i  ,j,k)+Ver_wm_8%m(k)*F_t(i  ,j,km)
               barzp = Ver_wp_8%m(k)*F_t(i+1,j,k)+Ver_wm_8%m(k)*F_t(i+1,j,km)
               t_interp = (barz+barzp)*half/Cstv_Tstr_8

               v_interp = 0.25d0*(F_v(i,j,k)+F_v(i,j-1,k)+F_v(i+1,j,k)+F_v(i+1,j-1,k))

               F_nu(i,j,k) = (t_interp-one) * ( F_q(i+1,j,k) - F_q(i,j,k) ) * geomh_invDX_8(j) &
                           - ( Cori_fcoru_8(i,j) + geomh_tyoa_8(j) * F_u(i,j,k) ) * v_interp

   !           Adding vertical coordinate metric terms
   !           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               F_nu(i,j,k) = F_nu(i,j,k) - (t_interp - one*isol_i) * mc_Jx(i,j,k) * ( &
                             Ver_wp_8%m(k)*half*( (F_q(i+1,j,k+1)-F_q(i+1,j,k ))*mc_iJz(i+1,j,k )   &
                                                 +(F_q(i  ,j,k+1)-F_q(i  ,j,k ))*mc_iJz(i  ,j,k ) ) &
                            +Ver_wm_8%m(k)*half*( (F_q(i+1,j,k  )-F_q(i+1,j,km))*mc_iJz(i+1,j,km)   &
                                                 +(F_q(i  ,j,k  )-F_q(i  ,j,km))*mc_iJz(i  ,j,km) ) )
            end do
         end do

!        Compute Nv
!        ~~~~~~~~~~
         do j = j0v, jnv
            do i = i0, in

               barz  = Ver_wp_8%m(k)*F_t(i,j  ,k)+Ver_wm_8%m(k)*F_t(i,j  ,km)
               barzp = Ver_wp_8%m(k)*F_t(i,j+1,k)+Ver_wm_8%m(k)*F_t(i,j+1,km)
               t_interp = (barz+barzp)*half/Cstv_Tstr_8

               u_interp = 0.25d0*(F_u(i,j,k)+F_u(i-1,j,k)+F_u(i,j+1,k)+F_u(i-1,j+1,k))

               F_nv(i,j,k) = (t_interp-one) * ( F_q(i,j+1,k) - F_q(i,j,k) ) * geomh_invDY_8 &
                           + ( Cori_fcorv_8(i,j) + geomh_tyoav_8(j) * u_interp ) * u_interp

   !           Adding vertical coordinate metric terms
   !           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               F_nv(i,j,k) = F_nv(i,j,k) - (t_interp - one*isol_i) * mc_Jy(i,j,k) * ( &
                             Ver_wp_8%m(k)*half*( (F_q(i,j+1,k+1)-F_q(i,j+1,k ))*mc_iJz(i,j+1,k )   &
                                                 +(F_q(i,j  ,k+1)-F_q(i,j  ,k ))*mc_iJz(i,j  ,k ) ) &
                            +Ver_wm_8%m(k)*half*( (F_q(i,j+1,k  )-F_q(i,j+1,km))*mc_iJz(i,j+1,km)   &
                                                 +(F_q(i,j  ,k  )-F_q(i,j  ,km))*mc_iJz(i,j  ,km) ) )
            end do
         end do

      end do
!$omp end do

      if (.not.Grd_yinyang_L) then

!        Set  Nu=0  on the east and west boundaries of the LAM grid
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         if (l_west) then
!$omp do
            do k=1,l_nk
               do j=j0,jn
                  F_nu(pil_w,j,k) = 0.
               end do
            end do
!$omp end do
         end if
         if (l_east) then
!$omp do
            do k=1,l_nk
               do j=j0,jn
                  F_nu(l_ni-pil_e,j,k) = 0.
               end do
            end do
!$omp end do
         end if

!        Set  Nv=0  on the north and south boundaries  of the LAM grid
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         if (l_south) then
!$omp do
            do k=1,l_nk
               do i=i0,in
                  F_nv(i,pil_s,k) = 0.
               end do
            end do
!$omp end do
         end if
         if (l_north) then
!$omp do
            do k=1,l_nk
               do i=i0,in
                  F_nv(i,l_nj-pil_n,k) = 0.
               end do
            end do
!$omp end do
         end if

      end if

!$omp do
      do k=1,l_nk

         km=max(k-1,1)
         do j= j0, jn
            do i= i0, in
               xtmp_8(i,j) = F_t(i,j,k)/real(Cstv_Tstr_8)
            end do
         end do
         call vlog ( ytmp_8, xtmp_8, nij )
         do j= j0, jn
            do i= i0, in

   !           Compute Nw
   !           ~~~~~~~~~~
               F_nw(i,j,k) = ((xtmp_8(i,j)-one*isol_i)*mc_iJz(i,j,k) - &
                              isol_d*Ver_idz_8%t(k))*(F_q(i,j,k+1)-F_q(i,j,k)) &
                           -  (xtmp_8(i,j)-one)*grav_8*(one-one/xtmp_8(i,j))
   !           Compute Nc
   !           ~~~~~~~~~~
               F_nc(i,j,k) = isol_d * ( half * ( mc_Ix(i,j,k)*(F_u(i,j,k)+F_u(i-1,j,k))   &
                                               + mc_Iy(i,j,k)*(F_v(i,j,k)+F_v(i,j-1,k)) ) &
                                               + mc_Iz(i,j,k)*(Ver_wpA_8(k)*F_zd(i,j,k) + &
                                                Ver_wmA_8(k)*Ver_onezero(k)*F_zd(i,j,km)) )
   !           Compute Nt
   !           ~~~~~~~~~~
               F_nt(i,j,k) = Cstv_invT_8*( ytmp_8(i,j) - (one-one/xtmp_8(i,j)) )

   !           Compute Nt'
   !           ~~~~~~~~~~~
               F_nt(i,j,k) = gama_8 * ( c1 * F_nt(i,j,k) + F_nw(i,j,k) )

            end do
         end do

      end do
!$omp end do

!$omp do
      do k=1,l_nk
         km=max(k-1,1)
         do j = j0, jn
            do i = i0, in
               w1= (Ver_idz_8%m(k) + (isol_i*mc_Iz(i,j,k) - epsi_8)*Ver_wpA_8(k))*Ver_zeronk(k)
               w2= (Ver_idz_8%m(k) - (isol_i*mc_Iz(i,j,k) - epsi_8)*Ver_wmA_8(k))*Ver_onezero(k)

   !           Compute Nc'
   !           ~~~~~~~~~~~
               div = (F_nu(i,j,k)-F_nu(i-1,j,k)) * geomh_invDXM_8(j) &
                   + (F_nv(i,j,k)*geomh_cyM_8(j)-F_nv(i,j-1,k)* &
                      geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                    + isol_i*(half * ( mc_Ix(i,j,k)*(F_nu(i,j,k)+F_nu(i-1,j,k)) &
                                     + mc_Iy(i,j,k)*(F_nv(i,j,k)+F_nv(i,j-1,k)) ))

               F_nc(i,j,k) = div - F_nc(i,j,k)*Cstv_invT_m_8

   !           Compute Nc"
   !           ~~~~~~~~~~~
               F_nc(i,j,k) = F_nc(i,j,k) + w1 * F_nt(i,j,k) - w2 * F_nt(i,j,km)

            end do
         end do
      end do
!$omp end do

!$omp do
      do k=1,l_nk
         do j= j0, jn
            do i= i0, in
               F_rhs(i,j,k) =  c0 * ( F_rc(i,j,k) - F_nc(i,j,k) )
            end do
         end do
      end do
!$omp end do

!     Apply bottom boundary conditions
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!$omp do
      do j= j0, jn
         do i= i0, in
            F_nt(i,j,l_nk) = F_nt(i,j,l_nk)-Ver_wmstar_8(G_nk)*F_nt(i,j,l_nk-1)
         end do
      end do
!$omp end do

!$omp end parallel

      if ( trim(Sol_type_S) == 'ITERATIVE_3D') then

         call  boundary ( F_rhs,F_rt,F_nt,Minx,Maxx,Miny,Maxy, &
                          Nk,ni,nj,i0,j0,in,jn,i0u,inu,j0v,jnv )
      end if


1000 format(/,5X,'COMPUTE NON-LINEAR RHS: (S/R NLI_H)')
!     __________________________________________________________________
!
      return
      end
