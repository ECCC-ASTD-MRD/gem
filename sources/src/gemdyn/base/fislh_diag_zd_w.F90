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

!** s/r diag_zd_w_H - Computes model vertical velocities zd and w diagnostically.
!                     Height-type vertical coordinate

      subroutine fislh_diag_zd_w (F_zd, F_w, F_u, F_v, F_t, F_q,  &
                                  Minx, Maxx, Miny, Maxy, Nk, F_zd_L, F_w_L )
      use dyn_fisl_options
      use gem_options
      use geomh
      use glb_ld
      use lun
      use metric
      use tdpack
      use ver
      implicit none
#include <arch_specific.hf>

      integer, intent(in) ::  Minx, Maxx, Miny, Maxy, Nk
      logical, intent(in) ::  F_zd_L, F_w_L
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(out)   :: F_zd, F_w
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(inout) :: F_u, F_v, F_t
      real, dimension(Minx:Maxx,Miny:Maxy,Nk+1),intent(inout) :: F_q

      integer :: i, j, k, kp, km, i0, in, j0, jn
      real, dimension(Minx:Maxx,Miny:Maxy,Nk) :: rau
      real, dimension(Minx:Maxx,Miny:Maxy)    :: rJzX, rJzY, rJzZ
!     ________________________________________________________________
!
      if(.not.(F_zd_L.or.F_w_L)) then
         F_zd = 0.
         F_w = 0.
         return
      end if

      if (Lun_debug_L.and.F_zd_L) write (Lun_out,1000)
      if (Lun_debug_L.and.F_w_L ) write (Lun_out,1001)
!
! Halo exchange needed because scope below goes one point in halo
!
      call rpn_comm_xch_halo (F_u,  l_minx,l_maxx,l_miny,l_maxy, l_niu,l_nj, Nk,   &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
      call rpn_comm_xch_halo (F_v,  l_minx,l_maxx,l_miny,l_maxy, l_ni,l_njv, Nk,   &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
      call rpn_comm_xch_halo( F_t,  l_minx,l_maxx,l_miny,l_maxy, l_ni, l_nj ,Nk,   &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
      call rpn_comm_xch_halo( F_q,  l_minx,l_maxx,l_miny,l_maxy, l_ni, l_nj ,Nk+1, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

!     Initializations

!     local grid setup for final results
      i0 = 1
      in = l_ni
      j0 = 1
      jn = l_nj
      if(l_west)  i0 = 2
      if(l_east)  in = in-1
      if(l_south) j0 = 2
      if(l_north) jn = jn-1

      if (F_zd_L) then

!     CALCULATION of rau=p/(Rd*Tv)

         rau=0.
         do k=1,Nk
            km=max(k-1,1)
            do j=j0-1,jn+1
               do i=i0-1,in+1
                  rau(i,j,k) = exp(F_q(i,j,k)/(rgasd_8*Cstv_Tstr_8)+lg_pstar(i,j,k)) &
                             / (rgasd_8*(Ver_wp_8%m(k)*F_t(i,j,k)+Ver_wm_8%m(k)*F_t(i,j,km)))
               end do
            end do
         end do

!        CALCULATION of Zdot

         F_zd=0.0
         rJzZ=0.0
         do k=1,Nk-1
            km=max(k-1,1)
            do j=j0,jn
               do i=i0-1,in
                  rJzX(i,j) = 0.5d0*(rau(i+1,j,k)*(ztht(i+1,j,k)-ztht(i+1,j,k-1)) &
                            + rau(i  ,j,k)*(ztht(i  ,j,k)-ztht(i  ,j,k-1)))*Ver_idz_8%m(k)
               end do
            end do
            do j=j0-1,jn
               do i=i0,in
                  rJzY(i,j) = 0.5d0*(rau(i,j+1,k)*(ztht(i,j+1,k)-ztht(i,j+1,k-1)) &
                            +rau(i,j  ,k)*(ztht(i,j  ,k)-ztht(i,j  ,k-1)))*Ver_idz_8%m(k)
               end do
            end do
            do j=j0,jn
               do i=i0,in
                  F_zd(i,j,k) = rJzZ(i,j)*F_zd(i,j,km) - Ver_dz_8%m(k)*( &
                                (rJzX(i  ,j)*F_u(i  ,j,k)                   &
                               -rJzX(i-1,j)*F_u(i-1,j,k))*geomh_invDX_8(j) &
                              +(rJzY(i,j  )*F_v(i,j  ,k)*geomh_cyv_8(j  )  &
                               -rJzY(i,j-1)*F_v(i,j-1,k)*geomh_cyv_8(j-1))*geomh_invDYM_8(j) )
                  rJzZ(i,j) = 0.5d0*(rau(i,j,k+1)+rau(i,j,k))* &
                              (zmom(i,j,k+1)-zmom(i,j,k))*Ver_idz_8%t(k)
                  F_zd(i,j,k) = F_zd(i,j,k)/rJzZ(i,j)
               end do
            end do
         end do
      else
         F_zd = 0.
      end if

      if(F_w_L) then

!        CALCULATION of W depending on Zdot

         do k=1,Nk
            km=max(k-1,1)
            kp=min(k+1,Nk)
            do j=j0,jn
               do i=i0,in
                  F_w(i,j,k) = 0.25d0* ( &
                              (F_u(i,j,kp)*mc_Jx(i,j,kp)+F_u(i-1,j,kp)*mc_Jx(i-1,j,kp))   &
                             +(F_u(i,j,k )*mc_Jx(i,j,k )+F_u(i-1,j,k )*mc_Jx(i-1,j,k ))   &
                             +(F_v(i,j,kp)*mc_Jy(i,j,kp)+F_v(i,j-1,kp)*mc_Jy(i,j-1,kp))   &
                             +(F_v(i,j,k )*mc_Jy(i,j,k )+F_v(i,j-1,k )*mc_Jy(i,j-1,k )) ) &
                             +(Ver_wpstar_8(k)*F_zd(i,j,k)+Ver_wmstar_8(k)*F_zd(i,j,km))  &
                             *(zmom(i,j,k+1)-zmom(i,j,k))*Ver_idz_8%t(k)
               end do
            end do
         end do
      else
         F_w = 0.
      end if

 1000 format(3X,'COMPUTE DIAGNOSTIC ZDT1: (S/R DIAG_ZD_W_H)')
 1001 format(3X,'COMPUTE DIAGNOSTIC WT1:  (S/R DIAG_ZD_W_H)')
!
!     ________________________________________________________________
!
      return
      end
