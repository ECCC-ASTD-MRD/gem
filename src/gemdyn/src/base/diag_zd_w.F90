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

!** s/r diag_zd_w - Computes model vertical velocities zd
!                   and w diagnostically.
      
      subroutine diag_zd_w ( F_zd, F_w, F_u, F_v, F_t, F_s, F_sls,&
                             Minx, Maxx, Miny, Maxy, Nk, F_zd_L, F_w_L )
      use gem_options
      use init_options
      use geomh
      use glb_ld
      use gmm_geof
      use lun
      use tdpack
      use ver
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(out) :: F_zd, F_w
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(inout) :: F_u, F_v
      real, dimension(Minx:Maxx,Miny:Maxy), intent(in   ) :: F_sls
      real, dimension(Minx:Maxx,Miny:Maxy), intent(inout) :: F_s
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(in) :: F_t
      logical, intent(in) ::  F_zd_L, F_w_L
!
!authors
!      C.Girard & A.Plante, August 2011, based on routine uv2psd
!
!arguments
!___________________________________________________________________
!        |                                             |           |
! NAME   |             DESCRIPTION                     | DIMENSION |
!--------|---------------------------------------------|-----------|
! F_zd   | coordinat vertical motion ( 1/s )           | 3D (Nk)   |
! F_w    | true vertical motion      ( m/s )           | 3D (Nk)   |
!--------|---------------------------------------------|-----------|
! F_u    | x component of velocity                     | 3D (Nk)   |
! F_v    | y component of velocity                     | 3D (Nk)   |
! F_t    | temperature                                 | 3D (Nk)   |
! F_s    | s log of surface pressure over constant     | 2D (1)    |
! F_zd_L | true to compute zdot                        | scal      |
! F_w_L  | true to compute w                           | scal      |
!________|_____________________________________________|___________|
!
      integer i, j, k, kp, i0, in, j0, jn
      real(kind=REAL64), parameter :: half=0.5d0
      real(kind=REAL64) c1,c2
      real lnpi,pbX,pbY,adv,advl,RoverG
      real, dimension(l_minx:l_maxx,l_miny:l_maxy)      :: UdpX,VdpY
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,4)    :: sbar
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,nk)   :: div,pi_t,&
                                                    lnpi_t,xd,lapse
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,0:nk) :: div_i
!     ________________________________________________________________
!
      call gem_xch_halo (F_s, l_minx, l_maxx, l_miny, l_maxy, 1)
      call gem_xch_halo (F_u, l_minx, l_maxx, l_miny, l_maxy, nk)
      call gem_xch_halo (F_v, l_minx, l_maxx, l_miny, l_maxy, nk)

      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox-1
            sbar(i,j,1) = (F_s(i,j) + F_s(i+1,j))*half
            sbar(i,j,2) = (F_sls(i,j) + F_sls(i+1,j))*half
         end do
      end do
      if (l_east)  then
         sbar(l_ni+G_halox,:,1)= F_s(l_ni+G_halox,:)
         sbar(l_ni+G_halox,:,2)= F_sls(l_ni+G_halox,:)
      endif
      
      do j=1-G_haloy,l_nj+G_haloy-1
         do i=1-G_halox,l_ni+G_halox
            sbar(i,j,3) = (F_s(i,j) + F_s(i,j+1))*half
            sbar(i,j,4) = (F_sls(i,j) + F_sls(i,j+1))*half
         end do
      end do
      if (l_north)  then
         sbar(:,l_nj+G_haloy,3)= F_s(:,l_nj+G_haloy)
         sbar(:,l_nj+G_haloy,4)= F_sls(:,l_nj+G_haloy)
      endif

      call gem_xch_halo (sbar, l_minx, l_maxx, l_miny, l_maxy, 4)
      
! local grid setup for final results
      i0 = 1
      in = l_ni
      j0 = 1
      jn = l_nj
      if (l_west)  i0 = 2 - G_halox
      if (l_east)  in = l_niu + G_halox
      if (l_south) j0 = 2 - G_haloy
      if (l_north) jn = l_njv + G_haloy

      do k = 1, Nk

! Compute U*dpi/dz and V*dpi/dz
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               lnpi = Ver_z_8%m(k)+Ver_b_8%m(k)*sbar(i,j,1)+Ver_c_8%m(k)*sbar(i,j,2)
               pbX = exp(lnpi)
               UdpX(i,j) = F_u(i,j,k)*pbX*(1.d0+Ver_dbdz_8%m(k)*sbar(i,j,1)+Ver_dcdz_8%m(k)*sbar(i,j,2))
            end do
         end do
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               lnpi = Ver_z_8%m(k)+Ver_b_8%m(k)*sbar(i,j,3)+Ver_c_8%m(k)*sbar(i,j,4)
               pbY = exp(lnpi)
               VdpY(i,j) = F_v(i,j,k)*geomh_cyv_8(j)*pbY*(1.d0+Ver_dbdz_8%m(k)*sbar(i,j,3)+Ver_dcdz_8%m(k)*sbar(i,j,4))
            end do
         end do

! Compute DIV(V*dpi/dz)
         do j=j0,jn
            do i=i0,in
               div(i,j,k) = (UdpX(i,j)-UdpX(i-1,j))*geomh_invDX_8(j) &
                          + (VdpY(i,j)-VdpY(i,j-1))*geomh_invcy_8(j)*geomh_invDY_8
            end do
         end do

! Compute lnpi_t and pi_t
         do j=j0,jn
            do i=i0,in
               lnpi_t(i,j,k) = Ver_z_8%t(k) + Ver_b_8%t(k)*F_s(i,j) &
                                            + Ver_c_8%t(k)*F_sls(i,j)
               pi_t(i,j,k) = exp(lnpi_t(i,j,k))
            end do
         end do

      end do

      if (Zdot_divHLM_L) then
! Compute lapse rate
         call verder (lapse,F_t,lnpi_t,2.0,2.0,&
                      l_minx,l_maxx,l_miny,l_maxy,Nk,i0,in,j0,jn)

! Scale divergence by the lapse rate
         do k=1,Nk
            do j=j0,jn
               do i=i0,in
                 div(i,j,k)=div(i,j,k)*(tanh(2.*lapse(i,j,k)/pi_8)+1.)/2.
               end do
            end do
         end do
      end if

!     Vertical integral of divergence
      div_i(:,:,0) = 0.
      do k=1,Nk
         do j=j0,jn
            do i=i0,in
               div_i(i,j,k)=div_i(i,j,k-1)+div(i,j,k)*Ver_dz_8%m(k)
            end do
         end do
      end do

! Compute ZDOT
      if (F_zd_L) then
         if (Lun_debug_L.and.F_zd_L) write (Lun_out,1000)
         do k=1,Nk-1
            do j=j0,jn
               do i=i0,in
! Ver_b_8%t(k) comes from the time tendency of ln(pi), therefore
! Ver_c_8 is not present since this part of ln(pi) Bl*sl is constant.
                  F_zd(i,j,k) = ( Ver_b_8%t(k) * div_i(i,j,Nk) /&
                                  exp(Ver_z_8%m(Nk+1)+F_s(i,j)) &
                                  - div_i(i,j,k)/pi_t(i,j,k) ) /&
                                  ( 1.+Ver_dbdz_8%t(k)*F_s(i,j) &
                                  +Ver_dcdz_8%t(k)*F_sls(i,j) )
               end do
            end do
         end do
      end if

! Compute W
      if (F_w_L) then
         if (Lun_debug_L.and.F_w_L ) write (Lun_out,1001)

!        Compute W=-omega/(g*ro)      ro=p/RT
         RoverG=rgasd_8/grav_8
         do k=1,Nk
            kp = min(k+1,Nk)
            do j=j0,jn
               c1 = geomh_cyv_8(j  ) ! pgi optimizer is bugged without those 2 lines
               c2 = geomh_cyv_8(j-1)
               do i=i0,in
                  !ADV = V*grad(s) = DIV(s*Vbarz)-s*DIV(Vbarz)
                  adv = 0.5 * ( geomh_invDX_8(j) *  &
                      ( (F_u(i  ,j,kp)+F_u(i  ,j,k))*(sbar(i  ,j,1)-(F_s(i,j)+0.d0))   &
                       -(F_u(i-1,j,kp)+F_u(i-1,j,k))*(sbar(i-1,j,1)-(F_s(i,j)+0.d0)) ) &
                      + geomh_invcy_8(j) * geomh_invDY_8 *  &
                      ( (F_v(i,j  ,kp)+F_v(i,j  ,k))*c1*(sbar(i,j,3  )-(F_s(i,j)+0.d0))   &
                       -(F_v(i,j-1,kp)+F_v(i,j-1,k))*c2*(sbar(i,j-1,3)-(F_s(i,j)+0.d0)) ) )
                  advl = 0.5 * ( geomh_invDX_8(j) *  &
                       ( (F_u(i  ,j,kp)+F_u(i  ,j,k))*(sbar(i  ,j,2)-(F_sls(i,j)+0.d0))   &
                        -(F_u(i-1,j,kp)+F_u(i-1,j,k))*(sbar(i-1,j,2)-(F_sls(i,j)+0.d0)) ) &
                       + geomh_invcy_8(j) * geomh_invDY_8 *  &
                       ( (F_v(i,j  ,kp)+F_v(i,j  ,k))*c1*(sbar(i,j,4  )-(F_sls(i,j)+0.d0))   &
                        -(F_v(i,j-1,kp)+F_v(i,j-1,k))*c2*(sbar(i,j-1,4)-(F_sls(i,j)+0.d0)) ) )
                  xd (i,j,k) = Ver_b_8%t(k)*adv + Ver_c_8%t(k)*advl - div_i(i,j,k)/pi_t(i,j,k)
                  F_w(i,j,k) = -RoverG*F_t(i,j,k)*xd(i,j,k)
               end do
            end do
         end do

         do j=j0,jn
            do i=i0,in
               F_w(i,j,Nk) = -RoverG*F_t(i,j,Nk) * ( Ver_wpstar_8(Nk) * &
                              xd(i,j,Nk)+ Ver_wmstar_8(Nk)*xd(i,j,Nk-1) )
            end do
         end do
      end if

1000  format(3X,'COMPUTE DIAGNOSTIC ZDT1: (S/R DIAG_ZD_W)')
1001  format(3X,'COMPUTE DIAGNOSTIC WT1:  (S/R DIAG_ZD_W)')
!
!     ________________________________________________________________
!
      return
      end subroutine diag_zd_w
