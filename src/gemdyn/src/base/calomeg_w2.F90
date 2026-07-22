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

!**s/p calomeg_w - compute vertical velocity in pressure coordinates
!                    from advection
!
      subroutine calomeg_w2 (F_ww,F_s,F_s0,F_sl,F_u,F_v,&
            F_zd,Minx,Maxx,Miny,Maxy,Nk)
      use dynkernel_options
      use tdpack
      use glb_ld
      use cstv
      use ver
      use type_mod
      use geomh
      use gmm_geof
      use HORgrid_options
      use gem_options
       
      implicit none
#include <arch_specific.hf>
!
      integer, intent(in) :: Minx,Maxx,Miny,Maxy, Nk
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(out) :: F_ww
      real, dimension(Minx:Maxx,Miny:Maxy),    intent(in)  :: F_s, F_sl, F_s0
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(in)  :: F_u, F_v, F_zd      

      
!objective
!	compute vertical velocity in hydrostatic pressure coordinates,
!       omega = dpi/dt
!       in GEM-P using the definition of the vertical coordinate, i.e.,
!       ln pi= zeta + B*s
!
!
!arguments
!  Name               Description
!---------------------------------------------------
! F_ww                dpi/dt (Pa/s)
! F_s                 s at time t1 
! F_s0                s at time t0
! F_u                 u-component of wind at time t1
! F_v                 v-component of wind at time t1
! F_zd                vertical motion in zeta coordinate at time t1
!

      integer i, j, k, kp, km, i0, in, j0, jn, istat
      real*8 c1, c2, adv, advl, pidot, w1, w2
      real pi_t(Minx:Maxx,Miny:Maxy,Nk)
      real  sbX(Minx:Maxx,Miny:Maxy),  sbY(Minx:Maxx,Miny:Maxy)
      real slbX(Minx:Maxx,Miny:Maxy), slbY(Minx:Maxx,Miny:Maxy)
      real*8, parameter :: half=0.5d0, one=1.0d0

!     __________________________________________________________________
!
!    Halo exchange needed because scope below goes one point in halo
!
      call rpn_comm_xch_halo (F_s, l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj , 1,   &
                              G_halox, G_haloy, G_periodx, G_periody, l_ni, 0)
      call rpn_comm_xch_halo (F_sl, l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj , 1,   &
                              G_halox, G_haloy, G_periodx, G_periody, l_ni, 0)
      call rpn_comm_xch_halo (F_u, l_minx, l_maxx, l_miny, l_maxy, l_niu,l_nj, Nk,   &
                              G_halox, G_haloy, G_periodx, G_periody, l_ni, 0)
      call rpn_comm_xch_halo (F_v, l_minx, l_maxx, l_miny, l_maxy, l_ni, l_njv, Nk,   &
                              G_halox, G_haloy, G_periodx, G_periody, l_ni, 0)
      call rpn_comm_xch_halo (F_s0, l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj , 1,   &
                              G_halox, G_haloy, G_periodx, G_periody, l_ni, 0)


!     Initializations

!     local grid setup for final results
      i0 = 1
      in = l_ni
      j0 = 1
      jn = l_nj

      if (.not. Grd_yinyang_L) then
         if (l_west)  i0 = 2
         if (l_east)  in = l_niu
         if (l_south) j0 = 2
         if (l_north) jn = l_njv
         F_ww=0.
      end if  
      F_ww=0.

!     CALCULATION of sbX, slbX, i.e. F_s and F_sl on U points
      do j=j0,jn
         do i=i0-1,in
            sbX(i,j) =(F_s(i,j)+F_s(i+1,j))*half
            slbX(i,j)=(F_sl(i,j)+F_sl(i+1,j))*half
         end do
      end do

!     CALCULATION of sbY, slbY, i.e. F_s and F_sl on V points
      do j=j0-1,jn
         do i=i0,in
            sbY(i,j) =(F_s(i,j)+F_s(i,j+1))*half
            slbY(i,j)=(F_sl(i,j)+F_sl(i,j+1))*half
         end do
      end do

      do k = 1, Nk
!        CALCULATION of pi_t (Pressure at thermodynamic levels)
         do j=j0,jn
            do i=i0,in
               pi_t(i,j,k) = exp(Ver_z_8%t(k) + Ver_b_8%t(k)*F_s(i,j) &
                                            + Ver_c_8%t(k)*F_sl(i,j))
            end do
         end do
      enddo 
    

      do k=1,Nk
         kp=min(k+1,Nk)
         km=max(k-1,1)
         do j=j0,jn
            c1=geomh_cyv_8(j  ) 
            c2=geomh_cyv_8(j-1)
            do i=i0,in
               !ADV = V*grad(s) = DIV(s*V)-s*DIV(V)
               adv = half * ( geomh_invDX_8(j) *  &
                     ( (F_u(i  ,j,kp)+F_u(i  ,j,k))*(sbX(i  ,j)-F_s(i,j))   &
                      -(F_u(i-1,j,kp)+F_u(i-1,j,k))*(sbX(i-1,j)-F_s(i,j)) ) &
                      + geomh_invcy_8(j) * geomh_invDY_8 *  &
                      ( (F_v(i,j  ,kp)+F_v(i,j  ,k))*c1*(sbY(i,j  )-F_s(i,j))  &
                      -(F_v(i,j-1,kp)+F_v(i,j-1,k))*c2*(sbY(i,j-1)-F_s(i,j)) ) )

               advl = half * ( geomh_invDX_8(j) *  &
                       ( (F_u(i  ,j,kp)+F_u(i  ,j,k))*(slbX(i  ,j)-F_sl(i,j))   &
                       -(F_u(i-1,j,kp)+F_u(i-1,j,k))*(slbX(i-1,j)-F_sl(i,j)) ) &
                       + geomh_invcy_8(j) * geomh_invDY_8 *  &
                       ( (F_v(i,j  ,kp)+F_v(i,j  ,k))*c1*(slbY(i,j  )-F_sl(i,j))   &
                       -(F_v(i,j-1,kp)+F_v(i,j-1,k))*c2*(slbY(i,j-1)-F_sl(i,j)) ) )
               
               ! Interpolating F_zd to momentum levels
               w1=Ver_wpstar_8(k)*F_zd(i,j,k)+Ver_wmstar_8(k)*F_zd(i,j,km)

               ! Computing (\partial s)/(\partial t)
               w2=(F_s(i,j)-F_s0(i,j))/Cstv_dt_8

               pidot=(w1*(one+F_s(i,j)*Ver_dbdz_8%t(k))+&
                  Ver_b_8%t(k)*adv + Ver_c_8%t(k)*advl + w1*F_sl(i,j)*Ver_dcdz_8%t(k))

               F_ww(i,j,k)=pi_t(i,j,k)*(pidot+w2*Ver_b_8%t(k))
            end do
         end do
      end do


!     __________________________________________________________________
!
      return
      end
