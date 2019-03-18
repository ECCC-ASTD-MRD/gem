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

!**s/r dcmip_baroclinic_wave_2016 - Setup for Baroclinic wave with Toy Terminal Chemistry (DCMIP 2016)

      subroutine dcmip_baroclinic_wave_2016 (F_u,F_v,F_w,F_tv,F_zd,F_s,F_topo,F_q,F_cl,F_cl2, &
                                             Mminx,Mmaxx,Mminy,Mmaxy,Nk,Deep,Moist,Pertt,X,F_stag_L)

      use baroclinic_wave
      use Terminator
      use geomh
      use glb_ld
      use cstv
      use lun
      use ver
      use ptopo

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer Mminx,Mmaxx,Mminy,Mmaxy,Nk
      real F_u   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Scalar u
           F_v   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Scalar v
           F_w   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_tv  (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Virtual temperature
           F_zd  (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_s   (Mminx:Mmaxx,Mminy:Mmaxy),    &
           F_topo(Mminx:Mmaxx,Mminy:Mmaxy),    &
           F_q   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_cl  (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_cl2 (Mminx:Mmaxx,Mminy:Mmaxy,Nk)

      integer  :: Deep,      & ! Deep (1) or Shallow (0) test case
                  Moist,     & ! Moist (1) or Dry (0) test case
                  Pertt        ! Perturbation type

      real(8)  :: X            ! Earth scaling parameter

      logical  :: F_stag_L     ! Staggered uv if .T. / Scalar uv if .F.

      !object
      !===================================================================
      !   Setup for Baroclinic wave with Terminator chemistry (DCMIP 2016)
      !===================================================================

      !----------------------------------------------------------

      integer i,j,kk

      real(8) x_a_8,y_a_8,utt_8,vtt_8,s_8(2,2),rlon_8

      real(8)  :: lon,     & ! Longitude (radians)
                  lat,     & ! Latitude (radians)
                  z,       & ! Altitude (m)
                  lon_d,   & ! Longitude (degrees)
                  lat_d      ! Latitude (degrees)

      real(8)  :: p          ! Pressure  (Pa)

      integer  :: zcoords    ! 0 if p coordinates are specified
                             ! 1 if z coordinates are specified

      real(8)  :: u,       & ! Zonal wind (m s^-1)
                  v,       & ! Meridional wind (m s^-1)
                  w,       & ! Vertical Velocity (m s^-1)
                  t,       & ! Temperature (K)
                  tv,      & ! Virtual Temperature (K)
                  thetav,  & ! Virtual potential temperature (K)
                  phis,    & ! Surface Geopotential (m^2 s^-2)
                  ps,      & ! Surface Pressure (Pa)
                  rho,     & ! density (kg m^-3)
                  q,       & ! water vapor mixing ratio (kg/kg)
                  cl,cl2     ! molar mixing ratio of cl and cl2

      real(8), parameter :: radians_to_degrees = 180.0_8/pi

      !----------------------------------------------------------

      if (Lun_out > 0) write (Lun_out,1000) X

      zcoords = 0  ! p coordinates are specified

      w = 0 !Vertical Velocity is zero

      !Initial conditions: T,ZD,W,Q,S,TOPO
      !-----------------------------------
      do kk = 1,Nk

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call baroclinic_wave_test (Deep,Moist,Pertt,X,lon,lat,p,z,zcoords,Ver_a_8%t(kk),Ver_b_8%t(kk), &
                                             Cstv_pref_8,u,v,t,tv,thetav,phis,ps,rho,q)

                  F_tv  (i,j,kk) = tv
                  F_q   (i,j,kk) = q
                  F_s   (i,j)    = log(ps/Cstv_pref_8)
                  F_topo(i,j)    = phis
                  F_zd  (i,j,kk) = w ! It is zero
                  F_w   (i,j,kk) = w ! It is zero

                  lon_d = lon * radians_to_degrees
                  lat_d = lat * radians_to_degrees

                  call initial_value_Terminator (lat_d,lon_d,cl,cl2)

                  !Chemical tracers
                  !----------------
                  F_cl (i,j,kk) = cl
                  F_cl2(i,j,kk) = cl2

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call baroclinic_wave_test (Deep,Moist,Pertt,X,lon,lat,p,z,zcoords,Ver_a_8%t(kk),Ver_b_8%t(kk), &
                                             Cstv_pref_8,utt_8,vtt_8,t,tv,thetav,phis,ps,rho,q)

                  F_tv  (i,j,kk) = tv
                  F_q   (i,j,kk) = q
                  F_s   (i,j)    = log(ps/Cstv_pref_8)
                  F_topo(i,j)    = phis
                  F_zd  (i,j,kk) = w ! It is zero
                  F_w   (i,j,kk) = w ! It is zero

                  lon_d = lon * radians_to_degrees
                  lat_d = lat * radians_to_degrees

                  call initial_value_Terminator (lat_d,lon_d,cl,cl2)

                  !Chemical tracers
                  !----------------
                  F_cl (i,j,kk) = cl
                  F_cl2(i,j,kk) = cl2

               end do

            end if

         end do

      end do

      !######################
      if (.not.F_stag_L) then
      !######################

      !Initial conditions: U True
      !--------------------------
      do kk = 1,Nk

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call baroclinic_wave_test (Deep,Moist,Pertt,X,lon,lat,p,z,zcoords,Ver_a_8%m(kk),Ver_b_8%m(kk), &
                                             Cstv_pref_8,u,v,t,tv,thetav,phis,ps,rho,q)

                  F_u(i,j,kk) = u

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call baroclinic_wave_test (Deep,Moist,Pertt,X,lon,lat,p,z,zcoords,Ver_a_8%m(kk),Ver_b_8%m(kk), &
                                             Cstv_pref_8,utt_8,vtt_8,t,tv,thetav,phis,ps,rho,q)

                  u = s_8(1,1)*utt_8 + s_8(1,2)*vtt_8

                  F_u(i,j,kk) = u

               end do

            end if

         end do

      end do

      !Initial conditions: V True
      !--------------------------
      do kk = 1,Nk

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call baroclinic_wave_test (Deep,Moist,Pertt,X,lon,lat,p,z,zcoords,Ver_a_8%m(kk),Ver_b_8%m(kk), &
                                             Cstv_pref_8,u,v,t,tv,thetav,phis,ps,rho,q)

                  F_v(i,j,kk) = v

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call baroclinic_wave_test (Deep,Moist,Pertt,X,lon,lat,p,z,zcoords,Ver_a_8%m(kk),Ver_b_8%m(kk), &
                                             Cstv_pref_8,utt_8,vtt_8,t,tv,thetav,phis,ps,rho,q)

                  v = s_8(2,1)*utt_8 + s_8(2,2)*vtt_8

                  F_v(i,j,kk) = v

               end do

            end if

         end do

      end do

      !######################
      else !F_stag_L
      !######################

      !Initial conditions: U True
      !--------------------------
      do kk = 1,Nk

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_niu

                  lon = geomh_xu_8(i)

                  call baroclinic_wave_test (Deep,Moist,Pertt,X,lon,lat,p,z,zcoords,Ver_a_8%m(kk),Ver_b_8%m(kk), &
                                             Cstv_pref_8,u,v,t,tv,thetav,phis,ps,rho,q)

                  F_u(i,j,kk) = u

               end do

            else

               do i = 1,l_niu

                  x_a_8 = geomh_xu_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call baroclinic_wave_test (Deep,Moist,Pertt,X,lon,lat,p,z,zcoords,Ver_a_8%m(kk),Ver_b_8%m(kk), &
                                             Cstv_pref_8,utt_8,vtt_8,t,tv,thetav,phis,ps,rho,q)

                  u = s_8(1,1)*utt_8 + s_8(1,2)*vtt_8

                  F_u(i,j,kk) = u

               end do

            end if

         end do

      end do

      !Initial conditions: V True
      !--------------------------
      do kk = 1,Nk

         do j = 1,l_njv

            lat   = geomh_yv_8(j)
            y_a_8 = geomh_yv_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call baroclinic_wave_test (Deep,Moist,Pertt,X,lon,lat,p,z,zcoords,Ver_a_8%m(kk),Ver_b_8%m(kk), &
                                             Cstv_pref_8,u,v,t,tv,thetav,phis,ps,rho,q)

                  F_v(i,j,kk) = v

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call baroclinic_wave_test (Deep,Moist,Pertt,X,lon,lat,p,z,zcoords,Ver_a_8%m(kk),Ver_b_8%m(kk), &
                                             Cstv_pref_8,utt_8,vtt_8,t,tv,thetav,phis,ps,rho,q)

                  v = s_8(2,1)*utt_8 + s_8(2,2)*vtt_8

                  F_v(i,j,kk) = v

               end do

            end if

         end do

      end do

      !######################
      end if !F_stag_L
      !######################

      return

 1000 format( &
      /,'USE INITIAL CONDITIONS FOR MOIST BAROCLINIC WAVE WITH TERMINATOR CHEMISTRY: (S/R DCMIP_BAROCLINIC_WAVE_2016)',   &
      /,'============================================================================================================',/, &
        ' X Scaling Factor = ',F5.0                                                                                   ,/, &
      /,'============================================================================================================',/,/)

      end subroutine dcmip_baroclinic_wave_2016
