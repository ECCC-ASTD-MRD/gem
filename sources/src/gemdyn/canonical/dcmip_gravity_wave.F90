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

!**s/r dcmip_gravity_wave - Setup for Gravity wave on a small planet along the equator (DCMIP 2012)

      subroutine dcmip_gravity_wave (F_u,F_v,F_zd,F_tv,F_q,F_topo,F_s,F_thbase, &
                                     Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L)

      use dcmip_2012_init_1_2_3
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
      real F_u     (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Scalar u
           F_v     (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Scalar v
           F_zd    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_tv    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Virtual temperature
           F_q     (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_s     (Mminx:Mmaxx,Mminy:Mmaxy),    &
           F_topo  (Mminx:Mmaxx,Mminy:Mmaxy),    &
           F_thbase(Mminx:Mmaxx,Mminy:Mmaxy,Nk)

      logical F_stag_L ! Staggered uv if .T. / Scalar uv if .F.

      !object
      !==========================================================================
      !   Setup for Gravity wave on a small planet along the equator (DCMIP 2012)
      !==========================================================================

      !-------------------------------------------------------------

      integer i,j,k

      real(8) x_a_8,y_a_8,utt_8,vtt_8,s_8(2,2),rlon_8

      real(8)  :: lon, &          ! Longitude (radians)
                  lat, &          ! Latitude (radians)
                  z               ! Height (m)

      real(8)  :: p               ! Pressure  (Pa)

      integer  :: zcoords         ! 0 or 1 see below

      real(8)  :: u, &            ! Zonal wind (m s^-1)
                  v, &            ! Meridional wind (m s^-1)
                  w, &            ! Vertical Velocity (m s^-1)
                  t, &            ! Temperature (K)
                  tv,&            ! Virtual Temperature (K)
                  phis, &         ! Surface Geopotential (m^2 s^-2)
                  ps, &           ! Surface Pressure (Pa)
                  rho, &          ! density (kg m^-3)
                  q               ! Specific Humidity (kg/kg)

      ! if zcoords = 1, then we use z and output z
      ! if zcoords = 0, then we use p

      real(8) :: theta_base      ! Base Potential temperature

      !-------------------------------------------------------------

      if (Lun_out > 0) write (Lun_out,1000)

      zcoords = 0

      !Initial conditions: T,ZD,Q,S,TOPO,THBASE
      !----------------------------------------
      do k = 1,Nk

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call test3_gravity_wave (lon,lat,p,z,zcoords,Ver_a_8%t(k),Ver_b_8%t(k),Cstv_pref_8,u,v,w,t,tv,phis,ps,rho,q,theta_base)

                  F_tv    (i,j,k) = tv
                  F_q     (i,j,k) = q
                  F_s     (i,j)   = log(ps/Cstv_pref_8)
                  F_topo  (i,j)   = phis
                  F_zd    (i,j,k) = w ! It is zero
                  F_thbase(i,j,k) = theta_base

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call test3_gravity_wave (lon,lat,p,z,zcoords,Ver_a_8%t(k),Ver_b_8%t(k),Cstv_pref_8,utt_8,vtt_8,w,t,tv,phis,ps,rho,q,theta_base)

                  F_tv     (i,j,k) = tv
                  F_q      (i,j,k) = q
                  F_s      (i,j)   = log(ps/Cstv_pref_8)
                  F_topo   (i,j)   = phis
                  F_zd     (i,j,k) = w ! It is zero
                  F_thbase(i,j,k) = theta_base

               end do

            end if

         end do

      end do

      !######################
      if (.not.F_stag_L) then
      !######################

      !Initial conditions: U True
      !--------------------------
      do k = 1,Nk

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call test3_gravity_wave (lon,lat,p,z,zcoords,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,u,v,w,t,tv,phis,ps,rho,q,theta_base)

                  F_u(i,j,k) = u

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call test3_gravity_wave (lon,lat,p,z,zcoords,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,utt_8,vtt_8,w,t,tv,phis,ps,rho,q,theta_base)

                  u = s_8(1,1)*utt_8 + s_8(1,2)*vtt_8

                  F_u(i,j,k) = u

               end do

            end if

         end do

      end do

      !Initial conditions: V True
      !--------------------------
      do k = 1,Nk

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call test3_gravity_wave (lon,lat,p,z,zcoords,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,u,v,w,t,tv,phis,ps,rho,q,theta_base)

                  F_v(i,j,k) = v

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call test3_gravity_wave (lon,lat,p,z,zcoords,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,utt_8,vtt_8,w,t,tv,phis,ps,rho,q,theta_base)

                  v = s_8(2,1)*utt_8 + s_8(2,2)*vtt_8

                  F_v(i,j,k) = v

               end do

            end if

         end do

      end do

      !######################
      else !F_stag_L
      !######################

      !Initial conditions: U True
      !--------------------------
      do k = 1,Nk

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_niu

                  lon = geomh_xu_8(i)

                  call test3_gravity_wave (lon,lat,p,z,zcoords,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,u,v,w,t,tv,phis,ps,rho,q,theta_base)

                  F_u(i,j,k) = u

               end do

            else

               do i = 1,l_niu

                  x_a_8 = geomh_xu_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call test3_gravity_wave (lon,lat,p,z,zcoords,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,utt_8,vtt_8,w,t,tv,phis,ps,rho,q,theta_base)

                  u = s_8(1,1)*utt_8 + s_8(1,2)*vtt_8

                  F_u(i,j,k) = u

               end do

            end if

         end do

      end do

      !Initial conditions: V True
      !--------------------------
      do k = 1,Nk

         do j = 1,l_njv

            lat   = geomh_yv_8(j)
            y_a_8 = geomh_yv_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call test3_gravity_wave (lon,lat,p,z,zcoords,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,u,v,w,t,tv,phis,ps,rho,q,theta_base)

                  F_v(i,j,k) = v

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call test3_gravity_wave (lon,lat,p,z,zcoords,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,utt_8,vtt_8,w,t,tv,phis,ps,rho,q,theta_base)

                  v = s_8(2,1)*utt_8 + s_8(2,2)*vtt_8

                  F_v(i,j,k) = v

               end do

            end if

         end do

      end do

      !######################
      end if !F_stag_L
      !######################

      return

 1000 format( &
      /,'USE INITIAL CONDITIONS FOR GRAVITY WAVE ON A SMALL PLANET ALONG THE EQUATOR : (S/R DCMIP_GRAVITY_WAVE)',   &
      /,'======================================================================================================',/,/)

      end subroutine dcmip_gravity_wave
