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

!**s/r dcmip_tropical_cyclone - Setup for Tropical cyclone (DCMIP 2016)

      subroutine dcmip_tropical_cyclone (F_u,F_v,F_w,F_zd,F_tv,F_qv,F_topo,F_s,F_ps, &
                                         Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L)

      use tropical_cyclone
      use geomh
      use glb_ld
      use cstv
      use lun
      use ver
      use ptopo
      use dynkernel_options

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer Mminx,Mmaxx,Mminy,Mmaxy,Nk

      real F_u    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_v    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_w    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_zd   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_tv   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Virtual temperature
           F_qv   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Specific humidity
           F_s    (Mminx:Mmaxx,Mminy:Mmaxy)   , &
           F_ps   (Mminx:Mmaxx,Mminy:Mmaxy)   , &
           F_topo (Mminx:Mmaxx,Mminy:Mmaxy)

      logical F_stag_L ! Staggered uv if .T. / Scalar uv if .F.

      !==========================================
      !   Setup for Tropical cyclone (DCMIP 2016)
      !==========================================

      !----------------------------------------------------------

      integer i,j,k

      real(kind=REAL64) x_a_8,y_a_8,utt_8,vtt_8,s_8(2,2),rlon_8

      real(kind=REAL64)  :: lon,     & ! Longitude (radians)
                  lat,     & ! Latitude (radians)
                  z          ! Altitude (m)

      real(kind=REAL64)  :: p          ! Pressure  (Pa)

      integer  :: zcoords    ! 0 if p coordinates are specified
                             ! 1 if z coordinates are specified

      real(kind=REAL64)  :: u,       & ! Zonal wind (m s^-1)
                  v,       & ! Meridional wind (m s^-1)
                  w,       & ! Vertical Velocity (m s^-1)
                  t,       & ! Temperature (K)
                  tv,      & ! Virtual Temperature (K)
                  thetav,  & ! Virtual potential temperature (K)
                  phis,    & ! Surface Geopotential (m^2 s^-2)
                  ps,      & ! Surface Pressure (Pa)
                  rho,     & ! Density (kg m^-3)
                  q          ! Water vapor mixing ratio (kg/kg)

      logical :: GEM_P_L

      !----------------------------------------------------------

      if (Lun_out > 0) write (Lun_out,1000)

      GEM_P_L = trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P'

      zcoords = 1
      if (GEM_P_L) zcoords = 0

      !Initial conditions: TV,S,PS,W,ZD,QV,TOPO
      !----------------------------------------
      do k = 1,Nk

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call tropical_cyclone_test (lon,lat,p,z,zcoords,Ver_z_8%t(k),Ver_a_8%t(k),Ver_b_8%t(k), &
                                              Cstv_pref_8,u,v,t,tv,thetav,phis,ps,rho,q)

                  F_tv(i,j,k) = tv
                  F_qv(i,j,k) = q

                  if (GEM_P_L) then
                     F_s (i,j) = log(ps/Cstv_pref_8)
                  else
                     F_ps(i,j) = ps
                  end if

                  F_topo(i,j) = phis

                  w = 0.

                  F_w (i,j,k) = w
                  F_zd(i,j,k) = w ! It is zero

                  if (k==Nk) F_zd(i,j,k) = 0.

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call tropical_cyclone_test (lon,lat,p,z,zcoords,Ver_z_8%t(k),Ver_a_8%t(k),Ver_b_8%t(k), &
                                              Cstv_pref_8,utt_8,vtt_8,t,tv,thetav,phis,ps,rho,q)

                  F_tv(i,j,k) = tv
                  F_qv(i,j,k) = q

                  if (GEM_P_L) then
                     F_s (i,j) = log(ps/Cstv_pref_8)
                  else
                     F_ps(i,j) = ps
                  end if

                  F_topo(i,j) = phis

                  w = 0.

                  F_w (i,j,k) = w
                  F_zd(i,j,k) = w ! It is zero

                  if (k==Nk) F_zd(i,j,k) = 0.

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

                  call tropical_cyclone_test (lon,lat,p,z,zcoords,Ver_z_8%m(k),Ver_a_8%m(k),Ver_b_8%m(k), &
                                              Cstv_pref_8,u,v,t,tv,thetav,phis,ps,rho,q)

                  F_u(i,j,k) = u

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call tropical_cyclone_test (lon,lat,p,z,zcoords,Ver_z_8%m(k),Ver_a_8%m(k),Ver_b_8%m(k), &
                                              Cstv_pref_8,utt_8,vtt_8,t,tv,thetav,phis,ps,rho,q)

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

                  call tropical_cyclone_test (lon,lat,p,z,zcoords,Ver_z_8%m(k),Ver_a_8%m(k),Ver_b_8%m(k), &
                                              Cstv_pref_8,u,v,t,tv,thetav,phis,ps,rho,q)

                  F_v(i,j,k) = v

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call tropical_cyclone_test (lon,lat,p,z,zcoords,Ver_z_8%m(k),Ver_a_8%m(k),Ver_b_8%m(k), &
                                              Cstv_pref_8,utt_8,vtt_8,t,tv,thetav,phis,ps,rho,q)

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

                  call tropical_cyclone_test (lon,lat,p,z,zcoords,Ver_z_8%m(k),Ver_a_8%m(k),Ver_b_8%m(k), &
                                              Cstv_pref_8,u,v,t,tv,thetav,phis,ps,rho,q)

                  F_u(i,j,k) = u

               end do

            else

               do i = 1,l_niu

                  x_a_8 = geomh_xu_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call tropical_cyclone_test (lon,lat,p,z,zcoords,Ver_z_8%m(k),Ver_a_8%m(k),Ver_b_8%m(k), &
                                              Cstv_pref_8,utt_8,vtt_8,t,tv,thetav,phis,ps,rho,q)

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

                  call tropical_cyclone_test (lon,lat,p,z,zcoords,Ver_z_8%m(k),Ver_a_8%m(k),Ver_b_8%m(k), &
                                              Cstv_pref_8,u,v,t,tv,thetav,phis,ps,rho,q)

                  F_v(i,j,k) = v

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call tropical_cyclone_test (lon,lat,p,z,zcoords,Ver_z_8%m(k),Ver_a_8%m(k),Ver_b_8%m(k), &
                                              Cstv_pref_8,utt_8,vtt_8,t,tv,thetav,phis,ps,rho,q)

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
      /,'USE INITIAL CONDITIONS FOR IDEALIZED TROPICAL CYCLONE: (S/R DCMIP_TROPICAL_CYCLONE)', &
      /,'===================================================================================',/,/)

      end subroutine dcmip_tropical_cyclone
