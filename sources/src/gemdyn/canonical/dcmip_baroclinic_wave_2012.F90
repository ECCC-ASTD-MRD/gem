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

!**s/r dcmip_baroclinic_wave_2012 - Setup for Baroclinic Instability on a Small Planet (DCMIP 2012)

      subroutine dcmip_baroclinic_wave_2012 (F_u,F_v,F_zd,F_tv,F_q,F_topo,F_s,F_q1,F_q2, &
                                             Mminx,Mmaxx,Mminy,Mmaxy,Nk,Moist,X,Tracers,F_stag_L)

      use canonical
      use dcmip_initial_conditions_test_4
      use gem_options
      use geomh
      use glb_ld
      use cstv
      use lun
      use ver
      use gmm_itf_mod
      use ptopo

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer Mminx,Mmaxx,Mminy,Mmaxy,Nk
      real F_u    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Scalar u
           F_v    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Scalar v
           F_zd   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_tv   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Virtual temperature
           F_q    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_s    (Mminx:Mmaxx,Mminy:Mmaxy),    &
           F_topo (Mminx:Mmaxx,Mminy:Mmaxy),    &
           F_q1   (Mminx:Mmaxx,Mminy:Mmaxy,*),  &
           F_q2   (Mminx:Mmaxx,Mminy:Mmaxy,*)

      integer  :: Moist      ! Moist (1) or non-moist (0) test case
      real(8)  :: X          ! Scale factor
      integer  :: Tracers    ! Tracers (1) or non-tracers (0) test case

      logical  :: F_stag_L   ! Staggered uv if .T. / Scalar uv if .F.

      !object
      !==================================================================
      !   Setup for Baroclinic Instability on a Small Planet (DCMIP 2012)
      !==================================================================

      !-----------------------------------------------------------------------

      integer i,j,k

      real(8) x_a_8,y_a_8,utt_8,vtt_8,s_8(2,2),rlon_8

      real(8)  :: &
                  lon,     & ! Longitude (radians)
                  lat,     & ! Latitude (radians)
                  z,       & ! Height (m)
                  eta_GEM    ! eta_GEM

      real(8)  :: p          ! Pressure  (Pa)

      integer  :: zcoords    ! 0 or 1 see below

      real(8)  :: &
                  u,       & ! Zonal wind (m s^-1)
                  v,       & ! Meridional wind (m s^-1)
                  w,       & ! Vertical Velocity (m s^-1)
                  t,       & ! Temperature (K)
                  tv,      & ! Virtual Temperature (K)
                  phis,    & ! Surface Geopotential (m^2 s^-2)
                  ps,      & ! Surface Pressure (Pa)
                  rho,     & ! density (kg m^-3)
                  q,       & ! Specific Humidity (kg/kg)
                  q1,      & ! Tracer q1 - Potential temperature (kg/kg)
                  q2         ! Tracer q2 - Ertel's potential vorticity (kg/kg)

      !-----------------------------------------------------------------------

      if (Lun_out > 0) write (Lun_out,1000) Moist,X

      zcoords = 0

      !Initial conditions: T,ZD,Q,Q1,Q2,S,TOPO
      !---------------------------------------
      do k = 1,Nk

         !Obtain Eta GEM based on: Zeta - Zeta_s = ln(Eta)
         !------------------------------------------------
         eta_GEM = exp(Ver_z_8%t(k) - Cstv_zsrf_8)

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call test4_baroclinic_wave (Moist,X,lon,lat,p,z,zcoords,eta_GEM,u,v,w,t,tv,phis,ps,rho,q,q1,q2)

                  F_tv  (i,j,k) = tv
                  F_q   (i,j,k) = q
                  F_s   (i,j)   = log(ps/Cstv_pref_8)
                  F_topo(i,j)   = phis
                  F_zd  (i,j,k) = w ! It is zero

                  if (Tracers == 1) then
                  F_q1  (i,j,k) = q1
                  F_q2  (i,j,k) = q2
                  end if

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call test4_baroclinic_wave (Moist,X,lon,lat,p,z,zcoords,eta_GEM,utt_8,vtt_8,w,t,tv,phis,ps,rho,q,q1,q2)

                  F_tv  (i,j,k) = tv
                  F_q   (i,j,k) = q
                  F_s   (i,j)   = log(ps/Cstv_pref_8)
                  F_topo(i,j)   = phis
                  F_zd  (i,j,k) = w ! It is zero

                  if (Tracers == 1) then
                  F_q1  (i,j,k) = q1
                  F_q2  (i,j,k) = q2
                  end if

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

         !Obtain Eta GEM based on: Zeta - Zeta_s = ln(Eta)
         !------------------------------------------------
         eta_GEM = exp(Ver_z_8%m(k) - Cstv_zsrf_8)

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call test4_baroclinic_wave (Moist,X,lon,lat,p,z,zcoords,eta_GEM,u,v,w,t,tv,phis,ps,rho,q,q1,q2)

                  F_u(i,j,k) = u

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call test4_baroclinic_wave (Moist,X,lon,lat,p,z,zcoords,eta_GEM,utt_8,vtt_8,w,t,tv,phis,ps,rho,q,q1,q2)

                  u = s_8(1,1)*utt_8 + s_8(1,2)*vtt_8

                  F_u(i,j,k) = u

               end do

            end if

         end do

      end do

      !Initial conditions: V True
      !--------------------------
      do k = 1,Nk

         !Obtain Eta GEM based on: Zeta - Zeta_s = ln(Eta)
         !------------------------------------------------
         eta_GEM = exp(Ver_z_8%m(k) - Cstv_zsrf_8)

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call test4_baroclinic_wave (Moist,X,lon,lat,p,z,zcoords,eta_GEM,u,v,w,t,tv,phis,ps,rho,q,q1,q2)

                  F_v(i,j,k) = v

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call test4_baroclinic_wave (Moist,X,lon,lat,p,z,zcoords,eta_GEM,utt_8,vtt_8,w,t,tv,phis,ps,rho,q,q1,q2)

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

         !Obtain Eta GEM based on: Zeta - Zeta_s = ln(Eta)
         !------------------------------------------------
         eta_GEM = exp(Ver_z_8%m(k) - Cstv_zsrf_8)

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_niu

                  lon = geomh_xu_8(i)

                  call test4_baroclinic_wave (Moist,X,lon,lat,p,z,zcoords,eta_GEM,u,v,w,t,tv,phis,ps,rho,q,q1,q2)

                  F_u(i,j,k) = u

               end do

            else

               do i = 1,l_niu

                  x_a_8 = geomh_xu_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call test4_baroclinic_wave (Moist,X,lon,lat,p,z,zcoords,eta_GEM,utt_8,vtt_8,w,t,tv,phis,ps,rho,q,q1,q2)

                  u = s_8(1,1)*utt_8 + s_8(1,2)*vtt_8

                  F_u(i,j,k) = u

               end do

            end if

         end do

      end do

      !Initial conditions: V True
      !--------------------------
      do k = 1,Nk

         !Obtain Eta GEM based on: Zeta - Zeta_s = ln(Eta)
         !------------------------------------------------
         eta_GEM = exp(Ver_z_8%m(k) - Cstv_zsrf_8)

         do j = 1,l_njv

            lat   = geomh_yv_8(j)
            y_a_8 = geomh_yv_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call test4_baroclinic_wave (Moist,X,lon,lat,p,z,zcoords,eta_GEM,u,v,w,t,tv,phis,ps,rho,q,q1,q2)

                  F_v(i,j,k) = v

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call test4_baroclinic_wave (Moist,X,lon,lat,p,z,zcoords,eta_GEM,utt_8,vtt_8,w,t,tv,phis,ps,rho,q,q1,q2)

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
      /,'USE INITIAL CONDITIONS FOR BAROCLINIC INSTABILITY WITH DYNAMIC TRACERS: (S/R DCMIP_BAROCLINIC_WAVE_2012)',   &
      /,'========================================================================================================',/, &
        ' MOIST            = ',I1                                                                                 ,/, &
        ' X Scaling Factor = ',F5.0                                                                               ,/, &
      /,'========================================================================================================',/,/)

      end subroutine dcmip_baroclinic_wave_2012
