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

!**s/r dcmip_Schaer_mountain - Setup for Mountain waves over a Schaer-type mountain on a small planet (DCMIP 2012)

      subroutine dcmip_Schaer_mountain (F_u,F_v,F_w,F_zd,F_tv,F_qv,F_topo,F_orls,F_s,F_ps, &
                                        Mminx,Mmaxx,Mminy,Mmaxy,Nk,Shear,Set_topo_L,F_stag_L)

      use dcmip_2012_init_1_2_3
      use dyn_fisl_options
      use geomh
      use glb_ld
      use cstv
      use lun
      use ver
      use ptopo
      use dynkernel_options
      use gem_options

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
           F_topo (Mminx:Mmaxx,Mminy:Mmaxy)   , & !S grid
           F_orls (Mminx:Mmaxx,Mminy:Mmaxy)       !S grid

      integer Shear      ! If 0 then we use constant u
                         ! If 0 then we use shear flow

      logical Set_topo_L ! If .T.: Set F_topo to initialize  Topo_High
                         ! If .F.: Use F_topo initialized as Topo_low

      logical F_stag_L   ! Staggered uv if .T. / Scalar uv if .F.

      !object
      !======================================================================================
      !   Setup for Mountain waves over a Schaer-type mountain on a small planet (DCMIP 2012)
      !======================================================================================

      !--------------------------------------------------------

      integer i,j,k,i0,in,j0,jn,inu,jnv

      real(kind=REAL64) x_a_8,y_a_8,utt_8,vtt_8,s_8(2,2),rlon_8

      real(kind=REAL64)  :: lon,     & ! Longitude (radians)
                  lat,     & ! Latitude (radians)
                  z          ! Height (m)

      real(kind=REAL64)  :: p          ! Pressure  (Pa)

      integer  :: zcoords    ! 0 if p coordinates are specified
                             ! 1 if z coordinates are specified

      real(kind=REAL64)  :: u,       & ! Zonal wind (m s^-1)
                  v,       & ! Meridional wind (m s^-1)
                  w,       & ! Vertical Velocity (m s^-1)
                  t,       & ! Temperature (K)
                  tv,      & ! Virtual Temperature (K)
                  phis,    & ! Surface Geopotential (m^2 s^-2)
                  orls,    & ! Surface Geopotential (m^2 s^-2) Large-Scale
                  ps,      & ! Surface Pressure (Pa)
                  rho,     & ! Density (kg m^-3)
                  q          ! Specific Humidity (kg/kg)

      logical :: GEM_P_L

      !--------------------------------------------------------

      if (Lun_out > 0.and.Cstv_tstr_8 /= 300.0) call handle_error(-1,'DCMIP_SCHAER_MOUNTAIN','SET TSTR AS T AT EQUATOR')

      if (Lun_out > 0) write (Lun_out,1000) Shear,Set_topo_L

      i0= 1-G_halox ; in= l_ni+G_halox ; inu= l_niu+G_halox
      j0= 1-G_haloy ; jn= l_nj+G_haloy ; jnv= l_njv+G_haloy

      GEM_P_L = trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P'

      zcoords = 1
      if (GEM_P_L) zcoords = 0

      !Initial conditions: TV,S,PS,W,ZD,QV,TOPO
      !----------------------------------------
      do k = 1,Nk

         do j = j0,jn

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = i0,in

                  lon = geomh_x_8(i)

                  if (.NOT.Set_topo_L) phis = F_topo(i,j)
                  if (.NOT.Set_topo_L) orls = F_orls(i,j)

                  call test2_schaer_mountain (lon,lat,p,z,zcoords,Ver_z_8%t(k),Ver_a_8%t(k),Ver_b_8%t(k),Ver_c_8%t(k),Cstv_pref_8, &
                                              Shear,u,v,w,t,tv,phis,orls,ps,rho,q,Set_topo_L)

                  F_tv(i,j,k) = tv
                  F_qv(i,j,k) = q

                  if (GEM_P_L) then
                     F_s (i,j) = log(ps/Cstv_pref_8)
                  else
                     F_ps(i,j) = ps
                  end if

                  if (Set_topo_L) F_topo(i,j) = phis
                  if (Set_topo_L) F_orls(i,j) = orls 

                  F_w (i,j,k) = w
                  F_zd(i,j,k) = w ! It is zero

                  if (k==Nk) F_zd(i,j,k) = 0.

               end do

            else

               do i = i0,in

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  if (.NOT.Set_topo_L) phis = F_topo(i,j)
                  if (.NOT.Set_topo_L) orls = F_orls(i,j)

                  call test2_schaer_mountain (lon,lat,p,z,zcoords,Ver_z_8%t(k),Ver_a_8%t(k),Ver_b_8%t(k),Ver_c_8%t(k),Cstv_pref_8, &
                                              Shear,utt_8,vtt_8,w,t,tv,phis,orls,ps,rho,q,Set_topo_L)

                  F_tv(i,j,k) = tv
                  F_qv(i,j,k) = q

                  if (GEM_P_L) then
                     F_s (i,j) = log(ps/Cstv_pref_8)
                  else
                     F_ps(i,j) = ps
                  end if

                  if (Set_topo_L) F_topo(i,j) = phis
                  if (Set_topo_L) F_orls(i,j) = orls 

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

         do j = j0,jn

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = i0,in

                  lon = geomh_x_8(i)

                  if (.NOT.Set_topo_L) phis = F_topo(i,j)
                  if (.NOT.Set_topo_L) orls = F_orls(i,j)

                  call test2_schaer_mountain (lon,lat,p,z,zcoords,Ver_z_8%m(k),Ver_a_8%m(k),Ver_b_8%m(k),Ver_c_8%m(k),Cstv_pref_8, &
                                              Shear,u,v,w,t,tv,phis,orls,ps,rho,q,Set_topo_L)

                  F_u(i,j,k) = u

               end do

            else

               do i = i0,in

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  if (.NOT.Set_topo_L) phis = F_topo(i,j)
                  if (.NOT.Set_topo_L) orls = F_orls(i,j)

                  call test2_schaer_mountain (lon,lat,p,z,zcoords,Ver_z_8%m(k),Ver_a_8%m(k),Ver_b_8%m(k),Ver_c_8%m(k),Cstv_pref_8, &
                                              Shear,utt_8,vtt_8,w,t,tv,phis,orls,ps,rho,q,Set_topo_L)

                  u = s_8(1,1)*utt_8 + s_8(1,2)*vtt_8

                  F_u(i,j,k) = u

               end do

            end if

         end do

      end do

      !Initial conditions: V True
      !--------------------------
      do k = 1,Nk

         do j = j0,jn

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = i0,in

                  lon = geomh_x_8(i)

                  if (.NOT.Set_topo_L) phis = F_topo(i,j)
                  if (.NOT.Set_topo_L) orls = F_orls(i,j)

                  call test2_schaer_mountain (lon,lat,p,z,zcoords,Ver_z_8%m(k),Ver_a_8%m(k),Ver_b_8%m(k),Ver_c_8%m(k),Cstv_pref_8, &
                                              Shear,u,v,w,t,tv,phis,orls,ps,rho,q,Set_topo_L)

                  F_v(i,j,k) = v

               end do

            else

               do i = i0,in

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  if (.NOT.Set_topo_L) phis = F_topo(i,j)
                  if (.NOT.Set_topo_L) orls = F_orls(i,j)

                  call test2_schaer_mountain (lon,lat,p,z,zcoords,Ver_z_8%m(k),Ver_a_8%m(k),Ver_b_8%m(k),Ver_c_8%m(k),Cstv_pref_8, &
                                              Shear,utt_8,vtt_8,w,t,tv,phis,orls,ps,rho,q,Set_topo_L)

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

         do j = j0,jn

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = i0,inu

                  lon = geomh_xu_8(i)

                  if (.NOT.Set_topo_L) phis = F_topo(i,j) !ZERO
                  if (.NOT.Set_topo_L) orls = F_orls(i,j) !ZERO

                  call test2_schaer_mountain (lon,lat,p,z,zcoords,Ver_z_8%m(k),Ver_a_8%m(k),Ver_b_8%m(k),Ver_c_8%m(k),Cstv_pref_8, &
                                              Shear,u,v,w,t,tv,phis,orls,ps,rho,q,Set_topo_L)

                  F_u(i,j,k) = u

               end do

            else

               do i = i0,inu

                  x_a_8 = geomh_xu_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  if (.NOT.Set_topo_L) phis = F_topo(i,j) !ZERO
                  if (.NOT.Set_topo_L) orls = F_orls(i,j) !ZERO

                  call test2_schaer_mountain (lon,lat,p,z,zcoords,Ver_z_8%m(k),Ver_a_8%m(k),Ver_b_8%m(k),Ver_c_8%m(k),Cstv_pref_8, &
                                              Shear,utt_8,vtt_8,w,t,tv,phis,orls,ps,rho,q,Set_topo_L)

                  u = s_8(1,1)*utt_8 + s_8(1,2)*vtt_8

                  F_u(i,j,k) = u

               end do

            end if

         end do

      end do

      !Initial conditions: V True
      !--------------------------
      do k = 1,Nk

         do j = j0,jnv

            lat   = geomh_yv_8(j)
            y_a_8 = geomh_yv_8(j)

            if (Ptopo_couleur == 0) then

               do i = i0,in

                  lon = geomh_x_8(i)

                  if (.NOT.Set_topo_L) phis = F_topo(i,j) !ZERO
                  if (.NOT.Set_topo_L) orls = F_orls(i,j) !ZERO

                  call test2_schaer_mountain (lon,lat,p,z,zcoords,Ver_z_8%m(k),Ver_a_8%m(k),Ver_b_8%m(k),Ver_c_8%m(k),Cstv_pref_8, &
                                              Shear,u,v,w,t,tv,phis,orls,ps,rho,q,Set_topo_L)

                  F_v(i,j,k) = v

               end do

            else

               do i = i0,in

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  if (.NOT.Set_topo_L) phis = F_topo(i,j) !ZERO
                  if (.NOT.Set_topo_L) orls = F_orls(i,j) !ZERO

                  call test2_schaer_mountain (lon,lat,p,z,zcoords,Ver_z_8%m(k),Ver_a_8%m(k),Ver_b_8%m(k),Ver_c_8%m(k),Cstv_pref_8, &
                                              Shear,utt_8,vtt_8,w,t,tv,phis,orls,ps,rho,q,Set_topo_L)

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
      /,'USE INITIAL CONDITIONS FOR MOUNTAIN WAVES OVER A SCHAER-TYPE MOUNTAIN ON A SMALL PLANET : (S/R DCMIP_SCHAER_MOUNTAIN)',   &
      /,'=====================================================================================================================',/, &
        ' Shear wind = ',I1                                                                                                    ,   &
        ' Set Topo   = ',L2                                                                                                    ,   &
      /,'=====================================================================================================================',/,/)

      end subroutine dcmip_Schaer_mountain
