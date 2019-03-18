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

!**s/r dcmip_tracers13_transport - Setup for Horizontal advection of thin cloud-like tracers
!                                  in the presence of orography (DCMIP 2012)

      subroutine dcmip_tracers13_transport (F_u,F_v,F_zd,F_tv,F_q,F_topo,F_s,F_q1,F_q2,F_q3,F_q4, &
                                            Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L)

      use dcmip_2012_init_1_2_3
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

      real F_u    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_v    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_zd   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_tv   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Virtual temperature
           F_q    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_s    (Mminx:Mmaxx,Mminy:Mmaxy)   , &
           F_topo (Mminx:Mmaxx,Mminy:Mmaxy)   , &
           F_q1   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_q2   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_q3   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_q4   (Mminx:Mmaxx,Mminy:Mmaxy,Nk)

      logical F_stag_L ! Staggered uv if .T. / Scalar uv if .F.

      !object
      !======================================================================================================
      !   Setup for Horizontal advection of thin cloud-like tracers in the presence of orography (DCMIP 2012)
      !======================================================================================================

      !-----------------------------------------------------------------------------------------------------------

      integer i,j,k

      real(8) x_a_8,y_a_8,utt_8,vtt_8,s_8(2,2),rlon_8

      real(8)  :: &
                  lon,  &    ! Longitude (radians)
                  lat,  &    ! Latitude (radians)
                  z, &       ! Height (m)
               !!!hyam, &    ! A coefficient for hybrid-eta coordinate, at model level midpoint
               !!!hybm, &    ! B coefficient for hybrid-eta coordinate, at model level midpoint
                  gc         ! bar{z} for Gal-Chen coordinate

      logical  :: hybrid_eta ! flag to indicate whether the hybrid sigma-p (eta) coordinate is used
                             ! if set to .true., then the pressure will be computed via the
                          !!!!    hybrid coefficients hyam and hybm, they need to be initialized
                             !    hybrid coefficients GEM Ver_a and Ver_b, they need to be initialized
                             ! if set to .false. (for pressure-based models): the pressure is already pre-computed
                             !    and is an input value for this routine
                             ! for height-based models: pressure will always be computed based on the height and
                             !    hybrid_eta is not used

      real(8)  :: p          ! Pressure  (Pa)

      integer  :: zcoords    ! 0 or 1 see below
      integer  :: cfv        ! 0, 1 or 2 see below

      real(8)  :: &
                  u,       & ! Zonal wind (m s^-1)
                  v,       & ! Meridional wind (m s^-1)
                  w,       & ! Vertical Velocity (m s^-1)
                  zd,      & ! Zdot GEM
                  t,       & ! Temperature (K)
                  tv,      & ! Virtual Temperature (K)
                  phis,    & ! Surface Geopotential (m^2 s^-2)
                  ps,      & ! Surface Pressure (Pa)
                  rho,     & ! density (kg m^-3)
                  q,       & ! Specific Humidity (kg/kg)
                  q1,      & ! Tracer q1 (kg/kg)
                  q2,      & ! Tracer q2 (kg/kg)
                  q3,      & ! Tracer q3 (kg/kg)
                  q4         ! Tracer q4 (kg/kg)

      ! if zcoords = 1, then we use z and output p
      ! if zcoords = 0, then we use p

      ! if cfv = 0 we assume that our horizontal velocities are not coordinate following
   !!!! if cfv = 1 then our velocities follow hybrid eta coordinates and we need to specify w
      ! if cfv = 1 then our velocities follow hybrid eta GEM coordinates and we need to specify w
      ! if cfv = 2 then our velocities follow Gal-Chen coordinates and we need to specify w

   !!!! In hybrid-eta coords    : p = hyam p0 + hybm ps
      ! In hybrid-eta GEM coords: p = exp(Ver_a p0 + Ver_b ps)
      ! In Gal-Chen coords: z = zs + (gc/ztop)*(ztop - zs)

      ! if other orography-following coordinates are used, the w wind needs to be newly derived for them

      !-----------------------------------------------------------------------------------------------------------

      if (Lun_out > 0) write (Lun_out,1000)

      zcoords    = 0
      cfv        = 1
      hybrid_eta = .TRUE. ! as in GEM

      !Initial conditions: T,ZD,Q,Q1,Q2,Q3,Q4,S,TOPO
      !---------------------------------------------
      do k = 1,Nk

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call test1_advection_orography (lon,lat,p,z,zcoords,cfv,Ver_a_8%t(k),Ver_b_8%t(k),Cstv_pref_8,hybrid_eta,gc, &
                                                  u,v,w,zd,t,tv,phis,ps,rho,q,q1,q2,q3,q4)

                  F_tv  (i,j,k) = tv
                  F_q   (i,j,k) = q
                  F_s   (i,j)   = log(ps/Cstv_pref_8)
                  F_topo(i,j)   = phis
                  F_zd  (i,j,k) = zd

                  !Tracers
                  !-------
                  F_q1  (i,j,k) = q1
                  F_q2  (i,j,k) = q2
                  F_q3  (i,j,k) = q3
                  F_q4  (i,j,k) = q4

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call test1_advection_orography (lon,lat,p,z,zcoords,cfv,Ver_a_8%t(k),Ver_b_8%t(k),Cstv_pref_8,hybrid_eta,gc, &
                                                  utt_8,vtt_8,w,zd,t,tv,phis,ps,rho,q,q1,q2,q3,q4)

                  F_tv  (i,j,k) = tv
                  F_q   (i,j,k) = q
                  F_s   (i,j)   = log(ps/Cstv_pref_8)
                  F_topo(i,j)   = phis
                  F_zd  (i,j,k) = zd

                  !Tracers
                  !-------
                  F_q1  (i,j,k) = q1
                  F_q2  (i,j,k) = q2
                  F_q3  (i,j,k) = q3
                  F_q4  (i,j,k) = q4

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

                  call test1_advection_orography (lon,lat,p,z,zcoords,cfv,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,hybrid_eta,gc, &
                                                  u,v,w,zd,t,tv,phis,ps,rho,q,q1,q2,q3,q4)

                  F_u(i,j,k) = u

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call test1_advection_orography (lon,lat,p,z,zcoords,cfv,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,hybrid_eta,gc, &
                                                  utt_8,vtt_8,w,zd,t,tv,phis,ps,rho,q,q1,q2,q3,q4)

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

                  call test1_advection_orography (lon,lat,p,z,zcoords,cfv,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,hybrid_eta,gc, &
                                                  u,v,w,zd,t,tv,phis,ps,rho,q,q1,q2,q3,q4)

                  F_v(i,j,k) = v

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call test1_advection_orography (lon,lat,p,z,zcoords,cfv,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,hybrid_eta,gc, &
                                                  utt_8,vtt_8,w,zd,t,tv,phis,ps,rho,q,q1,q2,q3,q4)

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

                  call test1_advection_orography (lon,lat,p,z,zcoords,cfv,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,hybrid_eta,gc, &
                                                  u,v,w,zd,t,tv,phis,ps,rho,q,q1,q2,q3,q4)

                  F_u(i,j,k) = u

               end do

            else

               do i = 1,l_niu

                  x_a_8 = geomh_xu_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call test1_advection_orography (lon,lat,p,z,zcoords,cfv,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,hybrid_eta,gc, &
                                                  utt_8,vtt_8,w,zd,t,tv,phis,ps,rho,q,q1,q2,q3,q4)

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

                  call test1_advection_orography (lon,lat,p,z,zcoords,cfv,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,hybrid_eta,gc, &
                                                  u,v,w,zd,t,tv,phis,ps,rho,q,q1,q2,q3,q4)

                  F_v(i,j,k) = v

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.d0)

                  call test1_advection_orography (lon,lat,p,z,zcoords,cfv,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,hybrid_eta,gc, &
                                                  utt_8,vtt_8,w,zd,t,tv,phis,ps,rho,q,q1,q2,q3,q4)

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
      /,'PRESCRIBED CONDITIONS For Horizontal advection of thin cloud-like tracers in the presence of orography (DCMIP 2012 T13)',   &
      /,'=======================================================================================================================',/,/)

      end subroutine dcmip_tracers13_transport
