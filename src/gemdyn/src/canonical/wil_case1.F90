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

!**s/r wil_case1 - To setup Williamson/Lauritzen/Nair's cases (TRACER)

      subroutine wil_case1 (F_tr,F_minx,F_maxx,F_miny,F_maxy,F_nk,The_tracer,F_istep)

      use cstv
      use gem_options
      use glb_ld
      use lun
      use ptopo
      use tdpack
      use wil_options

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer, intent(in) :: F_minx,F_maxx,F_miny,F_maxy,F_nk,The_tracer,F_istep
      real, intent(out) :: F_tr(F_minx:F_maxx,F_miny:F_maxy,F_nk)

      !authors
      !     Abdessamad Qaddouri and Vivian Lee
      !
      !revision
      ! v5_00 - Tanguay M. - New Testcases
      !
      !object
      !=======================================================================
      ! Williamson_NAIR    |=0=Solid body rotation of a cosine bell          |
      !                    |   Williamson et al.,1992,JCP,102,211-224        |
      !                    |=1=Deformational Non-divergent winds             |
      !                    |   Lauritzen et al.,2012,GMD,5,887-901           |
      !                    |=2=Deformational divergent winds                 |
      !                    |   Lauritzen et al.,2012,GMD,5,887-901           |
      !                    |=3=Deformational Flow for Circular vortex        |
      !                    |   Nair and Machenhauer,2002,MWR,130,649-667     |
      !=======================================================================

      integer :: i,j,k,g_i0,g_in,g_j0,g_jn,i0,in,j0,jn,zlist
      real(kind=REAL64):: phi0_8,rlon_8,rlat_8,sint_8,cost_8,        &
                          s_8(2,2),x_a_8,y_a_8,sinl_8,cosl_8,        &
                          sinl1_8,cosl1_8,sinl2_8,cosl2_8,           &
                          sint1_8,cost1_8,sint2_8,cost2_8,           &
                          rlon1_8,rlat1_8,rlon2_8,rlat2_8,           &
                          x_8,y_8,z_8,x1_8,y1_8,z1_8,x2_8,y2_8,z2_8, &
                          radius_8,dist_8,rrlata_8,rrlona_8,         &
                          rlat0_8,rlon0_8,dist1_8,dist2_8,           &
                          rlonr_8,rlatr_8,rho_8,vt_8,wil_omega_8,    &
                          rlona_8,rlata_8,rlon_r_8,rlat_r_8,         &
                          rlon_nr_8,rlat_nr_8,clon0_8,clat0_8,       &
                          time_8,gamma_8,phi1_8,phi2_8

      real :: picll(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy), &
              trloc(F_minx:F_maxx,F_miny:F_maxy)

      logical, save :: done_L=.false.
!
!---------------------------------------------------------------------
!
      g_i0= 1-G_halox ; g_in= G_ni+G_halox
      g_j0= 1-G_haloy ; g_jn= G_nj+G_haloy

      i0= 1-G_halox ; in= l_ni+G_halox
      j0= 1-G_haloy ; jn= l_nj+G_haloy

      !----------------------------------------------------
      !Williamson et al.,1992,JCP,102,211-224 [Cosine bell]
      !----------------------------------------------------
      if (The_tracer==0) then

         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) ''
         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) '--------------------------------------------------------------'
         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON CASE1 (Tracer): Cosine Bell, Williamson et al.,1992'
         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) '--------------------------------------------------------------'

         !--------------------------------------------------------------------------
         !Williamson_lon_pole_r_8 and Williamson_lat_pole_r_8 are Longitude,Latitude
         !in the no rotated system of the north pole P' of the rotated system
         !--------------------------------------------------------------------------

         Williamson_lon_pole_r_8 = Williamson_rlon0
         Williamson_lat_pole_r_8 = 0.5 * pi_8 - Williamson_alpha

         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON_ALPHA      = ',Williamson_alpha
         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON_LON_POLE_R = ',Williamson_lon_pole_r_8
         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON_LAT_POLE_R = ',Williamson_lat_pole_r_8
         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) '--------------------------------------------------------------'

         !------------------------------------------------------------------------------------------
         !RLON and RLAT are Longitude,Latitude of the center of Cosine bell in the no rotated system
         !------------------------------------------------------------------------------------------
         clon0_8 = Williamson_clon0
         clat0_8 = Williamson_clat0

         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON_CENTER_LON = ',clon0_8
         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON CENTER_LAT = ',clat0_8

         !------------------------------------------
         !Prescribed velocity of Solid body rotation
         !------------------------------------------
         Williamson_ubar_8 = (2.0*pi_8*rayt_8)/Williamson_period

         wil_omega_8 = Williamson_ubar_8/rayt_8

         time_8 = F_istep*cstv_dt_8

         phi0_8 = Williamson_phi0

         radius_8 = rayt_8/Williamson_radius

         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON_PHI0       = ',phi0_8
         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON_RADIUS     = ',radius_8
         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) '--------------------------------------------------------------'

         !Find associated LON LAT rotated at TIME 0 of CENTER
         !---------------------------------------------------
         rlon_nr_8 = clon0_8
         rlat_nr_8 = clat0_8

         call wil_rotate (rlon_nr_8,rlat_nr_8,rlon_r_8,rlat_r_8,0)

         !Advection of LON LAT rotated from TIME 0 to TIME T
         !--------------------------------------------------
         rlona_8 = rlon_r_8 + wil_omega_8*time_8
         rlata_8 = rlat_r_8

         !Find associated LON LAT not rotated at TIME T
         !---------------------------------------------
         call wil_rotate (rlona_8,rlata_8,rrlona_8,rrlata_8,1)

         !Compute tracer for YIN
         !----------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jn

               rlat_8 = G_yg_8(j)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i)

                  dist_8 = rayt_8*acos(sin(rrlata_8)*sin(rlat_8) + cos(rrlata_8) &
                           *cos(rlat_8)*cos(rlon_8-rrlona_8))

                  if (dist_8 <= radius_8) then
                      picll(i,j) = phi0_8/2.0*(1.0 + cos(pi_8*dist_8/radius_8))
                  else
                      picll(i,j) = 0.0
                  end if

               end do

            end do

         !Compute tracer for YAN
         !----------------------
         else

            do j=g_j0,g_jn

               do i=g_i0,g_in

                  x_a_8  = G_xg_8(i)-acos(-1.d0)
                  y_a_8  = G_yg_8(j)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0)

                  dist_8 = rayt_8*acos(sin(rrlata_8)*sin(rlat_8) + cos(rrlata_8) &
                           *cos(rlat_8)*cos(rlon_8-rrlona_8))

                  if (dist_8 <= radius_8) then
                      picll(i,j) = phi0_8/2.0*(1.0 + cos(pi_8*dist_8/radius_8))
                  else
                      picll(i,j) = 0.0
                  end if

               end do

            end do

         end if

      !--------------------------------------------------------------------------
      !Lauritzen et al.,2012,GMD,5,887-901 [Cosine bells/Correlated Cosine bells]
      !--------------------------------------------------------------------------
      else if (The_tracer==2.or.The_tracer==4) then

         if (.NOT.done_L.and.The_tracer==2.and.Lun_out>0) &
             write(Lun_out,*) 'WILLIAMSON CASE1 (Tracer): Cosine bells, Lauritzen et al.,2012'

         if (.NOT.done_L.and.The_tracer==4.and.Lun_out>0) then
             write(Lun_out,*) 'WILLIAMSON CASE1 (Tracer): Correlated Cosine bells, Lauritzen et al.,2012'
             write(Lun_out,*) '--------------------------------------------------------------'
         end if

         rlon1_8 = 5.*pi_8/6.
         rlat1_8 = 0.0

         rlon2_8 = 7.*pi_8/6.
         rlat2_8 = 0.0

         phi0_8   = 1.0
         radius_8 = rayt_8/2.

         !Compute tracer for YIN
         !----------------------
         if (ptopo_couleur==0) then

            do j=g_j0,g_jn

               rlat_8 = G_yg_8(j)

               cost_8 = cos(rlat_8)
               sint_8 = sin(rlat_8)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i)

                  sinl_8 = sin(rlon_8)
                  cosl_8 = cos(rlon_8)

                  dist1_8 = rayt_8*acos(sin(rlat1_8)*sin(rlat_8) + cos(rlat1_8) &
                            *cos(rlat_8)*cos(rlon_8-rlon1_8))
                  dist2_8 = rayt_8*acos(sin(rlat2_8)*sin(rlat_8) + cos(rlat2_8) &
                            *cos(rlat_8)*cos(rlon_8-rlon2_8))

                  picll(i,j) = 0.1

                  if (dist1_8 < radius_8) picll(i,j) = 0.1 + 0.9*phi0_8/2.0*(1.0 + cos(pi_8*dist1_8/radius_8))
                  if (dist2_8 < radius_8) picll(i,j) = 0.1 + 0.9*phi0_8/2.0*(1.0 + cos(pi_8*dist2_8/radius_8))

                  if (The_tracer==4) picll(i,j) = -0.8*picll(i,j)**2 + 0.9

               end do

            end do

         !Compute tracer for YAN
         !----------------------
         else

            do j=g_j0,g_jn

               do i=g_i0,g_in

                  x_a_8  = G_xg_8(i)-acos(-1.d0)
                  y_a_8  = G_yg_8(j)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0)

                  cost_8 = cos(rlat_8)
                  sint_8 = sin(rlat_8)

                  sinl_8 = sin(rlon_8)
                  cosl_8 = cos(rlon_8)

                  dist1_8 = rayt_8*acos(sin(rlat1_8)*sin(rlat_8) + cos(rlat1_8) &
                           *cos(rlat_8)*cos(rlon_8-rlon1_8))
                  dist2_8 = rayt_8*acos(sin(rlat2_8)*sin(rlat_8) + cos(rlat2_8) &
                           *cos(rlat_8)*cos(rlon_8-rlon2_8))

                  picll(i,j) = 0.1

                  if (dist1_8 < radius_8) picll(i,j) = 0.1 + 0.9*phi0_8/2.0*(1.0 + cos(pi_8*dist1_8/radius_8))
                  if (dist2_8 < radius_8) picll(i,j) = 0.1 + 0.9*phi0_8/2.0*(1.0 + cos(pi_8*dist2_8/radius_8))

                  if (The_tracer==4) picll(i,j) = -0.8*picll(i,j)**2 + 0.9

               end do

            end do

         end if

      !----------------------------------------------------
      !Lauritzen et al.,2012,GMD,5,887-901 [Gaussian hills]
      !----------------------------------------------------
      else if (The_tracer==1) then

         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) ''
         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) '--------------------------------------------------------------'
         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON CASE1 (Tracer): Gaussian hills, Lauritzen et al.,2012'

         rlon1_8 = 5.*pi_8/6.
         rlat1_8 = 0.0

         rlon2_8 = 7.*pi_8/6.
         rlat2_8 = 0.0

         cost1_8 = cos(rlat1_8)
         sint1_8 = sin(rlat1_8)
         cost2_8 = cos(rlat2_8)
         sint2_8 = sin(rlat2_8)

         sinl1_8 = sin(rlon1_8)
         cosl1_8 = cos(rlon1_8)
         sinl2_8 = sin(rlon2_8)
         cosl2_8 = cos(rlon2_8)

         !Compute tracer for YIN
         !----------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jn

               rlat_8 = G_yg_8(j)

               cost_8 = cos(rlat_8)
               sint_8 = sin(rlat_8)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i)

                  sinl_8 = sin(rlon_8)
                  cosl_8 = cos(rlon_8)

                  x_8  = cost_8*cosl_8
                  y_8  = cost_8*sinl_8
                  z_8  = sint_8
                  x1_8 = cost1_8*cosl1_8
                  y1_8 = cost1_8*sinl1_8
                  z1_8 = sint1_8
                  x2_8 = cost2_8*cosl2_8
                  y2_8 = cost2_8*sinl2_8
                  z2_8 = sint2_8

                  phi1_8 = 0.95 * exp(-5. * ( (x_8-x1_8)**2 + (y_8-y1_8)**2 + (z_8-z1_8)**2 ) )
                  phi2_8 = 0.95 * exp(-5. * ( (x_8-x2_8)**2 + (y_8-y2_8)**2 + (z_8-z2_8)**2 ) )

                  picll(i,j) = phi1_8 + phi2_8

               end do

            end do

         !Compute tracer for YAN
         !----------------------
         else

            do j=g_j0,g_jn

               do i=g_i0,g_in

                  x_a_8 = G_xg_8(i)-acos(-1.d0)
                  y_a_8 = G_yg_8(j)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0)

                  cost_8 = cos(rlat_8)
                  sint_8 = sin(rlat_8)

                  sinl_8 = sin(rlon_8)
                  cosl_8 = cos(rlon_8)

                  x_8  = cost_8*cosl_8
                  y_8  = cost_8*sinl_8
                  z_8  = sint_8
                  x1_8 = cost1_8*cosl1_8
                  y1_8 = cost1_8*sinl1_8
                  z1_8 = sint1_8
                  x2_8 = cost2_8*cosl2_8
                  y2_8 = cost2_8*sinl2_8
                  z2_8 = sint2_8

                  phi1_8 = 0.95 * exp(-5. * ( (x_8-x1_8)**2 + (y_8-y1_8)**2 + (z_8-z1_8)**2 ) )
                  phi2_8 = 0.95 * exp(-5. * ( (x_8-x2_8)**2 + (y_8-y2_8)**2 + (z_8-z2_8)**2 ) )

                  picll(i,j) = phi1_8 + phi2_8

               end do

            end do

         end if

      !-------------------------------------------------------
      !Lauritzen et al.,2012,GMD,5,887-901 [Slotted cylinders]
      !-------------------------------------------------------
      else if (The_tracer==3) then

         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON CASE1 (Tracer): Slotted cylinders, Lauritzen et al.,2012'

         rlon1_8 = 5.*pi_8/6.
         rlat1_8 = 0.0

         rlon2_8 = 7.*pi_8/6.
         rlat2_8 = 0.0

         radius_8 = rayt_8/2.

         !Compute tracer for YIN
         !----------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jn

               rlat_8 = G_yg_8(j)

               cost_8 = cos(rlat_8)
               sint_8 = sin(rlat_8)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i)

                  sinl_8 = sin(rlon_8)
                  cosl_8 = cos(rlon_8)

                  dist1_8 = rayt_8*acos(sin(rlat1_8)*sin(rlat_8) + cos(rlat1_8) &
                            *cos(rlat_8)*cos(rlon_8-rlon1_8))
                  dist2_8 = rayt_8*acos(sin(rlat2_8)*sin(rlat_8) + cos(rlat2_8) &
                            *cos(rlat_8)*cos(rlon_8-rlon2_8))

                  picll(i,j) = 0.1

                  if (   dist1_8 <= radius_8 .and. abs(rlon_8-rlon1_8) >= radius_8/(6.*rayt_8)) then
                      picll(i,j)  = 1.0
                  else if(dist1_8 <= radius_8 .and. abs(rlon_8-rlon1_8) <  radius_8/(6.*rayt_8) .and.  &
                         rlat_8-rlat1_8 < (-(5./12.)*radius_8)/rayt_8) then
                      picll(i,j)  = 1.0
                  end if

                  if (   dist2_8 <= radius_8 .and. abs(rlon_8-rlon2_8) >= radius_8/(6.*rayt_8)) then
                      picll(i,j)  = 1.0
                  else if(dist2_8 <= radius_8 .and. abs(rlon_8-rlon2_8) <  radius_8/(6.*rayt_8) .and.  &
                         rlat_8-rlat2_8 > (+(5./12.)*radius_8)/rayt_8) then
                      picll(i,j)  = 1.0
                  end if

               end do

            end do

         !Compute tracer for YAN
         !----------------------
         else

            do j=g_j0,g_jn

               do i=g_i0,g_in

                  x_a_8 = G_xg_8(i)-acos(-1.d0)
                  y_a_8 = G_yg_8(j)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0)

                  cost_8 = cos(rlat_8)
                  sint_8 = sin(rlat_8)

                  sinl_8 = sin(rlon_8)
                  cosl_8 = cos(rlon_8)

                  dist1_8 = rayt_8*acos(sin(rlat1_8)*sin(rlat_8) + cos(rlat1_8) &
                            *cos(rlat_8)*cos(rlon_8-rlon1_8))
                  dist2_8 = rayt_8*acos(sin(rlat2_8)*sin(rlat_8) + cos(rlat2_8) &
                            *cos(rlat_8)*cos(rlon_8-rlon2_8))

                  picll(i,j) = 0.1

                  if (   dist1_8 <= radius_8 .and. abs(rlon_8-rlon1_8) >= radius_8/(6.*rayt_8)) then
                      picll(i,j)  = 1.0
                  else if(dist1_8 <= radius_8 .and. abs(rlon_8-rlon1_8) <  radius_8/(6.*rayt_8) .and.  &
                         rlat_8-rlat1_8 < (-(5./12.)*radius_8)/rayt_8) then
                      picll(i,j)  = 1.0
                  end if

                  if (   dist2_8 <= radius_8 .and. abs(rlon_8-rlon2_8) >= radius_8/(6.*rayt_8)) then
                      picll(i,j)  = 1.0
                  else if(dist2_8 <= radius_8 .and. abs(rlon_8-rlon2_8) <  radius_8/(6.*rayt_8) .and.  &
                         rlat_8-rlat2_8>(+(5./12.)*radius_8)/rayt_8) then
                      picll(i,j)  = 1.0
                  end if

               end do

            end do

         end if

      !-----------------------------------------------------------
      !Nair and Machenhauer,2002,MWR,130,649-667 [Circular vortex]
      !-----------------------------------------------------------
      else if (The_tracer==5) then

         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) ''
         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) '--------------------------------------------------------------'
         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON CASE1 (Tracer): Circular vortex, Nair and Machenhauer,2002'
         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) '--------------------------------------------------------------'

         Williamson_v0_8    = -(2.0*pi_8*rayt_8)/(12.0*24.0*3600.0)
         Williamson_rho_i_8 = 3.

         gamma_8 = 5.

         rlon0_8 = Williamson_rlon0
         rlat0_8 = Williamson_rlat0

         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON_LON_POLE_R = ',rlon0_8
         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON_LAT_POLE_R = ',rlat0_8
         if (.NOT.done_L.and.Lun_out>0) write(Lun_out,*) '--------------------------------------------------------------'

         time_8 = F_istep*cstv_dt_8

         !Compute tracer for YIN
         !----------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jn

               rlat_8 = G_yg_8(j)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i)

                  rlonr_8 = atan2( cos(rlat_8)*sin(rlon_8-rlon0_8), &
                                   cos(rlat_8)*sin(rlat0_8)*cos(rlon_8-rlon0_8)-cos(rlat0_8)*sin(rlat_8) )

                  if (rlonr_8 < 0.0d0) rlonr_8 = rlonr_8 + 2.*pi_8

                  rlatr_8 = asin( sin(rlat_8)*sin(rlat0_8) + cos(rlat_8)*cos(rlat0_8)*cos(rlon_8-rlon0_8) )

                  rho_8 = Williamson_rho_i_8*cos(rlatr_8)

                  vt_8 = Williamson_v0_8*(3.*sqrt(3.)/2) * (1.0d0/cosh(rho_8))**2 * tanh(rho_8)

                  wil_omega_8 = 0.0d0
                  if (rho_8 /= 0.0d0) wil_omega_8 = vt_8/(rayt_8*rho_8)

                  picll(i,j) = 1. - tanh( (rho_8/gamma_8) * sin(rlonr_8 - wil_omega_8*time_8) )

               end do

            end do

         !Compute tracer for YAN
         !----------------------
         else

            do j=g_j0,g_jn

               do i=g_i0,g_in

                  x_a_8 = G_xg_8(i)-acos(-1.d0)
                  y_a_8 = G_yg_8(j)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0)

                  rlonr_8 = atan2( cos(rlat_8)*sin(rlon_8-rlon0_8), &
                                   cos(rlat_8)*sin(rlat0_8)*cos(rlon_8-rlon0_8)-cos(rlat0_8)*sin(rlat_8) )

                  if (rlonr_8 < 0.0d0) rlonr_8 = rlonr_8 + 2.*pi_8

                  rlatr_8 = asin( sin(rlat_8)*sin(rlat0_8) + cos(rlat_8)*cos(rlat0_8)*cos(rlon_8-rlon0_8) )

                  rho_8 = Williamson_rho_i_8*cos(rlatr_8)

                  vt_8 = Williamson_v0_8*(3.*sqrt(3.)/2) * (1.0d0/cosh(rho_8))**2 * tanh(rho_8)

                  wil_omega_8 = 0.0d0

                  if (rho_8 /= 0.0d0) wil_omega_8 = vt_8/(rayt_8*rho_8)

                  picll(i,j) = 1. - tanh( (rho_8/gamma_8) * sin(rlonr_8 - wil_omega_8*time_8) )

               end do

            end do

         end if

      end if

      trloc = 0.

      zlist = 1

      call glbdist_os (picll, trloc,&
                       F_minx,F_maxx,F_miny,F_maxy,1,&
                       G_ni+G_halox,G_nj+G_haloy,zlist,1,1.0d0,0.d0)

      do k=1,F_nk
         F_tr(i0:in,j0:jn,k) = trloc(i0:in,j0:jn)
      end do

      if ( Williamson_Nair==0.or.Williamson_Nair==3                   ) done_L = .TRUE.
      if ((Williamson_Nair==1.or.Williamson_Nair==2).and.The_tracer==4) done_L = .TRUE.
!
!---------------------------------------------------------------------
!
      return
      end

