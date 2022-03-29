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

!**s/r wil_uvcase1 - To setup Williamson/Lauritzen/Nair's cases (WINDS)

      subroutine wil_uvcase1 (F_u,F_v,F_minx,F_maxx,F_miny,F_maxy,F_nk,F_stag_L,F_istep)

      use wil_options
      use gem_options
      use tdpack
      use glb_ld
      use cstv
      use lun
      use ptopo

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer F_minx,F_maxx,F_miny,F_maxy,F_nk,F_istep
      real    F_u(F_minx:F_maxx,F_miny:F_maxy,F_nk),F_v(F_minx:F_maxx,F_miny:F_maxy,F_nk)

      logical F_stag_L ! Staggered uv if .T. / Scalar uv if .F.

      !authors
      !     Abdessamad Qaddouri & Vivian Lee
      !
      !revision
      ! v5_00 - Tanguay M. - New Testcases
      !
      !object
      !----------------------------------------------------------------------|
      ! Williamson_NAIR    |=0=Solid body rotation of a cosine bell          |
      !                    |   Williamson et al.,1992,JCP,102,211-224        |
      !                    |=1=Deformational Non-divergent winds             |
      !                    |   Lauritzen et al.,2012,GMD,5,887-901           |
      !                    |=2=Deformational divergent winds                 |
      !                    |   Lauritzen et al.,2012,GMD,5,887-901           |
      !                    |=3=Deformational Flow for Circular vortex        |
      !                    |   Nair and Machenhauer,2002,MWR,130,649-667     |
      !----------------------------------------------------------------------|

      !-------------------------------------------------------

      integer i,j,k,g_i0,g_in,g_j0,g_jn,g_inu,g_jnv,i0,in,j0,jn,inu,jnv,zlist
      real(kind=REAL64)  sina_8,cosa_8,                              &
              ui_u_8(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy), &
              ui_v_8(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy), &
              vi_u_8(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy), &
              vi_v_8(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy), &
              rlon_8,rlat_8,sint_8,cost_8,time_frac_8,               &
              s_8(2,2),x_a_8,y_a_8,sinl_8,cosl_8,                    &
              rlon0_8,rlat0_8,rlonr_8,rlatr_8,rho_8,vt_8,            &
              wil_omega_8,theta_8,lambda_8,                          &
              xgu_8(1-G_halox:G_niu+G_halox),                        &
              ygv_8(1-G_haloy:G_njv+G_haloy)

      real    uloc(F_minx:F_maxx,F_miny:F_maxy),                     &
              vloc(F_minx:F_maxx,F_miny:F_maxy),                     &
              uicll(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy),  &
              vicll(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy)

      real(kind=REAL64), parameter ::  DIX_8 = 10.0d0
      real(kind=REAL64), parameter :: FIVE_8 =  5.0d0
      real(kind=REAL64), parameter ::  TWO_8 =  2.0d0

      !-------------------------------------------------------

      g_i0= 1-G_halox ; g_in= G_ni+G_halox ; g_inu= G_niu+G_halox
      g_j0= 1-G_haloy ; g_jn= G_nj+G_haloy ; g_jnv= G_njv+G_haloy

      i0= 1-G_halox ; in= l_ni+G_halox ; inu= l_niu+G_halox
      j0= 1-G_haloy ; jn= l_nj+G_haloy ; jnv= l_njv+G_haloy

      !U grid
      !------
      do i=g_i0,g_inu
         xgu_8(i) = (G_xg_8(i+1)+G_xg_8(i))*.5
      end do

      !V grid
      !------
      do j=g_j0,g_jnv
         ygv_8(j) = (G_yg_8(j+1)+G_yg_8(j))*.5
      end do

      !------------------------------------------------------------
      !Williamson et al.,1992,JCP,102,211-224 [Solid body rotation]
      !------------------------------------------------------------
      if (Williamson_NAIR==0) then

         if (F_istep==0.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON CASE1 (Winds) : Solid body rotation, Williamson et al.,1992'
         if (F_istep==0.and.Lun_out>0) write(Lun_out,*) '--------------------------------------------------------------'
         if (F_istep==0.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON_UBAR (m/s) = ',Williamson_ubar_8

         theta_8  = Williamson_lat_pole_r_8
         lambda_8 = Williamson_lon_pole_r_8

         sina_8 = sin(theta_8)
         cosa_8 = cos(theta_8)

         !######################
         if (.not.F_stag_L) then
         !######################

         uicll=0.; vicll=0.

         !Compute U vector for YIN
         !------------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jn

               rlat_8 = G_yg_8(j)

               sint_8 = sin(rlat_8)
               cost_8 = cos(rlat_8)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i)

                  sinl_8 = sin(rlon_8-lambda_8)
                  cosl_8 = cos(rlon_8-lambda_8)

                  uicll(i,j) = Williamson_ubar_8*(cost_8*sina_8 - cosl_8*sint_8*cosa_8)

               end do

            end do

         !Compute U vector for YAN
         !------------------------
         else

            do j=g_j0,g_jn

               y_a_8 = G_yg_8(j)

               do i=g_i0,g_in

                  x_a_8 = G_xg_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0)

                  sint_8 = sin(rlat_8)
                  cost_8 = cos(rlat_8)

                  sinl_8 = sin(rlon_8-lambda_8)
                  cosl_8 = cos(rlon_8-lambda_8)

                  ui_u_8(i,j) = Williamson_ubar_8*(cost_8*sina_8 - cosl_8*sint_8*cosa_8)

                  vi_u_8(i,j) = Williamson_ubar_8*sinl_8*cosa_8

                   uicll(i,j) = s_8(1,1)*ui_u_8(i,j) + s_8(1,2)*vi_u_8(i,j)

               end do

            end do

         end if

         !Compute V vector for YIN
         !------------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jn

               rlat_8 = G_yg_8(j)

               sint_8 = sin(rlat_8)
               cost_8 = cos(rlat_8)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i)

                  sinl_8 = sin(rlon_8-lambda_8)
                  cosl_8 = cos(rlon_8-lambda_8)

                  vicll(i,j) = Williamson_ubar_8*sinl_8*cosa_8

               end do

            end do

         !Compute V vector for YAN
         !------------------------
         else

            do j=g_j0,g_jn

               y_a_8 = G_yg_8(j)

               do i=g_i0,g_in

                  x_a_8 = G_xg_8(i)-acos(-1.d0)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0)

                  sint_8 = sin(rlat_8)
                  cost_8 = cos(rlat_8)

                  sinl_8 = sin(rlon_8-lambda_8)
                  cosl_8 = cos(rlon_8-lambda_8)

                  ui_v_8(i,j) = Williamson_ubar_8*(cost_8*sina_8 - cosl_8*sint_8*cosa_8)

                  vi_v_8(i,j) = Williamson_ubar_8*sinl_8*cosa_8

                   vicll(i,j) = s_8(2,1)*ui_v_8(i,j) + s_8(2,2)*vi_v_8(i,j)

               end do

            end do

         end if

         !######################
         else !F_stag_L
         !######################

         uicll=0.; vicll=0.

         !Compute U vector for YIN
         !------------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jn

               rlat_8 = G_yg_8(j)

               sint_8 = sin(rlat_8)
               cost_8 = cos(rlat_8)

               do i=g_i0,g_inu

                  rlon_8 = xgu_8(i)

                  sinl_8 = sin(rlon_8-lambda_8)
                  cosl_8 = cos(rlon_8-lambda_8)

                  uicll(i,j) = Williamson_ubar_8*(cost_8*sina_8 - cosl_8*sint_8*cosa_8)

               end do

            end do

         !Compute U vector for YAN
         !------------------------
         else

            do j=g_j0,g_jn

               y_a_8 = G_yg_8(j)

               do i=g_i0,g_inu

                  x_a_8 = xgu_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0)

                  sint_8 = sin(rlat_8)
                  cost_8 = cos(rlat_8)

                  sinl_8 = sin(rlon_8-lambda_8)
                  cosl_8 = cos(rlon_8-lambda_8)

                  ui_u_8(i,j) = Williamson_ubar_8*(cost_8*sina_8 - cosl_8*sint_8*cosa_8)

                  vi_u_8(i,j) = Williamson_ubar_8*sinl_8*cosa_8

                   uicll(i,j) = s_8(1,1)*ui_u_8(i,j) + s_8(1,2)*vi_u_8(i,j)

               end do

            end do

         end if

         !Compute V vector for YIN
         !------------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jnv

               rlat_8 = ygv_8(j)

               sint_8 = sin(rlat_8)
               cost_8 = cos(rlat_8)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i)

                  sinl_8 = sin(rlon_8-lambda_8)
                  cosl_8 = cos(rlon_8-lambda_8)

                  vicll(i,j) = Williamson_ubar_8*sinl_8*cosa_8

               end do

            end do

         !Compute V vector for YAN
         !------------------------
         else

            do j=g_j0,g_jnv

               y_a_8 = ygv_8(j)

               do i=g_i0,g_in

                  x_a_8 = G_xg_8(i)-acos(-1.d0)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0)

                  sint_8 = sin(rlat_8)
                  cost_8 = cos(rlat_8)

                  sinl_8 = sin(rlon_8-lambda_8)
                  cosl_8 = cos(rlon_8-lambda_8)

                  ui_v_8(i,j) = Williamson_ubar_8*(cost_8*sina_8 - cosl_8*sint_8*cosa_8)

                  vi_v_8(i,j) = Williamson_ubar_8*sinl_8*cosa_8

                   vicll(i,j) = s_8(2,1)*ui_v_8(i,j) + s_8(2,2)*vi_v_8(i,j)

               end do

            end do

         end if

         !######################
         end if !F_stag_L
         !######################

      !-----------------------------------------------------------------------
      !Lauritzen et al.,2012,GMD,5,887-901 [Deformational Non-divergent winds]
      !-----------------------------------------------------------------------
      else if (Williamson_NAIR==1) then

         if (F_istep==0.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON CASE1 (Winds) : Deformational no-div. winds, Lauritzen et al.,2012'

         Williamson_ubar_8 = rayt_8/(12.0*24.0*3600.0) !Earth's radius/time to complete one revolution

         time_frac_8 = F_istep*cstv_dt_8/(12.0*24.0*3600.0)

         !######################
         if (.not.F_stag_L) then
         !######################

         uicll=0.; vicll=0.

         !Compute U vector for YIN
         !------------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_in

               rlat_8 = G_yg_8(j)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i)-2.*pi_8*time_frac_8

                  uicll(i,j) = DIX_8*Williamson_ubar_8*sin(rlon_8)**2*sin(2*rlat_8)* &
                               cos(pi_8*time_frac_8) +               &
                               2.*pi_8*Williamson_ubar_8*cos(rlat_8)

               end do

            end do

         !Compute U vector for YAN
         !------------------------
         else

            do j=g_j0,g_jn

               y_a_8 = G_yg_8(j)

               do i=g_i0,g_in

                  x_a_8 = G_xg_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0)-2.*pi_8*time_frac_8

                  ui_u_8(i,j) = DIX_8*Williamson_ubar_8*sin(rlon_8)**2*sin(2*rlat_8)* &
                                cos(pi_8*time_frac_8) + &
                                2.*pi_8*Williamson_ubar_8*cos(rlat_8)

                  vi_u_8(i,j) = DIX_8*Williamson_ubar_8*sin(2.*rlon_8)*cos(rlat_8)*cos(pi_8*time_frac_8)

                  uicll(i,j) = s_8(1,1)*ui_u_8(i,j) + s_8(1,2)*vi_u_8(i,j)

               end do

            end do

         end if

         !Compute V vector for YIN
         !------------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jn

               rlat_8 = G_yg_8(j)

               sint_8 = sin(rlat_8)
               cost_8 = cos(rlat_8)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i)-2.*pi_8*time_frac_8

                  vicll(i,j) = DIX_8*Williamson_ubar_8*sin(2.*rlon_8)*cos(rlat_8)*cos(pi_8*time_frac_8)

               end do

            end do

         !Compute V vector for YAN
         !------------------------
         else

            do j=g_j0,g_jn

               y_a_8 = G_yg_8(j)

               do i=g_i0,g_in

                  x_a_8 = G_xg_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 - 2.*pi_8*time_frac_8

                  ui_v_8(i,j) = DIX_8*Williamson_ubar_8*sin(rlon_8)**2*sin(2*rlat_8)* &
                                cos(pi_8*time_frac_8) +               &
                                2.*pi_8*Williamson_ubar_8*cos(rlat_8)

                  vi_v_8(i,j) = DIX_8*Williamson_ubar_8*sin(2.*rlon_8)*cos(rlat_8)*cos(pi_8*time_frac_8)

                  vicll(i,j) = s_8(2,1)*ui_v_8(i,j) + s_8(2,2)*vi_v_8(i,j)

               end do

            end do

         end if

         !######################
         else !F_stag_L
         !######################

         uicll=0.; vicll=0.

         !Compute U vector for YIN
         !------------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jn

               rlat_8 = G_yg_8(j)

               do i=g_i0,g_inu

                  rlon_8 = xgu_8(i)-2.*pi_8*time_frac_8

                  uicll(i,j) = DIX_8*Williamson_ubar_8*sin(rlon_8)**2*sin(2*rlat_8)* &
                               cos(pi_8*time_frac_8) +               &
                               2.*pi_8*Williamson_ubar_8*cos(rlat_8)

               end do

            end do

         !Compute U vector for YAN
         !------------------------
         else

            do j=g_j0,g_jn

               y_a_8 = G_yg_8(j)

               do i=g_i0,g_inu

                  x_a_8 = xgu_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0)-2.*pi_8*time_frac_8

                  ui_u_8(i,j) = DIX_8*Williamson_ubar_8*sin(rlon_8)**2*sin(2*rlat_8)* &
                                cos(pi_8*time_frac_8) + &
                                2.*pi_8*Williamson_ubar_8*cos(rlat_8)

                  vi_u_8(i,j) = DIX_8*Williamson_ubar_8*sin(2.*rlon_8)*cos(rlat_8)*cos(pi_8*time_frac_8)

                  uicll(i,j) = s_8(1,1)*ui_u_8(i,j) + s_8(1,2)*vi_u_8(i,j)

               end do

            end do

         end if

         !Compute V vector for YIN
         !------------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jnv

               rlat_8 = ygv_8(j)

               sint_8 = sin(rlat_8)
               cost_8 = cos(rlat_8)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i)-2.*pi_8*time_frac_8

                  vicll(i,j) = DIX_8*Williamson_ubar_8*sin(2.*rlon_8)*cos(rlat_8)*cos(pi_8*time_frac_8)

               end do

            end do

         !Compute V vector for YAN
         !------------------------
         else

            do j=g_j0,g_jnv

               y_a_8 = ygv_8(j)

               do i=g_i0,g_in

                  x_a_8 = G_xg_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 - 2.*pi_8*time_frac_8

                  ui_v_8(i,j) = DIX_8*Williamson_ubar_8*sin(rlon_8)**2*sin(2*rlat_8)* &
                                cos(pi_8*time_frac_8) +               &
                                2.*pi_8*Williamson_ubar_8*cos(rlat_8)

                  vi_v_8(i,j) = DIX_8*Williamson_ubar_8*sin(2.*rlon_8)*cos(rlat_8)*cos(pi_8*time_frac_8)

                  vicll(i,j) = s_8(2,1)*ui_v_8(i,j) + s_8(2,2)*vi_v_8(i,j)

               end do

            end do

         end if

         !######################
         end if !F_stag_L
         !######################

      !-------------------------------------------------------------------
      !Lauritzen et al.,2012,GMD,5,887-901 [Deformational divergent winds]
      !-------------------------------------------------------------------
      else if (Williamson_NAIR==2) then

         if (F_istep==0.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON CASE1 (Winds) : Deformational div. winds, Lauritzen et al.,2012'

         Williamson_ubar_8 = rayt_8/(12.0*24.0*3600.0) !Earth's radius/time to complete one revolution

         time_frac_8 = F_istep*cstv_dt_8/(12.0*24.0*3600.0)

         !######################
         if (.not.F_stag_L) then
         !######################

         uicll=0.; vicll=0.

         !Compute U vector for YIN
         !------------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jn

               rlat_8 = G_yg_8(j)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i) - 2.*pi_8*time_frac_8

                  uicll(i,j) = -FIVE_8*Williamson_ubar_8*sin(0.5*rlon_8)**2*sin(2*rlat_8)* &
                               cos(rlat_8)**2*cos(pi_8*time_frac_8) + &
                               2.*pi_8*Williamson_ubar_8*cos(rlat_8)

               end do

            end do

         !Compute U vector for YAN
         !------------------------
         else

            do j=g_j0,g_jn

               y_a_8 = G_yg_8(j)

               do i=g_i0,g_in

                  x_a_8 = G_xg_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0) - 2.*pi_8*time_frac_8

                  ui_u_8(i,j) = -FIVE_8*Williamson_ubar_8*sin(0.5*rlon_8)**2*sin(2*rlat_8)* &
                                cos(rlat_8)**2*cos(pi_8*time_frac_8) + &
                                2.*pi_8*Williamson_ubar_8*cos(rlat_8)

                  vi_u_8(i,j) = (FIVE_8/TWO_8)*Williamson_ubar_8*sin(rlon_8)*cos(rlat_8)**3*cos(pi_8*time_frac_8)

                  uicll(i,j) = s_8(1,1)*ui_u_8(i,j) + s_8(1,2)*vi_u_8(i,j)

               end do

            end do

         end if

         !Compute V vector for YIN
         !------------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jn

               rlat_8 = G_yg_8(j)

               sint_8 = sin(rlat_8)
               cost_8 = cos(rlat_8)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i) - 2.*pi_8*time_frac_8

                  vicll(i,j) = (FIVE_8/TWO_8)*Williamson_ubar_8*sin(rlon_8)*cos(rlat_8)**3*cos(pi_8*time_frac_8)

               end do

            end do

         !Compute V vector for YAN
         !------------------------
         else

            do j=g_j0,g_jn

               y_a_8 = G_yg_8(j)

               do i=g_i0,g_in

                  x_a_8 = G_xg_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 - 2.*pi_8*time_frac_8

                  ui_v_8(i,j) = -FIVE_8*Williamson_ubar_8*sin(0.5*rlon_8)**2*sin(2*rlat_8)* &
                                cos(rlat_8)**2*cos(pi_8*time_frac_8) + &
                                2.*pi_8*Williamson_ubar_8*cos(rlat_8)

                  vi_v_8(i,j) = (FIVE_8/TWO_8)*Williamson_ubar_8*sin(rlon_8)*cos(rlat_8)**3*cos(pi_8*time_frac_8)

                  vicll(i,j)= s_8(2,1)*ui_v_8(i,j) + s_8(2,2)*vi_v_8(i,j)

               end do

            end do

         end if

         !######################
         else !F_stag_L
         !######################

         uicll=0.; vicll=0.

         !Compute U vector for YIN
         !------------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jn

               rlat_8 = G_yg_8(j)

               do i=g_i0,g_inu

                  rlon_8 = xgu_8(i) - 2.*pi_8*time_frac_8

                  uicll(i,j) = -FIVE_8*Williamson_ubar_8*sin(0.5*rlon_8)**2*sin(2*rlat_8)* &
                               cos(rlat_8)**2*cos(pi_8*time_frac_8) + &
                               2.*pi_8*Williamson_ubar_8*cos(rlat_8)

               end do

            end do

         !Compute U vector for YAN
         !------------------------
         else

            do j=g_j0,g_jn

               y_a_8 = G_yg_8(j)

               do i=g_i0,g_inu

                  x_a_8 = xgu_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0) - 2.*pi_8*time_frac_8

                  ui_u_8(i,j) = -FIVE_8*Williamson_ubar_8*sin(0.5*rlon_8)**2*sin(2*rlat_8)* &
                                cos(rlat_8)**2*cos(pi_8*time_frac_8) + &
                                2.*pi_8*Williamson_ubar_8*cos(rlat_8)

                  vi_u_8(i,j) = (FIVE_8/TWO_8)*Williamson_ubar_8*sin(rlon_8)*cos(rlat_8)**3*cos(pi_8*time_frac_8)

                  uicll(i,j) = s_8(1,1)*ui_u_8(i,j) + s_8(1,2)*vi_u_8(i,j)

               end do

            end do

         end if

         !Compute V vector for YIN
         !------------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jnv

               rlat_8 = ygv_8(j)

               sint_8 = sin(rlat_8)
               cost_8 = cos(rlat_8)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i) - 2.*pi_8*time_frac_8

                  vicll(i,j) = (FIVE_8/TWO_8)*Williamson_ubar_8*sin(rlon_8)*cos(rlat_8)**3*cos(pi_8*time_frac_8)

               end do

            end do

         !Compute V vector for YAN
         !------------------------
         else

            do j=g_j0,g_jnv

               y_a_8 = ygv_8(j)

               do i=g_i0,g_in

                  x_a_8 = G_xg_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 - 2.*pi_8*time_frac_8

                  ui_v_8(i,j) = -FIVE_8*Williamson_ubar_8*sin(0.5*rlon_8)**2*sin(2*rlat_8)* &
                                cos(rlat_8)**2*cos(pi_8*time_frac_8) + &
                                2.*pi_8*Williamson_ubar_8*cos(rlat_8)

                  vi_v_8(i,j) = (FIVE_8/TWO_8)*Williamson_ubar_8*sin(rlon_8)*cos(rlat_8)**3*cos(pi_8*time_frac_8)

                  vicll(i,j)= s_8(2,1)*ui_v_8(i,j) + s_8(2,2)*vi_v_8(i,j)

               end do

            end do

         end if

         !######################
         end if !F_stag_L
         !######################

      !--------------------------------------------------------------
      !Nair and Machenhauer,2002,MWR,130,649-667 [Deformational flow]
      !--------------------------------------------------------------
      else if (Williamson_NAIR==3) then

         if (F_istep==0.and.Lun_out>0) write(Lun_out,*) 'WILLIAMSON CASE1 (Winds) : Circular vortex, Nair and Machenhauer,2002'

         rlon0_8 = Williamson_rlon0
         rlat0_8 = Williamson_rlat0

         !######################
         if (.not.F_stag_L) then
         !######################

         uicll=0.; vicll=0.

         !Compute U vector for YIN
         !------------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jn

               rlat_8 = G_yg_8(j)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i)

                  rlonr_8 = atan2( cos(rlat_8)*sin(rlon_8-rlon0_8), cos(rlat_8)*sin(rlat0_8)*cos(rlon_8-rlon0_8)-cos(rlat0_8)*sin(rlat_8) )
                  if (rlonr_8 < 0.0d0) rlonr_8 = rlonr_8 + 2.*pi_8

                  rlatr_8 = asin( sin(rlat_8)*sin(rlat0_8) + cos(rlat_8)*cos(rlat0_8)*cos(rlon_8-rlon0_8) )

                  rho_8 = Williamson_rho_i_8*cos(rlatr_8)

                  vt_8  = Williamson_v0_8*(3.*sqrt(3.)/2) * (1.0d0/cosh(rho_8))**2 * tanh(rho_8)

                  wil_omega_8 = 0.0d0
                  if (rho_8 /= 0.0d0) wil_omega_8 = vt_8/(rayt_8*rho_8)

                  uicll(i,j) = rayt_8*wil_omega_8 * (sin(rlat0_8)*cos(rlat_8) - cos(rlat0_8)*cos(rlon_8-rlon0_8)*sin(rlat_8))

               end do

            end do

         !Compute U vector for YAN
         !------------------------
         else

            do j=g_j0,g_jn

               y_a_8 = G_yg_8(j)

               do i=g_i0,g_in

                  x_a_8 = G_xg_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0)

                  rlonr_8 = atan2( cos(rlat_8)*sin(rlon_8-rlon0_8), cos(rlat_8)*sin(rlat0_8)*cos(rlon_8-rlon0_8)-cos(rlat0_8)*sin(rlat_8) )
                  if (rlonr_8 < 0.0d0) rlonr_8 = rlonr_8 + 2.*pi_8

                  rlatr_8 = asin( sin(rlat_8)*sin(rlat0_8) + cos(rlat_8)*cos(rlat0_8)*cos(rlon_8-rlon0_8) )

                  rho_8 = Williamson_rho_i_8*cos(rlatr_8)

                  vt_8  = Williamson_v0_8*(3.*sqrt(3.)/2) * (1.0d0/cosh(rho_8))**2 * tanh(rho_8)

                  wil_omega_8 = 0.0d0
                  if (rho_8 /= 0.0d0) wil_omega_8 = vt_8/(rayt_8*rho_8)

                  ui_u_8(i,j) = rayt_8*wil_omega_8 * (sin(rlat0_8)*cos(rlat_8) - cos(rlat0_8)*cos(rlon_8-rlon0_8)*sin(rlat_8))

                  vi_u_8(i,j) = rayt_8*wil_omega_8 * cos(rlat0_8) * sin(rlon_8-rlon0_8)

                   uicll(i,j) = s_8(1,1)*ui_u_8(i,j) + s_8(1,2)*vi_u_8(i,j)

               end do

            end do

         end if

         !Compute V vector for YIN
         !------------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jn

               rlat_8 = G_yg_8(j)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i)

                  rlonr_8 = atan2( cos(rlat_8)*sin(rlon_8-rlon0_8), cos(rlat_8)*sin(rlat0_8)*cos(rlon_8-rlon0_8)-cos(rlat0_8)*sin(rlat_8) )
                  if (rlonr_8 < 0.0d0) rlonr_8 = rlonr_8 + 2.*pi_8

                  rlatr_8 = asin( sin(rlat_8)*sin(rlat0_8) + cos(rlat_8)*cos(rlat0_8)*cos(rlon_8-rlon0_8) )

                  rho_8 = Williamson_rho_i_8*cos(rlatr_8)

                  vt_8 = Williamson_v0_8*(3.*sqrt(3.)/2) * (1.0d0/cosh(rho_8))**2 * tanh(rho_8)

                  wil_omega_8 = 0.0d0
                  if (rho_8 /= 0.0d0) wil_omega_8 = vt_8/(rayt_8*rho_8)

                  vicll(i,j) = rayt_8*wil_omega_8 * cos(rlat0_8) * sin(rlon_8-rlon0_8)

               end do

            end do

         !Compute V vector for YAN
         !------------------------
         else

            do j=g_j0,g_jn

               y_a_8 = G_yg_8(j)

               do i=g_i0,g_in

                  x_a_8 = G_xg_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0)

                  rlonr_8 = atan2( cos(rlat_8)*sin(rlon_8-rlon0_8), cos(rlat_8)*sin(rlat0_8)*cos(rlon_8-rlon0_8)-cos(rlat0_8)*sin(rlat_8) )
                  if (rlonr_8 < 0.0d0) rlonr_8 = rlonr_8 + 2.*pi_8

                  rlatr_8 = asin( sin(rlat_8)*sin(rlat0_8) + cos(rlat_8)*cos(rlat0_8)*cos(rlon_8-rlon0_8) )

                  rho_8 = Williamson_rho_i_8*cos(rlatr_8)

                  vt_8  = Williamson_v0_8*(3.*sqrt(3.)/2) * (1.0d0/cosh(rho_8))**2 * tanh(rho_8)

                  wil_omega_8 = 0.0d0
                  if (rho_8 /= 0.0d0) wil_omega_8 = vt_8/(rayt_8*rho_8)

                  ui_v_8(i,j) = rayt_8*wil_omega_8 * (sin(rlat0_8)*cos(rlat_8) - cos(rlat0_8)*cos(rlon_8-rlon0_8)*sin(rlat_8))

                  vi_v_8(i,j) = rayt_8*wil_omega_8 * cos(rlat0_8) * sin(rlon_8-rlon0_8)

                   vicll(i,j) = s_8(2,1)*ui_v_8(i,j) + s_8(2,2)*vi_v_8(i,j)

               end do

            end do

         end if

         !######################
         else !F_stag_L
         !######################

         uicll=0.; vicll=0.

         !Compute U vector for YIN
         !------------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jn

               rlat_8 = G_yg_8(j)

               do i=g_i0,g_inu

                  rlon_8 = xgu_8(i)

                  rlonr_8 = atan2( cos(rlat_8)*sin(rlon_8-rlon0_8), cos(rlat_8)*sin(rlat0_8)*cos(rlon_8-rlon0_8)-cos(rlat0_8)*sin(rlat_8) )
                  if (rlonr_8 < 0.0d0) rlonr_8 = rlonr_8 + 2.*pi_8

                  rlatr_8 = asin( sin(rlat_8)*sin(rlat0_8) + cos(rlat_8)*cos(rlat0_8)*cos(rlon_8-rlon0_8) )

                  rho_8 = Williamson_rho_i_8*cos(rlatr_8)

                  vt_8  = Williamson_v0_8*(3.*sqrt(3.)/2) * (1.0d0/cosh(rho_8))**2 * tanh(rho_8)

                  wil_omega_8 = 0.0d0
                  if (rho_8 /= 0.0d0) wil_omega_8 = vt_8/(rayt_8*rho_8)

                  uicll(i,j) = rayt_8*wil_omega_8 * (sin(rlat0_8)*cos(rlat_8) - cos(rlat0_8)*cos(rlon_8-rlon0_8)*sin(rlat_8))

               end do

            end do

         !Compute U vector for YAN
         !------------------------
         else

            do j=g_j0,g_jn

               y_a_8 = G_yg_8(j)

               do i=g_i0,g_inu

                  x_a_8 = xgu_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0)

                  rlonr_8 = atan2( cos(rlat_8)*sin(rlon_8-rlon0_8), cos(rlat_8)*sin(rlat0_8)*cos(rlon_8-rlon0_8)-cos(rlat0_8)*sin(rlat_8) )
                  if (rlonr_8 < 0.0d0) rlonr_8 = rlonr_8 + 2.*pi_8

                  rlatr_8 = asin( sin(rlat_8)*sin(rlat0_8) + cos(rlat_8)*cos(rlat0_8)*cos(rlon_8-rlon0_8) )

                  rho_8 = Williamson_rho_i_8*cos(rlatr_8)

                  vt_8  = Williamson_v0_8*(3.*sqrt(3.)/2) * (1.0d0/cosh(rho_8))**2 * tanh(rho_8)

                  wil_omega_8 = 0.0d0
                  if (rho_8 /= 0.0d0) wil_omega_8 = vt_8/(rayt_8*rho_8)

                  ui_u_8(i,j) = rayt_8*wil_omega_8 * (sin(rlat0_8)*cos(rlat_8) - cos(rlat0_8)*cos(rlon_8-rlon0_8)*sin(rlat_8))

                  vi_u_8(i,j) = rayt_8*wil_omega_8 * cos(rlat0_8) * sin(rlon_8-rlon0_8)

                   uicll(i,j) = s_8(1,1)*ui_u_8(i,j) + s_8(1,2)*vi_u_8(i,j)

               end do

            end do

         end if

         !Compute V vector for YIN
         !------------------------
         if (Ptopo_couleur==0) then

            do j=g_j0,g_jnv

               rlat_8 = ygv_8(j)

               do i=g_i0,g_in

                  rlon_8 = G_xg_8(i)

                  rlonr_8 = atan2( cos(rlat_8)*sin(rlon_8-rlon0_8), cos(rlat_8)*sin(rlat0_8)*cos(rlon_8-rlon0_8)-cos(rlat0_8)*sin(rlat_8) )
                  if (rlonr_8 < 0.0d0) rlonr_8 = rlonr_8 + 2.*pi_8

                  rlatr_8 = asin( sin(rlat_8)*sin(rlat0_8) + cos(rlat_8)*cos(rlat0_8)*cos(rlon_8-rlon0_8) )

                  rho_8 = Williamson_rho_i_8*cos(rlatr_8)

                  vt_8 = Williamson_v0_8*(3.*sqrt(3.)/2) * (1.0d0/cosh(rho_8))**2 * tanh(rho_8)

                  wil_omega_8 = 0.0d0
                  if (rho_8 /= 0.0d0) wil_omega_8 = vt_8/(rayt_8*rho_8)

                  vicll(i,j) = rayt_8*wil_omega_8 * cos(rlat0_8) * sin(rlon_8-rlon0_8)

               end do

            end do

         !Compute V vector for YAN
         !------------------------
         else

            do j=g_j0,g_jnv

               y_a_8 = ygv_8(j)

               do i=g_i0,g_in

                  x_a_8 = G_xg_8(i) - acos(-1.d0)

                  call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

                  rlon_8 = rlon_8 + acos(-1.d0)

                  rlonr_8 = atan2( cos(rlat_8)*sin(rlon_8-rlon0_8), cos(rlat_8)*sin(rlat0_8)*cos(rlon_8-rlon0_8)-cos(rlat0_8)*sin(rlat_8) )
                  if (rlonr_8 < 0.0d0) rlonr_8 = rlonr_8 + 2.*pi_8

                  rlatr_8 = asin( sin(rlat_8)*sin(rlat0_8) + cos(rlat_8)*cos(rlat0_8)*cos(rlon_8-rlon0_8) )

                  rho_8 = Williamson_rho_i_8*cos(rlatr_8)

                  vt_8  = Williamson_v0_8*(3.*sqrt(3.)/2) * (1.0d0/cosh(rho_8))**2 * tanh(rho_8)

                  wil_omega_8 = 0.0d0
                  if (rho_8 /= 0.0d0) wil_omega_8 = vt_8/(rayt_8*rho_8)

                  ui_v_8(i,j) = rayt_8*wil_omega_8 * (sin(rlat0_8)*cos(rlat_8) - cos(rlat0_8)*cos(rlon_8-rlon0_8)*sin(rlat_8))

                  vi_v_8(i,j) = rayt_8*wil_omega_8 * cos(rlat0_8) * sin(rlon_8-rlon0_8)

                   vicll(i,j) = s_8(2,1)*ui_v_8(i,j) + s_8(2,2)*vi_v_8(i,j)

               end do

            end do

         end if

         !######################
         end if !F_stag_L
         !######################

      end if

      zlist = 1

      call glbdist_os (uicll,uloc,&
                       F_minx,F_maxx,F_miny,F_maxy,1,&
                       G_ni+G_halox,G_nj+G_haloy,zlist,1,1.0d0,0.0d0)

      call glbdist_os (vicll,vloc,&
                       F_minx,F_maxx,F_miny,F_maxy,1,&
                       G_ni+G_halox,G_nj+G_haloy,zlist,1,1.0d0,0.0d0)

      !######################
      if (.not.F_stag_L) then
      !######################

      do k=1,F_nk
         F_u(i0:in,j0:jn,k) = uloc(i0:in,j0:jn)
      end do

      do k=1,F_nk
         F_v(i0:in,j0:jn,k) = vloc(i0:in,j0:jn)
      end do

      !######################
      else !F_stag_L
      !######################

      do k=1,F_nk
         F_u(i0:inu,j0:jn,k) = uloc(i0:inu,j0:jn)
      end do

      do k=1,F_nk
         F_v(i0:in,j0:jnv,k) = vloc(i0:in,j0:jnv)
      end do

      !######################
      end if !F_stag_L
      !######################

      if (F_istep==0.and.Lun_out>0) write(Lun_out,*) '--------------------------------------------------------------'

      return
      end
