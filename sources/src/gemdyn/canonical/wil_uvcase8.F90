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

!**s/r wil_uvcase8 - To setup Williamson Case 8 == Galewsky's Case: Barotropic wave (WINDS)

      subroutine wil_uvcase8 (F_u,F_v,F_minx,F_maxx,F_miny,F_maxy,F_nk,F_stag_L)

      use gem_options
      use glb_ld
      use ptopo

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer F_minx,F_maxx,F_miny,F_maxy,F_nk
      real    F_u(F_minx:F_maxx,F_miny:F_maxy,F_nk),F_v(F_minx:F_maxx,F_miny:F_maxy,F_nk)

      logical F_stag_L ! Staggered uv if .T. / Scalar uv if .F.

      !authors
      !     Abdessamad Qaddouri & Vivian Lee
      !
      !revision
      ! v5_00 - Tanguay M. - Clean Up
      !
      !object
      !===========================================================================
      !     To setup Williamson Case 8 == Galewsky's Case: Barotropic wave (WINDS)
      !     Galewsky et al.,2004,Tellus,56A,429-440
      !===========================================================================

      !---------------------------------------------------------------

      integer i,j,k
      real*8 s_8(2,2),x_a_8,y_a_8,sinl_8,cosl_8,sint_8,cost_8, &
             rlon_8,rlat_8,wil_galewski_wind_8,                &
             ui_u_8(G_ni,G_nj),ui_v_8(G_ni,G_nj),              &
             vi_u_8(G_ni,G_nj),vi_v_8(G_ni,G_nj),              &
             xgu_8(G_niu),ygv_8(G_njv)
      external wil_galewski_wind_8
      real    uloc(F_minx:F_maxx,F_miny:F_maxy), &
              vloc(F_minx:F_maxx,F_miny:F_maxy), &
              uicll(G_ni,G_nj),vicll(G_ni,G_nj)
      real*8 ONE_8, CLXXX_8
      parameter( ONE_8  = 1.0,CLXXX_8 = 180.0 )

      !---------------------------------------------------------------

      !U grid
      !------
      do i=1,G_niu
         xgu_8(i) = (G_xg_8(i+1)+G_xg_8(i))*.5
      end do

      !V grid
      !------
      do j=1,G_njv
         ygv_8(j) = (G_yg_8(j+1)+G_yg_8(j))*.5
      end do

      uicll = 0.; vicll = 0.

      !######################
      if (.not.F_stag_L) then
      !######################

      !Compute U vector for YIN
      !------------------------
      if (Ptopo_couleur==0) then

         do j=1,G_nj

            rlat_8 = G_yg_8(j)

            sint_8 = sin(rlat_8)
            cost_8 = cos(rlat_8)

            do i=1,G_ni

               rlon_8 = G_xg_8(i)

               sinl_8 = sin(rlon_8)
               cosl_8 = cos(rlon_8)

               uicll(i,j) = wil_galewski_wind_8(rlat_8,1)

           end do

         end do

      !Compute U vector for YAN
      !------------------------
      else

         do j=1,G_nj

            y_a_8 = G_yg_8(j)

            do i=1,G_ni

               x_a_8 = G_xg_8(i) - acos(-1.D0)

               call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

               rlon_8 = rlon_8 + acos(-1.D0)

               sint_8 = sin(rlat_8)
               cost_8 = cos(rlat_8)

               sinl_8 = sin(rlon_8)
               cosl_8 = cos(rlon_8)

               ui_u_8(i,j) = wil_galewski_wind_8(rlat_8,1)
               vi_u_8(i,j) = 0.0

                uicll(i,j) = s_8(1,1)*ui_u_8(i,j) + s_8(1,2)*vi_u_8(i,j)

           end do

         end do

      end if

      !Compute V vector for YIN
      !------------------------
      if (Ptopo_couleur==0) then

         do j=1,G_nj

            rlat_8 = G_yg_8(j)

            sint_8 = sin(rlat_8)
            cost_8 = cos(rlat_8)

            do i=1,G_ni

               rlon_8 = G_xg_8(i)

               sinl_8 = sin(rlon_8)
               cosl_8 = cos(rlon_8)

               vicll(i,j) = 0.0

           end do

         end do

      !Compute V vector for YAN
      !------------------------
      else

         do j=1,G_nj

            y_a_8 = G_yg_8(j)

            do i=1,G_ni

               x_a_8 = G_xg_8(i) - acos(-1.D0)

               call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

               rlon_8 = rlon_8 + acos(-1.D0)

               sint_8 = sin(rlat_8)
               cost_8 = cos(rlat_8)

               sinl_8 = sin(rlon_8)
               cosl_8 = cos(rlon_8)

               ui_v_8(i,j) = wil_galewski_wind_8(rlat_8,1)
               vi_v_8(i,j) = 0.0

                vicll(i,j) = s_8(2,1)*ui_v_8(i,j) + s_8(2,2)*vi_v_8(i,j)

           end do

         end do

      end if

      !######################
      else !F_stag_L
      !######################

      !Compute U vector for YIN
      !------------------------
      if (Ptopo_couleur==0) then

         do j=1,G_nj

            rlat_8 = G_yg_8(j)

            sint_8 = sin(rlat_8)
            cost_8 = cos(rlat_8)

            do i=1,G_niu

               rlon_8 = xgu_8(i)

               sinl_8 = sin(rlon_8)
               cosl_8 = cos(rlon_8)

               uicll(i,j) = wil_galewski_wind_8(rlat_8,1)

           end do

         end do

      !Compute U vector for YAN
      !------------------------
      else

         do j=1,G_nj

            y_a_8 = G_yg_8(j)

            do i=1,G_niu

               x_a_8 = xgu_8(i) - acos(-1.D0)

               call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

               rlon_8 = rlon_8 + acos(-1.D0)

               sint_8 = sin(rlat_8)
               cost_8 = cos(rlat_8)

               sinl_8 = sin(rlon_8)
               cosl_8 = cos(rlon_8)

               ui_u_8(i,j) = wil_galewski_wind_8(rlat_8,1)
               vi_u_8(i,j) = 0.0

                uicll(i,j) = s_8(1,1)*ui_u_8(i,j) + s_8(1,2)*vi_u_8(i,j)

           end do

         end do

      end if

      !Compute V vector for YIN
      !------------------------
      if (Ptopo_couleur==0) then

         do j=1,G_njv

            rlat_8 = ygv_8(j)

            sint_8 = sin(rlat_8)
            cost_8 = cos(rlat_8)

            do i=1,G_ni

               rlon_8 = G_xg_8(i)

               sinl_8 = sin(rlon_8)
               cosl_8 = cos(rlon_8)

               vicll(i,j) = 0.0

           end do

         end do

      !Compute V vector for YAN
      !------------------------
      else

         do j=1,G_njv

            y_a_8 = ygv_8(j)

            do i=1,G_ni

               x_a_8 = G_xg_8(i) - acos(-1.D0)

               call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

               rlon_8 = rlon_8 + acos(-1.D0)

               sint_8 = sin(rlat_8)
               cost_8 = cos(rlat_8)

               sinl_8 = sin(rlon_8)
               cosl_8 = cos(rlon_8)

               ui_v_8(i,j) = wil_galewski_wind_8(rlat_8,1)
               vi_v_8(i,j) = 0.0

                vicll(i,j) = s_8(2,1)*ui_v_8(i,j) + s_8(2,2)*vi_v_8(i,j)

           end do

         end do

      end if

      !######################
      end if !F_stag_L
      !######################

      call glbdist (uicll,G_ni,G_nj,uloc,F_minx,F_maxx,F_miny,F_maxy,1,G_halox,G_haloy)

      call glbdist (vicll,G_ni,G_nj,vloc,F_minx,F_maxx,F_miny,F_maxy,1,G_halox,G_haloy)

      !######################
      if (.not.F_stag_L) then
      !######################

      do k=1,F_nk
         F_u(1:l_ni,1:l_nj,k) = uloc(1:l_ni,1:l_nj)
      end do

      do k=1,F_nk
         F_v(1:l_ni,1:l_nj,k) = vloc(1:l_ni,1:l_nj)
      end do

      !######################
      else !F_stag_L
      !######################

      do k=1,F_nk
         F_u(1:l_niu,1:l_nj,k) = uloc(1:l_niu,1:l_nj)
      end do

      do k=1,F_nk
         F_v(1:l_ni,1:l_njv,k) = vloc(1:l_ni,1:l_njv)
      end do

      !######################
      end if !F_stag_L
      !######################

      return
      end
