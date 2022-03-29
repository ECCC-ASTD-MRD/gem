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

      use, intrinsic :: iso_fortran_env
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

      integer i,j,k,g_i0,g_in,g_j0,g_jn,g_inu,g_jnv,i0,in,j0,jn,inu,jnv,zlist
      real(kind=REAL64) s_8(2,2),x_a_8,y_a_8,sinl_8,cosl_8,sint_8,cost_8, &
             rlon_8,rlat_8,wil_galewski_wind_8,                     &
             ui_u_8(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy), &
             ui_v_8(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy), &
             vi_u_8(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy), &
             vi_v_8(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy), &
             xgu_8(1-G_halox:G_niu+G_halox),                        &
             ygv_8(1-G_haloy:G_njv+G_haloy)
      external wil_galewski_wind_8
      real    uloc(F_minx:F_maxx,F_miny:F_maxy),                    &
              vloc(F_minx:F_maxx,F_miny:F_maxy),                    &
              uicll(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy), &
              vicll(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy)
      real(kind=REAL64) ONE_8, CLXXX_8
      parameter( ONE_8  = 1.0,CLXXX_8 = 180.0 )

      !---------------------------------------------------------------

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

      uicll = 0.; vicll = 0.

      !######################
      if (.not.F_stag_L) then
      !######################

      !Compute U vector for YIN
      !------------------------
      if (Ptopo_couleur==0) then

         do j=g_j0,g_jn

            rlat_8 = G_yg_8(j)

            sint_8 = sin(rlat_8)
            cost_8 = cos(rlat_8)

            do i=g_i0,g_in

               rlon_8 = G_xg_8(i)

               sinl_8 = sin(rlon_8)
               cosl_8 = cos(rlon_8)

               uicll(i,j) = wil_galewski_wind_8(rlat_8,1)

           end do

         end do

      !Compute U vector for YAN
      !------------------------
      else

         do j=g_j0,g_jn

            y_a_8 = G_yg_8(j)

            do i=g_i0,g_in

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

         do j=g_j0,g_jn

            rlat_8 = G_yg_8(j)

            sint_8 = sin(rlat_8)
            cost_8 = cos(rlat_8)

            do i=g_i0,g_in

               rlon_8 = G_xg_8(i)

               sinl_8 = sin(rlon_8)
               cosl_8 = cos(rlon_8)

               vicll(i,j) = 0.0

           end do

         end do

      !Compute V vector for YAN
      !------------------------
      else

         do j=g_j0,g_jn

            y_a_8 = G_yg_8(j)

            do i=g_i0,g_in

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

         do j=g_j0,g_jn

            rlat_8 = G_yg_8(j)

            sint_8 = sin(rlat_8)
            cost_8 = cos(rlat_8)

            do i=g_i0,g_inu

               rlon_8 = xgu_8(i)

               sinl_8 = sin(rlon_8)
               cosl_8 = cos(rlon_8)

               uicll(i,j) = wil_galewski_wind_8(rlat_8,1)

           end do

         end do

      !Compute U vector for YAN
      !------------------------
      else

         do j=g_j0,g_jn

            y_a_8 = G_yg_8(j)

            do i=g_i0,g_inu

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

         do j=g_j0,g_jnv

            rlat_8 = ygv_8(j)

            sint_8 = sin(rlat_8)
            cost_8 = cos(rlat_8)

            do i=g_i0,g_in

               rlon_8 = G_xg_8(i)

               sinl_8 = sin(rlon_8)
               cosl_8 = cos(rlon_8)

               vicll(i,j) = 0.0

           end do

         end do

      !Compute V vector for YAN
      !------------------------
      else

         do j=g_j0,g_jnv

            y_a_8 = ygv_8(j)

            do i=g_i0,g_in

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

      return
      end
