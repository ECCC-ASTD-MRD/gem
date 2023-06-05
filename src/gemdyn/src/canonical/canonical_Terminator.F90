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

!**s/r canonical_terminator_0

      subroutine canonical_terminator_0 (F_cnt)

      use canonical

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      integer, intent(out) :: F_cnt

      logical, external :: wil_check_Terminator,dcmip_check_Terminator
      logical :: Terminator_L
!
!-------------------------------------------------------------------
!
      Terminator_L = wil_check_Terminator() .or. dcmip_check_Terminator()

      if ( Terminator_L ) then

         F_cnt = 0
         cly = 0.

      end if
!
!-------------------------------------------------------------------
!
      return
      end subroutine canonical_terminator_0

!**s/r canonical_terminator_1

      subroutine canonical_terminator_1 (F_src,F_name_S,F_cnt,minx,maxx,miny,maxy,F_ni,F_nj,F_nk)

      use canonical

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      character(len=*), intent(in) :: F_name_S
      integer, intent(in) :: minx,maxx,miny,maxy,F_ni,F_nj,F_nk
      integer, intent(inout) :: F_cnt
      real, intent(in) :: F_src(minx:maxx,miny:maxy,F_nk)

      logical, external :: wil_check_Terminator,dcmip_check_Terminator
      logical :: Terminator_L
!
!-------------------------------------------------------------------
!
      Terminator_L = wil_check_Terminator() .or. dcmip_check_Terminator()

      if ( Terminator_L ) then

         if (F_name_S(1:6)=='TR/CL:' ) cly(1:F_ni,1:F_nj,1:F_nk) = cly(1:F_ni,1:F_nj,1:F_nk) +     F_src(1:F_ni,1:F_nj,1:F_nk)
         if (F_name_S(1:6)=='TR/CL:' ) F_cnt = F_cnt + 1
         if (F_name_S(1:7)=='TR/CL2:') cly(1:F_ni,1:F_nj,1:F_nk) = cly(1:F_ni,1:F_nj,1:F_nk) + 2.0*F_src(1:F_ni,1:F_nj,1:F_nk)
         if (F_name_S(1:7)=='TR/CL2:') F_cnt = F_cnt + 1

      end if
!
!-------------------------------------------------------------------
!
      return
      end subroutine canonical_terminator_1

!**s/r canonical_terminator_2

      subroutine canonical_terminator_2 (F_airmass,F_tracer_8,F_cnt,minx,maxx,miny,maxy,F_nk,F_k0, &
                                         F_unout,F_couleur,F_type_S,F_time_S,F_comment_S)

      use adz_options
      use canonical
      use glb_ld

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      character(len=*), intent(in) :: F_type_S,F_time_S,F_comment_S
      integer, intent(in) :: minx,maxx,miny,maxy,F_nk,F_k0,F_unout,F_couleur
      integer, intent(inout) :: F_cnt
      real(kind=REAL64), intent(out) :: F_tracer_8
      real, intent(in) :: F_airmass(minx:maxx,miny:maxy,F_nk)

      logical,external :: wil_check_Terminator,dcmip_check_Terminator
      logical :: Terminator_L
      integer :: i0,in,j0,jn
!
!-------------------------------------------------------------------
!
      Terminator_L = wil_check_Terminator() .or. dcmip_check_Terminator()

      i0 = 1+pil_w ; j0 = 1+pil_s ; in = l_ni-pil_e ; jn = l_nj-pil_n

      if ( Terminator_L .and. F_cnt==2 ) then

         !Print the mass of CLY scaled by area (CORE)
         !-------------------------------------------
         call mass_tr (F_tracer_8,cly,F_airmass,minx,maxx,miny,maxy,F_nk,i0,in,j0,jn,F_k0)

         if (F_unout>0.and.F_couleur==0) then
            write(F_unout,1002) 'TRACERS: ', F_type_S, F_time_S,'  C= ',F_tracer_8/Adz_gc_area_8,"CLY ",F_comment_S
         end if

         F_cnt = 0

      end if

1002  format(1X,A9,A21,1X,A7,A5,E19.12,1X,A4,1X,A16)
!
!-------------------------------------------------------------------
!
      return
      end subroutine canonical_terminator_2

!**s/r canonical_Terminator - Terminator "Toy" chemistry module for 2D/3D model

      subroutine canonical_Terminator ()

      use canonical
      use cstv
      use dcmip_options
      use geomh
      use glb_ld
      use lun
      use mem_tracers
      use ptopo
      use Terminator
      use wil_options

      use, intrinsic :: iso_fortran_env
      implicit none

      !object
      !===================================================
      !  Terminator "Toy" chemistry module for 2D/3D model
      !  Lauritzen et al.,2015,GMD,8,1299-1313
      !===================================================

      integer :: i,j,k,istat,lower_value

      real(kind=REAL64) :: x_a_8,y_a_8,s_8(2,2),rlon_8

      real(kind=REAL64) :: lon,             &  ! Longitude (radians)
                           lat,             &  ! Latitude (radians)
                           lon_d,           &  ! Longitude (degrees)
                           lat_d,           &  ! Latitude (degrees)
                           cl_f_8, cl2_f_8, &  ! time rate of change of cl and cl2
                           cl_p_8, cl2_p_8     ! molar mixing ratio of cl and cl2 (TIME P)

      real(kind=REAL64), parameter :: pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164_8

      real(kind=REAL64), parameter :: radians_to_degrees = 180.0_8/pi

      real, pointer, dimension (:,:,:) :: cl_p,cl2_p

      real :: lower_value_4, check_4
!
!-------------------------------------------------------------------
!
      !Set admissible lowest value
      !---------------------------
      if (Canonical_williamson_L) lower_value = Williamson_lower_value
      if (Canonical_dcmip_L     ) lower_value = Dcmip_lower_value

      if (lower_value == 0) lower_value_4 = - huge(check_4)
      if (lower_value == 1) lower_value_4 = 0.
      if (lower_value == 2) lower_value_4 = 1.0e-15

      if (Lun_out>0) write (Lun_out,1000) lower_value_4

      istat = tr_get('CL:P',cl_p)
      istat = tr_get('CL2:P',cl2_p)

      do k = 1,l_nk

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  lon_d = lon * radians_to_degrees
                  lat_d = lat * radians_to_degrees

                  cl_p (i,j,k) = max(cl_p (i,j,k),lower_value_4)
                  cl2_p(i,j,k) = max(cl2_p(i,j,k),lower_value_4)

                  cl_p_8  = cl_p (i,j,k)
                  cl2_p_8 = cl2_p(i,j,k)

                  call tendency_Terminator( lat_d, lon_d, cl_p_8, cl2_p_8, Cstv_dt_8, cl_f_8, cl2_f_8 )

                  cl_p (i,j,k) = cl_p_8  + Cstv_dt_8 * cl_f_8
                  cl2_p(i,j,k) = cl2_p_8 + Cstv_dt_8 * cl2_f_8

                  cly  (i,j,k) = cl_p(i,j,k) + 2.0d0 * cl2_p(i,j,k)

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.D0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.D0)

                  lon_d = lon * radians_to_degrees
                  lat_d = lat * radians_to_degrees

                  cl_p (i,j,k) = max(cl_p (i,j,k),lower_value_4)
                  cl2_p(i,j,k) = max(cl2_p(i,j,k),lower_value_4)

                  cl_p_8  = cl_p (i,j,k)
                  cl2_p_8 = cl2_p(i,j,k)

                  call tendency_Terminator( lat_d, lon_d, cl_p_8, cl2_p_8, Cstv_dt_8, cl_f_8, cl2_f_8 )

                  cl_p (i,j,k) = cl_p_8  + Cstv_dt_8 * cl_f_8
                  cl2_p(i,j,k) = cl2_p_8 + Cstv_dt_8 * cl2_f_8

                  cly  (i,j,k) = cl_p(i,j,k) + 2.0d0 * cl2_p(i,j,k)

               end do

            end if

         end do

      end do
!
!-------------------------------------------------------------------
!
      return

 1000 format( &
      /,'USE TOY CHEMISTRY : (S/R CANONICAL_TERMINATOR) with LOWEST_VALUE = ',e9.2,/, &
        '============================================================================')

      end subroutine canonical_Terminator
