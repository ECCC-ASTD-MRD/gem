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

!**s/r wil_case9 - To setup Williamson Case 9 == The Matsuno baroclinic wave test case (HEIGHT)

      subroutine wil_case9 (F_gz,F_minx,F_maxx,F_miny,F_maxy,F_nk,F_istep)

      use cstv
      use gem_options
      use glb_ld
      use lun
      use ptopo
      use tdpack, only: grav_8
      use wil_options
      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer, intent(in)  :: F_minx,F_maxx,F_miny,F_maxy,F_nk,F_istep
      real,    intent(out) :: F_gz(F_minx:F_maxx,F_miny:F_maxy,F_nk)

      !object
      !=================================================================================
      !     To setup Williamson Case 9 == The Matsuno baroclinic wave test case (HEIGHT)
      !     ----------------------------------------------------------------------------
      !     Shamir et al.,2019,GMD,12,2181-2193
      !     ----------------------------------------------------------------------------
      !     Notes: This function supports k>=1 and n>=1 inputs only.
      !     Special treatments are required for k=0 and n=-1,0/-.
      !=================================================================================

      integer :: i,j,k

      real :: picll(G_ni,G_nj),gzloc(F_minx:F_maxx,F_miny:F_maxy)

      real(kind=REAL64) :: rlon_8,rlat_8,s_8(2,2),x_a_8,y_a_8,time_8, &
                           wil_eval_fld,wil_eval_omega

      external :: wil_eval_fld,wil_eval_omega
!
!---------------------------------------------------------------------
!
      if (Lun_out>0) write(Lun_out,*) ''
      if (Lun_out>0) write(Lun_out,*) '--------------------------------------------'
      if (Lun_out>0) write(Lun_out,*) 'WILLIAMSON CASE9, Shamir et al. (2019)      '
      if (Lun_out>0) write(Lun_out,*) 'The Matsuno baroclinic wave test case       '
      if (Lun_out>0) write(Lun_out,*) '--------------------------------------------'

      !Evaluate time (sec)
      !-------------------
      time_8 = F_istep*cstv_dt_8

      !Unpack Earth dictionary and Evaluate wave frequency
      !---------------------------------------------------
      Williamson_omega_8 = wil_eval_omega (Williamson_k, Williamson_n)

      !Compute tracer for YIN
      !----------------------
      if (Ptopo_couleur==0) then

         do j=1,G_nj

            rlat_8 = G_yg_8(j)

            do i=1,G_ni

               rlon_8 = G_xg_8(i)

               picll(i,j) = wil_eval_fld (rlat_8, rlon_8, time_8, Williamson_k, Williamson_n, Williamson_amp_8, 'phi')/grav_8 + &
                                          Williamson_mean_depth_8

            end do

         end do

      !Compute tracer for YAN
      !----------------------
      else

         do j=1,G_nj

            y_a_8 = G_yg_8(j)

            do i=1,G_ni

               x_a_8 = G_xg_8(i)-acos(-1.D0)

               call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

               rlon_8 = rlon_8+acos(-1.D0)

               picll(i,j) = wil_eval_fld (rlat_8, rlon_8, time_8, Williamson_k, Williamson_n, Williamson_amp_8, 'phi')/grav_8 + &
                                          Williamson_mean_depth_8

            end do

         end do

      end if

      call glbdist (picll,G_ni,G_nj,gzloc,F_minx,F_maxx,F_miny,F_maxy,1,G_halox,G_haloy)

      do k=1,F_nk
         F_gz(1:l_ni,1:l_nj,k) = gzloc(1:l_ni,1:l_nj)
      end do
!
!---------------------------------------------------------------------
!
      return
      end
