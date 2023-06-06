!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2019 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------

!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet/Project : GEM-MACH
! Fichier/File   : chm_set_dbg_point.ftn90
! Creation       : Deji Akingunola, Fall 2019
! Description    : Locate and set debug grid/PE point from input geographic
!                  coordinates.
!
!
!==============================================================================
!!if_on
integer function chm_set_dbg_point()
!!if_off
   use chm_nml_mod,   only: dbg_lat, dbg_lon
   use chm_utils_mod, only: dbg_jd, dbg_id, chm_lun_out
   use phygridmap,    only: mapmod2phy, phy_glb_gid, phy_lcl_ni, phydim_ni, &
                            phy_lcl_gid
   use ezgrid_mod,    only: ezgrid_find_ij0
   use rpn_comm_itf_mod

   implicit none
   integer(kind=4)               :: my_pe, numproc
   integer(kind=4)               :: i0, j0, lni, lnj
   integer(kind=4)               :: il, int_x, int_y, ier, ij(2)
   integer(kind=4)               :: lx_lbound, lx_ubound, ly_lbound, ly_ubound
   real(kind=4),    dimension(1) :: dbg_lat2x, dbg_lon2y
   integer(kind=4), dimension(3) :: dbg_loc
   integer(kind=4), dimension(:,:), allocatable :: pi0_j0, rpi0_j0
!
!  Declaration of external functions and subroutines
   external rpn_comm_reduce, rpn_comm_bcast
   integer(kind=4), external     :: gdxyfll

   if (abs(dbg_lat) > 90.0) then
      if (chm_lun_out > 0) &
         write(chm_lun_out, *) "CHM debug coordinate not set"
      chm_set_dbg_point = 1
      return
   end if

   chm_set_dbg_point = - 1

   call rpn_comm_rank(RPN_COMM_GRID, my_pe, ier)
   call rpn_comm_size(RPN_COMM_GRID, numproc, ier)
!
   allocate (pi0_j0(5, numproc), rpi0_j0(5, numproc))
   pi0_j0 = 0
   rpi0_j0 = 0
!
   ier = ezgrid_find_ij0(phy_lcl_gid, phy_glb_gid, i0, j0, lni, lnj)
!
   pi0_j0(1, my_pe+1) = i0
   pi0_j0(2, my_pe+1) = j0
   pi0_j0(3, my_pe+1) = lni
   pi0_j0(4, my_pe+1) = lnj
   pi0_j0(5, my_pe+1) = my_pe
   call rpn_comm_reduce(pi0_j0, rpi0_j0, 5*numproc, "MPI_INTEGER", &
                        "MPI_SUM", 0, RPN_COMM_GRID, ier)

   if (my_pe == 0) then

      if (dbg_lon < 0.0) dbg_lon = 360.0 + dbg_lon

      ier = gdxyfll(phy_glb_gid, dbg_lat2x, dbg_lon2y, [dbg_lat], [dbg_lon], 1)
      ! Round to the nearest i-j point
      int_x = nint(dbg_lat2x(1))
      int_y = nint(dbg_lon2y(1))
!
      do il = 1, numproc
         lx_lbound = rpi0_j0(1, il)                        ! i0
         lx_ubound = rpi0_j0(1, il) + rpi0_j0(3, il) - 1   ! in
         ly_lbound = rpi0_j0(2, il)                        ! j0
         ly_ubound = rpi0_j0(2, il) + rpi0_j0(4, il) - 1   ! jn

         if (int_x >= lx_lbound .and. int_x <= lx_ubound .and. &
             int_y >= ly_lbound .and. int_y <= ly_ubound) then
            dbg_loc(1) = int_x
            dbg_loc(2) = int_y
            dbg_loc(3) = rpi0_j0(5, il) ! pe
!            write(*, *) 'dbg info: ', il, int_x, int_y, rpi0_j0(5, il), &
!                        rpi0_j0(1, il), rpi0_j0(2, il), rpi0_j0(3, il), &
!                        rpi0_j0(4, il), int_x - lx_lbound + 1, &
!                        int_y - ly_lbound + 1
            exit
         end if
      end do

      if (dbg_loc(3) < 0) then
         write(*, *) 'WARNING !!! The requested GEM-MACH debug locations ', &
                     dbg_lon, ' and ', dbg_lat, ' are outside the model domain '
      else
         write(*, *) 'GEM-MACH debug co-ordinate (dbg_id, dbg_jd) at (',   &
                      int_x - lx_lbound + 1, ',',  int_y - ly_lbound + 1, &
                      ') on PE=', dbg_loc(3)
      end if
   end if
!
!***  Assign the debug coordinates id, and jd in the appropriate PE
   call RPN_COMM_bcast(dbg_loc, 3, "MPI_INTEGER", 0, RPN_COMM_GRID, ier)
   if (my_pe == dbg_loc(3)) then
      int_x = dbg_loc(1) - i0 + 1
      int_y = dbg_loc(2) - j0 + 1
      ij = mapmod2phy(int_x, int_y, phy_lcl_ni, phydim_ni)
      dbg_id = ij(1)
      dbg_jd = ij(2)
      write(*, *) 'debug location on this PE ', dbg_id, dbg_jd, int_x, int_y
      write(*, *) 'debug location on entire model grid ', dbg_loc(1), dbg_loc(2)
   else
      dbg_id = -1
      dbg_jd = -1
   end if

   chm_set_dbg_point = 1
   return
end function chm_set_dbg_point
