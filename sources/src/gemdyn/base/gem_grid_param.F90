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

!**s/r gem_grid_extension - returns the grid parameters

      subroutine gem_grid_param (F_bsc_base, F_bsc_ext1, F_extension     ,&
                                 F_maxcfl, F_iref, F_jref, F_lonr, F_latr,&
                                 F_ni, F_nj, F_dx, F_dy, F_x0_8, F_y0_8  ,&
                                 F_xl_8, F_yl_8, F_overlap, F_yinyang_L  ,&
                                 F_uout, F_err)
      implicit none
#include <arch_specific.hf>

      logical, intent(IN)  :: F_yinyang_L ! is this a Yin-Yang grid?
      integer, intent(IN)  :: F_maxcfl    ! Max Supported Courrant number at global boundaries
      integer, intent(OUT) :: F_bsc_base  ! Basic global lateral boundary conditions width
      integer, intent(OUT) :: F_bsc_ext1  ! Added extension for proper de-staggering of u,v at physics interface
      integer, intent(IN)  :: F_iref, F_jref
      integer, intent(OUT) :: F_extension ! Total extension to user specified grid configuration

      integer, intent(IN)   :: F_uout         ! Write error messages to unit F_uout
      integer, intent(OUT)  :: F_err          ! Error code (0=OK ; -1=error)
      integer, intent(INOUT):: F_ni, F_nj     ! user  grid dimensions as input
                                              ! final grid dimensions as output
      real, intent(INOUT) :: F_dx, F_dy       ! grid spacing in x and y directions
      real, intent(IN)    :: F_overlap        ! Overlap extent for Yin-Yang grid
      real, intent(IN   ) :: F_lonr, F_latr   ! coordinates of reference point in rotated space
      real*8, intent(OUT) :: F_x0_8, F_y0_8   ! coordinates of lower left  corner in rotated space
      real*8, intent(OUT) :: F_xl_8, F_yl_8   ! coordinates of upper right corner in rotated space

      integer iref, jref
      real*8 delta_8, lonr, latr
!
!-------------------------------------------------------------------
!
      ! basic global lateral boundary conditions width
      F_bsc_base = 5
      if(F_yinyang_L) F_bsc_base=F_bsc_base+1

      ! added points for proper de-staggering of u,v at physics interface
      F_bsc_ext1 = 2

      ! total extension to user specified grid configuration
      F_extension= F_maxcfl + F_bsc_base + F_bsc_ext1

      if (F_yinyang_L) then

         F_x0_8 =   45.0 - 3.0*F_overlap
         F_xl_8 =  315.0 + 3.0*F_overlap
         F_y0_8 = -45.0  -     F_overlap
         F_yl_8 =  45.0  +     F_overlap

         delta_8 = (F_xl_8-F_x0_8)/(F_ni-1)
         F_dx    = delta_8
         F_x0_8  = F_x0_8 - F_extension*delta_8
         F_xl_8  = F_xl_8 + F_extension*delta_8

         delta_8 = (F_yl_8-F_y0_8)/(F_nj-1)
         F_dy    = delta_8
         F_y0_8  = F_y0_8 - F_extension*delta_8
         F_yl_8  = F_yl_8 + F_extension*delta_8

         F_ni = F_ni + 2  *F_extension
         F_nj = F_nj + 2  *F_extension

      else

         if ((F_iref==-1) .and. (F_jref==-1)) then
            iref = F_ni / 2 + F_extension
            if (mod(F_ni,2)==0) then
               lonr = dble(F_lonr) - dble(F_dx)/2.d0
            else
               iref = iref + 1
               lonr = dble(F_lonr)
            end if
            jref = F_nj / 2 + F_extension
            if (mod(F_nj,2)==0) then
               latr = dble(F_latr) - dble(F_dy)/2.d0
            else
               jref = F_nj / 2 + F_extension + 1
               latr = dble(F_latr)
            end if
         else
            iref = F_iref + F_extension
            jref = F_jref + F_extension
            lonr = dble(F_lonr)
            latr = dble(F_latr)
         end if

         F_ni   = F_ni + 2*F_extension
         F_nj   = F_nj + 2*F_extension
         F_x0_8 = lonr - dble(iref-1) * dble(F_dx)
         F_y0_8 = latr - dble(jref-1) * dble(F_dy)
         F_xl_8 = F_x0_8 + dble(F_ni  -1) * dble(F_dx)
         F_yl_8 = F_y0_8 + dble(F_nj  -1) * dble(F_dy)
         if (F_x0_8 < 0.) F_x0_8=F_x0_8+360.
         if (F_xl_8 < 0.) F_xl_8=F_xl_8+360.

         F_err = 0
         if (F_x0_8 < 0.) then
            if (F_uout > 0) write (F_uout,1201) 'Longitude of WEST',F_x0_8,'< 0.'
            F_err = -1
         end if
         if (F_y0_8 < -90.) then
            if (F_uout > 0) write (F_uout,1201) 'Latitude of SOUTH',F_y0_8,'< -90.'
            F_err = -1
         end if
         if (F_xl_8 > 360.) then
            if (F_uout > 0) write (F_uout,1201) 'Longitude of EAST',F_xl_8,'> 360.'
            F_err = -1
         end if
         if (F_yl_8 > 90.) then
            if (F_uout > 0) write (F_uout,1201) 'Latitude of NORTH',F_yl_8,'> 90.'
            F_err = -1
         end if

      end if

 1201 format(/,' WRONG LAM GRID CONFIGURATION --- ABORT ---'/, &
               x,a,' boundary= ',f10.4,x,a/)
!
!-------------------------------------------------------------------
!
      return
      end subroutine gem_grid_param
