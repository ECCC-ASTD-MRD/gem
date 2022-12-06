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

!**s/r wil_Terminator_0 - Initial conditions for Terminator "Toy" chemistry module in Williamson's case

      subroutine wil_Terminator_0 (F_cl,F_cl2,Mminx,Mmaxx,Mminy,Mmaxy,Nk)

      use geomh
      use Terminator

      use glb_ld
      use lun
      use ptopo
      use gem_options
      implicit none

      integer Mminx,Mmaxx,Mminy,Mmaxy,Nk
      real F_cl(Mminx:Mmaxx,Mminy:Mmaxy,Nk),F_cl2(Mminx:Mmaxx,Mminy:Mmaxy,Nk)

      !object
      !===============================================================================
      !  Initial conditions for Terminator "Toy" chemistry module in Williamson's case
      !===============================================================================


      !-----------------------------------------------------------------------

      integer i,j,k,i0,in,j0,jn

      real(kind=REAL64) x_a_8,y_a_8,s_8(2,2),rlon_8

      real(kind=REAL64)  :: lon,     & ! Longitude (radians)
                  lat,     & ! Latitude (radians)
                  lon_d,   & ! Longitude (degrees)
                  lat_d      ! Latitude (degrees)

      real(kind=REAL64)  :: cl,cl2     ! molar mixing ratio of cl and cl2

      real(kind=REAL64), parameter :: pi = 3.14159265358979d0 ! pi

      real(kind=REAL64), parameter :: radians_to_degrees = 180.0_8/pi

      !-----------------------------------------------------------------------

      if (Lun_out>0) write (Lun_out,1000)

      i0= 1-G_halox ; in= l_ni+G_halox
      j0= 1-G_haloy ; jn= l_nj+G_haloy

      !Initial conditions: CL,CL2
      !--------------------------
      do k = 1,Nk

         do j = j0,jn

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = i0,in

                  lon = geomh_x_8(i)

                  lon_d = lon * radians_to_degrees
                  lat_d = lat * radians_to_degrees

                  call initial_value_Terminator (lat_d,lon_d,cl,cl2)

                  !Chemical tracers
                  !----------------
                  F_cl (i,j,k) = cl
                  F_cl2(i,j,k) = cl2

               end do

            else

               do i = i0,in

                  x_a_8 = geomh_x_8(i) - acos(-1.D0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.D0)

                  lon_d = lon * radians_to_degrees
                  lat_d = lat * radians_to_degrees

                  call initial_value_Terminator (lat_d,lon_d,cl,cl2)

                  !Chemical tracers
                  !----------------
                  F_cl (i,j,k) = cl
                  F_cl2(i,j,k) = cl2

               end do

            end if

         end do

      end do

      !---------------------------------------------------------------

      return

 1000 format( &
      /,'USE INITIAL CONDITIONS FOR TERMINATOR CHEMISTRY IN WILLIAMSON: (S/R WIL_TERMINATOR_0)', &
      /,'=====================================================================================')

      end subroutine wil_Terminator_0
