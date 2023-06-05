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

!**s/r canonical_coriolis - compute coriolis FACTOR for AUTOBAROTROPE

      subroutine canonical_coriolis ( F_u_8, F_v_8, F_x_8, F_y_8,  &
                                      F_xu_8, F_yv_8, F_rot_8, Minx,Maxx,Miny,Maxy)
      use wil_options
      use gem_options
      use HORgrid_options
      use tdpack, only : omega_8

      use glb_ld
      use lun
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer Minx,Maxx,Miny,Maxy
      real(kind=REAL64) F_u_8(Minx:Maxx,Miny:Maxy), F_v_8(Minx:Maxx,Miny:Maxy),  &
             F_x_8 (Minx:Maxx), F_y_8 (Miny:Maxy), &
             F_xu_8(Minx:Maxx), F_yv_8(Miny:Maxy), &
             F_rot_8(3,3)
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_u_8        O    - coriolis FACTOR on U grid
! F_v_8        O    - coriolis FACTOR on V grid
! F_x_8        I    - longitudes in radians PHI grid
! F_y_8        I    - latitudes in radians PHI grid
! F_xu_8       I    - longitudes in radians U grid
! F_yv_8       I    - latitudes in radians V grid
! F_rot_8      I    - rotation matrix of the grid


      real(kind=REAL64) ZERO, ONE, TWO
      parameter( ZERO = 0.0 )
      parameter( ONE  = 1.0 )
      parameter( TWO  = 2.0 )

      integer i, j
      real(kind=REAL64)  c0, sa, ca, s0, ang, coef_8, s(2,2), Y_lat, Y_lon, y_a, x_a
!
!-------------------------------------------------------------------
!
      if ( Williamson_case==2 ) goto 100

!     ---------------------------------------
!     Set rotation parameters (standard mode)
!     ---------------------------------------
      s0 = F_rot_8(3,3)

      if ( abs( (abs(s0)-ONE) ) > 1.0e-10 ) then
         ang = atan2( F_rot_8(2,3), F_rot_8(1,3) )
      else
         s0 = sign( ONE, s0 )
         ang = ZERO
      endif

      c0 = sqrt( max( ZERO, ONE - s0 ** 2 ) )

!	processing coriolis FACTOR on V grid
!       ____________________________________

      do j=1-G_haloy,l_nj+G_haloy
         sa = ( TWO * omega_8 ) * s0 * sin(F_yv_8(j))
         ca = ( TWO * omega_8 ) * c0 * cos(F_yv_8(j))
         do i=1-G_halox,l_ni+G_halox
            F_v_8(i,j) = ca * cos(F_x_8(i)-ang) + sa
         enddo
      enddo

!	processing coriolis FACTOR on U grid
!       ____________________________________

      do j=1-G_haloy,l_nj+G_haloy
         sa = ( TWO * omega_8 ) * s0 * sin(F_y_8(j))
         ca = ( TWO * omega_8 ) * c0 * cos(F_y_8(j))
         do i=1-G_halox,l_ni+G_halox
            F_u_8(i,j) = ca * cos(F_xu_8(i) - ang) + sa
         enddo
      enddo

      return

  100 continue

!     ---------------------------------------------------------------------
!     Use Coriolis based on alpha when AUTOBAROTROPE and Williamson_case==2
!     ---------------------------------------------------------------------

      coef_8 = TWO * omega_8

      if (trim(Grd_yinyang_S) == 'YAN') then

!       processing coriolis FACTOR on V grid
!       ____________________________________

        do j=1-G_haloy,l_nj+G_haloy

           y_a = F_yv_8(j)

           do i=1-G_halox,l_ni+G_halox

              x_a = F_x_8(i)-acos(-1.D0)

              call smat (S,Y_lon,Y_lat,x_a,y_a)

              Y_lon = Y_lon+acos(-1.D0)

              F_v_8(i,j) = coef_8 * (-cos(Y_lon)*cos(Y_lat)*sin(Williamson_alpha)+ &
                                                 sin(Y_lat)*cos(Williamson_alpha))

           enddo

        enddo

!       processing coriolis FACTOR on U grid
!       ____________________________________

        do j=1-G_haloy,l_nj+G_haloy

           y_a = F_y_8(j)

           do i=1-G_halox,l_ni+G_halox

              x_a = F_xu_8(i)-acos(-1.D0)

              call smat (S,Y_lon,Y_lat,x_a,y_a)

              Y_lon = Y_lon+acos(-1.D0)

              F_u_8(i,j) = coef_8 * (-cos(Y_lon)*cos(Y_lat)*sin(Williamson_alpha)+ &
                                                 sin(Y_lat)*cos(Williamson_alpha))

           enddo

        enddo

      else

        if (Lun_out>0) then
           write(Lun_out,*) ''
           write(Lun_out,*) 'Coriolis evaluation using alpha when Williamson_case=2'
           write(Lun_out,*) ''
        endif

!       processing coriolis FACTOR on V grid
!       ____________________________________

        do j=1-G_haloy,l_nj+G_haloy
        do i=1-G_halox,l_ni+G_halox
           F_v_8(i,j) = coef_8 * (-cos( F_x_8(i))* cos(F_yv_8(j))*sin(Williamson_alpha)+ &
                                                   sin(F_yv_8(j))*cos(Williamson_alpha))
        enddo
        enddo

!       processing coriolis FACTOR on U grid
!       ____________________________________

        do j=1-G_haloy,l_nj+G_haloy
        do i=1-G_halox,l_ni+G_halox
           F_u_8(i,j) = coef_8 * (-cos( F_xu_8(i))* cos(F_y_8(j))*sin(Williamson_alpha)+ &
                                                    sin(F_y_8(j))*cos(Williamson_alpha))
        enddo
        enddo

      endif
!
!-------------------------------------------------------------------
!
      return
      end
