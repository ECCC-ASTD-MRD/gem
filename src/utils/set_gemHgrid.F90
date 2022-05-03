!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

subroutine set_gemHgrid4(F_xgi_8, F_ygi_8, F_Grd_ni, F_Grd_nj      , &
     F_Grd_dx, F_Grd_dy, F_Grd_x0_8, F_Grd_xl_8, &
     F_Grd_y0_8, F_Grd_yl_8, F_Grd_yinyang_L)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
   include "rmnlib_basics.inc"
   !@object Compute model horizontal grid
   !@arguments
   logical F_Grd_yinyang_L
   integer F_Grd_ni, F_Grd_nj
   real   F_Grd_dx, F_Grd_dy
   real(REAL64) :: F_xgi_8(F_Grd_ni), F_ygi_8(F_Grd_nj), &
        F_Grd_x0_8, F_Grd_xl_8, F_Grd_y0_8, F_Grd_yl_8

   !     ---------------------------------------------------------------

   call UOG_parpos ( F_xgi_8, F_ygi_8,                               &
        F_Grd_x0_8, F_Grd_y0_8, F_Grd_xl_8, F_Grd_yl_8, &
        F_Grd_ni, F_Grd_nj )

   if (F_Grd_yinyang_L) then

      F_Grd_dx   = abs(F_xgi_8(2)-F_xgi_8(1))
      F_Grd_dy   = abs(F_ygi_8(2)-F_ygi_8(1))

   endif

   !     ---------------------------------------------------------------

   return
end subroutine set_gemHgrid4


subroutine set_gemHgrid3(F_xgi_8, F_ygi_8, F_Grd_ni, F_Grd_nj, F_Grd_dx, F_Grd_dy, &
     F_Grd_x0_8, F_Grd_xl_8, F_Grd_left,                       &
     F_Grd_y0_8, F_Grd_yl_8, F_Grd_belo,                       &
     F_Grd_nila, F_Grd_njla, F_Grd_dxmax, F_Grd_dymax,         &
     F_Grd_yinyang_L, F_Grd_gauss_L, F_LAM, F_Grd_uniform_L,   &
     F_errx, F_erry, F_print_L)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
   include "rmnlib_basics.inc"
   !@arguments
   logical F_Grd_yinyang_L,F_Grd_gauss_L, F_LAM, F_Grd_uniform_L, F_print_L
   integer F_Grd_ni, F_Grd_nj, F_Grd_left, F_Grd_belo, F_Grd_nila,  &
        F_Grd_njla, F_errx, F_erry
   real   F_Grd_dx, F_Grd_dy, F_Grd_dxmax, F_Grd_dymax
   real(REAL64) :: F_xgi_8(F_Grd_ni), F_ygi_8(F_Grd_nj), &
        F_Grd_x0_8, F_Grd_xl_8, F_Grd_y0_8, F_Grd_yl_8
   !@author M. Desgagne -  Fall 2011
   !@revision
   ! v4_40 - Desgagne     - initial version
   ! v4_50 - Desgagne     - call to UOG_parpos for all grids but GV and Gaussian

   integer, external :: stretch_axis2

   integer nimax,njmax,i,j,ni,nila
   real    r1,s1,x0,y0,xl,yl
   real(REAL64) :: halfres, y0_8, yl_8

   !     ---------------------------------------------------------------

   F_Grd_uniform_L = .true.
   F_errx=0 ; F_erry=0

   IF_YY: if (F_Grd_yinyang_L) then

      call UOG_parpos ( F_xgi_8, F_ygi_8,                               &
           F_Grd_x0_8, F_Grd_y0_8, F_Grd_xl_8, F_Grd_yl_8, &
           F_Grd_ni, F_Grd_nj )
      F_Grd_nila = F_Grd_ni
      F_Grd_njla = F_Grd_nj
      F_Grd_left = 0
      F_Grd_belo = 0
      F_Grd_dx   = abs(F_xgi_8(2)-F_xgi_8(1))
      F_Grd_dy   = abs(F_ygi_8(2)-F_ygi_8(1))

   else if (F_LAM) then

      call UOG_parpos ( F_xgi_8, F_ygi_8,                               &
           F_Grd_x0_8, F_Grd_y0_8, F_Grd_xl_8, F_Grd_yl_8, &
           F_Grd_ni, F_Grd_nj )

   else

      IF_GU: if ((F_Grd_ni.eq.F_Grd_nila).and.(.not. F_Grd_gauss_L)) then

         halfres = (F_Grd_yl_8-F_Grd_y0_8)/dble((F_Grd_nj))/2.d0
         y0_8    = F_Grd_y0_8 + halfres
         yl_8    = F_Grd_yl_8 - halfres

         call UOG_parpos ( F_xgi_8, F_ygi_8,                   &
              F_Grd_x0_8, y0_8, F_Grd_xl_8, yl_8, &
              F_Grd_ni+1, F_Grd_nj )

      else ! Gaussian .or. Variable resolution grid

         ni   = F_Grd_ni+1
         nila = F_Grd_nila
         if (F_Grd_ni.eq.F_Grd_nila)  nila = F_Grd_nila+1
         x0 = F_Grd_x0_8
         y0 = F_Grd_y0_8
         xl = F_Grd_xl_8
         yl = F_Grd_yl_8
         F_errx= stretch_axis2( F_xgi_8, F_Grd_dx, x0, xl, F_Grd_left     , &
              F_Grd_ni+1,  nila   , r1, .false., .false., &
              F_Grd_dxmax, nimax, F_Grd_gauss_L )
         F_erry= stretch_axis2( F_ygi_8, F_Grd_dy, y0, yl, F_Grd_belo     , &
              F_Grd_nj, F_Grd_njla, s1, .true. , .false., &
              F_Grd_dymax, njmax, F_Grd_gauss_L )
         F_Grd_uniform_L = .false.

         IF_PRINT: if (F_print_L) then
            x0 = F_xgi_8(1)
            y0 = F_ygi_8(1)
            xl = F_xgi_8(ni  ) !# Warning: out of bound (ni=F_Grd_ni+1)
            yl = F_ygi_8(F_Grd_nj)

            write(6,1020) ni,x0,xl, F_Grd_nj,y0,yl
            write(6,1025) nila,F_Grd_dx,1+F_Grd_left,1+F_Grd_left+nila-1, &
                 F_Grd_njla,F_Grd_dy,1+F_Grd_belo,1+F_Grd_belo+F_Grd_njla-1

            i = ni-nila-F_Grd_left
            j = F_Grd_nj-F_Grd_njla-F_Grd_belo
            write(6,1030) F_Grd_left,i,r1,F_xgi_8(2)-F_xgi_8(1), &
                 F_Grd_belo,j,s1,F_ygi_8(2)-F_ygi_8(1)

            if ( nimax .gt. 0 ) write(6,1035) F_Grd_dxmax, nimax, 'X','X'
            if ( njmax .gt. 0 ) write(6,1035) F_Grd_dymax, njmax, 'Y','Y'
            write(6,1031)
         endif IF_PRINT

      endif IF_GU

   endif IF_YY


1020 FORMAT (/1X,'FINAL HORIZONTAL GRID CONFIGURATION:' &
        /1X,' NI=',I4,' FROM x0=',F9.3,' TO xl=',F9.3,' DEGREES' &
        /1X,' NJ=',I4,' FROM y0=',F9.3,' TO yl=',F9.3,' DEGREES' )
1025 FORMAT(/2X,'THE CONSTANT RESOLUTION AREA HAS:', &
        /1X,' NILA=',I4,' OF GRID-LENGTH=',F9.4,' DEGREES', &
        1x,'(',i4,',',i4,' )', &
        /1X,' NJLA=',I4,' OF GRID-LENGTH=',F9.4,' DEGREES', &
        1x,'(',i4,',',i4,' )')
1030 FORMAT(/2X,'THE VARIABLE RESOLUTION AREA HAS:' &
        /1X,i3,' POINTS TO THE WEST  AND ',i3,' POINTS TO THE EAST' &
        /4x,'WITH STRETCHING FACTOR=',F8.4, &
        ' AND MINIMUM RESOLUTION=',F8.4, &
        /1X,i3,' POINTS ON THE SOUTH AND ',i3,' POINTS ON THE NORTH' &
        /4x,'WITH STRETCHING FACTOR=',F8.4, &
        ' AND MINIMUM RESOLUTION=',F8.4)
1031 FORMAT(1x,66('='))
1035 FORMAT(2x,'RESOLUTION IS LIMITED TO ',F9.4,1x, &
        'DEGREES OVER LAST',I4,' DELTA-',a1,' AT ', &
        'EACH ENDS OF THE ',a1,' AXIS.')
   !     ---------------------------------------------------------------
   return
end subroutine set_gemHgrid3
