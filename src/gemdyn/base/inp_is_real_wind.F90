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

! Remove this function by 2021 since there will be no more archive left with UT1 in m/s

!**s/r inp_is_real_wind - Detect if UT1 input is real wind

integer function inp_is_real_wind ( F_u, F_n, F_var_S) result(status)

   implicit none
#include <arch_specific.hf>

   integer :: F_n
   real, dimension(F_n) :: F_u
   character(len=*) :: F_var_S

   ! Local variables

   integer :: i
   ! Largest real wind in the atmosphere is of the order of 1000 m/s
   ! image wind = real wind * cos(lat) / rayt
   ! Therefore lastest image wind is 1.0E+03 / 1.0E+06 = 1.0E-03
   real, parameter :: max_image_wind = 1.0E-03

   status = -1

   if( trim(F_var_S) /= 'UT1' .and. trim(F_var_S) /= 'VT1' ) return

   do i = 1, F_n
      if( abs( F_u(i) ) > max_image_wind )then
         status = 1
         return
      end if
   end do

   return

end function inp_is_real_wind
