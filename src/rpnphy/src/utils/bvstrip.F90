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

!/@*
function bvstrip(F_name) result(F_norm)
  use phymem, only: PHY_NAMELEN
  implicit none
  !@objective Strip bus variable names of extraneous information
  !@arguments
  character(len=*), intent(in) :: F_name        !Bus variable name
  !@return
  character(len=PHY_NAMELEN) :: F_norm          !Normalized name
  !@author Ron McTaggart-Cowan, 2014-03
  !@description
  integer :: i
  character(len=1), dimension(2), parameter :: ext = (/'w','W'/)
  !*@/
  F_norm = F_name
  if (len_trim(F_norm) < 2) return  !minimum length for separator
  do i=1,size(ext)
     if (F_name(len_trim(F_name)-1:) == ','//ext(i)) &
          F_norm = F_name(:len_trim(F_name)-2)
  enddo
  return
end function bvstrip
