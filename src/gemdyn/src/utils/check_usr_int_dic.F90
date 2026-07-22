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
subroutine check_usr_int_dic()

  use utl
  
  implicit none
  character(len=16), dimension(:), pointer :: list_S
  character(len=1000) :: file_S
  type(UTL_usr_int_info), dimension(:), pointer :: ff
  logical :: verbose_L=.true.
  nullify(list_S,ff)
  
  file_S="/home/apm000/ords/datafiles/constants/user_output/debug/interpolation_dictionary.txt"

  ! Try to read user interpolation dirctonary
  if(.not.UTL_read_usr_int_dic(list_S,file_S,F_verbose_L=verbose_L))then
     print*,'Error reading user interpolation dirctonary'
     error stop 1
  endif

  allocate(ff(size(list_S)))

  ! Try to crack user interpolation dirctonary
  if(.not.UTL_crack_usr_int_dic(ff,list_S,file_S,F_verbose_L=verbose_L))then
     print*,'Error cracking user interpolation dirctonary'
     error stop 1
  endif
  print*,'Interpolation dictionary is OK'
 
end subroutine check_usr_int_dic
