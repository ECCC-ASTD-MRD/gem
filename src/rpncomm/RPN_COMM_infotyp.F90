!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 2021  Division de Recherche en Prevision Numerique
! !                     Environnement Canada
! !
! ! This library is free software; you can redistribute it and/or
! ! modify it under the terms of the GNU Lesser General Public
! ! License as published by the Free Software Foundation,
! ! version 2.1 of the License.
! !
! ! This library is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! ! Lesser General Public License for more details.
! !/
integer function RPN_COMM_infotyp(info)                 !InTf!
  use rpn_comm
  implicit none
  character(len=*), intent(IN) :: info                  !InTf!

  if(trim(info) == 'MPI_INFO_NULL') then
    RPN_COMM_infotyp = MPI_INFO_NULL
  else
    write(rpn_u,*) 'ERROR: Unknown info type ',info
  endif

  return
end function RPN_COMM_infotyp                           !InTf!
