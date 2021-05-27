!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 2020  Division de Recherche en Prevision Numerique
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
! !
! ! You should have received a copy of the GNU Lesser General Public
! ! License along with this library; if not, write to the
! ! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! ! Boston, MA 02111-1307, USA.
! !/
!! end interface                                                !InTf!
!! interface RPN_COMM_get_my_domain                             !InTf!
function RPN_COMM_get_my_domain_1(call_back) result (domain)    !InTf!
  implicit none
  external :: call_back                                         !InTf!
  integer :: domain                                             !InTf!
  integer :: temp
  call RPN_COMM_mydomain (call_back, temp)
  domain = temp
  return
end function RPN_COMM_get_my_domain_1                           !InTf!

function RPN_COMM_get_my_domain_2(ndomains, offset) result (domain) !InTf!
  implicit none
  integer, intent(IN) :: ndomains, offset                       !InTf!
  integer :: domain                                             !InTf!
  integer, save :: local_ndomains, local_offset
  integer :: temp
  local_ndomains = ndomains
  local_offset = offset
  call RPN_COMM_mydomain (internal, temp)
  domain = temp
  return
contains
  subroutine internal(ndomains, offset, err)
    integer, intent(OUT) :: ndomains, offset, err
    ndomains = local_ndomains
    offset = local_offset
    err = 0
    return
  end subroutine internal
end function RPN_COMM_get_my_domain_2                           !InTf!
!! end interface RPN_COMM_get_my_domain                         !InTf!
!! interface                                                    !InTf!
