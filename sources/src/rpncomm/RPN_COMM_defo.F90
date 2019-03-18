!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! !                          Environnement Canada
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
!InTf!
!
! this routine sets the default communicator for some collective operations
!   halo exchanges, adj halo exchanges, RPN_COMM_diag, RPN_COMM_bcst_world,
!   RPN_COMM_dist, RPN_COMM_coll, RPN_COMM_globalsum,
!   RPN_COMM_tmg_wrt, RPN_COMM_move
!
SUBROUTINE RPN_COMM_defo(com)             !InTf!
  use rpncomm_com
  implicit none                           !InTf!
  character(len=*), intent(IN) ::  com    !InTf!
  integer comm
  integer, external :: rpn_comm_comm

  if(.not. associated(com_tab)) call init_com_tab   ! in case communicator table is not initialized yet
  comm=rpn_comm_comm(com)               ! get communicator
  pe_defcomm = comm                     ! set default communicator
  com_tab(defcom_index)%number = comm   ! do not forget to update the communicator symbol table

  return
end SUBROUTINE RPN_COMM_defo              !InTf!
