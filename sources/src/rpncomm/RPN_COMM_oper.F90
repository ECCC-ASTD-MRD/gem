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
! */
!InTf!
        integer function RPN_COMM_oper(op)              !InTf!
	use rpn_comm
!	Luc Corbeil, 2000-11-20
!	lien entre datatype et MPI_datatype

        implicit none                                   !InTf!
!        include 'mpif.h'
        character(len=*), intent(IN) :: op              !InTf!
        character(len=32) operation
        integer :: i

        call rpn_comm_low2up(op,operation)

        RPN_COMM_oper = -999999  ! precondition to error return

        do i = 1 , size(op_tab)
          if(op_tab(i)%string == operation) then
            RPN_COMM_oper = op_tab(i)%number
            return
          endif
        enddo
#if defined(OBSOLETE)
        goto 777

        if (operation(1:11).eq.'MPI_OP_NULL') then
           RPN_COMM_oper=MPI_OP_NULL
           return
        endif
        if (operation(1:7).eq.'MPI_MAX') then
           RPN_COMM_oper=MPI_MAX
           return
        endif
        if (operation(1:7).eq.'MPI_MIN') then
           RPN_COMM_oper=MPI_MIN
           return
        endif
        if (operation(1:7).eq.'MPI_SUM') then
           RPN_COMM_oper=MPI_SUM
           return
        endif
        if (operation(1:8).eq.'MPI_PROD') then
           RPN_COMM_oper=MPI_PROD
           return
        endif
        if (operation(1:8).eq.'MPI_LAND') then
           RPN_COMM_oper=MPI_LAND
           return
        endif
        if (operation(1:8).eq.'MPI_BAND') then
           RPN_COMM_oper=MPI_BAND
           return
        endif
        if (operation(1:7).eq.'MPI_LOR') then
           RPN_COMM_oper=MPI_LOR
           return
        endif
        if (operation(1:7).eq.'MPI_BOR') then
           RPN_COMM_oper=MPI_BOR
           return
        endif
        if (operation(1:8).eq.'MPI_LXOR') then
           RPN_COMM_oper=MPI_LXOR
           return
        endif
        if (operation(1:8).eq.'MPI_BXOR') then
           RPN_COMM_oper=MPI_BXOR
           return
        endif
        if (operation(1:10).eq.'MPI_MAXLOC') then
           RPN_COMM_oper=MPI_MAXLOC
           return
        endif
        if (operation(1:10).eq.'MPI_MINLOC') then
           RPN_COMM_oper=MPI_MINLOC
           return
        endif


777     continue
#endif
        write(rpn_u,*) 'Unknown operation ',op,' aborting'
        stop
          
        return
        end function RPN_COMM_oper                  !InTf!
subroutine RPN_COMM_i_oper(op,r_oper)             !InTf!
  use rpn_comm
  implicit none
!! import :: rpncomm_operator                     !InTf!
  character(len=*), intent(IN) :: op              !InTf!
  type(rpncomm_operator), intent(OUT) :: r_oper   !InTf!

  integer, external :: RPN_COMM_oper

  r_oper%t2 = RPN_COMM_oper(op)
  r_oper%t1 = ieor(r_oper%t2,RPN_COMM_MAGIC)
  r_oper%p  = C_LOC(WORLD_COMM_MPI)            ! signature

  return
end subroutine RPN_COMM_i_oper                    !InTf!
function RPN_COMM_i_valid_oper(r_oper) result (is_valid)  !InTf!
  use rpn_comm
!! import :: rpncomm_operator                     !InTf!
  implicit none
  type(rpncomm_operator), intent(IN) :: r_oper    !InTf!
  logical :: is_valid                             !InTf!

  type(C_PTR) :: temp

  temp = C_LOC(WORLD_COMM_MPI)
  is_valid = c_associated(temp,r_oper%p)
  if(.not. is_valid ) then
    write(rpn_u,*) 'ERROR: bad signature for rpncomm_operator'
    return
  endif
  is_valid = ieor(r_oper%t1,RPN_COMM_MAGIC) == r_oper%t2
  if(.not. is_valid ) then
    write(rpn_u,*) 'ERROR: bad checksum for rpncomm_operator'
  endif
  
  return
end function RPN_COMM_i_valid_oper                !InTf!
