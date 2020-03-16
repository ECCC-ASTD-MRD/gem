
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
module rpncomm_com
  use rpn_comm
  implicit none
  integer, parameter :: MAX_COMM_TAB=128
  type(symtab), dimension(:), pointer, save :: com_tab => NULL()  ! communicator translation table
  integer, save :: defcom_index=-1                                ! index for RPN_COMM_DEFAULT
  integer, save :: max_com_index=0                    
contains
!
! allocate and initialize com_tab, the communicator table
!
  subroutine init_com_tab
    implicit none
    integer :: i

    if( associated(com_tab) ) return  ! job already done

    allocate(com_tab(MAX_COMM_TAB+1))   ! allocate table and initialize "the usual suspects"
    com_tab( 1) = symtab(pe_grid_peers,RPN_COMM_GRIDPEERS)
    com_tab( 2) = symtab(pe_indomm,RPN_COMM_GRID)
    com_tab( 3) = symtab(pe_indomm,RPN_COMM_DOMM)
    com_tab( 4) = symtab(WORLD_COMM_MPI,RPN_COMM_WORLD)
    com_tab( 5) = symtab(WORLD_COMM_MPI,RPN_COMM_ALLDOMAINS)
    com_tab( 6) = symtab(pe_a_domain,RPN_COMM_ALLGRIDS)
    com_tab( 7) = symtab(pe_indomms,RPN_COMM_MULTIGRID)
    com_tab( 8) = symtab(pe_wcomm,RPN_COMM_ALL)
    com_tab( 9) = symtab(pe_defcomm,RPN_COMM_DEFAULT)  ; defcom_index = 9  ! this item might get updated by rpn_comm_defo
    com_tab(10) = symtab(pe_myrow,RPN_COMM_EW)
    com_tab(11) = symtab(pe_mycol,RPN_COMM_NS)
    com_tab(12) = symtab(pe_blocmaster,RPN_COMM_BLOCMASTER)
    com_tab(13) = symtab(pe_bloc,RPN_COMM_BLOCK)
    com_tab(14) = symtab(MPI_COMM_WORLD,RPN_COMM_UNIVERSE)
    com_tab(15) = symtab(MPI_COMM_NULL,RPN_COMM_NULL)
    max_com_index = 15
    do i = 16,MAX_COMM_TAB+1
      com_tab(i) = symtab(MPI_COMM_NULL,"")
    enddo
!print *,'DEBUG: com_tab initialized'
!print *,com_tab(1:15)
  end subroutine init_com_tab
!
! add new communicator to table
!
  function new_com(com,mpicom) result(indx)
    implicit none
    character(len=*), intent(IN) :: com      ! rpn_comm character style communicator
    integer, intent(IN) :: mpicom            ! MPI integer communicator
    integer :: indx

    if(.not. associated(com_tab)) call init_com_tab
    if(max_com_index < MAX_COMM_TAB) then
      max_com_index = max_com_index + 1
      com_tab(max_com_index)%number = mpicom
      com_tab(max_com_index)%string = trim(com)
      indx = max_com_index
    else
      indx = -1
    endif
    return
  end function new_com
!
! get index of string com in com_tab
!
  function indx_com_tab(com) result(indx)
    implicit none
    character(len=*), intent(IN) :: com
    integer :: indx
    integer :: i

    if(.not. associated(com_tab)) call init_com_tab
    indx = MAX_COMM_TAB+1  ! will be returned if com is not found
    do i = 1,max_com_index
      if( trim(com_tab(i)%string) == trim(com) ) then
        indx = i
!print *,'DEBUG: index of "'//trim(com)//'" is',indx
        return
      endif
    enddo
  end function indx_com_tab
end module rpncomm_com
!InTf!
      integer function RPN_COMM_comm(com)                    !InTf!
!	Luc Corbeil, 2000-11-21
!
!	lien entre chaine de caractere de communicateur
!	GRID, EW et NS et leur numero assigne par
!	MPI.
!
      use rpncomm_com
      implicit none                                 !InTf!
!      include mpif.h
!        include rpn_comm.h
      character(len=*), intent(IN) :: com           !InTf!
      character(len=32) comm
      integer :: i, indx

      if(.not. associated(com_tab)) call init_com_tab
      call rpn_comm_low2up(com,comm)
!print *,'DEBUG: avant indx_com_tab'
      indx = indx_com_tab(comm)
!print *,'DEBUG: indx_com_tab=',indx
      RPN_COMM_comm = com_tab(indx_com_tab(comm))%number
!print *,'DEBUG: number=',RPN_COMM_comm

      return
#if defined(DEPRECATED)
      RPN_COMM_comm = MPI_COMM_NULL

      if (trim(comm) == RPN_COMM_GRIDPEERS) then
         RPN_COMM_comm=pe_grid_peers
         return
      endif
      if (trim(comm) == RPN_COMM_GRID) then
         RPN_COMM_comm=pe_indomm  ! alias pe_grid
         return
      endif
      if(trim(comm) == RPN_COMM_DOMM) then
         RPN_COMM_comm=pe_indomm  ! alias pe_grid
         return
      endif
      if (trim(comm) == RPN_COMM_WORLD) then
         RPN_COMM_comm=WORLD_COMM_MPI  ! alias pe_all_domains
         return
      endif
      if(trim(comm) == RPN_COMM_ALLDOMAINS) then
         RPN_COMM_comm=WORLD_COMM_MPI  ! alias pe_all_domains
         return
      endif
      if (trim(comm) == RPN_COMM_ALLGRIDS) then
         RPN_COMM_comm=pe_a_domain
         return
      endif
      if (trim(comm) == RPN_COMM_MULTIGRID) then
         RPN_COMM_comm=pe_indomms ! alias pe_multi_grid
         return
      endif
      if (trim(comm) == RPN_COMM_ALL) then
        RPN_COMM_comm=pe_wcomm  ! alias pe_grid
        return
      endif
      if(trim(comm) == RPN_COMM_DEFAULT) then
        RPN_COMM_comm=pe_defcomm
        return
      endif
      if (trim(comm) == RPN_COMM_EW) then
         RPN_COMM_comm=pe_myrow
         return
      endif
      if (trim(comm) == RPN_COMM_NS) then
         RPN_COMM_comm=pe_mycol
         return
      endif
      if (trim(comm) == RPN_COMM_BLOCMASTER) then
         RPN_COMM_comm=pe_blocmaster
         return
      endif
      if (trim(comm) == RPN_COMM_BLOCK) then
         RPN_COMM_comm=pe_bloc
         return
      endif
      if (trim(comm) == RPN_COMM_UNIVERSE) then
         RPN_COMM_comm=MPI_COMM_WORLD
         return
      endif
      if (trim(comm) == RPN_COMM_NULL) then
         RPN_COMM_comm=MPI_COMM_NULL
         return
      endif

      write(rpn_u,*) 'Unknown communicator ',comm,', aborting'
      stop
#endif
      end function RPN_COMM_comm                                  !InTf!
!InTf!
      integer function RPN_COMM_custom_comm(mpicom,name,mode)     !InTf!
      use rpncomm_com
      implicit none                                               !InTf!
!     lookup, create, or delete a custom communicator with a rpn_comm style name
      character(len=*), intent(IN) :: name                        !InTf!
      integer, intent(IN) :: mpicom                               !InTf!
      integer, intent(IN) :: mode                                 !InTf!
!
      integer :: i
      character (len=32) :: name2
      integer, save :: base = 0

      if(.not. associated(com_tab)) call init_com_tab
      if(base == 0) base = max_com_index

      RPN_COMM_custom_comm = MPI_COMM_NULL
      call rpn_comm_low2up(name,name2)

      if(mode==RPN_COMM_GET) then                   ! look for rpn_comm communicator named "name"
        RPN_COMM_custom_comm = com_tab(indx_com_tab(name2))%number
      else if(mode==RPN_COMM_SET) then             ! add "name" and com to the rpn_comm communicators
        i = new_com(name2,mpicom)
        if(i>0) RPN_COMM_custom_comm = mpicom      ! if create was successful
      else if(mode==RPN_COMM_DEL) then             ! delete "name" and com from rpn_comm communicators
        i = indx_com_tab(name2)                    ! find name
        if(i>0 .and. i>base) then                  ! if found and not a permanent communicator
          com_tab(i)%string = ''
          com_tab(i)%number = MPI_COMM_NULL
        endif
        RPN_COMM_custom_comm = MPI_COMM_NULL
      else
        write(rpn_u,*) 'ERROR: RPN_COMM_custom_comm illegal mode'
        RPN_COMM_custom_comm = MPI_COMM_NULL
      endif
      return

#if defined(DEPRECATED)
      integer, parameter :: MAX_NAMES=128
      type(SYMTAB), save, dimension(:), pointer :: names => NULL()
      integer, save :: entries=0
!
      if(.not. associated(names)) allocate(names(MAX_NAMES))

      RPN_COMM_custom_comm=MPI_COMM_NULL
      name2 = trim(name)
!
      if(mode==RPN_COMM_GET) then                 ! look for rpn_comm communicator named "name"
         do i = 1 , entries
            if(trim(name2)==trim(names(i)%string)) then
               RPN_COMM_custom_comm = names(i)%number
               return
            endif
         enddo
      else if(mode==RPN_COMM_SET) then             ! add "name" and com to the rpn_comm communicators
         if(entries<MAX_NAMES) then
            entries = entries + 1
            names(entries)%string = trim(name2)
            names(entries)%number = com
            RPN_COMM_custom_comm=com
         else
            write(rpn_u,*) 'ERROR: communicator table full'
         endif
      else if(mode==RPN_COMM_DEL) then              ! delete "name" and com from rpn_comm communicators
      else
         write(rpn_u,*) 'ERROR: RPN_COMM_custom_comm illegal mode'
      endif
      return
#endif
      end function RPN_COMM_custom_comm                      !InTf!
!
!       fill an entity of type rpncomm_communicator from type string
!       ctyp_c : character version of communicator (rpn comm char style) (see RPN_COMM_constants.inc)
!       ctyp   : new item of type rpncomm_communicator
!InTf!
        subroutine RPN_COMM_i_comm(ctyp_c,ctyp,mcom)         !InTf!
        use rpn_comm
!!      import :: rpncomm_communicator                       !InTf!
        implicit none
        type(rpncomm_communicator), intent(OUT) :: ctyp      !InTf!
        character(len=*), intent(IN) :: ctyp_c               !InTf!
        integer, optional, intent(IN) :: mcom                !InTf!

        integer, external :: RPN_COMM_comm
        integer :: temp, siz, ierr

        if( present(mcom) ) then
!print *,'DEBUG: using optional mcom'
          temp = mcom
        else
!print *,'DEBUG: using "'//trim(ctyp_c)//'"'
          temp = RPN_COMM_comm(ctyp_c)           ! converted communicator value
        endif
        call MPI_comm_size(temp,siz,ierr)        ! test communicator validity
        if( ierr .ne. MPI_SUCCESS ) temp = MPI_COMM_NULL
        ctyp%p = C_LOC(WORLD_COMM_MPI)            ! signature
        ctyp%t2 = temp
        ctyp%t1 = ieor(temp,RPN_COMM_MAGIC)     ! xor with magic token
        end subroutine RPN_COMM_i_comm                    !InTf!

