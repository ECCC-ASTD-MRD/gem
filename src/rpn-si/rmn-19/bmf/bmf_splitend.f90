!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2005  Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
subroutine bmf_splitend
! SUBROUTINE bmf_splitend, L. Corbeil
!
! ARGUMENTS
!    
! DESCRIPTION
! Subroutine which stops bmf_split mode. This mode is used in order
! to split and write data for later use in a parallel context.
! Files opened by bmf_splitinit will be closed.  
!
  use bmf_modsplit
  implicit none
  integer i,fnom
  integer troucoor(4)
  external fnom

  type(bmf_holelist), pointer ::  current_hole

! Write holes if any...

  current_hole => holelist
  do
     if(.not.associated(current_hole)) EXIT
     troucoor(1)=current_hole%bmf_trou%bmf_holei0
     troucoor(2)=current_hole%bmf_trou%bmf_holeisize
     troucoor(3)=current_hole%bmf_trou%bmf_holej0
     troucoor(4)=current_hole%bmf_trou%bmf_holejsize
     call bmf_splitwrall(current_hole%bmf_trou%bmf_holename,&
                         4,1,1,current_hole%bmf_trou%bmf_holet1,&
                         current_hole%bmf_trou%bmf_holet2, &
                         0,0,40,0,troucoor)
     current_hole=>current_hole%trou_suivant
  enddo
  if(.not.allocated(split_files)) then
      write(*,*) 'BMF_SPLITEND: split mode not started yet: use SPLITINIT'
      stop
  endif
  do i=1,size(split_unit)
    if(split_unit(i).ne.0) then
       call fclos(split_unit(i))
    endif
  enddo
  deallocate(split_files)
  deallocate(split_unit)

  current_hole => holelist
  do
     if(.not.associated(current_hole)) EXIT
     holelist => current_hole%trou_suivant
!     write(*,*) 'avant' ,associated(holelist)
     deallocate(current_hole)
     current_hole => holelist
!           write(*,*) associated(current_hole),associated(holelist)
  enddo



  return
end subroutine bmf_splitend
