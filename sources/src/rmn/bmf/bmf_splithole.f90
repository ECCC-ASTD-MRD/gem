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
subroutine bmf_splithole(holename,time1,time2,i0,isize,j0,jsize)
  use bmf_modsplit
!
! BMF_SPLITHOLE, L. Corbeil, 31/08/2002
!ARGUMENTS
! holename,       self-explanatory
! time1, time2    timestamps
! i0,j0           origin of the hole
! isize,jsize     size of the hole

  implicit none
  character(len=4), intent(IN) :: holename
  integer, intent(IN) :: i0,isize,j0,jsize,time1,time2
  integer ier

  type(bmf_holelist), pointer :: current_hole

  nullify(current_hole)
 
  allocate(current_hole, STAT=ier)
  if(ier/=0) then
     write(*,*) 'BMF_SPLITHOLE: Error memory allocation, abort'
  endif
  current_hole%bmf_trou%bmf_holename=holename
  current_hole%bmf_trou%bmf_holet1=time1
  current_hole%bmf_trou%bmf_holet2=time2
  current_hole%bmf_trou%bmf_holei0=i0
  current_hole%bmf_trou%bmf_holeisize=isize
  current_hole%bmf_trou%bmf_holej0=j0
  current_hole%bmf_trou%bmf_holejsize=jsize
  current_hole%trou_suivant => holelist

  holelist => current_hole
  
  return 
end subroutine bmf_splithole

