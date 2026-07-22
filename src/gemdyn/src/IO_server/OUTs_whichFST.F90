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

      subroutine OUTs_whichFST (F_unf, F_prefix_S)
      use iso_c_binding
      use IOs
      use OUTs
      use rmn_fst24
      implicit none

      character(len=*), intent(IN ) :: F_prefix_S
      integer         , intent(OUT) :: F_unf
      
      integer :: i
      character(len=2) :: prefix_S
!     
!--------------------------------------------------------------------
!
      do i=1,nout_files
         if(OUTs_hor_int_L(F_prefix_S))then
            ! User horizontal grid
            prefix_S=F_prefix_S(4:4)//F_prefix_S(2:2)
         else
            prefix_S=F_prefix_S(1:2)
         endif
         if (trim(prefix_S)==trim(OUTs_out(i)%type))then
            F_unf = Out_list_files(i)%get_unit()
            Out_file = Out_list_files(i)
         endif
      end do
!     
!--------------------------------------------------------------------
!
      return
      end subroutine OUTs_whichFST
