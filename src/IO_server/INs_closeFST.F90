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

!**s/r INs_closeFST - Close all fst input files already opened

      subroutine INs_closeFST (F_err)
      use ios
      use INs
      use rmn_fst24
      implicit none

      integer, intent(OUT) :: F_err

      logical :: success
      integer i, j, unit, err_code(size(Nes_fst))
!
!--------------------------------------------------------------------
!
      f_err= 0 ; err_code=0
      do i=1, size(Nes_fst)
         unit= Nes_fst(i)%get_unit()
         if (fst(i)%fnt>0) then
            if (fst(i)%type/='R') then
               success= Nes_fst(i)%close()
               if ( success .and. (lun_out>0)) then
                  write (6,'("Fortran unit:",i6," is closed")') unit
               else
                  err_code(i)= -1
               end if
            else
               do j= 1, Inp_nfiles
                  unit= Inp_list_files(j)%get_unit()
                  success= Inp_list_files(j)%close()
                  if ( success .and. (lun_out>0)) then
                     if (lun_out>0) &
                     write (6,'("Fortran unit:",i6," is closed")') unit
                  else
                     err_code(i)= -1
                  endif
               end do
            endif
         endif
      end do
      F_err= minval(err_code)

      call INs_deallocate ()
!
!--------------------------------------------------------------------
!
      return
      end subroutine INs_closeFST
