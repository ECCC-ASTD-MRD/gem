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
      subroutine OUTs_closFST (F_step)
      use iso_c_binding
      use OUTs
      use INs
      use clib_itf_mod
      use rmn_fst24
      implicit none

      integer, intent(IN) :: F_step
      
#include <rmnlib_basics.hf>

      character(len=2048) :: filen, link, step
      integer :: i,err
      logical :: success
!     
!--------------------------------------------------------------------
!
      if (OUTs_1o1_L) print*, 'Closing FST files unit', &
                             (OUTs_out(i)%fnt,i=1,nout_files)
      do i=1,nout_files
         success = Out_list_files(i)%close()
         OUTs_out(i)%fnt = 0
         OUTs_out(i)%name= ''
         OUTs_out(i)%type= ''
      end do

      if (OUTs_1o1_L) then
         if (F_step>-9999) then
            write(step,'(i10.10)') F_step
            filen= trim(Out_dirname_S)//'/../'//'output_ready_server_'//trim(step)
            link = trim(Out_dirname_S)//'/../'//'output_ready'
            err = clib_symlink ( trim(filen), trim(link) )
            print*, 'Launching Sortie job on step ',trim(step) 
         endif
      endif
!     
!--------------------------------------------------------------------
!
      return
      end subroutine OUTs_closFST
