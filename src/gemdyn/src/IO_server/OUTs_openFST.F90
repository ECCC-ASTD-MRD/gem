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

      subroutine OUTs_openFST (F_me_S)
      use iso_c_binding
      use OUTs
      use clib_itf_mod
      use rmn_fst24
      implicit none
      
#include <rmnlib_basics.hf>
      logical success
      character(len=*) :: F_me_S
      integer :: i,indx,err

!     
!--------------------------------------------------------------------
!
!      err = fstopc('MSGLVL','INFORM',RMN_OPT_SET)

      OUTs_out(:)%fnt = 0
      if (OUTs_1o1_L) print*, nout_files,'FST files for Out_npas= ',Out_npas

      err = clib_mkdir(trim(Out_dirname_S)//'/'//trim(F_me_S))
      do i=1,nout_files
         indx = index(Out_filenames_S(i),"_")
         if (indx>5) then
            OUTs_out(i)%type= Out_filenames_S(i)(1:2)
         else
            OUTs_out(i)%type= Out_filenames_S(i)(1:indx-1)
         endif
         ! Put user grid output in specific directory without prefix.
         if(  OUTs_out(i)%type(1:2) == 'dm' .or. &
              OUTs_out(i)%type(1:2) == 'dp' .or. &
              OUTs_out(i)%type(1:2) == 'pm' .or. &
              OUTs_out(i)%type(1:2) == 'pp' .or. &
              OUTs_out(i)%type(1:2) == 'dh' .or. &
              OUTs_out(i)%type(1:2) == 'ph' ) then
            OUTs_out(i)%name= trim(Out_dirname_S)//'/'//trim(F_me_S)//'/'//trim(Out_filenames_S(i))            
         else
            err = clib_mkdir(trim(Out_dirname_S)//'/usr'//Out_filenames_S(i)(2:2))
            err = clib_mkdir(trim(Out_dirname_S)//'/usr'//Out_filenames_S(i)(2:2)//'/'//trim(F_me_S))
            OUTs_out(i)%name=trim(Out_dirname_S)//'/usr'//Out_filenames_S(i)(2:2)//'/'//trim(F_me_S)//'/'//trim(Out_filenames_S(i)(3:))
         endif
         success=Out_list_files(i)%open(trim(OUTs_out(i)%name),'STD+RND+R/W')
         OUTs_out(i)%fnt=Out_list_files(i)%get_unit()
         if (OUTs_1o1_L) print*, 'FST: ', OUTs_out(i)%fnt,trim(OUTs_out(i)%name)
      end do
!     
!--------------------------------------------------------------------
!
      return
      end subroutine OUTs_openFST
