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

!**s/r inp_close - Close all fst input files already opened

      subroutine inp_close ()
      use inp_mod
      implicit none
#include <arch_specific.hf>

#include <rmnlib_basics.hf>

      integer i, err_code
!
!--------------------------------------------------------------------
!
      err_code= 0

      if (Inp_iome >= 0) then
         do i= 1, Inp_nfiles
            if (.not. Inp_list_files(i)%close()) then
               err_code= -1
            end if
         end do
         deallocate (Inp_list_files) ; nullify (Inp_list_files)
      end if
      err_code = vgd_free(Inp_vgd_src)

      call gem_error ( err_code, 'inp_close', &
                       'Problems opening input files' )
!
!--------------------------------------------------------------------
!
      return
      end
