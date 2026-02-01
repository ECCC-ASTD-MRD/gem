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

!**s/r INs_deallocate - 

      subroutine INs_deallocate ()
      use INs
      use ios
      implicit none

      integer err
!
!--------------------------------------------------------------------
!
      fst(:)%fnt= -1 ; fst(:)%type= '0'

      if (associated(Inp_list_files)) then
         deallocate (Inp_list_files)
         nullify (Inp_list_files)
      endif
      if (associated(GZ)) then
         call MPI_Win_free (GZ_win  ,err)
         nullify(GZ) ; deallocate (GZIP1)
      endif
      if (associated(NES_DATA)) then
         call MPI_Win_free (NEST_win,err)
         nullify(NES_DATA) ; deallocate (DIP1)
      endif
      if (associated(IAU_DATA)) then
         call MPI_Win_free (IAU_win,err)
         nullify(IAU_DATA) ; deallocate (DIP1)
      endif
      if (associated(SPN_DATA)) then
         call MPI_Win_free (SPN_win,err)
         nullify(SPN_DATA) ; deallocate (DIP1)
      endif
      INs_nia= -1 ; INs_nja= -1 ; INs_nka= -1

      if (INs_vgd_L) then
         err= vgd_free(Inp_vgd_src)
         if ((err==0) .and. (lun_out>0)) &
         print*, 'Inp_vgd_src is now free'
         INs_vgd_L= .false.
      endif
!
!--------------------------------------------------------------------
!
      return
      end subroutine INs_deallocate
