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
      subroutine INs_Iau_Data ( F_datev )
      use, intrinsic :: iso_fortran_env
      use iso_c_binding
      use omp_timing
      use IOs
      use INs
      implicit none

      character(len=*), intent(IN ) :: F_datev

      logical :: success
      integer :: un,dim,ierr
!     
!--------------------------------------------------------------------
!
      call clock ( Lun_out, 'IAU_FST', .false. )
      call gtmg_start ( 21, 'IAU_FST', 1)
      call INs_Iau_openFST ( F_datev, ierr )
      call gtmg_stop ( 21 )
      
      if (Lun_out>0) call clock ( Lun_out, 'IAU Read+Hint', .false. )
      
      call gtmg_start ( 22, 'IAU_Read', 1)
      dim= ubound(IAU_DATA,1)
      SRL(:)%vname(2) = ''
      INs_nplans=0
      if (IOS_couleur==0) then
         call INs_read (IAU_DATA(1,1),Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4, dim )
      endif
      if (IOS_couleur==1) then
         call INs_read (IAU_DATA(1,2),RotY_ig1, RotY_ig2, RotY_ig3, RotY_ig4, dim )
      endif
      call gtmg_stop ( 22 )
      
      call gtmg_start ( 23, 'IAU_SEND', 1)
      call INs_Iau_send (IAU_DATA(1,IOS_couleur+1),1-G_halox,G_ni+G_halox,1-G_haloy,&
                         G_nj+G_haloy,ierr)
      call gtmg_stop ( 23 )

      un= IAU_file%get_unit()
      success = IAU_file%close()
      if (Lun_out>0) then
         write (Lun_out,'(" Fortran unit:",i6," is closed")') un
      endif
      
      call INs_deallocate ()
!     
!--------------------------------------------------------------------
!
      return
      end subroutine INs_Iau_Data
