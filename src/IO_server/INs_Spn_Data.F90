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
      subroutine INs_Spn_Data ( F_datev )
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
      call clock ( Lun_out, 'SPN_FST', .false. )
      call gtmg_start ( 41, 'SPN_FST', 1)
      call INs_Spn_openFST ( F_datev, ierr )
      call gtmg_stop ( 41 )
      if (ierr < 0) then
         if (Lun_out>0) then
            write (Lun_out,'(/" FOUND NO SPN data valid at ",&
                             a/," ABORT ABORT ABORT"//)') F_datev
         endif
         stop
      endif

      if (Lun_out>0) call clock ( Lun_out, 'SPN Read+Hint', .false. )
      
      call gtmg_start ( 42, 'SPN_Read', 1)
      dim= ubound(SPN_DATA,1)
      SRL(:)%vname(2) = ''
      INs_nplans=0

      if (IOS_couleur==0) then
         call INs_read (SPN_DATA(1,1),Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4, dim )
      endif
      if (IOS_couleur==1) then
         call INs_read (SPN_DATA(1,2),RotY_ig1, RotY_ig2, RotY_ig3, RotY_ig4, dim )
      endif
      call gtmg_stop ( 42 )
      
      call gtmg_start ( 43, 'SPN_SEND', 1)
      call INs_Spn_send (SPN_DATA(1,IOS_couleur+1),1-G_halox,G_ni+G_halox,1-G_haloy,&
                         G_nj+G_haloy,ierr)
      call gtmg_stop ( 43 )

      un= SPN_file%get_unit()
      success = SPN_file%close()
      if (Lun_out>0) then
         write (Lun_out,'(" Fortran unit:",i6," is closed")') un
      endif
      
      call INs_deallocate ()

 9000 format(/,' TREATING INPUT ',a,' DATA VALID AT: ',a,&
             /,' ===============================================')
!     
!--------------------------------------------------------------------
!
      return
      end subroutine INs_Spn_Data
