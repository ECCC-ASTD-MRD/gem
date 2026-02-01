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
      subroutine INs_Nes_Data ( F_datev )
      use, intrinsic :: iso_fortran_env
      use iso_c_binding
      use omp_timing
      use IOs
      use INs
      implicit none

      character(len=*), intent(IN ) :: F_datev

      logical, external :: INs_open
      
      character(len=256) :: component_S
      character(len=16 ) :: dateV_S
      logical :: process_L
      integer :: colors(2),COMMs(2),ierr,grid_info(100),tag,client,dim
!     
!--------------------------------------------------------------------
!
      call clock ( Lun_out, 'NES_FST+GZ', .false. )
      
      call gtmg_start ( 11, 'NES_FST+GZ', 1)
      call INs_Nes_openFST ( F_datev, ierr )
      call gtmg_stop ( 11 )
      
      if (Lun_out>0) call clock ( Lun_out, 'NES Read+Hint', .false. )
      
      call gtmg_start ( 12, 'NES_Read', 1)
      dim= ubound(NES_DATA,1)
      SRL(:)%vname(2) = ''
      INs_nplans=0
      call INs_read (NES_DATA, Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4, dim )
      call gtmg_stop ( 12 )
      
      call gtmg_start ( 13, 'NES_SEND', 1)
      call INs_Nes_send (GZ,NES_DATA,1-G_halox,G_ni+G_halox,1-G_haloy,&
                         G_nj+G_haloy,ierr)
      call gtmg_stop ( 13 )
      
      call INs_closeFST (ierr)
!     
!--------------------------------------------------------------------
!
      return
      end subroutine INs_Nes_Data
