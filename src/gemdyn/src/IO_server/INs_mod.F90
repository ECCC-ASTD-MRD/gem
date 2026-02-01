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

module INs
   use iso_c_binding
   use, intrinsic :: iso_fortran_env
   use vGrid_Descriptors
   use rmn_fst24
   implicit none
   public
   save

      character(len=2048) :: Path_input_S
      character(len=64), dimension(:), allocatable :: INs_list_S
      logical :: INs_Iau_listfst_L, INs_Spn_listfst_L
      integer :: INs_GEM_COMM,INs_me1o1,INs_gem1o1,INs_host
      integer :: GZ_win,NEST_win,IAU_win,SPN_win,INs_hostmyproc
      integer :: INs_CMCdate0
      logical :: INs_server_L=.false. , INs_1o1_L
      integer :: INs_nia,INs_nja,INs_nka,INs_nid,INs_njd,INs_hord
      integer :: INs_Nes_nkd,INs_Iau_nkd,INs_Spn_nkd
       
      real, dimension(:,:,:), pointer :: GZ
      real, dimension(:,:  ), pointer :: IAU_DATA, SPN_DATA
      real, dimension(:    ), pointer :: NES_data
      real, dimension(:,:  ), allocatable :: GZbuf, &
                               Nesbuf, Iaubuf, Spnbuf
      integer, dimension(:), allocatable :: DIP1, GZIP1, &
                              Nes_iBUF,Iau_iBUF,Spn_iBUF,&
                 INs_Nes_isend,INs_Iau_isend,INs_Spn_isend
      character(len=32), dimension(:), allocatable :: &
                             Nes_cBUF,Iau_cBUF,Spn_cBUF
      integer :: GZcaract(5),INs_maxreqs,INs_maxNKA 
      real(kind=REAL64), pointer :: vtbl_8(:,:,:)
      real(kind=REAL64), pointer :: Nes_VGD_tbl(:),&
                     Iau_VGD_tbl(:),Spn_VGD_tbl(:)
      
      character(len=16) :: Inp_datev, INs_runstrt_S
      logical Inp_src_hauteur_L, SRC_GZ_L, INs_vgd_L
      
      integer Inp_nfiles, Inp_kind, Inp_vgdkind, Ins_handle, &
              Inp_cmcdate, INs_nreq, Inp_rtag, INs_nplans, INs_n123(3)
      type(fst_file), dimension(:), contiguous,pointer :: Inp_list_files => null()
      type(vgrid_descriptor) :: Inp_vgd_src
      type(fst_file) :: Inp_file,IAU_file,SPN_file

      type :: REQ
         character(len=32) :: vname(2)
         character(len=4 ) :: stag
         character(len=1 ) :: src
         integer :: unf,nk,deb
      end type REQ
      type(REQ),dimension (:), allocatable :: SRL
      
      type :: STD
         character(len=1) :: type
         integer :: fnt
      end type STD
      type(STD) :: fst(10)

contains

      subroutine INs_getvgd ()
      use ios
      implicit none
      
      integer err
!
!-----------------------------------------------------------------------
!
      if (associated(vtbl_8)) deallocate (vtbl_8)
      INs_vgd_L= .false. ; nullify (vtbl_8)
      if (lun_out>0) print*, 'Attempting vgd_new on unit: ',Ins_handle
      if (vgd_new ( Inp_vgd_src, unit=Ins_handle, &
                    format='fst', ip1=-1, ip2=-1, &
                    quiet=.true. )==VGD_OK) then       
         INs_vgd_L = vgd_get ( Inp_vgd_src, key='KIND',&
                value=Inp_vgdkind,quiet=.true. )==VGD_OK
         err= vgd_get ( Inp_vgd_src, 'VTBL', vtbl_8)
         INs_n123(1:3) = ubound(vtbl_8)
      endif
      if (INs_vgd_L .and. (lun_out>0)) &
      print*, 'Inp_vgd_src and vtbl_8 from unit: ',Ins_handle
!
!-----------------------------------------------------------------------
!
      return
      end subroutine INs_getvgd
      
end module INs
