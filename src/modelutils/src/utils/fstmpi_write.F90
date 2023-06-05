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

!/@
module fstmpi_write_mod
   use iso_c_binding
   use rpn_comm_itf_mod
   use fst_mod
   use ezgrid_mod
   use ptopo_utils
   use hgrid_wb
   implicit none
   private
   !@objective 
   !@author  Stephane Chamberland, 2012-02
   !@description
   ! Public functions
   public :: fstmpi_write
   ! Public constants
   public :: FST_NPAK_DEFAULT,FST_NPAK_FULL32
!@/

#include <rmn/msg.h>
#include <rmnlib_basics.hf>
!!!#include <arch_specific.hf>

   interface fstmpi_write
      module procedure fstmpi_write_2d_r4
      module procedure fstmpi_write_2d_r4_s
      module procedure fstmpi_write_3d_r4_ip1
      module procedure fstmpi_write_3d_r4_ip1_s
      module procedure fstmpi_write_3d_r4_vgd
      module procedure fstmpi_write_3d_r4_vgd_s
   end interface

   logical,parameter :: BCAST2ALL = .true.

contains

   !/@
   function fstmpi_write_2d_r4_s(F_fileid,F_nomvar_S,F_data,F_gridid_S,F_ip1, &
        F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3,F_typvar_S,F_etiket_S, &
        F_rewrite_L,F_writegrid_L,F_lni,F_lnj,F_ci0,F_cj0,F_cin,F_cjn,F_cdij) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      real,pointer :: F_data(:,:)
      integer,intent(in) :: F_fileid,F_ip1
      character(len=*),intent(in) :: F_nomvar_S,F_gridid_S
      integer,intent(in),optional :: F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3
      integer,intent(in),optional :: F_lni,F_lnj   !# Subset dims of F_data
      integer,intent(in),optional :: F_ci0,F_cj0,F_cin,F_cjn,F_cdij
      !# Crop the grid (collected data) to that sub domain and aggregate over dij points
      character(len=*),intent(in),optional :: F_etiket_S,F_typvar_S
      logical,intent(in),optional :: F_rewrite_L,F_writegrid_L
      !@author
      !@return
      integer :: F_istat
      !@/
      logical :: rewrite_L,writegrid_L
      integer :: ip1list(1),dateo,deet,npas,npak,dtype,ip3,nijk(3),lijk(2),istat,ci0,cj0,cin,cjn,cdij
      real,pointer :: data3d(:,:,:)
      real,target :: dummy3d(1,1,1)
      character(len=RMN_ETK_LEN) :: etiket_S
      character(len=RMN_VARTYPE_LEN) :: typvar_S
      ! ---------------------------------------------------------------------
      call msg_toall(MSG_DEBUG,'(fstmpi) write_2d_r4_s [BGN]')
      nullify(data3d)
      if (associated(F_data)) data3d => dummy3d
      F_istat = priv_init(data3d,RMN_OK,dateo,deet,npas,npak,dtype,ip3,nijk, &
           typvar_S,etiket_S,rewrite_L,writegrid_L)
      if (.not.RMN_IS_OK(F_istat)) return

      ci0=-1 ; cj0=-1 ; cin=-1 ; cjn=-1 ; cdij=1
      if (present(F_dateo)) dateo = F_dateo
      if (present(F_deet)) deet = F_deet
      if (present(F_npas)) npas = F_npas
      if (present(F_npak)) npak = F_npak
      if (present(F_dtype)) dtype = F_dtype
      if (present(F_ip3)) ip3 = F_ip3
      if (present(F_typvar_S)) typvar_S = F_typvar_S
      if (present(F_etiket_S)) etiket_S = F_etiket_S
      if (present(F_rewrite_L)) rewrite_L = F_rewrite_L
      if (present(F_writegrid_L)) writegrid_L = F_writegrid_L
      nijk(1:2) = ubound(F_data) ; lijk = lbound(F_data)
      if (present(F_lni)) nijk(1) = min(F_lni,nijk(1))
      if (present(F_lnj)) nijk(2) = min(F_lnj,nijk(2))
      if (present(F_ci0)) ci0 = F_ci0
      if (present(F_cj0)) cj0 = F_cj0
      if (present(F_cin)) cin = F_cin
      if (present(F_cjn)) cjn = F_cjn
      if (present(F_cdij)) cdij = F_cdij

      allocate(data3d(lijk(1):nijk(1),lijk(2):nijk(2),1))
      data3d(:,:,1) = F_data(:,:)
      ip1list(1) = F_ip1
      F_istat = fstmpi_write(F_fileid,F_nomvar_S,data3d,F_gridid_S,ip1list,dateo, &
           deet,npas,npak,dtype,ip3,typvar_S,etiket_S,rewrite_L,writegrid_L, &
           nijk(1),nijk(2),ci0,cj0,cin,cjn,cdij)
      deallocate(data3d,stat=istat)

      call msg_toall(MSG_DEBUG,'(fstmpi) write_2d_r4_s [END]')
      ! ---------------------------------------------------------------------
      return
   end function fstmpi_write_2d_r4_s


   !/@
   function fstmpi_write_2d_r4(F_fileid,F_nomvar_S,F_data,F_gridid,F_ip1,F_dateo, &
        F_deet,F_npas,F_npak,F_dtype,F_ip3,F_typvar_S,F_etiket_S,F_rewrite_L, &
        F_writegrid_L,F_lni,F_lnj) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      real,pointer :: F_data(:,:)
      integer,intent(in) :: F_fileid,F_gridid,F_ip1 !TODO: F_ip1 never used
      character(len=*),intent(in) :: F_nomvar_S
      integer,intent(in),optional :: F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3
      integer,intent(in),optional :: F_lni,F_lnj   !# Subset dims of F_data 
      character(len=*),intent(in),optional :: F_etiket_S,F_typvar_S
      logical,intent(in),optional :: F_rewrite_L,F_writegrid_L
      !@author
      !@return
      integer :: F_istat
      !@/
      logical :: rewrite_L,writegrid_L
      integer :: ip1list(1),dateo,deet,npas,npak,dtype,ip3,nijk(3),lijk(2),istat
      real,pointer :: data3d(:,:,:)
      real,target :: dummy3d(1,1,1)
      character(len=RMN_ETK_LEN) :: etiket_S
      character(len=RMN_VARTYPE_LEN) :: typvar_S
      ! ---------------------------------------------------------------------
      call msg_toall(MSG_DEBUG,'(fstmpi) write_2d_r4 [BGN]')
      nullify(data3d)
      if (associated(F_data)) data3d => dummy3d
      F_istat = priv_init(data3d,F_gridid,dateo,deet,npas,npak,dtype,ip3,nijk, &
           typvar_S,etiket_S,rewrite_L,writegrid_L)
      if (.not.RMN_IS_OK(F_istat)) return

      if (present(F_dateo)) dateo = F_dateo
      if (present(F_deet)) deet = F_deet
      if (present(F_npas)) npas = F_npas
      if (present(F_npak)) npak = F_npak
      if (present(F_dtype)) dtype = F_dtype
      if (present(F_ip3)) ip3 = F_ip3
      if (present(F_typvar_S)) typvar_S = F_typvar_S
      if (present(F_etiket_S)) etiket_S = F_etiket_S
      if (present(F_rewrite_L)) rewrite_L = F_rewrite_L
      if (present(F_writegrid_L)) writegrid_L = F_writegrid_L
      nijk(1:2) = ubound(F_data) ; lijk = lbound(F_data)
      if (present(F_lni)) nijk(1) = min(F_lni,nijk(1))
      if (present(F_lnj)) nijk(2) = min(F_lnj,nijk(2))

      allocate(data3d(lijk(1):nijk(1),lijk(2):nijk(2),1))
      data3d(:,:,1) = F_data(:,:)
      ip1list(1) = F_ip1
      F_istat = fstmpi_write(F_fileid,F_nomvar_S,data3d,F_gridid,ip1list,dateo, &
           deet,npas,npak,dtype,ip3,typvar_S,etiket_S,rewrite_L,writegrid_L,nijk(1),nijk(2))
      deallocate(data3d,stat=istat)

      call msg_toall(MSG_DEBUG,'(fstmpi) write_2d_r4 [END]')
      ! ---------------------------------------------------------------------
      return
   end function fstmpi_write_2d_r4


   !/@
   function fstmpi_write_3d_r4_ip1_s(F_fileid,F_nomvar_S,F_data,F_gridid_S, &
        F_ip1list,F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3,F_typvar_S, &
        F_etiket_S,F_rewrite_L,F_writegrid_L,F_lni,F_lnj,F_ci0,F_cj0, &
        F_cin,F_cjn,F_cdij) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      real,pointer :: F_data(:,:,:)
      !TODO-later: F_ip1list should be a pointer so lbound is preseved
      integer,intent(in) :: F_fileid,F_ip1list(:)
      character(len=*),intent(in) :: F_nomvar_S,F_gridid_S
      integer,intent(in),optional :: F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3
      integer,intent(in),optional :: F_lni,F_lnj   !# Subset dims of F_data
      integer,intent(in),optional :: F_ci0,F_cj0,F_cin,F_cjn,F_cdij
      !# Crop the grid (collected data) to that sub domain and aggregate over dij points
      character(len=*),intent(in),optional :: F_etiket_S,F_typvar_S
      logical,intent(in),optional :: F_rewrite_L,F_writegrid_L
      !@author
      !@return
      integer :: F_istat
      !@/
      logical :: rewrite_L,writegrid_L
      integer :: istat,dateo,deet,npas,npak,dtype,ip3,nijk(3),ci0,cj0,cin,cjn,cdij
      real,pointer :: data(:,:,:)
      character(len=RMN_ETK_LEN) :: etiket_S
      character(len=RMN_VARTYPE_LEN) :: typvar_S
      character(len=32) :: gridid_S
      ! ---------------------------------------------------------------------
      call msg_toall(MSG_DEBUG,'(fstmpi) write_3d_r4_ip1_s [BGN]')
      F_istat = priv_init(F_data,RMN_OK,dateo,deet,npas,npak,dtype,ip3,nijk, &
           typvar_S,etiket_S,rewrite_L,writegrid_L)
      if (.not.RMN_IS_OK(F_istat)) return

      ci0=-1 ; cj0=-1 ; cin=-1 ; cjn=-1 ; cdij=1
      if (present(F_dateo)) dateo = F_dateo
      if (present(F_deet)) deet = F_deet
      if (present(F_npas)) npas = F_npas
      if (present(F_npak)) npak = F_npak
      if (present(F_dtype)) dtype = F_dtype
      if (present(F_ip3)) ip3 = F_ip3
      if (present(F_typvar_S)) typvar_S = F_typvar_S
      if (present(F_etiket_S)) etiket_S = F_etiket_S
      if (present(F_rewrite_L)) rewrite_L = F_rewrite_L
      if (present(F_writegrid_L)) writegrid_L = F_writegrid_L
      if (present(F_lni)) nijk(1) = min(F_lni,nijk(1))
      if (present(F_lnj)) nijk(2) = min(F_lnj,nijk(2))
      if (present(F_ci0)) ci0 = F_ci0
      if (present(F_cj0)) cj0 = F_cj0
      if (present(F_cin)) cin = F_cin
      if (present(F_cjn)) cjn = F_cjn
      if (present(F_cdij)) cdij = F_cdij


      F_istat = priv_prep_s(data,gridid_S,F_data,F_gridid_S,nijk,ci0,cj0,cin,cjn,cdij)
      if (.not.RMN_IS_OK(F_istat)) return

      if (ptopo_isblocmaster_L) then
         F_istat = RMN_ERR
         if (associated(data).and.len_trim(gridid_S)>0) then
            F_istat = RMN_OK
            F_istat = fst_write(F_fileid,F_nomvar_S,data,gridid_S,F_ip1list,dateo, &
                 deet,npas,npak,dtype,ip3,typvar_S,etiket_S,rewrite_L,writegrid_L)
         endif
      endif
      call collect_error(F_istat) !#TODO: comment

      if (associated(data)) deallocate(data,stat=istat)
      call msg_toall(MSG_DEBUG,'(fstmpi) write_3d_r4_ip1_s [END]')
      ! ---------------------------------------------------------------------
      return
   end function fstmpi_write_3d_r4_ip1_s


   !/@
   function fstmpi_write_3d_r4_ip1(F_fileid,F_nomvar_S,F_data,F_gridid,F_ip1list, &
        F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3,F_typvar_S,F_etiket_S, &
        F_rewrite_L,F_writegrid_L,F_lni,F_lnj) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      real,pointer :: F_data(:,:,:)
      !TODO-later: F_ip1list should be a pointer so lbound is preseved
      integer,intent(in) :: F_fileid,F_gridid,F_ip1list(:)
      character(len=*),intent(in) :: F_nomvar_S
      integer,intent(in),optional :: F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3,F_lni,F_lnj
      character(len=*),intent(in),optional :: F_etiket_S,F_typvar_S
      logical,intent(in),optional :: F_rewrite_L,F_writegrid_L
      !@author
      !@return
      integer :: F_istat
      !@/
      logical :: rewrite_L,writegrid_L
      integer :: istat,dateo,deet,npas,npak,dtype,ip3,gridid,nijk(3),i0,j0
      real,pointer :: data(:,:,:)
      character(len=RMN_ETK_LEN) :: etiket_S
      character(len=RMN_VARTYPE_LEN) :: typvar_S
      ! ---------------------------------------------------------------------
      call msg_toall(MSG_DEBUG,'(fstmpi) write_3d_r4_ip1 a[BGN]')
      F_istat = priv_init(F_data,F_gridid,dateo,deet,npas,npak,dtype,ip3,nijk, &
           typvar_S,etiket_S,rewrite_L,writegrid_L)
      if (.not.RMN_IS_OK(F_istat)) return

      if (present(F_dateo)) dateo = F_dateo
      if (present(F_deet)) deet = F_deet
      if (present(F_npas)) npas = F_npas
      if (present(F_npak)) npak = F_npak
      if (present(F_dtype)) dtype = F_dtype
      if (present(F_ip3)) ip3 = F_ip3
      if (present(F_typvar_S)) typvar_S = F_typvar_S
      if (present(F_etiket_S)) etiket_S = F_etiket_S
      if (present(F_rewrite_L)) rewrite_L = F_rewrite_L
      if (present(F_writegrid_L)) writegrid_L = F_writegrid_L
      if (present(F_lni)) nijk(1) = min(F_lni,nijk(1))
      if (present(F_lnj)) nijk(2) = min(F_lnj,nijk(2))


      F_istat = priv_prep(data,gridid,F_data,F_gridid,nijk,i0,j0)
      if (.not.RMN_IS_OK(F_istat)) return

      if (ptopo_isblocmaster_L) then
         F_istat = RMN_ERR
         if (associated(data).and.RMN_IS_OK(gridid)) then
            F_istat = fst_write(F_fileid,F_nomvar_S,data,gridid,F_ip1list,dateo, &
                 deet,npas,npak,dtype,ip3,typvar_S,etiket_S,rewrite_L,writegrid_L)
         endif
      endif
      !call collect_error(F_istat)

      if (associated(data)) deallocate(data,stat=istat)
      call msg_toall(MSG_DEBUG,'(fstmpi) write_3d_r4_ip1 [END]')
      ! ---------------------------------------------------------------------
      return
   end function fstmpi_write_3d_r4_ip1


   !/@
   function fstmpi_write_3d_r4_vgd_s(F_fileid,F_nomvar_S,F_data,F_gridid_S, &
        F_vgridid,F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3,F_typvar_S, &
        F_etiket_S,F_rewrite_L,F_writegrid_L,F_lni,F_lnj,F_ci0,F_cj0, &
        F_cin,F_cjn,F_cdij) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      real,pointer :: F_data(:,:,:)
      integer,intent(in) :: F_fileid,F_vgridid
      character(len=*),intent(in) :: F_nomvar_S,F_gridid_S
      integer,intent(in),optional :: F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3
      integer,intent(in),optional :: F_lni,F_lnj   !# Subset dims of F_data
      integer,intent(in),optional :: F_ci0,F_cj0,F_cin,F_cjn,F_cdij
      !# Crop the grid (collected data) to that sub domain and aggregate over dij points
      character(len=*),intent(in),optional :: F_etiket_S,F_typvar_S
      logical,intent(in),optional :: F_rewrite_L,F_writegrid_L
      !@author
      !@return
      integer :: F_istat
      !@/
      logical :: rewrite_L,writegrid_L
      integer :: istat,dateo,deet,npas,npak,dtype,ip3,nijk(3),ci0,cj0,cin,cjn,cdij
      real,pointer :: data(:,:,:)
      character(len=RMN_ETK_LEN) :: etiket_S
      character(len=RMN_VARTYPE_LEN) :: typvar_S
      character(len=32) :: gridid_S
      ! ---------------------------------------------------------------------
      call msg_toall(MSG_DEBUG,'(fstmpi) write_3d_r4_vgd_s [BGN]')
      F_istat = priv_init(F_data,RMN_OK,dateo,deet,npas,npak,dtype,ip3,nijk, &
           typvar_S,etiket_S,rewrite_L,writegrid_L)
      if (.not.RMN_IS_OK(F_istat)) return

      ci0=-1 ; cj0=-1 ; cin=-1 ; cjn=-1 ; cdij=1
      if (present(F_dateo)) dateo = F_dateo
      if (present(F_deet)) deet = F_deet
      if (present(F_npas)) npas = F_npas
      if (present(F_npak)) npak = F_npak
      if (present(F_dtype)) dtype = F_dtype
      if (present(F_ip3)) ip3 = F_ip3
      if (present(F_typvar_S)) typvar_S = F_typvar_S
      if (present(F_etiket_S)) etiket_S = F_etiket_S
      if (present(F_rewrite_L)) rewrite_L = F_rewrite_L
      if (present(F_writegrid_L)) writegrid_L = F_writegrid_L
      if (present(F_lni)) nijk(1) = min(F_lni,nijk(1))
      if (present(F_lnj)) nijk(2) = min(F_lnj,nijk(2))
      if (present(F_ci0)) ci0 = F_ci0
      if (present(F_cj0)) cj0 = F_cj0
      if (present(F_cin)) cin = F_cin
      if (present(F_cjn)) cjn = F_cjn
      if (present(F_cdij)) cdij = F_cdij

      F_istat = priv_prep_s(data,gridid_S,F_data,F_gridid_S,nijk,ci0,cj0,cin,cjn,cdij)
      if (.not.RMN_IS_OK(F_istat)) return

      if (ptopo_isblocmaster_L) then
         F_istat = RMN_ERR
         if (associated(data).and.len_trim(gridid_S)>0) then
            F_istat = fst_write(F_fileid,F_nomvar_S,data,gridid_S,F_vgridid, &
                 dateo,deet,npas,npak,dtype,ip3,typvar_S,etiket_S,rewrite_L,writegrid_L)
         endif
      endif
      !call collect_error(F_istat)

      if (associated(data)) deallocate(data,stat=istat)
      call msg_toall(MSG_DEBUG,'(fstmpi) write_3d_r4_vgd_s [END]')
      ! ---------------------------------------------------------------------
      return
   end function fstmpi_write_3d_r4_vgd_s


   !/@
   function fstmpi_write_3d_r4_vgd(F_fileid,F_nomvar_S,F_data,F_gridid, &
        F_vgridid,F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3,F_typvar_S, &
        F_etiket_S,F_rewrite_L,F_writegrid_L,F_lni,F_lnj) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      real,pointer :: F_data(:,:,:)
      integer,intent(in) :: F_fileid,F_gridid,F_vgridid
      character(len=*),intent(in) :: F_nomvar_S
      integer,intent(in),optional :: F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3,F_lni,F_lnj
      character(len=*),intent(in),optional :: F_etiket_S,F_typvar_S
      logical,intent(in),optional :: F_rewrite_L,F_writegrid_L
      !@author
      !@return
      integer :: F_istat
      !@/
      logical :: rewrite_L,writegrid_L
      integer :: istat,dateo,deet,npas,npak,dtype,ip3,gridid,nijk(3),i0,j0
      real,pointer :: data(:,:,:)
      character(len=RMN_ETK_LEN) :: etiket_S
      character(len=RMN_VARTYPE_LEN) :: typvar_S
      ! ---------------------------------------------------------------------
      call msg_toall(MSG_DEBUG,'(fstmpi) fstmpi_write_3d_r4_vgd [BGN]')
      F_istat = priv_init(F_data,F_gridid,dateo,deet,npas,npak,dtype,ip3,nijk, &
           typvar_S,etiket_S,rewrite_L,writegrid_L)
      if (.not.RMN_IS_OK(F_istat)) return

      if (present(F_dateo)) dateo = F_dateo
      if (present(F_deet)) deet = F_deet
      if (present(F_npas)) npas = F_npas
      if (present(F_npak)) npak = F_npak
      if (present(F_dtype)) dtype = F_dtype
      if (present(F_ip3)) ip3 = F_ip3
      if (present(F_typvar_S)) typvar_S = F_typvar_S
      if (present(F_etiket_S)) etiket_S = F_etiket_S
      if (present(F_rewrite_L)) rewrite_L = F_rewrite_L
      if (present(F_writegrid_L)) writegrid_L = F_writegrid_L
      if (present(F_lni)) nijk(1) = min(F_lni,nijk(1))
      if (present(F_lnj)) nijk(2) = min(F_lnj,nijk(2))

      F_istat = priv_prep(data,gridid,F_data,F_gridid,nijk,i0,j0)
      if (.not.RMN_IS_OK(F_istat)) return

      if (ptopo_isblocmaster_L) then
         F_istat = RMN_ERR
         if (associated(data).and.RMN_IS_OK(gridid)) then
            F_istat = fst_write(F_fileid,F_nomvar_S,data,gridid,F_vgridid,dateo, &
                 deet,npas,npak,dtype,ip3,typvar_S,etiket_S,rewrite_L,writegrid_L)
         endif
      endif
      !call collect_error(F_istat)

      if (associated(data)) deallocate(data,stat=istat)
      call msg_toall(MSG_DEBUG,'(fstmpi) fstmpi_write_3d_r4_vgd [END]')
      ! ---------------------------------------------------------------------
      return
   end function fstmpi_write_3d_r4_vgd


   !==== Private Functions =================================================


   function priv_init(F_data,F_gridid,dateo,deet,npas,npak,dtype,ip3,nijk, &
        typvar_S,etiket_S,rewrite_L,writegrid_L) result(F_istat)
      implicit none
      real,pointer :: F_data(:,:,:)
      integer,intent(in) :: F_gridid
      integer,intent(out) :: dateo,deet,npas,npak,dtype,ip3,nijk(:)
      character(len=*),intent(out) :: etiket_S,typvar_S
      logical,intent(out) :: rewrite_L,writegrid_L
      integer :: F_istat,istat
      !----------------------------------------------------------------------
      call ptopo_init_var()
      F_istat = RMN_OK
      dateo = 0 ; deet = 0 ; npas =0
      npak = FST_NPAK_DEFAULT
      dtype = RMN_DTYPE_IEEE
      ip3 = ptopo_grid_ibloc
      typvar_S = 'P'
      etiket_S = ' '
      rewrite_L = RMN_APPEND
      writegrid_L = .true.
      if (.not.(associated(F_data).and.RMN_IS_OK(F_gridid))) F_istat = RMN_ERR
      call collect_error(F_istat)
      if (.not.RMN_IS_OK(F_istat)) return
      nijk = ubound(F_data)
      istat = fst_write_opt(FST_OPT_OUTGRID_TYPE,FST_DIEZEGRID)
      if (ptopo_grid_nbloc == 1) &
           istat = fst_write_opt(FST_OPT_OUTGRID_TYPE,FST_ZGRID)
      !----------------------------------------------------------------------
      return
   end function priv_init


   function priv_prep(data,gridid,F_data,F_gridid,nijk,i0,j0) result(F_istat)
      implicit none
      real,pointer :: data(:,:,:),F_data(:,:,:)
      integer :: gridid,F_gridid,nijk(:),i0,j0,F_istat
      integer :: istat,nib,njb,hx,hy
      !----------------------------------------------------------------------
      call msg_toall(MSG_DEBUG,'(fstmpi) priv_prep [BGN]')
      F_istat = ptopo_collect_dims(RPN_COMM_BLOC_COMM,nijk(1),nijk(2),nib,njb,i0,j0)
      nullify(data)
      if (ptopo_isblocmaster_L) then
         allocate(data(nib,njb,nijk(3)),stat=istat)
         if (istat /= 0) then
            F_istat = RMN_ERR
            call msg(MSG_ERROR,'(fst_mpi_write) Problem allocating mem')
         endif
      endif
      F_istat = min(ptopo_collect(data,F_data,RPN_COMM_BLOC_COMM,i0,j0,nijk(1),nijk(2)),F_istat)
      call collect_error(F_istat)
      if (.not.RMN_IS_OK(F_istat)) return

      hx = max(0,1-lbound(F_data,1)) !TODO: need better way to define halo
      hy = max(0,1-lbound(F_data,2))
      gridid = ezgrid_merge(F_gridid,RPN_COMM_BLOC_COMM,.not.BCAST2ALL,1+hx,1+hy,nijk(1),nijk(2))
      call msg_toall(MSG_DEBUG,'(fstmpi) priv_prep [END]')
      !----------------------------------------------------------------------
      return
   end function priv_prep


   function priv_prep_s(F_dataout,F_grididout_S,F_datain,F_grididin_S,F_nijk, &
        F_ci0,F_cj0,F_cin,F_cjn,F_cdij) result(F_istat)
      implicit none
      integer, parameter :: N_ACHAR = 256
      integer, parameter :: N_SKIP  = 32
      integer, parameter :: N_ACHAR2 = N_ACHAR - N_SKIP
      real,pointer :: F_dataout(:,:,:),dataout(:,:,:),F_datain(:,:,:)
      integer :: F_nijk(:),F_ci0,F_cj0,F_cin,F_cjn,F_cdij,F_istat
      character(len=*) :: F_grididout_S,F_grididin_S
      integer :: grid_id,grid_gi0,grid_gj0,grid_lni,grid_lnj,grid_hx,grid_hy,&
           grid_dieze,grid_periodx,grid_periody,copyx,copyy
      integer :: istat,i0,j0,grid_id2,nib,njb,gi0,gj0,gnij(2),gi0c,gj0c,ginc, &
           gjnc,nibc,njbc,gi0bc,gj0bc,li0bc,lj0bc
      logical :: periodx_L, periody_L
      !----------------------------------------------------------------------
      call msg_toall(MSG_DEBUG,'(fstmpi) priv_prep_s [BGN]')
      nullify(F_dataout)
      F_istat = hgrid_wb_get(F_grididin_S,grid_id,grid_gi0,grid_gj0,grid_lni, &
           grid_lnj,grid_hx,grid_hy,grid_dieze,grid_periodx,grid_periody)
      call collect_error(F_istat)
      if (.not.RMN_IS_OK(F_istat)) return

      !TODO: remove halo from grid if grid_hx,grid_hy

      !# Output grid name
      F_grididout_S = trim(F_grididin_S)//'/bloc'
      if (any((/F_ci0,F_cj0,F_cin,F_cjn/) /= -1).or.F_cdij/=1) &
           F_grididout_S = trim(F_grididin_S)//'/b/'// &
           achar(N_SKIP+mod(abs(F_ci0),N_ACHAR2))// &
           achar(N_SKIP+mod(abs(F_cj0),N_ACHAR2))// &
           achar(N_SKIP+mod(abs(F_cin),N_ACHAR2))// &
           achar(N_SKIP+mod(abs(F_cjn),N_ACHAR2))// &
           achar(N_SKIP+mod(abs(F_cdij),N_ACHAR2))

      !# revert back to older interface if not dieze grid
      if (grid_dieze /= HGRID_DIEZE) then
         F_istat = priv_prep(F_dataout,grid_id2,F_datain,grid_id,F_nijk,i0,j0)
         !TODO: crop
         istat = hgrid_wb_put(F_grididout_S,grid_id2)
         return
      endif

      !# Bloc Collect
      F_istat = ptopo_collect_dims_ij0(RPN_COMM_BLOC_COMM,F_nijk(1),F_nijk(2), &
           grid_gi0,grid_gj0,nib,njb,i0,j0,gi0,gj0)
      istat = ezgrid_params(grid_id,gnij)

      nullify(F_dataout,dataout)
      if (ptopo_isblocmaster_L) then
         allocate(dataout(nib,njb,F_nijk(3)),stat=istat)
         if (.not.(associated(dataout) .and. istat == 0)) then
            F_istat = RMN_ERR
            call msg(MSG_ERROR,'(fst_mpi_write) Problem allocating mem')
         endif
         F_dataout => dataout
      endif
      F_istat = min(ptopo_collect(dataout,F_datain,RPN_COMM_BLOC_COMM,i0,j0,F_nijk(1),F_nijk(2)),F_istat)
      call collect_error(F_istat)
      if (.not.RMN_IS_OK(F_istat)) return

      !# Crop & periodic row/col
      F_ci0 = max(F_ci0,1) ; F_cj0 = max(F_cj0,1)
      if (F_cin < 1) F_cin = gnij(1)
      if (F_cjn < 1) F_cjn = gnij(2)
      F_cin = min(F_cin,gnij(1)) ; F_cjn = min(F_cjn,gnij(2))

      if (F_cdij /= 1) then
         call msg(MSG_WARNING,'(fstmpi_write) output aggregation not yet supported, writing none aggregated data')
         F_cdij = 1 !TODO: Allow agregation
      endif

      gi0c = max(F_ci0,gi0)           ; gj0c = max(F_cj0,gj0)
      ginc = min(gi0+nib-1,F_cin)     ; gjnc = min(gj0+njb-1,F_cjn)
      nibc = ginc - gi0c + 1          ; njbc = gjnc - gj0c + 1
      gi0bc = gi0c -F_ci0+1           ; gj0bc = gj0c - F_cj0+1
      li0bc = gi0c-gi0+1              ; lj0bc = gj0c-gj0+1
      !TODO: revise if F_cdij > 1

!!$      print *,'8',':p',ptopo_grid_ipe,ptopo_grid_ibloc,':g  ',1,gnij(1),':',gnij(1)
!!$      print *,'8',':p',ptopo_grid_ipe,ptopo_grid_ibloc,':gl ',grid_gi0,grid_gi0+grid_lni-1,':',grid_lni
!!$      print *,'8',':p',ptopo_grid_ipe,ptopo_grid_ibloc,':gb ',gi0,gi0+nib-1,':',nib
!!$
!!$      print *,'8',':p',ptopo_grid_ipe,ptopo_grid_ibloc,':guc',F_ci0,F_cin,':',F_cin-F_ci0+1
!!$      print *,'8',':p',ptopo_grid_ipe,ptopo_grid_ibloc,':gc ',gi0c,ginc,':',nibc
!!$      print *,'8',':p',ptopo_grid_ipe,ptopo_grid_ibloc,':gbc',gi0bc
!!$      print *,'8',':p',ptopo_grid_ipe,ptopo_grid_ibloc,':lbc',li0bc,ginc-gi0+1,':',(ginc-gi0+1)-li0bc+1

      periodx_L = (F_ci0==1 .and. F_cin==gnij(1) .and. grid_periodx==HGRID_PERIODIC)
      periody_L = (F_cj0==1 .and. F_cjn==gnij(2) .and. grid_periody==HGRID_PERIODIC)
      grid_periodx = 0 ; grid_periody = 0 ; 
      if (periodx_L .and. ginc==gnij(1)) grid_periodx = 1
      if (periody_L .and. gjnc==gnij(2)) grid_periody = 1
  
      IS_MASTER: if (ptopo_isblocmaster_L .and. (nib /= nibc .or. periodx_L .or. periody_L)) then
         nullify(F_dataout)
         if (nibc > 0 .and. njbc > 0) then
            allocate(F_dataout(nibc+grid_periodx,njbc+grid_periody,F_nijk(3)),stat=istat)
            !TODO: revise if F_cdij > 1
            F_dataout(1:nibc,1:njbc,:) = dataout(li0bc:li0bc+nibc-1,lj0bc:lj0bc+njbc-1,:)
            copyx = 0 ; copyy = 0
            if (periodx_L) copyx = 1
            if (periody_L) copyy = 1
            istat = ptopo_copyfirst2last(F_dataout,copyx,copyy,.true.)
         endif
         deallocate(dataout,stat=istat)
      endif IS_MASTER

      grid_id2 = ezgrid_sub(grid_id,F_ci0,F_cj0,F_cin,F_cjn,F_cdij,F_cdij)
      if (periodx_L .or. periody_L) &
           grid_id2 = ezgrid_addperiod(grid_id2,periodx_L,periody_L)
      istat = hgrid_wb_put(F_grididout_S,grid_id2,gi0bc,gj0bc,nibc+grid_periodx, &
           njbc+grid_periody,F_quiet_L=.true.)

      call msg_toall(MSG_DEBUG,'(fstmpi) priv_prep_s [END]')
      !----------------------------------------------------------------------
      return
   end function priv_prep_s


end module fstmpi_write_mod
