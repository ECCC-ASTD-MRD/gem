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
#include <rmn/msg.h>

!/@
module fst_write_mod
   use vGrid_Descriptors
   use ezgrid_mod
   use hgrid_wb
   use vgrid_wb
   implicit none
   private
   !@objective 
   !@author  Stephane Chamberland, 2011-04
   !@description
   ! Public functions
   public :: fst_write,fst_write_opt
   ! Public constants
   integer,parameter,public :: FST_NPAK_DEFAULT = -16
   integer,parameter,public :: FST_NPAK_FULL32 = -32
   integer,parameter,public :: FST_OPT_OUTGRID_TYPE = 1
   integer,parameter,public :: FST_ZGRID = 1
   integer,parameter,public :: FST_DIEZEGRID = 2
!@/

#include <rmnlib_basics.hf>
!!!#include <arch_specific.hf>

   interface fst_write
      module procedure fst_write_2d_r4
      module procedure fst_write_2d_r4_s
      module procedure fst_write_3d_r4_ip1
      module procedure fst_write_3d_r4_ip1_s
      module procedure fst_write_3d_r4_vgd
      module procedure fst_write_3d_r4_vgd_s
   end interface


   integer,save :: m_outgrid_type = FST_DIEZEGRID

contains

   !/@
   function fst_write_opt(F_opt,F_value) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      integer,intent(in) :: F_opt
      integer,intent(in) :: F_value
      !@author
      !@return
      integer :: F_istat
      !@/
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR
      select case(F_opt)
      case(FST_OPT_OUTGRID_TYPE)
         if (F_value >= 0) m_outgrid_type = F_value
         F_istat = m_outgrid_type
      end select
      ! ---------------------------------------------------------------------
      return
   end function fst_write_opt


   !/@
   function fst_write_2d_r4_s(F_fileid,F_nomvar_S,F_data,F_gridid_S,F_ip1, &
        F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3,F_typvar_S,F_etiket_S, &
        F_rewrite_L,F_writegrid_L) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      integer,intent(in) :: F_fileid
      real,pointer :: F_data(:,:)
      character(len=*),intent(in) :: F_gridid_S
      character(len=*),intent(in) :: F_nomvar_S
      integer,intent(in),optional :: F_dateo,F_deet,F_npas,F_ip1,F_npak,F_dtype,F_ip3
      character(len=*),intent(in),optional :: F_typvar_S,F_etiket_S
      logical,intent(in),optional :: F_rewrite_L,F_writegrid_L
      !@author
      !@return
      integer :: F_istat
      !@/
      integer :: npak,dtype,dateo,deet,npas,ip1,ip2,ip3,nijk(3),ig14(4),lnij(2),ij0(2)
      logical :: writegrid_L,rewrite_L
      character(len=2) :: grtyp_S
      character(len=RMN_ETK_LEN) :: etiket_S
      character(len=RMN_VARTYPE_LEN) :: typvar_S
      real,pointer :: data3d(:,:,:)
      real,target :: dummy3d(1,1,1)
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(fst) fst_write_2d_r4 [BGN]')
      nullify(data3d)
      if (associated(F_data)) data3d => dummy3d
      F_istat = priv_init(data3d,RMN_OK,dateo,deet,npas,npak,dtype,ip3,nijk, &
           typvar_S,etiket_S,rewrite_L,writegrid_L)
      ip1 = 0
      if (present(F_npak)) npak = F_npak
      if (present(F_dtype)) dtype = F_dtype
      if (present(F_dateo)) dateo = F_dateo
      if (present(F_deet)) deet = F_deet
      if (present(F_npas)) npas = F_npas
      if (present(F_ip1)) ip1 = F_ip1
      if (present(F_ip3)) ip3 = F_ip3
      if (present(F_etiket_S)) etiket_S = F_etiket_S
      if (present(F_typvar_S)) typvar_S = F_typvar_S
      if (present(F_rewrite_L)) rewrite_L = F_rewrite_L
      if (present(F_writegrid_L)) writegrid_L = F_writegrid_L

      F_istat = priv_grid_s(F_fileid,F_gridid_S,F_nomvar_S,F_etiket_S, &
           writegrid_L,grtyp_S,ig14,lnij,ij0)
      if (.not.RMN_IS_OK(F_istat)) return

      nijk(1:2) = ubound(F_data)
      if (any(lnij /= nijk(1:2))) then
         call msg(MSG_ERROR,'(fst_write) Inconsistent grid/data dims')
         F_istat = RMN_ERR
         return
      endif

      ip2 = nint(dble(deet) * dble(npas) / 3600.D0)
      F_istat = fstecr(F_data,F_data,npak,F_fileid,dateo,deet,npas, &
           nijk(1),nijk(2),1,ip1,ip2,ip3,typvar_S,F_nomvar_S(1:4),etiket_S, &
           grtyp_S(1:1),ig14(1),ig14(2),ig14(3),ig14(4),dtype,rewrite_L)

      call msg(MSG_DEBUG,'(fst) fst_write_2d_r4 [END]')
      ! ---------------------------------------------------------------------
      return
   end function fst_write_2d_r4_s


   !/@
   function fst_write_2d_r4(F_fileid,F_nomvar_S,F_data,F_gridid,F_ip1,F_dateo, &
        F_deet,F_npas,F_npak,F_dtype,F_ip3,F_typvar_S,F_etiket_S,F_rewrite_L, &
        F_writegrid_L) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      integer,intent(in) :: F_fileid
      real,pointer :: F_data(:,:)
      integer,intent(in) :: F_gridid
      character(len=*),intent(in) :: F_nomvar_S
      integer,intent(in),optional :: F_dateo,F_deet,F_npas,F_ip1,F_npak,F_dtype,F_ip3
      character(len=*),intent(in),optional :: F_typvar_S,F_etiket_S
      logical,intent(in),optional :: F_rewrite_L,F_writegrid_L
      !@author
      !@return
      integer :: F_istat
      !@/
      integer :: npak,dtype,dateo,deet,npas,ip1,ip2,ip3,nijk(3),ig14(4),lnij(2),ij0(2),lnij0(2)
      logical :: writegrid_L,rewrite_L
      character(len=2) :: grtyp_S
      character(len=RMN_ETK_LEN) :: etiket_S
      character(len=RMN_VARTYPE_LEN) :: typvar_S
      real,pointer :: data3d(:,:,:)
      real,target :: dummy3d(1,1,1)
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(fst) fst_write_2d_r4 [BGN]')
      nullify(data3d)
      if (associated(F_data)) data3d => dummy3d
      F_istat = priv_init(data3d,F_gridid,dateo,deet,npas,npak,dtype,ip3,nijk, &
           typvar_S,etiket_S,rewrite_L,writegrid_L)
      ip1 = 0
      if (present(F_npak)) npak = F_npak
      if (present(F_dtype)) dtype = F_dtype
      if (present(F_dateo)) dateo = F_dateo
      if (present(F_deet)) deet = F_deet
      if (present(F_npas)) npas = F_npas
      if (present(F_ip1)) ip1 = F_ip1
      if (present(F_ip3)) ip3 = F_ip3
      if (present(F_etiket_S)) etiket_S = F_etiket_S
      if (present(F_typvar_S)) typvar_S = F_typvar_S
      if (present(F_rewrite_L)) rewrite_L = F_rewrite_L
      if (present(F_writegrid_L)) writegrid_L = F_writegrid_L

      ij0 = (/1,1/)
      lnij0 = (/huge(1),huge(1)/)
      F_istat = priv_grid(F_fileid,F_gridid,ij0,lnij0,F_nomvar_S,F_etiket_S, &
           writegrid_L,grtyp_S,ig14,lnij)
      if (.not.RMN_IS_OK(F_istat)) return

      nijk(1:2) = ubound(F_data)
      if (any(lnij /= nijk(1:2))) then
         call msg(MSG_ERROR,'(fst_write) Inconsistent grid/data dims')
         F_istat = RMN_ERR
         return
      endif

      ip2 = nint(dble(deet) * dble(npas) / 3600.D0)
      F_istat = fstecr(F_data,F_data,npak,F_fileid,dateo,deet,npas, &
           nijk(1),nijk(2),1,ip1,ip2,ip3,typvar_S,F_nomvar_S(1:4),etiket_S, &
           grtyp_S(1:1),ig14(1),ig14(2),ig14(3),ig14(4),dtype,rewrite_L)

      call msg(MSG_DEBUG,'(fst) fst_write_2d_r4 [END]')
      ! ---------------------------------------------------------------------
      return
   end function fst_write_2d_r4


   !/@
   function fst_write_3d_r4_vgd_s(F_fileid,F_nomvar_S,F_data,F_gridid_S, &
        F_vgridid,F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3,F_typvar_S, &
        F_etiket_S,F_rewrite_L,F_writegrid_L) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      integer,intent(in) :: F_fileid
      real,pointer :: F_data(:,:,:)
      integer,intent(in) :: F_vgridid
      character(len=*),intent(in) :: F_nomvar_S,F_gridid_S
      integer,intent(in),optional :: F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3
      character(len=*),intent(in),optional :: F_typvar_S,F_etiket_S
      logical,intent(in),optional :: F_rewrite_L,F_writegrid_L
      !@author
      !@return
      integer :: F_istat
      !@/
      integer,pointer :: ip1list(:)
      integer :: npak,dtype,dateo,deet,npas,ip3,istat,nijk(3)
      logical :: writegrid_L,rewrite_L
      character(len=RMN_ETK_LEN) :: etiket_S
      character(len=RMN_VARTYPE_LEN) :: typvar_S
      type(vgrid_descriptor),pointer :: myvgrid  !#TODO: pointer?
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(fst) fst_write_3d_r4_vgd [BGN]')
      F_istat = priv_init(F_data,RMN_OK,dateo,deet,npas,npak,dtype,ip3,nijk, &
           typvar_S,etiket_S,rewrite_L,writegrid_L)
      if (present(F_npak)) npak = F_npak
      if (present(F_dtype)) dtype = F_dtype
      if (present(F_dateo)) dateo = F_dateo
      if (present(F_deet)) deet = F_deet
      if (present(F_npas)) npas = F_npas
      if (present(F_ip3)) ip3 = F_ip3
      if (present(F_etiket_S)) etiket_S = F_etiket_S
      if (present(F_typvar_S)) typvar_S = F_typvar_S
      if (present(F_rewrite_L)) rewrite_L = F_rewrite_L
      if (present(F_writegrid_L)) writegrid_L = F_writegrid_L

      nullify(myvgrid,ip1list)
      istat = vgrid_wb_get(F_vgridid,myvgrid,ip1list)
      if (.not.(RMN_IS_OK(istat).and.associated(ip1list))) then
          call msg(MSG_ERROR,'(fst_write) Problem getting vgrid description')
          return
      endif

      F_istat = fst_write_3d_r4_ip1_s(F_fileid,F_nomvar_S,F_data,F_gridid_S, &
           ip1list,dateo,deet,npas,npak,dtype,ip3,typvar_S,etiket_S,rewrite_L,writegrid_L)

      if (associated(ip1list)) deallocate(ip1list,stat=istat)

      if (writegrid_L) then
         istat = vgd_write(myvgrid,F_fileid,'fst')
      endif
      istat = vgd_free(myvgrid)

      call msg(MSG_DEBUG,'(fst) fst_write_3d_r4_vgd [END]')
      ! ---------------------------------------------------------------------
      return
   end function fst_write_3d_r4_vgd_s


   !/@
   function fst_write_3d_r4_vgd(F_fileid,F_nomvar_S,F_data,F_gridid,F_vgridid, &
        F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3,F_typvar_S,F_etiket_S, &
        F_rewrite_L,F_writegrid_L) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      integer,intent(in) :: F_fileid
      real,pointer :: F_data(:,:,:)
      integer,intent(in) :: F_gridid,F_vgridid
      character(len=*),intent(in) :: F_nomvar_S
      integer,intent(in),optional :: F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3
      character(len=*),intent(in),optional :: F_typvar_S,F_etiket_S
      logical,intent(in),optional :: F_rewrite_L,F_writegrid_L
      !@author
      !@return
      integer :: F_istat
      !@/
      integer,pointer :: ip1list(:)
      integer :: npak,dtype,dateo,deet,npas,ip3,istat,nijk(3)
      logical :: writegrid_L,rewrite_L
      character(len=RMN_ETK_LEN) :: etiket_S
      character(len=RMN_VARTYPE_LEN) :: typvar_S
      type(vgrid_descriptor),pointer :: myvgrid  !#TODO: pointer?
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(fst) fst_write_3d_r4_vgd [BGN]')
      F_istat = priv_init(F_data,F_gridid,dateo,deet,npas,npak,dtype,ip3,nijk, &
           typvar_S,etiket_S,rewrite_L,writegrid_L)
      if (present(F_npak)) npak = F_npak
      if (present(F_dtype)) dtype = F_dtype
      if (present(F_dateo)) dateo = F_dateo
      if (present(F_deet)) deet = F_deet
      if (present(F_npas)) npas = F_npas
      if (present(F_ip3)) ip3 = F_ip3
      if (present(F_etiket_S)) etiket_S = F_etiket_S
      if (present(F_typvar_S)) typvar_S = F_typvar_S
      if (present(F_rewrite_L)) rewrite_L = F_rewrite_L
      if (present(F_writegrid_L)) writegrid_L = F_writegrid_L

      nullify(myvgrid,ip1list)
      istat = vgrid_wb_get(F_vgridid,myvgrid,ip1list)
      if (.not.(RMN_IS_OK(istat).and.associated(ip1list))) then
          call msg(MSG_ERROR,'(fst_write) Problem getting vgrid description')
          return
      endif

      F_istat = fst_write_3d_r4_ip1(F_fileid,F_nomvar_S,F_data,F_gridid,ip1list, &
           dateo,deet,npas,npak,dtype,ip3,typvar_S,etiket_S,rewrite_L,writegrid_L)

      if (associated(ip1list)) deallocate(ip1list,stat=istat)

      if (writegrid_L) then
         istat = vgd_write(myvgrid,F_fileid,'fst')
      endif
      istat = vgd_free(myvgrid)

      call msg(MSG_DEBUG,'(fst) fst_write_3d_r4_vgd [END]')
      ! ---------------------------------------------------------------------
      return
   end function fst_write_3d_r4_vgd


   !/@
   function fst_write_3d_r4_ip1_s(F_fileid,F_nomvar_S,F_data,F_gridid_S, &
        F_ip1list,F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3,F_typvar_S, &
        F_etiket_S,F_rewrite_L,F_writegrid_L) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      integer,intent(in) :: F_fileid
      real,pointer :: F_data(:,:,:)
      !TODO-later: F_ip1list should be a pointer so lbound is preseved
      integer,intent(in) :: F_ip1list(:)
      character(len=*),intent(in) :: F_nomvar_S,F_gridid_S
      integer,intent(in),optional :: F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3
      character(len=*),intent(in),optional :: F_typvar_S,F_etiket_S
      logical,intent(in),optional :: F_rewrite_L,F_writegrid_L
      !@author
      !@return
      integer :: F_istat
      !@/
      integer :: npak,dtype,dateo,deet,npas,ip1,ip2,ip3,l_ijk(3),u_ijk(3), &
           nijk(3),ig14(4),k,arbitrarykind,lnij(2),ij0(2)
      real :: zp1
      logical :: writegrid_L,rewrite_L
      character(len=2) :: grtyp_S
      character(len=RMN_ETK_LEN) :: etiket_S
      character(len=RMN_VARTYPE_LEN) :: typvar_S
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(fst) fst_write_3d_r4 ip1 [BGN]')

      F_istat = priv_init(F_data,RMN_OK,dateo,deet,npas,npak,dtype,ip3,u_ijk, &
           typvar_S,etiket_S,rewrite_L,writegrid_L)
      if (present(F_npak)) npak = F_npak
      if (present(F_dtype)) dtype = F_dtype
      if (present(F_dateo)) dateo = F_dateo
      if (present(F_deet)) deet = F_deet
      if (present(F_npas)) npas = F_npas
      if (present(F_ip3)) ip3 = F_ip3
      if (present(F_etiket_S)) etiket_S = F_etiket_S
      if (present(F_typvar_S)) typvar_S = F_typvar_S
      if (present(F_rewrite_L)) rewrite_L = F_rewrite_L
      if (present(F_writegrid_L)) writegrid_L = F_writegrid_L

      F_istat = priv_grid_s(F_fileid,F_gridid_S,F_nomvar_S,F_etiket_S,writegrid_L, &
           grtyp_S,ig14,lnij,ij0)
      if (.not.RMN_IS_OK(F_istat)) return

      l_ijk = lbound(F_data)
      if (l_ijk(3) < 1) then
         call msg(MSG_WARNING,'(fst_write) k index <= 0 for ip1list not yet supported, will use default convip values')
      endif
      nijk  = shape(F_data)

      if (any(lnij /= nijk(1:2))) then
         call msg(MSG_ERROR,'(fst_write) Inconsistent grid/data dims')
         F_istat = RMN_ERR
         return
      endif
      if (size(F_ip1list) < u_ijk(3) .or. any(F_ip1list(1:u_ijk(3)) < 0))  then
         call msg(MSG_ERROR,'(fst_write) Missing ip1 values')
         F_istat = RMN_ERR
         return
      endif

      ip2 = nint(dble(deet) * dble(npas) / 3600.D0)
      arbitrarykind = RMN_CONV_ARBITRARY
      F_istat = RMN_OK
      do k=l_ijk(3),u_ijk(3)
         if (k > 0) then
            ip1 = F_ip1list(k)
         else
            zp1 = abs(float(k))
            call convip_plus(ip1, zp1, arbitrarykind, RMN_CONV_P2IPNEW, ' ', .not.RMN_CONV_USEFORMAT_L)
         endif
         F_istat = min(fstecr(F_data(:,:,k),F_data(:,:,k),npak,F_fileid,&
              dateo,deet,npas, &
              nijk(1),nijk(2),1,ip1,ip2,ip3,typvar_S,F_nomvar_S(1:4),etiket_S, &
              grtyp_S(1:1),ig14(1),ig14(2),ig14(3),ig14(4),dtype,rewrite_L),F_istat)
      enddo
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_ERROR,'(fst_write) problem writing: '//trim(F_nomvar_S))
      else
         call msg(MSG_DEBUG,'(fst_write) wrote: '//trim(F_nomvar_S))
      endif
      call msg(MSG_DEBUG,'(fst) fst_write_3d_r4 ip1 [END]')
      ! ---------------------------------------------------------------------
      return
   end function fst_write_3d_r4_ip1_s


   !/@
   function fst_write_3d_r4_ip1(F_fileid,F_nomvar_S,F_data,F_gridid,F_ip1list, &
        F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3,F_typvar_S,F_etiket_S, &
        F_rewrite_L,F_writegrid_L) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      integer,intent(in) :: F_fileid
      real,pointer :: F_data(:,:,:)
      !TODO-later: F_ip1list should be a pointer so lbound is preseved
      integer,intent(in) :: F_gridid,F_ip1list(:)
      character(len=*),intent(in) :: F_nomvar_S
      integer,intent(in),optional :: F_dateo,F_deet,F_npas,F_npak,F_dtype,F_ip3
      character(len=*),intent(in),optional :: F_typvar_S,F_etiket_S
      logical,intent(in),optional :: F_rewrite_L,F_writegrid_L
      !@author
      !@return
      integer :: F_istat
      !@/
      integer :: npak,dtype,dateo,deet,npas,ip1,ip2,ip3,l_ijk(3),u_ijk(3),nijk(3), &
           ig14(4),k,arbitrarykind,lnij(2),ij0(2),lnij0(2)
      real :: zp1
      logical :: writegrid_L,rewrite_L
      character(len=2) :: grtyp_S
      character(len=RMN_ETK_LEN) :: etiket_S
      character(len=RMN_VARTYPE_LEN) :: typvar_S
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(fst) fst_write_3d_r4 ip1 [BGN]')
      F_istat = priv_init(F_data,F_gridid,dateo,deet,npas,npak,dtype,ip3,u_ijk, &
           typvar_S,etiket_S,rewrite_L,writegrid_L)
      if (present(F_npak)) npak = F_npak
      if (present(F_dtype)) dtype = F_dtype
      if (present(F_dateo)) dateo = F_dateo
      if (present(F_deet)) deet = F_deet
      if (present(F_npas)) npas = F_npas
      if (present(F_ip3)) ip3 = F_ip3
      if (present(F_etiket_S)) etiket_S = F_etiket_S
      if (present(F_typvar_S)) typvar_S = F_typvar_S
      if (present(F_rewrite_L)) rewrite_L = F_rewrite_L
      if (present(F_writegrid_L)) writegrid_L = F_writegrid_L

      ij0 = (/1,1/)
      lnij0 = (/huge(1),huge(1)/)
      F_istat = priv_grid(F_fileid,F_gridid,ij0,lnij0,F_nomvar_S,F_etiket_S, &
           writegrid_L,grtyp_S,ig14,lnij)
      if (.not.RMN_IS_OK(F_istat)) return

      l_ijk = lbound(F_data)
      if (l_ijk(3) < 1) then
         call msg(MSG_WARNING,'(fst_write) k index <= 0 for ip1list not yet supported, will use default convip values')
      endif
      nijk  = shape(F_data)

      if (any(lnij /= nijk(1:2))) then
         call msg(MSG_ERROR,'(fst_write) Inconsistent grid/data dims')
         F_istat = RMN_ERR
         return
      endif
      if (size(F_ip1list) < u_ijk(3) .or. any(F_ip1list(1:u_ijk(3)) < 0))  then
         call msg(MSG_ERROR,'(fst_write) Missing ip1 values')
         F_istat = RMN_ERR
         return
      endif

      ip2 = nint(dble(deet) * dble(npas) / 3600.D0)
      arbitrarykind = RMN_CONV_ARBITRARY
      F_istat = RMN_OK
      do k=l_ijk(3),u_ijk(3)
         if (k > 0) then
            ip1 = F_ip1list(k)
         else
            zp1 = abs(float(k))
            call convip_plus(ip1, zp1, arbitrarykind, RMN_CONV_P2IPNEW, ' ', .not.RMN_CONV_USEFORMAT_L)
         endif
         F_istat = min(fstecr(F_data(:,:,k),F_data(:,:,k),npak,F_fileid,&
              dateo,deet,npas, &
              nijk(1),nijk(2),1,ip1,ip2,ip3,typvar_S,F_nomvar_S(1:4),etiket_S, &
              grtyp_S(1:1),ig14(1),ig14(2),ig14(3),ig14(4),dtype,rewrite_L),F_istat)
      enddo
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_ERROR,'(fst_write) problem writing: '//trim(F_nomvar_S))
      else
         call msg(MSG_DEBUG,'(fst_write) wrote: '//trim(F_nomvar_S))
      endif
      call msg(MSG_DEBUG,'(fst) fst_write_3d_r4 ip1 [END]')
      ! ---------------------------------------------------------------------
      return
   end function fst_write_3d_r4_ip1


   !==== Private Functions =================================================


   function priv_init(F_data,F_gridid,dateo,deet,npas,npak,dtype,ip3,nijk,typvar_S, &
        etiket_S,rewrite_L,writegrid_L) result(F_istat)
      implicit none
      real,pointer :: F_data(:,:,:)
      integer,intent(in) :: F_gridid
      integer,intent(out) :: dateo,deet,npas,npak,dtype,ip3,nijk(:)
      character(len=*),intent(out) :: etiket_S,typvar_S
      logical,intent(out) :: rewrite_L,writegrid_L
      integer :: F_istat
      !----------------------------------------------------------------------
      F_istat = RMN_OK
      dateo = 0 ; deet = 0 ; npas =0
      npak = FST_NPAK_DEFAULT
      dtype = RMN_DTYPE_IEEE
      ip3 = 0
      typvar_S = 'P'
      etiket_S = ' '
      rewrite_L = RMN_APPEND
      writegrid_L = .true.
      if (.not.(associated(F_data).and.RMN_IS_OK(F_gridid))) F_istat = RMN_ERR
      if (.not.RMN_IS_OK(F_istat)) return
      nijk = ubound(F_data)
      !----------------------------------------------------------------------
      return
   end function priv_init


   function priv_grid_S(F_fileid,F_gridid_S,F_nomvar_S,F_etiket_S,F_writegrid_L, &
        F_grtyp_S,F_ig14,F_lnij,F_ij0) result(F_istat)
      implicit none
      integer,intent(in) :: F_fileid
      logical,intent(in) :: F_writegrid_L
      character(len=*),intent(in) :: F_gridid_S,F_nomvar_S,F_etiket_S
      character(len=*),intent(out) :: F_grtyp_S
      integer,intent(out) :: F_ig14(:),F_lnij(:),F_ij0(:)
      integer :: F_istat
      integer :: gridid,hx,hy,nij(2),isdieze,isperiodx,isperiody
      !----------------------------------------------------------------------
      F_istat = hgrid_wb_get(F_gridid_S,gridid,F_ij0(1),F_ij0(2), &
           F_lnij(1),F_lnij(2),hx,hy,isdieze,isperiodx,isperiody)
      F_grtyp_S = 'Z'
      if (RMN_IS_OK(F_istat)) &
           F_istat = priv_grid(F_fileid,gridid,F_ij0,F_lnij,F_nomvar_S,F_etiket_S, &
           &                   F_writegrid_L,F_grtyp_S,F_ig14,nij)
     if (any(F_grtyp_S(1:1) == EZGRID_REF_TYPES)) then
        if (m_outgrid_type == FST_DIEZEGRID) then
           if (any(F_grtyp_S(1:1) == (/'z','Z'/))) F_grtyp_S = '#'
           F_ig14(3:4) = F_ij0
        else
           if (F_grtyp_S(1:1) == '#') F_grtyp_S = 'Z'
        endif
     endif
      !----------------------------------------------------------------------
      return
   end function priv_grid_S


   function priv_grid(F_fileid,F_gridid,F_ij0,F_lnij,F_nomvar_S,F_etiket_S, &
        F_writegrid_L,F_grtyp_S,F_ig14,F_nij) result(F_istat)
      implicit none
      integer,intent(in) :: F_fileid,F_gridid,F_lnij(:)
      integer,intent(inout) :: F_ij0(:)
      logical,intent(in) :: F_writegrid_L
      character(len=*),intent(in) :: F_nomvar_S,F_etiket_S
      character(len=*),intent(out) :: F_grtyp_S
      integer,intent(out) :: F_ig14(:),F_nij(2)
      integer :: F_istat
      integer :: istat,nij(2),ij0(2),igp14(4),gridid
      character(len=2) :: grref_S
      real,pointer :: ax(:,:),ay(:,:)
      !----------------------------------------------------------------------
      nullify(ax,ay)
      if (F_writegrid_L) then
         F_istat = ezgrid_params(F_gridid,F_nij,F_grtyp_S,grref_S,F_ig14,ij0,ax,ay,igp14)
        if (m_outgrid_type /= FST_DIEZEGRID) then
           if (any(F_ij0 > 1) .or. any(F_lnij < F_nij)) then
              gridid = ezgrid_sub(F_gridid,F_ij0(1),F_ij0(2), &
                   min(F_ij0(1)+F_lnij(1)-1,F_nij(1)), &
                   min(F_ij0(2)+F_lnij(2)-1,F_nij(2)))
              if (associated(ax)) deallocate(ax,stat=istat)
              if (associated(ay)) deallocate(ay,stat=istat)
              nullify(ax,ay)
              F_istat = ezgrid_params(gridid,F_nij,F_grtyp_S,grref_S,F_ig14,ij0,ax,ay,igp14)
              F_ij0 = (/1,1/)
           endif
        endif
      else
         F_istat = ezgrid_params(F_gridid,F_nij,F_grtyp_S,grref_S,F_ig14,ij0,F_igp14=igp14)
      endif
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_ERROR,'(fst_write) Problem getting grid params')
         deallocate(ax,ay,stat=istat)
         return
      endif

      if (associated(ax) .and. associated(ay)) then
         nij = shape(ax)
         F_istat = fstecr(ax,ax,FST_NPAK_FULL32,F_fileid,0,0,0, &
              nij(1),nij(2),1,igp14(1),igp14(2),igp14(3),'X', &
              '>>',F_etiket_S,grref_S(1:1),F_ig14(1),F_ig14(2),F_ig14(3),F_ig14(4), &
              RMN_DTYPE_IEEE,RMN_REWRITE)
         nij = shape(ay)
         F_istat = min(fstecr(ay,ay,FST_NPAK_FULL32,F_fileid,0,0,0, &
              nij(1),nij(2),1,igp14(1),igp14(2),igp14(3),'X', &
              '^^',F_etiket_S,grref_S(1:1),F_ig14(1),F_ig14(2),F_ig14(3),F_ig14(4), &
              RMN_DTYPE_IEEE,RMN_REWRITE),F_istat)
         if (.not.RMN_IS_OK(F_istat)) then
            call msg(MSG_ERROR,'(fst_write) problem writing grid for: '//trim(F_nomvar_S))
         else
            call msg(MSG_DEBUG,'(fst_write) wrote grid for: '//trim(F_nomvar_S))
         endif
      endif
      deallocate(ax,ay,stat=istat)
      if (any(F_grtyp_S(1:1) == EZGRID_REF_TYPES)) F_ig14 = igp14
      !----------------------------------------------------------------------
      return
   end function priv_grid

end module fst_write_mod
