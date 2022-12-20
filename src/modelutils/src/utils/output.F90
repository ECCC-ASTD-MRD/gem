!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

!/@
module output_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use wb_itf_mod
   use outcfg_mod
   use output_files_mod
   use fstmpi_mod
   use ezgrid_mod
   use hgrid_wb
   use vgrid_wb
   use vGrid_Descriptors
   use ptr_store
   use rmn_gmm
   implicit none
   private
   !@objective 
   !@author Stephane Chamberland, 2012-01
   !@description
   ! Public functions
   public :: output_new, output_getlist, output_close, output_set_basedir, &
        output_get_basedir, output_set_diag_level, output_get_diag_level, &
        output_writevar, output_writestep, output_set_postproc
   ! Public constants
   !
!@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   interface output_new
      module procedure outcfg_new_4
      module procedure outcfg_new_8
   end interface

   interface output_getlist
      module procedure outcfg_getlist
   end interface

   interface output_writevar
      module procedure output_writevar_r4_2d
      module procedure output_writevar_r4_3d
   end interface

   interface output_close
      module procedure output_files_close
   end interface

   interface output_set_basedir
      module procedure output_files_set_basedir
   end interface

   interface output_get_basedir
      module procedure output_files_get_basedir
   end interface

   interface output_set_postproc
      module procedure output_files_set_postproc
   end interface

   integer, external :: output_writestep !TODO: include full interface

   integer,parameter :: NLEVELS_MAX = 32768
   real,parameter :: PRES_MB2PA = 100.

contains

   !/@*
   function output_set_diag_level(F_id,F_diag_level) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_diag_level
      !@return
      integer :: F_istat
      !*@/
      character(len=WB_MAXSTRINGLENGTH) :: name_S
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG,'(output) set_diag_level [BEGIN]')
      write(name_S,OUTPUT_FILES_WBFMT) F_id,'diaglvl'
      F_istat = wb_put(name_S,F_diag_level,WB_REWRITE_MANY)
      write(name_S,'(i6)') F_diag_level
      if (RMN_IS_OK(F_istat)) then
         call msg(MSG_INFO,'(output) set_diag_level: '//trim(name_S))
      else
         call msg(MSG_WARNING,'(output) Probleme saving diag_level: '//trim(name_S))
      endif
      call msg(MSG_DEBUG,'(output) set_diag_level [END]')
      !----------------------------------------------------------------------
      return
   end function output_set_diag_level


   !/@*
   function output_get_diag_level(F_id) result(F_diag_level)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id
      !@return
      integer :: F_diag_level
      !*@/
      integer :: istat
      character(len=WB_MAXSTRINGLENGTH) :: name_S
      !----------------------------------------------------------------------
      write(name_S,OUTPUT_FILES_WBFMT) F_id,'diaglvl'
      istat = wb_get(name_S,F_diag_level)
      if (.not.RMN_IS_OK(istat)) then
         F_diag_level = -999
      endif
      !----------------------------------------------------------------------
      return
   end function output_get_diag_level


   !/@*
   function output_writevar_r4_2d(F_id,F_step,F_varname_S,F_data,F_hgrid_S,F_vgrid_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_step
      character(len=*),intent(in) :: F_varname_S
      real,pointer :: F_data(:,:)
      character(len=*),intent(in) :: F_hgrid_S,F_vgrid_S
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2012-02
      !*@/
      character(len=128) :: msg_S
      integer :: l_ijk(3),u_ijk(3)
      real,pointer :: data3d(:,:,:)
      !----------------------------------------------------------------------
      !TODO-later: try to avoid "data3d = data2d" with a simple pointer association
      F_istat = RMN_OK
      write(msg_S,'(a,I5.5,a)') '(Output) Step=',F_step,' ('//trim(F_varname_S)//')'
      nullify(data3d)
      if (associated(F_data)) then
         l_ijk = 1
         u_ijk = 1
         l_ijk(1:2) = lbound(F_data)
         u_ijk(1:2) = ubound(F_data)
         call ptr_store_get(data3d,l_ijk,u_ijk)
         if (.not.associated(data3d)) then
            call msg(MSG_ERROR,trim(msg_S)//' unable to allocate workspace')
            F_istat = RMN_ERR
         endif
      endif
      call collect_error(F_istat)
      if (.not.RMN_IS_OK(F_istat)) then
         call ptr_store_free(data3d)
         return
      endif
      data3d(:,:,1) = F_data(:,:)
      F_istat = output_writevar_r4_3d(F_id,F_step,F_varname_S,data3d,F_hgrid_S,F_vgrid_S)
      call ptr_store_free(data3d)
      !----------------------------------------------------------------------
      return
   end function output_writevar_r4_2d


   !/@*
   function output_writevar_r4_3d(F_id,F_step,F_varname_S,F_data,F_hgrid_S,F_vgrid_S,F_sfc_field_ref) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_step
      character(len=*),intent(in) :: F_varname_S
      real,pointer :: F_data(:,:,:)
      character(len=*),intent(in) :: F_hgrid_S,F_vgrid_S
      real,pointer,optional :: F_sfc_field_ref(:,:) !- sfc field for vgd_levels
      !@return
      integer :: F_istat
      !@author Stephane Chamberland, 2012-01
      !*@/
      integer,parameter :: NLIST_MAX = 128

!!$      logical :: periodx_L,periody_L
!!$      integer :: gridid_out,lni,lnj,gi0,gj0,li0,lin,lj0,ljn,l_ijk2(3),u_ijk2(3)
      integer :: istat, istat2, datev, dateo, deet, l_ijk(3), u_ijk(3), &
           mystep, nidxlist, idxlist(NLIST_MAX), nn, idx, dtype, nbits, &
           filt_pass, ip3, nlevels, nlinbot, dij, fileid, level_type_in, &
           surf_level_idx, crop_gi0, crop_gj0, crop_gin, crop_gjn, &
           src_grid_id, src_grid_gi0, src_grid_gj0,src_lni, src_lnj
      logical :: rewrite_L,flip_L,fullplane_L
      integer,target :: ip1list_out(NLEVELS_MAX)
      character(len=16) :: etk_S,grid_etk_S,lvltype_S,v_int_S
      character(len=128) :: msg_S
      type(vgrid_descriptor) :: vgrid_in
      real :: levels(NLEVELS_MAX),filt_coef
      integer,pointer :: ip1list_in(:),ip1list_out_p(:)
      real,pointer :: data2(:,:,:),data3(:,:,:),sfc_field(:,:)
      !----------------------------------------------------------------------
      nullify(ip1list_in,ip1list_out_p,data2,data3,sfc_field)
      F_istat = RMN_OK
      write(msg_S,'(a,I5.5,a)') '(Output) Step=',F_step,' ('//trim(F_varname_S)//')'
      if (.not.associated(F_data)) then
         call msg(MSG_ERROR,trim(msg_S)//' no data')
         F_istat = RMN_ERR
      endif
      istat = vgrid_wb_get(F_vgrid_S,vgrid_in,ip1list_in,level_type_in)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_ERROR,trim(msg_S)//' invalid vgrid')
         F_istat = RMN_ERR
      endif
      istat = outcfg_time(F_id,F_step,datev,dateo,deet)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_ERROR,trim(msg_S)//' unable to get time info')
         F_istat = RMN_ERR
      endif

      if (associated(F_data)) then
         l_ijk= lbound(F_data)
         u_ijk = ubound(F_data)
         call ptr_store_get(data2,l_ijk,u_ijk)
         if (istat/=0) then
            call msg(MSG_ERROR,trim(msg_S)//' unable to allocate workspace')
            F_istat = RMN_ERR
         endif
      endif

      call collect_error(F_istat)
      if (.not.RMN_IS_OK(F_istat)) then
         call ptr_store_free(data2)
        return
      endif

      surf_level_idx = output_get_diag_level(F_id)
      mystep = output_files_get_closestep(F_id,F_step)
      nidxlist = outcfg_get_idxlist(F_id,F_step,F_varname_S,idxlist)
      DO_IDX: do nn=1,nidxlist
         istat = RMN_OK

         idx = idxlist(nn)
         istat = outcfg_var_meta(F_id,idx,F_varname_S,&
              dtype,nbits,filt_pass,filt_coef,rewrite_L,flip_L,etk_S,ip3,fullplane_L)
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_ERROR,trim(msg_S)//' unable to get var meta')
         endif

         !- Vertical
         istat2 = outcfg_level_meta(F_id,idx,lvltype_S,v_int_S,levels,nlevels,nlinbot)
         if (.not.RMN_IS_OK(istat2) .or. nlevels<=0) then
            call msg(MSG_ERROR,trim(msg_S)//' unable to get level meta')
            istat = RMN_ERR
         endif
         if (nlevels > u_ijk(3)) then
            call msg(MSG_ERROR,trim(msg_S)//' Not enough workspace levels')
            istat = RMN_ERR
         endif

         ip1list_out_p => ip1list_out
         IF_OK: if (RMN_IS_OK(istat)) then
            IF_P_OR_M: if (lvltype_S == 'p') then

               nullify(sfc_field)
               if (present(F_sfc_field_ref)) sfc_field => F_sfc_field_ref
               nlevels = priv_vinterp(sfc_field,vgrid_in,levels(1:nlevels), &
                    nlinbot,v_int_S,F_varname_S,F_data,ip1list_in,data2, &
                    ip1list_out_p)
               if (.not.RMN_IS_OK(nlevels)) then
                  call msg(MSG_ERROR, trim(msg_S)// &
                       ' Problem computing vertical interpolation ; lvl='// &
                       trim(F_vgrid_S)//' ==> P')
                  istat = RMN_ERR
               endif

            else !IF_P_OR_M

               nlevels = priv_level_subset(levels(1:nlevels),F_data,level_type_in,ip1list_in,surf_level_idx,data2,ip1list_out)
               if (.not.RMN_IS_OK(nlevels)) then
                  call msg(MSG_ERROR,trim(msg_S)//' Problem sampling levels ; lvl='//trim(F_vgrid_S))
                  istat = RMN_ERR
               endif

            endif IF_P_OR_M
         end if IF_OK

         !- Horizontal
         istat = outcfg_grid_meta(F_id,idx,crop_gi0,crop_gj0,crop_gin,crop_gjn,dij,grid_etk_S)
         istat = hgrid_wb_get(F_hgrid_S,src_grid_id,src_grid_gi0,src_grid_gj0,src_lni,src_lnj)
         data3 => data2

!!$         gridid_out = priv_outgrid_d(F_id,idx,F_hgrid_S,gi0,gj0,li0,lj0,lni,lnj,dij,periodx_L,periody_L)
!!$         lin = li0+lni-1 ; ljn = lj0+lnj-1
!!$
!!$         !#TODO: what tile is outside requested domain?
!!$         ! All tiles need to go through fstmpi_write (mpi hang otherwise)...
!!$         ! but some w/o data!!!
!!$
!!$         if (.not.RMN_IS_OK(gridid_out)) then
!!$            call msg(MSG_ERROR,trim(msg_S)//' Problem defining the output grid ; hgrid='//trim(F_hgrid_S))
!!$            istat = RMN_ERR
!!$         endif
!!$         !TODO-later: try pointer association instead of mem copy
!!$         l_ijk2 = 1
!!$         u_ijk2 = (/lni,lnj,nlevels/)
!!$         call ptr_store_get(data3,l_ijk2,u_ijk2)
!!$         if (.not.(associated(data2).and.associated(data3))) then
!!$            istat = RMN_ERR
!!$         else
!!$            !# data3 => data2(li0:lin,lj0:ljn,1:nlevels)
!!$            data3(1:u_ijk2(1),1:u_ijk2(2),1:u_ijk2(3)) = data2(li0:lin,lj0:ljn,1:nlevels)
!!$         endif

         if (RMN_IS_OK(istat) .and. filt_pass > 0) then
            !TODO-later: implement filter
            call msg(MSG_WARNING,trim(msg_S)//' Field filtering not yet supported for:'//trim(F_varname_S))
         endif
         if (RMN_IS_OK(istat) .and. dij > 1) then
            !TODO-later: implement data aggreg specified by reduc
            call msg(MSG_WARNING,trim(msg_S)//' Data aggregation not yet supported')
         endif
         if (RMN_IS_OK(istat) .and. flip_L) then
            !TODO-later:  implement data flip (x-z plan)
            call msg(MSG_WARNING,trim(msg_S)//' Field fliping not yet supported')
         endif
         if (RMN_IS_OK(istat) .and. fullplane_L) then
            !TODO-later:  implement data fullplane_L
            call msg(MSG_WARNING,trim(msg_S)//' fullplane_L not yet supported')
         endif

         fileid = output_files_open(F_id,mystep,lvltype_S)
         if (.not.RMN_IS_OK(fileid)) then
            call msg(MSG_ERROR,trim(msg_S)//' Cannot write field, file not open')
            istat = RMN_ERR
         endif

         call collect_error(istat)
         if (.not.RMN_IS_OK(istat)) then
            F_istat = RMN_ERR
            call ptr_store_free(data3)
            cycle 
            !- WARNING: make call collect_error before any cycle or return to avoid fstmpi_write MPI hang/deadlock
         endif

!!$         print *,trim(F_varname_S),minval(data3),maxval(data3),nlevels, &
!!$              gridid_out,ip1list_out_p(1:nlevels),nbits,ip3 ; call flush(6)

         !TODO: pass vgrid_S (updated with ip1list_out) to have vgrid in file
         !TODO: when model levels, need to write F_sfc_field_ref on same grid
!!$         F_istat = min(fstmpi_write(fileid,F_varname_S,data3,&
!!$              gridid_out,ip1list_out_p(1:nlevels),dateo,deet,F_step, &
!!$              -abs(nbits),dtype,ip3,'P',etk_S,rewrite_L,F_writegrid_L=.true., &
!!$              F_periodx_L=periodx_L,F_periodx_L=periody_L),F_istat)
         F_istat = min(&
              & fstmpi_write( &
              &   fileid,F_varname_S,data3, &
              &    F_hgrid_S,ip1list_out_p(1:nlevels),dateo,deet,F_step, &
              &    -abs(nbits),dtype,ip3,'P',etk_S,rewrite_L,.true., &
              &    src_lni,src_lnj,crop_gi0,crop_gj0,crop_gin,crop_gjn,dij), &
              & F_istat)
         if (RMN_IS_OK(F_istat)) then
            call msg(MSG_INFO,trim(msg_S)//' Written on lvltype='//lvltype_S(1:1))
         else
            call msg(MSG_ERROR,trim(msg_S)//' Writing Probleme')
         endif
         call ptr_store_free(data3)
      enddo DO_IDX
      call ptr_store_free(data2)
      !----------------------------------------------------------------------
      return
   end function output_writevar_r4_3d


   !==== Private Functions =================================================

!!$   !/@*
!!$   function priv_outgrid_d(F_id,F_idx,F_hgrid_S,F_gi0,F_gj0,F_li0,F_lj0,F_lni, &
!!$        F_lnj,F_dij,F_periodx_L,F_periody_L) result(F_gridid_out)
!!$      implicit none
!!$      !@objective Return a #-grid corped to user's request with i0,j0,ln,lni 
!!$      !           of local-mpi tile, li0,lj0 are relative to the croped grid
!!$      !@arguments
!!$      integer,intent(in) :: F_id,F_idx
!!$      character(len=*),intent(in) :: F_hgrid_S
!!$      integer,intent(out) :: F_gi0,F_gj0 !#position 1,1 of the local data crop in the dst dieze grid
!!$      integer,intent(out) :: F_li0,F_lj0,F_lni,F_lnj,F_dij !#local data crop limit
!!$      logical,intent(out) :: F_periodx_L,F_periody_L
!!$      !@return
!!$      integer :: F_gridid_out
!!$      !@description
!!$      !- Note: since ezsinct grid axes cannot have negative indexes (halo)
!!$      !  grid indexes are  1:ni+2*halo,   1:nj+2*halo
!!$      !  while for data    1-halo:ni+halo,1-halo:nj+halo
!!$      !  thus src_i0/j0 are relative to grid indexes (=data_index+halo)
!!$      !       dst_i0/j0 are relatives to data indexes
!!$      !
!!$      !  grid:   1     1+hg                        gni+hg    gni+2hg 
!!$      !          |......|-------|=========|-----------|........|
!!$      !  data:  1-hg    1       gi0    gi0+lni-1     gni     gni+hg
!!$      !
!!$      !  where:
!!$      !    hg : global halo or BCS region
!!$      !    gi0: index i=1 uof local domain in global domain
!!$      !    gni: global domain size (w/o halos)
!!$      !    ...: global halo or  BCS region
!!$      !    ---: computational domain outside of the local mpi tile
!!$      !    ===: computational domain in the local mpi tile
!!$      !
!!$      !         1-hg    1       gi0    gi0+lni-1         gni     gni+hg
!!$      !  src:    |......|---|----|=========|------|-------|........|
!!$      !                    di0                   din
!!$      !  dst:               |----|=========|------|
!!$      !                     1   dgi0  dgi0+lni-1 dni
!!$      !  local:                  |=========|
!!$      !                         li0       lin
!!$      !
!!$      !  where:
!!$      ! 
!!$      !*@/
!!$      character(len=16) :: grid_etk_S
!!$      integer :: istat,src_grid_id,src_grid_gnij(2),src_data_gni,src_data_gnj
!!$      integer :: src_grid_gi0, src_grid_gj0, src_lni, src_lnj, src_hx, src_hy
!!$      integer :: src_crop_gi0, src_crop_gj0, src_crop_gin, src_crop_gjn
!!$      integer :: src_crop_gni, src_crop_gnj
!!$      integer :: src_cropg_gi0,src_cropg_gj0,src_cropg_gin,src_cropg_gjn
!!$      logical :: src_dieze_L,src_periodx_L,src_periody_L
!!$      !----------------------------------------------------------------------
!!$      call msg_toall(MSG_DEBUG,'(Output) priv_outgrid_d [BGN]')
!!$      F_gridid_out = RMN_ERR
!!$      F_gi0 = 1 ; F_gj0 = 1
!!$      F_li0 = 1 ; F_lj0 = 1
!!$      F_lni = 0 ; F_lnj = 0 
!!$      F_dij = 1
!!$      F_periodx_L = .false. ; F_periody_L = .false.
!!$
!!$      !TODO-later: Make sure provided grid is on full domain with i0,j0,lni,lnj defining the local mpi tile (# style)
!!$
!!$      !# Get src grid info
!!$      istat = hgrid_wb_get(F_hgrid_S,src_grid_id,src_grid_gi0,src_grid_gj0, &
!!$           src_lni,src_lnj,src_hx,src_hy,src_dieze,src_periodx,src_periody)
!!$    if (.not.RMN_IS_OK(istat)) then
!!$         call msg(MSG_ERROR,'(Output) invalid hgrid')
!!$         return
!!$      endif
!!$      F_lni = src_lni ; F_lnj = src_lnj
!!$
!!$      istat = ezgrid_params(src_grid_id,F_nij=src_grid_gnij)
!!$      src_data_gni = src_grid_gnij(1)-2*src_hx
!!$      src_data_gnj = src_grid_gnij(2)-2*src_hy
!!$
!!$      !# Crop to user resquested sub-set (removing halos)
!!$      istat = outcfg_grid_meta(F_id,F_idx,src_crop_gi0,src_crop_gj0,src_crop_gin,src_crop_gjn,F_dij,grid_etk_S)
!!$      if (src_crop_gin < 0 .or. src_crop_gjn < 0) then
!!$         if (src_crop_gin < 0) then
!!$            F_periodx_L = src_periodx_L
!!$            src_crop_gin = src_data_gni
!!$         endif
!!$         if (src_crop_gjn < 0) then
!!$            F_periody_L = src_periody_L
!!$            src_crop_gjn = src_data_gnj
!!$         endif
!!$      endif
!!$
!!$      !# Make sure the desired grid is inside src grid limits (w/o halos)
!!$      src_crop_gi0 = max(1,src_crop_gi0)
!!$      src_crop_gj0 = max(1,src_crop_gj0)
!!$      src_crop_gin = min(src_crop_gin,src_data_gni)
!!$      src_crop_gjn = min(src_crop_gjn,src_data_gnj)
!!$      src_crop_gni = src_crop_gin - src_crop_gi0 + 1
!!$      src_crop_gnj = src_crop_gjn - src_crop_gj0 + 1
!!$
!!$      !# Shift indexes (for halos) from data to grid
!!$      src_cropg_gi0 = src_crop_gi0 + src_hx
!!$      src_cropg_gj0 = src_crop_gj0 + src_hy
!!$      src_cropg_gin = src_crop_gin + src_hx
!!$      src_cropg_gjn = src_crop_gjn + src_hy
!!$
!!$      !# Create output grid
!!$      F_gridid_out = ezgrid_sub(src_grid_id,src_cropg_gi0,src_cropg_gj0,src_cropg_gin,src_cropg_gjn)
!!$
!!$      !# Compute Local ij0,ijn on output subgrid
!!$      F_gi0 = src_grid_gi0 - (src_crop_gi0-1)
!!$      F_gj0 = src_grid_gj0 - (src_crop_gi0-1)
!!$      if (F_gi0 < 1) F_li0 = abs(F_gi0-1)+1
!!$      if (F_gj0 < 1) F_lj0 = abs(F_gj0-1)+1
!!$      F_lni = min((src_lni-F_li0+1),(src_crop_gni-F_gi0+1))
!!$      F_lnj = min((src_lnj-F_lj0+1),(src_crop_gnj-F_gj0+1))
!!$
!!$      call msg_toall(MSG_DEBUG,'(Output) priv_outgrid_d [END]')
!!$     !----------------------------------------------------------------------
!!$      return
!!$   end function priv_outgrid_d


!TODO: priv_outgrid_z needs to be reviewed
!!$   !/@*
!!$   function priv_outgrid_z(F_id,F_idx,F_hgrid_S,F_li0,F_lj0,F_lin,F_ljn, &
!!$        F_dij,F_periodx_L,F_periody_L) result(F_gridid_out)
!!$      implicit none
!!$      !@objective Return a Z-grid corped to user's request and to local 
!!$      !           local-mpi tile, li0,lj0 are relative to the full source grid
!!$      !@arguments
!!$      integer,intent(in) :: F_id,F_idx
!!$      character(len=*),intent(in) :: F_hgrid_S
!!$      integer,intent(out) :: F_li0,F_lj0,F_lin,F_ljn,F_dij !#local data crop limit
!!$      logical,intent(out) :: F_periodx_L,F_periody_L
!!$      !@return
!!$      integer :: F_gridid_out
!!$      !@description
!!$      !- Note: since ezsinct grid axes cannot have negative indexes (halo)
!!$      !  grid indexes are  1:ni+2*halo,   1:nj+2*halo
!!$      !  while for data    1-halo:ni+halo,1-halo:nj+halo
!!$      !  thus src_i0/j0 are relative to grid indexes (=data_index+halo)
!!$      !       dst_i0/j0 are relatives to data indexes
!!$      !*@/
!!$      character(len=16) :: grid_etk_S
!!$      integer :: istat,src_grid_id,src_grid_gnij(2),src_data_gnij(2)
!!$      integer :: src_grid_gi0,src_grid_gj0,src_lni,src_lnj,src_hx,src_hy
!!$      integer :: src_data_gi0,src_data_gj0,src_data_gin,src_data_gjn
!!$      integer :: dst_data_gi0,dst_data_gj0,dst_data_gin,dst_data_gjn
!!$      integer :: dst_grid_gi0,dst_grid_gj0,dst_grid_gin,dst_grid_gjn
!!$      logical :: src_dieze_L,src_periodx_L,src_periody_L
!!$      !----------------------------------------------------------------------
!!$      call msg_toall(MSG_DEBUG,'(Output) priv_outgrid_z [BGN]')
!!$      F_gridid_out = RMN_ERR
!!$      F_periodx_L = .false. ; F_periody_L = .false.
!!$
!!$      !TODO-later: Make sure provided grid is on full domain with i0,j0,lni,lnj defining the local mpi tile (# style)
!!$
!!$      istat = hgrid_wb_get(F_hgrid_S,src_grid_id,src_grid_gi0,src_grid_gj0, &
!!$           src_lni,src_lnj,src_hx,src_hy,src_dieze_L,src_periodx_L,src_periody_L)
!!$    if (.not.RMN_IS_OK(istat)) then
!!$         call msg(MSG_ERROR,'(Output) invalid hgrid')
!!$         return
!!$      endif
!!$
!!$      !# Shift indexes (for halos) from grid to data
!!$   print '(a,7i6)','(Outgrid) '//trim(F_hgrid_S),src_grid_id,src_grid_gi0, &
!!$        src_grid_gj0,src_lni,src_lnj,src_hx,src_hy ; call flush(6)
!!$      src_data_gi0 = src_grid_gi0 - src_hx
!!$      src_data_gj0 = src_grid_gj0 - src_hy
!!$      src_data_gin = src_data_gi0 + src_lni - 1 !TODO: double check if hx != 0
!!$      src_data_gjn = src_data_gj0 + src_lnj - 1 !TODO: double check if hy != 0
!!$
!!$      !# Crop to user resquested sub-set
!!$      istat = outcfg_grid_meta(F_id,F_idx,dst_data_gi0,dst_data_gj0,dst_data_gin,dst_data_gjn,F_dij,grid_etk_S)
!!$      if (dst_data_gin < 0 .or. dst_data_gjn < 0) then
!!$         istat = ezgrid_params(src_grid_id,F_nij=src_grid_gnij)
!!$         src_data_gnij(1) = src_grid_gnij(1) - 2*src_hx
!!$         src_data_gnij(2) = src_grid_gnij(2) - 2*src_hy
!!$         if (dst_data_gin < 0) then
!!$            if (src_periodx_L) F_periodx_L = .true.
!!$            dst_data_gin = src_data_gnij(1) + dst_data_gin + 1
!!$         endif
!!$         if (dst_data_gjn < 0) then
!!$            if (src_periody_L) F_periody_L = .true.
!!$            dst_data_gjn = src_data_gnij(2) + dst_data_gjn + 1
!!$         endif
!!$      endif
!!$
!!$      dst_data_gi0 = max(src_data_gi0,dst_data_gi0)
!!$      dst_data_gj0 = max(src_data_gj0,dst_data_gj0)
!!$      dst_data_gin = min(dst_data_gin,src_data_gin)
!!$      dst_data_gjn = min(dst_data_gjn,src_data_gjn)
!!$
!!$ !!$      if (src_periodx_L .and. &
!!$ !!$           dst_data_gi0 == src_data_gi0 .and. &
!!$ !!$           dst_data_gin == dst_data_gin)  F_addx = 1
!!$ !!$      if (src_periody_L .and. &
!!$ !!$           dst_data_gj0 == src_data_gj0 .and. &
!!$ !!$           dst_data_gjn == dst_data_gjn)  F_addy = 1
!!$
!!$      !# Shift back indexes (for halos) from data to grid
!!$      dst_grid_gi0 = dst_data_gi0 + src_hx
!!$      dst_grid_gj0 = dst_data_gj0 + src_hy
!!$      dst_grid_gin = dst_data_gin + src_hx
!!$      dst_grid_gjn = dst_data_gjn + src_hy
!!$
!!$      !# Create output grid
!!$      F_gridid_out = ezgrid_sub(src_grid_id,dst_grid_gi0,dst_grid_gj0,dst_grid_gin,dst_grid_gjn)
!!$ !!$      if (RMN_IS_OK(F_gridid_out) .and. (F_addx>0 .or.F_addy>0)) then
!!$ !!$         src_grid_gi0 = F_gridid_out
!!$ !!$         istat = ezgrid_addperiod(F_gridid_out,src_grid_gi0, &
!!$ !!$              src_grid_gj0,src_lni,src_lnj,src_hx,src_hy, &
!!$ !!$              src_dieze_L,(F_addx>0),(F_addy>0))
!!$ !!$        endif
!!$
!!$      F_li0 = dst_data_gi0 - src_data_gi0 + 1
!!$      F_lin = dst_data_gin - src_data_gi0 + 1
!!$      F_lj0 = dst_data_gj0 - src_data_gj0 + 1
!!$      F_ljn = dst_data_gjn - src_data_gj0 + 1
!!$
!!$      call msg_toall(MSG_DEBUG,'(Output) priv_outgrid_z [END]')
!!$     !----------------------------------------------------------------------
!!$      return
!!$   end function priv_outgrid_z


   !/@*
   function priv_level_subset(F_levels,F_data_in,F_level_type_in,F_ip1list_in, &
        F_surf_level_idx,F_data_out,F_ip1list_out) result(F_nip1_out)
      implicit none
      !@arguments
      real,intent(in) :: F_levels(:)
      integer :: F_level_type_in,F_surf_level_idx,F_ip1list_in(:),F_ip1list_out(:)
      real,pointer :: F_data_in(:,:,:),F_data_out(:,:,:)
      !@return
      integer :: F_nip1_out
      !*@/
      integer :: nlevels,nlevels2,ii,ii2
      logical :: all_levels_L
      !----------------------------------------------------------------------
      F_nip1_out = RMN_ERR

      nlevels = size(F_levels)

      all_levels_L = (nint(F_levels(1)) == -1)
      if (.not.all_levels_L .and. F_levels(1) < 1.e-5) then
         if (any(F_level_type_in == (/VGRID_SURF_TYPE,VGRID_GROUND_TYPE/))) then
            all_levels_L = .true.
         else
            if (F_surf_level_idx >= lbound(F_ip1list_in,1) .and. &
                 F_surf_level_idx <= ubound(F_ip1list_in,1)) then
!!$            if (F_surf_level_idx >= lbound(F_ip1list_in)) then
               ii2 = min(F_surf_level_idx,ubound(F_data_in,3))
               F_ip1list_out(1) = F_ip1list_in(ii2)
               F_data_out(:,:,1) = F_data_in(:,:,ii2)
               F_nip1_out = 1
            else
               call msg(MSG_ERROR,'(Output) surf/diag level for up air fields requested, surface/diag level not defined')
            endif
            return
         endif
      endif

      IF_LEVEL: if (all_levels_L) then
         nlevels2 =size(F_data_in,3)
         if (ubound(F_ip1list_in,1) < nlevels2) then
            call msg(MSG_ERROR,'(Output) Inconsistent vertical description')
!!$            print *,'F_ip1list_in=',F_ip1list_in
!!$            print *,'size(F_ip1list_in),nlevels2=',ubound(F_ip1list_in,1),nlevels2
            return
         endif
         if (size(F_data_out,3) < nlevels2) then
            call msg(MSG_ERROR,'(Output) Not enough workspace levels')
            return
         endif
         F_ip1list_out(1:nlevels2) = F_ip1list_in(1:nlevels2)
         F_data_out(:,:,1:nlevels2) = F_data_in(:,:,1:nlevels2)
      else if (all(F_levels(1:nlevels) > 0.)) then
         nlevels2 = 0
         do ii = 1,nlevels
            ii2 = nint(F_levels(ii))
            if (all(ii2 <= (/size(F_ip1list_in),size(F_data_in,3)/))) then
               nlevels2 = nlevels2 + 1
               F_ip1list_out(nlevels2) = F_ip1list_in(ii2)
               F_data_out(:,:,nlevels2) = F_data_in(:,:,ii2)
            endif
         enddo
      else
         call msg(MSG_ERROR,'(Output) Level list not valid (<0)')
         return
      endif IF_LEVEL

      F_nip1_out = nlevels2
      !----------------------------------------------------------------------
      return
   end function priv_level_subset


   !/@*
   function priv_vinterp(F_sfc_field_ref,F_vgrid_in,F_levels,F_nlinbot, &
        F_v_int_S,F_varname_S,F_data_in,F_ip1list_in,F_data_out,F_ip1list_out) &
        result(F_nip1_out)
      implicit none
      !@arguments
      type(vgrid_descriptor),intent(in) :: F_vgrid_in
      real,intent(in) :: F_levels(:)
      integer,intent(in) :: F_nlinbot
      character(len=*),intent(in) :: F_v_int_S,F_varname_S
      integer,pointer :: F_ip1list_in(:),F_ip1list_out(:)
      real,pointer :: F_sfc_field_ref(:,:),F_data_in(:,:,:),F_data_out(:,:,:)
      !@return
      integer :: F_nip1_out
      !*@/
      integer :: istat,nlevels,kk,nlinbot,ikind
      real :: levels(NLEVELS_MAX)
      real,pointer :: levels_3d_in(:,:,:)
      type(gmm_metadata) :: mymeta
      character(len=16) :: sfc_field_S,tmp_S
      !----------------------------------------------------------------------
      F_nip1_out = RMN_ERR
 
      nlevels = size(F_levels)
      if (any(F_levels(1:nlevels) <= 0.)) then
         call msg(MSG_ERROR,'(Output) Pres level list not valid (<0)')
         return
      endif

      if (.not.associated(F_sfc_field_ref)) then
         istat = vgd_get(F_vgrid_in,key='RFLD',value=sfc_field_S)
         if (istat /= VGD_OK) sfc_field_S = ''
         istat = gmm_get(sfc_field_S,F_sfc_field_ref,mymeta)
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_ERROR,'(Output) Sfc ref field not found ('//trim(sfc_field_S)//') Cannot vertically interpolate')
            return
         endif
      endif
      !TODO-later: make sure F_sfc_field_ref is in PA
      istat = vgd_levels(F_vgrid_in,F_ip1list_in,levels_3d_in,F_sfc_field_ref,in_log=.false.)
      if (.not.(RMN_IS_OK(istat) .and. associated(levels_3d_in))) then
         call msg(MSG_ERROR,'(Output) Problem getting input 3d levels for vertical interpolation')
         return
      endif

      ikind = RMN_CONV_PRESS
      tmp_S = ''
      levels(1:nlevels) = F_levels(1:nlevels)
      do kk=1,nlevels
         call convip_plus(F_ip1list_out(kk),levels(kk),ikind,RMN_CONV_P2IPNEW,tmp_S,.not.RMN_CONV_USEFORMAT_L)
      enddo

      nlinbot = F_nlinbot
      if (any(F_v_int_S(1:1) == (/'l','L'/))) nlinbot = nlevels
      levels(1:nlevels) = levels(1:nlevels)*PRES_MB2PA
      call vte_intvertx_isodst(F_data_out,F_data_in,levels_3d_in,levels, &
           size(F_data_in,1) * size(F_data_in,2),size(F_data_in,3),nlevels, &
           F_varname_S,nlinbot)

      if (associated(levels_3d_in)) deallocate(levels_3d_in,stat=istat)

      F_nip1_out = nlevels
      !----------------------------------------------------------------------
      return
   end function priv_vinterp

end module output_mod
