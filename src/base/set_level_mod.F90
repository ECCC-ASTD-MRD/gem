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

!**s/r set_level_mod - Utility functions for output set_level routine
!
module set_level_mod

  implicit none
  
  private
  
  public :: set_level_usr_vgrid, set_level_usr_val

  integer, public :: set_level_OK = 0, set_level_ERROR = -1

  contains

    !===========================================================================
    integer function set_level_usr_vgrid(F_vu, F_file_S) result(istat)
      use vGrid_Descriptors, only: vgrid_descriptor, vgd_new, vgd_get, &
           vgd_print, VGD_OK, VGD_ERROR
      use levels, only: vgrid_user
      use ptopo, only: Ptopo_myproc
      use iso_c_binding
      use, intrinsic :: iso_fortran_env
      use lun, only: Lun_debug_L

      implicit none

      type(vgrid_user), intent(inout) :: F_vu
      character(len=*), intent(in) :: F_file_S
      ! Local variables
      integer, dimension(3) :: n123
      integer :: ier
      real(kind=REAL64), dimension(:,:,:), pointer :: vtbl_8
      include 'mpif.h'
      include "rpn_comm.inc"

      nullify(vtbl_8)

      istat = set_level_ERROR
      if( Ptopo_myproc == 0 )then
         if( set_level_read_vgrid(F_vu%vgd, F_file_S) == set_level_ERROR) return
         if( vgd_get(F_vu%vgd,"VTBL - vgrid_descriptor vtbl",vtbl_8) == &
              VGD_ERROR )then
            ! TODO handle error (only PE 0!!)
            print*,'ERROR in set_level_user with vgd_get for key VTBL'
         endif
         n123 = ubound(vtbl_8)         
      endif
      call rpn_comm_bcast(n123, size(n123), RPN_COMM_INTEGER, 0, "grid", ier)
      ! TODO handle error rpn_comm_bcast?
      if( Ptopo_myproc /= 0 )then
         allocate(vtbl_8(n123(1), n123(2), n123(3)), stat=ier)
         if( ier /= 0 )then
            print*,'ERROR in set_level_user, cannot allocate vtbl_8 of size', &
                 n123
            return
         endif
      endif
      call rpn_comm_bcast(vtbl_8, size(vtbl_8), RPN_COMM_REAL8, 0, "grid", ier)
      ! TODO handle error rpn_comm_bcast?
      if( Ptopo_myproc /= 0 )then
         if( vgd_new(F_vu%vgd,vtbl_8) == VGD_ERROR )then
            ! TODO handle error (all but PE 0!!)
            print*,'ERROR in set_level_user with vgd_new with table'
         endif
      endif
      !if( vgd_print(F_vu%vgd) == VGD_ERROR )then
      !   print*,'Warning vgd_print returned error code in ', &
      !        'set_level_usr_vgrid on file ',trim(F_file_S)
      !endif
      deallocate(vtbl_8)
      if( vgd_get(F_vu%vgd, "VCOD", value = F_vu%vcode ) == VGD_ERROR ) return
      F_vu%usr_grid_L=.true.
      if (Lun_debug_L) ier = vgd_print(F_vu%vgd)
      istat = set_level_OK

    end function set_level_usr_vgrid
    
    !===========================================================================

    Integer function set_level_read_vgrid(F_vgd, F_file_S) result(istat)
      use vGrid_Descriptors, only: vgrid_descriptor, vgd_new, vgd_get, VGD_OK, &
           VGD_ERROR, vgd_print
      use path
      implicit none
      type(vgrid_descriptor), intent(inout) :: F_vgd
      character(len=*), intent(in) :: F_file_S
      ! Local variables
      ! TODO get intelligent way to set unit number.
      integer, save :: unit = 62
      integer :: fnom, fstouv, fstfrm
      unit = unit + 1
      istat = set_level_ERROR

      if( fnom(unit,trim(Path_ind_S)//'/'//trim(F_file_S),'RND+OLD+R/O',0) &
           < 0 )then
         print*,'Error in set_level_read_vgrid with fnom on file ',&
              trim(Path_ind_S)//'/'//trim(F_file_S)
         return
      endif
      if( fstouv(unit,'RND') < 0 )then
         print*,'Error in set_level_read_vgrid with fstouv on file ',&
              trim(Path_ind_S)//'/'//trim(F_file_S)
         return
      endif
      if( vgd_new(F_vgd, unit) == VGD_ERROR )then
         print*,'ERROR with vgd_new in set_level_read_vgrid on file ', &
              trim(Path_ind_S)//'/'//trim(F_file_S)
         return
      endif
      !if( vgd_print(F_vgd) == VGD_ERROR )then
      !   print*,'Warning vgd_print returned error code in ', &
      !   'set_level_read_vgrid on file ', trim(Path_ind_S)//'/'//trim(F_file_S)
      !endif
      if( fstfrm(unit) < 0 )return
      call fclos(unit)
      istat = set_level_OK
    end function set_level_read_vgrid

    !===========================================================================

    integer function set_level_usr_eta_lnp(F_vu, F_lnp, F_nk, F_eta,&
         F_kind, F_in_lnp_L,F_level,F_level_max,F_minx,F_maxx,F_miny,F_maxy, F_p0) result(status)
      use vGrid_Descriptors, only: vgrid_descriptor, vgd_print, vgd_new,&
           vgd_get, vgd_levels, VGD_NO_REF_NOMVAR, VGD_OK, VGD_ERROR
      use levels, only: vgrid_user
      implicit none

      type(vgrid_user), intent(inout) :: F_vu
      real, dimension(:,:,:), pointer, intent(inout) :: F_lnp
      integer, intent(out) :: F_nk
      real, dimension(:), pointer, intent(inout) :: F_eta
      integer, intent(out) :: F_kind
      logical, intent(out) :: F_in_lnp_L
      integer :: F_level_max
      real, dimension(F_level_max) :: F_level
      integer, intent(in) :: F_minx,F_maxx,F_miny,F_maxy
      real, dimension(F_minx:F_maxx,F_miny:F_maxy), optional, intent(in) :: F_p0
      ! Local varibales
      integer :: ier, version, kind, k, ind, isize
      integer, dimension(:), pointer :: ip1s,ip1sub
      real, dimension(:), pointer :: lnp_1D
      character(len=1) :: dummy_L
      character(len=4) :: rfld_S
      logical :: alloc_ip1sub_L
      logical, save :: done_L=.false.
      real, save, dimension(:), pointer :: lnp_target,eta_target

      if(.not. done_L)then
         done_L=.true.
         nullify(lnp_target,eta_target)
      endif
      nullify(ip1s,ip1sub,lnp_1D)

      status = set_level_ERROR
      if( vgd_get(F_vu%vgd, "KIND - vertical coordinate ip1 kind",&
           value = F_kind )  == VGD_ERROR ) return
      if( vgd_get(F_vu%vgd, "VERS - vertical coordinate version",&
           value = version ) == VGD_ERROR ) return
      if( vgd_get(F_vu%vgd, "RFLD", value = rfld_S, quiet = .true. ) == VGD_ERROR )then
         ! vgd_get will be in error if there is no RFLD, like for pressure
         ! We can disredard this errror. There will be no error message du to
         ! the quiet = .true. option
      endif
      if(F_vu%vcode == 1002)then
         F_in_lnp_L = .true.
      elseif(F_vu%vcode == 2001)then
         F_in_lnp_L = .true.
      elseif(F_vu%vcode == 4001)then
         F_in_lnp_L = .false.
      else
         print*,'Expecting one of the following vertical descriptor:'
         print*,'   kind = 1 and version = 2 (eta)'
         print*,'   kind = 2 and version = 1 (pressure)'
         print*,'   kind = 4 and version = 1 (Height above surface)'
         print*,'   but got kind = ', F_kind,' version = ',version
         return
      endif
      if( vgd_get(F_vu%vgd, "NL_M - number of momentum levels", value = F_nk )&
           == VGD_ERROR ) return    
      
      allocate( ip1s(F_nk), stat = ier)      
      if(ier /= 0 )then
         print*,'Allocation problem in set_level_usr_eta_lnp on ip1s of size =>'&
              , F_nk
         return
      endif
      if( vgd_get(F_vu%vgd, "VIPM - level ip1 list (m)", value = ip1s )&
           == VGD_ERROR ) return

      if(F_nk == F_level_max)then
         alloc_ip1sub_L=.false.
         ip1sub => ip1s
      else
         alloc_ip1sub_L=.true.
         allocate(ip1sub(F_level_max))
         do k=1,F_level_max
            ind=nint(F_level(k))
            if(ind < 1 .or. ind > F_nk)then
               print*,'Problem with level list in set_level_usr_eta_lnp, got value of',ind
               print*,'but expecting values between 1 and ',F_nk
               deallocate(ip1sub)
               return
            endif
            ip1sub(k)=ip1s(ind)
         end do 
         F_nk=F_level_max
      endif

      isize=(F_maxx-F_minx+1)*(F_maxy-F_miny+1)*F_nk
      if(associated(lnp_target))then
         if(size(lnp_target) < isize)deallocate(lnp_target)
      endif
      if(.not.associated(lnp_target))then         
         allocate(lnp_target(isize), stat = ier)
         if(ier /= 0 )then
            print*,'Allocation problem in set_level_usr_eta_lnp on lnp_target of size =>'
            print*,'(F_maxx-F_minx+1)*(F_maxy-F_miny+1)*F_nk=',isize
            return
         endif
      endif
      F_lnp(F_minx:F_maxx,F_miny:F_maxy,1:F_nk) => lnp_target(1:isize)
      
      if( rfld_S == VGD_NO_REF_NOMVAR)then
         if( vgd_levels(F_vu%vgd, ip1sub, lnp_1D, in_log=F_in_lnp_L) == VGD_ERROR )&
              return       
         do k=1, F_nk
            F_lnp(F_minx:F_maxx, F_miny:F_maxy, k) = lnp_1D(k)
         enddo
      else
         if(.not. present(F_p0))then
            print*,'ERROR in set_level_usr_eta_lnp, argument F_p0 must be present'
            return
         endif
         if( vgd_levels(F_vu%vgd, ip1sub, F_lnp, F_p0, in_log=F_in_lnp_L) == VGD_ERROR )&
              return
      endif

      if(associated(eta_target))then
         if(size(eta_target) < F_nk)deallocate(eta_target)
      endif
      if(.not.associated(eta_target))then
         allocate(eta_target(F_nk), stat = ier)
         if(ier /= 0 )then
            print*,'Allocation problem in set_level_usr_eta_lnp on eta_target of size =>'
            print*,'F_nk=',F_nk
            return
         endif
      endif
      F_eta(1:F_nk) => eta_target(1:F_nk)
      do k=1,F_nk
         call convip( ip1sub(k), F_eta(k), kind, -1, dummy_L, .false.)
      end do
      deallocate(ip1s)
      if(associated(lnp_1D))deallocate(lnp_1D)
      if(alloc_ip1sub_L)deallocate(ip1sub)
      status = set_level_OK
      return
    end function set_level_usr_eta_lnp
    
    !===================================================================

  integer function set_level_usr_val(F_val, F_indo, F_rf, F_kind, F_nko, F_usr_src, &
       F_levset,F_level,F_level_max, F_lev_src_name_S, &
       F_minx, F_maxx, F_miny, F_maxy, F_nk, F_stag_S,&
       F_level_typ_S, F_hgrid_index, F_is_lnp_L) result(istat)
    use glb_ld
    use gmm_pw, only: gmmk_pw_p0_plus_s, pw_p0_plus
    use gmm_geof, only : gmmk_fis0_s, fis0
    use rmn_gmm, only: gmm_get
    use levels, only: Level_typ_S, Level_vgrid_usr,MAXLEV,MAXSET
    use outgrid, only: OutGrid_hgrid_usr
    use vGrid_Descriptors, only: vgrid_descriptor, vgd_get, vgd_free, VGD_OK, &
         VGD_ERROR
    use tdpack, only : grav_8

    implicit none

    real, dimension(:,:,:), pointer, intent(inout)  :: F_val, F_usr_src
    integer, dimension(:), pointer, intent(inout) :: F_indo
    real, dimension(:), pointer, intent(inout) :: F_rf
    integer, intent(out) :: F_kind
    integer, intent(out) :: F_nko
    integer, intent(in) :: F_levset,F_hgrid_index
    real, dimension(MAXLEV,MAXSET) :: F_level
    integer, dimension(MAXSET) :: F_level_max
    integer, intent(in) :: F_minx, F_maxx, F_miny, F_maxy, F_nk
    character(len=*), intent(in) :: F_lev_src_name_S,F_level_typ_S
    character(len=*), intent(out) :: F_stag_S
    logical, intent(out), optional :: F_is_lnp_L
    
    ! Local variables
    integer :: ier, k
    real p0(F_minx:F_maxx, F_miny:F_maxy)
    logical :: is_lnp_L
    logical, save :: done_L=.false.
    integer, save, dimension(:), pointer :: indo_target

    if(.not. done_L)then
       done_L=.true.
       nullify(indo_target)
    endif
    
    istat = set_level_ERROR

    if( trim(F_lev_src_name_S) /= 'MOMENTUM' .and. trim(F_lev_src_name_S) /= 'THERMO' )then
       print*,'ERROR in set_level_usr_val, argument F_lev_src_name_S must be set to '
       print*,'"MOMENTUM" or "THERMO", got ',trim(F_lev_src_name_S)
       return
    endif

    if( set_level_get_usr_src(F_usr_src, F_levset, F_lev_src_name_S, &
         F_minx, F_maxx, F_miny, F_maxy, F_nk) == set_level_ERROR )then
       print*,'ERROR in set_level_usr_val, with set_level_get_usr_src'
       return
    endif

    if( Level_vgrid_usr(F_levset)%vcode == 1002 )then
       ! eta
       F_stag_S(2:4)=F_level_typ_S(1:1)//" "//OutGrid_hgrid_usr(F_hgrid_index)%usr_grid_index_S
       ier = gmm_get(gmmk_pw_p0_plus_s, pw_p0_plus)
       p0 = pw_p0_plus
       call out_padbuf(p0,F_minx,F_maxx,F_miny,F_maxy,1)
       if( set_level_usr_eta_lnp(Level_vgrid_usr(F_levset), F_val, F_nko, F_rf, &
            F_kind, is_lnp_L,F_level(1,F_levset),F_level_max(F_levset),&
            F_minx, F_maxx, F_miny, F_maxy, F_p0=p0) == SET_LEVEL_ERROR )then
          ! TODO handle error gracefully
          print*,'Problem in set_level_usr_val with levels of Vcode' ,Level_vgrid_usr(F_levset)%vcode, ' ', &
               trim(Level_typ_S(F_levset))
          stop
       endif
       if(present(F_is_lnp_L))F_is_lnp_L=is_lnp_L
    elseif(Level_vgrid_usr(F_levset)%vcode == 2001)then
       ! pressure
       F_stag_S(2:4)="P "//OutGrid_hgrid_usr(F_hgrid_index)%usr_grid_index_S
       if( set_level_usr_eta_lnp(Level_vgrid_usr(F_levset), F_val, F_nko, F_rf, &
            F_kind, is_lnp_L,F_level(1,F_levset),F_level_max(F_levset),&
            F_minx, F_maxx, F_miny, F_maxy) == SET_LEVEL_ERROR )then
          ! TODO handle error gracefully
          print*,'Problem in set_level_usr_val with levels of Vcode ',Level_vgrid_usr(F_levset)%vcode, ' ', &
               trim(Level_typ_S(F_levset))
          stop
       endif
       if(present(F_is_lnp_L))F_is_lnp_L=is_lnp_L
    elseif(Level_vgrid_usr(F_levset)%vcode == 4001)then
       ! Heights above ground level
       F_stag_S(2:4)="H "//OutGrid_hgrid_usr(F_hgrid_index)%usr_grid_index_S
       if( set_level_usr_eta_lnp(Level_vgrid_usr(F_levset), F_val, F_nko, F_rf, &
            F_kind, is_lnp_L,F_level(1,F_levset),F_level_max(F_levset),&
            F_minx, F_maxx, F_miny, F_maxy) == SET_LEVEL_ERROR )then
          ! TODO handle error gracefully
          print*,'Problem in set_level_usr_val with levels of Vcode ',Level_vgrid_usr(F_levset)%vcode, ' ', &
               trim(Level_typ_S(F_levset))
          stop
       endif
       ier = gmm_get(gmmk_fis0_s, fis0)
       do k=1,F_nko
          F_val(1:l_ni,1:l_nj,k)=F_val(1:l_ni,1:l_nj,k)*grav_8+fis0(1:l_ni,1:l_nj)
       end do
       if(present(F_is_lnp_L))F_is_lnp_L=is_lnp_L
    else
       ! TODO handle error gracefully
       print*,'Problem in set_level_usr_val unsuported Vcode ',Level_vgrid_usr(F_levset)%vcode
       stop
    endif
    
    if(associated(indo_target))then
       if(size(indo_target) < F_nko)deallocate(indo_target)
    endif
    if(.not.associated(indo_target))then
       allocate(indo_target(F_nko), stat = ier)
       if(ier /= 0 )then
          print*,'Allocation problem in set_level_usr_val on indo_target of size =>'
          print*,'F_nko=',F_nko
          return
       endif
    endif
    F_indo(1:F_nko) => indo_target(1:F_nko)
    do k=1,F_nko
       F_indo(k) = k
    end do
    
    istat = set_level_OK
    return    
  end function set_level_usr_val

  integer function set_level_get_usr_src(F_usr_src, F_levset, F_lev_src_name_S, &
       F_minx, F_maxx, F_miny, F_maxy, F_nk) result(istat)
    use glb_ld
    use rmn_gmm, only: gmm_get
    use gmm_pw, only: gmmk_pw_log_pm_s, gmmk_pw_log_pt_s, gmmk_pw_gz_plus_s, pw_gz_plus
    use gmm_geof, only : gmmk_fis0_s, fis0
    use levels, only: Level_vgrid_usr
    use vGrid_Descriptors, only: vgd_get, vgd_print,vgd_associated, VGD_ERROR

    use tdpack, only: grav_8

    implicit none

    real, dimension(:,:,:), pointer, intent(inout)  :: F_usr_src
    integer, intent(in) :: F_levset
    character(len=*), intent(in) :: F_lev_src_name_S
    integer, intent(in) :: F_minx, F_maxx, F_miny, F_maxy, F_nk
    
    ! Local variables
    integer :: ier, k, isize
    real, save, dimension(:), pointer :: usr_src_target
    logical, save :: done_L=.false.    

    istat = set_level_ERROR

    if(.not.done_L)then
       done_L=.true.
       nullify(usr_src_target)
    endif

    if( trim(F_lev_src_name_S) /= 'MOMENTUM' .and. trim(F_lev_src_name_S) /= 'THERMO' )then
       print*,'ERROR in set_level_get_usr_src, argument F_lev_src_name_S must be set'
       print*,'to "MOMENTUM" or "THERMO", got ',trim(F_lev_src_name_S)
       return
    endif
    !ier = vgd_print(Level_vgrid_usr(F_levset)%vgd)
    
    if( Level_vgrid_usr(F_levset)%vcode == 1002 )then
       if( trim(F_lev_src_name_S) == 'MOMENTUM' )then
          ier = gmm_get(gmmk_pw_log_pm_s, F_usr_src)
       elseif( trim(F_lev_src_name_S) == 'THERMO' )then
          ier = gmm_get(gmmk_pw_log_pt_s, F_usr_src)
       endif
    elseif(Level_vgrid_usr(F_levset)%vcode == 2001)then
       if( trim(F_lev_src_name_S) == 'MOMENTUM' )then
          ier = gmm_get(gmmk_pw_log_pm_s, F_usr_src)
       elseif( trim(F_lev_src_name_S) == 'THERMO' )then
          ier = gmm_get(gmmk_pw_log_pt_s, F_usr_src)
       endif
    elseif(Level_vgrid_usr(F_levset)%vcode == 4001)then
       isize=(F_maxx-F_minx+1)*(F_maxy-F_miny+1)*(F_nk+1)
       if(associated(usr_src_target))then
          if(size(usr_src_target) < isize)deallocate(usr_src_target)
       endif
       if(.not.associated(usr_src_target))then         
          allocate(usr_src_target(isize), stat = ier)
          if(ier /= 0 )then
             print*,'Allocation problem in set_level_get_usr_src on usr_src_target of size =>'
             print*,'(F_maxx-F_minx+1)*(F_maxy-F_miny+1)*(F_nk+1)=',isize
             return
          endif
       endif
       F_usr_src(F_minx:F_maxx,F_miny:F_maxy,1:F_nk+1) => usr_src_target(1:isize)
       if( trim(F_lev_src_name_S) == 'MOMENTUM' )then
          ier = gmm_get(gmmk_fis0_s, fis0)
          ier = gmm_get(gmmk_pw_gz_plus_s, pw_gz_plus)
          do k = 1, F_nk
             F_usr_src(1:l_ni,1:l_nj,k) = pw_gz_plus(1:l_ni,1:l_nj,k)
          enddo
          F_usr_src(1:l_ni,1:l_nj,F_nk+1) = fis0(1:l_ni,1:l_nj)+10*grav_8
       elseif( trim(F_lev_src_name_S) == 'THERMO' )then
          ier = gmm_get(gmmk_fis0_s, fis0)
          ier = gmm_get(gmmk_pw_gz_plus_s, pw_gz_plus)
          do k = 1, F_nk-1
             F_usr_src(1:l_ni,1:l_nj,k) = .5*(pw_gz_plus(1:l_ni,1:l_nj,k)+pw_gz_plus(1:l_ni,1:l_nj,k+1))
          enddo
          F_usr_src(1:l_ni,1:l_nj,F_nk)   = .5*(pw_gz_plus(1:l_ni,1:l_nj,F_nk)+fis0(:,:))
          F_usr_src(1:l_ni,1:l_nj,F_nk+1) = fis0(1:l_ni,1:l_nj)
       endif
    else
       ! TODO handle error gracefully
       print*,'Problem in set_level_get_usr_src unsuported Vcode ',Level_vgrid_usr(F_levset)%vcode
       stop
    endif
    istat = set_level_OK
    return

  end function set_level_get_usr_src

  
end module set_level_mod
