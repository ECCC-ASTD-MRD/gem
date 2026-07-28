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

module OUTs
   use vGrid_Descriptors
   use IOs
   use UTL
   use rmn_fst24
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   private :: OUTs_get_source_ezgrid_id
   save

      type :: sfc_var
         sequence
         character(len=4) :: stag
         character(len=4) :: nv
         integer :: knd, nbits, indx, k0, skip
         real :: lvl
         type(UTL_usr_int_info) :: usr_int_info
      end type sfc_var
      type :: fst
         character(len=32  ) :: type
         character(len=1024) :: name
         integer :: fnt
      end type fst
      type(fst) :: OUTs_out(20)
      type(fst_file) :: Out_file
      type(fst_file) :: Out_list_files(20)
      type(fst_record) :: Out_rec
      type :: vgrid_usr
         real(kind=REAL64), pointer :: vtbl_8(:,:,:)
         type(vgrid_descriptor) :: vgd
         integer :: stag,vcode,kind,ni,nj,nk
         real, dimension(:), pointer, contiguous :: levels
      end type vgrid_usr
      type(vgrid_usr), pointer, dimension(:) :: vgd_usr
      type :: hgrid_usr
         integer :: ni, nj, ip1, ip2, ip3, ig1, ig2, ig3, ig4, usr_grid_index
         real, pointer, dimension(:) :: tic,tac
         integer :: ezgdid
         logical :: same_rotation_L
      end type hgrid_usr
      type :: hgrid_source
         character(len=1) :: h_stag_grid_S
         integer :: i0,in,ezgdid         
      end type hgrid_source
      integer, parameter :: Out_max_hgd_src = 100
      integer :: Out_n_hgd_src=0
      type(hgrid_source), dimension(Out_max_hgd_src) :: hgd_src
      type(hgrid_usr), pointer, dimension(:) :: hgd_usr
      type :: wind_rec
         logical :: got_UU_L
         real, dimension(:,:,:), pointer :: data
         real, dimension(:), pointer :: levels
         integer :: dateo,npas,nis,wk_njs,k0,kn
      end type wind_rec
      character(len=2048) :: Out_path_input_S
      character(len=1024) :: Out_dirname_S
      character(len=32  ) :: Out_filenames_S(20)
      character(len=12  ) :: Out_etik_S
      character(len=4   ) :: Out_nomvar
      character(len=1   ) :: Out_typvar_S,chac1,Out_previous_usr_grid_S
      logical :: OUTs_server_L=.false.
      logical :: OUTs_1o1_L, Out_reduc_L, Out_rewrit_L
      integer :: OUTs_GEM_COMM,OUTs_me1o1,OUTs_gem1o1
      integer :: clients,data_wm,IOS_events=0
      integer :: nout_files, Out_kind, Out_nbit
      integer :: Out_i0,Out_in,Out_j0,Out_jn
      integer :: Out_hgd_usr_index, OUTs_hostmyproc
      integer :: Out_ip2,Out_ip3,Out_npas,Out_dateo,Out_deet
      integer :: Out_ig1,Out_ig2,Out_ig3,Out_ig4,Out_unf,Out_previous_src_ezgdid
      integer :: SHARED_WIN
      integer, dimension (:    ), allocatable :: metaG
      real, dimension(:), pointer :: levels
      type(UTL_usr_int_info), dimension(:), pointer :: OUT_usr_int_info
    contains
      
      logical function OUTs_vgrid_usr_L(F_stag) result(status_L)
        character(len=4) :: F_stag
        status_L = (&
             F_stag(2:2)=='1' .or. F_stag(2:2)=='2' .or.&
             F_stag(2:2)=='3' .or. F_stag(2:2)=='4' .or.&
             F_stag(2:2)=='4' .or. F_stag(2:2)=='6' .or.&
             F_stag(2:2)=='7' .or. F_stag(2:2)=='8' .or.&
             F_stag(2:2)=='9' )
        return        
      end function OUTs_vgrid_usr_L
      
      logical function OUTs_hor_int_L(F_stag) result(status_L)
        character(len=4) :: F_stag
        status_L = (&
             F_stag(4:4)=='1' .or. F_stag(4:4)=='2' .or.&
             F_stag(4:4)=='3' .or. F_stag(4:4)=='4' .or.&
             F_stag(4:4)=='4' .or. F_stag(4:4)=='6' .or.&
             F_stag(4:4)=='7' .or. F_stag(4:4)=='8' .or.&
             F_stag(4:4)=='9' )
        return        
      end function OUTs_hor_int_L

      logical function OUTS_read_hgrid_usr(F_hgd_usr) result(status_L)
        use fst_utils
        implicit none
        type(hgrid_usr), target, dimension(:), intent(inout)  :: F_hgd_usr
        !Local variables
        integer  :: i,unit,ier,key,ni,nj,nk
        integer, external :: fnom,fstouv,fstfrm,fclos,fstinf,fstluk
        character(len=1) :: index_S
        character(len=2064) :: file_S
        type(hgrid_usr), pointer :: fh
        type(fst_rpn) :: rec
        status_L=.false.
        ! TODO get a better way to get unit no
        unit=62
        do i=1,size(F_hgd_usr)
           fh => F_hgd_usr(i) ! just to simplify code
           write(index_S,'(i1)')fh%usr_grid_index
           file_S=trim(Out_path_input_S)//'/MODEL_INPUT/user_hgrid'//index_S
           if( fnom(unit,file_S,'RND+OLD+R/O',0) < 0 )then
              print*,'Error in OUTS_read_hgrid_usr with fnom on file ',trim(file_S)
              return
           endif
           if( fstouv(unit,'RND') < 0 )then
              print*,'Error in OUTS_read_hgrid_usr with fstouv on file ',trim(file_S)
              return
           endif
           key=fstinf(unit,fh%ni,nj,nk,-1,' ',-1,-1,-1,' ','>>')
           if(key < 0)then
              print*,'Error in OUTS_read_hgrid_usr, no >> found in file ',trim(file_S)
              return
           endif
           if(fst_fstprm(key,rec) == FST_ERROR)return         
           allocate(fh%tic(fh%ni))
           if( fstluk(fh%tic, key, fh%ni, nj, nk) < 0 )then
              print*,'Error in OUTS_read_hgrid_usr allocating user grid tic'
              return
           endif
           key=fstinf(unit,ni,fh%nj,nk,-1,' ',-1,-1,-1,' ','^^')
           if(key < 0)then
              print*,'Error in OUTS_read_hgrid_usr no ^^ found in ',trim(file_S)
              return
           endif
           if( fst_fstprm(key,rec) == FST_ERROR)return
           allocate( fh%tac(fh%nj))
           if( fstluk(fh%tac, key, ni, fh%nj, nk) < 0 )then
              print*,'Error in OUTS_read_hgrid_usr allocating user grid tac'
              return
           endif
           fh%ip1=rec%ip1; fh%ip2=rec%ip2; fh%ip3=rec%ip3
           fh%ig1=rec%ig1; fh%ig2=rec%ig2; fh%ig3=rec%ig3; fh%ig4=rec%ig4
           ier=fstfrm(unit); ier=fclos(unit)
           unit=unit+1
        end do
        status_L=.true.
      end function OUTS_read_hgrid_usr

      logical function OUTS_read_vgrid_usr(F_vgd_usr) result(status_L)
        use vGrid_Descriptors
        implicit none        
        type(vgrid_usr), target, dimension(:), intent(inout)  :: F_vgd_usr        
        !Local variables
        integer, external :: fnom,fstouv,fstfrm,fclos
        integer :: i,unit,ier
        character(len=1) :: index_S
        character(len=2064) :: file_S
        type(vgrid_usr), pointer :: fv
        status_L=.false.
        ! TODO get a better way to get unit no
        unit=72
        do i=1,size(F_vgd_usr)
           fv => F_vgd_usr(i) ! just to simplify code           
           write(index_S,'(i1)')fv%stag
           file_S=trim(Out_path_input_S)//'/MODEL_INPUT/user_vgrid'//index_S
           if( fnom(unit,file_S,'RND+OLD+R/O',0) < 0 )then
              print*,'Error in OUTS_read_vgrid_usr with fnom on file ',trim(file_S)
              return
           endif
           if( fstouv(unit,'RND') < 0 )then
              print*,'Error in OUTS_read_vgrid_usr with fstouv on file ',trim(file_S)
              return
           endif
           if( vgd_new(fv%vgd,unit) == VGD_ERROR )then
              print*,'Error in OUTS_read_vgrid_usr with vgd_new on file ',trim(file_S)
              return
           endif
           !ier = vgd_print(fv%vgd)
           if( vgd_get (fv%vgd,'VTBL',fv%vtbl_8) == VGD_ERROR )return
           fv%ni=size(fv%vtbl_8,1); fv%nj=size(fv%vtbl_8,2); fv%nk=size(fv%vtbl_8,2);
           ier=fstfrm(unit); ier=fclos(unit)
           unit=unit+1
        end do
        status_L=.true.
      end function OUTS_read_vgrid_usr
      
      logical function OUTs_set_horizontal_interpolation(F_stag_S, F_src_ezgdid) result(found_L)
        implicit none
        character(len=4), intent(in) ::F_stag_S
        integer, intent(in) :: F_src_ezgdid
        ! Local variables
#include <rmnlib_basics.hf>
        integer :: i,h_stag,gdid
        !Find user horizontal grid
        found_L=.false.
        read(F_stag_S(4:4),*)h_stag
        do i=1,size(hgd_usr)
           if(hgd_usr(i)%usr_grid_index == h_stag)then
              Out_hgd_usr_index=i
              found_L=.true.
              exit
           endif
        end do
        if(.not.found_L)then
           write(Lun_out,'("OUTs_set_horizontal_interpolation, internale error could not find matching user horizontal grid for stag ",s)')F_stag_S
           return
        endif
        gdid = ezdefset(hgd_usr(Out_hgd_usr_index)%ezgdid,F_src_ezgdid)
      end function OUTs_set_horizontal_interpolation
      
      logical function OUTs_get_source_ezgrid_id(F_ezgdid,F_stag_S) result(status_L)
        implicit none
        integer, intent(out) :: F_ezgdid
        character(len=1), intent(in) ::F_stag_S
#include <rmnlib_basics.hf>
        integer :: nis,njs,i
        integer, dimension(2) :: YY_gdid
        real, dimension(:), pointer :: posx,posy
        ! Try to find existing grid having the same definition
        status_L=.true.
        do i=1,Out_n_hgd_src
           if( F_stag_S(1:1) == hgd_src(i)%h_stag_grid_S .and. &
                Out_i0 == hgd_src(i)%i0 .and. &
                Out_in == hgd_src(i)%in ) then
              F_ezgdid=hgd_src(i)%ezgdid
              return
           endif
        end do
        ! Grid not defined yet, define it
        if (F_stag_S(1:1) =='M') then
           posx => geomh_longs           
           posy => geomh_latgs
        end if
        if (F_stag_S(1:1) =='U') then
           ! Note that UU is on the model grid, not on stag U
           posx => geomh_longu
           posy => geomh_latgs
        end if
        if (F_stag_S(1:1) =='V') then
           ! Note that VV is on the model grid, not on stag V
           posx => geomh_longs
           posy => geomh_latgv
        end if
        if (F_stag_S(1:1) =='F') then
           posx => geomh_longu
           posy => geomh_latgv
        end if
        nis=Out_in-Out_i0+1; njs=Out_jn-Out_j0+1
        ! Lam (model) grid or YIN grid
        YY_gdid(1) = ezgdef_fmem(nis,njs,'Z','E',Rot_ig1,Rot_ig2,Rot_ig3,Rot_ig4,posx(Out_i0),posy(Out_j0))     
        if ( Grd_yinyang_L ) then !Yin-Yang
           ! Yan grid
           ! model
           YY_gdid(2) = ezgdef_fmem(nis, njs, 'Z', 'E', RotY_ig1, RotY_ig2, RotY_ig3, RotY_ig4, posx(Out_i0), posy(Out_j0))
           F_ezgdid = ezgdef_supergrid(nis, njs*2,'U','F',1,2,YY_gdid)
        else
           F_ezgdid = YY_gdid(1)
        endif
        if(Out_n_hgd_src+1 > Out_max_hgd_src)then
           if (Lun_out>0)then
              write(Lun_out,'("WARNING!! Too many source grid definitions, skipping current and next ones.")')
              status_L=.false.
              return
           endif
        endif
        Out_n_hgd_src=Out_n_hgd_src+1
        hgd_src(Out_n_hgd_src)%h_stag_grid_S=F_stag_S(1:1)
        hgd_src(Out_n_hgd_src)%i0=Out_i0; hgd_src(Out_n_hgd_src)%in=Out_in
        hgd_src(Out_n_hgd_src)%ezgdid=F_ezgdid
      end function OUTs_get_source_ezgrid_id

      subroutine OUTs_wrtref_usr(F_stag_S,F_f)        
        implicit none
        character(len=4), intent(IN) :: F_stag_S
        type(hgrid_usr) :: F_f
        ! Local variables
#include <rmnlib_basics.hf>
        integer :: err
        real :: wk
        err=fstecr(F_f%tic,wk,-32,Out_unf,Out_dateo   ,&
             0,0,F_f%ni,1,1,F_f%ip1,F_f%ip2,F_f%ip3,'X','>>',&
             Out_etik_S,'E',&
             F_f%ig1,F_f%ig2,F_f%ig3,F_f%ig4,5,.true.)
        err=fstecr(F_f%tac,wk,-32,Out_unf,Out_dateo   ,&
                    0,0,1,F_f%nj,1,F_f%ip1,F_f%ip2,F_f%ip3,'X', '^^'    ,&
                    Out_etik_S,'E',&
                    F_f%ig1,F_f%ig2,F_f%ig3,F_f%ig4,5, .true.)
        call OUTs_wrtvref(F_stag_S,F_f%ip1,F_f%ip2)
      end subroutine OUTs_wrtref_usr
      
      subroutine OUTs_wrtvref(F_stag_S,F_ig1,F_ig2)
        ! Write vertical grid descritor !!
        ! To respect the CMC legacy, the !! is not writen in eta files
        implicit none
        character(len=4), intent(IN) :: F_stag_S
        integer, intent(IN) :: F_ig1,F_ig2
        ! Local variables
        integer :: c1,c2,k,err,vcode
        integer, dimension(:), allocatable :: ip1s
        real(kind=REAL64), dimension(:), allocatable :: zero
        type(vgrid_descriptor) :: vgd
        err = vgd_get(gem_vgd,'VCOD',vcode)
        if(vcode == 1002)return
        if ( F_stag_S(2:2) == 'M' ) then
           err = vgd_put  (gem_vgd,'IP_1 - record ip1',F_ig1)
           err = vgd_put  (gem_vgd,'IP_2 - record ip2',F_ig2)
           err = vgd_write(gem_vgd,unit=Out_unf,format='fst')
        else if ( F_stag_S(2:2) == 'P') then
           allocate(ip1s(size(Level_allpres)))
           allocate(zero(size(Level_allpres)))
           zero = 0.d0 ; c1=2 ; c2=1
            do k=1,size(Level_allpres)
               call convip(ip1s(k),Level_allpres(k),c1,c2,'',.false.)
            end do
            err = vgd_new(vgd,            &
                 kind     = 2,              &
                 version  = 1,              &
                 nk       = size(Level_allpres),     &
                 ip1      = F_ig1,        &
                 ip2      = F_ig2,        &
                 a_m_8    = dble(Level_allpres*100.),&
                 b_m_8    = zero,           &
                 ip1_m    = ip1s)
            err = vgd_write(vgd,unit=Out_unf,format='fst')
            err = vgd_free(vgd)
         else if ( F_stag_S(2:2) == 'H' ) then
            allocate(ip1s(size(Level_allheights)))
            allocate(zero(size(Level_allheights)))
            zero = 0.d0 ; c1=4 ; c2=1
            do k=1,size(Level_allheights)
               call convip(ip1s(k),Level_allheights(k),c1,c2,'',.false.)
            end do
            err = vgd_new(vgd,&
                 kind     = 4,&
                 version  = 1,&
                 nk       = size(Level_allheights),&
                 ip1      = F_ig1,&
                 ip2      = F_ig2,&
                 a_m_8    = dble(Level_allheights),&
                 b_m_8    = zero,&
                 ip1_m    = ip1s)
            err = vgd_write(vgd,unit=Out_unf,format='fst')
            err = vgd_free(vgd)
         endif

       end subroutine OUTs_wrtvref
      
      subroutine OUTS_init_wind_rec(F_wind_rec)
        type(wind_rec) :: F_wind_rec
        nullify(F_wind_rec%data)
        nullify(F_wind_rec%levels)
        F_wind_rec%got_UU_L=.false.
        F_wind_rec%dateo=-1; F_wind_rec%npas=-1; F_wind_rec%nis=-1; F_wind_rec%wk_njs=-1
        F_wind_rec%k0=-999;F_wind_rec%kn=-999
      end subroutine OUTS_init_wind_rec
      
      logical function OUTs_hor_int(wk2,wk2_uu,F_wk,F_nomvar_S,F_stag_S,F_windu,F_nis,F_wk_njs,&
           F_k0,F_kn, F_levels,F_gtyp,F_ig1,F_ig2,F_ig3,F_ig4,F_etiket_S,F_vector_int_L, &
           F_wait_for_v_L,F_skip_wind_L) result(status_L)
        implicit none
        integer, intent(IN) :: F_k0,F_kn
        integer, intent(INOUT) :: F_nis,F_wk_njs,F_ig1,F_ig2,F_ig3,F_ig4
        real, intent(INOUT), dimension(F_nis,F_wk_njs,F_k0:F_kn) :: F_wk
        real, intent(INOUT), dimension(:,:,:), pointer :: wk2,wk2_uu
        real, intent(IN) :: F_levels(F_k0:F_kn)
        character(len=1), intent(INOUT) :: F_gtyp
        character(len=4), intent(IN) :: F_nomvar_S, F_stag_S
        character(len=12), intent(INOUT) ::F_etiket_S
        Logical :: F_vector_int_L,F_wait_for_v_L
        logical, optional :: F_skip_wind_L
        type(wind_rec), intent(INOUT) :: F_windu
        ! Local variables
        integer, external :: ezuvint,ezsint,ezsetopt
        integer :: k,ier
        logical :: skip_wind_L
        logical, save :: done_L = .false.
        integer :: current_src_ezgdid
        real, save, dimension(:), pointer :: target_wk2,target_wk2_uu,target_uu,target_levels
        type(UTL_usr_int_info) :: usr_int_info
        if(.not.done_L)then
           done_L=.true.
           nullify(target_wk2,target_wk2_uu,target_uu,target_levels)
        endif
        skip_wind_L=.false.
        ! TODO remove F_skip_wind_L
        if(present(F_skip_wind_L))skip_wind_L=F_skip_wind_L
        status_L=.false.        
        F_vector_int_L=.false.
        F_wait_for_v_L=.false.
        if(.not.UTL_get_usr_int_info(usr_int_info,F_nomvar_S,OUT_usr_int_info))then
           print*,'Error in OUTs_process_data with, UTL_get_usr_int_info'
           ! TODO handle error gracefully
        endif
        
        ! Note, user grid is fully itentified by F_stag_S(4:4)
        !       However, source grid having the same F_stag_S(1:1) may differ due to their i0 and in.
        !       So we have to compare their ez grid id.
        if(.not.OUTs_get_source_ezgrid_id(current_src_ezgdid,F_stag_S))return
        if(  current_src_ezgdid /= Out_previous_src_ezgdid .or. &
             F_stag_S(4:4) /= Out_previous_usr_grid_S )then
           if(.not. OUTs_set_horizontal_interpolation(F_stag_S,current_src_ezgdid))return
           Out_previous_src_ezgdid=current_src_ezgdid
           Out_previous_usr_grid_S=F_stag_S(4:4)
        endif
        if(trim(F_nomvar_S) == "UU" .and. &
             (.not. hgd_usr(Out_hgd_usr_index)%same_rotation_L) .and. &
             .not. skip_wind_L)then
            ! For the wind interpolation on rotated grid we suppose that UU appears in the stack before VV
            ! otherwise the wind will not be output.
            ! Save UU wind component for vector interpolation
           call OUTs_check_allocation_3D(F_windu%data,target_uu,F_nis,F_wk_njs,F_k0,F_kn)
           call OUTs_check_allocation_1D(F_windu%levels,target_levels,F_k0,F_kn)
           
            ! Note k range is specified in copy in case F_windu%data and F_windu%levels are reused and
            ! have a larger k scope than wk
            F_windu%data(:,:,F_k0:F_kn)=F_wk(:,:,F_k0:F_kn)
            F_windu%levels(F_k0:F_kn)=F_levels(F_k0:F_kn)
            F_windu%nis=F_nis; F_windu%wk_njs=F_wk_njs
            F_windu%got_UU_L=.true.
            F_windu%dateo=Out_dateo; F_windu%npas=Out_npas; F_windu%nis=F_nis; F_windu%wk_njs=F_wk_njs
            F_windu%k0=F_k0; F_windu%kn=F_kn
            ! Now that record UU is loaded routine must return
            F_wait_for_v_L=.true.
            status_L=.true.
            return
         endif
         if(trim(F_nomvar_S) == "VV" .and. &
              (.not. hgd_usr(Out_hgd_usr_index)%same_rotation_L) .and. &
              .not. skip_wind_L)then
            F_vector_int_L=.true.
            ! Check if UU was saved before
            if(F_windu%got_UU_L)then
               ! Can procede with vector interpolation
               ! Test if UU parameters fit VV's
               if(F_windu%nis /= F_nis .or. F_windu%wk_njs /= F_wk_njs)then
                  print*,'WARNING OUTs_hor_int, UU and VV paires do not have the same horizontal size, skipping these records'
                  return
               endif
               if(F_windu%k0 /= F_k0 .or. F_windu%kn /= F_kn)then
                  print*,'WARNING OUTs_hor_int, levels scopes are not equivalent for UU and VV, skipping these records'
                  return
               endif
               if(.not.OUTs_levels_are_equivalent(F_windu,F_levels,F_k0,F_kn))then
                  print*,'WARNING OUTs_hor_int, levels are not equivalent for UU and VV, skipping these records'
                  return
               endif
               if(F_windu%dateo /= Out_dateo .or. F_windu%npas /= Out_npas) then
                  print*,'WARNING OUTs_hor_int, UU and VV paires do not have the same valid time, skipping these records'
                  return
               endif
            else
               print*,'WARNING OUTs_hor_int, got VV but did not see UU, make sure UU appears first in output request'
               return
            endif
         endif            
         ! Do interpolation
         ier = ezsetopt ('INTERP_DEGREE', usr_int_info%hor%int_deg_S)
         call OUTs_check_allocation_3D(wk2,target_wk2,hgd_usr(Out_hgd_usr_index)%ni,hgd_usr(Out_hgd_usr_index)%nj, &
              F_k0,F_kn)
         
         if(F_vector_int_L)then
            call OUTs_check_allocation_3D(wk2_uu,target_wk2_uu,hgd_usr(Out_hgd_usr_index)%ni,&
                 hgd_usr(Out_hgd_usr_index)%nj, F_k0,F_kn)
            do k=F_k0, F_kn
               if(ezuvint(wk2_uu(1,1,k), wk2(1,1,k), F_windu%data(1,1,k), F_wk(1,1,k)) < 0)then
                  print*,'ERROR in horizontal vector interpolation in OUTs_hor_int'
                  ! TODO handle error gracefully
                  ! TODO trap EXTRAP ABORT gracefully when version of ezscint returns error value
                  stop
               endif
            end do
            ! Reset F_windu for next UU,VV output
            ! Don't deallocate in order to try reuse arrays
            F_windu%got_UU_L=.false.
            F_windu%dateo=-1; F_windu%npas=-1; F_windu%nis=-1; F_windu%wk_njs=-1
            F_windu%k0=-999;F_windu%kn=-999
         else
            do k=F_k0, F_kn
               if( ezsint(wk2(1,1,k), F_wk(1,1,k)) < 0)then
                  print*,'ERROR in horizontal interpolation in OUTs_hor_int'
                  ! TODO handle error gracefully
                  ! TODO trap EXTRAP ABORT gracefully when version of ezscint returns error value
                  stop
               endif
            end do
            call apply_limits_3D(wk2,F_k0,F_kn,usr_int_info)
         endif
         F_nis=hgd_usr(Out_hgd_usr_index)%ni
         F_wk_njs=hgd_usr(Out_hgd_usr_index)%nj
         F_gtyp='Z'
         F_ig1=hgd_usr(Out_hgd_usr_index)%ip1; F_ig2=hgd_usr(Out_hgd_usr_index)%ip2;&
              F_ig3=hgd_usr(Out_hgd_usr_index)%ip3; F_ig4=0
         status_L=.true.
       end function OUTs_hor_int

       logical function OUTs_hor_int_sfc(F_wk2,F_wk,F_vsfc,F_nis,F_wk_njs,&
            F_gtyp,F_ig1,F_ig2,F_ig3,F_ig4,F_etiket_S,wind_cycle_L) result(status_L)         
         implicit none
         integer, intent(INOUT) :: F_nis,F_wk_njs,F_ig1,F_ig2,F_ig3,F_ig4
         real, intent(INOUT), dimension(F_nis,F_wk_njs) :: F_wk
         real, intent(INOUT), dimension(:,:), pointer :: F_wk2
         type(sfc_var) :: F_vsfc
         character(len=1) :: F_gtyp         
         character(len=12), intent(INOUT) ::F_etiket_S
         logical, intent(INOUT) :: wind_cycle_L
         ! Local variables
         integer :: current_src_ezgdid,ier
         integer, external :: ezsint, ezsetopt
         logical, save :: done_L=.false., wind_warning_done_L=.false.
         real, save, dimension(:), pointer :: target_wk2         
         if(.not.done_L)then
            done_L=.true.
            nullify(target_wk2)
         endif
         status_L=.false.
         if(.not.OUTs_get_source_ezgrid_id(current_src_ezgdid,F_vsfc%stag))return
         if(  current_src_ezgdid /= Out_previous_src_ezgdid .or. &
              F_vsfc%stag(4:4) /= Out_previous_usr_grid_S )then
            if(.not. OUTs_set_horizontal_interpolation(F_vsfc%stag,current_src_ezgdid))return
            Out_previous_src_ezgdid=current_src_ezgdid
            Out_previous_usr_grid_S=F_vsfc%stag(4:4)
         endif
         wind_cycle_L=.false.
         if(.not.hgd_usr(Out_hgd_usr_index)%same_rotation_L .and. &
              (trim(F_vsfc%nv) == 'UU' .or. trim(F_vsfc%nv) == 'VV')) then
            ! Note: with the actual surface output parallelisation, it is not possible to
            ! be certain that UU and VV will be processed by the same PE. Therefore we
            ! cannot do the vector wind horizontal interpolation that are needed for user
            ! grid not having the same rotation as source grid. For such cases we print
            ! a warning
            wind_warning_done_L=.true.
            if( Lun_out>0.and. (.not.wind_warning_done_L) )then
               write(Lun_out,'("Skipping diag level wind output since user and model grids have a different rotation")')
               status_L=.true.
            endif
            wind_cycle_L=.true.
            return
         endif         
         call OUTs_check_allocation_2D(F_wk2,target_wk2,hgd_usr(Out_hgd_usr_index)%ni,hgd_usr(Out_hgd_usr_index)%nj)

         ier = ezsetopt ('INTERP_DEGREE', F_vsfc%usr_int_info%hor%int_deg_S)
         if( ezsint(F_wk2, F_wk) < 0)then
            print*,'ERROR in horizontal interpolation in OUTs_hor_int'
            ! TODO handle error gracefully
            ! TODO trap EXTRAP ABORT gracefully when version of ezscint returns error value
            stop
         endif
         call apply_limits_2D(F_wk2,F_vsfc%usr_int_info)
         F_nis=hgd_usr(Out_hgd_usr_index)%ni
         F_wk_njs=hgd_usr(Out_hgd_usr_index)%nj
         F_gtyp='Z'
         F_ig1=hgd_usr(Out_hgd_usr_index)%ip1; F_ig2=hgd_usr(Out_hgd_usr_index)%ip2;&
              F_ig3=hgd_usr(Out_hgd_usr_index)%ip3; F_ig4=0
         status_L=.true.
       end function OUTs_hor_int_sfc

       subroutine apply_limits_2D(F_f,F_usr_int_info)
         implicit none
         real, dimension(:,:), pointer :: F_f
         type(UTL_usr_int_info) :: F_usr_int_info
         ! Local variables
         integer :: i,j
         if(F_usr_int_info%hor%low_lim_L .and. F_usr_int_info%hor%high_lim_L )then
            do j=1,size(F_f,2)
               do i=1,size(F_f,1)
                  F_f(i,j)=max(F_usr_int_info%hor%low_lim ,F_f(i,j))
                  F_f(i,j)=min(F_usr_int_info%hor%high_lim,F_f(i,j))
               end do
            end do
         else
            if(F_usr_int_info%hor%low_lim_L)then
               do j=1,size(F_f,2)
                  do i=1,size(F_f,1)
                     F_f(i,j)=max(F_usr_int_info%hor%low_lim ,F_f(i,j))
                  end do
               end do
            endif
            if(F_usr_int_info%hor%high_lim_L)then
               do j=1,size(F_f,2)
                  do i=1,size(F_f,1)
                     F_f(i,j)=min(F_usr_int_info%hor%high_lim,F_f(i,j))
                  end do
               end do
            endif
         endif
       end subroutine apply_limits_2D
       
       subroutine apply_limits_3D(F_f,F_k0,F_kn,F_usr_int_info)
         implicit none
         real, dimension(:,:,:), pointer :: F_f
         integer :: F_k0,F_kn
         type(UTL_usr_int_info) :: F_usr_int_info
         ! Local variables
         integer :: i,j,k
         if(F_usr_int_info%hor%low_lim_L .and. F_usr_int_info%hor%high_lim_L )then
            do k=F_k0,F_kn
               do j=1,size(F_f,2)
                  do i=1,size(F_f,1)
                     F_f(i,j,k)=max(F_usr_int_info%hor%low_lim ,F_f(i,j,k))
                     F_f(i,j,k)=min(F_usr_int_info%hor%high_lim,F_f(i,j,k))
                  end do
               end do
            end do
         else
            if(F_usr_int_info%hor%low_lim_L)then
               do k=F_k0,F_kn
                  do j=1,size(F_f,2)
                     do i=1,size(F_f,1)
                        F_f(i,j,k)=max(F_usr_int_info%hor%low_lim ,F_f(i,j,k))
                     end do
                  end do
               end do
            endif
            if(F_usr_int_info%hor%high_lim_L)then
               do k=F_k0,F_kn
                  do j=1,size(F_f,2)
                     do i=1,size(F_f,1)
                        F_f(i,j,k)=min(F_usr_int_info%hor%high_lim,F_f(i,j,k))
                     end do
                  end do
               end do
            endif
         endif
       end subroutine apply_limits_3D
       
    logical function OUTs_levels_are_equivalent(windu,lv2,k0,kn) result(equal_L)
      implicit none      
      integer, intent(IN) :: k0,kn
      real, intent(IN), dimension(k0:kn) :: lv2
      type(wind_rec), intent(IN) :: windu
      ! Local variables
      integer :: k
      
      equal_L=.false.
      do k=k0,kn
         if(windu%levels(k) == 0.)then
            if(abs(windu%levels(k)-lv2(k)) > 1.e-5)return
         else
            if(abs(windu%levels(k)-lv2(k))/lv2(k) > 1.e-5)return
         end if
      end do
      equal_L=.true.
      return
    end function OUTs_levels_are_equivalent

    subroutine OUTs_check_allocation_1D(F_data,F_target,F_k0,F_kn)
      implicit none
      integer, intent(IN) :: F_k0,F_kn
      real, dimension(:), pointer, intent(INOUT) :: F_data,F_target
      ! Local variables
      integer :: nk, n_target
      nk=F_kn-F_k0+1
      if(associated(F_target))then
         n_target=size(F_target)
         if(n_target < nk)then
            deallocate(F_target)
            allocate(F_target(nk))
         endif
      else
         allocate(F_target(nk))
      endif
      F_data(F_k0:F_kn) => F_target(1:nk)
    end subroutine OUTs_check_allocation_1D
    
    subroutine OUTs_check_allocation_2D(F_data_2D,F_target,F_ni,F_nj)
      implicit none
      integer, intent(IN) :: F_ni,F_nj
      real, dimension(:,:), pointer, intent(INOUT) :: F_data_2D
      real, dimension(:), pointer, intent(INOUT) :: F_target
      ! Local variables
      integer :: nij, n_target
      nij=F_ni*F_nj
      if(associated(F_target))then
         n_target=size(F_target)
         if(n_target < nij)then
            deallocate(F_target)
            allocate(F_target(nij))
         endif
      else
         allocate(F_target(nij))
      endif
      F_data_2D(1:F_ni,1:F_nj) => F_target(1:nij)
    end subroutine OUTs_check_allocation_2D
      
    subroutine OUTs_check_allocation_3D(F_data_3D,F_target,F_ni,F_nj,F_k0,F_kn)
      implicit none
      integer, intent(IN) :: F_ni,F_nj,F_k0,F_kn
      real, dimension(:,:,:), pointer, intent(INOUT) :: F_data_3D
      real, dimension(:), pointer, intent(INOUT) :: F_target
      ! Local variables
      integer :: nijk, n_target
      nijk=F_ni*F_nj*(F_kn-F_k0+1)
      if(associated(F_target))then
         n_target=size(F_target)
         if(n_target < nijk)then
            deallocate(F_target)
            allocate(F_target(nijk))
         endif
      else
         allocate(F_target(nijk))
      endif
      F_data_3D(1:F_ni,1:F_nj,F_k0:F_kn) => F_target(1:nijk)
    end subroutine OUTs_check_allocation_3D
    
end module OUTs
