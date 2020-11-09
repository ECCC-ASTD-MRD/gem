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

!**s/r itf_phy_init - Initializes physics parameterization package
!
      subroutine itf_phy_init
      use vGrid_Descriptors, only: vgrid_descriptor,vgd_get,vgd_put,vgd_free,VGD_OK,VGD_ERROR
      use vgrid_wb, only: vgrid_wb_get, vgrid_wb_put
      use phy_itf, only: phy_init, phymeta,phy_getmeta
      use itf_phy_filter, only: ipf_init, sfcflxfilt_o, sfcflxfilt_i, nsurfag
      use step_options
      use inp_options
      use init_options
      use dyn_fisl_options
      use adz_options, only: Adz_slt_winds
      use ctrl
      use glb_ld
      use cstv
      use lun
      use out3
      use levels
      use outp
      use ver
      use gmm_itf_mod
      use rstr
      use var_gmm
      use gmm_pw, only: gmmk_pw_uslt_s, gmmk_pw_vslt_s
      use path
      use clib_itf_mod
      use wb_itf_mod
      use ptopo_utils, only: ptopo_io_set !#TODO: should the phyics define its own?
      use dcmip_options, only: dcmip_case
      implicit none
#include <arch_specific.hf>


#include <msg.h>
#include <rmnlib_basics.hf>

      character(len=32), parameter  :: VGRID_M_S = 'ref-m'
      character(len=32), parameter  :: VGRID_T_S = 'ref-t'

      type(vgrid_descriptor) :: vcoord, vcoordt
      integer err,zuip,ztip,vtype
      integer, dimension(:), pointer :: ip1m, ip1t
      real :: zu,zt,cond_sig,gwd_sig,lhn_sig
      real, dimension(:,:), pointer :: ptr2d

      integer,parameter :: tlift= 0
      integer,parameter :: NBUS = 4
      character(len=9) :: BUS_LIST_S(NBUS) = &
                  (/'Entry    ', 'Permanent', 'Volatile ', 'Dynamics '/)
      character(len=GMM_MAXNAMELENGTH) :: diag_prefix
      logical :: have_slt_winds

!For sorting the output
      integer istat, iverb
      character(len=32) :: varname_S,outname_S,bus0_S,refp0_S,refp0ls_S
      character(len=32) :: varlist_S(MAXELEM*MAXSET), input_type_S
      integer k,j,n,ibus,multxmosaic
      type(phymeta) :: pmeta
!
!     ---------------------------------------------------------------
!
! Create space for diagnostic level values
      diag_prefix = 'diag/'
      gmmk_diag_tt_s = trim(diag_prefix)//'TT'
      gmmk_diag_hu_s = trim(diag_prefix)//'HU'
      gmmk_diag_uu_s = trim(diag_prefix)//'UU'
      gmmk_diag_vv_s = trim(diag_prefix)//'VV'
      istat = GMM_OK
      nullify(ptr2d)
      istat = min(gmm_create(gmmk_diag_tt_s, ptr2d ,meta2d),istat)
      nullify(ptr2d)
      istat = min(gmm_create(gmmk_diag_hu_s, ptr2d ,meta2d),istat)
      nullify(ptr2d)
      istat = min(gmm_create(gmmk_diag_uu_s, ptr2d ,meta2d),istat)
      nullify(ptr2d)
      istat = min(gmm_create(gmmk_diag_vv_s, ptr2d ,meta2d),istat)
      if (GMM_IS_ERROR(istat)) &
           call msg(MSG_ERROR,'itf_phy_init ERROR at gmm_create('//trim(diag_prefix)//'*)')

      Out3_sfcdiag_L= .false.

      ! Continue only if the physics is being run
      if (.not.Ctrl_phyms_L) return

      if (Lun_out > 0) write(Lun_out,1000)

! Collect the list of potentialy requested output physics vars
      n = 0
      varlist_S(:) = ' '
      do k = 1, Outp_sets
         do j = 1, Outp_var_max(k)
            varname_S = Outp_varnm_S(j,k)
            istat = clib_tolower(varname_S)
            if (.not.any(varname_S == varlist_S(1:n))) then
               n = n + 1
               varlist_S(n) = varname_S
            end if
         end do
      end do
      !#TODO: trim the output list to the ones actually requested within the run
      if (n > 0) then
         err = wb_put('itf_phy/PHYOUT', varlist_S(1:n), WB_REWRITE_AT_RESTART)
      end if

! We put mandatory variables in the WhiteBoard

      err= 0
      err= min(wb_put('itf_phy/TLIFT'       , tlift        , WB_REWRITE_AT_RESTART), err)
      err= min(wb_put('itf_phy/DYNOUT'      , Out3_accavg_L, WB_REWRITE_AT_RESTART), err)
      err= min(wb_put('itf_phy/slt_winds'   , Adz_slt_winds, WB_REWRITE_AT_RESTART), err)

! Complete physics initialization (see phy_init for interface content)

      istat = ptopo_io_set(Inp_npes) !#TODO mv this in phy... pass npes as arg
      call timing_start2(41, 'PHY_init', 40 )
      err= phy_init ( Path_phy_S, Step_CMCdate0, real(Cstv_dt_8), &
                      'model/Hgrid/lclphy', 'model/Hgrid/lclcore', &
                      'model/Hgrid/global', 'model/Hgrid/local'  , &
                      'model/Hgrid/glbphy', 'model/Hgrid/glbphycore', &
                      G_nk+1, Ver_std_p_prof%m)
      call timing_stop(41)

! Option consistency check
      if (WB_IS_OK(wb_get('phy/input_type', input_type_S))) then
         istat = clib_toupper(input_type_S)
         if (Iau_interval >0. .and. input_type_S /= 'DIST') then
            err = RMN_ERR
            call msg(MSG_ERROR, '(itf_phy_init) with IAU the physics input type should be DIST.')
         endif
      end if

! Initialize filter weights for smoothing
      if (.not.WB_IS_OK(wb_get('phy/cond_infilter',cond_sig))) cond_sig=-1.
      if (.not.WB_IS_OK(wb_get('phy/sgo_tdfilter',gwd_sig))) gwd_sig=-1.
      if (.not.WB_IS_OK(wb_get('phy/lhn_filter',lhn_sig))) lhn_sig=-1.
      if (.not.WB_IS_OK(wb_get('phy/sfcflx_filter_order',sfcflxfilt_o))) sfcflxfilt_o=-1
      if (.not.WB_IS_OK(wb_get('phy/sfcflx_filter_iter',sfcflxfilt_i))) sfcflxfilt_i=1
      if (.not.WB_IS_OK(wb_get('phy/nsurfag',nsurfag))) nsurfag=1

      err = min(ipf_init(F_sig=cond_sig, F_sig2=gwd_sig, F_sig3=lhn_sig), err)
      err = min(wb_put('dyn/cond_infilter',cond_sig,WB_REWRITE_AT_RESTART), err)
      err = min(wb_put('dyn/sgo_tdfilter',gwd_sig,WB_REWRITE_AT_RESTART), err)
      err = min(wb_put('dyn/lhn_filter',lhn_sig,WB_REWRITE_AT_RESTART), err)
      err = min(wb_put('dyn/sfcflx_filter_order',sfcflxfilt_o,WB_REWRITE_AT_RESTART), err)

! Retrieve the heights of the diagnostic levels (thermodynamic
! and momentum) from the physics ( zero means NO diagnostic level)

      iverb = wb_verbosity(WB_MSG_FATAL)
      err= min(wb_get('phy/zu', zu), err)
      err= min(wb_get('phy/zt', zt), err)
      iverb = wb_verbosity(iverb)

      call gem_error ( err,'itf_phy_init','phy_init or WB_get' )

      err = VGD_OK
      if ((zu > 0.) .and. (zt > 0.) ) then
         nullify(ip1m, ip1t)
         Level_kind_diag=4
         err = min ( vgrid_wb_get(VGRID_M_S,vcoord, ip1m,vtype,refp0_S,refp0ls_S), err)
         err = min ( vgrid_wb_get(VGRID_T_S,vcoordt,ip1t), err)
         deallocate(ip1m, ip1t,stat=err) ; nullify(ip1m, ip1t)
         call convip(zuip,zu,Level_kind_diag,+2,'',.true.)
         call convip(ztip,zt,Level_kind_diag,+2,'',.true.)
         err = min(vgd_put(vcoord, 'DIPM - IP1 of diagnostic level (m)',zuip), err)
         err = min(vgd_put(vcoord, 'DIPT - IP1 of diagnostic level (t)',ztip), err)
         err = min(vgd_put(vcoordt,'DIPM - IP1 of diagnostic level (m)',zuip), err)
         err = min(vgd_put(vcoordt,'DIPT - IP1 of diagnostic level (t)',ztip), err)
         if (vgd_get(vcoord ,'VIPM - level ip1 list (m)',ip1m,quiet=.true.) /= VGD_OK) err = -1
         if (vgd_get(vcoordt,'VIPT - level ip1 list (t)',ip1t,quiet=.true.) /= VGD_OK) err = -1
         out3_sfcdiag_L= .true.
         err = min(vgrid_wb_put(VGRID_M_S,vcoord, ip1m,refp0_S,refp0ls_S,F_overwrite_L=.true.), err)
         err = min(vgrid_wb_put(VGRID_T_S,vcoordt,ip1t,refp0_S,refp0ls_S,F_overwrite_L=.true.), err)
         err = vgd_free(vcoord)
         err = vgd_free(vcoordt)
         if (associated(ip1m)) deallocate(ip1m,stat=err)
         if (associated(ip1t)) deallocate(ip1t,stat=err)
         nullify(ip1m, ip1t)
      end if
      call gem_error ( err,'itf_phy_init','setting diagnostic level in vertical descriptor' )

! Determine whether winds on the lowest thermodynamic level are available for advection
      gmmk_pw_uslt_s     = 'PW_USLT'
      gmmk_pw_vslt_s     = 'PW_VSLT'
      have_slt_winds = (phy_getmeta(pmeta,gmmk_pw_uslt_s,F_npath='V',F_bpath='V') > 0 .and. &
           phy_getmeta(pmeta,gmmk_pw_vslt_s,F_npath='V',F_bpath='V') > 0)
      if (Adz_slt_winds .and. .not.have_slt_winds) then
           call gem_error ( -1,'itf_phy_init','use of surface layer winds requeseted without physics support' )
      endif
! Print table of variables requested for output

      if (Lun_out > 0) write(Lun_out,1001)
      multxmosaic = 0

      do ibus = 1,NBUS
         bus0_S = BUS_LIST_S(ibus)
         if (Lun_out > 0)  then
            write(Lun_out,1006)
            write(Lun_out,1002) bus0_S
            write(Lun_out,1006)
            write(Lun_out,1003)
         end if
         do k=1, Outp_sets
            do j=1,Outp_var_max(k)
               istat = phy_getmeta(pmeta,Outp_varnm_S(j,k), &
                       F_npath='VO',F_bpath=bus0_S(1:1),F_quiet=.true.)
               if (istat <= 0) then
                  cycle
               end if
               varname_S = pmeta%vname
               outname_S = pmeta%oname
               istat = clib_toupper(varname_S)
               istat = clib_toupper(outname_S)
               Outp_var_S(j,k) = outname_S(1:4)
               multxmosaic = max(multxmosaic,pmeta%fmul*(pmeta%mosaic+1))
               if (Lun_out > 0) write(Lun_out,1007) &
                    outname_S(1:4),varname_S(1:16),Outp_nbit(j,k), &
                    Outp_filtpass(j,k),Outp_filtcoef(j,k), &
                    Level_typ_S(Outp_lev(k))
            end do
         end do
      end do
!     maximum size of slices with one given field that is multiple+mosaic
      Outp_multxmosaic = multxmosaic+10
      if (Lun_out > 0)  write(Lun_out,1006)

      call heap_paint ()
!     ---------------------------------------------------------------
 1000 format(/,'INITIALIZATION OF PHYSICS PACKAGE (S/R itf_phy_init)', &
             /,'====================================================')
 1001 format(/'+',35('-'),'+',17('-'),'+',5('-'),'+'/'| PHYSICS VARIABLES REQUESTED FOR OUTPUT              |',5x,'|')
 1002 format('|',5X,a9,' Bus ',40x, '|')
 1003 format('|',1x,'OUTPUT',1x,'|',2x,'PHYSIC NAME ',2x,'|',2x,' BITS  |','FILTPASS|FILTCOEF| LEV |')
 1006 format('+',8('-'),'+',16('-'),'+',9('-'),'+',8('-'),'+',8('-'),'+',5('-'))
 1007 format('|',2x,a4,2x,'|',a16,'|',i5,'    |',i8,'|',f8.3,'|',a4,' |')
!     ---------------------------------------------------------------
!
      return
      end
