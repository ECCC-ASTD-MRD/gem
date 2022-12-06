!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

module phy_init_mod
   use iso_c_binding
   use rpn_comm_itf_mod
   use wb_itf_mod
   use mu_jdate_mod, only: jdate_from_cmc
   use ptopo_utils, only: ptopo_grid_ipe
   use timestr_mod, only: timestr2step, timestr2sec

   use cnv_options
   use phybudget, only: pb_init
   use phy_status, only: PHY_OK, PHY_ERROR, PHY_NONE, PHY_CTRL_NML_OK, PHY_CTRL_INI_OK, phy_init_ctrl, phy_error_L
   use phy_options
   use phygridmap
   use series_mod, only: series_init
   use sfcexch_options, only: sfcexch_options3
   use ens_perturb, only: ptp_L, ptp_nc, spp_L, spp_nc, ptpenvu, ptpenvb, ptpcape, ptpcritw, &
        ptpfacreduc, ens_nc2d, ens_spp_init, ENS_OK
   private
   public :: phy_init

   interface phy_init
      module procedure phy_init_2grids
      module procedure phy_init_4grids
      module procedure phy_init_6grids
   end interface phy_init

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   include "tables.cdk"
   include "surface.cdk"

   integer, parameter :: HALO = 0
   integer, parameter :: RELAX_MODE = 0
   logical, parameter :: ALONGX = .true.
   logical, parameter :: FILL = .true.
   logical, parameter :: DO_ABORT = .true.
   logical, parameter :: PBL_FLUX_CONSISTENCY = .true.

   integer, external :: chm_init

contains

   !/@*
   function phy_init_2grids(F_path_S, F_dateo, F_dt, F_phygrid_S, &
        F_lclcore_S, F_nk, F_std_pres, F_vgrid_M_S, F_vgrid_T_S) result(F_istat)
      implicit none

      character(len=*), intent(in) :: F_path_S
      character(len=*), intent(in) :: F_phygrid_S,F_lclcore_S
      integer, intent(in) :: F_dateo
      integer, intent(in) :: F_nk
      real,    intent(in) :: F_dt, F_std_pres(F_nk)
      character(len=*), intent(in), optional :: F_vgrid_M_S, F_vgrid_T_S
      integer :: F_istat !Return status
      !*@/
      ! --------------------------------------------------------------------
      if (present(F_vgrid_M_S) .and. present(F_vgrid_T_S)) then
         F_istat = phy_init_4grids(F_path_S, F_dateo, F_dt, F_phygrid_S, &
              F_lclcore_S, 'NULL' ,'NULL', F_nk, F_std_pres, &
              F_vgrid_M_S, F_vgrid_T_S)
      else
         F_istat = phy_init_4grids(F_path_S, F_dateo, F_dt, F_phygrid_S, &
              F_lclcore_S, 'NULL' ,'NULL', F_nk, F_std_pres)
      endif
      ! --------------------------------------------------------------------
      return
   end function phy_init_2grids


   !/@*
   function phy_init_6grids(F_path_S, F_dateo, F_dt, F_phygrid_S , &
        F_lclcore_S, F_drv_glb_S, F_drv_lcl_S, &
        F_phy_glb_S, F_phy_glbcore_S, &
        F_nk, F_std_pres, F_vgrid_M_S, F_vgrid_T_S) result(F_istat)
      use hgrid_wb, only: hgrid_wb_get
      implicit none
      character(len=*), intent(in) :: F_path_S
      character(len=*), intent(in) :: F_phygrid_S, F_lclcore_S, F_drv_glb_S, &
           F_drv_lcl_S, F_phy_glb_S, F_phy_glbcore_S
      integer, intent(in) :: F_dateo
      integer, intent(in) :: F_nk
      real,    intent(in) :: F_dt, F_std_pres(F_nk)
      character(len=*), intent(in), optional :: F_vgrid_M_S, F_vgrid_T_S
      integer :: F_istat  !Return status
      !*@/
      ! --------------------------------------------------------------------
      !# Store phy global and blg core grid information
      F_istat = hgrid_wb_get(F_phy_glb_S, phy_glb_gid, &
           F_lni=phy_glb_ni, F_lnj=phy_glb_nj)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg_toall(MSG_ERROR, '(phy_init) Unable to retrieve grid info for '// &
              trim(F_phy_glb_S))
      endif
      if (RMN_IS_OK(F_istat)) &
           F_istat = hgrid_wb_get(F_phy_glbcore_S, phy_glbcore_gid)
      call collect_error(F_istat)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_ERROR, '(phy_init) Unable to retrieve grid info for '// &
              trim(F_phy_glbcore_S))
         return
      endif

      if (present(F_vgrid_M_S) .and. present(F_vgrid_T_S)) then
         F_istat = phy_init_4grids(F_path_S, F_dateo, F_dt, F_phygrid_S , &
              F_lclcore_S, F_drv_glb_S, F_drv_lcl_S, F_nk, F_std_pres, &
              F_vgrid_M_S, F_vgrid_T_S)
      else
         F_istat = phy_init_4grids(F_path_S, F_dateo, F_dt, F_phygrid_S , &
              F_lclcore_S, F_drv_glb_S, F_drv_lcl_S, F_nk, F_std_pres)
      endif
     ! --------------------------------------------------------------------
      return
   end function phy_init_6grids


   !/@*
   function phy_init_4grids(F_path_S, F_dateo, F_dt, F_phygrid_S , &
        F_lclcore_S, F_drv_glb_S, F_drv_lcl_S, &
        F_nk, F_std_pres, F_vgrid_M_S, F_vgrid_T_S) result(F_istat)
      use hgrid_wb, only: hgrid_wb_get
      implicit none
      character(len=*), intent(in) :: F_path_S
      character(len=*), intent(in) :: F_phygrid_S,F_lclcore_S,F_drv_glb_S, F_drv_lcl_S
      integer, intent(in) :: F_dateo
      integer, intent(in) :: F_nk
      real,    intent(in) :: F_dt, F_std_pres(F_nk)
      character(len=*), intent(in), optional :: F_vgrid_M_S, F_vgrid_T_S
      integer :: F_istat !Return status
      !@authors Desgagne, Chamberland, McTaggart-Cowan, Spacek -- Spring 2014
      !*@/
      integer, external :: msg_getUnit, phydebu2, sfc_init1, itf_cpl_init

      logical :: print_L
      integer :: unout, options, itype, isizeof, ntr, nsurf
      integer :: ier, p_ni, p_nj, master_pe
      integer :: type1, sizeof1, options1
      integer :: mini, maxi, lni, lnimax, li0
      integer :: minj, maxj, lnj, lnjmax, lj0
      ! --------------------------------------------------------------------
      F_istat = RMN_ERR

      call fpe_setup()

      if (phy_init_ctrl /= PHY_CTRL_NML_OK) then
         if (phy_init_ctrl == PHY_NONE .or. &
             phy_init_ctrl == PHY_CTRL_INI_OK) then
            F_istat = PHY_NONE
         else
            call msg(MSG_ERROR, '(phy_init) Must call phy_nml first!')
         endif
         return
      endif
      phy_init_ctrl = PHY_ERROR

      unout   = msg_getUnit(MSG_INFO)
      print_L = (unout > 0)

      jdateo = jdate_from_cmc(F_dateo)

      delt   = F_dt
      ier    = timestr2step(kntrad, kntrad_S, dble(delt))
      kntrad = max(1, kntrad)
      if (kntraduv_S /= '') then
         ier    = timestr2step(kntraduv, kntraduv_S, dble(delt))
         kntraduv = max(1, kntraduv)
      endif
      ier       = timestr2step(lhn_ramp, lhn_ramp_S, dble(delt))
      lhn_ramp  = max(0, lhn_ramp)
      ier       = timestr2step(lhn_start, lhn_start_S, dble(delt))
      lhn_start = max(0, lhn_start)
      ier       = timestr2step(lhn_stop, lhn_stop_S, dble(delt))
      lhn_stop  = max(0, lhn_stop)

      !# Update convective timescales to include timestep support
      if (deep_timeent_sec > 0.) ier = timestr2sec(deep_timeent_sec,deep_timeent,dble(F_dt))
      if (deep_timeconv_sec > 0.) ier = timestr2sec(deep_timeconv_sec,deep_timeconv,dble(F_dt))
      if (deep_timerefresh_sec > 0.) ier = timestr2sec(deep_timerefresh_sec,deep_timerefresh,dble(F_dt))
      if (shal_timeconv_sec > 0.) ier = timestr2sec(shal_timeconv_sec,shal_timeconv,dble(F_dt))

      nphyoutlist = 0
      ier = wb_get_meta('itf_phy/PHYOUT', type1, sizeof1, nphyoutlist, options1)
      if (.not.WB_IS_OK(ier)) nphyoutlist = 0
      allocate(phyoutlist_S(max(1,nphyoutlist)))
      phyoutlist_S(:) = ' '
      if (nphyoutlist > 0) then
         ier = wb_get('itf_phy/PHYOUT', phyoutlist_S, nphyoutlist)
         if (.not.WB_IS_OK(ier)) nphyoutlist = 0
      endif

      ier = wb_get('itf_phy/DYNOUT', dynout)
      if (.not.WB_IS_OK(ier)) dynout = .false.

      ier = wb_get('itf_phy/TLIFT', Tlift)
      if (.not.WB_IS_OK(ier)) Tlift = 0

      ier = wb_get('itf_phy/slt_winds', slt_winds)
      if (.not.WB_IS_OK(ier)) slt_winds = .false.
      
      !# Store local core grid information
      ier = hgrid_wb_get(F_lclcore_S, phy_lclcore_gid)
      if (.not.RMN_IS_OK(ier)) then
         call msg_toall(MSG_ERROR, '(phy_init) Unable to retrieve grid info for '&
              //trim(F_lclcore_S))
         return
      endif

      !# Establish physics v-grid information
      if (present(F_vgrid_M_S)) vgrid_M_S = F_vgrid_M_S
      if (present(F_vgrid_T_S)) vgrid_T_S = F_vgrid_T_S
      
      !# Establish physics h-grid information
      ier = hgrid_wb_get(F_phygrid_S, phy_lcl_gid, F_i0=phy_lcl_i0, &
           F_j0=phy_lcl_j0, F_lni=phy_lcl_ni, F_lnj=phy_lcl_nj)
      if (.not.RMN_IS_OK(ier)) then
         call msg_toall(MSG_ERROR, '(phy_init) Unable to retrieve grid info for '&
              //trim(F_phygrid_S))
         return
      endif

      phy_lcl_in = phy_lcl_i0 + phy_lcl_ni - 1
      phy_lcl_jn = phy_lcl_j0 + phy_lcl_nj - 1

      if (F_drv_glb_S /= 'NULL') then
         ier = hgrid_wb_get(F_drv_glb_S, drv_glb_gid, &
              F_lni=drv_glb_ni, F_lnj=drv_glb_nj)
         if (.not.RMN_IS_OK(ier)) then
            call msg_toall(MSG_ERROR, '(phy_init) Unable to retrieve grid info for '&
                 //trim(F_drv_glb_S))
            return
         endif
      else
         drv_glb_gid = -99
      endif

      if (F_drv_lcl_S /= 'NULL') then
         ier = hgrid_wb_get(F_drv_lcl_S, drv_lcl_gid, F_i0=drv_lcl_i0, &
              F_j0=drv_lcl_j0, F_lni=drv_lcl_ni, F_lnj=drv_lcl_nj)
         if (.not.RMN_IS_OK(ier)) then
            call msg_toall(MSG_ERROR, '(phy_init) Unable to retrieve grid info for '&
                 //trim(F_drv_lcl_S))
            return
         endif
         drv_lcl_in = drv_lcl_i0 + drv_lcl_ni - 1
         drv_lcl_jn = drv_lcl_j0 + drv_lcl_nj - 1
      else
         drv_lcl_gid = -99
      endif

      if (p_runlgt <= 0) p_runlgt = phy_lcl_ni
      p_runlgt = min(phy_lcl_ni*phy_lcl_nj,max(1,p_runlgt))
      p_ni = p_runlgt
      p_nj = phy_lcl_ni*phy_lcl_nj/p_ni
      if (p_ni*p_nj < phy_lcl_ni*phy_lcl_nj) p_nj = p_nj + 1

      phydim_ni = p_ni
      phydim_nj = p_nj

      call mapping2drivergrid()

      !# Establish physics v-grid information
      phydim_nk = F_nk

      allocate(std_p_prof(F_nk))
      std_p_prof = F_std_pres

      !# Share params with other components
      options = WB_REWRITE_NONE+WB_IS_LOCAL
      ier = WB_OK
      ier = min(wb_put('phy/atm_external', (fluvert == 'SURFACE'), options), ier)
      ier = min(wb_put('phy/jdateo'  ,jdateo   , options),ier)
      ier = min(wb_put('phy/climat'  ,climat   , options),ier)
      ier = min(wb_put('phy/sgo_tdfilter', sgo_tdfilter, options), ier)
      ier = min(wb_put('phy/lhn_filter', lhn_filter, options), ier)
      ier = min(wb_put('phy/sfcflx_filter_order', sfcflx_filter_order, options), ier)
      ier = min(wb_put('phy/sfcflx_filter_iter', sfcflx_filter_iter, options), ier)
      ier = min(wb_put('phy/convec'  ,convec   , options),ier)
      ier = min(wb_put('phy/delt'    ,delt     , options),ier)
      if (fluvert /= 'SURFACE') &
           ier = min(wb_put('phy/flux_consist', PBL_FLUX_CONSISTENCY, options), ier)
      ier = min(wb_put('phy/rad_off', (radia == 'NIL'), options),ier)
      ier = min(wb_put('phy/radslope',radslope , options),ier)
      ier = min(wb_put('phy/update_alwater', (radia == 'CCCMARAD2'), options),ier)
      ier = min(wb_put('phy/z0veg_only', (tofd /= 'NIL'), options),ier)
      ier = min(wb_put('phy/test_phy',test_phy , options),ier)
      ier = min(wb_put('phy/deep_cloudobj', deep_cloudobj, options), ier)
      ier = min(wb_put('phy/timings', timings_L, options), ier)
      ier = min(wb_put('phy/nphyoutlist', nphyoutlist, options), ier)
      if (nphyoutlist > 0) ier = min(wb_put('phy/phyoutlist', phyoutlist_S, options), ier)
      ier = min(wb_put('phy/debug_alldiag', debug_alldiag_L, options), ier)
      ier = min(wb_put('phy/input_type', input_type, options), ier)
      if (.not.WB_IS_OK(ier)) then
         call msg_toall(MSG_ERROR,'(phy_init) Problem with WB_put')
      endif
      ier = min(sfcexch_options3(), ier)
      call collect_error(ier)
      if (.not.WB_IS_OK(ier)) then
         call msg(MSG_ERROR,'(phy_init) Problem with wb_put/get for surface')
         return
      endif

      !# Stochastic physic init
      if (ens_spp_init() /= ENS_OK) then
         call msg(MSG_ERROR, '(phy_init) Problem initializing SPP perturbations')
         return
      endif
      ptp_nc = 0
      if (WB_IS_OK(ier)) &
           ier = wb_get_meta('ens/PTP', itype, isizeof, ntr, options)
      if (WB_IS_OK(ier)) then
         ier = min(wb_get('ens/PTP', ptp_L), ier)
      else
         ptp_L = .false.
         ier = WB_OK
      endif
      if (WB_IS_OK(ier) .and. ptp_L) then
         ier = min(wb_get('ens/PTP_NC'     , ptp_nc)     ,ier)
         ier = min(wb_get('ens/PTPENVU'    , ptpenvu)    ,ier)
         ier = min(wb_get('ens/PTPENVB'    , ptpenvb)    ,ier)
         ier = min(wb_get('ens/PTPCAPE'    , ptpcape)    ,ier)
         ier = min(wb_get('ens/PTPCRITW'   , ptpcritw)   ,ier)
         ier = min(wb_get('ens/PTPFACREDUC', ptpfacreduc),ier)
         if (.not.WB_IS_OK(ier)) then
            call msg_toall(MSG_ERROR,'(phydebu) Problem with wb_get #3')
         endif
      endif
      call collect_error(ier)
      if (.not.WB_IS_OK(ier)) then
         call msg(MSG_ERROR,'(phy_init) Problem initializing stochastic physics')
         return
      endif
      ens_nc2d = max(ptp_nc + spp_nc, 1)

      !# Main physic init
      ier = phydebu2(p_ni, p_nj, F_nk, F_path_S)
      if (phy_error_L) ier = RMN_ERR

      if ((moyhr>0 .or. acchr>0) .and. dynout) then
         call msg_toall(MSG_ERROR,'(phy_init) cannot use MOYHR nor ACCHR with sortie_p(average/accum)')
         ier = RMN_ERR
      endif

      call collect_error(ier)
      if (.not.RMN_IS_OK(ier)) then
         call msg(MSG_ERROR,'(phy_init) Problem in phydebu')
         return
      endif

      !# Surface init
      ier = sfc_init1()
      call collect_error(ier)
      if (.not.RMN_IS_OK(ier)) then
         call msg(MSG_ERROR,'(phy_init) Problem in sfc_debu')
         return
      endif

      !# Chemistry init
      ier = chm_init(F_path_S)
      call collect_error(ier)
      if (.not.RMN_IS_OK(ier)) then
         call msg(MSG_ERROR,'(phy_init) Problem in chm_init')
         return
      endif

      ier = WB_OK
      ier = min(wb_get('sfc/beta'        ,beta        ),ier)
      ier = min(wb_get('sfc/bh91_a'      ,bh91_a      ),ier)
      ier = min(wb_get('sfc/bh91_b'      ,bh91_b      ),ier)
      ier = min(wb_get('sfc/bh91_c'      ,bh91_c      ),ier)
      ier = min(wb_get('sfc/bh91_d'      ,bh91_d      ),ier)
      ier = min(wb_get('sfc/critlac'     ,critlac     ),ier)
      ier = min(wb_get('sfc/critmask'    ,critmask    ),ier)
      ier = min(wb_get('sfc/critsnow'    ,critsnow    ),ier)
      ier = min(wb_get('sfc/d97_as'      ,d97_as      ),ier)
      ier = min(wb_get('sfc/dg92_ci'     ,dg92_ci     ),ier)
      ier = min(wb_get('sfc/impflx'      ,impflx      ),ier)
      ier = min(wb_get('sfc/indx_soil'   ,indx_soil   ),ier)
      ier = min(wb_get('sfc/indx_glacier',indx_glacier),ier)
      ier = min(wb_get('sfc/indx_water'  ,indx_water  ),ier)
      ier = min(wb_get('sfc/indx_ice'    ,indx_ice    ),ier)
      ier = min(wb_get('sfc/indx_urb'    ,indx_urb    ),ier)
      ier = min(wb_get('sfc/indx_lake'   ,indx_lake   ),ier)
      ier = min(wb_get('sfc/indx_river'  ,indx_river  ),ier)
      ier = min(wb_get('sfc/indx_agrege' ,indx_agrege ),ier)
      ier = min(wb_get('sfc/l07_ah'      ,l07_ah      ),ier)
      ier = min(wb_get('sfc/l07_am'      ,l07_am      ),ier)
      ier = min(wb_get('sfc/leadfrac'    ,leadfrac    ),ier)
      ier = min(wb_get('sfc/n0rib'       ,n0rib       ),ier)
      ier = min(wb_get('sfc/nsurf'       ,nsurf       ),ier)
      ier = min(wb_get('sfc/sl_func_stab',sl_func_stab),ier)
      ier = min(wb_get('sfc/sl_func_unstab',sl_func_unstab),ier)
      ier = min(wb_get('sfc/tdiaglim'    ,tdiaglim    ),ier)
      ier = min(wb_get('sfc/vamin'       ,vamin       ),ier)
      ier = min(wb_get('sfc/veg_rs_mult' ,veg_rs_mult ),ier)
      ier = min(wb_get('sfc/z0dir'       ,z0dir       ),ier)
      ier = min(wb_get('sfc/zt'          ,zt          ),ier)
      ier = min(wb_get('sfc/zu'          ,zu          ),ier)
      if (.not.WB_IS_OK(ier)) then
         call msg_toall(MSG_ERROR,'(phy_init) Problem with WB_get #2')
         ier = RMN_ERR
      endif

      if (.not.WB_IS_OK(wb_get('itf_phy/zta', zta))) zta = -1.
      if (.not.WB_IS_OK(wb_get('itf_phy/zua', zua))) zua = -1.

      ier = WB_OK
      ier = min(wb_put('phy/zu',zu,options),ier)
      ier = min(wb_put('phy/zt',zt,options),ier)
      ier = min(wb_put('phy/nsurfag',nsurf+1,options),ier)

      call collect_error(ier)
      if (.not.WB_IS_OK(ier)) then
         call msg(MSG_ERROR,'(phy_init) Problem with WB_put #2')
         return
      endif

      !# Print list of physics var
      if (print_L) then
         call printbus('E')
         call printbus('D')
         call printbus('P')
         call printbus('V')
      endif

      !# Init other components
      ier = wb_get('model/outout/pe_master', master_pe)
      F_istat = series_init(phy_lcl_gid, drv_glb_gid, &
           phydim_nk, phy_lcl_ni, phy_lcl_nj, &
           moyhr, delt, jdateo, satuco, master_pe)
      if (.not.RMN_IS_OK(F_istat)) then
         call physeterror('phy_init', 'Cannot initialize time series module')
         return
      endif
      if (pb_init() /= PHY_OK) then
         call physeterror('phy_init', 'Cannot initialize physics budget module')
         return
      endif

      !# Init input I/O ezshufl distribution
      IF_INPUTIO: if (drv_glb_gid >= 0 .and. input_type == 'DIST' ) then

         ier = rpn_comm_topo_2(drv_glb_ni, mini, maxi, lni, lnimax, HALO, li0, &
              ALONGX, .not.FILL, RELAX_MODE, .not.DO_ABORT)
         if (RMN_IS_OK(ier)) then
            ier = rpn_comm_topo_2(drv_glb_nj, minj, maxj, lnj, lnjmax, HALO, lj0, &
                 .not.ALONGX, .not.FILL, RELAX_MODE, .not.DO_ABORT)
         endif
         if (.not.(RMN_IS_OK(ier) .and. &
              all((/mini, minj/) == (/1, 1/)) .and. &
              maxi >=phy_lcl_ni .and. maxj >= phy_lcl_nj)) then
            call msg_toall(MSG_ERROR, '(phy_init) Problem decomposing the domain')
            print *, ptopo_grid_ipe, '(phy_init) grd0: ', drv_glb_ni, drv_glb_nj, ":", phy_lcl_ni, phy_lcl_nj
            print *, ptopo_grid_ipe, '(phy_init) grd1: ', mini, minj, lni, lnj
            call flush(6)
            F_istat = RMN_ERR
         endif

         if (RMN_IS_OK(F_istat)) then
            phy_comm_io_id = rpn_comm_create_2dgrid(drv_glb_ni, drv_glb_nj, &
                 mini, maxi, minj, maxj)
            if (.not.RMN_IS_OK(phy_comm_io_id)) then
               call msg_toall(MSG_ERROR, '(phy_init) Unable to create an rpn_comm_2d grid')
               print *, ptopo_grid_ipe, '(phy_init) rpn_comm_create_2dgrid: ', phy_comm_io_id
               print *, ptopo_grid_ipe, '(phy_init) glb grd: ', drv_glb_ni, drv_glb_nj, ":", phy_lcl_ni, phy_lcl_nj
               print *, ptopo_grid_ipe, '(phy_init) drv grid i=', drv_lcl_ni, ":", drv_lcl_i0, drv_lcl_in
               print *, ptopo_grid_ipe, '(phy_init) drv grid j=', drv_lcl_nj, ":", drv_lcl_j0, drv_lcl_jn
               call flush(6)
               F_istat = RMN_ERR
            endif
         endif
         if (.not.RMN_IS_OK(F_istat)) then
            call msg_toall(MSG_ERROR,'(phy_init) Problem with I/O initialization')
         endif
      endif IF_INPUTIO
      call collect_error(F_istat)

      if (RMN_IS_OK(F_istat)) &
           F_istat = itf_cpl_init(F_path_S, print_L, unout, F_dateo, F_dt)

!!$      print *,'(phy_init) printbus' ; call flush(6)
!!$      call printbus('P')

      phy_init_ctrl = PHY_CTRL_INI_OK
      ! --------------------------------------------------------------------
      return
   end function phy_init_4grids

end module phy_init_mod
