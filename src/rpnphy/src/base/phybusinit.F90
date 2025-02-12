
!/@*
subroutine phybusinit(ni,nk)
   use bus_builder, only: bb_keylist, bb_n
   use wb_itf_mod
   use cnv_options
   use phy_options
   use phy_status, only: phy_error_L, PHY_OK
   use phybusidx
   use ens_perturb, only: ens_nc2d
   use microphy_utils, only: mp_phybusinit
   use phymem, only: phymem_init, phymem_add, phymem_find, phymem_alloc
   implicit none
!!!#include <arch_specific.hf>
   !@Object Establishes requirements in terms of variables in the 4 main buses
   !        (busent, busdyn, busper and busvol) for the entire physics.
   !@Arguments
   integer, intent(in) ::  ni,nk !# horiz and vert dimensions
   !@Author M. Desgagne (Oct 1995)

   !@Revision
   ! 001      L. Spacek  (Aug 2010) - Complete rewrite
   ! 002      B. Dugas   (Oct 2010) - A few small corrections
   ! 003      L. Spacek  (Sep 2011) - Eliminate obsolete convection options
   !*@/
   
#include <rmnlib_basics.hf>
#include <rmn/msg.h>
   include "surface.cdk"

   character(len=2), parameter :: E1 = 'e1'
   character(len=2), parameter :: D1 = 'd1'
   character(len=2), parameter :: P0 = 'p0'
   character(len=2), parameter :: P1 = 'p1'
   character(len=2), parameter :: V0 = 'v0'
   character(len=2), parameter :: U0 = 'u0'
   
   character(len=4), parameter :: LVLA = 'A'
   character(len=4), parameter :: LVLA4= 'A*4'
   character(len=4), parameter :: LVLE = 'E'
   character(len=4), parameter :: LVLM = 'M'
   character(len=4), parameter :: LVLM2= 'M*2'
   character(len=4), parameter :: LVLT = 'T'

   character(len=6)  :: nag, nmar, dwwz, nuv, psss
   integer :: ier, iverb, nsurf, i
   logical :: lbourg3d, lbourg
   logical :: lkfbe, lshal, lshbkf, lmid
   logical :: lmoistke, lrpnint
   logical :: lmoyhr, lmoyhrkf, lmoykfsh, lmoymid
   logical :: lgwdsm
   logical :: lccc2
   logical :: lghg, ltrigtau
   logical :: liuv
   logical :: lmoyhroz, lmoyhrgh, llinozout, llinghout, llinozage
   logical :: lmoycons
   logical :: lhn_init, lsfcflx
   logical :: lsurfonly, lwindgust
   logical :: lpcp_frac, ladvzn, ls2, lmp
   !---------------------------------------------------------------------

   ier = phymem_init()
   if (.not.RMN_IS_OK(ier)) then
      call physeterror('phybusinit', 'Problem with phymem module init')
      return
   endif

   iverb = wb_verbosity(WB_MSG_INFO)

   ier = WB_OK
   ier = min(wb_get('sfc/nsurf',   nsurf),   ier)
   ier = min(wb_get('sfc/schmsol', schmsol), ier)
   nagrege = nsurf + 1

   iverb = wb_verbosity(iverb)

   !# nagg is the dimension of aggregrated variables
   write(nag,'(a,i2)') 'A*', nagrege

   write(nuv,'(a,i2)') 'A*', RAD_NUVBRANDS

   !# nmar is the number of 2d Markov fields
   write(nmar,'(a,i2)') 'A*', ens_nc2d

   ! Retrieve bus requirements for microphysics scheme
   if (mp_phybusinit() /= PHY_OK) then
      call physeterror('phybusinit', &
           'Cannot retrieve microphysics bus information')
      return
   endif
   
   lbourg3d= (pcptype == 'BOURGE3D')
   lgwdsm  = (sgo_tdfilter > 0.)
   lmoyhr  = (moyhr > 0 .or. dynout)
   lkfbe   = any(convec == (/ &
        'BECHTOLD', &
        'KFC     ', &
        'KFC2    ' &
        /))
   lshal   = (conv_shal /= 'NIL')
   lshbkf  = (conv_shal == 'BECHTOLD')
   lmid    = (conv_mid /= 'NIL')
   lmoyhrkf= (lmoyhr .and. lkfbe)
   lmoymid= (lmoyhr .and. lmid)
   lmoykfsh= (lmoyhr .and. lshbkf .and. bkf_lshalm)
   lbourg  = any(pcptype == (/&
        'BOURGE ', &
        'NIL    ', &
        'SPS_W19', &
        'SPS_FRC', &
        'SPS_H13'  &
        /))
   lrpnint = (fluvert == 'RPNINT')
   ladvzn  = (advectke .and. lrpnint)
   lccc2   = (radia == 'CCCMARAD2')
   lghg    = (lccc2 .and. radghg_L)
   ls2     = (stcond == 'S2')
   lmp     = (stcond(1:3) == 'MP_')

   ! Compute linoz diags only on demand
   do i=1,nphyoutlist
      out_linoz = any(phyoutlist_S(i) == (/ &
           'AO3      ', 'AO3C     ', 'ACH4     ', 'AN2O     ', 'AF11     ','AF12     ', &
           'ZCH4     ', 'ZN2O     ', 'ZF11     ', 'ZF12     ', &
           'YCH4     ', 'YN2O     ', 'YF11     ', 'YF12     ', &
           'AZO3     ', 'AZOC     ', 'AZCH     ', 'AZN2     ', 'AZF1     ','AZF2     ', &
           'AYO3     ', 'AYOC     ', 'AYCH     ', 'AYN2     ', 'AYF1     ','AYF2     ', &
           'ADO3     ', 'ADO1     ', 'ADO4     ', 'ADO6     ', 'ADO7     ', &
           'ADCH     ', 'ADN2     ', 'ADF1     ', 'ADF2     ',  &
           'O3AVG    ', 'CH4AVG   ', 'N2OAVG   ', 'F11AVG   ', 'F12AVG   ', &
           'CH4COL   ', 'N2OCOL   ', 'F11COL   ', 'F12COL   ', &
           'CH4TC    ', 'N2OTC    ', 'F11TC    ', 'F12TC    ', &
           'O3COLM   ', 'O3CCOLM  ', 'CH4COLM  ', 'N2OCOLM  ', 'F11COLM  ', 'F12COLM  ', &
           'O3TCM    ', 'O3CTCM   ', 'CH4TCM   ', 'N2OTCM   ', 'F11TCM   ', 'F12TCM   ', &
           'O3CHMTDM ', 'O1CHMTDM ', 'O4CHMTDM ', 'O6CHMTDM ', 'O7CHMTDM ', &
           'CH4CHMTDM', 'N2OCHMTDM', 'F11CHMTDM', 'F12CHMTDM' &
      /))
      if (out_linoz) exit
   enddo
   out_linoz = (out_linoz .or. debug_alldiag_L)
   
   llinozage = (llinoz .and. age_linoz)              ! age of air tracer off 
   llinozout = (llinoz .and. out_linoz)
   llinghout = (llingh .and. out_linoz)
   lmoyhroz =(lmoyhr .and. llinoz .and. out_linoz)
   lmoyhrgh =(lmoyhr .and. llingh .and. out_linoz)
   ltrigtau = (kfctrigtau > 0.)
   liuv    = (any(radia == (/&
        'CCCMARAD ', &
        'CCCMARAD2'  &
        /)) .and. kntraduv_S /= '')      

   
   dwwz = 'd1'
   lsurfonly = (fluvert == 'SURFACE')
   if (lsurfonly) dwwz = 'd0'
   psss = 'p0'
   if (tofd /= 'NIL') psss = 'p1'
   lpcp_frac = lsurfonly .and. (pcptype == 'SPS_FRC')


   ! Activate energy budget diagnostics only if outputs are requested by the user
   i = 1
   do while (.not.ebdiag .and. i <= nphyoutlist)
      if (any(phyoutlist_S(i) == &
           (/ &
           'T2I   ', 'T2IM  ', 'TII   ', 'TIIM  ', &
           'Q1I   ', 'Q1IM  ', 'Q2I   ', 'Q2IM  ', &
           'Q1APP ', 'Q1APPM', 'Q2APP ', 'Q2APPM'  &
           /) &
           )) ebdiag = .true.
      i = i+1
   enddo
   ebdiag = (ebdiag .or. debug_alldiag_L)

   ! Activate ECMWF diagnostics only if outputs are requested by the user
   i = 1
   do while (.not.ecdiag .and. i <= nphyoutlist)
      if (any(phyoutlist_S(i) == (/ &
           'DQEC    ', 'TDEC    ', 'TJEC    ', 'UDEC    ', 'VDEC    ', &
           'QDIAGEC ', 'TDDIAGEC', 'TDIAGEC ', 'UDIAGEC ', 'VDIAGEC '  &
           /))) ecdiag = .true.
      i = i+1
   enddo
   ecdiag = (ecdiag .or. debug_alldiag_L)

   ! Activate lightning diagnostics only if outputs are requested by the user
   llight = .false.
   if (any(bb_keylist(:bb_n) == "LIGHTNING")) then
      i = 1
      do while (.not.llight .and. i <= nphyoutlist)
         if (any(phyoutlist_S(i) == (/ &
              'FDAC   ', 'FDRE   ', &
              'AFOUDRE', 'FOUDRE '  &
              /))) llight = .true.
         i = i+1
      enddo
      llight = (llight .or. debug_alldiag_L)
   endif
   
   ! Activate refractivity diagnostics only if outputs are requested by the user
   lrefract = .false.
   i = 1
   do while (.not.lrefract .and. i <= nphyoutlist)
      if (any(phyoutlist_S(i) == (/ &
           'DCBH      ', 'DCNB      ', 'DCLL      ', &
           'DC1M      ', 'DC1I      ', 'DCMR      ', &
           'DC2M      ', 'DC2I      ', 'DCST      ', &
           'DCTH      ', 'DC3M      ', 'DC3I      ', &
           'DCT_BH    ', 'DCT_NB    ', 'DCT_LVL   ', &
           'DCT_LVLMAX', 'DCT_LVLMIN', 'DCT_MOREF ', &
           'DCT_SNDMAX', 'DCT_SNDMIN', 'DCT_STR   ', &
           'DCT_THICK ', 'DCT_TRDMAX', 'DCT_TRDMIN'  &
           /))) lrefract = .true.
      i = i+1
   enddo
   lrefract = (lrefract .or. debug_alldiag_L)

   ! Activate wind gust estimate only if outputs are requested by the user
   lwindgust = .false.
   i = 1
   do while (.not.lwindgust .and. i <= nphyoutlist)
      if (any(phyoutlist_S(i) == (/ &
           'WGE   ', 'WGX   ', 'WGN   ', &
           'SDWD  ', 'SDWS  ', &
           'WGMAX ', 'WGMIN ', &
           'SDTSWD', 'SDTSWS'  &
           /))) lwindgust = .true.
      i = i+1
   enddo
   lwindgust = (lwindgust .or. debug_alldiag_L)
   
   ! Activate energy budget diagnostics only if outputs are requested by the user
   lcons = .false.
   i = 1
   do while (.not.lcons .and. i <= nphyoutlist)
      if (any(phyoutlist_S(i) == (/ &
           'CEC     ', 'CECM    ', &
           'CED     ', 'CEDM    ', 'CEP     ', &
           'CEPM    ', 'CEM     ', 'CEMM    ', &
           'CER     ', 'CERM    ', 'CES     ', &
           'CESM    ', 'CQC     ', &
           'CQCM    ', 'CQD     ', 'CQDM    ', &
           'CQP     ', 'CQPM    ', 'CQM     ', &
           'CQMM    ', 'CQR     ', 'CQRM    ', &
           'CQS     ', 'CQSM    ', 'CEY     ', &
           'CEYM    ', 'CEF     ', 'CEFM    ', &
           'CQY     ', 'CQYM    ', 'CQF     ', &
           'CQFM    ', &
           'CONECND ', 'CONECNDM', &
           'CONEDC  ', 'CONEDCM ', 'CONEPBL ', &
           'CONEPBLM', 'CONEMO  ', 'CONEMOM ', &
           'CONERAD ', 'CONERADM', &
           'CONESC  ', 'CONESCM ', 'CONQCND ', &
           'CONQCNDM', 'CONQDC  ', 'CONQDCM ', &
           'CONQPBL ', 'CONQPBLM', 'CONQMO  ', &
           'CONQMOM ', 'CONQRAD ', 'CONQRADM', &
           'CONQSC  ', 'CONQSCM ', 'CONEDYN ', &
           'CONEDYNM', 'CONEPHY ', 'CONEPHYM', &
           'CONQDYN ', 'CONQDYNM', 'CONQPHY ', &
           'CONQPHYM'  &
           /))) lcons = .true.
      i = i+1
   enddo
   lcons = (lcons .or. debug_alldiag_L)
   lmoycons = (lcons .and. lmoyhr)

   etccdiag = .false.
   i = 1
   do while (.not.etccdiag .and. i <= nphyoutlist)
      if (any(phyoutlist_S(i) == (/ &
           'TCCM', 'TCSH', 'TSHM', 'TCSL', 'TSLM', 'TCSM', 'TSMM', &
           'TCZH', 'TZHM', 'TCZL', 'TZLM', 'TCZM', 'TZMM' &
           /))) etccdiag = .true.
      i = i+1
   enddo
   etccdiag = (etccdiag .or. debug_alldiag_L)
   
   etccdiagout = .false.
   i = 1
   do while (.not.etccdiagout .and. i <= nphyoutlist)
      if (any(phyoutlist_S(i) == (/ &
           'TCCM', 'NF  ', 'TSHM', 'TSMM', 'TSLM', 'TZHM', 'TZMM', 'TZLM', &
           'NTAF' &
           /))) etccdiagout = .true.
      i = i+1
   enddo
   etccdiagout = (etccdiagout .or. debug_alldiag_L)

   
   lhn_init = (lhn /= 'NIL')
   lsfcflx = (sfcflx_filter_order > 0)

   cmt_comp_diag = .false.
   if (any(convec == (/'KFC2', 'KFC3'/)) .and. cmt_type_i /= CMT_NONE) then
      i = 1
      do while (.not.cmt_comp_diag .and. i <= nphyoutlist)
         if (any(phyoutlist_S(i) == (/ &
              'U6A   ', 'U6B   ', 'U6C   ', 'U6AM  ', 'U6BM  ', 'U6CM  ', &
              'V7A   ', 'V7B   ', 'V7C   ', 'V7AM  ', 'V7BM  ', 'V7CM  ', &
              'U6D   ', 'U6DM  ', 'V7D   ', 'V7DM  ', &
              'UFCP1 ', 'UFCP2 ', 'UFCP3 ', 'UFCP1M', 'UFCP2M', 'UFCP3M', &
              'VFCP1 ', 'VFCP2 ', 'VFCP3 ', 'VFCP1M', 'VFCP2M', 'VFCP3M', &
              'SUFCP ', 'SUFCPM', 'SVFCP ', 'SVFCPM'  &
              /))) cmt_comp_diag = .true.
         i = i+1
      enddo
   endif
   cmt_comp_diag = (cmt_comp_diag .or. debug_alldiag_L)

#include "phymkptr.hf"
#include "phyvar.hf"
   if (phy_error_L) return

   call sfc_businit(moyhr,ni,nk)
   if (phy_error_L) return

#ifdef HAVE_MACH
   call chm_businit(ni,nk)
#endif

   !#NOTE: phymem_alloc must be done before any call to phymem_find
   if (debug_alldiag_L) then
      ier = phymem_alloc(debug_mem_L, (/'*'/))
   else
      ier = phymem_alloc(debug_mem_L, phyoutlist_S)
   endif
   if (.not.RMN_IS_OK(ier)) &
        call physeterror('phybusinit', 'problem in phymem_alloc')
   if (phy_error_L .or. .not.RMN_IS_OK(ier)) return

#undef PHYMKPTR
#define PHYPTRGETIDX
#include "phymkptr.hf"
#include "phyvar.hf"
   if (phy_error_L) return   

   sigw = sigt

   if (lbourg3d) then
      fip = fip3d
      fneige = fneige3d
   endif

   if (qcmoinsmp > 0) qcmoins = qcmoinsmp
   if (qcphytdmp > 0) qcphytd = qcphytdmp
   if (qcplusmp > 0) qcplus = qcplusmp
   if (qrphytdmp > 0) qrphytd = qrphytdmp
   
   !-------------------------------------------------------------------
   return
end subroutine phybusinit

