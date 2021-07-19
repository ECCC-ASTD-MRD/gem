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

!/@*
subroutine phybusinit(ni,nk)
   use wb_itf_mod
   use cnv_options
   use phy_options
   use phy_status, only: phy_error_L
   use phybus
   use ens_perturb, only: ens_nc2d 
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

   include "surface.cdk"
   include "mcica.cdk"
   include "clefcon.cdk"

   character(len=6)  :: nag, ntp, nmar, wwz, nuv, isss
   integer :: ier, iverb, nsurf, i
   logical :: lcn_mpx, lcn_my2, lcn_p3i1, lcn_p3i2, lcn_p3i3, lcn_p3i4, lcn_none
   logical :: lbourg3d, lbourg
   logical :: lrslp
   logical :: lkfbe, lshal, lshbkf, lmid
   logical :: ladvtke, nadvtke, lmoistke
   logical :: lmoyhr, lmoyhrkf, lmoykfsh, lmoymid
   logical :: lgwdsm
   logical :: lccc2
   logical :: lghg, ltrigtau
   logical :: liuv
   logical :: lmoyhroz, lmoyhrgh, llinozout, llinghout, llinozage
   logical :: lcndsm
   logical :: lcons, lmoycons
   logical :: lhn_init, lsfcflx
   logical :: lsurfonly

   !---------------------------------------------------------------------

   iverb = wb_verbosity(WB_MSG_INFO)

   ier = WB_OK
   ier = min(wb_get('sfc/nsurf',   nsurf),   ier)
   ier = min(wb_get('sfc/schmsol', schmsol), ier)
   nagrege = nsurf + 1

   iverb = wb_verbosity(iverb)

   !# nagg is the dimension of aggregrated variables
   write(nag,'(a,i2)') 'A*', nagrege

   !# nipt is the number of tau/cloud top pressure bins in ISCCP histograms
   write(ntp,'(a,i2)') 'A*', ntau*nptop
   write(nuv,'(a,i2)') 'A*', RAD_NUVBRANDS

   !# nmar is the number of 2d Markov fields
   write(nmar,'(a,i2)') 'A*', ens_nc2d

   lcn_mpx  = (stcond(1:2) == 'MP')
   lcn_none = .not.lcn_mpx
   lcn_my2  = (stcond(1:6) == 'MP_MY2')
   lcn_p3i1 = (stcond == 'MP_P3' .and. p3_ncat >= 1)
   lcn_p3i2 = (stcond == 'MP_P3' .and. p3_ncat >= 2)
   lcn_p3i3 = (stcond == 'MP_P3' .and. p3_ncat >= 3)
   lcn_p3i4 = (stcond == 'MP_P3' .and. p3_ncat == 4)

   lbourg3d= (pcptype == 'BOURGE3D')
   lgwdsm  = (sgo_tdfilter > 0.)
   lrslp   = radslope
   ladvtke = advectke
   nadvtke = .not.ladvtke
   lmoyhr  = (moyhr > 0 .or. dynout)
   lkfbe   = any(convec == (/ &
        'BECHTOLD', &
        'KFC     ', &
        'KFC2    ', &
        'KFC3    '  &
        /))
   lshal   = (conv_shal /= 'NIL')
   lshbkf  = (conv_shal == 'BECHTOLD')
   lmid    = (conv_mid /= 'NIL')
   lmoyhrkf= (lmoyhr .and. lkfbe)
   lmoymid= (lmoyhr .and. lmid)
   lmoykfsh= (lmoyhr .and. lshbkf .and. bkf_lshalm)
   lbourg  = any(pcptype == (/&
        'BOURGE', &
        'NIL   '  &
        /))
   lmoistke= (fluvert == 'MOISTKE')
   lccc2   = (radia == 'CCCMARAD2')
   lghg    = (lccc2 .and. radghg_L)

   ! Compute linoz diags only on demand
   do i=1,nphyoutlist
      out_linoz = any(phyoutlist_S(i) == (/ &
           'ao3 ','ao3c','ach4','an2o','af11','af12',&
           'zch4','zn2o','zf11','zf12',&
           'ych4','yn2o','yf11','yf12',&
           'azo3','azoc','azch','azn2','azf1','azf2',&
           'ayo3','ayoc','aych','ayn2','ayf1','ayf2',&
           'ado3','ado1','ado4','ado6','ado7', &
           'adch','adn2','adf1','adf2' &
           /))
      if (out_linoz) exit
   enddo
   out_linoz = (out_linoz .or. debug_alldiag_L)
   
   llinozage = (llinoz .and. age_linoz)              ! age of air tracer off 
   llinozout = (llinoz .and. out_linoz)
   llinghout = (llingh .and. out_linoz)
   lmoyhroz =(lmoyhr .and. llinoz .and. out_linoz)
   lmoyhrgh =(lmoyhr .and. llingh .and. out_linoz)
   lcndsm  = (cond_infilter > 0.)
   ltrigtau = (kfctrigtau > 0.)
   liuv    = (any(radia == (/&
        'CCCMARAD ', &
        'CCCMARAD2'  &
        /)) .and. kntraduv_S /= '')      

   
   wwz = '1'
   lsurfonly = (fluvert == 'SURFACE')
   if (lsurfonly) wwz = '0'
   isss = '0'
   if (tofd /= 'NIL') isss = '1'

   ! Activate energy budget diagnostics only if outputs are requested by the user
   i = 1
   do while (.not.ebdiag .and. i <= nphyoutlist)
      if (any(phyoutlist_S(i) == &
           (/'t2i ', 't2im', 'tii ', 'tiim', 'q1i ', 'q1im', 'q2i ', 'q2im'/) &
           )) ebdiag = .true.
      i = i+1
   enddo
   ebdiag = (ebdiag .or. debug_alldiag_L)

   ! Activate ECMWF diagnostics only if outputs are requested by the user
   i = 1
   do while (.not.ecdiag .and. i <= nphyoutlist)
      if (any(phyoutlist_S(i) == (/'dqec','tdec','tjec','udec','vdec'/))) ecdiag = .true.
      i = i+1
   enddo
   ecdiag = (ecdiag .or. debug_alldiag_L)

   ! Activate lightning diagnostics only if outputs are requested by the user
   llight = .false.
   if (any(stcond(1:5) == (/'MP_P3 ', 'MP_MY2'/))) then
      i = 1
      do while (.not.llight .and. i <= nphyoutlist)
         if (any(phyoutlist_S(i) == (/'fdac', 'fdre'/))) llight = .true.
         i = i+1
      enddo
      llight = (llight .or. debug_alldiag_L)
   endif
   
   ! Activate refractivity diagnostics only if outputs are requested by the user
   lrefract = .false.
   i = 1
   do while (.not.lrefract .and. i <= nphyoutlist)
      if (any(phyoutlist_S(i) == (/ &
           'dcbh', 'dcnb', 'dcll', &
           'dc1m', 'dc1i', 'dcmr', &
           'dc2m', 'dc2i', 'dcst', &
           'dcth', 'dc3m', 'dc3i' &
           /))) lrefract = .true.
      i = i+1
   enddo
   lrefract = (lrefract .or. debug_alldiag_L)

   ! Activate energy budget diagnostics only if outputs are requested by the user
   lcons = .false.
   i = 1
   do while (.not.lcons .and. i <= nphyoutlist)
      if (any(phyoutlist_S(i) == (/ &
           'cec ', 'cecm', &
           'ced ', 'cedm', 'cep ', &
           'cepm', 'cem ', 'cemm', &
           'cer ', 'cerm', 'ces ', &
           'cesm', 'cqc ', &
           'cqcm', 'cqd ', 'cqdm', &
           'cqp ', 'cqpm', 'cqm ', &
           'cqmm', 'cqr ', 'cqrm', &
           'cqs ', 'cqsm', 'cey ', &
           'ceym', 'cef ', 'cefm', &
           'cqy ', 'cqym', 'cqf ', &
           'cqfm' &
           /))) lcons = .true.
      i = i+1
   enddo
   lcons = (lcons .or. debug_alldiag_L)
   lmoycons = (lcons .and. lmoyhr)
   
   ! Compute tendencies from P3 microphysics only on demand
   do i=1,nphyoutlist
      p3_comptend = any(phyoutlist_S(i) == (/ &
           'ste ','sqe ','sqce','sqre','mtqi', &
           'w7  ','w9  ','w5  ','ta  ' &
           /))
      if (p3_comptend) exit
   enddo
   p3_comptend = (p3_comptend .or. debug_alldiag_L)
   !#TODO: Check if some alloc can be avoided when p3_comptend == .false.

   etccdiag = .false.
   i = 1
   do while (.not.etccdiag .and. i <= nphyoutlist)
      if (any(phyoutlist_S(i) == (/ &
           'tccm', 'tcsh', 'tshm', 'tcsl', 'tslm', 'tcsm', 'tsmm', &
           'tczh', 'tzhm', 'tczl', 'tzlm', 'tczm', 'tzmm' &
           /))) etccdiag = .true.
      i = i+1
   enddo
   etccdiag = (etccdiag .or. debug_alldiag_L)

   
   lhn_init = (lhn /= 'NIL')
   lsfcflx = (sfcflx_filter_order > 0)

   cmt_comp_diag = .false.
   if (any(convec == (/'KFC2', 'KFC3'/)) .and. cmt_type_i /= CMT_NONE) then
      i = 1
      do while (.not.cmt_comp_diag .and. i <= nphyoutlist)
         if (any(phyoutlist_S(i) == (/ &
              'u6a ', 'u6b ', 'u6c ', 'u6am', 'u6bm', 'u6cm', &
              'v7a ', 'v7b ', 'v7c ', 'v7am', 'v7bm', 'v7cm', &
              'u6d ', 'u6dm', 'v7d ', 'v7dm' &
              /))) cmt_comp_diag = .true.
         i = i+1
      enddo
   endif
   cmt_comp_diag = (cmt_comp_diag .or. debug_alldiag_L)

#include "phyvar.hf"
   if (phy_error_L) return

   if (lbourg3d) then
      fip = fip3d
      fneige = fneige3d
   endif
   if (lcn_mpx) then
      qcmoins = qcmoinsmp
      qcphytd = qcphytdmp
      qcplus = qcplusmp
      qrphytd = qrphytdmp
   endif

   call sfc_businit(moyhr,ni,nk)
   if (phy_error_L) return

   call chm_businit(ni,nk)
   !-------------------------------------------------------------------
   return
end subroutine phybusinit


