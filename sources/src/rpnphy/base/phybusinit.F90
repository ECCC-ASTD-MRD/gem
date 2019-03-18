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
   use cnv_options
   use phy_options
   use phy_status, only: phy_error_L 
   use phybus
   implicit none
#include <arch_specific.hf>
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

#include <WhiteBoard.hf>
   include "ens.cdk"
   include "surface.cdk"
   include "mcica.cdk"
   include "clefcon.cdk"
   include "buses.cdk"

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
   logical :: llight, lcndsm
   logical :: lcons, lmoycons

   !---------------------------------------------------------------------

   buslck = .false.

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
   write(nmar,'(a,i2)') 'A*', imrkv2

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
   llight  = lightning_diag
   lcndsm  = (cond_infilter > 0.)
   ltrigtau = (kfctrigtau > 0.)

   wwz = '1'
   if (offline) wwz = '0'
   isss = '0'
   if (tofd /= 'NIL') isss = '1'


   ! Activate energy budget diagnostics only if outputs are requested by the user
   lcons = .false.
   i = 1
   do while (.not.lcons .and. i <= nphyoutlist)
      if (any(phyoutlist_S(i) == (/ &
           'clse', 'cec ', 'cecm', &
           'ced ', 'cedm', 'cep ', &
           'cepm', 'cem ', 'cemm', &
           'cer ', 'cerm', 'ces ', &
           'cesm', 'cqt ', 'cqc ', &
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
   lmoycons = (lcons .and. lmoyhr)

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


