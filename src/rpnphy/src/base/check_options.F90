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
function check_options2() result(F_istat)
   use phy_options
   use cnv_options
   implicit none
!!!#include <arch_specific.hf>
   !@Object Check the consistency of some options
   !@Note
   !  Option validation and derived options setting are now done in
   !  phy_nml (phy_nml_check and phy_nml_post_init)
   !@returns
   integer :: F_istat
   !@Author B. Bilodeau (Feb 2007)
   !@Revisions
   ! 001      J. Milbrandt (Oct 2007) - added options for experimental and full
   !                                    versions of M-Y scheme; changed STCOND
   !                                    keys to 'my_sm', 'my_exp1', and 'my_full'
   ! 002      B. Dugas (Dec 2008)     - Add Bechtold-Kain-Fritsch schemes
   !                                    (Bechtold)
   ! 003      J. Milbrandt (Apr 2009) - Replaced 'my_exp1' with 'my_dm'
   ! 004      J. Toviessi (July 2009) - radslope will not run with oldrad
   ! 005      J. Milbrandt (Aug 2009) - Replaced 'my_full' with 'my_tm'
   ! 006      A-M. Leduc   (Mar 2010) - Fix the bechtold shallow condition
   !                                    ISHLCVT(2) = 3.
   !                                    Allow Bechtold with moistke.
   ! 007      L. Spacek    (Sep 2011)   Eliminate obsolete options
   ! 008      J. Milbrandt (Mar 2015) - Added options for MP_MY2 (new) and MP_P3
   ! 009      A. Glazer   (July 2015) - Added options for lightning diagnostics
   !*@/

#include <rmn/msg.h>
#include <rmnlib_basics.hf>

   integer, parameter :: MAXSLOFLUX = 20
   character(len=512) :: str512
   !-------------------------------------------------------------------
   F_istat = RMN_ERR

   IF_CONV_STDC: if (convec == 'SEC' .and. stcond /= 'NIL') then
      call msg(MSG_ERROR,'(check_options) stcond must be NIL when convec='// &
           trim(convec))
      return

   else if (stcond == 'CONSUN' .and. &
        .not.any(convec == (/ &
        'KFC     ', &
        'KFC2    ', &
        'BECHTOLD', &
        'NIL     '  &
        /))) then
      call msg(MSG_ERROR,'(check_options) option mismatch: stcond='// &
           trim(stcond)//' and convec='//trim(convec))
      return

   endif IF_CONV_STDC

   if (fluvert == 'MOISTKE' .and. &
        .not.any(convec == (/ &
        'KFC     ', &
        'KFC2    ', &
        'BECHTOLD', &
        'NIL     '  &
        /))) then
      call msg(MSG_ERROR,'(check_options) option mismatch: fluvert='// &
           trim(fluvert)//' and convec='//trim(convec))
      return
   endif

   if (fluvert == 'CLEF' .and. &
        .not.any(longmel == (/ &
        'TURBOUJO', &
        'BOUJO   ', &
        'BLAC62  '  &
        /))) then
      call msg(MSG_ERROR,'(check_options) option mismatch: fluvert='//&
           trim(fluvert)//' and longmel='//trim(longmel))
      return
   endif

  if (fluvert == 'CLEF' .and. &
        .not.(pbl_nonloc == 'NIL')) then
      call msg(MSG_ERROR,'(check_options) option mismatch: fluvert='//&
           trim(fluvert)//' and pbl_nonloc='//trim(pbl_nonloc))
      return
   endif

   if (NSLOFLUX > MAXSLOFLUX) then
      write(str512, '(a,i3)') &
           '(check_options) NSLOFLUX CANNOT EXCEED A VALUE OF ', MAXSLOFLUX
      call msg(MSG_ERROR, str512)
      return
   endif

   if ( PCPTYPE == 'BOURGE3D')then
      if (STCOND /= 'CONSUN'    .or. &
           (CONVEC /= 'KFC' .and. CONVEC /= 'KFC2' .and. CONVEC /= 'BECHTOLD'))   then
         call msg(MSG_ERROR,'(check_options) option mismatch: with ' // &
              'PCPTYPE=BOURGE3D you can only use STCOND=CONSUN or KFC')
         return
      endif
   endif

   if (stcond == 'MP_P3' .and. (p3_ncat < 1 .or. p3_ncat > 4)) then
      write(str512, '(a,i2,a)') '(check_options) STCOND = MP_P3 option ' // &
           'P3_NCAT=',p3_ncat,' (MUST BE BETWEEN 1 AND 4)'
      call msg(MSG_ERROR, str512)
      return
   endif

   if (RADSLOPE .and. RADIA(1:8) /= 'CCCMARAD') then
      call msg(MSG_ERROR,'(check_options) option mismatch: ' // &
           'RADSLOPE=.true. and RADIA='//trim(RADIA))
      return
   endif

   if (shal == 'BECHTOLD' .and. bkf_closures == 'EQUILIBRIUM' .and. &
        any((/bkf_entrains, bkf_detrains/) == 'BECHTOLD01')) then
      call msg(MSG_ERROR,'(check_options) Bechtold shallow ' // &
           'bkf_closures="equilibrium" closure is incompatible '// &
           'with "bechtold01" entrainment/detrainment rates '// &
           '(bkf_entrains, bkf_detrains)')
      return
   endif

   ! Check for an invalid flux enhancement modification in conres
!!$      if (ifluvert == 2) then
!!$         if (fnnmod <= 1.) then
!!$            write(options_character(1),'(e16.8)') fnnmod
!!$            if (prout) write(6,1030) 'FNNMOD',options_character(1),'> 1 ONLY'
!!$            return
!!$         endif
!!$      endif
!!$1030   format (  2(/58('*'))/'*',56x,'*'/4('****** ABORT '),'******'/ &
!!$              / '* FUNCTION CHECK_OPTIONS: ILLEGAL VALUE', &
!!$              / '*   FOR OPTION ',a,': ',a, &
!!$              / '*   ALLOWED: ',A//58('*'))

   F_istat = RMN_OK
   !-------------------------------------------------------------------
   return
end function check_options2
