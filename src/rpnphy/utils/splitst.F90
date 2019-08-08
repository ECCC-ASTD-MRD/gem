!-------------------------------------- LICENCE BEGIN --------------------------
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
!-------------------------------------- LICENCE END ----------------------------

module splitst
   implicit none
   private
   public :: splitst4

contains

   function splitst4(cvn, con, cin, csn, cvd1, cvd2, cvs, fmosaik, fmul,  &
        cvb, dynini, stagg, F_vmin, F_vmax, F_wload, F_hzd, F_monot, F_massc, &
        F_string_S) result(F_istat)
      use str_mod
      use tracers_attributes_mod, only: tracers_attributes
      use phy_options
      implicit none
!!!#include <arch_specific.hf>
      character(len=*), intent(in)  :: F_string_S
      character(len=*), intent(out) :: con,cvn,cin,csn,cvd1,cvd2,cvb,cvs
      integer, intent(out) ::  fmul,fmosaik,dynini,stagg
      real, intent(out) :: F_vmin, F_vmax
      integer, intent(out) :: F_wload, F_hzd, F_monot, F_massc
      integer :: F_istat
      !
      !Author
      !          M. Desgagne (Oct 1995)
      !
      !Revision
      ! 001      B. Bilodeau (Sept 1996) - Add 2-letter names
      ! 002      B. Bilodeau (Aug  1998) - Add staggered levels
      ! 003      B. Bilodeau (Jun  2005) - Add mosaic capability for CLASS
      !                                    and remove fadd
      ! 004      V. Lee (Mar 2011) - fmosaik = 0 if no mosaic tiles.
      !
      !Object
      !
      !Arguments
      !            - Output -
      ! cvn       formal name (VN)
      ! con       output name (ON)
      ! cin       input name (IN)
      ! cvd1      formal description (VD)
      ! cvd2      complete shape (VS)
      ! cvs       shape --A, M, T, E-- (VS)
      ! fmosaik   mosaic factor (number of types of soil surfaces for CLASS)
      ! fmul      multiplicative factor
      ! cvb       bus identification (VB)
      ! dynini    flag for initialysation by the dynamics (1=yes)
      ! stagg     flag for staggered levels (0=non staggered; 1=staggered)
      !
      !            - Input -
      ! string    input description string including all tokens (IN is optional)

#include <rmnlib_basics.hf>
#include <msg.h>
      include "buses.cdk"

      integer, external :: str_split2keyval

      integer,parameter :: IDX_VN = 1
      integer,parameter :: IDX_VS = 2
      integer,parameter :: IDX_VB = 3
      integer,parameter :: IDX_VD = 4
      integer,parameter :: IDX_ON = 5
      integer,parameter :: IDX_IN = 6
      integer,parameter :: IDX_SN = 7
      integer,parameter :: IDX_MIN = 8
      integer,parameter :: IDX_MAX = 9
      integer,parameter :: IDX_WLOAD = 10
      integer,parameter :: IDX_HZD = 11
      integer,parameter :: IDX_MONOT = 12
      integer,parameter :: IDX_MASSC = 13

      integer,parameter :: NIDX = 13
      integer,parameter :: NIDX_MIN = 5
      integer,parameter :: NIDX_EXTRA = 32 - NIDX
      integer,parameter :: NIDX_MAX = NIDX + NIDX_EXTRA
      integer,parameter :: KEY = 1
      integer,parameter :: VAL = 2

      character(len=32),parameter :: KNOWN_KEYS(NIDX) = (/&
           'vn     ', & !# mandatory
           'vs     ', & !# mandatory
           'vb     ', & !# mandatory
           'vd     ', & !# mandatory
           'on     ', & !# mandatory
           'in     ', &
           'sn     ', &
           'min    ', &
           'max    ', &
           'wload  ', &
           'hzd    ', &
           'monot  ', &
           'massc  '  &
           /)

      character(len=1024) :: string_S,kv_S(2,NIDX_MAX),s1_S,s2_S,attr_S
      integer :: nkeys, istat, ind,iwload, ihzd, imonot, imassc
      real :: rvmin
      !-------------------------------------------------------------------
      F_istat = RMN_ERR

      string_S = F_string_S
      call str_tab2space(string_S)
      string_S = adjustl(string_S)

      kv_S = ' '
      kv_S(KEY,1:NIDX) = KNOWN_KEYS(1:NIDX)
      nkeys = str_split2keyval(kv_S,string_S,NIDX_MAX)
      if (nkeys < NIDX_MIN .or. any(kv_S(VAL,1:NIDX_MIN) == ' ')) then
         call msg(MSG_ERROR,'(gesdict) Mandatory params not all provided for: '//trim(string_S))
         return
      endif

      cvn  = kv_S(VAL,IDX_VN)
      attr_S = ' '
      ind = index(cvn,",")
      if (ind /= 0) then
         attr_S = cvn(ind+1: )
         cvn  = cvn(1:ind-1)
      endif
      istat = tracers_attributes(attr_S, iwload, ihzd, imonot, imassc, rvmin, F_ignore_L=.true.)
      !#TODO: should we allow use of this syntax anymore? if yes, apply these... warning, default values can be changed in tracers_attributes

      con  = kv_S(VAL,IDX_ON)
      cin = con
      if (kv_S(VAL,IDX_IN) /= '') cin  = kv_S(VAL,IDX_IN)
      csn = con
      if (kv_S(VAL,IDX_SN) /= '') csn  = kv_S(VAL,IDX_SN)
      cvd1 = kv_S(VAL,IDX_VD)
      cvd2 = kv_S(VAL,IDX_VS)
      cvs  = kv_S(VAL,IDX_VS)(1:1)
      cvb  = kv_S(VAL,IDX_VB)(1:1)

      select case(cvs)
      case("A")  !# arbitrary
         stagg = BUSPAR_STAG_SFC
      case("M")  !# momentum
         stagg = BUSPAR_STAG_MOM
      case("T")  !# thermo
         stagg = BUSPAR_STAG_THERMO
      case("E")  !# energy
         stagg = BUSPAR_STAG_ENERGY
      case default
         call msg(MSG_ERROR,'(gesdict) VS=(SHAPE) NOT ALLOWED: '//trim(string_S))
         return
      end select

      dynini = 0
      if (kv_S(VAL,IDX_VB)(2:2) == '1') dynini = 1

      if (.not.any(cvb == (/'E','D','P','V'/)))  then
         call msg(MSG_ERROR,'(gesdict) VB=(BUS) NOT ALLOWED: '//trim(string_S))
         return
      endif

      string_S = kv_S(VAL,IDX_VS)(2:len_trim(kv_S(VAL,IDX_VS)))
      call str_split(s1_S,s2_S,string_S,'@')
      fmosaik = 0
      if (s2_S /= '') then
         istat = str_toint(fmosaik, s2_S)
         if (.not.RMN_IS_OK(istat)) return
      endif

      string_S =  s1_S
      call str_split(s1_S,s2_S,string_S,'*')
      fmul = 1
      if (s2_S /= '') then
         istat = str_toint(fmul, s2_S)
         if (.not.RMN_IS_OK(istat)) return
      endif

      F_vmin = -1.*huge(F_vmin)
      if (kv_S(VAL,IDX_MIN) /= '') then
         istat = str_toreal(F_vmin, kv_S(VAL,IDX_MIN))
         if (.not.RMN_IS_OK(istat)) return
      endif

      F_vmax = huge(F_vmax)
      if (kv_S(VAL,IDX_MAX) /= '') then
         istat = str_toreal(F_vmax, kv_S(VAL,IDX_MAX))
         if (.not.RMN_IS_OK(istat)) return
      endif

      F_wload = 0
      if (kv_S(VAL,IDX_WLOAD) /= '') then
         istat = str_toint(F_wload, kv_S(VAL,IDX_WLOAD))
         if (.not.RMN_IS_OK(istat)) return
      endif

      F_hzd = 0
      if (kv_S(VAL,IDX_HZD) /= '') then
         istat = str_toint(F_hzd, kv_S(VAL,IDX_HZD))
         if (.not.RMN_IS_OK(istat)) return
      endif

      F_monot = 1
      if (kv_S(VAL,IDX_MONOT) /= '') then
         istat = str_toint(F_monot, kv_S(VAL,IDX_MONOT))
         if (.not.RMN_IS_OK(istat)) return
      endif

      F_massc = 0
      if (kv_S(VAL,IDX_MASSC) /= '') then
         istat = str_toint(F_massc, kv_S(VAL,IDX_MASSC))
         if (.not.RMN_IS_OK(istat)) return
      endif

      F_istat = RMN_OK
      !-------------------------------------------------------------------
      return
   end function splitst4

end module splitst
