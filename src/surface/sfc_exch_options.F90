!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------

module sfcexch_options
   private
   public :: sfcexch_options3

contains

   !/@*
   function sfcexch_options3() result(F_istat)
      use wb_itf_mod
      use sfc_options
      use sfcbus_mod
      use timestr_mod, only: timestr2step
      implicit none
!!!#include <arch_specific.hf>
      integer :: F_istat
      !@Object initialization of the surface parameters at the beginning
      !        of each execution of the model
      !@Author  L. Spacek (Spring 2013)
      !@Revisions
      ! 001 K. Winger/M. Mackay    (Feb 2017/Sep 2022) - Add 'indx_lake' and 'indx_river' (M.A.)
      !*@/

#include <rmn/msg.h>
#include <rmnlib_basics.hf>
      include "isbapar.cdk"
      include "tebcst.cdk"

      integer :: ier, options, iverb, i
      logical :: debug_alldiag_L, tmp_atm_tplus
      !---------------------------------------------------------------------
      F_istat = RMN_ERR

      options = WB_REWRITE_NONE+WB_IS_LOCAL

      iverb = wb_verbosity(WB_MSG_FATAL)
      ier = wb_get('phy/flux_consist',tmp_atm_tplus)
      if (RMN_IS_OK(ier)) atm_tplus = tmp_atm_tplus
      ier = wb_verbosity(WB_MSG_INFO)
      
      ier = WB_OK

      ier = min(wb_get('phy/jdateo',jdateo),ier)
      ier = min(wb_get('phy/climat',climat),ier)
      ier = min(wb_get('phy/delt',delt),ier)
      ier = min(wb_get('phy/rad_off',rad_off),ier)
      ier = min(wb_get('phy/radslope',radslope),ier)
      ier = min(wb_get('phy/atm_external',atm_external),ier)
      ier = min(wb_get('phy/timings',timings_L),ier)
      ier = min(wb_get('phy/update_alwater',update_alwater),ier)
      ier = min(wb_get('phy/z0veg_only',z0veg_only),ier)
      ier = min(wb_get('phy/nphyoutlist',nphyoutlist),ier)
      if (nphyoutlist > 0) then
         allocate(phyoutlist_S(nphyoutlist))
         ier = min(wb_get('phy/phyoutlist',phyoutlist_S,nphyoutlist),ier)
      endif
      ier = min(wb_get('phy/debug_alldiag',debug_alldiag_L),ier)
            
      ier = min(wb_put('sfc/beta',beta,options),ier)
      ier = min(wb_put('sfc/bh91_a',bh91_a,options),ier)
      ier = min(wb_put('sfc/bh91_b',bh91_b,options),ier)
      ier = min(wb_put('sfc/bh91_c',bh91_c,options),ier)
      ier = min(wb_put('sfc/bh91_d',bh91_d,options),ier)
      ier = min(wb_put('sfc/critlac',critlac,options),ier)
      ier = min(wb_put('sfc/critmask',critmask,options),ier)
      ier = min(wb_put('sfc/critsnow',critsnow,options),ier)
      ier = min(wb_put('sfc/d97_as',d97_as,options),ier)
      ier = min(wb_put('sfc/dg92_ci',dg92_ci,options),ier)
      ier = min(wb_put('sfc/impflx',impflx,options),ier)
      ier = min(wb_put('sfc/indx_soil',indx_soil,options),ier)
      ier = min(wb_put('sfc/indx_glacier',indx_glacier,options),ier)
      ier = min(wb_put('sfc/indx_water',indx_water,options),ier)
      ier = min(wb_put('sfc/indx_ice',indx_ice,options),ier)
      ier = min(wb_put('sfc/indx_urb',indx_urb,options),ier)
      ier = min(wb_put('sfc/indx_lake',indx_lake,options),ier)
      ier = min(wb_put('sfc/indx_river',indx_river,options),ier)
      ier = min(wb_put('sfc/indx_agrege',indx_agrege,options),ier)
      ier = min(wb_put('sfc/l07_ah',l07_ah,options),ier)
      ier = min(wb_put('sfc/l07_am',l07_am,options),ier)
      ier = min(wb_put('sfc/leadfrac',leadfrac,options),ier)
      ier = min(wb_put('sfc/n0rib',n0rib,options),ier)
      ier = min(wb_put('sfc/sl_Lmin_soil',sl_Lmin_soil,options),ier)
      ier = min(wb_put('sfc/sl_Lmin_glacier',sl_Lmin_glacier,options),ier)
      ier = min(wb_put('sfc/sl_Lmin_water',sl_Lmin_water,options),ier)
      ier = min(wb_put('sfc/sl_Lmin_seaice',sl_Lmin_seaice,options),ier)
      ier = min(wb_put('sfc/sl_Lmin_town',sl_Lmin_town,options),ier)
      ier = min(wb_put('sfc/sl_func_stab',sl_func_stab,options),ier)
      ier = min(wb_put('sfc/sl_func_unstab',sl_func_unstab,options),ier)
      ier = min(wb_put('sfc/tdiaglim',tdiaglim,options),ier)
      ier = min(wb_put('sfc/vamin',vamin,options),ier)
      ier = min(wb_put('sfc/veg_rs_mult',veg_rs_mult,options),ier)
      ier = min(wb_put('sfc/z0dir',z0dir,options),ier)
      ier = min(wb_put('sfc/zt',zt,options),ier)
      ier = min(wb_put('sfc/zu',zu,options),ier)

      iverb = wb_verbosity(iverb)

      if (.not.RMN_IS_OK(ier)) then
         call msg(MSG_ERROR,'(sfc_exch_options) probleme in wb_put/get')
         return
      endif

      if (kntveg_S /= '') then
         ier = timestr2step(kntveg, kntveg_S, dble(delt))
         if (.not.RMN_IS_OK(ier)) then
            call msg(MSG_ERROR,'(sfc_exch_options) Problem converting kntveg_S='//trim(kntveg_S))
            return
         endif
      endif

      IF_PHYOUTLIST: if (associated(phyoutlist_S)) then
         ! Trigger calculation of on-demand diagnostics. Modular computation to save comp. time

         ! List for UTCI (algorithm costs a bit) (4)
         i = 1
         do while (i <= nphyoutlist .and. .not.thermal_stress_utci)
            if (any(phyoutlist_S(i) == (/'DXSU','DXHD', &
                 'DXFS','DXFD'/))) then
               thermal_stress_utci = .true.
               thermal_stress = .true.
            endif
            i = i+1
         enddo

         ! List for var in shade
         if (associated(phyoutlist_S)) then
            i = 1
            do while (i <= nphyoutlist .and. .not.thermal_stress_shade)
               if (any(phyoutlist_S(i) == (/'GXHD','RTHD','GTHD','DXHD','DXFD'/))) then
                  thermal_stress_shade = .true.
                  thermal_stress = .true.
               endif
               i = i+1
            enddo
         endif

         ! List for var over roof (17)
         if (associated(phyoutlist_S)) then
            i = 1
            do while (i <= nphyoutlist .and. .not.thermal_stress_roof)
               if (any(phyoutlist_S(i) == (/'DXFS','DXFD', &
                    'GXFS','GXFD','RTFS','RTFD','GTFS','GTFD',  &
                    'QSRF','QSFK','QLRF','QLFK','QSFW','QLFW',  &
                    'DCFS','DCFD','WBRF'  &
                    /))) then
                  thermal_stress_roof = .true.
                  thermal_stress = .true.
               endif
               i = i+1
            enddo
         endif

         ! List for main var
         if (associated(phyoutlist_S)) then
            i = 1
            do while (i <= nphyoutlist .and. .not.thermal_stress)
               if (any(phyoutlist_S(i) == (/'GXSU','RTSU', &
                    'GTSU','QSSU','QSSK','QLSK',    &
                    'QSWL','QLWL','QSRD','QLRD'/))) thermal_stress = .true.
               i = i+1
            enddo
         endif
         
      else !IF_PHYOUTLIST
         
         iverb = wb_verbosity(WB_MSG_FATAL)
         if (.not.RMN_IS_OK(wb_get('sps_cfgs/ThStress_L',       thermal_stress)))       thermal_stress       = .false.
         if (.not.RMN_IS_OK(wb_get('sps_cfgs/ThStress_utci_L',  thermal_stress_utci)))  thermal_stress_utci  = .false.
         if (.not.RMN_IS_OK(wb_get('sps_cfgs/ThStress_shade_L', thermal_stress_shade))) thermal_stress_shade = .false.
         if (.not.RMN_IS_OK(wb_get('sps_cfgs/ThStress_roof_L',  thermal_stress_roof)))  thermal_stress_roof  = .false.
         iverb = wb_verbosity(iverb)

      endif IF_PHYOUTLIST
      
      thermal_stress_utci  = (thermal_stress_utci  .or. debug_alldiag_L)
      thermal_stress_shade = (thermal_stress_shade .or. debug_alldiag_L)
      thermal_stress_roof  = (thermal_stress_roof  .or. debug_alldiag_L)
      thermal_stress       = (thermal_stress       .or. debug_alldiag_L .or. &
           thermal_stress_utci .or. thermal_stress_shade .or. thermal_stress_roof)

      F_istat = RMN_OK
      !----------------------------------------------------------------------
      return
   end function sfcexch_options3

end module sfcexch_options
