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

subroutine lacs4(F_climat_L, ni, trnch)
   use tdpack_const, only: TRPL
   use sfc_options
   use sfcbus_mod
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   logical F_climat_L
   integer ni, trnch

   !@Author Bernard Bilodeau (June 2001)
   !@Object define sea ice fraction and thickness over lakes
   !@Arguments
   !          - Input -
   ! F_climat_L climate mode logical key
   ! ni         horizontal length

   include "sfcinput.cdk"

   integer i
   real twater_ini

   real, pointer, dimension(:) :: zglsea, zicedp, ziceline, zml, ztwater, &
        zvegf

#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni) => busptr(vd%NAME2%i)%ptr(:,trnch)

   MKPTR1D(zglsea,glsea)
   MKPTR1D(zicedp,icedp)
   MKPTR1D(ziceline,iceline)
   MKPTR1D(zml,ml)
   MKPTR1D(ztwater,twater)
   MKPTR1D(zvegf,vegf)

!VDIR NODEP
   DO_I: do i = 1, ni

      ! keep a copy of the original water temperature
      twater_ini = ztwater(i)

      ! Array ML is the fraction of lakes.

      IF_LAC: if (zml(i) < critlac) then  !# Not over lakes

         ! Leadfrac is the climatological value of % of leads in
         ! marine ice. The ice-covered lakes remain untouched.
         ! glsea0 contains the original value of the sea-ice
         ! analysis (updated if needed with daily increments
         ! if switch "climat" is true in subroutine climphs).

         !CC  Most current analysis already incorporate this feature.
         !CC  This approach MAY have been appropriate for older datasets,
         !CC  but it clearly deteriorates results today (BD, Feb '08).

         if (.not.F_climat_L .and. &
              (any('glsea0'==phyinread_list_s(1:phyinread_n)) .or. &
              any('vegf'==phyinread_list_s(1:phyinread_n)))) then
            zglsea(i) = zglsea(i) * (1. - leadfrac)
         endif

         ! Set minimal ice depth if zero ice thickness in
         ! analysis (or climatology) while ice fraction is non zero
         if (any('glsea0'==phyinread_list_s(1:phyinread_n)) .or. &
              any('icedp'==phyinread_list_s(1:phyinread_n)) .or. &
              any('vegf'==phyinread_list_s(1:phyinread_n))) then
            if (zglsea(i) >= 0..and.zicedp(i) <= 0.) then
               zicedp(i) = max(minicedp, 0.25 * zglsea(i))
            endif
         endif

      else !# IF_LAC - Over lakes

         ! Water temperature of ice-covered lakes is 0 C. only
         ! if point is "north" of ice line.
         if (any('iceline'==phyinread_list_s(1:phyinread_n)) .or. &
              any('glsea0'==phyinread_list_s(1:phyinread_n)) .or. &
              any('vegf'==phyinread_list_s(1:phyinread_n))) then
            if (zglsea(i) > 10.*critmask .and. ziceline(i) >= 0.5) then
               ztwater(i) = TRPL
            endif
         endif

         ! Over lakes, ice thickness is set to a minimum of 30 cm
         ! if GL > 50%. GL is then changed to 100% since the lakes
         ! are usually fully frozen. Otherwise, linear interpolation
         ! of ice thickness (varying between 0 and 30 cm) is performed
         ! when GL ranges from 0 to 50%.

         if (any('glsea0'==phyinread_list_s(1:phyinread_n)) .or. &
              any('vegf'==phyinread_list_s(1:phyinread_n))) then
            if (zglsea(i) > 0.50) then
               zicedp (i)  = max(0.30, zicedp(i))
               zglsea (i)  = 1.00 - lake_leadfrac
            else
               zicedp (i)  = 0.60 * zglsea(i)
            end if
         endif

      end if IF_LAC

      ! Restore original water temperature over salt water points if gl < 5%
      ! this fixes a bug over georgia strait in reg15 model.
      if (zvegf(i) > critmask.and.zglsea(i) < 0.05) &
           ztwater(i) = twater_ini

   end do DO_I

   return
end subroutine lacs4
