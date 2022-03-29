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

subroutine ccc1_raylev(taur, ig, dp, rmu3, il1, il2, ilg, lay)
   implicit none
!!!#include <arch_specific.hf>

   integer ilg, lay, il1, il2, k, i, ig
   real ri0(6), ri2(3), taur(ilg,lay), dp(ilg,lay), rmu3(ilg)

   !Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !        JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)

   !Object
   !        Rayleigh scattering for each sub-band in bands1, visible region
   !        taur is the optical depth rayleigh scattering for a given layer
   !        for uvc (35700 - 50000 cm^-1), since the optical depth of o3 and
   !        o2 are very large, rayleigh scattering effect is neglected, it
   !        is shown even for 10% o3 amount of the standard atmo, the
   !        rayleigh scattering for uvc still can be neglected.
   !        For par and uva, since their spectral ranges are very wide, small
   !        errors could occur for large zenith angle, slightly adjustment
   !        is needed, this does mean the rayleigh optical depth is related
   !        solar zenith angle for multiple scattering process in swtran.

   !Arguments
   ! taur   rayleigh optical depth
   ! dp     air mass path for a model layer (explained in raddriv).
   ! rmu3   a factor of solar zenith angle, given in raddriv

   data ri0 / .75996e-04, .20144e-03, .44813e-03, .80206e-03, &
        .10430e-02, .12756e-02 /
   data ri2 / .38900e-05, .11500e-04, .19100e-04 /

   if (ig .le. 3) then
      do k = 1, lay
         do i = il1, il2
            taur(i,k) = (ri0(ig) - ri2(ig) * rmu3(i)) * dp(i,k)
         enddo
      enddo
   else
      do k = 1, lay
         do i = il1, il2
            taur(i,k) =  ri0(ig) * dp(i,k)
         enddo
      enddo
   endif

   return
end subroutine ccc1_raylev
