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

subroutine compz0_1(optz0, z0, z0loc, z0t, fm, va, fcor, n)
   use tdpack_const, only: GRAV, KARMAN
   use sfc_options
   implicit none
!!!#include <arch_specific.hf>
   integer n, optz0
   real va(n),fm(n)
   real z0(n),z0t(n), z0loc(n)
   real fcor(n)

   !@Author  S.Leroyer (Apr 2011)
   !@Revision
   !  001     M. Abrahamowicz (May 2013) - Shift options by 1., so that option 0. is a RETURN
   !  002     M.A (June 2020) - Add new input: z0loc for option 9
   !@Object compute aerodynamic and scalar roughness lenghts
   ! 
   !          options :
   ! 0. RETURN        : no calculations done
   ! 1. sea ice       : z0 fixed and z0t=z0   (SEAICE1 classic)
   ! 2. land surface  : z0 fixed (related to surface covers) and z0t/z0 fixed ratio  (ISBA / FORCE-RESTORE classic)
   ! 3. water surface : z0 from Charnock (1955) and z0t dependant on latitude (WATER classic)
   ! 4. glaciers      : z0 from geophy and z0t classic (GLACIERS1 classic)
   ! 5. water surface : Deacu et al. ...
   ! 6. water surface : Fairall et al. 2003
   ! 7. urban surface : z0 from geophy and z0t from Brutsaert (1982)
   ! 8. urban surface : z0 from geophy and z0t from Kanda (2007)
   ! 9. land surface  : LOCAL z0 from geophy and z0t from Zilitinkevich (1995) 
   ! 10.urban surface : z0 fixed (related to surface covers) and z0t/z0 fixed ratio=200 (Mascart et al. 1995) 
   !                 
   !@Arguments
   !          - Output -
   ! z0       aerodynamic roughness length
   ! z0loc    local aerodynamic roughness length (no orography)
   ! z0t      scalar roughness lenght for heat
   !          - Input -
   ! fm       momentum stability function
   ! va       wind speed at ZU
   ! n        horizontal dimension
   ! optz0    option for z0 formulations
   ! fcor     Coriolis factor
   !@Variables
   !          - Local -
   ! ue       friction velocity
   ! trrn     turbulent roughness Reynolds number
   ! lat      latitude (computed with fcor)
   !          - Constants -
   ! VISCOSITY air viscosity
   ! Z0MAX     upper limit for z0 over water
   ! Z0HCON    fixed value for deep convection zone
   ! ac        Charnock's formulae constant
   ! z0tz0     z0t / z0 constant ratio over land

   integer j
   real ac, Z0MAX, VISCOSITY
   real ue(n)
   real trrn(n)
   real lat(n)
   real z0tz0

   data VISCOSITY / 1.461E-5/
   data Z0MAX / 5.E-3 /
   ! DATA Z0MIN          /1.5E-5/        ! in sfc_opt.ftn90
   !       DATA Z0HCON/ 4.E-5/           ! in sfc_opt.ftn90
   ac = 0.018
   z0tz0=0.2
 
   ! z0 and z0t computation depending on option optz0
 
   ! 0.  RETURN: no calculations done
 
   if (optz0 .eq. 0 ) then

      return


      ! 1. sea ice : z0 fixed and z0t=z0   (SEAICE classic)
   elseif (optz0 .eq. 1 ) then

      do j=1,n
         z0t(j)= z0(j)
      enddo

      ! 2. land surface : z0 fixed (related to surface covers) and z0t/z0 fixed ratio  (ISBA / FORCE-RESTORE classic)
      !         z0t still capped at 0.2
   elseif (optz0 .eq. 2 ) then

      do j=1,n
         z0t(j)=min(  z0(j) * z0tz0, 0.2)
      enddo

      ! 3. water surface : z0 from Charnock (1955) and z0t dependant on latitude (WATER classic)
      !     Note:  For |lat| >= Z0TLAT(2)  Charnock's relation is used
      !            For |lat| <= Z0TLAT(1)  Z0HCON is used.
      !            For Z0TLAT(1) < |lat| < Z0TLAT(2)
      !            we do a linear interpolation between Charnock and Z0HCON.

   elseif (optz0 .eq. 3 ) then

      do j=1,n
         ue(j)= karman/fm(j) *va(j)
         z0(j) = max( min( ac*ue(j)**2/GRAV,Z0MAX ) , Z0MIN )
         lat(j)= asin(fcor(j) / 1.45441e-4)
         if (abs(lat(j)) .ge. Z0TLAT(2)) then
            z0t(j)=z0(j)
         else if (abs(lat(j)) .le. Z0TLAT(1)) then
            z0t(j) = Z0HCON
         else
            z0t(j)=( ((abs(lat(j))-Z0TLAT(1))/(Z0TLAT(2)-Z0TLAT(1)))   &
                 *(Z0(j)-Z0HCON) ) + Z0HCON
         endif
      enddo

      ! 4. glaciers : z0 from geophy and z0t classic
   elseif (optz0 .eq. 4 ) then
      do j=1,n
         z0t(j) = max(   0.2 * z0(j)  ,  Z0GLA   )
      enddo

      ! 5. water surface : Deacu and Fortin (2012)
   elseif (optz0 .eq. 5 ) then
      do j=1,n
         ue(j)= karman/fm(j) * va(j)
         z0(j) = max( min( ac*ue(j)**2/GRAV,Z0MAX ) , Z0MIN )
         if (ue(j) .gt. 0.2) then
            z0t(j) = 1.333 * VISCOSITY /ue(j)
         else
            z0t(j) = 1.0E-4
         endif
      enddo


      ! 6. water surface : Fairall et al. 2003
      !            z0 from  Smith (1988),
      !            but different ac (between 0.011 et 0.2 in original paper)
   elseif (optz0 .eq. 6 ) then
      do j=1,n
         ue(j)= karman/fm(j) *va(j)
         z0(j) = max( min( (ac*ue(j)**2/GRAV + 0.11*         &
              VISCOSITY /ue(j)), Z0MAX ) , Z0MIN )
         trrn(j)=z0(j) * ue(j) /VISCOSITY
         z0t(j) = min(1.1e-04 ,5.5e-05*trrn(j)**(-0.6))
      enddo

      ! 7. urban surface: z0 from geophy and z0t from Brutsaert (1982)
   elseif (optz0 .eq. 7 ) then

      do j=1,n
         ue(j)= karman/fm(j) *va(j)
         trrn(j)=z0(j) * ue(j) /VISCOSITY
         z0t(j) = z0(j) * 7.4 * exp( - 2.46 * trrn(j)**0.25)
      enddo

      ! 8. urban surface: z0 from geophy and z0t from Kanda (2007)
   elseif (optz0 .eq. 8 ) then

      do j=1,n
         ue(j)= karman/fm(j) *va(j)
         trrn(j)=z0(j) * ue(j) /VISCOSITY
         z0t(j) = z0(j) * 7.4 * exp( - 1.29 * trrn(j)**0.25)
      enddo

      ! 9. land surface : z0 from geophy and z0t from Zilitinkevich (1995)
      !    here use local mom. roughness when available 
   elseif (optz0 .eq. 9 ) then
      do j=1,n
         ue(j)= karman/fm(j) *va(j)
         trrn(j)= z0loc(j) * ue(j) /VISCOSITY
         z0t(j) = max( z0loc(j) * exp( - KARMAN  * 0.1 * trrn(j)**0.5) , z0min )
      enddo


      ! 10. urban surface : z0 fixed (related to surface covers)
      !    and z0t/z0 fixed ratio=200 (Mascart et al. 1995)
   elseif (optz0 .eq. 10 ) then

      do j=1,n
         z0t(j)=z0(j) / 200.
      enddo

   endif

   return
end subroutine compz0_1
