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


module sfclayer_compz0
   use clib_itf_mod, only: clib_toupper
   implicit none
   private

#include <rmnlib_basics.hf>

   public :: compz0_set  !#, compz0

   integer, parameter, public :: Z0SCFTYPE_MINVAL = 1
   integer, parameter, public :: Z0SCFTYPE_SEAICE = 1
   integer, parameter, public :: Z0SCFTYPE_LAND = 2
   integer, parameter, public :: Z0SCFTYPE_WATER = 3
   integer, parameter, public :: Z0SCFTYPE_GLACIER = 4
   integer, parameter, public :: Z0SCFTYPE_WATER_DF12 = 5
   integer, parameter, public :: Z0SCFTYPE_WATER_FA03 = 6
   integer, parameter, public :: Z0SCFTYPE_URB_BR82 = 7
   integer, parameter, public :: Z0SCFTYPE_URB_KA07 = 8
   integer, parameter, public :: Z0SCFTYPE_LAND_ZI95 = 9
   integer, parameter, public :: Z0SCFTYPE_URB_MA95 = 10
   integer, parameter, public :: Z0SCFTYPE_MAXVAL = 10

   !#TODO: rm public
   real, public :: z0min = 1.5e-5
   real, public :: z0gla = 0.0003
   real, public :: z0tlat(2) = 0.
   real, public :: z0hcon = 4.0e-5

!!$   interface compz0
!!$      module procedure compz0_s
!!$      module procedure compz0_a
!!$   end interface compz0
   
contains

   function compz0_set(F_name_S, F_val, F_val2) result(F_istat)
      implicit none
      character(len=*), intent(in) :: F_name_s
      real, intent(in) :: F_val
      real, intent(in), optional :: F_val2
      integer :: F_istat
      
      character(len=32) :: name_S
      
      name_S = F_name_S
      F_istat = clib_toupper(name_S)
      select case(name_S)
      case ('Z0MIN')
         z0min = F_val
         return
      case ('Z0GLA')
         z0gla = F_val
         return
      case ('Z0HCON')
         z0hcon = F_val
         return
      case ('Z0LAT')
         if (present(F_val2)) then
            z0tlat(1) = F_val
            z0tlat(2) = F_val2
            return
         endif
      end select
      F_istat = RMN_ERR
      return   
   end function compz0_set

!!$   subroutine compz0_s(optz0, z0, z0loc, z0t, fm, va, fcor)
!!$      use tdpack_const, only: GRAV, KARMAN
!!$      implicit none
!!$      integer, intent(in) :: optz0
!!$      real, intent(in) :: z0loc, fm, va, fcor
!!$      real, intent(inout) :: z0
!!$      real, intent(inout) :: z0t !# Actually only out, but sometime not computed, hence inout
!!$
!!$      real, dimension(1) :: lz0loc, lfm, lva, lfcor, lz0, lz0t
!!$
!!$
!!$      if (optz0 < Z0SCFTYPE_MINVAL .or. optz0 > Z0SCFTYPE_MAXVAL) return
!!$      
!!$      lz0loc = z0loc
!!$      lfm = fm
!!$      lva = va
!!$      lfcor = fcor
!!$      lz0 = z0
!!$      lz0t = z0t
!!$      call compz0_a(optz0, lz0, lz0loc, lz0t, lfm, lva, lfcor, 1)
!!$      z0 = lz0(1)
!!$      z0t = lz0t(1)
!!$      
!!$      return
!!$   end subroutine compz0_s
!!$
!!$   subroutine compz0_a(optz0, z0, z0loc, z0t, fm, va, fcor, n)
!!$      use tdpack_const, only: GRAV, KARMAN
!!$      implicit none
!!$      integer, intent(in) :: n, optz0
!!$      real, intent(in) :: z0loc(n), fm(n), va(n), fcor(n)
!!$      real, intent(inout) :: z0(n)
!!$      real, intent(inout) :: z0t(n) !# Actually only out, but sometime not computed, hence inout
!!$
!!$      call compz0...(optz0, z0, z0loc, z0t, fm, va, fcor, n)
!!$      
!!$      return
!!$   end subroutine compz0_a
   
end module sfclayer_compz0


   subroutine compz0_a(optz0, z0, z0loc, z0t, fm, va, fcor, n)
      use tdpack_const, only: GRAV, KARMAN
      use sfclayer_compz0
      implicit none
!!!#include <arch_specific.hf>
      integer, intent(in) :: n, optz0
      real, intent(in) :: z0loc(n), fm(n), va(n), fcor(n)
      real, intent(inout) :: z0(n)
      real, intent(inout) :: z0t(n) !# Actually only out, but sometime not computed, hence inout

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
      ! z0t      scalar roughness lenght for heat
      !          - Input -
      ! fm       momentum stability function
      ! va       wind speed at ZU
      ! n        horizontal dimension
      ! optz0    option for z0 formulations
      ! fcor     Coriolis factor
      ! z0loc    local aerodynamic roughness length (no orography)
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

      real, parameter :: VISCOSITY = 1.461E-5
      real, parameter :: Z0MAX = 5.E-3
      real, parameter :: AC = 0.018
      real, parameter :: Z0TZ0 = 0.2

      integer :: j
      real :: ue(n), trrn(n), lat(n)

      ! z0 and z0t computation depending on option optz0

      !#TODO: check if this it is safe not to always return a z0 or z0t value
      if (optz0 < Z0SCFTYPE_MINVAL .or. optz0 > Z0SCFTYPE_MAXVAL) return

      select case (optz0)

      case (Z0SCFTYPE_SEAICE)
         ! 1. sea ice : z0 fixed and z0t=z0   (SEAICE classic)

         z0t(1:n) = z0(1:n)

      case (Z0SCFTYPE_LAND)
         ! 2. land surface : z0 fixed (related to surface covers) and z0t/z0 fixed ratio  (ISBA / FORCE-RESTORE classic)
         !         z0t still capped at 0.2

         do j=1,n
            z0t(j) = min(z0(j) * z0tz0, 0.2)
         enddo

      case (Z0SCFTYPE_WATER)
         ! 3. water surface : z0 from Charnock (1955) and z0t dependant on latitude (WATER classic)
         !     Note:  For |lat| >= Z0TLAT(2)  Charnock's relation is used
         !            For |lat| <= Z0TLAT(1)  Z0HCON is used.
         !            For Z0TLAT(1) < |lat| < Z0TLAT(2)
         !            we do a linear interpolation between Charnock and Z0HCON.

         do j=1,n
            ue(j) = karman / fm(j) * va(j)
            z0(j) = max(Z0MIN, min(Z0MAX, ac*ue(j)**2/GRAV))
            lat(j) = asin(fcor(j) / 1.45441e-4)
            if (abs(lat(j)) >= Z0TLAT(2)) then
               z0t(j) = z0(j)
            else if (abs(lat(j)) <= Z0TLAT(1)) then
               z0t(j) = Z0HCON
            else
               z0t(j) = ( ((abs(lat(j))-Z0TLAT(1))/(Z0TLAT(2)-Z0TLAT(1)))   &
                    *(Z0(j)-Z0HCON) ) + Z0HCON
            endif
         enddo

      case (Z0SCFTYPE_GLACIER)
         ! 4. glaciers : z0 from geophy and z0t classic

         do j=1,n
            z0t(j) = max(0.2*z0(j), Z0GLA)
         enddo

      case (Z0SCFTYPE_WATER_DF12)
         ! 5. water surface : Deacu and Fortin (2012)

         do j=1,n
            ue(j) = karman / fm(j) * va(j)
            z0(j) = max(Z0MIN, min(Z0MAX, ac*ue(j)**2/GRAV))
            if (ue(j) > 0.2) then
               z0t(j) = 1.333 * VISCOSITY / ue(j)
            else
               z0t(j) = 1.0E-4
            endif
         enddo

      case (Z0SCFTYPE_WATER_FA03)
         ! 6. water surface : Fairall et al. 2003
         !            z0 from  Smith (1988),
         !            but different ac (between 0.011 et 0.2 in original paper)

         do j=1,n
            ue(j) = karman / fm(j) * va(j)
            z0(j) = max(Z0MIN, min(Z0MAX, &
                 ac*ue(j)**2/GRAV + 0.11*VISCOSITY/ue(j) &
                 ))
            trrn(j) = z0(j) * ue(j) / VISCOSITY
            z0t(j) = min(1.1e-04, 5.5e-05*trrn(j)**(-0.6))
         enddo

      case (Z0SCFTYPE_URB_BR82)
         ! 7. urban surface: z0 from geophy and z0t from Brutsaert (1982)

         do j=1,n
            ue(j) = karman / fm(j) * va(j)
            trrn(j) = z0(j) * ue(j) / VISCOSITY
            z0t(j) = z0(j) * 7.4 * exp(-2.46 * trrn(j)**0.25)
         enddo

      case (Z0SCFTYPE_URB_KA07)
         ! 8. urban surface: z0 from geophy and z0t from Kanda (2007)

         do j=1,n
            ue(j) = karman / fm(j) * va(j)
            trrn(j) = z0(j) * ue(j) / VISCOSITY
            z0t(j) = z0(j) * 7.4 * exp(-1.29 * trrn(j)**0.25)
         enddo

      case (Z0SCFTYPE_LAND_ZI95)
         ! 9. land surface : z0 from geophy and z0t from Zilitinkevich (1995)
         !    here use local mom. roughness when available 

         do j=1,n
            ue(j) = karman / fm(j) * va(j)
            trrn(j) = z0loc(j) * ue(j) / VISCOSITY
            z0t(j) = max(Z0MIN, z0loc(j) * exp(-KARMAN * 0.1 * trrn(j)**0.5))
         enddo

      case (Z0SCFTYPE_URB_MA95)
         ! 10. urban surface : z0 fixed (related to surface covers)
         !    and z0t/z0 fixed ratio=200 (Mascart et al. 1995)

         z0t(1:n) = z0(1:n) / 200.

      end select

      return
   end subroutine compz0_a
