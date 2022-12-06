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
 
module diagnosurf
   implicit none
   private
   public :: diagnosurf5

contains

   !/@*
   subroutine diagnosurf5(ni, trnch)
      use series_mod, only: series_xst, series_isstep, series_isvar
      use sfc_options
      use sfcbus_mod
      implicit none
!!!#include <arch_specific.hf>

      integer, intent(in) :: ni, trnch

      !@Author  B. Bilodeau (Sept 1999)
      !@Revisions
      ! 001      B. Bilodeau (Nov 2000) - New comdeck phybus.cdk
      ! 002      B. Dugas               - Add S5
      ! 003      L. Spacek   (Oct 2010) - Replace UD/VD by UDWE/VDSN
      !                                   Add WSPD/WD (speed/diresction)
      !@Object Time-series and zonal diagnostics extractions
      !       of surface-related variables
      !@Arguments
      !          - Input -
      ! NI       horizontal dimensions of fields
      !          number of elements processed in the horizontal
      !          (could be a subset of M in future versions)
      !          (not used for the moment)
      ! TRNCH    row number
      !*@/
#include <rmn/msg.h>

#define PTR1D(NAME2,IDX) busptr(vd%NAME2%i)%ptr(1+IDX:1+IDX+ni-1,trnch)

      if (.not.series_isstep()) return

      call msg_toall(MSG_DEBUG, 'diagnosurf [BEGIN]')
      if (timings_L) call timing_start_omp(490, 'diagnosurf', 46)

      ! inverse of Monin-Obukhov length
      call series_xst(PTR1D(ILMO,(indx_agrege-1)*NI), 'IL', trnch)
      ! specific humidity at the surface
      call series_xst(PTR1D(QSURF,(indx_agrege-1)*NI), 'QG', trnch)

      ! glacier surface temperature
      call series_xst(PTR1D(TGLACIER,NI), '9I', trnch)

      ! marine ice temperature (3 levels)
      call series_xst(PTR1D(TMICE,NI),  '7I', trnch)
      call series_xst(PTR1D(TMICE,2*NI),'8I', trnch)

      ! snow thickness (soil, glacier and marine ice)
      call series_xst(PTR1D(SNODP, (indx_soil-1)*NI), 'S1', trnch)
      call series_xst(PTR1D(SNODP, (indx_glacier-1)*NI), 'S2', trnch)
      call series_xst(PTR1D(SNODP, (indx_ice-1)*NI), 'S4', trnch)
      call series_xst(PTR1D(SNODP, (indx_agrege-1)*NI), 'S5', trnch)

      ! surface albedo (soil, glacier, marine ice, water and average)
      call series_xst(PTR1D(ALVIS, (indx_soil-1)*NI), 'XS', trnch)
      call series_xst(PTR1D(ALVIS, (indx_glacier-1)*NI), 'XG', trnch)
      call series_xst(PTR1D(ALVIS, (indx_water-1)*NI), 'XW', trnch)
      call series_xst(PTR1D(ALVIS, (indx_ice-1)*NI), 'XI', trnch)
      call series_xst(PTR1D(ALVIS, (indx_agrege-1)*NI), 'AL', trnch)

      ! soil temperature
      call series_xst(PTR1D(TSURF,(indx_agrege-1)*NI), 'TS', trnch)
      call series_xst(PTR1D(TSOIL, 0), 'TS1', trnch)

      ! radiative surface temperature
      call series_xst(PTR1D(TSRAD, 0), 'G3', trnch) !TG

      ! diagnostic U component of the wind at screen level
      call series_xst(PTR1D(UDIAG, 0), 'UDWE', trnch) !UD
      call series_xst(PTR1D(UDIAG, 0), 'WSPD', trnch) !UD

      ! diagnostic V component of the wind at screen level
      call series_xst(PTR1D(VDIAG, 0), 'VDSN', trnch) !VD
      call series_xst(PTR1D(VDIAG, 0), 'WD', trnch) !VD


      IF_ISBA: if (schmsol == 'ISBA') then              ! isba

         ! deep soil temperature
         call series_xst(PTR1D(TSOIL, NI), 'TP', trnch)

         ! soil moisture content
         call series_xst(PTR1D(WSOIL, 0), 'WG', trnch) !I1

         ! deep soil moisture content
         call series_xst(PTR1D(WSOIL, NI), 'WR', trnch)

         ! liquid water stored on canopy
         call series_xst(PTR1D(WVEG, 0), 'C5', trnch) !I3

         ! mass of snow cover
         call series_xst(PTR1D(SNOMA, 0), 'C6', trnch) !I5

         ! snow albedo
         call series_xst(PTR1D(SNOAL, 0), 'C7', trnch) !I6

         ! snow density
         call series_xst(PTR1D(SNORO, 0), 'C8', trnch) !7S

         ! net radiation
         call series_xst(PTR1D(RNET_S, 0), 'C9', trnch) !NR

         ! latent heat flux over bare ground
         call series_xst(PTR1D(LEG, 0), 'D3', trnch) !L2

         ! latent heat flux over vegetation
         call series_xst(PTR1D(LEV, 0), 'D4', trnch) !LV

         ! latent heat flux over snow
         call series_xst(PTR1D(LES, 0), 'D5', trnch) !LS

         ! direct latent heat flux from vegetation leaves
         call series_xst(PTR1D(LER, 0), 'D6', trnch) !LR

         ! latent heat of evapotranspiration
         call series_xst(PTR1D(LETR, 0), 'D7', trnch) !LT

         ! runoff
         call series_xst(PTR1D(OVERFL, 0), 'E2', trnch) !RO

         ! drainage
         call series_xst(PTR1D(DRAIN, 0), 'E3', trnch) !DR

         ! fraction of the grid covered by snow
         call series_xst(PTR1D(PSN, 0), 'E5', trnch) !5P

         ! fraction of bare ground covered by snow
         call series_xst(PTR1D(PSNG, 0), 'E6', trnch) !3P

         ! fraction of vegetation covered by snow
         call series_xst(PTR1D(PSNV, 0), 'E7', trnch) !4P

         ! stomatal resistance
         call series_xst(PTR1D(RST, 0), 'E8', trnch) !R1

         ! specific humidity of the surface
         call series_xst(PTR1D(HUSURF, 0), 'E9', trnch) !FH

         ! Halstead coefficient (relative humidity of veg. canopy)
         call series_xst(PTR1D(HV, 0), 'G1', trnch) !HV

         ! soil volumetric ice content
         call series_xst(PTR1D(ISOIL, 0), 'G4', trnch) !I2

         ! liquid water in snow
         call series_xst(PTR1D(WSNOW, 0), 'G5', trnch) !I4

         ! liquid precip. rate
         call series_xst(PTR1D(RAINRATE, 0), 'G6', trnch) !U1

         ! solid precip. rate
         call series_xst(PTR1D(SNOWRATE, 0), 'G7', trnch) !U3

      endif IF_ISBA

      if (timings_L) call timing_stop_omp(490)
      call msg_toall(MSG_DEBUG, 'diagnosurf [END]')

      return
   end subroutine diagnosurf5

end module diagnosurf
