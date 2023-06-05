!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2013 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------

!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet/Project : GEM-MACH
! Fichier/File   : mach_canopy_levels.ftn90
! Creation       : Paul Makar - 2015
! Description    : Uses canopy height (hc), and fractions of canopy height to
!                  gather meteorological and geophysical data into vertical
!                  slices containing no canopy (these have the same vertical
!                  dimension as the original ungathered model) and combined
!                  canopy + no-canopy columns (these have the original,
!                  "resolved scale" model layers interleaved with nkc canopy
!                  layers.  For the canopy layers, various algorithms are used
!                  to fill in the missing meteorological data.
!
! Extra info     : On the first timestep, the subroutine also initializes the
!                  concentrations of the canopy layers using the value of the
!                  lowest resolved scale model layer.
!
! Arguments:
! Array dimensions:
!  nkc                        :  number of canopy levels
!  nkt                        :  number of levels in gathered canopy + resolved scale columns = chm_nk+nkc
!  ni_can                     :  number of columns containing a forest canopy
!  busper                     :  permanent bus
!  can_frac(nkc)              :  fractions of canopy height for canopy layers;
!                                can_frac * hc = height above sfc for canopy layer
!  metvar2d(chm_ni, SIZE_MV2D):  met 2d variables, resolved model layers and unscattered horizontal coordinate
!  metvar3d(chm_ni, chm_nk, SIZE_MV3D):  met 3d variables, resolved model layers and unscattered horizontal coordinate
!
! Output variables, gathered canopy columns:
!
!  metvar3dcan(ni_can, nkt, SIZE_MV3D)   :  met 3d variables, gathered canopy + resolved scale columns
!
!  Output variables, original horizontal coordinate
!   busper                   :  permanent bus
!   zcan(ni_can,nkc)         :  Height above ground of nkc canopy layers
!   kmod(ni_can,chm_nk)      :  Vertical index location of original ungathered model layer in combined
!                               canopy + resolved scale column
!   kcan(ni_can,nkc)         :  Vertical index location of canopy ungathered model layer in combined
!                               canopy + resolved scale column
!
!  Local variables:
!
!    hc(ni_can)              :  canopy height in m
!    klower_can(ni_can, nkc) :  Index of the resolved scale layer just below the given canopy layer
!    del                     :  Minimum allowable distance between a resolved model layer and a
!                               canopy layer (fraction of canopy layer height)
!  zt(ni_can, chm_nk)        :  (aka gz_chm) thermodynamic heights  (m asl), resolved model layers
!  sigmcan (ni_can, nkt)  :  Sigma coordinate for momentum heights in gathered canopy + resolved scale columns
!  sigtcan (ni_can, nkt)  :  Sigma coordinate for thermodynamic levels in gathered canopy + resolved scale columns
!  zmomcan (ni_can, nkt)  :  Momentum heights in gathered canopy + resolved scale columns
!  zthrmcan(ni_can, nkt)  :  Thermodynamic heights in gathered forest canopy columns
!=============================================================================
!
!!if_on
subroutine mach_canopy_levels(busper, metvar3dcan, metvar3dnocan, metvar3d, &
                            metvar2d, imod, imod2, kmod, kcan, ni_can, ni_nocan)
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk, nkt, nkc
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
!!if_off
   use chm_utils_mod,        only: ik, chm_error_l, global_debug
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_HC, sp_FRT
   use chm_consphychm_mod,   only: karman, pi, rgasd, delta
   use chm_metvar_mod,       only: MV2D_ILMO, MV2D_PPLUS, MV2D_QDIAG, MV2D_UE, &
                                   MV2D_TDIAG, MV3D_FTOT, MV3D_KT, MV3D_TPLUS, &
                                   MV3D_SIGM, MV3D_SIGT, MV3D_PEVP, MV3D_PPRO, &
                                   MV3D_ZMOM, MV3D_ZPLUS, MV3D_WS, MV3D_RHO,  &
                                   MV3D_QCPLUS, MV3D_HUPLUS, MV3D_CLDRAD,      &
                                   MV3D_RNFLX, MV3D_SNOFLX, MV3D_LWC

   implicit none
!!if_on
   integer(kind=4), intent   (in) :: ni_can, ni_nocan
   integer(kind=4), intent  (out) :: imod(ni_can)
   integer(kind=4), intent  (out) :: imod2(ni_nocan)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
   real(kind=4),    intent  (out) :: metvar3dcan(ni_can, nkt, SIZE_MV3D)
   real(kind=4),    intent  (out) :: metvar3dnocan(ni_nocan, chm_nk, SIZE_MV3D)
   integer(kind=4), intent  (out) :: kmod(ni_can, chm_nk)
   integer(kind=4), intent  (out) :: kcan(ni_can, nkc)

!!if_off
!
! Local variables:
   integer(kind=4) :: i, k, kk, npass, ii, ic, kc
   integer(kind=4) :: k2, n
   logical(kind=4) :: flag_error
   real(kind=4)    :: tmp, presscan
   real(kind=4)    :: zm, td, hd, ddel
   real (kind=4)   :: zr, uh, uspr, wndr, sigw, tl, ktr
   real(kind=4),    parameter :: del = 0.2
   real(kind=4),    parameter :: min_kt = 0.1
   real(kind=4), dimension(ni_can)           :: hc
   real(kind=4), dimension(ni_can)           :: safe_inv_mo_length
   real(kind=4), dimension(ni_can, chm_nk)   :: zt
   real(kind=4), dimension(ni_can, nkc)      :: zcan
   real(kind=4), dimension(ni_can, nkt)      :: sigtcan, sigmcan, zthrmcan, zmomcan
   real(kind=4), dimension(ni_can, chm_nk+1) :: sigt2, z2
   integer(kind=4), dimension(ni_can)        :: ka, kcan_top
   integer(kind=4), dimension(ni_can, nkc)   :: klower_can(ni_can, nkc)

   real(kind=4) :: zl, a1, b1, c1, rat
   real(kind=4), parameter :: THRESHOLD = 1.e06 ! MOL threshold, similar to mach_plumerise

!  Assign the fractional heights of the canopy layers (fraction of canopy height)
   real(kind=4), dimension(3), parameter   :: can_frac = (/1.0, 0.5, 0.2/)

   logical(kind=4)                         :: local_dbg

   local_dbg = (.false. .or. global_debug)
!
! Initializations: Define canopy and non-canopy grid indirect addresses
   ic = 0
   ii = 0
   do i = 1, chm_ni
      if (busper(sm(sp_FRT) % per_offset + i - 1) > 0.0) then
         ic = ic + 1
         imod(ic) = i
      else
         ii = ii + 1
         imod2(ii) = i
!  Re-pack the metvar3d fields in non-canopy grids
         do n = 1, SIZE_MV3D
            do k = 1, chm_nk
               metvar3dnocan(ii, k, n) = metvar3d(i, k, n)
            end do
         end do
      end if
   end do
!
   do ic = 1, ni_can
      ii = imod(ic)
!
! Generate initial canopy levels, as altitude above sea level
      hc(ic) = busper(sm(sp_HC) % per_offset + ii - 1)
!
      do kc = 1, nkc
         zcan(ic, kc) = hc(ic) * can_frac(kc)
      end do

      do k = 1, chm_nk
         zt(ic, k) =  metvar3d(ii, k, MV3D_ZPLUS)
      end do

! The model level above the tallest canopy in grid
      kcan_top(ic) = 2
      do k = chm_nk, 3, -1
         if (zt(ic, k) > hc(ic)) then
            kcan_top(ic) = k - 1
            exit
         end if
      end do
!
!  Setup of Monin-Obhukov Length similar to plumerise for upper limit:
      safe_inv_mo_length(ic) = metvar2d(ii, MV2D_ILMO)
      if (abs(metvar2d(ii, MV2D_ILMO)) > THRESHOLD) then
         safe_inv_mo_length(ic) = sign(THRESHOLD, metvar2d(ii, MV2D_ILMO))
      end if
!
! Adjust the canopy levels: we don't want canopy levels to get closer than del
! to the model levels to prevent possible differencing errors in the diffusion.
!  If zcan > zt but is too close to zt, move zcan up by ddel.  If zcan < zt but
!  is too close to zt, move zcan down by ddel.  The net result will be that the
!  canopy levels are never closer than del from the original model levels.
      do k = kcan_top(ic), chm_nk
         do kc = 1, nkc
            if (abs(zt(ic, k) - zcan(ic, kc)) < del) then
               ddel = max(0.0, del - abs(zcan(ic, kc) - zt(ic, k)))
               zcan(ic, kc) = zcan(ic, kc) + sign(ddel, zcan(ic, kc) - zt(ic, k))
            end if
         end do
      end do
   end do
!
!  Set the initial values of the combined height array:
!
! Note that here, zthrmcan is created, but the heights within each column have
! yet to be sorted to rearrange the layers in the correct order.
   do ic = 1, ni_can
      do k = 1, chm_nk
         zthrmcan(ic, k) = zt(ic, k)
      end do
! Add zcan additional thermo levels into zthrmcan array for later sorting
      do kc = 1, nkc
         zthrmcan(ic, chm_nk+kc) = zcan(ic, kc)
      end do
   end do
!
!  Determine locations of canopy and resolved model levels within
!  the combined array for the canopy columns:
   do ic = 1, ni_can
      if (zthrmcan(ic, chm_nk) < zthrmcan(ic, chm_nk+1)) then
!
!  Non-trivial case:  the ancilliary and original array levels intermingle.
!  Sort the combined height array to get the right order of the the heights:
!
!  zthrmcan is the height locations of the combined array, which needs to be
!  sorted: since there are only nkc levels in the canopy, and both zcan and z
!  decrease monotonically, only nkc+1 passes are needed to sort the combined
!  array:
         do npass = 1, nkc+1
            flag_error = .false.
            do k = nkt, 2, -1
               if (zthrmcan(ic, k) > zthrmcan(ic, k-1)) then
!  The combined array heights are out of order, sort them:
                  tmp = zthrmcan(ic, k-1)
                  zthrmcan(ic, k-1) = zthrmcan(ic,k)
                  zthrmcan(ic, k)   = tmp
                  flag_error = .true.
               end if
            end do
         end do
         if (flag_error) then
            write(6,*) 'nkc+1 passes insufficient to sort canopy array '
            write(6,*) 'in mach_canopy_level.ftn90.  Scream and die.'
            chm_error_l = .true.
            return
         end if
      end if
   end do

!
!  Heights in zthrmcan should now be monotonically decreasing.
!
!  Next, identify the locations of the vertical levels in the combined
!  array relative to the resolved model array and canopy array
   kcan = -999
   kmod = -999
   do ic = 1, ni_can
      do kc = 1, nkc
         do kk = nkt, 1, -1
            if (zthrmcan(ic, kk) == zcan(ic,kc)) then
               kcan(ic, kc) = kk
               exit
            endif
         end do
      end do
      do k = 1, chm_nk
         do kk = k, nkt
            if (zthrmcan(ic, kk) == zt(ic, k)) then
               kmod(ic, k) = kk
               exit
            endif
         end do
      end do

      if (local_dbg) then
      do kc = 1, nkc
         if (kcan(ic, kc) < 1) then
            write(6,*) 'mach_canopy_levels: kcan undefined: ', &
                       ic, kc, kcan(ic, kc)
            chm_error_l = .true.
            return
         end if
      end do
      do k = 1, chm_nk
         if (kmod(ic, k) < 1) then
            write(6,*) 'mach_canopy_levels: kmod undefined: ',ic, k, kmod(ic, k)
            chm_error_l = .true.
            return
         end if
      end do
      end if
   end do
!
!  Create the corresponding momentum height array
!  The original methodology adopted made use of the at2m array and the thermodynamic heights determined above.
!  However, this methodology resulted in momentum levels which did not match the original model levels
!  above the region modified for canopy layers.  Here, the thermodynamic layers will be used to
!  (1) Determine whether the original model and canopy thermodynamic layers coincide, and if so,
!  (2) Use the existing model layer values for the momentum layers, while if not,
!  (3) Assign the new momentum layers as being 1/2 way between the canopy layers
! Note that these changes only exist inside the chemistry part of GEM-MACH and do not affect the
! model physics
   do ic = 1, ni_can
      ii = imod(ic)

!      Default case:  all added canopy thermodynamic layers are
!      below the lowest resolved model thermodynamic layer
      do k = 1, kcan_top(ic) - 1
         zmomcan(ic, k) = metvar3d(ii, k, MV3D_ZMOM)
      end do
      ka(ic) = chm_nk
      inner0: do k = kcan_top(ic), chm_nk-1
!  Starting from the top, scan down through the original and combined thermo levels, to see when they
!  first deviate from each other
         if (zthrmcan(ic, k) == zt(ic, k) .and. zthrmcan(ic, k+1) == zt(ic, k+1) ) then
            zmomcan(ic, k) = metvar3d(ii, k, MV3D_ZMOM)
         else
            ka(ic) = k
            exit inner0
         end if
      end do inner0
!
! ka  is the last layer for which zmomcan = metvar3d(imod(ic),k,MV3D_ZMOM)
      zmomcan(ic, ka(ic)) = metvar3d(ii, ka(ic), MV3D_ZMOM)
      do k = ka(ic)+1, nkt
         zmomcan(ic, k) = (zthrmcan(ic, k-1) + zthrmcan(ic, k)) * 0.5
      end do
   end do
!
!  create original model arrays of z and sigma-t which include the surface, to
!  allow interpolation:
   sigtcan = 0.0
   do k = 1, chm_nk
      do ic = 1, ni_can
         sigt2(ic, k) = metvar3d(imod(ic), k, MV3D_SIGT)
         z2(ic, k)    = zt(ic, k)
!
!  Fill in the thermodynamic sigma levels (Pre-existing levels first):
         sigtcan(ic, kmod(ic, k)) = sigt2(ic, k)
      end do
   end do
   klower_can = -999
   do ic = 1, ni_can
      ii = imod(ic)
      sigt2(ic, chm_nk+1) = 1.0
      z2(ic, chm_nk+1)    = 0.0
!
!  fill in the remaining sigma levels by interpolating in z:
      do kc = 1, nkc
         do k2 = kcan_top(ic), chm_nk+1
            if (zcan(ic, kc) > z2(ic, k2) .and. zcan(ic, kc) <= z2(ic, k2-1)) then
! Interpolate in sigma
               sigtcan(ic, kcan(ic, kc)) = metvar3d(ii, k2-1, MV3D_SIGT)  +     &
                                           (sigt2(ic, k2) - sigt2(ic, k2-1)) /  &
                                           (z2(ic, k2) - z2(ic, k2-1)) *        &
                                           (zcan(ic, kc) - z2(ic, k2-1))
! Store grid locations for use in later interpolations
               klower_can(ic, kc) = k2
            end if
         end do
!
!
         if (klower_can(ic, kc) < 1) then
            write(6,*) 'mach_canopy_levels:  klower_can is unassigned at ii ic kc: ',ii ,ic, kc
            write(6,*) 'mach_canopy_levels:  zcan(ic, kc): ',zcan(ic,kc)
            do kk = kcan_top(ic), chm_nk+1
               write(6,*) 'mach_canopy_levels: kk z2(ic kk) which should bracket the above zcan: ',kk, z2(ic,kk)
            end do
            do kk = 1, chm_nk+1
               write(6,*) 'mach_canopy_levels:  kk z2(ic kk) full set of z2 values: ', kk, z2(ic,kk)
            end do
            do kk = 1, nkc
               write(6,*) 'mach_canopy_levels:  kc zcan(ic kc) hc(ic) fr(kc) for full set of zcan values: ',kk, zcan(ic, kk), hc(ic), can_frac(kk)
            end do
            chm_error_l = .true.
            return
         end if
      end do
   end do

!
   if (local_dbg) then
!  Check on klower_can for NaN or out of bounds:
   do kc = 1, nkc
      do ic = 1, ni_can
         if ((klower_can(ic, kc) /= klower_can(ic, kc)) .or. &
            (klower_can(ic, kc) <= 0)             .or. &
            (klower_can(ic, kc) > chm_nk + 1) ) then
            write(6,*) 'mach_canopy_levels klower_can after creation NaN or <=0 or >chm_nk+1 : ', &
                       ic, kc, klower_can(ic, kk)
            chm_error_l = .true.
            return
         end if
      end do
   end do
   end if
!
!  Create sigma coordinate  momentum levels:
!
!  As above, the existing momentum levels and the canopy values are used to create SIGM levels
!
!  (1) Determine whether the original model and canopy thermodynamic layers coincide, and if so,
!  (2) Use the existing model layer values for the momentum layers, while if not,
!  (3) Assign the new momentum layers as being 1/2 way between the canopy layers
! Note that these changes only exist inside the chemistry part of GEM-MACH and do not affect the
! model physics

   do ic = 1, ni_can
      ii = imod(ic)
!      Default case:  all added canopy thermodynamic layers are
!      below the lowest resolved model thermodynamic layer
      ka(ic) = chm_nk
      inner2:   do k = 1, chm_nk-1
         if (sigtcan(ic, k) == sigt2(ic, k) .and. sigtcan(ic, k+1) == sigt2(ic, k+1) ) then
            sigmcan(ic,k) = metvar3d(ii, k, MV3D_SIGM)
         else
            ka(ic) = k
            exit inner2
         end if
      end do inner2
! ka  is the last layer for which sigmcan= metvar3d(imod(ic),k,MV3D_SIGM)
      sigmcan(ic, ka(ic)) = metvar3d(ii, ka(ic), MV3D_SIGM)
      do k = ka(ic)+1, nkt
         sigmcan(ic, k) = (sigtcan(ic, k-1) + sigtcan(ic, k)) * 0.5
      end do
   end do
!
!  Next, do a sort of all of the variables in the original METV3D array into
!  canopy.  Note that the declaration of the met arrays for the new canopy
!  subdomain has occurred earlie in the code.
!  Three-D variables are a bit more complicated, in that one must make
!  decisions regarding the values of the met variables in the canopy region.
!  The code which follows is based on chm_load_metvar.ftn90
!
   metvar3dcan = -999.  !initialize for error checking
!
!  First, carry over original model values for the matching layers
   do n = 1, size_mv3d
      do k = 1, chm_nk
         do ic = 1, ni_can
            metvar3dcan(ic, kmod(ic,k), n) = metvar3d(imod(ic), k, n)
         end do
      end do
   end do
!  Assign pre-calculated fields; MV3D_ZPLUS, MV3D_ZMOM, MV3D_SIGM, and MV3D_SIGT
   do k = 1, nkt
      do ic = 1, ni_can
         metvar3dcan(ic, k, MV3D_ZPLUS) = zthrmcan(ic, k)
         metvar3dcan(ic, k, MV3D_ZMOM)  = zmomcan(ic, k)
         metvar3dcan(ic, k, MV3D_SIGM)  = sigmcan(ic, k)
         metvar3dcan(ic, k, MV3D_SIGT)  = sigtcan(ic, k)
      end do
   end do
!
!----------------------------------------------------------------------------
!  Canopy region:  next, go through each variable to work out canopy values.
!
!  (1) Do those variables for which special canopy formulae will NOT be used:
   do kc = 1, nkc
      do ic = 1, ni_can
!  Each of the following four variables have a screen height (2m) value
!  in the 2D met arrays
!          Temperature:         TPLUS, TDIAG
!          Specific humidity:  HUPLUS, QDIAG
         ii = imod(ic)
         kk = kcan(ic, kc)
         if (klower_can(ic, kc) <= chm_nk) then
            k2 = klower_can(ic, kc)
            zm = (zcan(ic, kc) - z2(ic, k2-1)) / (z2(ic, k2) - z2(ic, k2-1))
            td = (metvar3d(ii, k2, MV3D_TPLUS)  - metvar3d(ii, k2-1, MV3D_TPLUS))   * zm
            hd = (metvar3d(ii, k2, MV3D_HUPLUS) - metvar3d(ii, k2-1, MV3D_HUPLUS))  * zm
!
            metvar3dcan(ic, kk, MV3D_TPLUS)  = metvar3d(ii, k2-1, MV3D_TPLUS)  + td
            metvar3dcan(ic, kk, MV3D_HUPLUS) = metvar3d(ii, k2-1, MV3D_HUPLUS) + hd
         else
            if (zcan(ic, kc) - z2(ic, chm_nk+1) >= 2.0) then
!  Level is below first resolved model level but above screen height
               zm = (zcan(ic, kc) - z2(ic, chm_nk+1) - 2.0) / (z2(ic, chm_nk) - z2(ic, chm_nk+1) - 2.0)
               td = (metvar3d(ii, chm_nk, MV3D_TPLUS)  - metvar2d(ii, MV2D_TDIAG))  * zm
               hd = (metvar3d(ii, chm_nk, MV3D_HUPLUS) - metvar2d(ii, MV2D_QDIAG))  * zm
               metvar3dcan(ic, kk, MV3D_TPLUS)  = metvar2d(ii, MV2D_TDIAG) + td
               metvar3dcan(ic, kk, MV3D_HUPLUS) = metvar2d(ii, MV2D_QDIAG) + hd
            else
! Level in canopy is below screen height; assume constant values below screen height
               metvar3dcan(ic, kk, MV3D_TPLUS)  = metvar2d(ii, MV2D_TDIAG)
               metvar3dcan(ic, kk, MV3D_HUPLUS) = metvar2d(ii, MV2D_QDIAG)
            end if
         end if
!
!  The following variables are assumed to have uniform values throughout the
!  lowest resolved model layer:
!
!    Cloud liquid water mass mixing ratio (QCPLUS)
!    Total cloud fraction (FTOT)
!    Stratospheric cloud fraction (FXP)
!    Convective cloud fraction (FDC)
!    Total liquid water flux (RNFLX)
!    Total solid water flux (SNOFLX)
!    Precipitation evaporation (FEVP)
!    Cloud to rain collection tendency (PPRO)
!  Search over the original model layers (k).  Note that the outer loop above this
!  one is over the canopy layers kc:  we are looking for the values to assign the
!  canopy layers in the combined canopy+resolved scale space.  For these variables,
!  the resolved scale values will be used, hence the aim is to determine the
!  resolved scale layer in which the canopy layer resides, and assign the
!  corresponding values to the locations of the canopy layers in the combined
!  canopy + resolved scale space (kk).
         do k = kcan_top(ic), chm_nk
!  If the canopy layer kc falls within a given resolved scale layer:
            if (zcan(ic, kc) < metvar3d(ii, k-1, MV3D_ZMOM) .and. &
                zcan(ic, kc) >= metvar3d(ii, k, MV3D_ZMOM)) then
!  Assign those resolved scale values to the canopy layers in the resolved + canopy space (kk):
               metvar3dcan(ic, kk, MV3D_LWC)     = metvar3d(ii, k-1, MV3D_LWC)
               metvar3dcan(ic, kk, MV3D_CLDRAD)  = metvar3d(ii, k-1, MV3D_CLDRAD)
               metvar3dcan(ic, kk, MV3D_QCPLUS)  = metvar3d(ii, k-1, MV3D_QCPLUS)
               metvar3dcan(ic, kk, MV3D_FTOT)    = metvar3d(ii, k-1, MV3D_FTOT)
               metvar3dcan(ic, kk, MV3D_RNFLX)   = metvar3d(ii, k-1, MV3D_RNFLX)
               metvar3dcan(ic, kk, MV3D_SNOFLX)  = metvar3d(ii, k-1, MV3D_SNOFLX)
               metvar3dcan(ic, kk, MV3D_PEVP)    = metvar3d(ii, k-1, MV3D_PEVP)
               metvar3dcan(ic, kk, MV3D_PPRO)    = metvar3d(ii, k-1, MV3D_PPRO)
            end if
         end do
!  If the canopy layer kc falls below the lowest resolved scale layer...
         if (zcan(ic,kc) < metvar3d(ii, chm_nk, MV3D_ZMOM)) then
            metvar3dcan(ic, kk, MV3D_LWC)    = metvar3d(ii, chm_nk, MV3D_LWC)
            metvar3dcan(ic, kk, MV3D_CLDRAD) = metvar3d(ii, chm_nk, MV3D_CLDRAD)
            metvar3dcan(ic, kk, MV3D_QCPLUS) = metvar3d(ii, chm_nk, MV3D_QCPLUS)
            metvar3dcan(ic, kk, MV3D_FTOT)   = metvar3d(ii, chm_nk, MV3D_FTOT)
            metvar3dcan(ic, kk, MV3D_RNFLX)  = metvar3d(ii, chm_nk, MV3D_RNFLX)
            metvar3dcan(ic, kk, MV3D_SNOFLX) = metvar3d(ii, chm_nk, MV3D_SNOFLX)
            metvar3dcan(ic, kk, MV3D_PEVP)   = metvar3d(ii, chm_nk, MV3D_PEVP)
            metvar3dcan(ic, kk, MV3D_PPRO)   = metvar3d(ii, chm_nk, MV3D_PPRO)
         end if
!
!  Evaluate the air density in canopy columns using values determined above
         presscan        = sigtcan(ic, kk) * metvar2d(ii, MV2D_PPLUS)
         metvar3dcan(ic, kk, MV3D_RHO) = presscan / (rgasd * metvar3dcan(ic, kk, MV3D_TPLUS) * &
                              (1.0 + delta * metvar3dcan(ic, kk, MV3D_HUPLUS)) )
!
      end do !  chm_ni
   end do ! kc

   if (local_dbg) then
! Several checks for suspicious values:
   do kk = 1, nkt
      do ic = 1, ni_can
         if (metvar3dcan(ic, kk, MV3D_TPLUS) < 150.0) then
            write(6,*) 'mach_canopy_levels:  suspicious temperature detected in mach_canopy_levels after creation (ic kk value): ',&
                        ic, kk, metvar3dcan(ic, kk, MV3D_TPLUS)
            do kc = 1, nkc
               write(6,*) 'mach_canopy_levels: value of zcan(ic kc) z2(ic,chm_nk+1) and difference  at this value of ic for kk: ',&
                            kc,' are: ',zcan(ic,kc),z2(ic,chm_nk+1), zcan(ic,kc)-z2(ic,chm_nk+1)
            end do

            do k = 1, nkt
               write(6,*) 'mach_canopy_levels: value of zthrmcan for ic = ',ic,' at k = ',k,' is: ',zthrmcan(ic,k)
            end do

            do kc = 1, nkc
               write(6,*) 'mach_canopy_levels:  values of kcan zcan and original zcan for ic = ',ic,' at kc = ',kc,' are: ',&
                           kcan(ic, kc), zcan(ic, kc), hc(ic) * can_frac(kc)
            end do

            do k = 1,chm_nk
               write(6,*) 'mach_canopy_levels:  value of kmod and z for ic = ',ic,' at k = ',k,' are: ',kmod(ic,k), zt(ic, k)
            end do

            do kc = 1, nkc
               write(6,*) 'mach_canopy_levels: value of klower_can at this value of ic for kc: ',kc,' is: ',klower_can(ic,kc)
            end do

            chm_error_l = .true.
            return
         end if
      end do
   end do
   end if
!
!  (2) For the last few variables, the value at the lowest
!  resolved model layer and typical profiles for that variable
!  within the canopy will be used to create the canopy values:
   do kc = 1, nkc
      do ic = 1, ni_can
         ii = imod(ic)
         kk = kcan(ic, kc)
!  Ratio of lowest model level to canopy height:
!
         zr = (zt(ic, chm_nk) - z2(ic, chm_nk+1)) / hc(ic)
!
!  Horizontal wind and KT profiles are from Raupach, Quarterly Journal
!  of the Royal Meteorological Society, vol 115, pp 609-632, 1989, examples
!  from page 626, equations (48) through (51).
!
!  Wind speed (equation 51), assumed to scale similarly in each horizontal dimension:
!
!   U(z) = ustar/karman * ln((z - d) / z0), where
!   k = 0.4
!   d = 0.75 hc
!   z0 = 0.07530 hc
! The next few lines calculate the average value of u(z), v(z), Raupach's eqn 51,
! at the first resolved level model height
         uh = metvar2d(ii, MV2D_UE) * 3.0
         if (zr >= 1.0) then
            uspr = metvar2d(ii, MV2D_UE) / karman * &
                   alog((zt(ic, chm_nk) - z2(ic, chm_nk+1) - 0.75 * hc(ic)) / &
                   (0.07530 * hc(ic)))
         else
            uspr = uh * exp(- 2.0 * ( 1.0 - zr))
         end if
!  wndr is the ratio of the wind to Raupach's average us(), eqn 51.
!  This is used to scale the wind speed with height values from eqn 51 to the current grid square
         wndr = metvar3d(ii, chm_nk, MV3D_WS) / uspr
!  Using Raupach's formulae for wind speed, multiplied by the above ratio, for the canopy layers:
!
         zr = (zcan(ic, kc) - z2(ic, chm_nk+1)) / hc(ic)
         if (zr >= 1.0) then
            uspr = alog((zcan(ic, kc) - z2(ic, chm_nk+1) - 0.75 * hc(ic)) / &
                   (0.07530 * hc(ic))) * metvar2d(ii, MV2D_UE)
         else
            uspr = uh * exp(- 2.0 * (1.0 - (zcan(ic, kc) - z2(ic, chm_nk+1)) / hc(ic)))
         end if
         metvar3dcan(ic, kk, MV3D_WS) = wndr * uspr
!
!  Coefficients of diffusivity:
!  Find value of K at first model level from raupach's sigw and TL formulae (eqns 48, 49)
         zr = (zt(ic, chm_nk) - z2(ic, chm_nk+1)) / hc(ic)
!  Gradient in stability under the canopy is reduced for higher stability conditions
!  in accord with Shaw, den Hartog and Neumann, BLM 45, 391-409, 1988, Fig 16.
         zl = hc(ic) * safe_inv_mo_length(ic)
! Unstable:
         if(zl < -0.1) then
             a1 = 0.75
             b1 = 0.5
             c1 = 1.25
         end if
! Neutral:
         if(zl >= -0.1 .and. zl < 0.1) then
             a1 = 0.625
             b1 = 0.375
             c1 = 1.0
         end if
! Stable:
         if(zl >= 0.1 .and. zl < 0.9) then
             rat = 4.375 - 3.75 * zl
             a1 = 0.125 * rat + 0.125
             b1 = 0.125 * rat - 0.125
             c1 = 0.25 * rat
         end if
! Very stable (from extrapolation of Shaw et al's values at 0.1 and 0.5:
         if(zl >= 0.9 .or. metvar3d(ii, chm_nk, MV3D_KT) <= min_kt) then
             a1 = 0.25
             b1 = 0.0
             c1 = 0.25
         end if
!  Raupach's originals:
!         if (zr >= 1.0) then
!            sigw = metvar2d(ii,MV2D_UE) * 1.25
!         else
!            sigw = metvar2d(ii,MV2D_UE) * ( 0.75 + 0.5 * cos(pi * (1.0 - (zt(ic,chm_nk) - z2(ic,chm_nk+1))/hc(ii)) ) )
!         end if
!  Replace Raupach's originals with fit to Patton et al and Shaw et al 1988
         if(zr < 0.175) then
               sigw = metvar2d(ii, MV2D_UE) * 0.25
         else
           if(zr < 1.25) then
               sigw = metvar2d(ii, MV2D_UE) * ( a1 + b1 * cos(pi / 1.06818 * &
                      (1.25 - (zt(ic, chm_nk) - z2(ic, chm_nk+1)) / hc(ic))))
           else
               sigw = metvar2d(ii, MV2D_UE) * c1
           end if
         end if

         tl = hc(ic) / metvar2d(ii, MV2D_UE)  * &
              (0.256 * ((zt(ic, chm_nk) - z2(ic, chm_nk+1) - 0.75 * hc(ic)) / hc(ic)) + &
               0.492 * exp (-(0.256 * ((zt(ic, chm_nk) - z2(ic, chm_nk+1)) / hc(ic)) / 0.492)))
! ktr is the ratio of the resolved model diffusivity at the lowest resolved
! model level to that derived by Raupach's formula
!
         ktr =  metvar3d(ii, chm_nk, MV3D_KT)  / (sigw * sigw * tl)
!
!  Use Raupach's formulae for diffusivity, multiplied by the above ratio, for the canopy layers:
!
         zr = (zcan(ic, kc) - z2(ic, chm_nk+1)) / hc(ic)
!  Gradient in stability under the canopy is reduced for higher stability conditions
!  in accord with Shaw, den Hartog and Neumann, BLM 45, 391-409, 1988, Fig 16.
!  Raupach's original:
!         if (zr >= 1.0) then
!            sigw = metvar2d(ii,MV2D_UE) * 1.25
!         else
!            sigw = metvar2d(ii,MV2D_UE) * ( 0.75 + 0.5 * cos(pi * (1.0 - (zcan(ic,kc) - z2(ic,chm_nk+1))/hc(ii)) ) )
!         end if
         if(zr < 0.175) then
               sigw = metvar2d(ii, MV2D_UE) * 0.25
         else
           if(zr < 1.25) then
               sigw = metvar2d(ii, MV2D_UE) * ( a1 + b1 * cos(pi / 1.06818 * &
                      (1.25 - (zcan(ic, kc) - z2(ic, chm_nk+1))/hc(ic))))
           else
               sigw = metvar2d(ii, MV2D_UE) * c1
           end if
         end if
!
         tl = hc(ic) / metvar2d(ii, MV2D_UE)  *  &
              (0.256 * ( (zcan(ic, kc) - z2(ic, chm_nk+1) - 0.75 * hc(ic)) / hc(ic)) + &
              (0.492 * exp (-(0.256 * (zcan(ic, kc) - z2(ic, chm_nk+1)) / hc(ic)) / 0.492) ) )

         metvar3dcan(ic, kk, MV3D_KT)  = (sigw * sigw * tl) * ktr
      end do
   end do
!
   if (local_dbg) then
      do kc = 1, nkc
         do ic = 1, ni_can
            if (kcan(ic, kc) == 0) then
               write(6,*) 'kcan zero inside mach_canopy_levels at ic kc = ', &
                           ic, kc
               chm_error_l = .true.
               return
            end if
         end do
      end do
   end if
!
   return

end subroutine mach_canopy_levels
