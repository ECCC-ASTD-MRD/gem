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
! Fichier/File   : mach_canopy_transfer.ftn90
! Creation       : Paul Makar - 2015
! Description    : Transfers resolved scale model layers mass to canopy layers
!                   and vice versa.
!                  
! Extra info     : 
!
! Arguments:
! Input variables
!-----------------------------------------------------------------------------
! Array dimensions:
!  nkc                                :  number of canopy levels
!  nkt                                :  number of levels in gathered canopy
!                                        + resolved scale columns = chm_nk+nkc
!  ni_can                             :  number of columns containing a forest canopy
!  metvar2d(chm_ni, SIZE_MV2D)        :  met 2d variables, (full slice)
!  metvar3d(chm_ni, chm_nk, SIZE_M32D):  met 3d variables, (full slab)
!  metvar3dcan(ni_can, nkt, SIZE_MV3D):  met 3d variables, gathered canopy + resolved scale columns
!  kmod(ni_can,chm_nk)                :  Vertical index location of original ungathered model layer in combined
!                                        canopy + resolved scale column  
!  kcan(ni_can,nkc)                   :  Vertical index location of canopy ungathered model layer in combined 
!                                          canopy + resolved scale column
!  flag                               : 0 -> resolved_to_canopy
!                                       1 -> canopy_to_resolved
!
!  Input/Output variables, original horizontal coordinate
!  tracers_can(ni_can, nkc, nb_tracers):  Chemical tracers concentrations within canopy
!  chem_tr(chm_ni, chm_nk+1, nb_tracers)    :  Chemical tracers concentrations on model levels 
!
!  Local variables:
!  massaircan(ni_can, nkt)    :  mass of air in canopy layers (kg)
!  massairmod(ni_can, chm_nk) :  mass of air in model layers (kg)
!                                        (gathered canopy + resolved scale columns)
!   nfrct(nkt,ni_can)        :  Number of original model levels contributing to canopy level k
!   ifrct(nkt,2,ni_can)      :  Index of the original model level contributing to canopy level k
!   frctr2c(nkt,2,ni_can)    :  Fractional contribution of the original model level to canopy level k
!   frctc2r(nkt,2,ni_can)    :  Fractional contribution of the canopy level to the original model level
!
!=============================================================================
!
!!if_on
subroutine mach_canopy_transfer(chem_tr, tracers_can, metvar2d, metvar3d, &
                                metvar3dcan, kmod, kcan, imod, ni_can, flag)
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk, nkt, nkc
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D 
   use chm_species_info_mod, only: nb_dyn_tracers
!!if_off
   use chm_utils_mod,        only: ik, chm_error_l, global_debug, chm_lun_out
   use chm_metvar_mod,       only: MV2D_DXDY, MV3D_ZMOM, MV3D_ZPLUS, MV3D_RHO
   implicit none
!!if_on
   integer(kind=4),   intent(in)    :: ni_can
   real(kind=4),      intent(inout) :: chem_tr(chm_ni, chm_nk+1, nb_dyn_tracers)
   real(kind=4),      intent(inout) :: tracers_can(ni_can, nkc, nb_dyn_tracers)
   real(kind=4),      intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),      intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
   real(kind=4),      intent   (in) :: metvar3dcan(ni_can, nkt, SIZE_MV3D)
   integer(kind=4),   intent   (in) :: imod(ni_can)
   integer(kind=4),   intent   (in) :: kmod(ni_can, chm_nk)
   integer(kind=4),   intent   (in) :: kcan(ni_can, nkc)
   integer(kind=4),   intent   (in) :: flag
!!if_off
!
! Local variables:
   integer(kind=4) :: k, kk, ic, ii, isp, kc, k2
   real(kind=4), dimension(ni_can, nkt)       :: massaircan
   real(kind=4), dimension(ni_can, chm_nk)    :: massairmod
   real(kind=4), dimension(nkt)               :: mmr_canopy
   real(kind=4), dimension(chm_nk)            :: mass_resolved
   real(kind=4), dimension(nkt)               :: masscan, conc
   real(kind=4), dimension(ni_can, nkt+1)     :: zmomcan
   real(kind=4), dimension(ni_can, chm_nk+1)  :: zmom
   integer(kind=4), dimension(nkt, ni_can)    :: nfrct
   integer(kind=4), dimension(nkt, 2, ni_can) :: ifrct
   real(kind=4),    dimension(nkt, 2, ni_can) :: frctr2c
   real(kind=4),    dimension(nkt, 2, ni_can) :: frctc2r

! Diagnostic height is the assumed height above ground of the
! sampling for observations
   real(kind=4),    parameter              :: diag_hgt = 2.0

   logical(kind=4)                         :: local_dbg

   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))
!
   do ic = 1, ni_can
      ii = imod(ic)
! Calculate mass of air in canopy levels
      zmomcan(ic, nkt + 1) = 0.0
      do k = nkt, 1, -1
         zmomcan(ic, k)    = metvar3dcan(ic, k, MV3D_ZMOM)
         massaircan(ic, k) = metvar3dcan(ic, k, MV3D_RHO) * metvar2d(ii, MV2D_DXDY) * &
                             (zmomcan(ic, k) - zmomcan(ic, k + 1))
      end do
!  Calculate mass of air in model levels
      zmom(ic, chm_nk + 1) = 0.0
      do k = chm_nk, 1, -1
         zmom(ic, k)       = metvar3d(ii, k, MV3D_ZMOM)
         massairmod(ic, k) = metvar3d(ii, k, MV3D_RHO) * metvar2d(ii, MV2D_DXDY) * &
                             (zmom(ic, k) - zmom(ic, k + 1))
      end do
   end do
!
!  Next, we need a set of arrays which track mass transfer from resolved to model layers; 
!  how much of the original (aka "resolved") model layer mass goes into each canopy layer,
!  given the above level structure.  The three arrays are:
!  nfrct(k,ic) : the number of resolved model levels contributing to canopy level k
!  ifrct(k,n,ic) : the index of the resolved model level contributing to canopy level k (n is at most 2)
!  frctr2c(k,n,ic)  : the fractional contribution of the resolved model level to canopy level k
!  frctc2r(k,n,ic)  : the fractional contribution of the canopy model level to the resolved model level
!
!  Check for coincident layers first:
!
   do ic = 1, ni_can
      inner: do k = 1, chm_nk
!  If the following IF statement is true, then the canopy and resolved 
!  model layer upper and lower boundaries coincide, and the entire resolved model
!  model layer contributes to the combined model layer (trivial case).
         if (zmomcan(ic, k) == zmom(ic, k) .and. zmomcan(ic, k+1) == zmom(ic, k+1)) then
            nfrct(k, ic) = 1  
            ifrct(k, 1, ic) = k    
            frctr2c(k, 1, ic) = 1.0   
            frctc2r(k, 1, ic) = 1.0 
         else
            exit inner
         end if
      end do inner
!
!  "k" is the first layer where boundary levels do not match on output from the above loops.
! Determine fractions of original model layer structure contributing to canopy model layers.
      k2 = k
      do k = k2, nkt
         do kk = k2, chm_nk
! (1) Upper boundaries of combined and resolved model layers coincide,
! lower boundary of combined layer is within resolved layer, so canopy
! layer resides entirely within resolved layer, and shares an upper boundary
! with the resolved layer:
            if ((zmomcan(ic, k) == zmom(ic, kk) .and. zmomcan(ic, k+1) > zmom(ic, kk+1)) .or. &
! (2) Lower boundaries coincide, upper boundary of combined layer is within resolved layer,
! so canopy layer resides entirely within the resolved layer, and shares a lower boundary
! with the canopy layer.
                (zmomcan(ic, k+1) == zmom(ic, kk+1) .and. zmomcan(ic, k) > zmom(ic, kk)) .or. &
! (3) Both canopy layer boundaries exist inside a resolved layer, with no shared boundaries:
                (zmomcan(ic, k) < zmom(ic, kk) .and. zmomcan(ic, k+1) >= zmom(ic, kk+1))) then
               nfrct(k, ic) = 1   
               ifrct(k, 1, ic) = kk  
               frctr2c(k, 1, ic) = (zmomcan(ic, k) - zmomcan(ic, k+1)) / (zmom(ic, kk) - zmom(ic, kk+1))
               frctc2r(k, 1, ic) = 1.0  ! canopy layer resides within resolved model layer
            end if
!  Resolved layer boundary splits a combined canopy layer:
!  This case arises if, due to the use of the momentum levels in the canopy column
!  sometimes being half-way between the thermodynamic levels, a resolved model
!  momentum layer falls within the canopy layer.  Since the resolved model layers are 
!  defacto thicker than the canopy layers, this means that there can at most be two 
!  resolved model layers contributing to the canopy layer (only case where nfrct = 2).  
            if (zmomcan(ic, k+1) < zmom(ic, kk) .and. zmomcan(ic, k) > zmom(ic, kk)) then
               nfrct(k, ic) = 2
               ifrct(k, 1, ic) = kk
               ifrct(k, 2, ic) = kk-1
!  Fraction of resolved model layer contributing to canopy layer:
               frctr2c(k, 1, ic) = (zmom(ic, kk) - zmomcan(ic, k+1)) / (zmom(ic, kk) - zmom(ic, kk+1))
               frctr2c(k, 2, ic) = (zmomcan(ic, k) - zmom(ic, kk)) / (zmom(ic, kk-1) - zmom(ic, kk))
!  Fraction of canopy layer contributing to resolved model layer:
               frctc2r(k, 1, ic) = (zmom(ic, kk) - zmomcan(ic, k+1)) / (zmomcan(ic, k) - zmomcan(ic, k+1))
               frctc2r(k, 2, ic) = (zmomcan(ic, k) - zmom(ic, kk)) / (zmomcan(ic, k) - zmomcan(ic, k+1))
            end if
         end do
      end do
!      
   end do
    
!
!  massaircan thus contains the mass of air in the canopy layers in kg, while massairmod contains the
!  mass of air in the original model layers, at the canopy columns ic.
!
!  return tracers to resolved scale model layers:
!
   if (flag == 1) then  ! "canopy_to_resolved"
      do isp = 1, nb_dyn_tracers
!
!  At this point, the model mass is distributed over the combined layers,
!  and the tracer concentration arrays are both in the combined layer system.  
!
!  (1) We start off by converting these mass mixing ratios to mass in ug:
         do ic = 1, ni_can
            ii = imod(ic)
            do k = 1, chm_nk
               kc = kmod(ic, k)
               conc(kc) = chem_tr(ii, k, isp)
            end do
            do k = 1, nkc
               kc = kcan(ic, k)
               conc(kc) = tracers_can(ic, k, isp)
            end do
! 
! (2) Array "masscan" now holds the mass of the tracer in each of the combined levels.  
! This mass must be added back to the resolved levels:
            mass_resolved = 0.
            do k = 1, nkt
               masscan(k) = conc(k) * massaircan(ic, k)
               do kk = 1, nfrct(k, ic)
                  kc = ifrct(k, kk, ic)
                  mass_resolved(kc) = mass_resolved(kc) + masscan(k) * frctc2r(k,kk,ic)
               end do
            end do
!
!  Check:  total mass in the column should be the same
            if (local_dbg) then
               call canopy_mass_check(masscan, mass_resolved, ic, flag)
               if (chm_error_l) return
            end if
!
!  (3) The masses in ug need to be converted back to ug/kg
            do k = 1, chm_nk
!               
               chem_tr(ii, k, isp) = mass_resolved(k) / massairmod(ic, k)
            end do
!
!  (4) Evaluate the diagnostic level concentration
!  Find the bounding layers above and below the diagnostic height:
!  kk'th layer is the layer above the inlet height
            kk = nkt
            do k = nkt, nkt-8, -1
               if (diag_hgt <= metvar3dcan(ic, k-1, MV3D_ZPLUS) .and. &
                  diag_hgt > metvar3dcan(ic, k, MV3D_ZPLUS)) then
                  kk = k - 1
               end if
            end do
!  If the diagnostic height is less than the lowest level, then use that level 
!  for the concentration.
            if (kk == nkt) then
               chem_tr(ii, chm_nk + 1, isp) = conc(nkt)
            else
               chem_tr(ii, chm_nk + 1, isp) = conc(kk) + (conc(kk) - conc(kk + 1)) / &
                                       (metvar3dcan(ic, kk, MV3D_ZPLUS) -     &
                                        metvar3dcan(ic, kk+1, MV3D_ZPLUS)) *  &
                                 (diag_hgt - metvar3dcan(ic, kk+1, MV3D_ZPLUS))
            end if
!
         end do  ! loop over ni_can
      end do !species index loop
!  Done transfering from combined canopy + resolved scale back to resolved scale.  :)
!
   else !(can_transfer == "resolved_to_canopy") then
!
      do isp = 1, nb_dyn_tracers
!
!  (1) Convert the original model domain values in the current column to mass from mass mixing ratio:
!  mass_resolved = Mass mixing ratio * (density) / (volume of original model layer)  (ug)
         do ic = 1, ni_can
            ii = imod(ic)
            do k = 1, chm_nk
               mass_resolved(k) = chem_tr(ii, k, isp) * massairmod(ic, k)
            end do
!
!  (2) Use the array fractions defined earlier to divide the resolved layer masses into the canopy layers,
!  and convert back to mixing ratios.  Note that the frctr2c fractions are vertical extent of the
!  contribution of the resolved layer into the canopy layer, hence the mass/volume can be divided up
!  this way:
!  mmr_canopy = sum of masses contributed / (density * volume of canopy model layeri)
            mmr_canopy = 0.
            masscan = 0.
            do k = 1, nkt
               do kk = 1, nfrct(k, ic)
                  kc = ifrct(k, kk, ic)
                  masscan(k) = masscan(k) + mass_resolved(kc) * frctr2c(k, kk, ic)
               end do
            end do
!
!  Check:  total mass in the column should be the same
            if (local_dbg) then
               call canopy_mass_check(masscan, mass_resolved, ic, flag)
               if (chm_error_l) return
            end if
!
            do k = 1, nkt
               mmr_canopy(k) = masscan(k) / massaircan(ic, k)
            end do
!
!  (3) Replace the original model layer values with the corresponding canopy layer values, when
!  a canopy exists:
            do k = 1, chm_nk
               chem_tr(ii, k, isp) = mmr_canopy(kmod(ic, k))
            end do
!
!  (4) Fill the canopy layers with the new mass mixing ratios
            do k = 1, nkc
               tracers_can(ic, k, isp) = mmr_canopy(kcan(ic, k))
            end do
!
         end do  !  loop over canopy columns
!
      end do !species index loop
! 
   end if  ! transfer resolved scale to canopy layers

   return
   
 contains
   
   subroutine canopy_mass_check(mass_can, mass_model, ic, flag)
      implicit none
      integer(kind=4),   intent(in) :: flag, ic
      real(kind=4),      intent(in) :: mass_can(nkt), mass_model(chm_nk)

      character(len=18) :: mode_transfer
      real(kind=4) :: masstotcan, masstotres, massrat
      real(kind=4) :: sum2can(nkt), sum2res(nkt)
      
      masstotcan = 0.
      masstotres = 0.
      do k = 1, nkt
         masstotcan = masstotcan + mass_can(k)
      end do
      do k = 1,chm_nk
         masstotres = masstotres + mass_model(k)
      end do

      if (flag == 1) then  
         mode_transfer = "canopy_to_resolved"
      else
         mode_transfer = "resolved_to_canopy"
      end if

      if (masstotres > 0.0) then
         massrat = masstotcan / masstotres
         if (massrat > 1.001 .or. massrat < 0.999) then
            write(*, *) 'Conversion of mass in mach_canopy_transfer not conserved'
            write(*, *) 'during ', mode_transfer, 'evaluation.  Stopping '
            write(*, *) 'code with masstotcan = ',masstotcan,' and masstotres = ', & 
                              masstotres
            write(*, *) 'Values of mass_can: ',(mass_can(k), k=1, nkt)
            write(*, *) 'Values of mass_resolved: ',(mass_model(k), k=1, chm_nk)
            do k = 1, nkt
               write(*, *) 'canopy layer ',k,'has ',nfrct(k, ic),' contributions'
               do kk = 1, nfrct(k, ic)
                  write(*, *) 'Resolved # ',ifrct(k,kk,ic),' with mass: ',mass_model(ifrct(k,kk,ic)),&
                             ' contributes ',frctr2c(k,kk,ic),' to canopy layer ',k,&
                             ' with masscan ',mass_can(k)
               end do
            end do

            chm_error_l = .true.
            return
         end if
      end if
!
!  Check on the values of the fractions:  they should sum to unity across the number 
!  of original model levels!
!          
      sum2can = 0.
      sum2res = 0.
      do k = nkt, 1, -1
         do kk = 1, nfrct(k, ic)
            kc = ifrct(k, kk, ic)
            sum2can(kc) = sum2can(kc) + frctr2c(k, kk, ic)
            sum2res(k) = sum2res(k) + frctc2r(k, kk, ic)
         end do
      end do
      do k = chm_nk, 1, -1
         if (sum2can(k) < 0.999 .or. sum2can(k) > 1.001) then
            write(*, *) 'layer mismatch in canopy level setup in resolved to canopy indexing'
            write(*, 20) 'sum of non-zero contributions from column ',ic, &
               ' layer ',k,' is ',sum2can(k),' (should be unity).'
            chm_error_l = .true.
            return
         end if
      end do
      do k = nkt, 1, -1
         if (sum2res(k) < 0.999 .or. sum2res(k) > 1.001) then
            write(*, *) 'layer mismatch in canopy level setup in canopy to resolved indexing'
            write(*, 20) 'sum of non-zero contributions from column ',ic, &
               ' layer ',k,' is ',sum2res(k),' (should be unity).'
            write(*, *) 'k nfrct(k ic) frctc2r'
            write(*, *) k, nfrct(k,ic),(frctc2r(k,kk,ic), kk = 1,nfrct(k,ic))
            chm_error_l = .true.
            return
         end if
      end do
      
 20   format(a42, i6, a7, i3, a5, 1pe10.3, a18)
!
      return
   end subroutine canopy_mass_check
   
end subroutine mach_canopy_transfer
