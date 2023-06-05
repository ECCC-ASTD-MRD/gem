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
! Fichier/File   : mach_hetv_hetchem.ftn90
! Creation       : V. Bouchet, S. Menard, P. Makar, S. Gravel, B. Pabla
! Description    : Calling heterogeneous chemistry. Rearrange arrays and units conversion
!
! Extra info     : On input it is assumed that gaseous species and aerosol species
!                  concentrations are in kg/kg-of-air. Concentration units are all
!                  converted back before exiting this subroutine.
!
! Arguments  IN
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_hetv_hetchem(gascon, aerocon, tempk, zpres, aeronum, &
                             rhrow, rhorow, ibulk, jlat, kount, pni, pnk)
   use mach_cam_utils_mod,    only: isize, icom, maxns
!!if_off
   use mach_hetv_mod,         only: maxhet, maxghet, lolimit
   use mach_hetv_headers_mod, only: mach_hetv_hetvcall, mach_hetv_rebin
   use mach_cam_utils_mod,    only: igs_SO4, igs_HNO3, igs_NH3,  &
                                    iae_SU, iae_NI, iae_AM
   use chm_utils_mod,         only: chm_error_l, CHM_MSG_DEBUG
   use chm_nml_mod,           only: chm_hetchem_s, chm_timings_L
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: ibulk
   integer(kind=4), intent   (in) :: jlat
   integer(kind=4), intent   (in) :: kount
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent(inout) :: gascon (pni, pnk, maxns)
   real(kind=4),    intent(inout) :: aerocon(pni, pnk, icom, isize)
   real(kind=4),    intent   (in) :: tempk  (pni, pnk)
   real(kind=4),    intent   (in) :: zpres  (pni, pnk)
   real(kind=4),    intent   (in) :: aeronum(pni, pnk, isize)
   real(kind=4),    intent   (in) :: rhrow  (pni, pnk)
   real(kind=4),    intent   (in) :: rhorow (pni, pnk)
!!if_off
!
!  local variables
!
   integer(kind=4) :: nptsnz
   integer(kind=4) :: i, k, isz, nn, mm
   real(kind=8)    :: orgmass, diff
   real(kind=8)    :: so4_ratio, hno3_ratio, nh3_ratio
   real(kind=8)    :: tem      (pni * pnk)
   real(kind=8)    :: rho      (pni * pnk)
   real(kind=8)    :: rh       (pni * pnk)
!  for heterogeneous chemistry)
   real(kind=8)    :: het      (pni, pnk, maxhet, isize)
   real(kind=8)    :: ghet     (pni, pnk, maxghet)
   real(kind=8)    :: hetnew   (pni * pnk, maxhet, isize)
   real(kind=8)    :: ghetnew  (pni * pnk, maxghet)
   real(kind=8)    :: het_orig (pni, pnk, maxhet)
   real(kind=8)    :: dhet_chem(pni, pnk, maxhet)
   real(kind=8)    :: totav    (pni, pnk, 3)
   real(kind=8)    :: totap    (pni, pnk, 3)
! for debugging
   real(kind=8)    :: het_old  (pni, pnk, maxhet, isize)
   real(kind=8)    :: ghet_old (pni, pnk, maxghet)

!          |                                      | T |           |    |
!  NAME    |          DESCRIPTION                 | Y |DIMENSIONS |IN/ |
!          |                                      | P |           |OUT |
!          |                                      | E |           |    |
!-------------------------------------------------------------------------
! g        | gas/part species conc (ppm)          | R |pni, pnk,  | i/o|
!          |                                      |   |maxns      |    |
! aerocon  | Aqueous sp. conc(molar) in cloud     | R |pni, pnk,  | i/o|
!          | water & Ice/Snow                     |   |icom, isize|    |
! tempk    | Atmospheric Temperature (Kelvin)     | R |pni, pnk   | i  |
! zpres    | Pressure (mb)                        | R |pni, pnk   | i  |
! rhorow   | air density (kg/m3)                  | R |pni, pnk   | i  |
! ibulk    | Switch for bulk chemistry or bin     | I | scalar    | i  |
!          | resolved                             |   |           |    |
!          | ibulk=1(bulk chemistry activated)    |   |           |    |
!          | ibulk=0(chemistry for activated bins)|   |           |    |

!
!  external subroutines
!
   external msg_toall, timing_start_omp, timing_stop_omp

   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_hetv [BEGIN]')
   if (chm_timings_L) call timing_start_omp(360, 'mach_hetv', 340)

   het_orig  = 0.0d0
   dhet_chem = 0.0d0

   so4_ratio   = 96.0636d0 / 98.0795d0
   hno3_ratio  = 62.0049d0 / 63.0128d0
   nh3_ratio   = 18.0385d0 / 17.03056d0

!  Modif: hetchem is now always called from mach_cam_aerocld, and species
!  are coming in as kg/kg het. chem applies to all the grid points at this point
!  If hetchem is calling isorropia, information coming in is in kg/kg
!  (need unit conversion to ug/m3). If hetchem is calling hetv (1), no unit
!  conversion is needed as species will be passed to hetv in kg/kg

   do k = 1, pnk
      do i = 1, pni
         ghet(i, k, 1) = dble(gascon(i, k, igs_SO4))
         ghet(i, k, 2) = dble(gascon(i, k, igs_HNO3))
         ghet(i, k, 3) = dble(gascon(i, k, igs_NH3))
! debugging
         ghet_old(i, k, 1) = ghet(i, k, 1)
         ghet_old(i, k, 2) = ghet(i, k, 2)
         ghet_old(i, k, 3) = ghet(i, k, 3)
!
      end do
   end do

   do isz = 1, isize
      do k = 1, pnk
         do i = 1, pni
            het(i, k, 1, isz) = dble(aerocon(i, k, iae_SU, isz))
            het(i, k, 2, isz) = dble(aerocon(i, k, iae_NI, isz))
            het(i, k, 3, isz) = dble(aerocon(i, k, iae_AM, isz))
! debugging
            het_old(i, k, 1, isz) = het(i, k, 1, isz)
            het_old(i, k, 2, isz) = het(i, k, 2, isz)
            het_old(i, k, 3, isz) = het(i, k, 3, isz)
!
         end do
      end do
   end do

   do k = 1, pnk
      do i = 1, pni
!         if (chm_hetchem_s == 'ISO') then
!   unit conversion from kg/kg to ug/m3 order needs to be changed to match what
!   will be used by isocall regardless of where it comes from
!            rhom(i, k) = dble(rhorow(i, k)) * 1.0d9 !air density (in ug/m^3)
!            ghet(i, k, 1) = ghet(i, k, 1) * rhom(i, k)  !h2so4
!            ghet(i, k, 2) = ghet(i, k, 2) * rhom(i, k)  !hno3
!            ghet(i, k, 3) = ghet(i, k, 3) * rhom(i, k)  !nh3
!            do isz = 1, isize
!               het(i, k, 1, isz) = het(i, k, 1, isz) * rhom(i, k) !so4
!               het(i, k, 2, isz) = het(i, k, 2, isz) * rhom(i, k) !no3
!               het(i, k, 3, isz) = het(i, k, 3, isz) * rhom(i, k) !nh4
!            end do
!         end if
         totav(i, k, 1) = ghet(i, k, 1) * so4_ratio
         totav(i, k, 2) = ghet(i, k, 2) * hno3_ratio
         totav(i, k, 3) = ghet(i, k, 3) * nh3_ratio

         do isz = 1, isize
            totav(i, k, 1) = totav(i, k, 1) + het(i, k, 1, isz)
            totav(i, k, 2) = totav(i, k, 2) + het(i, k, 2, isz)
            totav(i, k, 3) = totav(i, k, 3) + het(i, k, 3, isz)
         end do
      end do
   end do

!  call the original interface to heterogeneous chem.
   if (chm_hetchem_s == 'HETV' .and. ibulk == 1) then
!  hetvcall subroutine replaced, see e-mail from Paul feb 25/03head -54

   !  reshape dimension from (pni, pnk) to (pni * pnk)
   !  similar operation for size-dependant arrays are done later
      nptsnz = pni * pnk
      nn = 0
      do k = 1, pnk
         do i = 1, pni
            nn = nn + 1
            do mm = 1, maxghet
               ghetnew(nn, mm) = ghet(i, k, mm)
            end do

!  Keep track of original bulk conc of het in het_orig
            do mm = 1, maxhet
               do isz = 1, isize
                  het_orig(i, k, mm) = het_orig(i, k, mm) + het(i, k, mm, isz)
                  hetnew(nn, mm, isz) = het(i, k, mm, isz)
               end do
            end do

            tem(nn)  = dble(tempk(i, k))
            rho(nn)  = dble(rhorow(i, k))
            rh(nn)   = dble(rhrow(i, k))
         end do
      end do

      call mach_hetv_hetvcall(nptsnz, ghetnew, hetnew, tem, rh, rho, &
                              jlat, kount)
      if (chm_error_l) return

      !  bring buffer back to original format
      nn = 0
      do k = 1, pnk
         do i = 1, pni
            nn = nn + 1
            do mm = 1, maxghet
               ghet(i, k, mm) = ghetnew(nn, mm)
            end do

      !  compile dhet
            do mm = 1, maxhet
               dhet_chem(i, k, mm) = hetnew(nn, mm, 1) - het_orig(i, k, mm)
            end do
         end do
      end do

      !  Rebinning of heterogeneous chemistry
      call mach_hetv_rebin(het, dhet_chem, tempk, zpres, aeronum, &
                           kount, jlat, pni, pnk)
      if (chm_error_l) return
!
   end if
!
   do k = 1, pnk
      do i = 1, pni
!         if (chm_hetchem_s == 'ISO') then
!            ghet(i, k, 1) = ghet(i, k, 1) / rhom(i, k)  !h2so4
!            ghet(i, k, 2) = ghet(i, k, 2) / rhom(i, k)  !hno3
!            ghet(i, k, 3) = ghet(i, k, 3) / rhom(i, k)  !nh3
!            do isz = 1, isize
!               het(i, k, 1, isz) = het(i, k, 1, isz) / rhom(i, k) !so4
!               het(i, k, 2, isz) = het(i, k, 2, isz) / rhom(i, k) !no3
!               het(i, k, 3, isz) = het(i, k, 3, isz) / rhom(i, k) !nh4
!            end do
!         end if

         totap(i, k, 1) = ghet(i, k, 1) * so4_ratio
         totap(i, k, 2) = ghet(i, k, 2) * hno3_ratio
         totap(i, k, 3) = ghet(i, k, 3) * nh3_ratio
         do isz = 1, isize
            totap(i, k, 1) = totap(i, k, 1) + het(i, k, 1, isz)
            totap(i, k, 2) = totap(i, k, 2) + het(i, k, 2, isz)
            totap(i, k, 3) = totap(i, k, 3) + het(i, k, 3, isz)
         end do
         do mm = 1, 3
            orgmass = totav(i, k, mm) * 1.0d-3
            diff = totap(i, k, mm) - totav(i, k, mm)
            if (diff > orgmass .or. diff < -1.0d0 * orgmass) then
               write(0, *) '### Error in mach_hetv_hetchem (hetv)  ###'
               write(0, *) '# mass balance pb ', i, k, mm, kount, jlat
               write(0, *) '#', totav(i, k, mm), totap(i, k, mm)
               write(0, *) '#', diff, orgmass
! debugging
               write(0, *) '#', 'T, RH : ', tempk(i,k), rhrow(i,k)
               write(0, *) '#', 'before hetv...'
               write(0, *) '#', ghet_old(i,k,1),ghet_old(i,k,2),ghet_old(i,k,3)
               write(0, *) '#', (het_old(i,k,1,isz),isz=1,isize)
               write(0, *) '#', (het_old(i,k,2,isz),isz=1,isize)
               write(0, *) '#', (het_old(i,k,3,isz),isz=1,isize)
               write(0, *) '#', 'het_orig :', het_orig(i,k,1),het_orig(i,k,2) &
                              , het_orig(i,k,3)
               write(0, *) '#', 'after hetv+rebining...'
               write(0, *) '#', ghet(i,k,1),ghet(i,k,2),ghet(i,k,3)
               write(0, *) '#', (het(i,k,1,isz),isz=1,isize)
               write(0, *) '#', (het(i,k,2,isz),isz=1,isize)
               write(0, *) '#', (het(i,k,3,isz),isz=1,isize)
!debuggin end
               write(0, *) '###         ABORT         ###'
               chm_error_l = .true.
               return
            endif
         end do
      end do
   end do

   do k = 1, pnk
      do i = 1, pni
         gascon(i, k, igs_SO4)  = real(ghet(i, k, 1))
         gascon(i, k, igs_HNO3) = real(ghet(i, k, 2))
         gascon(i, k, igs_NH3)  = real(ghet(i, k, 3))
      end do
   end do

   do isz = 1, isize
      do k = 1, pnk
         do i = 1, pni
            aerocon(i, k, iae_SU, isz) = real(het(i, k, 1, isz))
            aerocon(i, k, iae_NI, isz) = real(het(i, k, 2, isz))
            aerocon(i, k, iae_AM, isz) = real(het(i, k, 3, isz))
         end do
      end do
   end do

   call msg_toall(CHM_MSG_DEBUG, 'mach_hetv [END]')
   if (chm_timings_L) call timing_stop_omp(360)
   !-----------------------------------------------------------------

   return
end subroutine mach_hetv_hetchem
