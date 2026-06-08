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

!  This subroutine subdivides the original input bin structure
!  into sub-bins, in order to reduce linear interpolation dispersion
!  and bin resolution problems during the transfer of mass between bins
!  accompanying bin growth.   The routine calls mach_cam_intrsec1.ftn90
!  to do the actual mass transfer between sub-bins.
!
!  6 options are possible, making use of the user input
!  variable "chm_intrsec_ver" on the call list.  The code
!  splits the input particle mass, according to the version number,
!  calculates the corresponding aerosol number in each bin,
!  and then splits the condensation rate by version number.
!
!  The options (chm_intrsec_ver = ) are:
!   1, 2 = mass and delta-mass evenly distributed over sub-bins by radius '
!
!   2, 3 = mass and delta-mass evenly distributed over sub-bins by ln radius'
!
!   3, 6 = mass and delta-mass distributed according to locally scaled trimodal lognormal distribution'
!
!
!  The algorithm is similar to mach_cam_intrsec_outer; the main difference
!  is that here, the input xrow contains the aerosol mass distribution
!  after mass has been added due to aqueous phase or inorganic heterogeneous
!  chemistry processes, the total amount of mass added being daqshem.
!  both xrow and daqchem are split using the same choice of algorithms.
!  The 6 version numbers are to provide continuity and a similar method
!  of mass splitting as mach_cam_intrsec_outer, which also has 2 options
!  for redistributing the aerosol condensation rate.  The latter is unnecessary
!  here, but the version numbers are preserved so that the same method
!  for mass splitting is used in both routines.
!
!  P.A. Makar, W. Gong, S. Gong, February 2009
!
! Arguments:  IN
!               totmas   -> Total mass of aerosol in each bin
!               nn       -> Loop index on aerosol types
!               rhopd    -> dry aerosol density [ug / m^3]
!               aeronum  -> Number concentration of aerosols
!               aeronum_n-> Number concentration of aerosols in sub bins
!               xrow_n   -> Tracer concentration with sub-bins for the particles
!               q_bin    -> liquid water content in each size bin; needed for
!                           distributing daqchm in the call following
!                           mach_incld_main.ftn90
!               rcrit    -> critical radius expressed in real bin index
!                           (including fraction)
!               iswitch  -> 1 for the call following mach_incld_main
!                           0 otherwise
!
!             IN / OUT
!              XROW      -> Tracer concentration in each bin before / after intersection tranport
!
!             MODULE: CAM_UTILS
!               ntr      -> Total number of trace substances (gases and aerosols)
! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
!
!!if_on
 subroutine mach_cam_intrsec1_outer(xrow, rhopd, daqchm, aeronum, q_bin, &
                                    rcrit, iswitch, pni, pnk)
   use mach_cam_utils_mod,   only: isize, ntr
!!if_off
   use chm_nml_mod,          only: chm_intrsec_ver, chm_intrsec_ndiv
   use mach_cam_headers_mod, only: mach_cam_intrsec1_inner
   use mach_cam_utils_mod,   only: nbnd, iae1, icom, nsb, sub_pvol, dvn
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: iswitch
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent(inout) :: xrow   (pni, pnk, ntr)
   real(kind=4),    intent   (in) :: rhopd  (pni, pnk, isize)
   real(kind=4),    intent   (in) :: daqchm (pni, pnk, isize)
   real(kind=4),    intent   (in) :: aeronum(pni, pnk, isize)
   real(kind=4),    intent   (in) :: q_bin  (pni, pnk, isize)
   real(kind=4),    intent   (in) :: rcrit  (pni, pnk)
!!if_off
!
!  local variables
!
   integer(kind=4)                              :: k, np, pos, il, i
!  New variables:
   real(kind=4), dimension(pni, pnk, nbnd) :: aeronum_n, daqchm_n, rmass_n
   real(kind=4), dimension(pni, pnk, nsb)  :: xrow_n
   real(kind=4), dimension(pni, pnk, nbnd) :: q_bin_n
   real(kind=4)    :: net_vol
   real(kind=4)    :: aeronum_t, aeronum_t1, aeronum_t2, aeronum_tsub, dlt
   integer(kind=4) :: j, jsub, jk
!
!
   if (isize == 12) then
!  skip sub-bin division when the 12-bin structure is used
!
      do k = 1, isize
         do il = 1, pnk
            do i = 1, pni
               daqchm_n(i, il, k) = daqchm(i, il, k)
               rmass_n(i, il, k)  = sub_pvol(k) * rhopd(i, il, k)
               aeronum_n(i, il, k) = aeronum(i, il, k)
!
               do np = 1, icom
                  pos = (np - 1) * isize + k + iae1 - 1
                  jk = (np - 1) * isize + k
                  xrow_n(i, il, jk) = xrow(i, il, pos)
               end do
            end do
         end do
      end do
!
   else
!
!  Commence subdivision of original bins:
!
!  Particle density is assumed to be the same within the sub-bin
!  as within the original bin:
!
      do k = 1, isize
         do j = 1, chm_intrsec_ndiv
            jsub = (k - 1) * chm_intrsec_ndiv + j
            do il = 1, pnk
               do i = 1, pni
!  Particle density is assumed to be the same within the sub-bin
!  as within the original bin:
                  rmass_n(i, il, jsub) = sub_pvol(jsub) * rhopd(i, il, k)
! Net particle mass change:
                  daqchm_n(i, il, jsub) = daqchm(i, il, k) * dvn(jsub)
!
!  Assign new particle mass, and new particle number based on individual 
!  particle volume, density,:
                  net_vol = 0.0
                  do np = 1, icom
                     pos = (np - 1) * isize + k + iae1 - 1
                     jk = (np - 1) * nbnd + jsub
                     xrow_n(i, il, jk) = xrow(i, il, pos) * dvn(jsub)
                     net_vol = net_vol + xrow_n(i, il, jk) / rhopd(i, il, k)
                  end do
                  aeronum_n(i,il,jsub)  = net_vol / sub_pvol(jsub)
               end do
            end do
         end do
      end do
!
!
!  Distributing bulk paritcle mass change to sub-bins for the call after
!  in-cloud chemistry (weighted by liquid water content in sun-bins which is
!  obtained by distributing q_bin equally amongst activated particles in sub-
!  divided bins) - W. Gong
!
      if (iswitch == 1) then
!
         q_bin_n = 0.0
         daqchm_n = 0.0
!
!   q_bin contains the liquid water content in each bin before subdivision
!   which is being distributed to activated sub-bins to obtain q_bin_n in
!   the following. Q_bin_n is then used for distributing daqchm to daqchm_n
!
         do il = 1, pnk
            do i = 1, pni
               do k = int(rcrit(i, il)), isize
                  aeronum_t = 0.0
                  do j = 1, chm_intrsec_ndiv
                     jsub = (k - 1) * chm_intrsec_ndiv + j
                     aeronum_t = aeronum_t + aeronum_n(i, il, jsub)
                  end do
                  if (aeronum_t > 0.0 .and. q_bin(i, il, k) > 0.0) then
                     if (k > int(rcrit(i, il))) then
                        do j = 1, chm_intrsec_ndiv
                           jsub = (k - 1) * chm_intrsec_ndiv + j
                           q_bin_n(i, il, jsub) = q_bin(i, il, k)   &
                                        * aeronum_n(i, il, jsub) / aeronum_t
                           daqchm_n(i, il, jsub) = daqchm(i, il, k) &
                                        * q_bin_n(i, il, jsub) / q_bin(i, il, k)
                        end do
                     else
                        aeronum_tsub = aeronum_t * (rcrit(i, il) - int(rcrit(i, il)))
                        aeronum_t1 = 0.0
                        aeronum_t2 = 0.0
                        dlt = 1.0
                        do j = 1, chm_intrsec_ndiv
                           jsub = (k - 1) * chm_intrsec_ndiv + j
                           aeronum_t1 = aeronum_t1 + aeronum_n(i, il, jsub) * dlt
                           if (aeronum_t1 <=  aeronum_tsub) then
                              aeronum_t2 = aeronum_t1
                           else
                              q_bin_n(i, il, jsub) = q_bin(i, il, k)   &
                                        * aeronum_n(i, il, jsub) / (aeronum_t - aeronum_t2)
                              daqchm_n(i, il, jsub) = daqchm(i, il, k) &
                                        * q_bin_n(i, il, jsub) / q_bin(i, il, k)
                              dlt = 0.0
                           end if
                        end do
                     end if
                  end if
               end do
            end do
         end do
!
      end if
!
!  Finished subdivision of inputs for condensation.
  
   end if
!
!  Call the original intrsec code to determine the mass transfer:
!
   call mach_cam_intrsec1_inner(xrow_n, daqchm_n, aeronum_n, rmass_n, pni, pnk)
!  Reassign the mass back into the original xrow locations:

   if (isize == 12) then
      do k = 1, isize
         do np = 1, icom
            pos = (np - 1) * isize + k + iae1 - 1
            jk = (np - 1) * isize + k
            do il = 1, pnk
               do i = 1, pni
                  xrow(i, il, pos) = xrow_n(i, il, jk)
               end do
            end do
         end do
      end do
   else
      do k = 1, isize
         do np = 1, icom
            pos = (np - 1) * isize + k + iae1 - 1
            do il = 1, pnk
               do i = 1, pni
                  xrow(i, il, pos) = 0.0
                  do j = 1, chm_intrsec_ndiv
                     jsub = (k - 1) * chm_intrsec_ndiv + j
                     jk = (np - 1) * nbnd + jsub
                     xrow(i, il, pos) = xrow(i, il, pos) + xrow_n(i, il, jk)
                  end do
               end do
            end do
         end do
      end do
   endif
!
   return
 end subroutine mach_cam_intrsec1_outer