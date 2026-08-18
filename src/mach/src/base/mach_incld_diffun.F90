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
! Fichier/File   : mach_incld_diffun.ftn90
! Creation       : S. Menard,  GEM-MACH,  June 2008.
!:
! Description    : Compute mass transfer / chemistry derivatives in the form
!                  (dc / dt) = cd = p - d * c
!
! Extra info     : Based on ADOMIIB,  level: 06/09/89,  ENSR(AES)
!
! Arguments:
!           IN
!              gaz_conc  --> Gas/part species conc (ppm)
!              baq       --> Aqueous variable coefficients
!              tempk     --> Atmospheric Temperature (Kelvin)
!              tempi     --> Inverse of tempk (1.0 / tempk)
!              psacw     --> CW to snow collection rates (gm / m3 s):
!              ideriv
!           OUT
!              RAQ       --> Rate constants for aqueous-phase
!              PG
!              DG
!              CDG
!              AQ        --> Aqueous sp. conc(molar) in cloud
!              PAQ
!              DAQ
!              CDAQ
!              Q         --> Liquid water concentrations of cloudwater
!
!=============================================================================
!
!!if_on
subroutine mach_incld_diffun(gaz_conc, pg, dg, cdg, aq, paq, daq, cdaq, q, &
                             ideriv, baq, tempk, tempi, raq, psacw)
   use mach_cam_utils_mod,     only: maxnsg, maxnsaq
!!if_off
   use chm_utils_mod,          only: chm_lun_out
   use mach_incld_headers_mod, only: mach_incld_upaqr, mach_incld_rates, &
                                     mach_incld_intrqf
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: ideriv
   real(kind=4),    intent   (in) :: q(2)
   real(kind=4),    intent   (in) :: baq(5, 2)
   real(kind=4),    intent(inout) :: raq(25, 2)
   real(kind=4),    intent   (in) :: gaz_conc(maxnsg)
   real(kind=4),    intent  (out) :: pg      (maxnsg)
   real(kind=4),    intent  (out) :: dg      (maxnsg)
   real(kind=4),    intent  (out) :: cdg     (maxnsg)
   real(kind=4),    intent   (in) :: tempk
   real(kind=4),    intent   (in) :: tempi
   real(kind=4),    intent   (in) :: psacw
   real(kind=4),    intent  (out) :: paq     (maxnsaq, 2)
   real(kind=4),    intent  (out) :: daq     (maxnsaq, 2)
   real(kind=4),    intent  (out) :: cdaq    (maxnsaq, 2)
   real(kind=4),    intent   (in) :: aq      (maxnsaq, 2)
!!if_off
!
!  Local variables
!
   integer(kind=4) :: izero = -1, iprt = -1
   integer(kind=4) :: naqend, i, j
   real(kind=4)    :: zero = 0.e0
   real(kind=4)    :: ppaq(maxnsaq, 2), ddaq(maxnsaq, 2), ccdaq(maxnsaq, 2)
!
!  External subroutines
!
   external msg_toall, timing_start_omp, timing_stop_omp

   naqend = maxnsaq
   if (ideriv < 0) then
      naqend = 3
   end if

!  zero out

   if (izero >= 0) then
      do i = 1, maxnsg
         pg(i) = zero
         dg(i) = zero
         cdg(i) = zero
      end do
      do j = 1, 2
         if (q(j) <= 0.) cycle
         do i = 1, maxnsaq
            paq(i, j) = zero
            daq(i, j) = zero
            cdaq(i, j) = zero
         end do
      end do
   end if

!  Update rate constants that depend on [c] if ideriv > 0

   if (ideriv >= 0) then
      if (q(1) > 0.0) call mach_incld_upaqr(aq, tempi, raq)
   end if

! Compute derivatives for chemical kinetics and
! mass transfer (excluding intra-hydrometeor fluxes)

   call mach_incld_rates(gaz_conc, dg, pg, cdg, aq, daq, paq, cdaq, raq, baq, &
                         q, ideriv)

! Compute intra-hydrometeor flux terms - only riming terms
! computed in this version.

   if (q(1) > 0 .and. q(2) > 0 .and. psacw > 0.) then
      call mach_incld_intrqf(aq, ppaq, ddaq, ccdaq, q, psacw, ideriv, tempk)

! Combine terms

      do j = 1, 2
         do i = 1, naqend
            paq(i, j)  = paq(i, j)  + ppaq(i, j)
            daq(i, j)  = daq(i, j)  + ddaq(i, j)
            cdaq(i, j) = cdaq(i, j) + ccdaq(i, j)
         end do
      end do
   end if
!
   if (iprt < 0) return
   if (chm_lun_out < 0) return
   write (chm_lun_out, 290)(gaz_conc(i), i = 1, maxnsg)
   write (chm_lun_out, 300)(pg(i), i = 1, maxnsg)
   write (chm_lun_out, 301)(dg(i), i = 1, maxnsg)
   write (chm_lun_out, 302)(cdg(i), i = 1, maxnsg)
   do j = 1, 2
      if (q(j) <= 0.) cycle
      write (chm_lun_out, 306)j
      write (chm_lun_out, 295)(aq(i, j), i = 1, maxnsaq)
      write (chm_lun_out, 303)(paq(i, j), i = 1, maxnsaq)
      write (chm_lun_out, 304)(daq(i, j), i = 1, maxnsaq)
      write (chm_lun_out, 305)(cdaq(i, j), i = 1, maxnsaq)
   end do
return
290   format ('   g    ', 12e10.3)
295   format ('   aq   ', 12e10.3)
300   format ('   pg   ', 12e10.3)
301   format ('   dg   ', 12e10.3)
302   format ('   cdg: ', 12e10.3)
303   format ('   paq: ', 12e10.3)
304   format ('   daq: ', 12e10.3)
305   format ('   cdaq:', 12e10.3)
307   format (' diffun gas/part derivative terms --  level =', i3)
306   format (' diffun aqueous derivative terms hydrometeor class =', i2)
end subroutine mach_incld_diffun
