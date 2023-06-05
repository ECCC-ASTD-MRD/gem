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
! Fichier/File   : mach_incld_dochem.ftn90
! Creation       : B. Pabla, W. Gong, S. Menard, S. Gravel, GEM-MACH, June 2008.
!
! Description    : This routine solves the chemistry equations. For increased efficiency,
!                  the slowly-varying species are integrated with larger time steps than
!                  the rapidly-varying species.  This modified version also has the ability
!                  to integrate the chemistry analytically, if no cloudwater is present.
!
! Extra info     : Based on ADOMIIB,  level: 06/09/89,  ENSR(AES)
!
! Arguments:
!           IN
!             qvec       --> liquid water concentrations of cloudwater, ice/snow, & rainwater in air (gm-w/m**3-air)
!             tin        --> initial time (seconds)
!             tout       --> final time (seconds)
!             idrf       --> = 1 on first call from wetchem for each level
!                            = 2 on subsequent calls
!             bvec       -->
!             tempkvec   --> Atmospheric Temperature (Kelvin
!             tempivec   --> Inverse of tempk (1.0/tempk
!             raqvec     --> Rate constants for aqueous-phase
!             psacwvec   --> CW to snow collection rates (gm/m3 s)
!             nptsnz     --> Total number of grids to integrate
!
!           OUT
!           INOUT
!             Gvec        --> gas/part species concentrations (ppm)
!             AQvec       --> aqueous species concentrations (molar) in cloudwater & ice/snow
!=============================================================================
!
!!if_on
subroutine mach_incld_dochem(gvec, aqvec, qvec, tempkvec, tempivec, bvec, &
                             raqvec, psacwvec, tin, tout, idrf, nptsnz)
   use mach_cam_utils_mod,     only: maxnsg, maxnsaq
!!if_off
   use chm_utils_mod,          only: chm_error_l
   use mach_incld_headers_mod, only: mach_incld_fn_dtnew, mach_incld_intsca_il,&
                                     mach_incld_diffun, mach_incld_steady
   use mach_incld_mod,         only: dtsrt, dtmax, gmin, aqmin
   implicit none
!!if_on

   integer(kind=4), intent   (in)  :: idrf
   integer(kind=4), intent   (in)  :: nptsnz
   real(kind=4),    intent   (in)  :: tin, tout
   real(kind=4),    intent   (in)  :: qvec     (nptsnz, 2)
   real(kind=4),    intent   (in)  :: bvec     (nptsnz, 5, 2)
   real(kind=4),    intent   (in)  :: tempkvec (nptsnz)
   real(kind=4),    intent   (in)  :: tempivec (nptsnz)
   real(kind=4),    intent   (in)  :: raqvec   (nptsnz, 25, 2)
   real(kind=4),    intent   (in)  :: psacwvec (nptsnz)
   real(kind=4),   intent (inout)  :: gvec     (nptsnz, maxnsg)
   real(kind=4),   intent (inout)  :: aqvec    (nptsnz, maxnsaq, 2)
!!if_off
!
!  local variables
!
   integer(kind=4)            :: ii, nstep
   integer(kind=4), parameter :: nstmax = 500
   integer(kind=4)            :: igc, imaxnsg, imaxns2, imaxnsaq, iqcntr, &
                                 ibcntr1, ibcntr2, iraqc1, iraqc2
   integer(kind=4)            :: ndt, i, kiter
   integer(kind=4), parameter :: minus1 = -1, plus1 = 1
   integer(kind=4), parameter :: ng1 = 3, naq1 = 4, ng2 = 10, naq2 = 7, ng1p1 = 4

   real(kind=4), parameter    :: smf = 1.e-30
   real(kind=4)               :: t, dtmul, dtbig, dtsm, dtloc, dtsmmx, delt, dtn, &
                                 dtloc2, dtout, dtn1, dtn2

!  following arrays used in original scalr version has one to one
!  corresondance with above defined arrays,
!  note that BAQ is B (see below)
   real(kind=4)  :: dt, b(5, 2)
   real(kind=4)  :: raq(25, 2)
   real(kind=4)  :: tempk, tempi, psacw
   real(kind=4)  :: gaz_conc(maxnsg), aq(maxnsaq, 2)
   real(kind=4)  :: dg(maxnsg), pg(maxnsg), cdg(maxnsg)
   real(kind=4)  :: daq(maxnsaq, 2), paq(maxnsaq, 2), cdaq(maxnsaq, 2)
   real(kind=4)  :: q(2)
   real(kind=4)  :: edt(maxnsaq), xc(maxnsaq)

!  scalar version of dochem operates at each grid cell at a time for all 25 reactions
!  vector version operates at (ni*nz i.e 71*28) points at a time
!  Threfore in order to keep consistency with dta structures, arrays in argument
!  list are kept like, they were in vector version
!  dimension value "5" in following arrays stand for number of stoichmetric coefficents
!                  "25"                     .......  number of reactions
!
!  THIS IS THE DRIVER FOR THE CHEMISTRY SOLVER.
!
!
!  scalar version was designed to process each grid separately.
!  Subset the data from vector version to the format needed by scalar version.
!  Process the array only from n1(start of L/w) in 1D array of (71*28) to the end.

   do igc = 1, nptsnz
      if (qvec(igc, 1) <= 0.0 .and. qvec(igc, 2) <= 0.0) cycle
!  vector gvec(nptsnz, maxnsg)  ; scalar gaz_conc(maxnsg)
      do imaxnsg = 1, maxnsg
         gaz_conc(imaxnsg) = gvec(igc, imaxnsg)
      end do
!  vector aqvec(nptsnz, maxnsaq, 2) ; scalar aq(maxnsaq, 2)
      do imaxns2 = 1, 2
         do imaxnsaq = 1, maxnsaq
            aq(imaxnsaq, imaxns2) = aqvec(igc, imaxnsaq, imaxns2)
         end do
      end do
!  vector qvec(nptsnz, 2) ; scalar q(2)
      do iqcntr = 1, 2
         q(iqcntr) = qvec(igc, iqcntr)
      end do
!  vector bvec(nptsnz, 5, 2) ; scalar b(5, 2)
      do ibcntr1 = 1, 2
         do ibcntr2 = 1, 5
            b(ibcntr2, ibcntr1) = bvec(igc, ibcntr2, ibcntr1)
         end do
      end do
!  vector tempkvec(nptsnz) ; scalar tempk
      tempk = tempkvec(igc)
      tempi = tempivec(igc)
!  vector raqvec(nptsnz, 25, 2) ; scalar raq(25, 2)
      do iraqc1 = 1, 2
         do iraqc2 = 1, 25
            raq(iraqc2, iraqc1) = raqvec(igc, iraqc2, iraqc1)
         end do
      end do
!  vector psacwvec(nptsnz) ; scalar psacw
      psacw = psacwvec(igc)

!  end of data rearranging

      t = tin
      nstep = 0
!  initial ratio of time steps (large time step/small time step)
      ndt = 5
      dtmul = 1.

      if (q(1) <= 0.) then
         ndt = 1
         dtmul = 5.
      end if

      dtbig = dtsrt(idrf) * dtmul * real(ndt)
      dtsm = dtbig / real(ndt)
      dtbig = min(dtbig, dtmax)
      dtbig = min(dtbig, tout - tin)

!  start integration

!  if cloudwater = 0,  integrate aqueous species analytically
      if (q(1) > 0.0) then
!
         q1_loop: do
            nstep = nstep + 1
!  here nstmax=500
            if (nstep > nstmax) then
               write(0, *) '### Error in mach_incld_dochem ###'
               write(0, *) '# Maximum number of iterative steps exceeded'
               write(0, *) '###         ABORT         ###'
               chm_error_l = .true.
               return
            end if

!  integrate short term species
            dtloc = 0.
            dtsmmx = dtsm
            short_sp: do
               dtsm = min(dtsm, dtbig)

               call mach_incld_intsca_il(plus1, ng1, plus1, naq1, minus1, &
                                         kiter, dtsm, gaz_conc, aq, q, b, &
                                         tempk, tempi, raq, psacw)

!   exit:  integration fails
               if (chm_error_l) then
                  write(0, *) '### Error in mach_incld_dochem ###'
                  write(0, *) '# integration for fast species failed...'
                  write(0, *) '###         ABORT         ###'

                  write(0, *) 'SYL gaz_conc'
                  do ii = 1, maxnsg
                     write(0, *) 'gaz_conc(ii), ii ',gaz_conc(ii), ii
                  end do

                  write(0, *) 'SYL AQ'
                  do ii = 1, maxnsaq
                     write(0, *) 'aq(ii,1), aq(ii,2), ii ', aq(ii,1), aq(ii,2), ii
                  end do

                  write(0, *) 'SYL Q'
                  write(0, *) 'q(1), q(2) ',q(1), q(2)

                  return
               end if

!  increment time .....exit if dtbig reached or if nstep = 1
               dtloc = dtloc + dtsm
               if ((nstep == 1) .or. (dtloc + 1.e-2 * dtbig) >= dtbig) exit short_sp
               delt = dtbig - dtloc

!  select time step for next step
               dtn = mach_incld_fn_dtnew(dtsm, kiter)
               dtsm = min(dtn, dtmax * dtmul, delt)
               dtsmmx = max(dtsmmx, dtsm)
            end do short_sp

!  integrate long-lived species
            dt = dtloc
            dtloc2 = 0.

            long_sp: do
               call mach_incld_intsca_il(ng1p1, ng2, naq1, naq2, plus1, &
                                         kiter, dt, gaz_conc, aq, q, b, &
                                         tempk, tempi, raq, psacw)

               dtloc2 = dtloc2 + dt
!  exit: integration fails

               if (chm_error_l) then
                  write(0, *) '### Error in mach_incld_dochem ###'
                  write(0, *) '# integration for slow species failed...'
                  write(0, *) '###         ABORT         ###'
                  return
               end if

!  increment time .....exit if tout reached
               if ((dtloc2 + 1.e-3 * dtloc) >= dtloc) exit long_sp
               dt = dtloc - dtloc2
            end do long_sp

!  check for negative concentrations
            gaz_conc = max(gaz_conc, gmin)
            
            t = t + dtloc
            if ((t + 1.e-2) >= tout) exit q1_loop
            dtout = tout - t

!  select time step for next step
            dtn1 = mach_incld_fn_dtnew(dtloc, kiter)
            dtn2 = max(dtbig, dtn1)
            dtbig = min(dtn2, dtmax * dtmul, dtout)
            dtsm = max(dtsmmx, dtn2 / real(ndt))
            dtsmmx = dtsm
!  return to top of integration loop for next time step
         end do q1_loop
      
!  integrate once for only snow present
      else if (q(2) > 0.0) then
         q2_loop: do
            nstep = nstep + 1

            call mach_incld_diffun(gaz_conc, pg, dg, cdg, aq, paq, daq, cdaq, &
                                   q, plus1, b, tempk, tempi, raq, psacw)

            do i = 1, ng2
               edt(i) = exp(-dg(i) * dtbig)
               if (dg(i) < 1.e-15) then
                  xc(i) = gaz_conc(i) + cdg(i) * dtbig
               else
                  xc(i) = (pg(i) / (dg(i) + smf)) * &
                          (1.0 - edt(i)) + gaz_conc(i) * edt(i)
               endif
               if (pg(i) > 1.e-20 .or. gaz_conc(i) > 1.e-18 ) then
                  gaz_conc(i) = xc(i)
               endif
            end do
            do i = 1, naq2
               edt(i) = exp(-daq(i, 2) * dtbig)
               if (daq(i, 2)<1.e-15 ) then
                  xc(i) = aq(i, 2) + cdaq(i, 2) * dtbig
               else
                  xc(i) = (paq(i, 2) / (daq(i, 2) + smf)) * &
                          (1.0 - edt(i)) + aq(i, 2) * edt(i)
               endif
               if (paq(i, 2) > 1.e-20 .or. aq(i, 2) > 1.e-18 ) then
                  aq(i, 2) = xc(i)
               endif
            end do

!  steady state solution for hco3-
            aq(8, 2) = paq(8, 2) / daq(8, 2)
            call mach_incld_steady(aq, 2)
         
            do i = 1, maxnsaq
               aq(i, 2) = max(aq(i, 2), aqmin)
            end do

!  check for negative concentrations
            gaz_conc = max(gaz_conc, gmin)

            t = t + dtbig
!  exit out of main q(2) loop
            if ((t + 1.e-2) >= tout) exit q2_loop
            dtout = tout - t

!  select time step for next step
            dtbig = 10. * dtbig
            dtn = dtbig
            dtbig = min(dtn, dtmax * dtmul * real(ndt), dtout)     

! --- return to top of integration loop for next step
         end do q2_loop
      end if

!  pass back the output arrays.
!  now rearrange the output data from scalr format to vector format
      do imaxnsg = 1, maxnsg
         gvec(igc, imaxnsg) = gaz_conc(imaxnsg)
      end do
      do imaxns2 = 1, 2
         do imaxnsaq = 1, maxnsaq
            aqvec(igc, imaxnsaq, imaxns2) = aq(imaxnsaq, imaxns2)
         end do
      end do

   end do
!  end of grid cell loop (main loop)

end
