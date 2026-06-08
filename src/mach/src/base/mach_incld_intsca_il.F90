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
! Fichier/File   : mach_incld_intsca_il.ftn90
! Creation       : S. Menard, S. Gravel, B. Pabla, GEM-MACH, June 2008
! Description    : Checks stiffness of system of equations and solves the
!                  predictor-corrector equations, with iterations and convergence check.
!
! Extra info     : QA ADOM VERSION: ADOMIIB(CLEAN) LEVEL: 06/09/89 INTRQF ENSR(AES)
!
! Arguments  IN
!               n1     -->
!               n2     -->
!               n3     -->
!               n4     -->
!               ideriv -->
!               q      --> Cloudwater, ice/snow
!               b      --> Aqueous variable coefficients
!               tempk  --> Atmospheric Temperature (Kelvin)
!               tempi  --> Inverse of tempk (1.0/tempk)
!               psacw  -->
!
!            OUT
!               KITER  -->
!
!            INOUT
!               DT     --> Timestep
!               RAQ    --> Rate constants for aqueous-phase
!               G      --> Gas/part species conc (ppm)
!               AQ     --> Aqueous species concentrations (molar) in cloudwater
!=============================================================================
!
!!if_on
subroutine mach_incld_intsca_il(n1, n2, n3, n4, ideriv, kiter, dt, &
                                gaz_conc, aq, q, b, tempk, tempi, raq, psacw)
   use mach_cam_utils_mod,     only: maxnsg, maxnsaq
!!if_off
   use chm_utils_mod,          only: chm_error_l
   use mach_incld_headers_mod, only: mach_incld_diffun, mach_incld_steady
   use mach_incld_mod,         only: aqmin, incld_dtmin, gmin
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: n1
   integer(kind=4), intent   (in) :: n2
   integer(kind=4), intent   (in) :: n3
   integer(kind=4), intent   (in) :: n4
   integer(kind=4), intent   (in) :: ideriv
   integer(kind=4), intent  (out) :: kiter
   real(kind=4),    intent(inout) :: dt
   real(kind=4),    intent(inout) :: gaz_conc(maxnsg)
   real(kind=4),    intent(inout) :: aq(maxnsaq, 2)
   real(kind=4),    intent   (in) :: q(2)
   real(kind=4),    intent   (in) :: b(5, 2)
   real(kind=4),    intent   (in) :: tempk
   real(kind=4),    intent   (in) :: tempi
   real(kind=4),    intent(inout) :: raq(25, 2)
   real(kind=4),    intent   (in) :: psacw
!!if_off
!
! Local variables
!
   integer(kind=4)                    :: i, j, kneg, ineg, kloop, jj, iconv
   integer(kind=4), parameter         :: niter = 3
   real(kind=4)                       :: dthalf
   real(kind=4),    parameter         :: epsg = .003, epsq = .003, &
                                         tepsg = .03, tepsq = .03
   real(kind=4)                       :: xcrit, xneg, xtemi, xtempi, xsxi, eri,&
                                         ersumg, errorg, erraq, erraqm, dum2
   real(kind=4)                       :: smf = 1.e-30, one = 1.0
   real(kind=4), dimension(maxnsg)    :: dg, pg, cdg, sg, sdg, spg, scdg, &
                                         yc1g, yc2g, erg
   real(kind=4), dimension(maxnsaq, 2):: daq, paq, cdaq, saq, sdaq, spaq, &
                                         scdaq, yc1q, yc2q
   real(kind=4), dimension(maxnsaq)   :: xtemp, eraq
   real(kind=4), dimension(2)         :: ersumq, errorq
   logical(kind=4)                    :: lstifg(maxnsg), lstifq(maxnsaq, 2)
!
!  ygmin  = error control minimum concentrations for gas/aerosols
!  yaqmin = error control minimum concentrations for aqueous

   real(kind=4), dimension(maxnsg),  parameter :: ygmin =  &
                        (/1.e-4, 1.e-5, 1.e-5, 1.e-4, 1.e-4, 1.e-4, &
                          1.e-5, 1.e-5, 1.e-5, 1.e-5, 0.0,   0.0  /)
   real(kind=4), dimension(maxnsaq), parameter :: yaqmin = &
                        (/1.e-12, 1.e-6, 1.e-6, 1.e-6, 1.e-6, 1.e-6, 1.e-6, &
                           0.0  , 0.0  , 0.0  , 0.0  , 0.0 /)
 
   this_step: do 
      do i = 1, maxnsg
         gaz_conc(i) = max(gaz_conc(i), gmin)
         sg(i) = gaz_conc(i)
         yc1g(i) = gaz_conc(i)
         yc2g(i) = gaz_conc(i)
      end do
      do j = 1, 2
         if (q(j) <= 0.) cycle
         do i = 1, maxnsaq
            aq(i, j) = max(aq(i, j), aqmin)
            if ((i==9 .or. i==10) .and. aq(i, j)< 1.e-12 ) then
               aq(i, j) = 1.e-12
            end if
            saq(i, j) = aq(i, j)
            yc1q(i, j) = aq(i, j)
            yc2q(i, j) = aq(i, j)
         end do
      end do

!  get derivatives

      call mach_incld_diffun(gaz_conc, pg, dg, cdg, aq, paq, daq, cdaq, q, &
                             ideriv, b, tempk, tempi, raq, psacw)

      xcrit = 0.9 / dt
      do i = n1, n2
         if (abs(dg(i)) > xcrit) then
            xtemp(i) = 1.00001
         else
            xtemp(i) = 0.0
         endif
         lstifg(i) = xtemp(i) > 0.
      end do
      do j = 1, 2
         if (q(j) <= 0.0) cycle
         do i = n3, n4
            if (abs(daq(i, j)) > xcrit) then
               xtemp(i) = 1.00001
            else
               xtemp(i) = 0.0
            endif
            lstifq(i, j) = xtemp(i) > 0.
         end do
      end do

!  predictor step
      predictor: do kneg = 1, 3
         if (kneg > 1) dt = dt * .5

!  gases / particles

!  start of inlining pred
         xneg = 0.
         dum2 = dt + smf
         do i = n1, n2
            xtemi = 2.0 / (dg(i) + smf)
            if (lstifg(i)) then
               xsxi = (gaz_conc(i) * (xtemi - dt) + xtemi * dt * pg(i)) / &
                      (xtemi + dum2)
            else
               xsxi = gaz_conc(i) + dt * cdg(i)
            end if

            if (pg(i) <= 1.e-20 .and. gaz_conc(i) <= 1.e-18) then
               sg(i) = gaz_conc(i)
            else
               sg(i) = xsxi
            end if

            if (sg(i) <= 1.e-30) then
               xneg = xneg - 1.00001
            else
               xneg = xneg + 1.00001
            end if
         end do
         ineg = xneg / real(n2 - n1 + 1)
         if (ineg <= 0 .and. kneg /= 3) cycle predictor
!  end of inlining of pred.f for gaseous species

!  aqueous species
         do j = 1, 2
            if (q(j) <= 0.) cycle

!   start of inlining of pred.f for aqueous species
            xneg = 0.
            dum2 = dt + smf
            do i = n3, n4
               xtemi = 2.0 / (daq(i, j) + smf)

               if (lstifq(i, j)) then
                  xsxi = (aq(i, j) * (xtemi - dt) + xtemi * dt * paq(i, j)) /  &
                         (xtemi + dum2)
               else
                  xsxi = aq(i, j) + dt * cdaq(i, j)
               end if

               if (paq(i, j) <= 1.e-20 .and. aq(i, j) <= 1.e-18) then
                  saq(i, j) = aq(i, j)
               else
                  saq(i, j) = xsxi
               end if

               if (saq(i, j) <= 1.e-30) then
                  xneg = xneg - 1.00001
               else
                  xneg = xneg + 1.00001
               end if
            end do
            ineg = xneg / real(n4 - n3 + 1)
            if (ineg <= 0 .and. kneg /= 3) cycle predictor
!  end of inlining of pred.f for aqueous species
         end do

!  exit loop if no negative concentrations in predictor step
         exit predictor
      end do predictor

!  corrector step
      do i = n1, n2
         yc1g(i) = max(sg(i), gmin)
      end do
      do j = 1, 2
         if (q(j) <= 0.) cycle
         do i = n3, n4
            yc1q(i, j) = max(saq(i, j), aqmin)
         end do
      end do

!  iteration loop
      do kloop = 1, niter
         ersumg = 0.0
         ersumq = 0.0

!  get derivatives
         call mach_incld_diffun(yc1g, spg, sdg, scdg, yc1q, spaq, sdaq, scdaq, &
                                q, ideriv, b, tempk, tempi, raq, psacw)

         errorg = -1.0
         errorq = -1.0

!  inlining of "corr" routine for gaseous species
         dthalf = dt * 0.5
         dum2 = dt + smf

         do i = n1, n2
            xtempi = one / (sdg(i) + smf) + one / (dg(i) + smf)
!  getting rid of "cvmgt" routine
            if (lstifg(i)) then
               xsxi = (dthalf * xtempi * (pg(i) + spg(i)) + &
                       gaz_conc(i) * (xtempi - dt)) / (xtempi + dum2)
            else
               xsxi = gaz_conc(i) + dthalf * (cdg(i) + scdg(i))
            end if

            if ((pg(i) + spg(i)) <= 2.e-20 .and. gaz_conc(i) <= 1.e-18) then
               yc2g(i) = gaz_conc(i)
            else
               yc2g(i) = xsxi
            end if
!  evaluate error
            eri = abs((yc2g(i) - yc1g(i)) / &
                  (max(min(yc1g(i), yc2g(i)), ygmin(i)) + smf))
            ersumg = ersumg + eri
            erg(i) = eri
         end do

         do i = n1, n2
            errorg = max(errorg, erg(i))
         end do
!  this finishes inlining of "corr" routine for gaseous species

!  compute average error
         ersumg = ersumg / real(n2 - n1 + 1)

!  aqueous species
         do j = 1, 2
            if (q(j) <= 0.0) cycle

!  inlining of "corr" routine for aqueous species
            dthalf = dt * 0.5
            dum2 = dt + smf

            do i = n3, n4
               xtempi = one / (sdaq(i, j) + smf) + one / (daq(i, j) + smf)
!  getting rid of "cvmgt" routine
               if (lstifq(i, j)) then
                  xsxi = (dthalf * xtempi * (paq(i, j) + spaq(i, j)) +  &
                         aq(i, j) * (xtempi - dt)) / (xtempi + dum2)
               else
                  xsxi = aq(i, j) + dthalf * (cdaq(i, j) + scdaq(i, j))
               end if

               if ((paq(i, j)+spaq(i, j))<=2.e-20 .and. aq(i, j)<=1.e-18) then
                  yc2q(i, j) = aq(i, j)
               else
                  yc2q(i, j) = xsxi
               end if
!  evaluate error
               eri = abs((yc2q(i, j) - yc1q(i, j)) / &
                     (max(min(yc1q(i, j), yc2q(i, j)), yaqmin(i)) + smf))
               ersumq(j) = ersumq(j) + eri
               eraq(i) = eri
            end do

            do i = n3, n4
               errorq(j) = max(errorq(j), eraq(i))
            end do
!  this finishes inlining of "corr" routine for aqueous species

!  compute average error
            ersumq(j) = ersumq(j) / real(n4 - n3 + 1)
         end do

         erraq = 0.0
         erraqm = 0.0
         jj = 2
         do j = 1, 2
            if (q(j) <= 0.) jj = jj - 1
         end do
      
         if (jj /= 0) then
!  lump average and maximum aqueous phase errors together
            do j = 1, 2
               if (q(j) <= 0.) cycle
               erraq = erraq + ersumq(j) / jj
               erraqm = erraqm + errorq(j) / jj
            end do
         end if

!   convergence test
!   error control method:
!   max error .lt. 10*eps and average error .lt. eps
         kiter = kloop
         iconv = 1
         if (errorg > tepsg) iconv = -1
         if (ersumg > epsg)  iconv = -2
         if (erraqm > tepsq) iconv = -3
         if (erraq > epsq)   iconv = -4

!   exit iteration loop if convergence obtained
         if (iconv > 0) exit this_step
!
!  use solution from this iteration for initial in next iteration
         do i = n1, n2
            yc1g(i) = max(gmin, yc2g(i))
         end do
         do j = 1, 2
            if (q(j) <= 0.) cycle
            do i = n3, n4
               yc1q(i, j) = max(yc2q(i, j), aqmin)
            end do
         end do
      end do

!  no convergence after niter iterations, select smaller time step
!  retry step beginning at predictor equations

      dt = max(dt * 0.5, incld_dtmin)
      if (dt <= incld_dtmin) then
         chm_error_l = .true.
         do i = n1, n2
            gaz_conc(i) = yc2g(i)
         end do
         do j = 1, 2
            if (q(j) <= 0.) cycle
            do i = n3, n4
               aq(i, j) = yc2q(i, j)
            end do
         end do
!         
         return
      end if
   end do this_step
!  convergence obtained

   if (ideriv >= 0) then
!  calculate diagnostic species
      do j = 1, 2
         if (q(j) <= 0.) cycle
         yc2q(8, j) = (paq(8, j) + spaq(8, j)) / (daq(8, j) + sdaq(8, j))
         call mach_incld_steady(yc2q, j)
      end do
   end if
!  convergence obtained

!  check for negative concentrations
   do i = n1, n2
      gaz_conc(i) = max(gmin, yc2g(i))
   end do
   do j = 1, 2
      if (q(j) <= 0.) cycle
      do i = 1, maxnsaq
         aq(i, j) = max(aqmin, yc2q(i, j))
      end do
   end do
    
   return
end subroutine mach_incld_intsca_il
