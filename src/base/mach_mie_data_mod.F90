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
! Fichier/File   : mach_mie_data_mod.ftn90
! Creation       : P.A. Makar, and D. Akingunola 
! Description    : Module to run and load the mie scattering optical properties

module mach_mie_data_mod
   use chm_consphychm_mod, only: pi
   implicit none

   integer(kind=4), parameter :: mxxb = 1001
   integer(kind=4), parameter :: mxx = 101
   integer(kind=4), parameter :: mxnang = 1000
   integer(kind=4) :: nwl_aod
   real(kind=4),    dimension(mxxb,mxx) :: qext_t, ssca_t, asym_t
   real(kind=4) :: l10alf(mxxb), l10df(mxx)

   save

 contains

!============================================================================
! Name           : mach_mie_optical_properties
!
! Description    : This program creates arrays of coefficients needed for a 
!                  lookup table foroptical depth, single scattering albedo,
!                  and asymmetry factor to be used in the Mie code.
!
! Arguments:  IN/OUT
!                
!============================================================================
   subroutine mach_mie_optical_properties()
!     
      implicit none
      integer(kind=4) :: nang0
      real(kind=8)    :: l10alf_r8(mxxb)
      real(kind=4)    :: alf, dfrac, wfrac
      real(kind=4)    :: qext(mxxb, mxx), qsca(mxxb, mxx), qback(mxxb, mxx)
      complex :: refrel
! wavelenght of visible light [ m ]
      real(kind=4), parameter :: LAM = 0.55E-6
      integer(kind=4) :: i, j
      real(kind=4)    :: nr, ni
!
! CONSTL to get units of m**-1
      real(kind=4), parameter :: CONSTL = 2.* PI / LAM
!
      nang0 = 1
!
! real part of refractive index
      do j = 1, mxxb
         l10alf_r8(j) = -1.52287874528034D0 + dble(j-1)*0.005D0
         l10alf(j) = real(l10alf_r8(j))
      end do
      do i = 1,mxx
         l10df(i) = -5.0+real(i-1)*0.05
         dfrac = 10.0**(l10df(i))
         wfrac = 1.0 - dfrac
         nr = 1.5 - 0.17 * wfrac
         ni = 0.01 * (1.0 - wfrac)
         refrel = cmplx(nr, ni)
! *** set up Mie parameters 
         do j = 1,mxxb
            alf = real (10.D0**(l10alf_r8(j)))
! *** Call extinction routine
            call bhmie(alf, refrel, nang0, qext(j, i), qsca(j, i), qback(j, i))
         end do 
      end do
!
!  Convert backscatter to asymmetry factor, and
!  Convert total scattering to single scattering albedo:
      do i = 1, mxx
         do j = 1, mxxb
            qback(j,i) = (qsca(j, i) - qback(j, i)) / qsca(j, i)
            if (qback(j,i) < 0.0) then
               write(6,*) 'asym < 0 detected'
            end if
            qsca(j, i) =qsca(j, i) / qext(j, i)
            if (qsca(j, i) < 0.0) then
               write(6, *) 'ssa < 0 detected'
            end if
         end do
      end do
!  Output
      qext_t = alog10(qext)
      ssca_t = alog10(qsca)
      asym_t = alog10(qback)

      return
   end subroutine mach_mie_optical_properties
!
!***********************************************************************************
!
   subroutine mach_mie_lookup(wfrac, rw, numconc, bext, ssca_out, asym_out, pni, pnk)
!============================================================================
! Name           : mach_mie_lookup
!
! Description    : !  Contains lookup table of coefficient for Mie scattering, used
!                     Estimate aerosol optical properties optical depth, asymmetry
!                     factor, and the single-scattering albedo
!
!===============================================================================
      use chm_nml_mod,        only: aero_opt_wavel
      use mach_cam_utils_mod, only: isize
      implicit none

!!if_on
      integer(kind=4), intent (in) :: pni, pnk
      real(kind=4),    intent (in) :: wfrac   (pni, pnk, isize)
      real(kind=4),    intent (in) :: rw      (pni, pnk, isize)
      real(kind=4),    intent (in) :: numconc (pni, pnk, isize)

      real(kind=4),    intent(out) :: bext    (nwl_aod, pni, pnk)
      real(kind=4),    intent(out) :: ssca_out(nwl_aod, pni, pnk)
      real(kind=4),    intent(out) :: asym_out(nwl_aod, pni, pnk)
!
!  Local variables:
!
!  wavelength of light (m)
!  conversion for summations from efficiencies to total cross-sections:
      real(kind=8)    :: l10alf_in, l10df_in
      real(kind=4)    :: dfracin, dfracd, alfd, alfd2, dfracd2
      real(kind=8)    :: alfdd, dfracdd
      integer(kind=4) :: andx, bndx
      real(kind=4)    :: qsca, qback
      real(kind=4)    :: qext_p, ssca_p, asym_p
      integer(kind=4) :: ii, kk, mm, nn
      real(kind=4)    :: alfin, fac
      real(kind=4), dimension(nwl_aod) :: constl
!
      constl = 2.0 * pi / aero_opt_wavel

      bext = 0.0
      ssca_out = 0.0
      asym_out = 0.0

      do mm = 1, isize
         do kk = 1, pnk
            do ii = 1, pni
!  Only bother with cases where there are more than 100 particles/m3
               if (numconc(ii, kk, mm) > 100.0) then
               ! Evaluate wavelength independent variables first 
                  dfracin = min(1.0, max((1.0 - wfrac(ii, kk, mm)), 1.E-05))
                  l10df_in = min(max(dlog10(dble(dfracin)), -5.D+00), 0.D+00)
                  bndx = 1 + int((l10df_in + 5.D0) / 0.05D0)
                  bndx = max(min(bndx, mxx - 1), 1)
                  dfracdd = max(min((l10df_in - l10df(bndx)) / &
                                (l10df(bndx + 1) - l10df(bndx)), 1.0D0), 0.0D0)
                  dfracd = real(dfracdd)
                  dfracd2 = 1.0 - dfracd
                  fac = pi * rw(ii, kk, mm) * rw(ii, kk, mm) * numconc(ii, kk, mm)
               
                  do nn = 1, nwl_aod
                     alfin = constl(nn) * rw(ii, kk, mm)
!
                     l10alf_in = min(max(dlog10(dble(alfin)), -1.52787875D+00), &
                                         3.47712125D+00)
!
                     andx = 1 + int((l10alf_in + 1.52787875D+00) / 0.005D0)
                     andx = max(min(andx, mxxb - 1), 1)
! Double precision bilinear interpolant:
                     alfdd = max(min((l10alf_in - l10alf(andx)) / &
                                 (l10alf(andx + 1) - l10alf(andx)), 1.0D0), 0.0D0)
! Convert to single precision:
                     alfd = real(alfdd)
!
                     alfd2 = 1.0 - alfd
!
                     qext_p = 10.0**((qext_t(andx, bndx) * alfd2 +               &
                                      qext_t(andx + 1, bndx) * alfd) * dfracd2 + &
                                     (qext_t(andx, bndx + 1) * alfd2 +           &
                                      qext_t(andx+1, bndx + 1) * alfd) * dfracd)
!
                     ssca_p = 10.0**((ssca_t(andx, bndx) * alfd2 +               &
                                      ssca_t(andx + 1, bndx) * alfd) * dfracd2 + &
                                     (ssca_t(andx, bndx + 1) * alfd2 +           &
                                      ssca_t(andx + 1, bndx + 1) * alfd) * dfracd)
!
                     asym_p = 10.0**((asym_t(andx, bndx) * alfd2 +               &
                                      asym_t(andx + 1, bndx) * alfd) * dfracd2 + &
                                     (asym_t(andx, bndx + 1) * alfd2 +           &
                                      asym_t(andx + 1, bndx + 1) * alfd) * dfracd)
!  The above 3 values are for a single particle:  need to convert back 
!  to scattering efficiences, 
!  Recover net scattering from single scattering albedo:
                     qsca = ssca_p * qext_p
                     qback = qsca * (1.0 - asym_p)
!  Multiply by fac and add to sum:
                     bext(nn, ii, kk) = bext(nn, ii, kk) + qext_p * fac
!  Note that at this stage ssca_t holds sigma_sca, asym_t hold sigma_back:
                     ssca_out(nn, ii, kk) = ssca_out(nn, ii, kk) + qsca * fac
                     asym_out(nn, ii, kk) = asym_out(nn, ii, kk) + qback * fac
                  end do ! nn
               end if
            end do ! ii
         end do ! kk
      end do ! mm
!
!  Recalculate net ssca_t and asym_t based on the above:
!  Convert backscatter to asymmetry factor:
      where (ssca_out == 0.0)
         asym_out = 1.0
      elsewhere
         asym_out = max(min(1.0, (ssca_out - asym_out) / ssca_out), 0.0)
      endwhere
      
!  Convert total scattering to single scattering albedo:
      where (bext == 0.0)
         ssca_out = 1.0
      elsewhere
         ssca_out = max(min(1.0, ssca_out / bext), 0.0)
      endwhere
!
      return
   end subroutine mach_mie_lookup
!
!***********************************************************************************
  
!
   SUBROUTINE BHMIE(X,REFREL,NANG,QEXT,QSCA,QBACK)
! Declare parameters:
      INTEGER, PARAMETER :: NMXX = 400000
!C Arguments:
      INTEGER NANG
      REAL QBACK,QEXT,QSCA,X
      COMPLEX REFREL
      COMPLEX S1(2*MXNANG-1) !,S2(2*MXNANG-1)
!C Local variables:
      INTEGER J,JJ,N,NSTOP,NMX,NN
      REAL APSI,APSI1,CHI,CHI0,CHI1,DANG,FN,P,PII, &
           RN,THETA,XSTOP,YMOD
      REAL AMU(MXNANG),PI(MXNANG),PI0(MXNANG),PI1(MXNANG),TAU(MXNANG)
      DOUBLE PRECISION PSI0,PSI1,PSI,DN,DX
      COMPLEX AN,AN1,BN,BN1,XI,XI1,Y
      COMPLEX D(NMXX)
!C***********************************************************************
!C Subroutine BHMIE is the Bohren-Huffman Mie scattering subroutine
!C    to calculate scattering and absorption by a homogenous isotropic
!C    sphere.
!C Given:
!C    X = 2*pi*a/lambda
!C    REFREL = (complex refr. index of sphere)/(real index of medium)
!C    NANG = number of angles between 0 and 90 degrees
!C           (will calculate 2*NANG-1 directions from 0 to 180 deg.)
!C           if called with NANG<2, will set NANG=2 and will compute
!C           scattering for theta=0,90,180.
!C Returns:
!C    S1(1 - 2*NANG-1) = -i*f_22 (incid. E perp. to scatt. plane,
!C                                scatt. E perp. to scatt. plane)
!C    S2(1 - 2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
!C                                scatt. E parr. to scatt. plane)
!C    QEXT = C_ext/pi*a**2 = efficiency factor for extinction
!C    QSCA = C_sca/pi*a**2 = efficiency factor for scattering
!C    QBACK = (dC_sca/domega)/pi*a**2
!C          = backscattering efficiency
!C
!C Original program taken from Bohren and Huffman (1983), Appendix A
!C Modified by B.T.Draine, Princeton Univ. Obs., 90/10/26
!C in order to compute <cos(theta)>
!C 91/05/07 (BTD): Modified to allow NANG=1
!C 91/08/15 (BTD): Corrected error (failure to initialize P)
!C 91/08/15 (BTD): Modified to enhance vectorizability.
!C 91/08/15 (BTD): Modified to make NANG=2 if called with NANG=1
!C 91/08/15 (BTD): Changed definition of QBACK.
!C 92/01/08 (BTD): Note that this version has been superceded by
!C                 fully double precision version = bhmie.f which,
!C                 unfortunately, is not standard f77.
!C                 However, retain this in case standard f77 version
!C                 is required for porting to some other system.
!C***********************************************************************
!C*** Safety checks
      IF (NANG > MXNANG) THEN
         write(6,*) '***Error: NANG > MXNANG in bhmie'
         STOP
      END IF
      IF (NANG < 2) NANG=2
!C*** Obtain pi:
      PII=4.E0*ATAN(1.E0)
      DX=X
      Y=X*REFREL
      YMOD=ABS(Y)
!C
!C*** Series expansion terminated after NSTOP terms
!C    Logarithmic derivatives calculated from NMX on down
      XSTOP=X+4.E0*X**0.3333+2.0
!C*** Original code:
!C      NMX=AMAX1(XSTOP,YMOD)+15
!C      NSTOP=XSTOP
!C*** Experimental code:
      NMX=1.0*AMAX1(XSTOP,YMOD)+15
      NSTOP=1.0*XSTOP
!C
!      write(6,*) 'Value of X and REFREL: ',x,refrel
      IF(NMX > NMXX) THEN
          WRITE(6,*) 'Error: NMX (',NMX,') > NMXX=',NMXX,' for |m|x=',YMOD,&
          'and X and REFREL = ',X,REFREL
          STOP
      ENDIF
!C*** Require NANG.GE.1 in order to calculate scattering intensities
      DANG=0.
      IF (NANG > 1) DANG=.5E0*PII/FLOAT(NANG-1)
      DO J=1,NANG
          THETA=FLOAT(J-1)*DANG
          AMU(J)=COS(THETA)
      END DO
      DO J=1,NANG
          PI0(J)=0.E0
          PI1(J)=1.E0
      END DO
      NN=2*NANG-1
      DO J=1,NN
          S1(J)=(0.E0,0.E0)
!          S2(J)=(0.E0,0.E0)
      END DO
!C
!C*** Logarithmic derivative D(J) calculated by downward recurrence
!C    beginning with initial value (0.,0.) at J=NMX
!C
      D(NMX)=(0.E0,0.E0)
      NN=NMX-1
      DO N=1,NN
          RN=NMX-N+1
          D(NMX-N)=(RN/Y)-(1.E0/(D(NMX-N+1)+RN/Y))
      END DO
!C
!C*** Riccati-Bessel functions with real argument X
!C    calculated by upward recurrence
!C
      PSI0=DCOS(DX)
      PSI1=DSIN(DX)
      CHI0=-SIN(X)
      CHI1=COS(X)
!C APSI0 never used, so this line removed from program:
!C      APSI0=PSI0
      APSI1=PSI1
!C XI0 never used, so this line removed from program:
!C      XI0=CMPLX(APSI0,-CHI0)
      XI1=CMPLX(APSI1,-CHI1)
      QSCA=0.E0
      P=-1.
      DO N=1,NSTOP
          DN=N
          RN=N
          FN=(2.E0*RN+1.E0)/(RN*(RN+1.E0))
          PSI=(2.E0*DN-1.E0)*PSI1/DX-PSI0
          APSI=PSI
          CHI=(2.E0*RN-1.E0)*CHI1/X-CHI0
          XI=CMPLX(APSI,-CHI)
!C
!C*** Store previous values of AN and BN for use
!C    in computation of g=<cos(theta)>
          IF(N.GT.1)THEN
              AN1=AN
              BN1=BN
          ENDIF
!C
!C*** Compute AN and BN:
          AN=(D(N)/REFREL+RN/X)*APSI-APSI1
          AN=AN/((D(N)/REFREL+RN/X)*XI-XI1)
          BN=(REFREL*D(N)+RN/X)*APSI-APSI1
          BN=BN/((REFREL*D(N)+RN/X)*XI-XI1)
!C
!C*** Augment sums for Qsca and g=<cos(theta)>
          QSCA=QSCA+(2.*RN+1.)*(CABS(AN)**2+CABS(BN)**2)
!C
!C*** Now calculate scattering intensity pattern
!C    First do angles from 0 to 90
          DO J=1,NANG
              JJ=2*NANG-J
              PI(J)=PI1(J)
              TAU(J)=RN*AMU(J)*PI(J)-(RN+1.E0)*PI0(J)
              S1(J)=S1(J)+FN*(AN*PI(J)+BN*TAU(J))
!              S2(J)=S2(J)+FN*(AN*TAU(J)+BN*PI(J))
          END DO
!C
!C*** Now do angles greater than 90 using PI and TAU from
!C    angles less than 90.
!C    P=1 for N=1,3,...; P=-1 for N=2,4,...
          P=-P
          DO J=1,NANG-1
              JJ=2*NANG-J
              S1(JJ)=S1(JJ)+FN*P*(AN*PI(J)-BN*TAU(J))
!              S2(JJ)=S2(JJ)+FN*P*(BN*PI(J)-AN*TAU(J))
          END DO
          PSI0=PSI1
          PSI1=PSI
          APSI1=PSI1
          CHI0=CHI1
          CHI1=CHI
          XI1=CMPLX(APSI1,-CHI1)
!C
!C*** Compute pi_n for next value of n
!C    For each angle J, compute pi_n+1
!C    from PI = pi_n , PI0 = pi_n-1
          DO J=1,NANG
              PI1(J)=((2.*RN+1.)*AMU(J)*PI(J)-(RN+1.)*PI0(J))/RN
              PI0(J)=PI(J)
          END DO
      END DO
!C
!C*** Have summed sufficient terms.
!C    Now compute QSCA,QEXT,QBACK
      QSCA=(2.E0/(X*X))*QSCA
      QEXT=(4.E0/(X*X))*REAL(S1(1))
!C  
!C*** Comment added by P.A. Makar, February 2013:  note that
!C    the following definition of Qback is what Bohren and Huffmann
!C    refer to as the "conventional" definition of the backscattering
!C    ratio (page 121 of their text):  the differential scattering cross section for 
!C    scattering into a unit solid angle around the backscattering 
!C    direction, which can lead to instances in which the 
!C    backscattering cross section for a sphere small compared with
!C    the wavelength is greater than the total scattering 
!C    cross-section.  Tests of the code suggest that this will NOT
!C    be the case down to particles <1.E-08 of the wavelength.  
!C    Also, Bohren and Huffman's
!C    original code would give values of QBACK that approach
!C    QSCA in magnitude, in turn giving asymmetry factors close
!C    to zero, whereas most texts such as Jacobson's assume
!C    forward scattering dominates, as does GEM's original internal
!C    code.  So the revision will be retained here.
!C    
      QBACK=CABS(S1(2*NANG-1))*CABS(S1(2*NANG-1))/(PII*X*X)
      RETURN
      END SUBROUTINE BHMIE

end module mach_mie_data_mod
