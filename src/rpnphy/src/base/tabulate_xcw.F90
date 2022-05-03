!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!-------------------------------------- LICENCE END --------------------------------------
!**S/P TABULATE_XCW - PART OF THE ISCCP CLOUD SIMULATOR PACKAGE
!
      SUBROUTINE TABULATE_XCW()

         print *,'TABULATE_XCW() called a stub - ABORT'
         call flush(6)
         stop
      end SUBROUTINE TABULATE_XCW

!!$!Author
!!$!        Jason Cole, MSC/Cloud Physics (Summer 2005)
!!$
!!$!Revisions
!!$! 001    ...
!!$
!!$!Object
!!$!        This subroutine tabulates
!!$!          XCW = ratio of cloud condensate mixing ratio (QC)
!!$!                to its mean value as a function of cumulative
!!$!          ***   cumulative probability (N1 points)
!!$!          ***   relative standrad deviation of CWC
!!$!               (N2 points, presently 0.100, 0.125, ....)
!!$!        Either a beta (DIST=1.0) or
!!$!               a gamma distribution (DIST=2.0) can be used
!!$
!!$
!!$      implicit none
!!$!!!#include <arch_specific.hf>
!!$#include "mcica.cdk"
!!$
!!$! Local variables
!!$      INTEGER I,J
!!$      REAL AVG, STD, A, B, BG, C, D, R, Q, BETA, ALPHA, PROB_MIN, &
!!$           PROB, DIST, A_GAMMA  !#, A_BETA
!!$
!!$! Loop over standard deviations
!!$
!!$      DO 10 J=1,N2
!!$
!!$! mean and std dev
!!$        AVG = 1.0
!!$        STD = 0.025*(J+3)
!!$!
!!$! upper and lower limits for beta dist (A and B).
!!$! upper limit for gamma dist (BG).
!!$! BG superficially optimized by P. Raisanen (July 2002?)
!!$!
!!$
!!$        A  = 0.0
!!$        B  = (5. + 5.*STD**2)*AVG
!!$        BG = (5. + 5.*STD**2)*AVG
!!$!
!!$! using mean and std dev, determine parameters of
!!$! beta dist and gamma dist
!!$!
!!$        C = (AVG - A) / (B - AVG)
!!$        D = ((B - A) / STD)**2
!!$        R = C * (D - 2.0 - C) / (C * (C**2 + 3.0 * C + 3.0) + 1.0)
!!$        Q = C * R
!!$
!!$        BETA  = AVG / STD**2
!!$        ALPHA = AVG * BETA
!!$
!!$        PROB_MIN = 0.0
!!$
!!$! - PROB = cumulative frequency
!!$! - A_BETA and A_GAMMA = returned value given PROB
!!$!
!!$        DO 20 I=1,N1
!!$          PROB = REAL(I-1.)/REAL(N1-1.)
!!$
!!$!         DIST = 1.0 ! <<<<< BETA distribution >>>>>
!!$!
!!$!         CALL ROOT_LIMIT (DIST, Q, R, A, B, ALPHA, BETA, PROB_MIN,
!!$!     +                    PROB, A_BETA)
!!$!         XCW(I,J) = A_BETA
!!$
!!$          DIST = 2.0 ! <<<<< GAMMA distribution >>>>>
!!$
!!$          CALL ROOT_LIMIT (DIST, Q, R, A, BG, ALPHA, BETA, PROB_MIN, &
!!$                          PROB, A_GAMMA)
!!$
!!$          XCW(I,J) = A_GAMMA
!!$ 20     CONTINUE
!!$ 10   CONTINUE
!!$
!!$      RETURN
!!$      END
!!$
!!$!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$      SUBROUTINE ROOT_LIMIT (DIST, Q, R, A, B, ALPHA, BETA, &
!!$                             A_PROB, C_PROB, VALUE)
!!$      implicit none
!!$      real :: DIST, Q, R, A, B, ALPHA, BETA, A_PROB, C_PROB, VALUE
!!$      real, parameter :: XACC = 0.00001, EPS = 1.0E-5
!!$      integer, parameter :: JMAX = 10000
!!$!      PARAMETER (XACC = 0.00001, JMAX = 10000, EPS = 1.0E-5)
!!$!     PARAMETER (XACC = 0.005, JMAX = 1000, EPS = 1.0E-5)
!!$
!!$      integer :: J
!!$      real :: DX, FF, FMID, RTBIS, X, X1, X2, XMID
!!$      !TODO: BETAI, GAMMP
!!$
!!$      IF (DIST .EQ. 1.0) THEN
!!$         X2   = (B - A) / (B - A)
!!$         FMID = BETAI(Q,R,X2) - A_PROB
!!$         X1   = (A+EPS - A) / (B - A)
!!$         FF   = BETAI(Q,R,X1) - A_PROB
!!$      ELSE IF (DIST .EQ. 2.0) THEN
!!$         X2   = B * BETA
!!$         FMID = GAMMP(ALPHA,X2) - A_PROB
!!$         X1   = (EPS) * BETA
!!$         FF   = GAMMP(ALPHA,X1) - A_PROB
!!$      END IF
!!$
!!$!      IF (FF * FMID .GE. 0.0) PAUSE
!!$!      IF (FF .LT. 0.0)THEN
!!$         RTBIS = X1
!!$         DX    = X2 - X1
!!$!      ELSE
!!$!         RTBIS = X2
!!$!         DX    = X1 - X2
!!$!      ENDIF
!!$      DO 11 J=1,JMAX
!!$         DX   = DX * .50
!!$         XMID = RTBIS + DX
!!$         X    = XMID
!!$         IF (DIST .EQ. 1.0) THEN
!!$            FMID = BETAI(Q,R,X) - C_PROB
!!$         ELSE IF (DIST .EQ. 2.0) THEN
!!$            FMID = GAMMP(ALPHA,X) - C_PROB
!!$         END IF
!!$         IF (FMID .LT. 0.0) RTBIS = XMID
!!$         IF (ABS(DX) .LT. XACC .OR. FMID .EQ. 0.0) go to 15
!!$11    CONTINUE
!!$      WRITE(*,*) 'too many bisections'
!!$15    CONTINUE
!!$
!!$      IF (DIST .EQ. 1.0) THEN
!!$         VALUE = RTBIS * (B - A) + A
!!$      ELSE IF (DIST .EQ. 2.0) THEN
!!$         VALUE = RTBIS / BETA
!!$      END IF
!!$
!!$      RETURN
!!$      END
!!$
!!$!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$      FUNCTION BETA_FTN(z,w)
!!$      REAL beta_ftn,w,z
!!$!U    USES gammln
!!$      REAL gammln
!!$      beta_ftn=exp(gammln(z)+gammln(w)-gammln(z+w))
!!$      return
!!$      END
!!$
!!$!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$      FUNCTION gammp(a,x)
!!$      REAL a,gammp,x
!!$!U    USES gcf,gser
!!$      REAL gammcf,gamser,gln
!!$      if(x.lt.0..or.a.le.0.) WRITE(*,*) 'bad arguments in gammp'
!!$      if(x.lt.a+1.)then
!!$        call gser(gamser,a,x,gln)
!!$        gammp=gamser
!!$      else
!!$        call gcf(gammcf,a,x,gln)
!!$        gammp=1.-gammcf
!!$      endif
!!$      return
!!$      END
!!$
!!$!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$      SUBROUTINE gcf(gammcf,a,x,gln)
!!$      INTEGER ITMAX
!!$      REAL a,gammcf,gln,x,EPS,FPMIN
!!$      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
!!$!U    USES gammln
!!$      INTEGER i
!!$      REAL an,b,c,d,del,h,gammln
!!$      gln=gammln(a)
!!$      b=x+1.-a
!!$      c=1./FPMIN
!!$      d=1./b
!!$      h=d
!!$      do 11 i=1,ITMAX
!!$        an=-i*(i-a)
!!$        b=b+2.
!!$        d=an*d+b
!!$        if(abs(d).lt.FPMIN)d=FPMIN
!!$        c=b+an/c
!!$        if(abs(c).lt.FPMIN)c=FPMIN
!!$        d=1./d
!!$        del=d*c
!!$        h=h*del
!!$        if(abs(del-1.).lt.EPS)goto 1
!!$11    continue
!!$      WRITE(*,*) 'a too large, ITMAX too small in gcf'
!!$
!!$1     gammcf=exp(-x+a*log(x)-gln)*h
!!$      return
!!$      END
!!$
!!$!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$      SUBROUTINE gser(gamser,a,x,gln)
!!$      INTEGER ITMAX
!!$      REAL a,gamser,gln,x,EPS
!!$      PARAMETER (ITMAX=100,EPS=3.e-7)
!!$!U    USES gammln
!!$      INTEGER n
!!$      REAL ap,del,sum,gammln
!!$      gln=gammln(a)
!!$      if(x.le.0.)then
!!$        if(x.lt.0.) WRITE(*,*) 'x < 0 in gser'
!!$        gamser=0.
!!$        return
!!$      endif
!!$      ap=a
!!$      sum=1./a
!!$      del=sum
!!$      do 11 n=1,ITMAX
!!$        ap=ap+1.
!!$        del=del*x/ap
!!$        sum=sum+del
!!$        if(abs(del).lt.abs(sum)*EPS)goto 1
!!$11    continue
!!$      WRITE(*,*) 'a too large, ITMAX too small in gser'
!!$1     gamser=sum*exp(-x+a*log(x)-gln)
!!$
!!$      return
!!$      END
!!$
!!$!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$      FUNCTION gammln(xx)
!!$      REAL gammln,xx
!!$      INTEGER j
!!$      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
!!$      SAVE cof,stp
!!$      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, &
!!$      24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
!!$      -.5395239384953d-5,2.5066282746310005d0/
!!$      x=xx
!!$      y=x
!!$      tmp=x+5.5d0
!!$      tmp=(x+0.5d0)*log(tmp)-tmp
!!$      ser=1.000000000190015d0
!!$      do 11 j=1,6
!!$        y=y+1.d0
!!$        ser=ser+cof(j)/y
!!$11    continue
!!$      gammln=tmp+log(stp*ser/x)
!!$
!!$      return
!!$      END
!!$
!!$!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$      FUNCTION betai(a,b,x)
!!$      REAL betai,a,b,x
!!$!U    USES betacf,gammln
!!$      REAL bt,betacf,gammln
!!$      if(x.lt.0..or.x.gt.1.)WRITE(*,*) 'bad argument x in betai'
!!$      if(x.eq.0..or.x.eq.1.)then
!!$        bt=0.
!!$      else
!!$        bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.-x))
!!$      endif
!!$      if(x.lt.(a+1.)/(a+b+2.))then
!!$        betai=bt*betacf(a,b,x)/a
!!$        return
!!$      else
!!$        betai=1.-bt*betacf(b,a,1.-x)/b
!!$        return
!!$      endif
!!$      END
!!$
!!$!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$      FUNCTION betacf(a,b,x)
!!$      INTEGER MAXIT
!!$      REAL betacf,a,b,x,EPS,FPMIN
!!$      PARAMETER (MAXIT=100,EPS=3.e-7,FPMIN=1.e-30)
!!$      INTEGER m,m2
!!$      REAL aa,c,d,del,h,qab,qam,qap
!!$      qab=a+b
!!$      qap=a+1.
!!$      qam=a-1.
!!$      c=1.
!!$      d=1.-qab*x/qap
!!$      if(abs(d).lt.FPMIN)d=FPMIN
!!$      d=1./d
!!$      h=d
!!$      do 11 m=1,MAXIT
!!$        m2=2*m
!!$        aa=m*(b-m)*x/((qam+m2)*(a+m2))
!!$        d=1.+aa*d
!!$        if(abs(d).lt.FPMIN)d=FPMIN
!!$        c=1.+aa/c
!!$        if(abs(c).lt.FPMIN)c=FPMIN
!!$        d=1./d
!!$        h=h*d*c
!!$
!!$        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
!!$        d=1.+aa*d
!!$        if(abs(d).lt.FPMIN)d=FPMIN
!!$        c=1.+aa/c
!!$        if(abs(c).lt.FPMIN)c=FPMIN
!!$        d=1./d
!!$        del=d*c
!!$        h=h*del
!!$        if(abs(del-1.).lt.EPS)goto 1
!!$11    continue
!!$      WRITE(*,*) 'a or b too big, or MAXIT too small in betacf'
!!$1     betacf=h
!!$      return
!!$      END
!!$!******************************************************************

