C Copyright 1981-2007 ECMWF
C 
C Licensed under the GNU Lesser General Public License which
C incorporates the terms and conditions of version 3 of the GNU
C General Public License.
C See LICENSE and gpl-3.0.txt for details.
C
C  this file contains a few modifications to the ECMWF original code 
C  introduced implicit none
C  appended _M8 to the original names to avoid accidents should the new vode
C    be used together with the original code
C  introduced a few top level entry points( setfft8, ffft8, etc ...) that needed
C     no work array and/or no explicit factor/trig constants initialization
C  made code work with 8 byte reals
C
C
	subroutine ffft_m8(a, n, inc, jump, lot, isign )
	implicit none
	integer n,inc, jump, lot, isign
	real*8  a(*)

	call setfft_M8( n )
	call fft_M8( a, inc, jump, lot, isign )
	return
	end

	subroutine setfft8( n )
	implicit none
	integer  n
	external set99_m8
	call setfft_M8( n )
	return
	end

	subroutine setfft_M8( n )
	implicit none
	integer  n
	external set99_m8

	integer npts
	real *8, pointer, dimension(:) :: trigs
	integer, parameter :: maxfac=20
	integer, dimension(maxfac) :: ifac
	common /QQQ_FFFT8_QQQ/ trigs,ifac,npts

	data npts /-1/

	if(n .eq. npts) return
	if(n .gt. npts) then
	  if(npts .gt. 0) deallocate(trigs)
	  allocate(trigs(n+2))
	endif
	npts = n
	ifac=0
	call set99_m8(trigs,ifac,npts)

	return
	end

	subroutine ffft8( a, inc, jump, lot, isign )
	implicit none
	integer inc, jump, lot, isign
	real*8  a(*)
	call fft_m8(a, inc, jump, lot, isign )
	return
	end

	subroutine fft_M8( a, inc, jump, lot, isign )
	implicit none
	integer inc, jump, lot, isign
	real*8  a(*)

	integer npts
	real *8, pointer, dimension(:) :: trigs
	integer, parameter :: maxfac=20
	integer, dimension(maxfac) :: ifac
	common /QQQ_FFFT8_QQQ/ trigs,ifac,npts

	external fft991_m8
C
C  CHANGE maxlot ON VECTOR MACHINES
C
	integer, parameter :: maxlot=16
	real *8 work(npts+2,maxlot)
	integer i,slice

	do i=1,lot,maxlot
	  slice=min(maxlot,lot+1-i)
	  call fft991_m8( a(1+(i-1)*jump), work,
     %                   trigs, ifac, inc, jump, npts, slice, isign)
	enddo
	return
	end

      SUBROUTINE FFT991_M8(A,WORK,TRIGS,IFAX,INC,JUMP,N,ILOT,ISIGN)
	implicit none
	integer INC,JUMP,N,ILOT,ISIGN
      REAL*8 A(*),WORK(*), TRIGS(*)
      INTEGER IFAX(*)
	integer nx,nblox,nvex,i,ia
	integer nfax,istart,nb,j,ila,igo,k,ifac,ierr
	integer ibase,jbase,jj,ii,ix,iz
C
C     SUBROUTINE 'FFT991' - MULTIPLE FAST REAL PERIODIC TRANSFORM
C     SUPERSEDES PREVIOUS ROUTINE 'FFT991'
C
C     REAL TRANSFORM OF LENGTH N PERFORMED BY REMOVING REDUNDANT
C     OPERATIONS FROM COMPLEX TRANSFORM OF LENGTH N
C
C     A IS THE ARRAY CONTAINING INPUT & OUTPUT DATA
C     WORK IS AN AREA OF SIZE (N+1)*MIN(ILOT,JP_LOT)
C     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
C     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N
C     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
C         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
C     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
C     N IS THE LENGTH OF THE DATA VECTORS
C     ILOT IS THE NUMBER OF DATA VECTORS
C     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
C           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
C
C     ORDERING OF COEFFICIENTS:
C         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
C         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
C
C     ORDERING OF DATA:
C         X(0),X(1),X(2),...,X(N-1), 0 , 0 ; (N+2) LOCATIONS REQUIRED
C
C     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS
C     IN PARALLEL
C
C     N MUST BE COMPOSED OF FACTORS 2,3 & 5 BUT DOES NOT HAVE TO BE EVEN
C
C     DEFINITION OF TRANSFORMS:
C     -------------------------
C
C     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
C         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
C
C     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
C               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
C
      real*8  zero, haf, two
      parameter( zero  = 0.0D0 )
      parameter( haf  = 0.5D0 )
      parameter( two   = 2.0D0 )
      IF(IFAX(10).NE.N) CALL SET99_M8(TRIGS,IFAX,N)
      NFAX=IFAX(1)
      NX=N+1
      IF (MOD(N,2).EQ.1) NX=N
C==cfse  NBLOX=1+(ILOT-1)/64
C==cfse  NVEX=ILOT-(NBLOX-1)*64
      NBLOX=1+(ILOT-1)/512
      NVEX=ILOT-(NBLOX-1)*512
      IF (ISIGN.EQ.-1) GO TO 300
C
C     ISIGN=+1, SPECTRAL TO GRIDPOINT TRANSFORM
C     -----------------------------------------
  100 CONTINUE
      ISTART=1
      DO 220 NB=1,NBLOX
      IA=ISTART
      I=ISTART
C==*vocl loop,novrec
      DO 110 J=1,NVEX
      A(I+INC)=haf*A(I)
      I=I+JUMP
  110 CONTINUE
      IF (MOD(N,2).EQ.1) GO TO 130
      I=ISTART+N*INC
C==*vocl loop,novrec
      DO 120 J=1,NVEX
      A(I)=haf*A(I)
      I=I+JUMP
  120 CONTINUE
  130 CONTINUE
      IA=ISTART+INC
      ILA=1
      IGO=+1
C
      DO 160 K=1,NFAX
      IFAC=IFAX(K+1)
      IERR=-1
      IF (IGO.EQ.-1) GO TO 140
      CALL RPASSM_M8(A(IA),A(IA+ILA*INC),WORK(1),WORK(IFAC*ILA+1),TRIGS,
     *    INC,1,JUMP,NX,NVEX,N,IFAC,ILA,IERR)
      GO TO 150
  140 CONTINUE
      CALL RPASSM_M8(WORK(1),WORK(ILA+1),A(IA),A(IA+IFAC*ILA*INC),TRIGS,
     *    1,INC,NX,JUMP,NVEX,N,IFAC,ILA,IERR)
  150 CONTINUE
      IF (IERR.NE.0) GO TO 500
      ILA=IFAC*ILA
      IGO=-IGO
      IA=ISTART
  160 CONTINUE
C
C     IF NECESSARY, COPY RESULTS BACK TO A
C     ------------------------------------
      IF (MOD(NFAX,2).EQ.0) GO TO 190
      IBASE=1
      JBASE=IA
      DO 180 JJ=1,NVEX
      I=IBASE
      J=JBASE
      DO 170 II=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
  170 CONTINUE
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
  180 CONTINUE
  190 CONTINUE
C
C     FILL IN ZEROS AT END
C     --------------------
      IX=ISTART+N*INC
C==*vocl loop,novrec
      DO 210 J=1,NVEX
      A(IX)=0.0
      A(IX+INC)=0.0
      IX=IX+JUMP
  210 CONTINUE
C
      ISTART=ISTART+NVEX*JUMP
C==cfse  NVEX=64
      NVEX=512
  220 CONTINUE
      RETURN
C
C     ISIGN=-1, GRIDPOINT TO SPECTRAL TRANSFORM
C     -----------------------------------------
  300 CONTINUE
      ISTART=1
      DO 410 NB=1,NBLOX
      IA=ISTART
      ILA=N
      IGO=+1
C
      DO 340 K=1,NFAX
      IFAC=IFAX(NFAX+2-K)
      ILA=ILA/IFAC
      IERR=-1
      IF (IGO.EQ.-1) GO TO 320
      CALL QPASSM_M8(A(IA),A(IA+IFAC*ILA*INC),WORK(1),WORK(ILA+1),TRIGS,
     *    INC,1,JUMP,NX,NVEX,N,IFAC,ILA,IERR)
      GO TO 330
  320 CONTINUE
      CALL QPASSM_M8(WORK(1),WORK(IFAC*ILA+1),A(IA),A(IA+ILA*INC),TRIGS,
     *    1,INC,NX,JUMP,NVEX,N,IFAC,ILA,IERR)
  330 CONTINUE
      IF (IERR.NE.0) GO TO 500
      IGO=-IGO
      IA=ISTART+INC
  340 CONTINUE
C
C     IF NECESSARY, COPY RESULTS BACK TO A
C     ------------------------------------
      IF (MOD(NFAX,2).EQ.0) GO TO 370
      IBASE=1
      JBASE=IA
      DO 360 JJ=1,NVEX
      I=IBASE
      J=JBASE
      DO 350 II=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
  350 CONTINUE
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
  360 CONTINUE
  370 CONTINUE
C
C     SHIFT A(0) & FILL IN ZERO IMAG PARTS
C     ------------------------------------
      IX=ISTART
      DO 380 J=1,NVEX
      A(IX)=A(IX+INC)
      A(IX+INC)=0.0
      IX=IX+JUMP
  380 CONTINUE
      IF (MOD(N,2).EQ.1) GO TO 400
      IZ=ISTART+(N+1)*INC
      DO 390 J=1,NVEX
      A(IZ)=0.0
      IZ=IZ+JUMP
  390 CONTINUE
  400 CONTINUE
C
      ISTART=ISTART+NVEX*JUMP
C==cfse  NVEX=64
      NVEX=512
  410 CONTINUE
      RETURN
C
C     ERROR MESSAGES
C     --------------
  500 CONTINUE
      IF(IERR.EQ.1) THEN
  520 FORMAT('VECTOR LENGTH =',I4,', GREATER THAN 64')
      ELSEIF(IERR.EQ.2) THEN
      WRITE(6,540) IFAC
  540 FORMAT( 'FACTOR =',I3,', NOT CATERED FOR')
      ELSE
      WRITE(6,560) IFAC
  560 FORMAT('FACTOR =',I3,', ONLY CATERED FOR IF ILA*IFAC=N')
      ENDIF
      END
C Copyright 1981-2007 ECMWF
C 
C Licensed under the GNU Lesser General Public License which
C incorporates the terms and conditions of version 3 of the GNU
C General Public License.
C See LICENSE and gpl-3.0.txt for details.
C

      SUBROUTINE QPASSM_M8(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,ILOT,N,
     *    IFAC,ILA,IERR)
	implicit none
      REAL *8  A(*),B(*),C(*),D(*),TRIGS(*)
	integer INC1,INC2,INC3,INC4,ILOT,N,IFAC,ILA,IERR
	integer ijump,jl,m,iink,jink,kstop,ibad,ibase,jbase,igo
	integer i,j,k,ijk,ia,ib,ja,jb,kb,ic,jc,kc,id,jd,kd
	integer ie,je,ke,if,jf,kf,ig,ih
	real *8  c1,s1,c2,s2,c3,s3,c4,s4,c5,s5,z
	real *8  zsin60,sin45,zqrt5,zsin36,zsin72,zsin45
	real *8  a1,b1,a2,b2,a3,b3,a0,b0,a4,b4,a5,b5,a6,b6,a10,b10
	real *8  a11,b11,a20,b20,a21,b21
C
C     SUBROUTINE 'QPASSM' - PERFORMS ONE PASS THROUGH DATA AS PART
C     OF MULTIPLE REAL FFT (FOURIER ANALYSIS) ROUTINE
C
C     A IS FIRST REAL INPUT VECTOR
C         EQUIVALENCE B(1) WITH A(IFAC*ILA*INC1+1)
C     C IS FIRST REAL OUTPUT VECTOR
C         EQUIVALENCE D(1) WITH C(ILA*INC2+1)
C     TRIGS IS A PRECALCULATED LIST OF SINES & COSINES
C     INC1 IS THE ADDRESSING INCREMENT FOR A
C     INC2 IS THE ADDRESSING INCREMENT FOR C
C     INC3 IS THE INCREMENT BETWEEN INPUT VECTORS A
C     INC4 IS THE INCREMENT BETWEEN OUTPUT VECTORS C
C     ILOT IS THE NUMBER OF VECTORS
C     N IS THE LENGTH OF THE VECTORS
C     IFAC IS THE CURRENT FACTOR OF N
C     ILA = N/(PRODUCT OF FACTORS USED SO FAR)
C     IERR IS AN ERROR INDICATOR:
C              0 - PASS COMPLETED WITHOUT ERROR
C              1 - ILOT GREATER THAN 64
C              2 - IFAC NOT CATERED FOR
C              3 - IFAC ONLY CATERED FOR IF ILA=N/IFAC
C
C-----------------------------------------------------------------------
C
	real *8  SIN36, SIN72, QRT5, SIN60
      parameter( sin36 = 0.587785252292473D0)
      parameter( sin72 = 0.951056516295154D0)
      parameter( qrt5  = 0.559016994374947D0)
      parameter( sin60 = 0.866025403784437D0)
      real*8  zero, haf, two
      parameter( zero  = 0.0D0 )
      parameter( haf  = 0.5D0 )
      parameter( two   = 2.0D0 )
C
      M=N/IFAC
      IINK=ILA*INC1
      JINK=ILA*INC2
      IJUMP=(IFAC-1)*IINK
      KSTOP=(N-IFAC)/(2*IFAC)
C
      IBAD=1
C     no longer necessary, temporary arrays have dimension ILOT now
C     IF (ILOT.GT.512) GO TO 910
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.EQ.7) IGO=6
      IBAD=2
      IF (IGO.LT.1.OR.IGO.GT.6) GO TO 910
      GO TO (200,300,400,500,600,800),IGO
C
C     CODING FOR FACTOR 2
C     -------------------
  200 CONTINUE
      IA=1
      IB=IA+IINK
      JA=1
      JB=JA+(2*M-ILA)*INC2
C
      IF (ILA.EQ.M) GO TO 290
C
      DO 220 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 210 IJK=1,ILOT
      C(JA+J)=A(IA+I)+A(IB+I)
      C(JB+J)=A(IA+I)-A(IB+I)
      I=I+INC3
      J=J+INC4
  210 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  220 CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JA.EQ.JB) GO TO 260
      DO 250 K=ILA,KSTOP,ILA
      KB=K+K
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      JBASE=0
      DO 240 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 230 IJK=1,ILOT
      C(JA+J)=A(IA+I)+(C1*A(IB+I)+S1*B(IB+I))
      C(JB+J)=A(IA+I)-(C1*A(IB+I)+S1*B(IB+I))
      D(JA+J)=(C1*B(IB+I)-S1*A(IB+I))+B(IA+I)
      D(JB+J)=(C1*B(IB+I)-S1*A(IB+I))-B(IA+I)
      I=I+INC3
      J=J+INC4
  230 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  240 CONTINUE
      IBASE=IBASE+IJUMP
      JA=JA+JINK
      JB=JB-JINK
  250 CONTINUE
      IF (JA.GT.JB) GO TO 900
  260 CONTINUE
      JBASE=0
      DO 280 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 270 IJK=1,ILOT
      C(JA+J)=A(IA+I)
      D(JA+J)=-A(IB+I)
      I=I+INC3
      J=J+INC4
  270 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  280 CONTINUE
      GO TO 900
C
  290 CONTINUE
      Z=1.0D0/FLOAT(N)
      DO 294 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 292 IJK=1,ILOT
      C(JA+J)=Z*(A(IA+I)+A(IB+I))
      C(JB+J)=Z*(A(IA+I)-A(IB+I))
      I=I+INC3
      J=J+INC4
  292 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  294 CONTINUE
      GO TO 900
C
C     CODING FOR FACTOR 3
C     -------------------
  300 CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      JA=1
      JB=JA+(2*M-ILA)*INC2
      JC=JB
C
      IF (ILA.EQ.M) GO TO 390
C
      DO 320 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 310 IJK=1,ILOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      C(JB+J)=A(IA+I)-haf*(A(IB+I)+A(IC+I))
      D(JB+J)=SIN60*(A(IC+I)-A(IB+I))
      I=I+INC3
      J=J+INC4
  310 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  320 CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JA.EQ.JC) GO TO 360
      DO 350 K=ILA,KSTOP,ILA
      KB=K+K
      KC=KB+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      JBASE=0
      DO 340 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 330 IJK=1,ILOT
      A1=(C1*A(IB+I)+S1*B(IB+I))+(C2*A(IC+I)+S2*B(IC+I))
      B1=(C1*B(IB+I)-S1*A(IB+I))+(C2*B(IC+I)-S2*A(IC+I))
      A2=A(IA+I)-haf*A1
      B2=B(IA+I)-haf*B1
      A3=SIN60*((C1*A(IB+I)+S1*B(IB+I))-(C2*A(IC+I)+S2*B(IC+I)))
      B3=SIN60*((C1*B(IB+I)-S1*A(IB+I))-(C2*B(IC+I)-S2*A(IC+I)))
      C(JA+J)=A(IA+I)+A1
      D(JA+J)=B(IA+I)+B1
      C(JB+J)=A2+B3
      D(JB+J)=B2-A3
      C(JC+J)=A2-B3
      D(JC+J)=-(B2+A3)
      I=I+INC3
      J=J+INC4
  330 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  340 CONTINUE
      IBASE=IBASE+IJUMP
      JA=JA+JINK
      JB=JB+JINK
      JC=JC-JINK
  350 CONTINUE
      IF (JA.GT.JC) GO TO 900
  360 CONTINUE
      JBASE=0
      DO 380 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 370 IJK=1,ILOT
      C(JA+J)=A(IA+I)+haf*(A(IB+I)-A(IC+I))
      D(JA+J)=-SIN60*(A(IB+I)+A(IC+I))
      C(JB+J)=A(IA+I)-(A(IB+I)-A(IC+I))
      I=I+INC3
      J=J+INC4
  370 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  380 CONTINUE
      GO TO 900
C
  390 CONTINUE
      Z=1.0D0/FLOAT(N)
      ZSIN60=Z*SIN60
      DO 394 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 392 IJK=1,ILOT
      C(JA+J)=Z*(A(IA+I)+(A(IB+I)+A(IC+I)))
      C(JB+J)=Z*(A(IA+I)-haf*(A(IB+I)+A(IC+I)))
      D(JB+J)=ZSIN60*(A(IC+I)-A(IB+I))
      I=I+INC3
      J=J+INC4
  392 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  394 CONTINUE
      GO TO 900
C
C     CODING FOR FACTOR 4
C     -------------------
  400 CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      JA=1
      JB=JA+(2*M-ILA)*INC2
      JC=JB+2*M*INC2
      JD=JB
C
      IF (ILA.EQ.M) GO TO 490
C
      DO 420 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 410 IJK=1,ILOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
      C(JB+J)=A(IA+I)-A(IC+I)
      D(JB+J)=A(ID+I)-A(IB+I)
      I=I+INC3
      J=J+INC4
  410 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  420 CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC-JINK
      JD=JD-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JB.EQ.JC) GO TO 460
      DO 450 K=ILA,KSTOP,ILA
      KB=K+K
      KC=KB+KB
      KD=KC+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      JBASE=0
      DO 440 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 430 IJK=1,ILOT
      A0=A(IA+I)+(C2*A(IC+I)+S2*B(IC+I))
      A2=A(IA+I)-(C2*A(IC+I)+S2*B(IC+I))
      A1=(C1*A(IB+I)+S1*B(IB+I))+(C3*A(ID+I)+S3*B(ID+I))
      A3=(C1*A(IB+I)+S1*B(IB+I))-(C3*A(ID+I)+S3*B(ID+I))
      B0=B(IA+I)+(C2*B(IC+I)-S2*A(IC+I))
      B2=B(IA+I)-(C2*B(IC+I)-S2*A(IC+I))
      B1=(C1*B(IB+I)-S1*A(IB+I))+(C3*B(ID+I)-S3*A(ID+I))
      B3=(C1*B(IB+I)-S1*A(IB+I))-(C3*B(ID+I)-S3*A(ID+I))
      C(JA+J)=A0+A1
      C(JC+J)=A0-A1
      D(JA+J)=B0+B1
      D(JC+J)=B1-B0
      C(JB+J)=A2+B3
      C(JD+J)=A2-B3
      D(JB+J)=B2-A3
      D(JD+J)=-(B2+A3)
      I=I+INC3
      J=J+INC4
  430 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  440 CONTINUE
      IBASE=IBASE+IJUMP
      JA=JA+JINK
      JB=JB+JINK
      JC=JC-JINK
      JD=JD-JINK
  450 CONTINUE
      IF (JB.GT.JC) GO TO 900
  460 CONTINUE
      SIN45=DSQRT(haf)
      JBASE=0
      DO 480 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 470 IJK=1,ILOT
      C(JA+J)=A(IA+I)+SIN45*(A(IB+I)-A(ID+I))
      C(JB+J)=A(IA+I)-SIN45*(A(IB+I)-A(ID+I))
      D(JA+J)=-A(IC+I)-SIN45*(A(IB+I)+A(ID+I))
      D(JB+J)=A(IC+I)-SIN45*(A(IB+I)+A(ID+I))
      I=I+INC3
      J=J+INC4
  470 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  480 CONTINUE
      GO TO 900
C
  490 CONTINUE
      Z=1.0D0/FLOAT(N)
      DO 494 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 492 IJK=1,ILOT
      C(JA+J)=Z*((A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I)))
      C(JC+J)=Z*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
      C(JB+J)=Z*(A(IA+I)-A(IC+I))
      D(JB+J)=Z*(A(ID+I)-A(IB+I))
      I=I+INC3
      J=J+INC4
  492 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  494 CONTINUE
      GO TO 900
C
C     CODING FOR FACTOR 5
C     -------------------
  500 CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      JA=1
      JB=JA+(2*M-ILA)*INC2
      JC=JB+2*M*INC2
      JD=JC
      JE=JB
C
      IF (ILA.EQ.M) GO TO 590
C
      DO 520 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 510 IJK=1,ILOT
      A1=A(IB+I)+A(IE+I)
      A3=A(IB+I)-A(IE+I)
      A2=A(IC+I)+A(ID+I)
      A4=A(IC+I)-A(ID+I)
      A5=A(IA+I)-0.25*(A1+A2)
      A6=QRT5*(A1-A2)
      C(JA+J)=A(IA+I)+(A1+A2)
      C(JB+J)=A5+A6
      C(JC+J)=A5-A6
      D(JB+J)=-SIN72*A3-SIN36*A4
      D(JC+J)=-SIN36*A3+SIN72*A4
      I=I+INC3
      J=J+INC4
  510 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  520 CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JB.EQ.JD) GO TO 560
      DO 550 K=ILA,KSTOP,ILA
      KB=K+K
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      JBASE=0
      DO 540 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 530 IJK=1,ILOT
      A1=(C1*A(IB+I)+S1*B(IB+I))+(C4*A(IE+I)+S4*B(IE+I))
      A3=(C1*A(IB+I)+S1*B(IB+I))-(C4*A(IE+I)+S4*B(IE+I))
      A2=(C2*A(IC+I)+S2*B(IC+I))+(C3*A(ID+I)+S3*B(ID+I))
      A4=(C2*A(IC+I)+S2*B(IC+I))-(C3*A(ID+I)+S3*B(ID+I))
      B1=(C1*B(IB+I)-S1*A(IB+I))+(C4*B(IE+I)-S4*A(IE+I))
      B3=(C1*B(IB+I)-S1*A(IB+I))-(C4*B(IE+I)-S4*A(IE+I))
      B2=(C2*B(IC+I)-S2*A(IC+I))+(C3*B(ID+I)-S3*A(ID+I))
      B4=(C2*B(IC+I)-S2*A(IC+I))-(C3*B(ID+I)-S3*A(ID+I))
      A5=A(IA+I)-0.25*(A1+A2)
      A6=QRT5*(A1-A2)
      B5=B(IA+I)-0.25*(B1+B2)
      B6=QRT5*(B1-B2)
      A10=A5+A6
      A20=A5-A6
      B10=B5+B6
      B20=B5-B6
      A11=SIN72*B3+SIN36*B4
      A21=SIN36*B3-SIN72*B4
      B11=SIN72*A3+SIN36*A4
      B21=SIN36*A3-SIN72*A4
      C(JA+J)=A(IA+I)+(A1+A2)
      C(JB+J)=A10+A11
      C(JE+J)=A10-A11
      C(JC+J)=A20+A21
      C(JD+J)=A20-A21
      D(JA+J)=B(IA+I)+(B1+B2)
      D(JB+J)=B10-B11
      D(JE+J)=-(B10+B11)
      D(JC+J)=B20-B21
      D(JD+J)=-(B20+B21)
      I=I+INC3
      J=J+INC4
  530 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  540 CONTINUE
      IBASE=IBASE+IJUMP
      JA=JA+JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
  550 CONTINUE
      IF (JB.GT.JD) GO TO 900
  560 CONTINUE
      JBASE=0
      DO 580 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 570 IJK=1,ILOT
      A1=A(IB+I)+A(IE+I)
      A3=A(IB+I)-A(IE+I)
      A2=A(IC+I)+A(ID+I)
      A4=A(IC+I)-A(ID+I)
      A5=A(IA+I)+0.25*(A3-A4)
      A6=QRT5*(A3+A4)
      C(JA+J)=A5+A6
      C(JB+J)=A5-A6
      C(JC+J)=A(IA+I)-(A3-A4)
      D(JA+J)=-SIN36*A1-SIN72*A2
      D(JB+J)=-SIN72*A1+SIN36*A2
      I=I+INC3
      J=J+INC4
  570 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  580 CONTINUE
      GO TO 900
C
  590 CONTINUE
      Z=1.0D0/FLOAT(N)
      ZQRT5=Z*QRT5
      ZSIN36=Z*SIN36
      ZSIN72=Z*SIN72
      DO 594 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 592 IJK=1,ILOT
      A1=A(IB+I)+A(IE+I)
      A3=A(IB+I)-A(IE+I)
      A2=A(IC+I)+A(ID+I)
      A4=A(IC+I)-A(ID+I)
      A5=Z*(A(IA+I)-0.25*(A1+A2))
      A6=ZQRT5*(A1-A2)
      C(JA+J)=Z*(A(IA+I)+(A1+A2))
      C(JB+J)=A5+A6
      C(JC+J)=A5-A6
      D(JB+J)=-ZSIN72*A3-ZSIN36*A4
      D(JC+J)=-ZSIN36*A3+ZSIN72*A4
      I=I+INC3
      J=J+INC4
  592 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  594 CONTINUE
      GO TO 900
C
C     CODING FOR FACTOR 6
C     -------------------
  600 CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      IF=IE+IINK
      JA=1
      JB=JA+(2*M-ILA)*INC2
      JC=JB+2*M*INC2
      JD=JC+2*M*INC2
      JE=JC
      JF=JB
C
      IF (ILA.EQ.M) GO TO 690
C
      DO 620 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 610 IJK=1,ILOT
      A11=(A(IC+I)+A(IF+I))+(A(IB+I)+A(IE+I))
      C(JA+J)=(A(IA+I)+A(ID+I))+A11
      C(JC+J)=(A(IA+I)+A(ID+I)-haf*A11)
      D(JC+J)=SIN60*((A(IC+I)+A(IF+I))-(A(IB+I)+A(IE+I)))
      A11=(A(IC+I)-A(IF+I))+(A(IE+I)-A(IB+I))
      C(JB+J)=(A(IA+I)-A(ID+I))-haf*A11
      D(JB+J)=SIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
      C(JD+J)=(A(IA+I)-A(ID+I))+A11
      I=I+INC3
      J=J+INC4
  610 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  620 CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
      JF=JF-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JC.EQ.JD) GO TO 660
      DO 650 K=ILA,KSTOP,ILA
      KB=K+K
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      KF=KE+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      C5=TRIGS(KF+1)
      S5=TRIGS(KF+2)
      JBASE=0
      DO 640 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 630 IJK=1,ILOT
      A1=C1*A(IB+I)+S1*B(IB+I)
      B1=C1*B(IB+I)-S1*A(IB+I)
      A2=C2*A(IC+I)+S2*B(IC+I)
      B2=C2*B(IC+I)-S2*A(IC+I)
      A3=C3*A(ID+I)+S3*B(ID+I)
      B3=C3*B(ID+I)-S3*A(ID+I)
      A4=C4*A(IE+I)+S4*B(IE+I)
      B4=C4*B(IE+I)-S4*A(IE+I)
      A5=C5*A(IF+I)+S5*B(IF+I)
      B5=C5*B(IF+I)-S5*A(IF+I)
      A11=(A2+A5)+(A1+A4)
      A20=(A(IA+I)+A3)-haf*A11
      A21=SIN60*((A2+A5)-(A1+A4))
      B11=(B2+B5)+(B1+B4)
      B20=(B(IA+I)+B3)-haf*B11
      B21=SIN60*((B2+B5)-(B1+B4))
      C(JA+J)=(A(IA+I)+A3)+A11
      D(JA+J)=(B(IA+I)+B3)+B11
      C(JC+J)=A20-B21
      D(JC+J)=A21+B20
      C(JE+J)=A20+B21
      D(JE+J)=A21-B20
      A11=(A2-A5)+(A4-A1)
      A20=(A(IA+I)-A3)-haf*A11
      A21=SIN60*((A4-A1)-(A2-A5))
      B11=(B5-B2)-(B4-B1)
      B20=(B3-B(IA+I))-haf*B11
      B21=SIN60*((B5-B2)+(B4-B1))
      C(JB+J)=A20-B21
      D(JB+J)=A21-B20
      C(JD+J)=A11+(A(IA+I)-A3)
      D(JD+J)=B11+(B3-B(IA+I))
      C(JF+J)=A20+B21
      D(JF+J)=A21+B20
      I=I+INC3
      J=J+INC4
  630 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  640 CONTINUE
      IBASE=IBASE+IJUMP
      JA=JA+JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
      JF=JF-JINK
  650 CONTINUE
      IF (JC.GT.JD) GO TO 900
  660 CONTINUE
      JBASE=0
      DO 680 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 670 IJK=1,ILOT
      C(JA+J)=(A(IA+I)+haf*(A(IC+I)-A(IE+I)))+ SIN60*(A(IB+I)-A(IF+I))
      D(JA+J)=-(A(ID+I)+haf*(A(IB+I)+A(IF+I)))-SIN60*(A(IC+I)+A(IE+I))
      C(JB+J)=A(IA+I)-(A(IC+I)-A(IE+I))
      D(JB+J)=A(ID+I)-(A(IB+I)+A(IF+I))
      C(JC+J)=(A(IA+I)+haf*(A(IC+I)-A(IE+I)))-SIN60*(A(IB+I)-A(IF+I))
      D(JC+J)=-(A(ID+I)+haf*(A(IB+I)+A(IF+I)))+SIN60*(A(IC+I)+A(IE+I))
      I=I+INC3
      J=J+INC4
  670 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  680 CONTINUE
      GO TO 900
C
  690 CONTINUE
      Z=1.0D0/FLOAT(N)
      ZSIN60=Z*SIN60
      DO 694 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 692 IJK=1,ILOT
      A11=(A(IC+I)+A(IF+I))+(A(IB+I)+A(IE+I))
      C(JA+J)=Z*((A(IA+I)+A(ID+I))+A11)
      C(JC+J)=Z*((A(IA+I)+A(ID+I))-haf*A11)
      D(JC+J)=ZSIN60*((A(IC+I)+A(IF+I))-(A(IB+I)+A(IE+I)))
      A11=(A(IC+I)-A(IF+I))+(A(IE+I)-A(IB+I))
      C(JB+J)=Z*((A(IA+I)-A(ID+I))-haf*A11)
      D(JB+J)=ZSIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
      C(JD+J)=Z*((A(IA+I)-A(ID+I))+A11)
      I=I+INC3
      J=J+INC4
  692 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  694 CONTINUE
      GO TO 900
C
C     CODING FOR FACTOR 8
C     -------------------
  800 CONTINUE
      IBAD=3
      IF (ILA.NE.M) GO TO 910
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      IF=IE+IINK
      IG=IF+IINK
      IH=IG+IINK
      JA=1
      JB=JA+ILA*INC2
      JC=JB+2*M*INC2
      JD=JC+2*M*INC2
      JE=JD+2*M*INC2
      Z=1.0D0/FLOAT(N)
      ZSIN45=Z*SQRT(haf)
C
      DO 820 JL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 810 IJK=1,ILOT
      C(JA+J)=Z*(((A(IA+I)+A(IE+I))+(A(IC+I)+A(IG+I)))+
     *    ((A(ID+I)+A(IH+I))+(A(IB+I)+A(IF+I))))
      C(JE+J)=Z*(((A(IA+I)+A(IE+I))+(A(IC+I)+A(IG+I)))-
     *    ((A(ID+I)+A(IH+I))+(A(IB+I)+A(IF+I))))
      C(JC+J)=Z*((A(IA+I)+A(IE+I))-(A(IC+I)+A(IG+I)))
      D(JC+J)=Z*((A(ID+I)+A(IH+I))-(A(IB+I)+A(IF+I)))
      C(JB+J)=Z*(A(IA+I)-A(IE+I))
     *    +ZSIN45*((A(IH+I)-A(ID+I))-(A(IF+I)-A(IB+I)))
      C(JD+J)=Z*(A(IA+I)-A(IE+I))
     *    -ZSIN45*((A(IH+I)-A(ID+I))-(A(IF+I)-A(IB+I)))
      D(JB+J)=ZSIN45*((A(IH+I)-A(ID+I))+(A(IF+I)-A(IB+I)))
     *    +Z*(A(IG+I)-A(IC+I))
      D(JD+J)=ZSIN45*((A(IH+I)-A(ID+I))+(A(IF+I)-A(IB+I)))
     *    -Z*(A(IG+I)-A(IC+I))
      I=I+INC3
      J=J+INC4
  810 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  820 CONTINUE
C
C     RETURN
C     ------
  900 CONTINUE
      IBAD=0
  910 CONTINUE
      IERR=IBAD
      RETURN
      END
C Copyright 1981-2007 ECMWF
C 
C Licensed under the GNU Lesser General Public License which
C incorporates the terms and conditions of version 3 of the GNU
C General Public License.
C See LICENSE and gpl-3.0.txt for details.
C

      SUBROUTINE RPASSM_M8(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,ILOT,N,
     *    IFAC,ILA,IERR)
	implicit none
      REAL *8 A(*),B(*),C(*),D(*),TRIGS(*)
	integer INC1,INC2,INC3,INC4,ILOT,N,IFAC,ILA,IERR
	integer m,iink,jink,jump,kstop,ibad,ibase,jbase,igo
	integer ia,ib,ja,jb,il,i,j,ijk,k,kb,ic,jc,kc
	real *8 c1,c2,s1,s2,ssin60,c3,s3,c4,s4,c5,s5,sin45
	real *8 qqrt5,ssin45,ssin36,ssin72
	integer id,jd,kd,ie,je,ke,if,jf,kf,jg,jh
C
C     SUBROUTINE 'RPASSM' - PERFORMS ONE PASS THROUGH DATA AS PART
C     OF MULTIPLE REAL FFT (FOURIER SYNTHESIS) ROUTINE
C
C     A IS FIRST REAL INPUT VECTOR
C         EQUIVALENCE B(1) WITH A (ILA*INC1+1)
C     C IS FIRST REAL OUTPUT VECTOR
C         EQUIVALENCE D(1) WITH C(IFAC*ILA*INC2+1)
C     TRIGS IS A PRECALCULATED LIST OF SINES & COSINES
C     INC1 IS THE ADDRESSING INCREMENT FOR A
C     INC2 IS THE ADDRESSING INCREMENT FOR C
C     INC3 IS THE INCREMENT BETWEEN INPUT VECTORS A
C     INC4 IS THE INCREMENT BETWEEN OUTPUT VECTORS C
C     ILOT IS THE NUMBER OF VECTORS
C     N IS THE LENGTH OF THE VECTORS
C     IFAC IS THE CURRENT FACTOR OF N
C     ILA IS THE PRODUCT OF PREVIOUS FACTORS
C     IERR IS AN ERROR INDICATOR:
C              0 - PASS COMPLETED WITHOUT ERROR
C              1 - ILOT GREATER THAN 64
C              2 - IFAC NOT CATERED FOR
C              3 - IFAC ONLY CATERED FOR IF ILA=N/IFAC
C
C-----------------------------------------------------------------------
C
      REAL *8 A10(ILOT),A11(ILOT),
     *     A20(ILOT),A21(ILOT),
     *     B10(ILOT),B11(ILOT),B20(ILOT),B21(ILOT)
	real *8  SIN36, SIN72, QRT5, SIN60
      parameter( sin36 = 0.587785252292473D0)
      parameter( sin72 = 0.951056516295154D0)
      parameter( qrt5  = 0.559016994374947D0)
      parameter( sin60 = 0.866025403784437D0)
      real*8  zero, haf, two
      parameter( zero  = 0.0D0 )
      parameter( haf  = 0.5D0 )
      parameter( two   = 2.0D0 )
C
      M=N/IFAC
      IINK=ILA*INC1
      JINK=ILA*INC2
      JUMP=(IFAC-1)*JINK
      KSTOP=(N-IFAC)/(2*IFAC)
C
      IBAD=1
C     no longer necessary, temporary arrays have dimension ILOT now
C     IF (ILOT.GT.512) GO TO 910
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.EQ.7) IGO=6
      IBAD=2
      IF (IGO.LT.1.OR.IGO.GT.6) GO TO 910
      GO TO (200,300,400,500,600,800),IGO
C
C     CODING FOR FACTOR 2
C     -------------------
  200 CONTINUE
      IA=1
      IB=IA+(2*M-ILA)*INC1
      JA=1
      JB=JA+JINK
C
      IF (ILA.EQ.M) GO TO 290
C
      DO 220 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 210 IJK=1,ILOT
      C(JA+J)=A(IA+I)+A(IB+I)
      C(JB+J)=A(IA+I)-A(IB+I)
      I=I+INC3
      J=J+INC4
  210 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  220 CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB-IINK
      IBASE=0
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IA.EQ.IB) GO TO 260
      DO 250 K=ILA,KSTOP,ILA
      KB=K+K
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      IBASE=0
      DO 240 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 230 IJK=1,ILOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)-B(IB+I)
      C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)+B(IB+I))
      D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)+B(IB+I))
      I=I+INC3
      J=J+INC4
  230 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  240 CONTINUE
      IA=IA+IINK
      IB=IB-IINK
      JBASE=JBASE+JUMP
  250 CONTINUE
      IF (IA.GT.IB) GO TO 900
  260 CONTINUE
      IBASE=0
      DO 280 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 270 IJK=1,ILOT
      C(JA+J)=A(IA+I)
      C(JB+J)=-B(IA+I)
      I=I+INC3
      J=J+INC4
  270 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  280 CONTINUE
      GO TO 900
C
  290 CONTINUE
      DO 294 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 292 IJK=1,ILOT
      C(JA+J)=two*(A(IA+I)+A(IB+I))
      C(JB+J)=two*(A(IA+I)-A(IB+I))
      I=I+INC3
      J=J+INC4
  292 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  294 CONTINUE
      GO TO 900
C
C     CODING FOR FACTOR 3
C     -------------------
  300 CONTINUE
      IA=1
      IB=IA+(2*M-ILA)*INC1
      IC=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
C
      IF (ILA.EQ.M) GO TO 390
C
      DO 320 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 310 IJK=1,ILOT
      C(JA+J)=A(IA+I)+A(IB+I)
      C(JB+J)=(A(IA+I)-haf*A(IB+I))-(SIN60*(B(IB+I)))
      C(JC+J)=(A(IA+I)-haf*A(IB+I))+(SIN60*(B(IB+I)))
      I=I+INC3
      J=J+INC4
  310 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  320 CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IA.EQ.IC) GO TO 360
      DO 350 K=ILA,KSTOP,ILA
      KB=K+K
      KC=KB+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      IBASE=0
      DO 340 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 330 IJK=1,ILOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)-B(IC+I))
      C(JB+J)=
     *    C1*((A(IA+I)-haf*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)+B(IC+
     *  I))))
     *   -S1*((B(IA+I)-haf*(B(IB+I)-B(IC+I)))+(SIN60*(A(IB+I)-A(IC+
     *  I))))
      D(JB+J)=
     *    S1*((A(IA+I)-haf*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)+B(IC+
     *  I))))
     *   +C1*((B(IA+I)-haf*(B(IB+I)-B(IC+I)))+(SIN60*(A(IB+I)-A(IC+
     *  I))))
      C(JC+J)=
     *    C2*((A(IA+I)-haf*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)+B(IC+
     *  I))))
     *   -S2*((B(IA+I)-haf*(B(IB+I)-B(IC+I)))-(SIN60*(A(IB+I)-A(IC+
     *  I))))
      D(JC+J)=
     *    S2*((A(IA+I)-haf*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)+B(IC+
     *  I))))
     *   +C2*((B(IA+I)-haf*(B(IB+I)-B(IC+I)))-(SIN60*(A(IB+I)-A(IC+
     *  I))))
      I=I+INC3
      J=J+INC4
  330 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  340 CONTINUE
      IA=IA+IINK
      IB=IB+IINK
      IC=IC-IINK
      JBASE=JBASE+JUMP
  350 CONTINUE
      IF (IA.GT.IC) GO TO 900
  360 CONTINUE
      IBASE=0
      DO 380 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 370 IJK=1,ILOT
      C(JA+J)=A(IA+I)+A(IB+I)
      C(JB+J)=(haf*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))
      C(JC+J)=-(haf*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))
      I=I+INC3
      J=J+INC4
  370 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  380 CONTINUE
      GO TO 900
C
  390 CONTINUE
      SSIN60=two*SIN60
      DO 394 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 392 IJK=1,ILOT
      C(JA+J)=two*(A(IA+I)+A(IB+I))
      C(JB+J)=(two*A(IA+I)-A(IB+I))-(SSIN60*B(IB+I))
      C(JC+J)=(two*A(IA+I)-A(IB+I))+(SSIN60*B(IB+I))
      I=I+INC3
      J=J+INC4
  392 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  394 CONTINUE
      GO TO 900
C
C     CODING FOR FACTOR 4
C     -------------------
  400 CONTINUE
      IA=1
      IB=IA+(2*M-ILA)*INC1
      IC=IB+2*M*INC1
      ID=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
C
      IF (ILA.EQ.M) GO TO 490
C
      DO 420 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 410 IJK=1,ILOT
      C(JA+J)=(A(IA+I)+A(IC+I))+A(IB+I)
      C(JB+J)=(A(IA+I)-A(IC+I))-B(IB+I)
      C(JC+J)=(A(IA+I)+A(IC+I))-A(IB+I)
      C(JD+J)=(A(IA+I)-A(IC+I))+B(IB+I)
      I=I+INC3
      J=J+INC4
  410 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  420 CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC-IINK
      ID=ID-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IB.EQ.IC) GO TO 460
      DO 450 K=ILA,KSTOP,ILA
      KB=K+K
      KC=KB+KB
      KD=KC+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      IBASE=0
      DO 440 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 430 IJK=1,ILOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)-B(IC+I))+(B(IB+I)-B(ID+I))
      C(JC+J)=
     *    C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     *   -S2*((B(IA+I)-B(IC+I))-(B(IB+I)-B(ID+I)))
      D(JC+J)=
     *    S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     *   +C2*((B(IA+I)-B(IC+I))-(B(IB+I)-B(ID+I)))
      C(JB+J)=
     *    C1*((A(IA+I)-A(IC+I))-(B(IB+I)+B(ID+I)))
     *   -S1*((B(IA+I)+B(IC+I))+(A(IB+I)-A(ID+I)))
      D(JB+J)=
     *    S1*((A(IA+I)-A(IC+I))-(B(IB+I)+B(ID+I)))
     *   +C1*((B(IA+I)+B(IC+I))+(A(IB+I)-A(ID+I)))
      C(JD+J)=
     *    C3*((A(IA+I)-A(IC+I))+(B(IB+I)+B(ID+I)))
     *   -S3*((B(IA+I)+B(IC+I))-(A(IB+I)-A(ID+I)))
      D(JD+J)=
     *    S3*((A(IA+I)-A(IC+I))+(B(IB+I)+B(ID+I)))
     *   +C3*((B(IA+I)+B(IC+I))-(A(IB+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  430 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  440 CONTINUE
      IA=IA+IINK
      IB=IB+IINK
      IC=IC-IINK
      ID=ID-IINK
      JBASE=JBASE+JUMP
  450 CONTINUE
      IF (IB.GT.IC) GO TO 900
  460 CONTINUE
      IBASE=0
      SIN45=SQRT(haf)
      DO 480 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 470 IJK=1,ILOT
      C(JA+J)=A(IA+I)+A(IB+I)
      C(JB+J)=SIN45*((A(IA+I)-A(IB+I))-(B(IA+I)+B(IB+I)))
      C(JC+J)=B(IB+I)-B(IA+I)
      C(JD+J)=-SIN45*((A(IA+I)-A(IB+I))+(B(IA+I)+B(IB+I)))
      I=I+INC3
      J=J+INC4
  470 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  480 CONTINUE
      GO TO 900
C
  490 CONTINUE
      DO 494 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 492 IJK=1,ILOT
      C(JA+J)=two*((A(IA+I)+A(IC+I))+A(IB+I))
      C(JB+J)=two*((A(IA+I)-A(IC+I))-B(IB+I))
      C(JC+J)=two*((A(IA+I)+A(IC+I))-A(IB+I))
      C(JD+J)=two*((A(IA+I)-A(IC+I))+B(IB+I))
      I=I+INC3
      J=J+INC4
  492 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  494 CONTINUE
      GO TO 900
C
C     CODING FOR FACTOR 5
C     -------------------
  500 CONTINUE
      IA=1
      IB=IA+(2*M-ILA)*INC1
      IC=IB+2*M*INC1
      ID=IC
      IE=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
      JE=JD+JINK
C
      IF (ILA.EQ.M) GO TO 590
C
      DO 520 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 510 IJK=1,ILOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      C(JB+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))+QRT5*(A(IB+I)-A(IC+
     *  I)))
     *    -(SIN72*B(IB+I)+SIN36*B(IC+I))
      C(JC+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))-QRT5*(A(IB+I)-A(IC+
     *  I)))
     *    -(SIN36*B(IB+I)-SIN72*B(IC+I))
      C(JD+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))-QRT5*(A(IB+I)-A(IC+
     *  I)))
     *    +(SIN36*B(IB+I)-SIN72*B(IC+I))
      C(JE+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))+QRT5*(A(IB+I)-A(IC+
     *  I)))
     *    +(SIN72*B(IB+I)+SIN36*B(IC+I))
      I=I+INC3
      J=J+INC4
  510 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  520 CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC+IINK
      ID=ID-IINK
      IE=IE-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IB.EQ.ID) GO TO 560
      DO 550 K=ILA,KSTOP,ILA
      KB=K+K
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      IBASE=0
      DO 540 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 530 IJK=1,ILOT
C
      A10(IJK)=(A(IA+I)-0.25*((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))))
     *    +QRT5*((A(IB+I)+A(IE+I))-(A(IC+I)+A(ID+I)))
      A20(IJK)=(A(IA+I)-0.25*((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))))
     *    -QRT5*((A(IB+I)+A(IE+I))-(A(IC+I)+A(ID+I)))
      B10(IJK)=(B(IA+I)-0.25*((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I))))
     *    +QRT5*((B(IB+I)-B(IE+I))-(B(IC+I)-B(ID+I)))
      B20(IJK)=(B(IA+I)-0.25*((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I))))
     *    -QRT5*((B(IB+I)-B(IE+I))-(B(IC+I)-B(ID+I)))
      A11(IJK)=SIN72*(B(IB+I)+B(IE+I))+SIN36*(B(IC+I)+B(ID+I))
      A21(IJK)=SIN36*(B(IB+I)+B(IE+I))-SIN72*(B(IC+I)+B(ID+I))
      B11(IJK)=SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))
      B21(IJK)=SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))
C
      C(JA+J)=A(IA+I)+((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I)))
      D(JA+J)=B(IA+I)+((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I)))
      C(JB+J)=C1*(A10(IJK)-A11(IJK))-S1*(B10(IJK)+B11(IJK))
      D(JB+J)=S1*(A10(IJK)-A11(IJK))+C1*(B10(IJK)+B11(IJK))
      C(JE+J)=C4*(A10(IJK)+A11(IJK))-S4*(B10(IJK)-B11(IJK))
      D(JE+J)=S4*(A10(IJK)+A11(IJK))+C4*(B10(IJK)-B11(IJK))
      C(JC+J)=C2*(A20(IJK)-A21(IJK))-S2*(B20(IJK)+B21(IJK))
      D(JC+J)=S2*(A20(IJK)-A21(IJK))+C2*(B20(IJK)+B21(IJK))
      C(JD+J)=C3*(A20(IJK)+A21(IJK))-S3*(B20(IJK)-B21(IJK))
      D(JD+J)=S3*(A20(IJK)+A21(IJK))+C3*(B20(IJK)-B21(IJK))
C
      I=I+INC3
      J=J+INC4
  530 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  540 CONTINUE
      IA=IA+IINK
      IB=IB+IINK
      IC=IC+IINK
      ID=ID-IINK
      IE=IE-IINK
      JBASE=JBASE+JUMP
  550 CONTINUE
      IF (IB.GT.ID) GO TO 900
  560 CONTINUE
      IBASE=0
      DO 580 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 570 IJK=1,ILOT
      C(JA+J)=(A(IA+I)+A(IB+I))+A(IC+I)
      C(JB+J)=(QRT5*(A(IA+I)-A(IB+I))+(0.25*(A(IA+I)+A(IB+I))-A(IC+
     *  I)))
     *    -(SIN36*B(IA+I)+SIN72*B(IB+I))
      C(JE+J)=-(QRT5*(A(IA+I)-A(IB+I))+(0.25*(A(IA+I)+A(IB+I))-A(IC+
     *  I)))
     *    -(SIN36*B(IA+I)+SIN72*B(IB+I))
      C(JC+J)=(QRT5*(A(IA+I)-A(IB+I))-(0.25*(A(IA+I)+A(IB+I))-A(IC+
     *  I)))
     *    -(SIN72*B(IA+I)-SIN36*B(IB+I))
      C(JD+J)=-(QRT5*(A(IA+I)-A(IB+I))-(0.25*(A(IA+I)+A(IB+I))-A(IC+
     *  I)))
     *    -(SIN72*B(IA+I)-SIN36*B(IB+I))
      I=I+INC3
      J=J+INC4
  570 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  580 CONTINUE
      GO TO 900
C
  590 CONTINUE
      QQRT5=two*QRT5
      SSIN36=two*SIN36
      SSIN72=two*SIN72
      DO 594 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 592 IJK=1,ILOT
      C(JA+J)=two*(A(IA+I)+(A(IB+I)+A(IC+I)))
      C(JB+J)=(two*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))
     *    +QQRT5*(A(IB+I)-A(IC+I)))-(SSIN72*B(IB+I)+SSIN36*B(IC+I))
      C(JC+J)=(two*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))
     *    -QQRT5*(A(IB+I)-A(IC+I)))-(SSIN36*B(IB+I)-SSIN72*B(IC+I))
      C(JD+J)=(two*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))
     *    -QQRT5*(A(IB+I)-A(IC+I)))+(SSIN36*B(IB+I)-SSIN72*B(IC+I))
      C(JE+J)=(two*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))
     *    +QQRT5*(A(IB+I)-A(IC+I)))+(SSIN72*B(IB+I)+SSIN36*B(IC+I))
      I=I+INC3
      J=J+INC4
  592 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  594 CONTINUE
      GO TO 900
C
C     CODING FOR FACTOR 6
C     -------------------
  600 CONTINUE
      IA=1
      IB=IA+(2*M-ILA)*INC1
      IC=IB+2*M*INC1
      ID=IC+2*M*INC1
      IE=IC
      IF=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
      JE=JD+JINK
      JF=JE+JINK
C
      IF (ILA.EQ.M) GO TO 690
C
      DO 620 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 610 IJK=1,ILOT
      C(JA+J)=(A(IA+I)+A(ID+I))+(A(IB+I)+A(IC+I))
      C(JD+J)=(A(IA+I)-A(ID+I))-(A(IB+I)-A(IC+I))
      C(JB+J)=((A(IA+I)-A(ID+I))+haf*(A(IB+I)-A(IC+I)))
     *    -(SIN60*(B(IB+I)+B(IC+I)))
      C(JF+J)=((A(IA+I)-A(ID+I))+haf*(A(IB+I)-A(IC+I)))
     *    +(SIN60*(B(IB+I)+B(IC+I)))
      C(JC+J)=((A(IA+I)+A(ID+I))-haf*(A(IB+I)+A(IC+I)))
     *    -(SIN60*(B(IB+I)-B(IC+I)))
      C(JE+J)=((A(IA+I)+A(ID+I))-haf*(A(IB+I)+A(IC+I)))
     *    +(SIN60*(B(IB+I)-B(IC+I)))
      I=I+INC3
      J=J+INC4
  610 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  620 CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC+IINK
      ID=ID-IINK
      IE=IE-IINK
      IF=IF-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IC.EQ.ID) GO TO 660
      DO 650 K=ILA,KSTOP,ILA
      KB=K+K
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      KF=KE+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      C5=TRIGS(KF+1)
      S5=TRIGS(KF+2)
      IBASE=0
      DO 640 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 630 IJK=1,ILOT
C
      A11(IJK)= (A(IE+I)+A(IB+I))+(A(IC+I)+A(IF+I))
      A20(IJK)=(A(IA+I)+A(ID+I))-haf*A11(IJK)
      A21(IJK)=SIN60*((A(IE+I)+A(IB+I))-(A(IC+I)+A(IF+I)))
      B11(IJK)= (B(IB+I)-B(IE+I))+(B(IC+I)-B(IF+I))
      B20(IJK)=(B(IA+I)-B(ID+I))-haf*B11(IJK)
      B21(IJK)=SIN60*((B(IB+I)-B(IE+I))-(B(IC+I)-B(IF+I)))
C
      C(JA+J)=(A(IA+I)+A(ID+I))+A11(IJK)
      D(JA+J)=(B(IA+I)-B(ID+I))+B11(IJK)
      C(JC+J)=C2*(A20(IJK)-B21(IJK))-S2*(B20(IJK)+A21(IJK))
      D(JC+J)=S2*(A20(IJK)-B21(IJK))+C2*(B20(IJK)+A21(IJK))
      C(JE+J)=C4*(A20(IJK)+B21(IJK))-S4*(B20(IJK)-A21(IJK))
      D(JE+J)=S4*(A20(IJK)+B21(IJK))+C4*(B20(IJK)-A21(IJK))
C
      A11(IJK)=(A(IE+I)-A(IB+I))+(A(IC+I)-A(IF+I))
      B11(IJK)=(B(IE+I)+B(IB+I))-(B(IC+I)+B(IF+I))
      A20(IJK)=(A(IA+I)-A(ID+I))-haf*A11(IJK)
      A21(IJK)=SIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
      B20(IJK)=(B(IA+I)+B(ID+I))+haf*B11(IJK)
      B21(IJK)=SIN60*((B(IE+I)+B(IB+I))+(B(IC+I)+B(IF+I)))
C
      C(JD+J)=
     *  C3*((A(IA+I)-A(ID+I))+A11(IJK))-S3*((B(IA+I)+B(ID+I))-B11(IJK))
      D(JD+J)=
     *  S3*((A(IA+I)-A(ID+I))+A11(IJK))+C3*((B(IA+I)+B(ID+I))-B11(IJK))
      C(JB+J)=C1*(A20(IJK)-B21(IJK))-S1*(B20(IJK)-A21(IJK))
      D(JB+J)=S1*(A20(IJK)-B21(IJK))+C1*(B20(IJK)-A21(IJK))
      C(JF+J)=C5*(A20(IJK)+B21(IJK))-S5*(B20(IJK)+A21(IJK))
      D(JF+J)=S5*(A20(IJK)+B21(IJK))+C5*(B20(IJK)+A21(IJK))
C
      I=I+INC3
      J=J+INC4
  630 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  640 CONTINUE
      IA=IA+IINK
      IB=IB+IINK
      IC=IC+IINK
      ID=ID-IINK
      IE=IE-IINK
      IF=IF-IINK
      JBASE=JBASE+JUMP
  650 CONTINUE
      IF (IC.GT.ID) GO TO 900
  660 CONTINUE
      IBASE=0
      DO 680 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 670 IJK=1,ILOT
      C(JA+J)=A(IB+I)+(A(IA+I)+A(IC+I))
      C(JD+J)=B(IB+I)-(B(IA+I)+B(IC+I))
      C(JB+J)=(SIN60*(A(IA+I)-A(IC+I)))-(haf*(B(IA+I)+B(IC+I))+B(IB+
     *  I))
      C(JF+J)=-(SIN60*(A(IA+I)-A(IC+I)))-(haf*(B(IA+I)+B(IC+I))+B(IB+
     *  I))
      C(JC+J)=SIN60*(B(IC+I)-B(IA+I))+(haf*(A(IA+I)+A(IC+I))-A(IB+I))
      C(JE+J)=SIN60*(B(IC+I)-B(IA+I))-(haf*(A(IA+I)+A(IC+I))-A(IB+I))
      I=I+INC3
      J=J+INC4
  670 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  680 CONTINUE
      GO TO 900
C
  690 CONTINUE
      SSIN60=two*SIN60
      DO 694 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 692 IJK=1,ILOT
      C(JA+J)=(two*(A(IA+I)+A(ID+I)))+(two*(A(IB+I)+A(IC+I)))
      C(JD+J)=(two*(A(IA+I)-A(ID+I)))-(two*(A(IB+I)-A(IC+I)))
      C(JB+J)=(two*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I)))
     *    -(SSIN60*(B(IB+I)+B(IC+I)))
      C(JF+J)=(two*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I)))
     *    +(SSIN60*(B(IB+I)+B(IC+I)))
      C(JC+J)=(two*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I)))
     *    -(SSIN60*(B(IB+I)-B(IC+I)))
      C(JE+J)=(two*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I)))
     *    +(SSIN60*(B(IB+I)-B(IC+I)))
      I=I+INC3
      J=J+INC4
  692 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  694 CONTINUE
      GO TO 900
C
C     CODING FOR FACTOR 8
C     -------------------
  800 CONTINUE
      IBAD=3
      IF (ILA.NE.M) GO TO 910
      IA=1
      IB=IA+ILA*INC1
      IC=IB+2*ILA*INC1
      ID=IC+2*ILA*INC1
      IE=ID+2*ILA*INC1
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
      JE=JD+JINK
      JF=JE+JINK
      JG=JF+JINK
      JH=JG+JINK
      SSIN45=SQRT(two)
C
      DO 820 IL=1,ILA
      I=IBASE
      J=JBASE
C==DIR$ IVDEP
*==VOCL LOOP,NOVREC
      DO 810 IJK=1,ILOT
      C(JA+J)=two*(((A(IA+I)+A(IE+I))+A(IC+I))+(A(IB+I)+A(ID+I)))
      C(JE+J)=two*(((A(IA+I)+A(IE+I))+A(IC+I))-(A(IB+I)+A(ID+I)))
      C(JC+J)=two*(((A(IA+I)+A(IE+I))-A(IC+I))-(B(IB+I)-B(ID+I)))
      C(JG+J)=two*(((A(IA+I)+A(IE+I))-A(IC+I))+(B(IB+I)-B(ID+I)))
      C(JB+J)=two*((A(IA+I)-A(IE+I))-B(IC+I))
     *    +SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
      C(JF+J)=two*((A(IA+I)-A(IE+I))-B(IC+I))
     *    -SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
      C(JD+J)=two*((A(IA+I)-A(IE+I))+B(IC+I))
     *    -SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
      C(JH+J)=two*((A(IA+I)-A(IE+I))+B(IC+I))
     *    +SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
      I=I+INC3
      J=J+INC4
  810 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  820 CONTINUE
C
C     RETURN
C     ------
  900 CONTINUE
      IBAD=0
  910 CONTINUE
      IERR=IBAD
      RETURN
      END
C Copyright 1981-2007 ECMWF
C 
C Licensed under the GNU Lesser General Public License which
C incorporates the terms and conditions of version 3 of the GNU
C General Public License.
C See LICENSE and gpl-3.0.txt for details.
C

      SUBROUTINE SET99_M8(TRIGS,IFAX,N)
	implicit none
	integer N, IFAX(N)
	real  *8 TRIGS(N)
      INTEGER JFAX(10),LFAX(7)
	integer ixxx, nil,nhl,k,nu,ifac,l,nfax,i
	real  *8 del, angle
C
C     SUBROUTINE 'SET99' - COMPUTES FACTORS OF N & TRIGONOMETRIC
C     FUNCTIONS REQUIRED BY FFT99 & FFT991
C
      DATA LFAX/6,8,5,4,3,2,1/
      IXXX=1
C
      DEL=4.0*ASIN(1.0D0)/FLOAT(N)
      NIL=0
      NHL=(N/2)-1
      DO 10 K=NIL,NHL
      ANGLE=FLOAT(K)*DEL
      TRIGS(2*K+1)=COS(ANGLE)
      TRIGS(2*K+2)=SIN(ANGLE)
   10 CONTINUE
C
C     FIND FACTORS OF N (8,6,5,4,3,2; ONLY ONE 8 ALLOWED)
C     LOOK FOR SIXES FIRST, STORE FACTORS IN DESCENDING ORDER
      NU=N
      IFAC=6
      K=0
      L=1
   20 CONTINUE
      IF (MOD(NU,IFAC).NE.0) GO TO 30
      K=K+1
      JFAX(K)=IFAC
      IF (IFAC.NE.8) GO TO 25
      IF (K.EQ.1) GO TO 25
      JFAX(1)=8
      JFAX(K)=6
   25 CONTINUE
      NU=NU/IFAC
      IF (NU.EQ.1) GO TO 50
      IF (IFAC.NE.8) GO TO 20
   30 CONTINUE
      L=L+1
      IFAC=LFAX(L)
      IF (IFAC.GT.1) GO TO 20
C
      WRITE(6,40) N
   40 FORMAT('1N =',I4,' - CONTAINS ILLEGAL FACTORS')
      RETURN
C
C     NOW REVERSE ORDER OF FACTORS
   50 CONTINUE
      NFAX=K
      IFAX(1)=NFAX
      DO 60 I=1,NFAX
      IFAX(NFAX+2-I)=JFAX(I)
   60 CONTINUE
      IFAX(10)=N
      RETURN
      END
