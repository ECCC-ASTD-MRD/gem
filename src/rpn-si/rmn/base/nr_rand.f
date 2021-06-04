*/* RMNLIB - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
* *                          Environnement Canada
* *
* * This library is free software; you can redistribute it and/or
* * modify it under the terms of the GNU Lesser General Public
* * License as published by the Free Software Foundation,
* * version 2.1 of the License.
* *
* * This library is distributed in the hope that it will be useful,
* * but WITHOUT ANY WARRANTY; without even the implied warranty of
* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* * Lesser General Public License for more details.
* *
* * You should have received a copy of the GNU Lesser General Public
* * License along with this library; if not, write to the
* * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
* * Boston, MA 02111-1307, USA.
* */
      INTEGER FUNCTION nr_rand_seed(idum)
      implicit none
      integer nr_rand_i
      integer idum
      SAVE iseed
      INTEGER iseed
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      DATA iseed /-1/

      iseed=idum
      nr_rand_seed=idum
      return

      entry nr_rand_i()

      if(iseed.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(iseed)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        iseed=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      nr_rand_i=mj
      return
      end
      real function nr_rand_r()
      implicit none
      integer nr_rand_i
      external nr_rand_i
      real *8 FAC,MBIG
      PARAMETER (MBIG=1000000000,FAC=1./MBIG)
      nr_rand_r=FAC*nr_rand_i()
      return
      END
      real*8 function nr_rand_d()
      implicit none
      integer nr_rand_i
      external nr_rand_i
      real *8 FAC,MBIG
      PARAMETER (MBIG=1000000000,FAC=1./MBIG)
      nr_rand_d=FAC*nr_rand_i()
      return
      END
