!-------------------------------------- LICENCE BEGIN -------------------------
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
!-------------------------------------- LICENCE END ---------------------------

subroutine qozon3 (oz, ozotoit,ozzx, press, s,lev,nk,np,nmax, &
     lref,pref, x, xm, y, ym, is,ism)
   use tdpack_const
   implicit none
!!!#include <arch_specific.hf>
   integer i,k,l,ii
   integer int,lev,lref,lstart,nk,nmax,np
   real oz(nmax,lev),ozzx(nmax,lref),press(np),s(np,nk+1)
   logical lo1
   real ac1, ac2, ac3, pl, slope, tercep, zz
   real x(np), xm(np), y(np), ym(np)
   integer is(np), ism(np)
   real pref(lref)
   real ozotoit(NMAX)

#include "tables.cdk"
!
!Author
!          L.Garand (June 1989)
!Revision
! 001      see version 5.5.0 for previous history
!
!Object
!          to vertically interpolate the climatological fields
!          OZZX from OZOREF2 to give the ozone mixing ratio (kg/kg)
!          at the centre of each layer.
!
!Arguments
!
!          - Output -
! oz       ozone mixing ratio in kg/kg at the center of each
!          layer
! ozotoit   total ozone (cm stp) above model roof
!
!          - Input -
! ozzx     ozone mixing ratio (kg/kg) at
!          the 37 levels of climatological files
! press    surface pressure (N/M**2)
! s        sigma levels at the centre of which OZ will be produced
! lev      number of sigma levels
! nk       number of layers (NK=LEV-1)
! np       number of points to calculate
! nmax     maximum number of points
! lref     number of climatological levels
! pref     ozone climatological pressures
! x        work field
! xm       work field
! y        work field
! ym       work field
! is       work field
! ism      work field
!
!Notes
!      Original code provided by J.P. Blanchet, CCRN
!      interpolation in mixing ratio now done; originally
!      interpolation made on integrated ozone amount
!*
!
      real fact
!
!----------------------------------------------------------------------
!                       >>> Interpolation <<<
!
!
!........ Locate indices for interpolation at full pressure levels.

      lstart = 2
      do 180 l = 1, lev
          do 120 i = 1,np
              is(i) = 2
  120     continue

          do 140 k = lstart, lref
              do 130 i = 1,np
              ac1=float(k)
              ac2=float(is(i))
              pl=s(i,l)*press(i)
              lo1=pref(k).ge.pl.and.is(i).le.2
              if (lo1) then
                 ac3=ac1
              else
                 ac3=ac2
              endif
              is(i)=int(ac3)
              ac1=float(lref)
              lo1=(k.eq.lref.and.pl.gt.pref(2)).and.is(i).eq.2
              if (lo1) then
                 ac2=ac1
              else
                 ac2=ac3
              endif
              is(i)=ac2
  130         continue
  140     continue

          lstart = lref
          do 150 i = 1,np
              ism(i) = is(i) - 1
              lstart = min0 (lstart, is(i))
  150     continue

!
       do 350 i=1,np
      x(i)=pref(is(i))
      xm(i)=pref(ism(i))
      y(i)=ozzx(i,is(i))
      ym(i)=ozzx(i,ism(i))
 350  continue

          do 170 i = 1,np
              slope   = (y(i) - ym(i)) / (x(i) - xm(i))
              tercep  = ym(i) - slope * xm(i)
              zz=amax1(s(i,l)*press(i), pref(ism(i)))
              oz(i,l) = slope * zz + tercep
  170     continue
  180 continue
!
!     approximate ozone (cm stp) above model roof for solar code
      fact = 1./2.144e-2
      ii=lref-1
      do 185 i = 1,np
         ozotoit(i) = 0.
         do 186 l=1,ii
            if (pref(l+1).lt.std_p_prof(1)) then
               zz= (pref(l+1)-pref(l))/grav* &
                0.5* (ozzx(i,l)+ozzx(i,l+1)) * fact
                ozotoit(i)=ozotoit(i)+zz
            endif
 186     continue
 185   continue
!
      return
      end
