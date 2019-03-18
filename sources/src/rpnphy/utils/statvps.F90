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
!**s/r sttvps - calcule la moyenne, la variance, le minimum et
!               le maximum de toutes les variables du bus busnom

      subroutine statvps (vp,kount,trnch,from,ni,nk,busnom)
      implicit none
#include <arch_specific.hf>
!
      integer ni,nk,kount,trnch
      character*(*) from
      character*1 busnom

      real vp(*)
!
!Author
!         Robert Benoit (Aug 93)
!
!Revision
! 001     B. Bilodeau (Feb 96) - Revised physics interface
! 002     B. Bilodeau (Nov 98) - Volatile bus diagnostics
! 003     B. Bilodeau (Feb 99) - Entry and dynamics buses diagnostics
! 004     B. Bilodeau (Nov 04) - Change format of output
! 005     B. Bilodeau (Jan 06) - Mosaic
! 006     V. Lee (Mar 2011) - mosaic pointer = dynpar(il,8)+1
!
!Object
!     calculates and  prints : the average  (moy)
!                              the variance (var)
!                              the minimum and the maximum of vp
!Arguments
!
!         - Input -
! vp      stack of the permanent variables of the physics
! no      counter
! from    name of the calling module
! ni      1st horizontal dimension of the grid
! nk      vertical dimension of the grid
! busnom  'P' : permanent bus
!         'V' : volatile  bus
!         'E' : entry     bus
!         'D' : dynamics  bus
!
!
!
!Implicits
!
!     to handle the list of vp names
!
#include "buses.cdk"
!
!*
      integer i,k,top
      real sum,moy,var,vpmin,vpmax
      integer imin,kmin,imax,kmax
      integer il, siz, esp, i0, m, mosaik, mul, stride
      character*1 busnomc
!
!--------------------------------------------------------------------
!
!     loop on the VP  list
!
!     conversion from lower case to upper case
      call low2up(busnom,busnomc)
!
      if      (busnomc.eq.'P') then
         top = pertop
      else if (busnomc.eq.'V') then
         top = voltop
      else if (busnomc.eq.'E') then
         top = enttop
      else if (busnomc.eq.'D') then
         top = dyntop
      endif
!
      do 100 il=1,top
!
         if      (busnomc.eq.'P') then
            siz   =perpar(il,5)
            mul   =perpar(il,6)
            mosaik=perpar(il,8)+1
         else if (busnomc.eq.'V') then
            siz   =volpar(il,5)
            mul   =volpar(il,6)
            mosaik=volpar(il,8)+1
         else if (busnomc.eq.'E') then
            siz   =entpar(il,5)
            mul   =entpar(il,6)
            mosaik=entpar(il,8)+1
         else if (busnomc.eq.'D') then
            siz   =dynpar(il,5)
            mul   =dynpar(il,6)
            mosaik=dynpar(il,8)+1
         endif
!
         if( siz.eq.ni) then
            esp = 1
         else
            esp = 2
         endif
!
         if      (busnomc.eq.'P') then
            i0=perpar(il,1)-1  !element "0"
         else if (busnomc.eq.'V') then
            i0=volpar(il,1)-1  !element "0"
         else if (busnomc.eq.'E') then
            i0=entpar(il,1)-1  !element "0"
         else if (busnomc.eq.'D') then
            i0=dynpar(il,1)-1  !element "0"
         endif
!
!
         do 110 m=1,mul*mosaik
!
!     ** On calcule la moyenne.
!
!        "stride" est utilise seulement si mul*mosaik > 1
         stride = (m-1)*siz
!
         sum = 0.0
         do 1 i=1,siz
            sum = sum + vp(i+i0+stride)
 1       continue
         moy = sum / float(siz)
!
!     ** On calcule la variance
!
            sum = 0.0
               do 2 i=1,siz
                  sum = sum + ((vp(i+i0+stride) - moy)*(vp(i+i0+stride) - moy))
 2             continue
               var = sqrt (sum / float(siz))
!
!     ** On identifie le minimum et le maximum.
!
               imin = 1
               kmin = 0
               imax = 1
               kmax = 0
               vpmax = vp(i0+1+stride)
               vpmin = vp(i0+1+stride)
!
                  do 3 i=1,siz
                     if (vp(i+i0+stride) .gt. vpmax) then
                        vpmax  = vp(i+i0+stride)
                        imax = i
!                       kmax = k
                     endif
                     if (vp(i+i0+stride) .lt. vpmin) then
                        vpmin  = vp(i+i0+stride)
                        imin = i
!                       kmin = k
                     endif
 3                continue
!     compute kmin/max if needed
                  if (esp.eq.2) then
!     min
                     i=mod(imin,ni)
                     if (i.eq.0) i=ni
                     k=1+(imin-i)/ni
                     imin=i
                     kmin=k
!     max
                     i=mod(imax,ni)
                     if (i.eq.0) i=ni
                     k=1+(imax-i)/ni
                     imax=i
                     kmax=k
                  else
                  endif
!
!     ** On imprime
!
         if      (busnomc.eq.'P') then
                  write(6,1000) kount,trnch,from,m,pernm(il,1),moy,var,imin,kmin,vpmin, &
                       imax,kmax,vpmax
         else if (busnomc.eq.'V') then
                  write(6,1000) kount,trnch,from,m,volnm(il,1),moy,var,imin,kmin,vpmin, &
                       imax,kmax,vpmax
         else if (busnomc.eq.'E') then
                  write(6,1000) kount,trnch,from,m,entnm(il,1),moy,var,imin,kmin,vpmin, &
                       imax,kmax,vpmax
         else if (busnomc.eq.'D') then
                  write(6,1000) kount,trnch,from,m,dynnm(il,1),moy,var,imin,kmin,vpmin, &
                       imax,kmax,vpmax
         endif
!
 110  continue
!
 100  continue
!
 1000             format (2i4,a10,i2,' ',a7,' Mean:',e15.8,'  Var:',e15.8, &
                       '  Min:[(',i3,',',i3,') ', &
                       e15.8,']',' Max:[(',i3,',',i3,') ', &
                       e15.8,']')
!
!----------------------------------------------------------------
      return
      end
