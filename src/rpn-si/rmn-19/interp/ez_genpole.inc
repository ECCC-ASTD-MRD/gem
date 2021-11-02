!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
!********************************************************************
!**
!**
      subroutine ez_genpole(vpolnor,vpolsud,fld,ni,nj,vecteur,grtyp,hem,      x,y,z,lat,lon,glat,ordint)
      implicit none
      real vpolnor,vpolsud
      integer ni,nj
      real fld(ni,nj)
      character*1 grtyp       
      logical vecteur   
      integer hem,ordint
      real x(2,ni),y(2,ni),z(2,ni),lat(2,ni),lon(2,ni),glat(*)

#include "pi.cdk"
#include "qqqpar.cdk"
#include "ez_qqqxtrp.cdk"
#include "ez_gaussgd.cdk"

      integer i,j
      real*8 x2,y2,xy,x3,y3,x2y,xy2
      real pi,pj,dlon,temp
      real*8 a(6,6),b(6),vpole,det,deta,detb,detc
      integer l(6),m(6)
      integer n,nglat,ier
      real lat1,lat2,dlat,xla0
      integer ilat1,ilat2

!      if (vecteur) then
!         return
!      endif
      
      
      if (grtyp.ne.'G'.and.grtyp.ne.'A') then
         print *
         print *,     ' <GENPOLE> BAD GRID TYPE'
         print *,     '           THIS FUNCTION IS DESIGNED ONLY FOR A AND G GRIDS'
         print *
         goto 999
      endif

      if (hem.ne. global) then
         nglat = nj*2
      else
         nglat = nj
      endif

      if (grtyp.eq.'G') then
         call ez_glat(glat,lroots,nglat,hem)
      else
         dlat = 180./float(nglat)
         xla0 = -90. + 0.5*dlat
         do j=1,nglat
            glat(j)=xla0+(j-1)*dlat
         enddo
      endif
      
      do n=1,2
         if (n.eq.nord) then
            lat2  = glat(nglat-1)
            lat1  = glat(nglat)
            if (hem.eq.global) then
               ilat1 = nj
               ilat2 = nj-1
            else if (hem.eq.nord) then
               ilat1 = nj
               ilat2 = nj-1
            else
               ilat1 = 1
               ilat2 = 2
            endif
         else
            lat1  = glat(1)
            lat2  = glat(2)
            if (hem.eq.global) then
               ilat1 = 1
               ilat2 = 2
            else if (hem.eq.nord) then
               ilat1 = nj
               ilat2 = nj-1
            else
               ilat1 = 1
               ilat2 = 2
            endif
         endif

         dlon = 360. / float(ni)
         do i=1,ni
            lat(1,i)=lat1
            lat(2,i)=lat2
            lon(1,i)=(i-1)*dlon
            lon(2,i)=lon(1,i)
            z(1,i)=fld(i,ilat1)
            z(2,i)=fld(i,ilat2)
         enddo

         pi = 0.0
         pj = 0.0
         call ez_vxyfll(x,y,lat,lon,2*ni,1000.0,0.0,pi,pj,n)
         
         do j=1,6
            b(j)=0.0
            do i=1,6
               a(i,j)=0.0
            enddo
         enddo
         
         do j=1,2
            do i=1,ni
               xy=x(j,i)*y(j,i)
               x2=x(j,i)*x(j,i)
               x3=x2*x(j,i)
               y2=y(j,i)*y(j,i)
               y3=y2*y(j,i)
               x2y=x2*y(j,i)
               xy2=x(j,i)*y2
               
               a(1,1)=a(1,1) + 1
               a(2,1)=a(2,1) + x(j,i)
               a(3,1)=a(3,1) + y(j,i)
               a(4,1)=a(4,1) + x2
               a(5,1)=a(5,1) + xy
               a(6,1)=a(6,1) + y2
               
               a(2,2)=a(2,2) + x2
               a(3,2)=a(3,2) + xy
               a(4,2)=a(4,2) + x3
               a(5,2)=a(5,2) + x2y
               a(6,2)=a(6,2) + xy2
               
               a(3,3)=a(3,3) + y2
               a(4,3)=a(4,3) + x2y
               a(5,3)=a(5,3) + xy2
               a(6,3)=a(6,3) + y3
               
               a(4,4)=a(4,4) + x2 * x2
               a(5,4)=a(5,4) + x3 * y(j,i)
               a(6,4)=a(6,4) + x2 * y2
               
               a(5,5)=a(5,5) + x2 * y2
               a(6,5)=a(6,5) + x(j,i) * y3
               
               a(6,6)=a(6,6) + y2*y2
               
               b(1) = b(1) + z(j,i)
               b(2) = b(2) + z(j,i)*x(j,i)
               b(3) = b(3) + z(j,i)*y(j,i)
               b(4) = b(4) + z(j,i)*x2
               b(5) = b(5) + z(j,i)*xy
               b(6) = b(6) + z(j,i)*y2
            enddo   
         enddo
         
         do i=1,6
            do j=i+1,6
               a(i,j)=a(j,i)
            enddo
         enddo
         
         if (ordint.eq.cubique) then
            call ez_mtxinv8(a,l,m,6)
            
            vpole = 0.0 
            do i=1,6
               vpole= vpole + a(i,1)*b(i)
            enddo
         endif

         if (ordint.eq.lineair) then
            det  = a(1,1)*(a(2,2)*a(3,3)-a(3,2)*a(2,3))            - a(2,1)*(a(1,2)*a(3,3)-a(3,2)*a(1,3))            + a(3,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
            deta = b(1)*(a(2,2)*a(3,3)-a(3,2)*a(2,3))            - a(2,1)*(b(2)*a(3,3)-a(3,2)*b(3))            + a(3,1)*(b(2)*a(2,3)-a(2,2)*b(3))
            detb = a(1,1)*(b(2)*a(3,3)-a(3,2)*b(3))            - b(1)*(a(1,2)*a(3,3)-a(3,2)*a(1,3))            + a(3,1)*(a(1,2)*b(3)-b(2)*a(1,3))
            detc = a(1,1)*(a(2,2)*b(3)-b(2)*a(2,3))            - a(2,1)*(a(1,2)*b(3)-b(2)*a(1,3))            + b(1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
            
            if (det.ne.0.0) then
               deta = deta/det
               detb = detb/det
               detc = detc/det
            endif
            
            if (det.ne.0.0) then
               vpole= deta
            else
               vpole = b(1)/a(1,1)
            endif
         endif
         
         if (ordint.eq.voisin) then
            vpole = b(1)/a(1,1)
         endif
         
         if (n.eq.nord) then
            vpolnor=real(vpole)
         else
            vpolsud=real(vpole)
         endif
      enddo

      if (ordint.ne.voisin.and.ordint.ne.lineair.and.ordint.ne.cubique) then
         print *, 'Tout un probleme mon gars'
      endif
      
 999  return
      end
