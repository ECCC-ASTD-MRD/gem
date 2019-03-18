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
!**s/r ez_vllfxy(i) - computes the grid co-ordinates of a point
!
      subroutine ez_vllfxy(dlat,dlon,x,y,ni,nj,d60,dgrw,pi,pj,nhem)
      implicit none

!arguments                                                              
!   out   - dlat - latitude in degrees (-90 to +90, positive n).        
!         - dlon - longitude in degrees (-180 to +180, positive e).     
!   in    - x    - x-co-ordinate of the point as measured with pole     
!                  as origin                                            
!         - y    - y-co-ordinate of the point as measured with pole     
!                  as origin                                            
!         - d60  - grid length (in metres) of the polar stereographic   
!                  grid at 60 degrees                                   
!         - dgrw - orientation of greenwich meridian with respect to    
!                  the grid (in degrees)                                
!         - nhem - 1 for northern hemisphere                            
!                  2 for southern hemisphere                            
!                                                                       
!notes    - the companion routine xyfll, which computes the grid        
!           co-ordinates given the latitude and longitude, is also      
!           available.                                                  
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
#include "pi.cdk"
#include "qqqpar.cdk"
      
      integer ni,nj, nhem
      real x(ni,nj), y(ni,nj), dlat(ni,nj), dlon(ni,nj)
      real d60, dgrw, pi, pj
      real*8 x1, y1
      real*8 re,re2,r2,rlat, rlon
      integer i,j
      
      re=1.866025d0*6.371e+6/d60
      re2=re**2
      
      do j=1,nj
         do i=1,ni
            x1 = x(i,j) - pi
            y1 = y(i,j) - pj
            
            if (x1 .eq. 0.  .and. y1 .eq. 0.)  then
               rlat=90.0d0                                   
               rlon=0.0d0
            endif
!     
!     calculate longitude in map coordinates.                         
!     
            if(x1 .eq. 0.0d0) rlon=sign(90.d0,y1)                      
            if(x1 .ne. 0.0d0) rlon=datan(y1/x1)*dble(rdtodg)
            if(x1.lt. 0.0d0) rlon=rlon+sign(180.d0,y1)  
!     
!     * adjust longitude for grid orientation.                          
!     
            rlon=rlon-dgrw                               
!     if(rlon.gt. 180.) rlon=rlon-360.
!     if(rlon.lt.-180.) rlon=rlon+360. 
            if(rlon.lt.0.0d0) rlon=rlon+3.6d2 
!     
!     * calculate latitude.                                             
!     
            r2=x1*x1+y1*y1
            rlat=(re2-r2)/(re2+r2)
            rlat=max(-1.0d0, min(rlat, 1.0d0))
            rlat= dasin(rlat)*dble(rdtodg)
!     change signs if in southern hemisphere.                         
!     
            if(nhem .eq. 2) then
               rlat=-rlat
               rlon=-rlon   
               if(rlon.lt.0.0d0) rlon=rlon+360.0d0
            endif
            dlat(i,j)=rlat
            dlon(i,j)=rlon
         enddo
      enddo
      
!-----------------------------------------------------------------------
!     
!     
      return
      end
!-----------------------------------------------------------------
      
