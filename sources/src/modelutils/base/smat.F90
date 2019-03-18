!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------
!**s/r smat -  permits to pass from one panel to another (Yin-Yang)
!

      Subroutine smat(s,t,p,x,y)
      implicit none

!author
!     A.Qaddouri  October 2009
!
!revision
!
!object
!      This subroutine gives:
!       1) panel polar coordinates in the other panel 
!       2) matrix to tranform vectors (like wind,coriolis...)
!        
!
!arguments
!        t,p: Lon. and Lat.  in other panel  o
!        S  : matrix for polar vectors trnsformation
!        x,y: Lon. and Lat.  in  panel       i 
! NT: atan uses and gives x between -pi and pi ( user must do shift if necessary)     
!
         real*8 x,y,t,p,s(2,2),atan2,asin,cos,sin
         real*8 det,det1

!
      t=atan2(sin(y),(-cos(y)*cos(x)))
      p=asin(cos(y)*sin(x))
      s(1,1)=-sin(x)*sin(t)
      s(1,2)=cos(x)/cos(p)
      s(2,1)=-cos(x)/cos(p)
      s(2,2)=-sin(x)*sin(t)
!      
      return
      end

