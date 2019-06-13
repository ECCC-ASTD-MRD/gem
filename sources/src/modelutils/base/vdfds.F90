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

!**s/r vdfds - Computes the first difference g of a function f
!              at unevenly spaced points (vectorized version)
!
      subroutine vdfds (F_g, F_f, F_hr, np, n, F_alfa, F_beta)
      implicit none
!!!#include <arch_specific.hf>
!
      integer n,np
      real    F_g(np,n), F_f(np,n), F_hr(n), F_alfa, F_beta
!
!AUTHOR
!
!REVISION
! v2_30   L. Corbeil             - Version parallele
!
!OBJECT
!        - GIVEN A FUNCTION F AT N UNEVENLY SPACED POINTS, THIS ROUTINE
!        - CALCULATES ITS FIRST DIFFERENCE G AT THESE POINTS.
!        - HR MUST CONTAIN THE INVERSE OF THE INTERVAL LENGTHS.
!        - BOUNDARY CONDITIONS SPECIFIED BY ALFA,BETA, (SEE BELOW). 
!
!ARGUMENTS
!  Name        I/O                 Description
!----------------------------------------------------------------
!   F_g         O         Result
!   F_f         I         Function to be differenced
!   F_hr        I         Inverse of the interval lengths 
!   F_alfa      I         Used for boundary conditions
!   F_beta      I         Used for boundary conditions
!----------------------------------------------------------------------
!
!*
      integer i, pt
      real    a(np)
!
!     ---------------------------------------------------------------
!
      do i=1,n-1
      do pt=1,np
         F_g(pt,i+1) = F_hr(i)*(F_f(pt,i+1)-F_f(pt,i))
      enddo
      enddo
!
      do pt=1,np
        a(pt) = F_g(pt,2)
      enddo
!
      do i=2,n-1
      do pt=1,np
         F_g(pt,i) = (F_hr(i)*F_g(pt,i+1)+F_hr(i-1)*F_g(pt,i)) &
                     /(F_hr(i)+F_hr(i-1))
      enddo
      enddo
!
!     BOUNDARIES
      do pt=1,np
         F_g(pt,1) = F_alfa*a(pt)     +(1.-F_alfa)*F_g(pt,2)
         F_g(pt,n) = F_beta*F_g(pt,n) +(1.-F_beta)*F_g(pt,n-1)
      enddo
!
!     ---------------------------------------------------------------
!
      return
      end 
