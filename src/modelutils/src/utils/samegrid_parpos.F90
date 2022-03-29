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

integer function samegrid_parpos (xps,nis, xpd,nid, ideb, epsi)
   implicit none
!!!#include <arch_specific.hf>

   !@arguments
   integer,intent(in) :: nis,nid,ideb
   real   ,intent(in) :: xps(nis), xpd(nid), epsi

   !@author   M. Desgagne  -- Winter 2014 --
   !@objective Check if grid points are colocated

   integer :: i,ii
   real :: r1,r2

   ! ---------------------------------------------------------------------
   samegrid_parpos = 0
   do i=1,nid
      ii = max(min((ideb+i-1),nis),1)
      r1 = max(abs(xps(ii)),epsi)
      r2 = max(abs(xpd( i)),epsi)
      if ( abs((r1-r2)/r1) .gt. epsi) &
           samegrid_parpos = samegrid_parpos + 1
   end do

   ! ---------------------------------------------------------------------
   return    
end function samegrid_parpos
