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

logical function samegrid_file ( unf, p1,p2,p3, g1,g2,g3,g4,xp,yp,ni,nj )
   implicit none
!!!#include <arch_specific.hf>

   integer,intent(in) :: unf            ! Source file unit number
   integer,intent(in) :: p1,p2,p3       ! Source grid search parameters
   integer,intent(in) :: g1,g2,g3,g4    ! Destination rotation
   integer,intent(in) :: ni,nj          ! Destination grid dimensions
   real   ,intent(in) :: xp(ni), yp(nj) ! Destination positions

   !@author   M. Desgagne  -- Winter 2014 --
   !@objective Compare positional parameters

   integer,external :: ezqkdef, fstinf, samegrid_gid
   integer :: src_gid, key, ni1,nj1,nk1,nis,njs

   ! ---------------------------------------------------------------------

   samegrid_file = .false.
   src_gid       = -1

   if ( (fstinf (unf, nis, nj1, nk1, -1, '', p1,p2,p3, '','>>').ge.0) .and.&
        (fstinf (unf, ni1, njs, nk1, -1, '', p1,p2,p3, '','^^').ge.0) ) &
      src_gid = ezqkdef(nis,njs,'Z', p1,p2,p3,-1, unf)

   if (src_gid < 0) then
      key = fstinf (unf, ni1, nj1, nk1, -1, '', p1,p2,p3, '','^>')
      if (key.ge.0) src_gid = ezqkdef(-1,-1,'U', p1,p2,p3,-1, unf)
   endif

   if (key < 0) return

   samegrid_file= samegrid_gid (src_gid, g1,g2,g3,g4,xp,yp,ni,nj) .ne. -1

   ! ---------------------------------------------------------------------
   return
end function samegrid_file
