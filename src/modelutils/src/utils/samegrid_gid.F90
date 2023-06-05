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

integer function samegrid_gid ( F_sgid, F_g1d,F_g2d,F_g3d,F_g4d, &
                                F_xpd,F_ypd,F_nid,F_njd )
   implicit none
!!!#include <arch_specific.hf>

   !@arguments
   integer,intent(in) :: F_sgid                     ! Source grid id (ezscint)
   integer,intent(in) :: F_g1d, F_g2d, F_g3d, F_g4d ! Destination rotation
   integer,intent(in) :: F_nid,F_njd                ! Destination grid dimensions
   real   ,intent(in) :: F_xpd(F_nid), F_ypd(F_njd) ! Destination positions

   !@author   M. Desgagne  -- Winter 2014 --
   !@objective Check if source grid had grid points colocated to a subgrid

   integer,external :: ezget_nsubgrids, ezget_subgridids, ezgxprm, &
                       gdgaxes, samegrid_parpos

   integer :: istat,nis,njs,g1,g2,g3,g4,g1s,g2s,g3s,g4s
   integer :: nsubgrids, igrid, i, j, ideb, jdeb, cnt
   integer, allocatable :: subgridsid(:)
   character(len=1)  :: grd_S, gref_S
   real, allocatable :: xps(:), yps(:)
   real, parameter   :: eps1 = 1.e-4
   real, parameter   :: eps2 = 1.e-5
   real :: moy

   ! ---------------------------------------------------------------------

   samegrid_gid = -1

   nsubgrids = ezget_nsubgrids(F_sgid)

   allocate(subgridsid(nsubgrids))
   istat = ezget_subgridids(F_sgid,subgridsid)

   if (istat < 0) return

   DOIGRID: do igrid= 1, nsubgrids

      istat = ezgxprm(subgridsid(igrid),nis,njs,grd_S,g1,g2,g3,g4,gref_S,g1s,g2s,g3s,g4s)

      if ((g1s.ne.F_g1d).or.(g2s.ne.F_g2d)  .or. &
          (g3s.ne.F_g3d).or.(g4s.ne.F_g4d)) cycle

      allocate (xps(nis),yps(njs))
      istat = gdgaxes(subgridsid(igrid),xps,yps)

      do i=1,nis
         ideb = i
         if ( xps(i) .gt. (F_xpd(1) + eps1) ) exit
      end do
!
!if ideb is last source point, check if xps(nis)>F_xpd(1)+eps
!before giving correct ideb to tst_parpo
!
!
   if (ideb .lt. nis) then
       ideb = max(1,ideb-1)
   else if ( xps(nis) .le. (F_xpd(1) + eps1) ) then
       ideb = nis
   else
       ideb = max(1,ideb-1)
   endif

      do j=1,njs
         jdeb = j
         if ( yps(j) .gt. (F_ypd(1) + eps1) ) exit
      end do
!
!if jdeb is last source point, check if yps(njs)>F_ypd(1)+eps
!before giving correct jdeb to tst_parpo
!
   if (jdeb .lt. njs) then
       jdeb = max(1,jdeb-1)
   else if ( yps(njs) .le. (F_ypd(1) + eps1) ) then
       jdeb = njs
   else
       jdeb = max(1,jdeb-1)
   endif


      cnt = samegrid_parpos (xps,nis, F_xpd,F_nid, ideb, eps2)
      moy = real(cnt)/real(F_nid)
      if (moy.gt.0.2) goto 667

      cnt = samegrid_parpos (xps,nis, F_xpd,F_nid, ideb, eps1)
      if (cnt.gt.0) goto 667

      cnt = samegrid_parpos (yps,njs, F_ypd,F_njd, jdeb, eps2)
      moy = real(cnt)/real(F_njd)
      if (moy.gt.0.2) goto 667

      cnt = samegrid_parpos (yps,njs, F_ypd,F_njd, jdeb, eps1)
      if (cnt.gt.0) goto 667

      samegrid_gid= subgridsid(igrid)

 667  deallocate(xps,yps,stat=istat)

      if (samegrid_gid.ge.0) then
         deallocate(subgridsid)
         exit DOIGRID
      endif

   enddo DOIGRID

   ! ---------------------------------------------------------------------
   return
end function samegrid_gid

integer function samegrid_rot ( F_sgid, F_g1d,F_g2d,F_g3d,F_g4d )
   implicit none
!!!#include <arch_specific.hf>

   !@arguments
   integer,intent(in) :: F_sgid                     ! Source grid id (ezscint)
   integer,intent(in) :: F_g1d, F_g2d, F_g3d, F_g4d ! Destination rotation

   !@author   M. Desgagne  -- Fall 2015 --
   !@objective Check if source grid has same rotation as target grid

   integer,external :: ezget_nsubgrids, ezget_subgridids, ezgxprm

   integer :: istat,g1,g2,g3,g4,g1s,g2s,g3s,g4s,nis,njs
   integer :: nsubgrids, igrid
   integer, allocatable :: subgridsid(:)
   character(len=1)  :: grd_S, gref_S

   ! ---------------------------------------------------------------------

   samegrid_rot= -1

   nsubgrids= max(1,ezget_nsubgrids(F_sgid))

   allocate ( subgridsid(nsubgrids) )
   istat= ezget_subgridids(F_sgid,subgridsid)

   if (istat < 0) goto 99

   DOIGRID: do igrid= 1, nsubgrids

      istat= ezgxprm( subgridsid(igrid),nis,njs,grd_S   ,&
                      g1,g2,g3,g4,gref_S,g1s,g2s,g3s,g4s )

      if ((g1s.ne.F_g1d).or.(g2s.ne.F_g2d)  .or. &
          (g3s.ne.F_g3d).or.(g4s.ne.F_g4d)) cycle

      samegrid_rot= 1

   enddo DOIGRID

99 deallocate (subgridsid)

   ! ---------------------------------------------------------------------
   return
end function samegrid_rot
