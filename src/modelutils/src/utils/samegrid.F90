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

!/@*
logical function samegrid(unf, ni,nj, p1,p2,p3,g1,g2,g3,g4,xp,yp)
   implicit none
!!!#include <arch_specific.hf>
   integer,intent(in) :: unf, ni,nj, p1,p2,p3, g1,g2,g3,g4
   real,intent(in) :: xp(*), yp(*)
   !@author M.Desgagne
   !@objective Compare positional parameters
   !*@/
   logical,external :: is_samegrid_sid
   integer,external :: ezqkdef
   integer :: src_gid
   ! ---------------------------------------------------------------------
   samegrid = .false.
   src_gid = ezqkdef(ni,nj,'Z', p1,p2,p3,-1, unf)
   if (src_gid < 0) src_gid = ezqkdef(-1,-1,'U', p1,p2,p3,-1, unf)
   !TODO: try with # and Y grids?
   if (src_gid < 0) then
      write(6,'(a)') 'ERROR: (samegrid) Can t find positional record describing grid -- ABORT --'
      return
   endif
   samegrid = is_samegrid_sid(src_gid, ni,nj, g1,g2,g3,g4, xp,yp)
   ! ---------------------------------------------------------------------
   return
end function samegrid


!/@*
function is_samegrid2(nis,njs, g1s, g2s, g3s, g4s, xps,yps, &
     nid,njd, g1d, g2d, g3d, g4d, xpd,ypd ) result(F_samegrid_L)
   implicit none
!!!#include <arch_specific.hf>
   !@arguments
   integer,intent(in) :: nis,njs, g1s, g2s, g3s, g4s
   integer,intent(in) :: nid,njd, g1d, g2d, g3d, g4d
   real   ,intent(in) :: xps(nis), yps(njs), xpd(nid), ypd(njd)
   !@return
   logical :: F_samegrid_L
   !@author   M. Desgagne	-- Summer 2010 --
   !@objective Check if grid points are colocated
   !*@/
   integer, external :: tst_parpo
   real,parameter :: eps1 = 1.e-4
   real,parameter :: eps2 = 1.e-5
   integer :: i,j,ideb,jdeb,cnt
   real :: moy
   ! ---------------------------------------------------------------------
   F_samegrid_L = .false.

   if ((g1s.ne.g1d).or.(g2s.ne.g2d).or.(g3s.ne.g3d).or.(g4s.ne.g4d)) return

   do i=1,nis
      ideb = i
      if ( xps(i) .gt. (xpd(1) + eps1) ) exit
   end do
!
!if ideb is last source point, check if xps(nis)>xpd(1)+eps
!before giving correct ideb to tst_parpo
!
!
   if (ideb .lt. nis) then
       ideb = max(1,ideb-1)
   else if ( xps(nis) .le. (xpd(1) + eps1) ) then
       ideb = nis
   else
       ideb = max(1,ideb-1)
   endif

   do j=1,njs
      jdeb = j
      if ( yps(j) .gt. (ypd(1) + eps1) ) exit
   end do
!
!if jdeb is last source point, check if yps(njs)>ypd(1)+eps
!before giving correct jdeb to tst_parpo
!
   if (jdeb .lt. njs) then
       jdeb = max(1,jdeb-1)
   else if ( yps(njs) .le. (ypd(1) + eps1) ) then
       jdeb = njs
   else
       jdeb = max(1,jdeb-1)
   endif

   cnt = tst_parpo (xps,nis, xpd,nid, ideb, eps2)
   moy = real(cnt)/real(nid)
   if (moy.gt.0.2) return

   cnt = tst_parpo (xps,nis, xpd,nid, ideb, eps1)
   if (cnt.gt.0) return

   cnt = tst_parpo (yps,njs, ypd,njd, jdeb, eps2)
   moy = real(cnt)/real(njd)
   if (moy.gt.0.2) return

   cnt = tst_parpo (yps,njs, ypd,njd, jdeb, eps1)
   if (cnt.gt.0) return

   F_samegrid_L = .true.
   ! ---------------------------------------------------------------------
   return    
end function is_samegrid2


!/@*
function samesubgrid(F_sgid, F_nid,F_njd, F_g1d, F_g2d, F_g3d, F_g4d, F_xpd,F_ypd) result(F_subgridid)
   implicit none
!!!#include <arch_specific.hf>
   !@arguments
   integer,intent(in) :: F_sgid !# Source grid id (ezscint)
   integer,intent(in) :: F_nid,F_njd, F_g1d, F_g2d, F_g3d, F_g4d
   real   ,intent(in) :: F_xpd(F_nid), F_ypd(F_njd)
   !@return
   integer :: F_subgridid
   !@objective Check if source grid had grid points colocated to a subgrid
   !*@/
   logical,external :: is_samegrid2
   integer,external :: ezget_nsubgrids,ezget_subgridids,ezgxprm,gdgaxes
   integer :: istat,nis,njs,g1,g2,g3,g4,g1s,g2s,g3s,g4s
   integer :: nsubgrids, igrid
   integer,allocatable :: subgridsid(:)
   character(len=1) :: grd_S, gref_S
   real,allocatable :: xps(:),yps(:)
   ! ---------------------------------------------------------------------
   F_subgridid = -1
   nsubgrids = ezget_nsubgrids(F_sgid)
   allocate(subgridsid(nsubgrids))
   istat = ezget_subgridids(F_sgid,subgridsid)
   if (istat < 0) return
   DOIGRID: do igrid=1,nsubgrids
      !# istat = ezgxprm(subgridsid(igrid),nis,njs,grd_S,g1s,g2s,g3s,g4s,gref_S,g1ref,g2ref,g3ref,g4ref)
      istat = ezgxprm(subgridsid(igrid),nis,njs,grd_S,g1,g2,g3,g4,gref_S,g1s,g2s,g3s,g4s)
      allocate(xps(nis),yps(njs),stat=istat)
      if (istat /= 0) cycle
      istat = gdgaxes(subgridsid(igrid),xps,yps)
      if (is_samegrid2(nis,njs, g1s, g2s, g3s, g4s, xps,yps, &
           F_nid,F_njd, F_g1d, F_g2d, F_g3d, F_g4d, F_xpd,F_ypd)) then
         F_subgridid = subgridsid(igrid)
      endif
      deallocate(xps,yps,stat=istat)
      if (F_subgridid.ge.0) then
          deallocate(subgridsid)
          exit DOIGRID
      endif
   enddo DOIGRID
   ! ---------------------------------------------------------------------
   return
end function samesubgrid


!/@*
function is_samegrid_sid(F_sgid, F_nid,F_njd, F_g1d, F_g2d, F_g3d, F_g4d, F_xpd,F_ypd) result(F_samegrid_L)
   implicit none
!!!#include <arch_specific.hf>
   !@arguments
   integer,intent(in) :: F_sgid !# Source grid id (ezscint)
   integer,intent(in) :: F_nid,F_njd, F_g1d, F_g2d, F_g3d, F_g4d
   real   ,intent(in) :: F_xpd(F_nid), F_ypd(F_njd)
   !@return
   logical :: F_samegrid_L
   !@objective Check if grid points are colocated
   !*@/
   integer,external :: samesubgrid
   integer :: igrid
   ! ---------------------------------------------------------------------
   igrid = samesubgrid(F_sgid, F_nid,F_njd, F_g1d, F_g2d, F_g3d, F_g4d, F_xpd,F_ypd)
   F_samegrid_L = (igrid.ne.-1)
   ! ---------------------------------------------------------------------
   return
end function is_samegrid_sid


!/@*
function tst_parpo(xps,nis, xpd,nid, ideb, epsi) result (F_cnt)
   implicit none
!!!#include <arch_specific.hf>
   !@arguments
   integer,intent(in) :: nis,nid,ideb
   real   ,intent(in) :: xps(nis), xpd(nid), epsi
   !@return
   integer :: F_cnt
   !@author   M. Desgagne	-- Spring 2012 --
   !@objective Check if grid points are colocated
   !*@/
   integer :: i,ii
   real :: r1,r2
   ! ---------------------------------------------------------------------
   F_cnt = 0
   do i=1,nid
      ii = max(min((ideb+i-1),nis),1)
      r1 = max(abs(xps(ii)),epsi)
      r2 = max(abs(xpd( i)),epsi)
      if ( abs((r1-r2)/r1) .gt. epsi) F_cnt = F_cnt + 1
   end do
   ! ---------------------------------------------------------------------
   return    
end function tst_parpo


!/@*
subroutine samegrid2(nis,njs, g1s, g2s, g3s, g4s, xps,yps, &
     nid,njd, g1d, g2d, g3d, g4d, xpd,ypd, inttyp)
   implicit none
!!!#include <arch_specific.hf>

   !@arguments
   character(len=*),intent(inout) ::  inttyp
   integer,intent(in) :: nis,njs, g1s, g2s, g3s, g4s
   integer,intent(in) :: nid,njd, g1d, g2d, g3d, g4d
   real   ,intent(in) :: xps(nis), yps(njs), xpd(nid), ypd(njd)

   !@author   M. Desgagne	-- Summer 2010 --
   !@objective Set interp_type=nearest if grid points are colocated
   !*@/
   logical,external :: is_samegrid2
   ! ---------------------------------------------------------------------
   if (inttyp == 'NEAREST' .or. inttyp == 'nearest') return
   if (is_samegrid2(nis,njs, g1s, g2s, g3s, g4s, xps,yps, &
        nid,njd, g1d, g2d, g3d, g4d, xpd,ypd)) &
        inttyp = 'NEAREST'
   ! ---------------------------------------------------------------------
   return
end subroutine samegrid2


!/@*
subroutine samegrid_sid(F_sgid, F_nid,F_njd, F_g1d, F_g2d, F_g3d, F_g4d, F_xpd,F_ypd, F_inttyp_S)
   implicit none
!!!#include <arch_specific.hf>
   !@arguments
   character(len=*),intent(inout) :: F_inttyp_S
   integer,intent(in) :: F_sgid !# Source grid id (ezscint)
   integer,intent(in) :: F_nid,F_njd, F_g1d, F_g2d, F_g3d, F_g4d
   real   ,intent(in) :: F_xpd(F_nid), F_ypd(F_njd)
   !@objective Set interp_type=nearest if grid points are colocated
   !*@/
   logical,external :: is_samegrid_sid
   ! ---------------------------------------------------------------------
   if (F_inttyp_S == 'NEAREST' .or. F_inttyp_S == 'nearest') return
   if (is_samegrid_sid(F_sgid,F_nid,F_njd,F_g1d,F_g2d,F_g3d,F_g4d,F_xpd,F_ypd)) &
        F_inttyp_S = 'NEAREST'
   ! ---------------------------------------------------------------------
   return
end subroutine samegrid_sid
