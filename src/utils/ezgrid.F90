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

!/@
module ezgrid_mod
   use iso_c_binding
   use rpn_comm_itf_mod
   use ptopo_utils
   use agg_filter_mod
   implicit none
   private
   !@objective Manipulates ezscint grids
   !@author  Stephane Chamberland, 2011-03
   !@revisions  
   !  2012-03, S.Chamberland: New functions for output_mod and offline
   !@description
   ! Public functions
   public :: ezgrid_bcast,ezgrid_serialize,ezgrid_unserialize,ezgrid_find_ij0,&
        ezgrid_params, ezgrid_sameproj, ezgrid_samegrid_params, ezgrid_samegrid, ezgrid_colocated, ezgrid_subcolocated, &
        ezgrid_sub, ezgrid_merge,ezgrid_latlon,ezgrid_addperiod
   !TODO-later: other grid fn: grid_split, grid_intersect, grid_comp...
   character(len=1),parameter,public :: EZGRID_REF_TYPES(5) = (/'Z','z','#','Y','y'/)
   integer,parameter,public :: EZGRID_180_180 = 1
   integer,parameter,public :: EZGRID_0_360 = 2
   !@/

#include <rmn/msg.h>
#include <rmnlib_basics.hf>

!!$   integer,external :: gdll,ezget_nsubgrids,ezget_subgridids,ezgdef_supergrid,msg_getUnit
   logical,external :: is_samegrid2,is_samegrid_sid

   !TODO-later: compute CHARPERBYTE value
   integer,parameter :: CHARPERBYTE = 4
   integer,parameter :: GRTYPLEN = 2
   integer,parameter :: NSUBGRIDS_MAX = 64

   type :: myprt2d
      sequence
      real,pointer :: p(:,:)
   end type myprt2d

contains


   !/@
   function ezgrid_params(F_gridid,F_nij,F_grtyp_S,F_grref_S,F_ig14,F_ij0,F_ax,F_ay,F_igp14) result(F_istat)
      implicit none
      !@objective Return grid parameters
      !@arguments
      integer,intent(in) :: F_gridid
      integer,intent(out),optional :: F_nij(2),F_ig14(4),F_ij0(2)
      character(len=*),intent(out),optional :: F_grtyp_S,F_grref_S
      real,pointer,optional :: F_ax(:,:),F_ay(:,:)
      integer,intent(out),optional :: F_igp14(4)  !- grid tag/id for ref grids
      !@author
      !@return
      integer :: F_istat
      !@description
      ! Grids at RPN in RPNStd files are encoded in 2 ways
      ! 1) Direct way in the field record description with 6 parameters
      !    nij    : horizontal dims of the grid (ni,nj)
      !    grtyp_S: Projection type = 'A', 'B', 'E', 'G', 'L', 'N', 'S'
      !    ig14   : encoded projection parameters (ig1,ig2,ig3,ig4)
      ! 2) An indirect way using additional grid records (ax,ay)
      !    that are actual grid position values in the projection
      !    (ax,ay units are projection type dependent)
      !    The full grid descrption then comes from 3 records.
      !    They come in 3 variants:
      ! 2.1) Z grids
      !    - the field data record
      !      nij    : horizontal dims of the grid (ni,nj)
      !      grtyp_S = 'Z'
      !      ig1-3  : 3 numbers for the grid id
      !      ig4    : ignored
      !    - the ax,ay records
      !      ip13   : grid id (as in ig1-3 of the field data record)
      !      nij    : ax record nij =(/ni,1/); ay record nij =(/1,nj/)
      !      grref_S: Projection type = 'A', 'B', 'E', 'G', 'L', 'N', 'S'
      !      ig14ref: encoded projection parameters (ig1,ig2,ig3,ig4)
      !      ax,ay values must be monotonuously increasing
      ! 2.2) # grids (as Z-grids but for a sub-grid)
      !    - the field data record
      !      ip3    : tile number (optional)
      !      nij    : horizontal dims of the fields/sub-grid (ni,nj)
      !      grtyp_S = '#'
      !      ig1-2  : 2 numbers for the grid id
      !      ig3-4  : i0,j0 of the sub-grid in the full grid described by ax,ay
      !    - the ax,ay records
      !      ip1-2  : grid id (as in ig1-2 of the field data record)
      !      nijref : ax (/niref,1/); ay (/1,njref/)
      !               note that we must have nijref >= nij+ij0
      !               it is not part of the interface but can be obtained with:
      !               niref = size(ax,1) ; njref = size(ay,2)
      !      grref_S: Projection type = 'A', 'B', 'E', 'G', 'L', 'N', 'S'
      !      ig14ref: encoded projection parameters (ig1,ig2,ig3,ig4)
      !      ax,ay values must be monotonuously increasing
      ! 2.3) Y grids - clouds of points (not a regular grid)
      !      Y grids are encoded as Z grids with 2 minor differences
      !      grtyp_S = 'Y'
      !      ax,ay records are 2d arrays of size/dims = (/ni,nj/)
      !    For the Z,#,Y grids, F_ig14ref are returned as F_ig14
      ! 2.4) U grids - multiple grids
      !      F means Yin-yang grid where it contains 2 Z grids with grref=E
      !@/
      integer :: istat, ig14ref(4),i,j,n
      integer :: nij(2),ig14(4),ij0(2),igp14(4)
      character(len=GRTYPLEN) :: grtyp_S,grref_S
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR
      if (.not.RMN_IS_OK(F_gridid)) then
         call msg(MSG_ERROR,'(ezgrid_params) invalid grid id')
         return
      endif

      ij0 = (/1,1/)
!!$      print *,'1a ij0=',ij0
      istat = ezgxprm(F_gridid, nij(1), nij(2), &
           grtyp_S(1:1), ig14(1),ig14(2),ig14(3),ig14(4), &
           grref_S(1:1), ig14ref(1),ig14ref(2),ig14ref(3),ig14ref(4))
!!$      print *,'1b ij0=',ij0,grtyp_S(1:1),grref_S(1:1)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_ERROR,'(ezgrid_params) ezgxprm error, connot get grid params')
         return
      endif

      igp14 = 0

      if (any(grtyp_S(1:1) == EZGRID_REF_TYPES)) then
!!$         igp14(1:3) = ig14(1:3)
         ig14 = ig14ref
         if (grtyp_S(1:1) == '#') then
            ij0 = (/ig14(3),ig14(4)/)
            igp14(3:4) = ij0(1:2)
         endif
!!$         if (present(F_igp14)) print *,'1 igp14=',igp14
      endif
!!$      print *,'2 ij0=',ij0

      IF_AXY: if (present(F_ax) .and. present(F_ay)) then

         if (associated(F_ax) .or. associated(F_ay)) return
         select case(grtyp_S(1:1))
         case('Y')
            allocate(F_ax(nij(1),nij(2)), &
                 F_ay(nij(1),nij(2)), &
                 stat=istat)
!!$         case('#') !TODO-later: does gdgaxes return full or trimed grid?
!!$            allocate(F_ax(nij(1),1),F_ay(1,nij(2)),stat=istat)
         case default !('Z' and '#')
            allocate(F_ax(nij(1),1), &
                 F_ay(1,nij(2)), &
                 stat=istat)
         end select
         if (istat /= 0) then
            call msg(MSG_ERROR,'(ezgrid_params) allocate error')
            return
         endif

         if (any(grtyp_S(1:1) == EZGRID_REF_TYPES)) then

            istat = gdgaxes(F_gridid, F_ax, F_ay)
            if (.not.RMN_IS_OK(istat)) then
               call msg(MSG_ERROR,'(ezgrid_params) gdgaxes error, connot get grid axes')
               if (associated(F_ax)) deallocate(F_ax,stat=istat)
               if (associated(F_ay)) deallocate(F_ay,stat=istat)
               return
            endif

            if (present(F_igp14)) then
               call set_igs2(igp14(1),igp14(2), &
                    F_ax, F_ay, nij(1), nij(2),      &
                    ig14(1),ig14(2),ig14(3),ig14(4), &
                    ij0(1), nij(1), 1,  ij0(2), nij(2), 1 )
!!$               print *,'2 igp14=',igp14
            endif

         else

            !TODO: should return values depending on grid type

            do j=1,size(F_ax,2)
               do i=1,size(F_ax,1)
                  F_ax(i,j) = float(i)
               enddo
            enddo
            do j=1,size(F_ay,2)
               do i=1,size(F_ay,1)
                  F_ay(i,j) = float(j)
               enddo
            enddo

         endif

      endif IF_AXY
!!$      print *,'3 ij0=',ij0

      if (present(F_grtyp_S)) F_grtyp_S = grtyp_S
      if (present(F_grref_S)) F_grref_S = grref_S
      if (present(F_nij)) then
         F_nij = 0
         n = min(max(1,size(F_nij)),size(nij))
         F_nij(1:n) = nij(1:n)
      endif
      if (present(F_ig14)) then
         F_ig14 = 0
         n = min(max(1,size(F_ig14)),size(ig14))
         F_ig14(1:n) = ig14(1:n)
      endif
      if (present(F_ij0)) then
         F_ij0  = 0
         n = min(max(1,size(F_ij0)),size(ij0))
         F_ij0(1:n) = ij0(1:n)
      endif
      if (present(F_igp14)) then
         F_igp14 = 0
         n = min(max(1,size(F_igp14)),size(igp14))
         F_igp14(1:n) = igp14(1:n)
!!$         print *,n,'F_igp14=',F_igp14(1:n)
      endif

      F_istat = RMN_OK
      ! ---------------------------------------------------------------------
      return
   end function ezgrid_params


   !/@
   function ezgrid_find_ij0(F_gridid_local,F_gridid_full,F_i0,F_j0,F_lni,F_lnj) result(F_istat)
      implicit none
      !@objective Find position F_gridid_local(1,1) in F_gridid_full
      !@arguments
      integer,intent(in) :: F_gridid_local,F_gridid_full
      integer,intent(out) :: F_i0,F_j0
      integer,intent(out),optional :: F_lni,F_lnj
      !@author
      !@return
      integer :: F_istat
      !@/
      logical :: coloc_L
      integer :: istat,lnij(2),gnij(2),minlocij(2)
      real,pointer :: lax(:,:),lay(:,:),gax(:,:),gay(:,:)
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR
      F_i0=0; F_j0=0; 
      if (present(F_lni)) F_lni = 0
      if (present(F_lnj)) F_lnj = 0
      coloc_L = ezgrid_colocated(F_gridid_local,F_gridid_full)
      if (.not.coloc_L) &
           coloc_L = ezgrid_colocated(F_gridid_full,F_gridid_local)
      if (.not.coloc_L) then
         call msg(MSG_ERROR,'(ezgrid) find_ij0 - Grids must be collocated')
         return
      endif
      nullify(lax,lay,gax,gay)
      F_istat = ezgrid_params(F_gridid_local,lnij,F_ax=lax,F_ay=lay)
      F_istat = min(ezgrid_params(F_gridid_full,gnij,F_ax=gax,F_ay=gay),F_istat)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_ERROR,'(ezgrid) find_ij0 - Probleme getting grid params')
         return
      endif
      if (present(F_lni)) F_lni = lnij(1)
      if (present(F_lnj)) F_lnj = lnij(2)
      gax = abs(gax - lax(1,1))
      gay = abs(gay - lay(1,1))
      minlocij = minloc(gax); F_i0 = minlocij(1)
      minlocij = minloc(gay); F_j0 = minlocij(2) 
      deallocate(lax,lay,gax,gay,stat=istat)
      F_istat = F_gridid_local
      ! ---------------------------------------------------------------------
      return
   end function ezgrid_find_ij0


   !/@
   function ezgrid_latlon(F_gridid,F_lat,F_lon,F_torad_L,F_lon_mode) result(F_istat)
      implicit none
      !@objective Compute Grid Lat Lon with simple transform
      !@arguments
      integer,intent(in) :: F_gridid
      real,pointer :: F_lat(:,:),F_lon(:,:)
      logical,intent(in),optional :: F_torad_L
      integer,intent(in),optional :: F_lon_mode
      !@author
      !@return
      integer :: F_istat
      !@/
      real :: deg2rad
      integer :: istat,nij(2),ig14(4),ij0(2)
      character(len=GRTYPLEN) :: grtyp_S,grref_S
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR
      istat  = ezgrid_params(F_gridid,nij,grtyp_S,grref_S,ig14,ij0)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_ERROR,'(ezgrid_latlon) Problem getting grid params')
         return
      endif
      if (.not.associated(F_lat)) allocate(F_lat(nij(1),nij(2)),stat=istat)
      if (.not.associated(F_lon)) allocate(F_lon(nij(1),nij(2)),stat=istat)
      if (any(shape(F_lat)/=nij) .or. any(shape(F_lat)/=nij)) then
         call msg(MSG_ERROR,'(ezgrid_latlon) Wrong array size')
         return
      endif
      F_istat = gdll(F_gridid,F_lat,F_lon)
      if (present(F_lon_mode)) then
         if (F_lon_mode == EZGRID_0_360) then
            where(F_lon < 0.)
               F_lon = F_lon + 360.
            endwhere
         elseif (F_lon_mode == EZGRID_180_180) then
            where(F_lon > 180.) F_lon = F_lon - 180.
         endif
      endif
      if (present(F_torad_L)) then
         if (F_torad_L) then
            deg2rad = acos(-1.)/180.
            F_lat = F_lat * deg2rad
            F_lon = F_lon * deg2rad
         endif
      endif
      ! ---------------------------------------------------------------------
      return
   end function ezgrid_latlon


   !/@
   function ezgrid_sameproj(F_gridid,F_gridid2) result(F_sameproj_L)
      implicit none
      !@objective Compare grids projection
      !@arguments
      integer,intent(in) :: F_gridid,F_gridid2
      !@author
      !@return
      logical :: F_sameproj_L
      !@/
      integer :: istat,istat2,nij(2),ig14(4),nij2(2),ig14_2(4),ij0(2),ij0_2(2)
      character(len=GRTYPLEN) :: grtyp_S,grref_S,grtyp2_S,grref2_S
      ! ---------------------------------------------------------------------
      if (F_gridid == F_gridid2) then
         F_sameproj_L = .true.
         return
      endif

      F_sameproj_L = .false.
      istat  = ezgrid_params(F_gridid,nij,grtyp_S,grref_S,ig14,ij0)
      istat2 = ezgrid_params(F_gridid2,nij2,grtyp2_S,grref2_S,ig14_2,ij0_2)
      if (.not.RMN_IS_OK(min(istat,istat2))) then
         call msg(MSG_ERROR,'(ezgrid_sameproj) problem getting grid params, cannot compare')
         return
      endif
      !TODO-later: for L,G,?  grids ig14 are not defining the proj
      F_sameproj_L = (all(ig14==ig14_2) .and. grtyp_S(1:1)==grtyp2_S(1:1) .and. grref_S(1:1)==grref2_S(1:1))
      ! ---------------------------------------------------------------------
      return
   end function ezgrid_sameproj


   !/@
   function ezgrid_samegrid(F_gridid,F_gridid2) result(F_samegrid_L)
      implicit none
      !@objective Compare grids projection, localisation and extent (grids must be identical)
      !@arguments
      integer,intent(in) :: F_gridid,F_gridid2
      !@author
      !@return
      logical :: F_samegrid_L
      !@/
      real,parameter :: EPSILON = 0.000001 !TODO: ok for LAT/LON but else?
      integer :: istat,istat2,msgunit,nij(2),ig14(4),nij2(2),ig14_2(4),ij0(2),ij0_2(2)
      character(len=GRTYPLEN) :: grtyp_S,grref_S,grtyp2_S,grref2_S
      real,pointer :: ax(:,:),ay(:,:), ax2(:,:),ay2(:,:)
      ! ---------------------------------------------------------------------
      if (F_gridid == F_gridid2) then
         F_samegrid_L = .true.
         return
      endif

      F_samegrid_L = .false.
      nullify(ax,ay, ax2,ay2)
      istat  = ezgrid_params(F_gridid,nij,grtyp_S,grref_S,ig14,ij0,ax,ay)
      istat2 = ezgrid_params(F_gridid2,nij2,grtyp2_S,grref2_S,ig14_2,ij0_2,ax2,ay2)
      if (.not.RMN_IS_OK(min(istat,istat2))) then
         call msg(MSG_ERROR,'(ezgrid_samegrid) problem getting grid params, cannot compare')         
         return
      endif

      if (all(ig14==ig14_2) .and. all(nij==nij2) .and. all(ij0==ij0_2) .and. &
           grtyp_S(1:1)==grtyp2_S(1:1) .and. grref_S(1:1)==grref2_S(1:1)) then
         if (any(grtyp_S(1:1) == EZGRID_REF_TYPES)) then
            if (associated(ax).and.associated(ax2).and.associated(ay).and.associated(ay2)) then
               !TODO: epsilon
               if (all(abs(ax-ax2)<EPSILON) .and. all(abs(ay-ay2)<EPSILON)) F_samegrid_L = .true.
            else
               call msg(MSG_INFO,'(ezgrid_samegrid) Cannot compare, missing src/dst axes')
            endif
         else
            F_samegrid_L = .true.
         endif
      endif
      msgunit = msg_getUnit(MSG_DEBUG)
      if (.not.F_samegrid_L .and. RMN_IS_OK(msgunit)) then
         print '(a,a,8i6)','ezgrid_samegrid 1: ',grtyp_S(1:1)//grref_S(1:1),nij,ig14,ij0
         print '(a,a,8i6)','ezgrid_samegrid 2: ',grtyp2_S(1:1)//grref2_S(1:1),nij2,ig14_2,ij0_2
      endif

      if (associated(ax)) deallocate(ax,stat=istat)
      if (associated(ay)) deallocate(ay,stat=istat)
      if (associated(ax2)) deallocate(ax2,stat=istat)
      if (associated(ay2)) deallocate(ay2,stat=istat)
      ! ---------------------------------------------------------------------
      return
   end function ezgrid_samegrid


   !/@
   function ezgrid_samegrid_params(F_gridid,nij2,grtyp2_S,grref2_S,ig14_2,ij0_2,ax2,ay2) result(F_samegrid_L)
      implicit none
      !@objective Compare grids projection, localisation and extent (grids must be identical)
      !@arguments
      integer,intent(in) :: F_gridid
      integer :: nij2(2),ig14_2(4),ij0_2(2)
      character(len=GRTYPLEN),intent(in) :: grtyp2_S,grref2_S
      real,pointer :: ax2(:,:),ay2(:,:)
      !@author
      !@return
      logical :: F_samegrid_L
      !@/
      real,parameter :: EPSILON = 0.000001 !TODO: ok for LAT/LON but else?
      integer :: istat,istat2,msgunit,nij(2),ig14(4),ij0(2)
      character(len=GRTYPLEN) :: grtyp_S,grref_S
      real,pointer :: ax(:,:),ay(:,:)
      ! ---------------------------------------------------------------------
      F_samegrid_L = .false.
      nullify(ax,ay)
      istat  = ezgrid_params(F_gridid,nij,grtyp_S,grref_S,ig14,ij0,ax,ay)
      if (.not.RMN_IS_OK(min(istat,istat2))) then
         call msg(MSG_ERROR,'(ezgrid_samegrid) problem getting grid params, cannot compare')         
         return
      endif

      if (all(ig14==ig14_2) .and. all(nij==nij2) .and. all(ij0==ij0_2) .and. &
           grtyp_S(1:1)==grtyp2_S(1:1) .and. grref_S(1:1)==grref2_S(1:1)) then
         if (any(grtyp_S(1:1) == EZGRID_REF_TYPES)) then
            if (associated(ax).and.associated(ax2).and.associated(ay).and.associated(ay2)) then
               if (all(abs(ax-ax2)<EPSILON) .and. all(abs(ay-ay2)<EPSILON)) F_samegrid_L = .true.
            endif
         else
            F_samegrid_L = .true.
         endif
      endif
      msgunit = msg_getUnit(MSG_DEBUG)
      if (.not.F_samegrid_L .and. RMN_IS_OK(msgunit)) then
         print '(a,a,8i6)','ezgrid_samegrid 1: ',grtyp_S(1:1)//grref_S(1:1),nij,ig14,ij0
         print '(a,a,8i6)','ezgrid_samegrid 2: ',grtyp2_S(1:1)//grref2_S(1:1),nij2,ig14_2,ij0_2
      endif

      if (associated(ax)) deallocate(ax,stat=istat)
      if (associated(ay)) deallocate(ay,stat=istat)
      ! ---------------------------------------------------------------------
      return
   end function ezgrid_samegrid_params


   !/@
   function ezgrid_colocated(F_src_gridid,F_dst_gridid) result(F_coloc_L)
      implicit none
      !@objective Compare grids projection, localisation and extent for point colacation (src_grid must be a sub-set of dst_grid)
      !@arguments
      integer,intent(in) :: F_src_gridid,F_dst_gridid
      !@author
      !@return
      logical :: F_coloc_L
      !@/
      integer :: gridid
      ! ---------------------------------------------------------------------
      if (F_src_gridid == F_dst_gridid) then
         F_coloc_L = .true.
         return
      endif
      gridid = ezgrid_subcolocated(F_src_gridid,F_dst_gridid)
      F_coloc_L = (gridid >= 0)
      ! ---------------------------------------------------------------------
      return
   end function ezgrid_colocated


   !/@
   function ezgrid_subcolocated(F_src_gridid,F_dst_gridid) result(F_subgrid)
      implicit none
      !@objective Compare grids projection, localisation and extent for point colacation (src_grid must be a sub-set of dst_grid)
      !@arguments
      integer,intent(in) :: F_src_gridid,F_dst_gridid
      integer,external :: samesubgrid
      !@author
      !@return
      integer :: F_subgrid
      !@/
      integer :: istat,istat2,nij(2),ig14(4),nij2(2),ig14_2(4),ij0(2),ij0_2(2)
      character(len=GRTYPLEN) :: grtyp_S,grref_S,grtyp2_S,grref2_S
      real,pointer :: ax(:,:),ay(:,:), ax2(:,:),ay2(:,:)
      ! ---------------------------------------------------------------------
      nullify(ax,ay, ax2,ay2)
      istat  = ezgrid_params(F_src_gridid,nij,grtyp_S,grref_S,ig14,ij0,ax,ay)
      istat2 = ezgrid_params(F_dst_gridid,nij2,grtyp2_S,grref2_S,ig14_2,ij0_2,ax2,ay2)
      F_subgrid=-1
      if (.not.RMN_IS_OK(min(istat,istat2))) then
         call msg(MSG_ERROR,'(ezgrid_subcolocated) problem getting grid params, cannot compare')         
         return
      endif
!     print *,'ezgrid_params:grtyp_S',grtyp_S(1:1),'grtyp2_S',grtyp2_S(1:1)
      if (any(grtyp_S(1:1) == (/'U','u'/)) .and. any(grtyp2_S(1:1) == EZGRID_REF_TYPES)) then

         if (associated(ax2).and.associated(ay2)) then
!           print *,'ez_subcoloc looks for subgridid'
            F_subgrid = samesubgrid(F_src_gridid, &
                 nij2(1),nij2(2), ig14_2(1),ig14_2(2),ig14_2(3),ig14_2(4), ax2,ay2)
         endif

      else if (F_src_gridid == F_dst_gridid) then
         F_subgrid = F_src_gridid
         return
      
      else if (any(grtyp_S(1:1) == EZGRID_REF_TYPES) .and. any(grtyp2_S(1:1) == EZGRID_REF_TYPES)) then

         !TODO-later: add ij0 subgrid comparison
         if (associated(ax).and.associated(ax2).and.associated(ay).and.associated(ay2)) then
            if( is_samegrid2(&
                 nij(1) ,nij(2) , ig14(1),ig14(2),ig14(3),ig14(4), ax,ay, &
                 nij2(1),nij2(2), ig14_2(1),ig14_2(2),ig14_2(3),ig14_2(4), ax2,ay2)) F_subgrid=F_src_gridid
         endif

      else
         if (ezgrid_samegrid(F_src_gridid,F_dst_gridid)) F_subgrid=F_src_gridid
      endif

      if (associated(ax)) deallocate(ax,stat=istat)
      if (associated(ay)) deallocate(ay,stat=istat)
      if (associated(ax2)) deallocate(ax2,stat=istat)
      if (associated(ay2)) deallocate(ay2,stat=istat)
      ! ---------------------------------------------------------------------
      return
   end function ezgrid_subcolocated

   !TODO-later: split out into ezgrid_mpi
   !/@
   function ezgrid_bcast(F_gridid,F_comm_S,F_from) result(F_grididout)
      implicit none
      !@objective MPI broadcast gridid
      !@arguments
      integer,intent(in) :: F_gridid
      character(len=*),intent(in) :: F_comm_S
      integer,intent(in),optional :: F_from
      !@author
      !@return
      integer :: F_grididout
      !@/
      logical :: ismaster_L
      integer :: istat,ns(2),me,master,mysize
      integer,pointer :: griddata(:)
      ! ---------------------------------------------------------------------
      call msg_toall(MSG_DEBUG,'(ezgrid_bcast) [BEGIN]')
      F_grididout = RMN_ERR

      master = RPN_COMM_MASTER
      if (present(F_from)) master = F_from
      call rpn_comm_rank(F_comm_S,me,istat)
      ismaster_L = (me == master)
      
      nullify(griddata)
      if (ismaster_L) then
         istat = ezgrid_serialize(F_gridid,griddata)
         ns(1) = size(griddata)
         ns(2) = istat
      endif
      mysize = size(ns)
      call rpn_comm_bcast(ns,mysize,RPN_COMM_INTEGER,master,F_comm_S,istat)
      istat = min(ns(2),istat)
      call collect_error(istat)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_ERROR,'(ezgrid_bcast) serialization error, cannot bcast')
         return
      endif
      istat = 0
      if (.not.ismaster_L) then
         allocate(griddata(ns(1)),stat=istat)
         griddata = 0
      endif
      call collect_error(istat)
      if (istat /= 0) then
         call msg(MSG_ERROR,'(ezgrid_bcast) allocate error')
         return
      endif

      mysize = ns(1)
      call rpn_comm_bcast(griddata,mysize,RPN_COMM_INTEGER,master,F_comm_S,istat)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_ERROR,'(ezgrid_bcast) rpn_comm_bcast error')
         return
      endif
      if (.not.ismaster_L) then
         F_grididout = ezgrid_unserialize(griddata)
      else
         F_grididout = F_gridid
      endif

      deallocate(griddata,stat=istat)
      call msg_toall(MSG_DEBUG,'(ezgrid_bcast) [END]')
      ! ---------------------------------------------------------------------
      return
   end function ezgrid_bcast


   !/@
   function ezgrid_serialize(F_gridid,F_griddata) result(F_istat)
      implicit none
      !@objective Transfer grid definition in a integer array
      !@arguments
      integer,intent(in) :: F_gridid
      integer,pointer :: F_griddata(:)
      !@return
      integer :: F_istat
      !@author 
      !@/
      integer,parameter :: NINTEGER = 14
      integer :: gridobjsize, istat, nij(2), ig14(4),ij0(2),i0,i00,i,j,axnij(2,NSUBGRIDS_MAX),aynij(2,NSUBGRIDS_MAX),nsubgrids,igrid
      type(myprt2d) :: axp(NSUBGRIDS_MAX),ayp(NSUBGRIDS_MAX)
      character(len=GRTYPLEN) :: grtyp_S,grref_S
      integer :: subgridsid(NSUBGRIDS_MAX)
      !---------------------------------------------------------------------
      F_istat = RMN_ERR
      if (F_gridid < 0) return

      nsubgrids = min(ezget_nsubgrids(F_gridid),NSUBGRIDS_MAX)
      istat = ezget_subgridids(F_gridid,subgridsid(1:nsubgrids))
      do igrid=1,nsubgrids
         nullify(axp(igrid)%p,ayp(igrid)%p)
         istat  = ezgrid_params(subgridsid(igrid),nij,grtyp_S,grref_S,ig14,ij0,axp(igrid)%p,ayp(igrid)%p)
         if (.not.RMN_IS_OK(istat)) then
            call msg_toall(MSG_ERROR,'(ezgrid_serialize) Problem getting grid params')
            return
         endif
      enddo

      !Size 
      ! 1 : N = nsubgrids
      ! N * 14 integers : ni,nj,ig1, ig2, ig3, ig4, i0,j0,axni,axnj,ayni,aynj,ichar(grtyp_S),ichar(grref_S)
      ! N * ?  real(REAL64) :: ax,xy
      gridobjsize = 1 + nsubgrids * NINTEGER
      axnij = 0
      aynij = 0
      do igrid=1,nsubgrids
         if (associated(axp(igrid)%p).and.associated(ayp(igrid)%p)) then
            axnij(:,igrid) = shape(axp(igrid)%p)
            aynij(:,igrid) = shape(ayp(igrid)%p)
            gridobjsize = gridobjsize + size(axp(igrid)%p) + size(ayp(igrid)%p)
         endif
      enddo

      if (.not.associated(F_griddata)) then
         allocate(F_griddata(gridobjsize),stat=istat)
         if (istat /= 0) then
            call msg_toall(MSG_ERROR,'(ezgrid_serialize) allocate error')
            return
         endif
      endif
      if (size(F_griddata) < gridobjsize) then
         call msg_toall(MSG_ERROR,'(ezgrid_serialize) provided data array too small')
         return
      endif

      F_griddata(1) = nsubgrids
      i0 = 2
      DOIGRID: do igrid=1,nsubgrids
         istat  = ezgrid_params(subgridsid(igrid),nij,grtyp_S,grref_S,ig14,ij0)
         i00 = i0
         F_griddata(i0) = nij(1) ; i0=i0+1
         F_griddata(i0) = nij(2) ; i0=i0+1
         F_griddata(i0) = ig14(1) ; i0=i0+1
         F_griddata(i0) = ig14(2) ; i0=i0+1
         F_griddata(i0) = ig14(3) ; i0=i0+1
         F_griddata(i0) = ig14(4) ; i0=i0+1
         F_griddata(i0) = ij0(1) ; i0=i0+1
         F_griddata(i0) = ij0(2) ; i0=i0+1
         F_griddata(i0) = axnij(1,igrid) ; i0=i0+1
         F_griddata(i0) = axnij(2,igrid) ; i0=i0+1
         F_griddata(i0) = aynij(1,igrid)  ; i0=i0+1
         F_griddata(i0) = aynij(2,igrid) ; i0=i0+1
         F_griddata(i0) = iachar(grtyp_S(1:1)); i0=i0+1
         F_griddata(i0) = iachar(grref_S(1:1)) ; i0=i0+1

         if (associated(axp(igrid)%p).and.associated(ayp(igrid)%p)) then
            do j=1,axnij(2,igrid)
               do i=1,axnij(1,igrid)
                  F_griddata(i0) = transfer(axp(igrid)%p(i,j),istat) ;i0 = i0+1 
               enddo
            enddo
            do j=1,aynij(2,igrid)
               do i=1,aynij(1,igrid)
                  F_griddata(i0) = transfer(ayp(igrid)%p(i,j),istat) ;i0 = i0+1
               enddo
            enddo
            deallocate(axp(igrid)%p,stat=istat)
            deallocate(ayp(igrid)%p,stat=istat)
         endif
      enddo DOIGRID

      if (i0-1 > gridobjsize) then
         call msg_toall(MSG_ERROR,'(ezgrid_serialize) provided data array too small')
         deallocate(F_griddata,stat=istat)
         return
      endif
      F_istat = RMN_OK
      !---------------------------------------------------------------------
      return
   end function ezgrid_serialize


   !/@
   function ezgrid_unserialize(F_griddata) result(F_gridid)
      implicit none
      !@objective Create a grid from definition in a integer array
      !@arguments
      integer,intent(in) :: F_griddata(:)
      !@return
      integer :: F_gridid
      !@author 
      !@/
      integer,parameter :: NINTEGER = 14
!!$      integer,parameter :: NMAXGRIDS = 128
!!$      integer,save :: ngrids = 0 
!!$      integer,save :: grididlist(NMAXGRIDS) = -1
!!$      integer :: ,n,nij(2),ig14(4)
      integer :: gridobjsize,istat, ni,nj, ig1,ig2,ig3,ig4,i0,i00,ij0(2),i,j, &
           axnij(2),aynij(2),nsubgrids,igrid,subgridsid(NSUBGRIDS_MAX),vercode
      real :: r4 = 0.
      real,pointer :: ax(:,:),ay(:,:)
      character(len=GRTYPLEN) :: grtyp_S,grref_S
!!$      logical :: same_L
      !---------------------------------------------------------------------
      F_gridid = RMN_ERR
      if (size(F_griddata) < NINTEGER) then
         call msg_toall(MSG_ERROR,'(ezgrid_unserialize) not enough data to specify grid')
         return
      endif

      subgridsid = -1
      nsubgrids = F_griddata(1)
      i0 = 2
      DOIGRID: do igrid=1,nsubgrids

         i00 = i0
         ni = F_griddata(i0) ; i0=i0+1
         nj = F_griddata(i0) ; i0=i0+1
         ig1 = F_griddata(i0) ; i0=i0+1
         ig2 = F_griddata(i0) ; i0=i0+1
         ig3 = F_griddata(i0) ; i0=i0+1
         ig4 = F_griddata(i0) ; i0=i0+1
         ij0(1) = F_griddata(i0) ; i0=i0+1
         ij0(2) = F_griddata(i0) ; i0=i0+1
         axnij(1) = F_griddata(i0) ; i0=i0+1
         axnij(2) = F_griddata(i0) ; i0=i0+1
         aynij(1) = F_griddata(i0) ; i0=i0+1
         aynij(2) = F_griddata(i0) ; i0=i0+1
         grtyp_S = achar(F_griddata(i0)) ; i0=i0+1
         grref_S = achar(F_griddata(i0)) ; i0=i0+1

         gridobjsize = i0 + axnij(1)*axnij(2) + aynij(1)*aynij(2) - 1
         if (size(F_griddata) < gridobjsize) then
            print *,'(ezgrid_unserialize) wrong data size0: ',igrid,size(F_griddata),'<',gridobjsize 
            print *,'(ezgrid_unserialize) wrong data size1: ',i0,axnij(1),axnij(2),aynij(1),aynij(2)
            call msg_toall(MSG_ERROR,'(ezgrid_unserialize) wrong data size')
            return
         endif

         nullify(ax,ay)
         if (all((/axnij(1),axnij(2),aynij(1),aynij(2)/) > 0)) then
            allocate(ax(axnij(1),axnij(2)),ay(aynij(1),aynij(2)),stat=istat)
            if (istat /= 0) then
               call msg(MSG_ERROR,'(ezgrid_unserialize) allocate error')
               return
            endif

            do j=1,axnij(2)
               do i=1,axnij(1)
                  ax(i,j) = transfer(F_griddata(i0),r4) ; i0=i0+1
               enddo
            enddo
            do j=1,aynij(2)
               do i=1,aynij(1)
                  ay(i,j) = transfer(F_griddata(i0),r4); i0=i0+1
               enddo
            enddo
         endif
         
!!$         nij = (/ni,nj/)
!!$         ig14 = (/ig1,ig2,ig3,ig4/)
!!$         do n=1,ngrids
!!$            same_L = ezgrid_samegrid_params(grididlist(n),nij,grtyp_S,grref_S,ig14,ij0,ax,ay)
!!$            if (same_L) then
!!$               subgridsid(igrid) = grididlist(n)
!!$               if (associated(ax) .and. associated(ay)) deallocate(ax,ay,stat=istat)
!!$               cycle DOIGRID
!!$            endif
!!$         enddo

         if (associated(ax) .and. associated(ay)) then
            subgridsid(igrid) = ezgdef_fmem(ni,nj, grtyp_S(1:1), grref_S(1:1), ig1,ig2,ig3,ig4, ax, ay)
            deallocate(ax,ay,stat=istat)
         else
            if (.not.any(grtyp_S(1:1)==EZGRID_REF_TYPES)) then
               subgridsid(igrid) = ezqkdef(ni,nj, grtyp_S(1:1), ig1, ig2, ig3, ig4, 0)
            else
               call msg_toall(MSG_ERROR,'(ezgrid_unserialize) Cannot Create '//grtyp_S(1:1)//' Grid, need positional params')
               return
            endif
         endif

!!$         if (ngrids < NMAXGRIDS .and. subgridsid(igrid) > 0) then
!!$            ngrids = ngrids + 1
!!$            grididlist(ngrids) = subgridsid(igrid)
!!$         endif

      enddo DOIGRID

      if (nsubgrids == 1) then
         F_gridid = subgridsid(1)
      else
         grtyp_S = 'U' !TODO: should not be hardcoded
         grref_S = 'F' !TODO: should not be hardcoded
         vercode = 1 !TODO: should not be hardcoded
         F_gridid = ezgdef_supergrid(ni,nj,grtyp_S,grref_S,vercode,nsubgrids,subgridsid(1:nsubgrids))
      endif

      if (F_gridid < 0) then
         call msg_toall(MSG_ERROR,'(ezgrid_unserialize) Cannot Create Grid')
         return
      endif

      if (i0-1 > size(F_griddata)) then
         call msg_toall(MSG_ERROR,'(ezgrid_unserialize) read past provided data')
         F_gridid = -1
         return
      endif
      !---------------------------------------------------------------------
      return
   end function ezgrid_unserialize


   !/@
   function ezgrid_sub(F_gridid,F_i0,F_j0,F_in,F_jn,F_di,F_dj) result(F_grididout)
      implicit none
      !@objective Create a new grid with from a subset of another one
      !@arguments
      integer,intent(in) :: F_gridid,F_i0,F_in,F_j0,F_jn
      integer,intent(in),optional :: F_di,F_dj !aggregation factor
      !@return
      integer :: F_grididout
      !@author 
      !@/
      integer :: istat, nij(2), ig14(4),ij0(2),i0,j0,in,jn,di,dj,ij0b,ijnb
      real,pointer :: ax(:,:),ay(:,:)
      character(len=GRTYPLEN) :: grtyp_S,grref_S
      !---------------------------------------------------------------------
      !TODO: abort if u-grid

      call msg_toall(MSG_DEBUG,'(ezgrid) ezgrid_sub [BGN]')
      F_grididout = RMN_ERR
      nullify(ax,ay)
      istat  = ezgrid_params(F_gridid,nij,grtyp_S,grref_S,ig14,ij0,ax,ay)
      if (.not.RMN_IS_OK(istat)) return
      if (.not.any(grtyp_S(1:1) == EZGRID_REF_TYPES)) then
         call msg(MSG_WARNING,'(ezgrid_sub) None referenced grid are not yet supported')
         return
      endif
      di = 1 ; dj = 1
      if (present(F_di)) di = max(1,F_di)
      if (present(F_dj)) dj = max(1,F_dj)
      if (associated(ax) .and. associated(ay)) then
         if (grtyp_S(1:1) == 'Z') then !TODO-later: support #,Y grids
            i0 = min(max(1,F_i0),nij(1))
            in = max(i0,min(F_in,nij(1)))
            j0 = min(max(1,F_j0),nij(2))
            jn = max(j0,min(F_jn,nij(2)))
            if (di+dj > 2) then
               ij0b = 1 ; ijnb = 1
               call aggreg(ax,i0,ij0b,in,ijnb,di,1)
               call aggreg(ay,ij0b,j0,ijnb,jn,1,dj)
            endif
            F_grididout = ezgdef_fmem(in-i0+1,jn-j0+1,grtyp_S(1:1),grref_S(1:1),&
                 ig14(1),ig14(2),ig14(3),ig14(4), ax(i0:in,:), ay(:,j0:jn))
         else
            call msg(MSG_ERROR,'(ezgrid_sub) #/Y grids not yet supported')
         endif
         deallocate(ax,ay,stat=istat)
      else
         call msg(MSG_ERROR,'(ezgrid_sub) Problem getting grid axes')
      endif
      call msg_toall(MSG_DEBUG,'(ezgrid) ezgrid_sub [END]')
      !---------------------------------------------------------------------
      return
   end function ezgrid_sub


   !/@
   function ezgrid_merge(F_gridid,F_comm_S,F_bcast2all_L,F_i0,F_j0,F_lni,F_lnj) result(F_grididout)
      implicit none
      !@objective Merge subgrid/tile from all MPI PE beloing to Communicator F_comm
      !@arguments
      integer,intent(in) :: F_gridid
      character(len=*) :: F_comm_S
      logical,optional :: F_bcast2all_L
      integer,optional :: F_i0,F_j0,F_lni,F_lnj !- # tile pos and size
      !@return
      integer :: F_grididout
      !@author 
      !@/
      logical :: bcast2all_L,ismaster_L
      integer :: istat,istat2,err,nij(2),nijb(2),ig14(4),ij0(2),ij0b(2),i0,j0,lni,lnj
      character(len=GRTYPLEN) :: grtyp_S,grref_S
      real,pointer :: ax(:,:),ay(:,:),ax2(:,:),ay2(:,:)
      !---------------------------------------------------------------------
      !TODO: abort if u-grid

      call msg_toall(MSG_DEBUG,'(ezgrid) ezgrid_merge [BGN]')
      F_grididout = RMN_ERR
      bcast2all_L = .false.
      if (present(F_bcast2all_L)) bcast2all_L = F_bcast2all_L
      nullify(ax,ay)
      istat = ezgrid_params(F_gridid,nij,grtyp_S,grref_S,ig14,ij0,ax,ay)
      i0 = 1 ; j0 = 1
      if (present(F_i0)) i0 = min(F_i0,nij(1))
      if (present(F_j0)) j0 = min(F_j0,nij(2))
      lni = nij(1) - i0 + 1
      lnj = nij(2) - j0 + 1
      if (present(F_lni)) lni = min(F_lni,lni)
      if (present(F_lnj)) lnj = min(F_lnj,lnj)
      istat = min(ptopo_collect_dims(F_comm_S,lni,lnj,nijb(1),nijb(2),ij0b(1),ij0b(2)),istat)
      if (.not.(associated(ax).and.associated(ay))) then
         !TODO-later: add support for not referenced grid
         call msg(MSG_ERROR,'(ezgrid) Merging non referenced grids is not supported yet')
         istat = RMN_ERR
      endif
      call rpn_comm_allreduce(istat,istat2,1,RPN_COMM_INTEGER,"MPI_MIN",F_comm_S,err)
      if (.not.RMN_IS_OK(istat2)) return

      ismaster_L = ptopo_ismaster_L(F_comm_S)
      !TODO-later: support 2d ax (Y grids)
      nullify(ax2,ay2)
      if (ismaster_L) allocate(ax2(nijb(1),1),ay2(1,nijb(2)),stat=istat)
      ax(1:lni,1) = ax(i0:i0+lni-1,1)
      ay(1,1:lnj) = ay(1,j0:j0+lnj-1)

      istat = ptopo_collect(ax2,ax,F_comm_S,ij0b(1),1,lni,1) 
      istat = min(ptopo_collect(ay2,ay,F_comm_S,1,ij0b(2),1,lnj),istat)

!!$      if (ismaster_L.and.associated(ay2)) then
!!$         do j=1,nijb(1)-1
!!$            if (ax2(j,1) > ax2(j+1,1)) istat=RMN_ERR
!!$         enddo
!!$         do j=1,nijb(2)-1
!!$            if (ay2(1,j) > ay2(1,j+1)) istat=RMN_ERR
!!$         enddo
!!$      endif
!!$      call handle_error(istat,'ezgrid_merge','wrong ax,ay')

      if (RMN_IS_OK(istat)) then
         if (ismaster_L) then
            if (associated(ax2).and.associated(ay2)) then
               F_grididout = ezgdef_fmem(nijb(1),nijb(2), grtyp_S, grref_S, ig14(1),ig14(2),ig14(3),ig14(4), ax2, ay2)
            else
               istat = RMN_ERR
            endif
         endif
      endif
      call rpn_comm_allreduce(istat,istat2,1,RPN_COMM_INTEGER,"MPI_MIN",F_comm_S,err)
      if (RMN_IS_OK(istat2)) then
         if (bcast2all_L) then
            F_grididout = ezgrid_bcast(F_grididout,F_comm_S,RPN_COMM_MASTER)
         else
            if (.not.ismaster_L) F_grididout = RMN_OK
         endif
      endif
      if (associated(ax)) deallocate(ax,stat=istat)
      if (associated(ay)) deallocate(ay,stat=istat)
      if (associated(ax2)) deallocate(ax2,stat=istat)
      if (associated(ay2)) deallocate(ay2,stat=istat)
      call msg_toall(MSG_DEBUG,'(ezgrid) ezgrid_merge [END]')
      !---------------------------------------------------------------------
      return
   end function ezgrid_merge


   !/@
   function ezgrid_addperiod(F_gridid,F_periodx_L,F_periody_L) result(F_gridid2)
      implicit none
      !@objective Return grid parameters
      !@arguments
      integer,intent(in) :: F_gridid
      logical,intent(in) :: F_periodx_L,F_periody_L
      !@author
      !@return
      integer :: F_gridid2
      !@description
      !@/
      integer :: istat,gnij(2),ig14(4)
      character(len=2) :: grtyp_S,grref_S
      real,pointer :: ax(:,:),ay(:,:)
      real,pointer :: ax2(:,:),ay2(:,:)
      !---------------------------------------------------------------------
      F_gridid2 = RMN_ERR
      nullify(ax,ay)
      istat = ezgrid_params(F_gridid,gnij,grtyp_S,grref_S,ig14,F_ax=ax,F_ay=ay)
      if (all(grtyp_S(1:1) /= EZGRID_REF_TYPES)) then
         F_gridid2 = F_gridid
         return
      endif
      if (.not.(associated(ax).and.associated(ay))) return
      if (all(grref_S(1:1) /= (/'E','G','L'/))) return
      ax2 => ax
      ay2 => ay
      if (F_periodx_L) then
         if (.not.(all(ax(1,:) == (ax(gnij(1),:)+360.)))) then
            allocate(ax2(gnij(1)+1,size(ax,2)))
            ax2(1:gnij(1),:) = ax(1:gnij(1),:)
            ax2(gnij(1)+1,:) = ax(1,:) + 360.
         endif
      endif
      if (F_periody_L) then
         call msg(MSG_WARNING,'(ezgrid_addperiod) periody not yet supported')
         return
      endif
      if (size(ax)/=size(ax2) .or. size(ay)/=size(ay2)) then
         F_gridid2 = ezgdef_fmem(size(ax2,1),size(ay2,2), grtyp_S(1:1), grref_S(1:1), ig14(1),ig14(2),ig14(3),ig14(4), ax2, ay2)
      endif
      if (.not.associated(ax2,ax)) deallocate(ax2,stat=istat)
      if (.not.associated(ay2,ay)) deallocate(ay2,stat=istat)
      if (associated(ax)) deallocate(ax,stat=istat)
      if (associated(ay)) deallocate(ay,stat=istat)
      !---------------------------------------------------------------------
      return
   end function ezgrid_addperiod

end module ezgrid_mod
