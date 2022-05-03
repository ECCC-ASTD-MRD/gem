!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------

module phyfoldmeta_mod
   use phygridmap, only: phy_lcl_ni, phy_lcl_nj, phydim_ni, phydim_nj, phydim_nk, ijdrv_phy
   use phy_typedef
   use phygetmetaplus_mod, only: phymetaplus, phygetmetaplus
   implicit none
   private
   public :: phyfoldmeta, phyfoldmeta2d, phyfoldmeta3d, phyfoldmeta3dicat, phyfoldmeta4d, phyfold

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>

   interface phyfoldmeta
      module procedure phyfoldmeta2d
      module procedure phyfoldmeta3d
      module procedure phyfoldmeta3dicat
      module procedure phyfoldmeta4d
   end interface phyfoldmeta

   interface phyfold
      module procedure phyfold2d
      module procedure phyfold3d
      module procedure phyfold3dicat
      module procedure phyfold4d
   end interface phyfold

contains

   !/@*
   function phyfoldmeta4d(F_src, F_ijk0, F_ijk1, F_metaplus) result(F_istat)
      implicit none
      !@object Transfer data to p_runlgt space
      !@params
      integer,intent(in) :: F_ijk0(4), F_ijk1(4)
      real,   pointer :: F_src(:,:,:,:)
      type(phymetaplus), intent(in) :: F_metaplus
      integer :: F_istat
      !@author Michel Desgagne  -   summer 2013
      !*@/
      integer :: i, j, k, ik, m
      type(phymeta) :: mymeta
      !---------------------------------------------------------------
      F_istat = RMN_ERR
      mymeta = F_metaplus%meta

      ! Bound checking
      if (any(F_ijk0(:) < 1) .or. &
           F_ijk1(1) > phy_lcl_ni .or. &
           F_ijk1(2) > phy_lcl_nj .or. &
           F_ijk1(3) > mymeta%n(3) .or. &
           (F_ijk1(4)-1)*phydim_nk > mymeta%n(3)) then
         call msg(MSG_WARNING,'(phyfold) Out of bounds for '//&
              trim(mymeta%vname)//' on '//trim(mymeta%bus))
         return
      endif

      if ( F_ijk0(1) /= 1 .or. F_ijk0(2) /= 1 .or. &
           F_ijk1(1) /= phy_lcl_ni .or. &
           F_ijk1(2) /= phy_lcl_nj) then
         call msg(MSG_WARNING,'(phyfold) Horizontal sub domaine Not yet supported')
         return
      endif

      ! Transfer from the 3D source grid into the physics folded space
!$omp parallel private(i,j,k,ik,m)
!$omp do
      do j = 1, phydim_nj
         do m = F_ijk0(4), F_ijk1(4)
            do k = F_ijk0(3), F_ijk1(3)
               do i = 1, phydim_ni
                  ik = (m-1)*phydim_ni*phydim_nk + (k-1)*phydim_ni + i
                  F_metaplus%vptr(ik,j) = F_src(ijdrv_phy(1,i,j),ijdrv_phy(2,i,j),k,m)
               end do
            end do
         end do
      end do
!$omp end do
!$omp end parallel

      F_istat = RMN_OK
      !---------------------------------------------------------------
      return
   end function phyfoldmeta4d

   !/@*
   function phyfoldmeta3dicat(F_src, F_ijk0, F_ijk1, F_icat, F_metaplus) result(F_istat)
      implicit none
      !@object Transfer data to p_runlgt space
      !@params
      integer,intent(in) :: F_ijk0(3), F_ijk1(3), F_icat
      real,   pointer :: F_src(:,:,:)
      type(phymetaplus), intent(in) :: F_metaplus
      integer :: F_istat
      !@author Michel Desgagne  -   summer 2013
      !*@/
      integer :: i, j, k, ik, m
      type(phymeta) :: mymeta
      !---------------------------------------------------------------
      F_istat = RMN_ERR
      mymeta = F_metaplus%meta

      ! Bound checking
      if (any(F_ijk0(:) < 1) .or. &
           F_ijk1(1) > phy_lcl_ni .or. &
           F_ijk1(2) > phy_lcl_nj .or. &
           F_ijk1(3) > mymeta%n(3) .or. &
           (F_icat-1)*phydim_nk > mymeta%n(3)) then
         call msg(MSG_WARNING,'(phyfold) Out of bounds for '//&
              trim(mymeta%vname)//' on '//trim(mymeta%bus))
         return
      endif

      if ( F_ijk0(1) /= 1 .or. F_ijk0(2) /= 1 .or. &
           F_ijk1(1) /= phy_lcl_ni .or. &
           F_ijk1(2) /= phy_lcl_nj) then
         call msg(MSG_WARNING,'(phyfold) Horizontal sub domaine Not yet supported')
         return
      endif

      ! Transfer from the 3D source grid into the physics folded space
      m = (F_icat-1)*phydim_ni*phydim_nk
!$omp parallel private(i,j,k,ik)
!$omp do
      do j = 1, phydim_nj
         do k = F_ijk0(3), F_ijk1(3)
            do i = 1, phydim_ni
               ik = m + (k-1)*phydim_ni + i
               F_metaplus%vptr(ik,j) = F_src(ijdrv_phy(1,i,j),ijdrv_phy(2,i,j),k)
            end do
         end do
      end do
!$omp end do
!$omp end parallel

      F_istat = RMN_OK
      !---------------------------------------------------------------
      return
   end function phyfoldmeta3dicat

   !/@*
   function phyfoldmeta3d(F_src, F_ijk0, F_ijk1, F_metaplus) result(F_istat)
      implicit none
      !@object Transfer data to p_runlgt space
      !@params
      integer,intent(in) :: F_ijk0(3), F_ijk1(3)
      real,   pointer :: F_src(:,:,:)
      type(phymetaplus), intent(in) :: F_metaplus
      integer :: F_istat
      !@author Michel Desgagne  -   summer 2013
      !*@/
      integer :: i, j, k, ik
      type(phymeta) :: mymeta
      !---------------------------------------------------------------
      F_istat = RMN_ERR
      mymeta = F_metaplus%meta

      ! Bound checking
      if (any(F_ijk0(:) < 1) .or. &
           F_ijk1(1) > phy_lcl_ni .or. &
           F_ijk1(2) > phy_lcl_nj .or. &
           F_ijk1(3) > mymeta%n(3)) then
         call msg(MSG_WARNING,'(phyfold) Out of bounds for '//&
              trim(mymeta%vname)//' on '//trim(mymeta%bus))
         return
      endif

      if ( F_ijk0(1) /= 1 .or. F_ijk0(2) /= 1 .or. &
           F_ijk1(1) /= phy_lcl_ni .or. &
           F_ijk1(2) /= phy_lcl_nj) then
         call msg(MSG_WARNING,'(phyfold) Horizontal sub domaine Not yet supported')
         return
      endif

      ! Transfer from the 3D source grid into the physics folded space
!$omp parallel private(i,j,k,ik)
!$omp do
      do j = 1, phydim_nj
         do k = F_ijk0(3), F_ijk1(3)
            do i = 1, phydim_ni
               ik = (k-1)*phydim_ni + i
               F_metaplus%vptr(ik,j) = F_src(ijdrv_phy(1,i,j),ijdrv_phy(2,i,j),k)
            end do
         end do
      end do
!$omp end do
!$omp end parallel

      F_istat = RMN_OK
      !---------------------------------------------------------------
      return
   end function phyfoldmeta3d


   !/@*
   function phyfoldmeta2d(F_src, F_ijk0, F_ijk1, F_metaplus) result(F_istat)
      implicit none
      !@object Transfer data to p_runlgt space
      !@params
      integer,intent(in) :: F_ijk0(2), F_ijk1(2)
      real,   intent(in) :: F_src(:,:)
      type(phymetaplus), intent(in) :: F_metaplus
      integer :: F_istat
      !@author Michel Desgagne  -   summer 2013
      !*@/
      integer :: i, j
      type(phymeta) :: mymeta
      !---------------------------------------------------------------
      F_istat = RMN_ERR
      mymeta = F_metaplus%meta

      ! Bound checking
      if (any(F_ijk0(:) < 1) .or. &
           F_ijk1(1) > phy_lcl_ni .or. &
           F_ijk1(2) > phy_lcl_nj) then
         call msg(MSG_WARNING,'(phyfold) Out of bounds for '//&
              trim(mymeta%vname)//' on '//trim(mymeta%bus))
         return
      endif

      if ( F_ijk0(1) /= 1 .or. F_ijk0(2) /= 1 .or. &
           F_ijk1(1) /= phy_lcl_ni .or. &
           F_ijk1(2) /= phy_lcl_nj) then
         call msg(MSG_WARNING,'(phyfold) Horizontal sub domaine Not yet supported')
         return
      endif

      ! Transfer from the 3D source grid into the physics folded space
!$omp parallel private(i,j)
!$omp do
      do j = 1, phydim_nj
         do i = 1, phydim_ni
            F_metaplus%vptr(i,j) = F_src(ijdrv_phy(1,i,j),ijdrv_phy(2,i,j))
         end do
      end do
!$omp end do
!$omp end parallel

      F_istat = RMN_OK
      !---------------------------------------------------------------
      return
   end function phyfoldmeta2d


   !/@*
   function phyfold2d(F_src, F_nomvar_S, F_bus_S, F_ijk0, F_ijk1) result(F_istat)
      implicit none
      !@object Transfer data to p_runlgt space
      !@params
      character(len=1),intent(in) :: F_bus_S
      character(len=*),intent(in) :: F_nomvar_S
      integer,intent(in) :: F_ijk0(2), F_ijk1(2)
      real,intent(in) :: F_src(:,:)
      integer :: F_istat
      !@author Michel Desgagne  -   summer 2013
      !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>
      type(phymetaplus) :: mymetaplus
      !---------------------------------------------------------------
      F_istat = phygetmetaplus(mymetaplus, F_nomvar_S, F_npath='V', &
           F_bpath=F_bus_S, F_quiet=.true., F_shortmatch=.false.)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(phyfold) No matching bus entry for '// &
              trim(F_nomvar_S)//' on '//trim(F_bus_S))
         return
      endif

      F_istat = phyfoldmeta2d(F_src, F_ijk0, F_ijk1, mymetaplus)
      !---------------------------------------------------------------
      return
   end function phyfold2d

   !/@*
   function phyfold3d(F_src, F_nomvar_S, F_bus_S, F_ijk0, F_ijk1) result(F_istat)
      implicit none
      !@object Transfer data to p_runlgt space
      !@params
      character(len=1),intent(in) :: F_bus_S
      character(len=*),intent(in) :: F_nomvar_S
      integer,intent(in) :: F_ijk0(3), F_ijk1(3)
      real,pointer :: F_src(:,:,:)
      integer :: F_istat
      !@author Michel Desgagne  -   summer 2013
      !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>
      type(phymetaplus) :: mymetaplus
      !---------------------------------------------------------------
      F_istat = phygetmetaplus(mymetaplus, F_nomvar_S, F_npath='V', &
           F_bpath=F_bus_S, F_quiet=.true., F_shortmatch=.false.)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(phyfold) No matching bus entry for '// &
              trim(F_nomvar_S)//' on '//trim(F_bus_S))
         return
      endif

      F_istat = phyfoldmeta3d(F_src, F_ijk0, F_ijk1, mymetaplus)
      !---------------------------------------------------------------
      return
   end function phyfold3d

   !/@*
   function phyfold3dicat(F_src, F_nomvar_S, F_bus_S, F_ijk0, F_ijk1, F_icat) result(F_istat)
      implicit none
      !@object Transfer data to p_runlgt space
      !@params
      character(len=1),intent(in) :: F_bus_S
      character(len=*),intent(in) :: F_nomvar_S
      integer,intent(in) :: F_ijk0(3), F_ijk1(3), F_icat
      real,pointer :: F_src(:,:,:)
      integer :: F_istat
      !@author Michel Desgagne  -   summer 2013
      !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>
      type(phymetaplus) :: mymetaplus
      !---------------------------------------------------------------
      F_istat = phygetmetaplus(mymetaplus, F_nomvar_S, F_npath='V', &
           F_bpath=F_bus_S, F_quiet=.true., F_shortmatch=.false.)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(phyfold) No matching bus entry for '// &
              trim(F_nomvar_S)//' on '//trim(F_bus_S))
         return
      endif

      F_istat = phyfoldmeta3dicat(F_src, F_ijk0, F_ijk1, F_icat, mymetaplus)
      !---------------------------------------------------------------
      return
   end function phyfold3dicat

   !/@*
   function phyfold4d(F_src, F_nomvar_S, F_bus_S, F_ijk0, F_ijk1) result(F_istat)
      implicit none
      !@object Transfer data to p_runlgt space
      !@params
      character(len=1),intent(in) :: F_bus_S
      character(len=*),intent(in) :: F_nomvar_S
      integer,intent(in) :: F_ijk0(4), F_ijk1(4)
      real,pointer :: F_src(:,:,:,:)
      integer :: F_istat
      !@author Michel Desgagne  -   summer 2013
      !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>
      type(phymetaplus) :: mymetaplus
      !---------------------------------------------------------------
      F_istat = phygetmetaplus(mymetaplus, F_nomvar_S, F_npath='V', &
           F_bpath=F_bus_S, F_quiet=.true., F_shortmatch=.false.)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(phyfold) No matching bus entry for '// &
              trim(F_nomvar_S)//' on '//trim(F_bus_S))
         return
      endif

      F_istat = phyfoldmeta4d(F_src, F_ijk0, F_ijk1, mymetaplus)
      !---------------------------------------------------------------
      return
   end function phyfold4d

end module phyfoldmeta_mod
