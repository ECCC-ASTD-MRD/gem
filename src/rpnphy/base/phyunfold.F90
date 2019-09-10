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

module phyunfoldmeta_mod
   use phygridmap, only: phy_lcl_ni, phy_lcl_nj, phydim_ni, phydim_nj, phydim_nk, ijdrv_phy
   use phy_typedef
   use phygetmetaplus_mod, only: phymetaplus, phygetmetaplus
   implicit none
   private
   public :: phyunfoldmeta, phyunfoldmeta2d, phyunfoldmeta3d, phyunfoldmeta4d, phyunfold

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>

   interface phyunfoldmeta
      module procedure phyunfoldmeta2d
      module procedure phyunfoldmeta3d
      module procedure phyunfoldmeta4d
   end interface phyunfoldmeta

   interface phyunfold
      module procedure phyunfold2d
      module procedure phyunfold3d
      module procedure phyunfold4d
   end interface phyunfold

contains

   !/@*
   function phyunfoldmeta4d(F_dest, F_ijk0, F_ijk1, F_metaplus) result(F_istat)
      implicit none
      !@object Transfer data from p_runlgt space
      !@params
      integer,intent(in) :: F_ijk0(4), F_ijk1(4)
      real,intent(out) :: F_dest(:,:,:,:)
      type(phymetaplus), intent(in) :: F_metaplus
      integer :: F_istat
      !@author Michel Desgagne  -   summer 2013
      !*@/
      integer :: i, j, k, ik, k1, m, m1
      type(phymeta) :: mymeta
      !---------------------------------------------------------------
      call msg(MSG_DEBUG, '(phyunfoldmeta3d) [BEGIN] '//trim(F_metaplus%meta%vname))
      F_istat = RMN_ERR
      mymeta = F_metaplus%meta

      ! Bound checking
      if (any(F_ijk0(:) < 1) .or. &
           F_ijk1(1) > phy_lcl_ni .or. &
           F_ijk1(2) > phy_lcl_nj .or. &
           F_ijk1(3)+(F_ijk1(4)-1)*phydim_nk > mymeta%n(3)) then
         call msg(MSG_WARNING,'(phyunfold) Out of bounds for '//&
              trim(mymeta%vname)//' on '//trim(mymeta%bus))
         print *,'Requested:',F_ijk1
         print *,'Phy      :',phy_lcl_ni,phy_lcl_nj,mymeta%n(3),':',mymeta%nk,mymeta%fmul
         return
      endif

      if ( F_ijk0(1) /= 1 .or. F_ijk0(2) /= 1 .or. &
           F_ijk1(1) /= phy_lcl_ni .or. &
           F_ijk1(2) /= phy_lcl_nj) then
         call msg(MSG_WARNING,'(phyunfold) Horizontal sub domaine Not yet supported')
         return
      endif

      ! Transfer from the physics folded space into the 3D dest grid
!$omp parallel private(i,j,k,ik,k1,m,m1)
!$omp do
      do j = 1, phydim_nj
         do m = F_ijk0(4), F_ijk1(4)
            do k = F_ijk0(3), F_ijk1(3)
               do i = 1, phydim_ni
                  ik = (m-1)*phydim_ni*phydim_nk + (k-1)*phydim_ni + i
                  k1 = 1+k-F_ijk0(3)
                  m1 = 1+m-F_ijk0(4)
                  F_dest(ijdrv_phy(1,i,j),ijdrv_phy(2,i,j),k1,m1) = F_metaplus%vptr(ik,j)
               end do
            end do
         end do
      end do
!$omp end do
!$omp end parallel

      F_istat = RMN_OK
      call msg(MSG_DEBUG, '(phyunfoldmeta3d) [END] '//trim(F_metaplus%meta%vname))
      !---------------------------------------------------------------
      return
   end function phyunfoldmeta4d


   !/@*
   function phyunfoldmeta3d(F_dest, F_ijk0, F_ijk1, F_metaplus) result(F_istat)
      implicit none
      !@object Transfer data from p_runlgt space
      !@params
      integer,intent(in) :: F_ijk0(3), F_ijk1(3)
      real,intent(out) :: F_dest(:,:,:)
      type(phymetaplus), intent(in) :: F_metaplus
      integer :: F_istat
      !@author Michel Desgagne  -   summer 2013
      !*@/
      integer :: i, j, k, ik, k1
      type(phymeta) :: mymeta
      !---------------------------------------------------------------
      call msg(MSG_DEBUG, '(phyunfoldmeta3d) [BEGIN] '//trim(F_metaplus%meta%vname))
      F_istat = RMN_ERR
      mymeta = F_metaplus%meta

      ! Bound checking
      if (any(F_ijk0(:) < 1) .or. &
           F_ijk1(1) > phy_lcl_ni .or. &
           F_ijk1(2) > phy_lcl_nj .or. &
           F_ijk1(3) > mymeta%n(3)) then
         call msg(MSG_WARNING,'(phyunfold) Out of bounds for '//&
              trim(mymeta%vname)//' on '//trim(mymeta%bus))
        return
      endif

      if ( F_ijk0(1) /= 1 .or. F_ijk0(2) /= 1 .or. &
           F_ijk1(1) /= phy_lcl_ni .or. &
           F_ijk1(2) /= phy_lcl_nj) then
         call msg(MSG_WARNING,'(phyunfold) Horizontal sub domaine Not yet supported')
         return
      endif

      ! Transfer from the physics folded space into the 3D dest grid
!$omp parallel private(i,j,k,ik)
!$omp do
      do j = 1, phydim_nj
         do k = F_ijk0(3), F_ijk1(3)
            do i = 1, phydim_ni
               ik = (k-1)*phydim_ni + i
               k1 = 1+k-F_ijk0(3)
               F_dest(ijdrv_phy(1,i,j),ijdrv_phy(2,i,j),k1) = F_metaplus%vptr(ik,j)
            end do
         end do
      end do
!$omp end do
!$omp end parallel

      F_istat = RMN_OK
      call msg(MSG_DEBUG, '(phyunfoldmeta3d) [END] '//trim(F_metaplus%meta%vname))
      !---------------------------------------------------------------
      return
   end function phyunfoldmeta3d

   !/@*
   function phyunfoldmeta2d(F_dest, F_ijk0, F_ijk1, F_metaplus) result(F_istat)
      implicit none
      !@object Transfer data from p_runlgt space
      !@params
      integer,intent(in) :: F_ijk0(2), F_ijk1(2)
      real,intent(out) :: F_dest(:,:)
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
         call msg(MSG_WARNING,'(phyunfold) Out of bounds for '//&
              trim(mymeta%vname)//' on '//trim(mymeta%bus))
        return
      endif

      if ( F_ijk0(1) /= 1 .or. F_ijk0(2) /= 1 .or. &
           F_ijk1(1) /= phy_lcl_ni .or. &
           F_ijk1(2) /= phy_lcl_nj) then
         call msg(MSG_WARNING,'(phyunfold) Horizontal sub domaine Not yet supported')
         return
      endif

      !#TODO: adapt for 4D vars (ni,nk*fmul,nj)
      ! Transfer from the physics folded space into the 3D dest grid
!$omp parallel private(i,j)
!$omp do
      do j = 1, phydim_nj
         do i = 1, phydim_ni
            F_dest(ijdrv_phy(1,i,j),ijdrv_phy(2,i,j)) = F_metaplus%vptr(i,j)
         end do
      end do
!$omp end do
!$omp end parallel

      F_istat = RMN_OK
      !---------------------------------------------------------------
      return
   end function phyunfoldmeta2d
 

   !/@*
   function phyunfold4d(F_dest, F_nomvar_S, F_bus_S, F_ijk0, F_ijk1) result(F_istat)
      implicit none
      !@object Transfer data from p_runlgt space
      !@params
      character(len=1),intent(in) :: F_bus_S
      character(len=*),intent(in) :: F_nomvar_S
      integer,intent(in) :: F_ijk0(4), F_ijk1(4)
      real,intent(out) :: F_dest(:,:,:,:)
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
         call msg(MSG_WARNING,'(phyunfold) No matching bus entry for '// &
              trim(F_nomvar_S)//' on '//trim(F_bus_S))
         return
      endif

      F_istat = phyunfoldmeta4d(F_dest, F_ijk0, F_ijk1, mymetaplus)
      !---------------------------------------------------------------
      return
   end function phyunfold4d


   !/@*
   function phyunfold3d(F_dest, F_nomvar_S, F_bus_S, F_ijk0, F_ijk1) result(F_istat)
      implicit none
      !@object Transfer data from p_runlgt space
      !@params
      character(len=1),intent(in) :: F_bus_S
      character(len=*),intent(in) :: F_nomvar_S
      integer,intent(in) :: F_ijk0(3), F_ijk1(3)
      real,intent(out) :: F_dest(:,:,:)
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
         call msg(MSG_WARNING,'(phyunfold) No matching bus entry for '// &
              trim(F_nomvar_S)//' on '//trim(F_bus_S))
         return
      endif

      F_istat = phyunfoldmeta3d(F_dest, F_ijk0, F_ijk1, mymetaplus)
      !---------------------------------------------------------------
      return
   end function phyunfold3d


   !/@*
   function phyunfold2d(F_dest, F_nomvar_S, F_bus_S, F_ijk0, F_ijk1) result(F_istat)
      implicit none
      !@object Transfer data from p_runlgt space
      !@params
      character(len=1),intent(in) :: F_bus_S
      character(len=*),intent(in) :: F_nomvar_S
      integer,intent(in) :: F_ijk0(2), F_ijk1(2)
      real,intent(out) :: F_dest(:,:)
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
         call msg(MSG_WARNING,'(phyunfold) No matching bus entry for '// &
              trim(F_nomvar_S)//' on '//trim(F_bus_S))
         return
      endif

      F_istat = phyunfoldmeta2d(F_dest, F_ijk0, F_ijk1, mymetaplus)
      !---------------------------------------------------------------
      return
   end function phyunfold2d

end module phyunfoldmeta_mod
