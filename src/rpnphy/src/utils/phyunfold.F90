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

module phyunfold
   use phygridmap, only: phy_lcl_ni, phy_lcl_nj, phydim_ni, phydim_nj, phydim_nk, ijdrv_phy
   use phymem, only: phymeta, phyvar, phymem_find
   implicit none
   private
   public :: phyunfoldmeta1, phyunfoldmeta2d, phyunfoldmeta3d, phyunfoldmeta4d, phyunfold1

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   interface phyunfoldmeta1
      module procedure phyunfoldmeta2d
      module procedure phyunfoldmeta3d
      module procedure phyunfoldmeta4d
   end interface phyunfoldmeta1

   interface phyunfold1
      module procedure phyunfold2d
      module procedure phyunfold3d
      module procedure phyunfold4d
   end interface phyunfold1

contains

   !/@*
   function phyunfoldmeta4d(F_dest, F_ijk0, F_ijk1, F_meta) result(F_istat)
      implicit none
      !@object Transfer data from p_runlgt space
      !@params
      integer,intent(in) :: F_ijk0(4), F_ijk1(4)
      real,intent(out) :: F_dest(:,:,:,:)
      type(phymeta), intent(in) :: F_meta
      integer :: F_istat
      !@author Michel Desgagne  -   summer 2013
      !*@/
      integer :: i, j, k, ik, k1, m, m1
      real, pointer, contiguous :: vptr(:)
      !---------------------------------------------------------------
      call msg(MSG_DEBUG, '(phyunfoldmeta3d) [BEGIN] '//trim(F_meta%vname))
      F_istat = RMN_ERR

      ! Bound checking
      if (any(F_ijk0(:) < 1) .or. &
           F_ijk1(1) > phy_lcl_ni .or. &
           F_ijk1(2) > phy_lcl_nj .or. &
           F_ijk1(3)+(F_ijk1(4)-1)*phydim_nk > F_meta%nlcl(3)) then
         call msg(MSG_WARNING,'(phyunfold) Out of bounds for '//&
              trim(F_meta%vname)//' on '//trim(F_meta%bus))
         print *,'Requested:',F_ijk1
         print *,'Phy      :',phy_lcl_ni,phy_lcl_nj,F_meta%nlcl(3),':',F_meta%nk,F_meta%fmul
         return
      endif

      if ( F_ijk0(1) /= 1 .or. F_ijk0(2) /= 1 .or. &
           F_ijk1(1) /= phy_lcl_ni .or. &
           F_ijk1(2) /= phy_lcl_nj) then
         call msg(MSG_WARNING,'(phyunfold) Horizontal sub domaine Not yet supported')
         return
      endif

      ! Transfer from the physics folded space into the 3D dest grid
!$omp parallel private(i,j,k,ik,k1,m,m1,vptr)
!$omp do
      do j = 1, phydim_nj
         vptr(1:F_meta%size) => F_meta%bptr(F_meta%i0:F_meta%in,j)
         do m = F_ijk0(4), F_ijk1(4)
            do k = F_ijk0(3), F_ijk1(3)
               do i = 1, phydim_ni
                  ik = (m-1)*phydim_ni*phydim_nk + (k-1)*phydim_ni + i
                  k1 = 1+k-F_ijk0(3)
                  m1 = 1+m-F_ijk0(4)
                  F_dest(ijdrv_phy(1,i,j),ijdrv_phy(2,i,j),k1,m1) = vptr(ik)
               end do
            end do
         end do
      end do
!$omp end do
!$omp end parallel

      F_istat = RMN_OK
      call msg(MSG_DEBUG, '(phyunfoldmeta3d) [END] '//trim(F_meta%vname))
      !---------------------------------------------------------------
      return
   end function phyunfoldmeta4d


   !/@*
   function phyunfoldmeta3d(F_dest, F_ijk0, F_ijk1, F_meta) result(F_istat)
      implicit none
      !@object Transfer data from p_runlgt space
      !@params
      integer,intent(in) :: F_ijk0(3), F_ijk1(3)
      real,intent(out) :: F_dest(:,:,:)
      type(phymeta), intent(in) :: F_meta
      integer :: F_istat
      !@author Michel Desgagne  -   summer 2013
      !*@/
      integer :: i, j, k, ik, k1
      real, pointer, contiguous :: vptr(:)
      !---------------------------------------------------------------
      call msg(MSG_DEBUG, '(phyunfoldmeta3d) [BEGIN] '//trim(F_meta%vname))
      F_istat = RMN_ERR

      ! Bound checking
      if (any(F_ijk0(:) < 1) .or. &
           F_ijk1(1) > phy_lcl_ni .or. &
           F_ijk1(2) > phy_lcl_nj .or. &
           F_ijk1(3) > F_meta%nlcl(3)) then
         call msg(MSG_WARNING,'(phyunfold) Out of bounds for '//&
              trim(F_meta%vname)//' on '//trim(F_meta%bus))
         print *,'Requested:',F_ijk1
         print *,'Phy      :',phy_lcl_ni,phy_lcl_nj,F_meta%nlcl(3),':',F_meta%nk,F_meta%fmul
         return
      endif

      if ( F_ijk0(1) /= 1 .or. F_ijk0(2) /= 1 .or. &
           F_ijk1(1) /= phy_lcl_ni .or. &
           F_ijk1(2) /= phy_lcl_nj) then
         call msg(MSG_WARNING,'(phyunfold) Horizontal sub domaine Not yet supported')
         return
      endif

      ! Transfer from the physics folded space into the 3D dest grid
!$omp parallel private(i,j,k,ik,vptr)
!$omp do
      do j = 1, phydim_nj
         vptr(1:F_meta%size) => F_meta%bptr(F_meta%i0:F_meta%in,j)
         do k = F_ijk0(3), F_ijk1(3)
            do i = 1, phydim_ni
               ik = (k-1)*phydim_ni + i
               k1 = 1+k-F_ijk0(3)
               F_dest(ijdrv_phy(1,i,j),ijdrv_phy(2,i,j),k1) = vptr(ik)
            end do
         end do
      end do
!$omp end do
!$omp end parallel

      F_istat = RMN_OK
      call msg(MSG_DEBUG, '(phyunfoldmeta3d) [END] '//trim(F_meta%vname))
      !---------------------------------------------------------------
      return
   end function phyunfoldmeta3d

   !/@*
   function phyunfoldmeta2d(F_dest, F_ijk0, F_ijk1, F_meta) result(F_istat)
      implicit none
      !@object Transfer data from p_runlgt space
      !@params
      integer,intent(in) :: F_ijk0(2), F_ijk1(2)
      real,intent(out) :: F_dest(:,:)
      type(phymeta), intent(in) :: F_meta
      integer :: F_istat
      !@author Michel Desgagne  -   summer 2013
      !*@/
      integer :: i, j
      real, pointer, contiguous :: vptr(:)
      !---------------------------------------------------------------
      F_istat = RMN_ERR

      ! Bound checking
      if (any(F_ijk0(:) < 1) .or. &
           F_ijk1(1) > phy_lcl_ni .or. &
           F_ijk1(2) > phy_lcl_nj) then
         call msg(MSG_WARNING,'(phyunfold) Out of bounds for '//&
              trim(F_meta%vname)//' on '//trim(F_meta%bus))
         print *,'Requested:',F_ijk1
         print *,'Phy      :',phy_lcl_ni,phy_lcl_nj,F_meta%nlcl(3),':',F_meta%nk,F_meta%fmul
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
!$omp parallel private(i,j,vptr)
!$omp do
      do j = 1, phydim_nj
         vptr(1:F_meta%size) => F_meta%bptr(F_meta%i0:F_meta%in,j)
         do i = 1, phydim_ni
            F_dest(ijdrv_phy(1,i,j),ijdrv_phy(2,i,j)) = vptr(i)
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
#include <rmn/msg.h>
      type(phyvar) :: myphyvar(1)
      !---------------------------------------------------------------
      F_istat = phymem_find(myphyvar, F_nomvar_S, F_npath='V', &
           F_bpath=F_bus_S, F_quiet=.true., F_shortmatch=.false.)
      if (F_istat <= 0) then
         call msg(MSG_WARNING,'(phyunfold) No matching bus entry for '// &
              trim(F_nomvar_S)//' on '//trim(F_bus_S))
         F_istat = RMN_ERR
         return
      endif

      F_istat = phyunfoldmeta4d(F_dest, F_ijk0, F_ijk1, myphyvar(1)%meta)
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
#include <rmn/msg.h>
      type(phyvar) :: myphyvar(1)
      !---------------------------------------------------------------
      F_istat = phymem_find(myphyvar, F_nomvar_S, F_npath='V', &
           F_bpath=F_bus_S, F_quiet=.true., F_shortmatch=.false.)
      if (F_istat <= 0) then
        call msg(MSG_WARNING,'(phyunfold) No matching bus entry for '// &
              trim(F_nomvar_S)//' on '//trim(F_bus_S))
         F_istat = RMN_ERR
         return
      endif

      F_istat = phyunfoldmeta3d(F_dest, F_ijk0, F_ijk1, myphyvar(1)%meta)
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
#include <rmn/msg.h>
      type(phyvar) :: myphyvar(1)
      !---------------------------------------------------------------
      F_istat = phymem_find(myphyvar, F_nomvar_S, F_npath='V', &
           F_bpath=F_bus_S, F_quiet=.true., F_shortmatch=.false.)
      if (F_istat <= 0) then
         call msg(MSG_WARNING,'(phyunfold) No matching bus entry for '// &
              trim(F_nomvar_S)//' on '//trim(F_bus_S))
         F_istat = RMN_ERR
         return
      endif

      F_istat = phyunfoldmeta2d(F_dest, F_ijk0, F_ijk1, myphyvar(1)%meta)
      !---------------------------------------------------------------
      return
   end function phyunfold2d

end module phyunfold
