!-------------------------------------- LICENCE BEGIN ------------------------------
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
!-------------------------------------- LICENCE END -------------------------------

module phygridmap
   implicit none
   private
   public :: mapping2drivergrid, mapphy2mod, mapmod2phy
!!!#include <arch_specific.hf>

   integer, save, public :: phy_lcl_gid = 0
   integer, save, public :: phy_lcl_ni = 0
   integer, save, public :: phy_lcl_nj = 0
   integer, save, public :: phy_lcl_i0 = 0
   integer, save, public :: phy_lcl_in = 0
   integer, save, public :: phy_lcl_j0 = 0
   integer, save, public :: phy_lcl_jn = 0
   integer, save, public :: drv_glb_gid = 0
   integer, save, public :: drv_glb_ni = 0
   integer, save, public :: drv_glb_nj = 0

   integer, save, public :: drv_lcl_gid = 0
   integer, save, public :: drv_lcl_ni = 0
   integer, save, public :: drv_lcl_nj = 0

   integer, save, public :: drv_lcl_i0 = 0
   integer, save, public :: drv_lcl_in = 0
   integer, save, public :: drv_lcl_j0 = 0
   integer, save, public :: drv_lcl_jn = 0
   integer, save, public :: phydim_ni = 0
   integer, save, public :: phydim_nj = 0
   integer, save, public :: phydim_nk = 0
   integer, save, public :: phy_lclcore_gid = 0

   integer, save, public :: phy_glb_ni = 0
   integer, save, public :: phy_glb_nj = 0
   integer, save, public :: phy_glb_gid = -1
   integer, save, public :: phy_glbcore_gid = -1
   integer, save, public :: phy_comm_io_id = -1

   integer, pointer, save, public :: ijdrv_mod(:,:,:) => NULL()
   integer, pointer, save, public :: ijdrv_phy(:,:,:)  => NULL()
   integer, pointer, save, public :: ijphy(:,:,:) => NULL()

contains

   !/@*
   subroutine mapping2drivergrid
      !*@/

      integer :: i, j, p_offi, p_offj, ij(2)
      !---------------------------------------------------------------
      allocate( &
           ijdrv_mod(2, phydim_ni, phydim_nj), &
           ijdrv_phy(2, phydim_ni, phydim_nj), &
           ijphy(2, phy_lcl_ni, phy_lcl_nj) )

      p_offi = phy_lcl_i0 - 1
      p_offj = phy_lcl_j0 - 1

      do j= 1, phydim_nj
         do i= 1, phydim_ni
            ij = mapphy2mod(i, j, p_offi, p_offj, phy_lcl_ni, phy_lcl_nj, phydim_ni)
            ijdrv_mod(1:2,i,j) = ij(1:2)
         end do
      end do

      p_offi = 0
      p_offj = 0

      do j= 1, phydim_nj
         do i= 1, phydim_ni
            ij = mapphy2mod(i, j, p_offi, p_offj, phy_lcl_ni, phy_lcl_nj, phydim_ni)
            ijdrv_phy(1:2,i,j) = ij(1:2)
         end do
      end do

      do j= 1, phy_lcl_nj
         do i= 1, phy_lcl_ni
            ij = mapmod2phy(i, j, phy_lcl_ni, phydim_ni)
            ijphy(1:2,i,j) = ij(1:2)
         end do
      end do
      !---------------------------------------------------------------
      return
   end subroutine mapping2drivergrid

   !/@*
   function mapphy2mod(F_iphy, F_jphy, F_offi, F_offj, F_lcl_ni, F_lcl_nj, F_dim_ni)
      implicit none
      !@Arguments
      integer, intent(in)  :: F_iphy, F_jphy, F_offi, F_offj, F_lcl_ni, F_lcl_nj, F_dim_ni
      !@Return
      integer :: mapphy2mod(2)
      !*@/

      integer :: nphy, ijp, jmod
      !---------------------------------------------------------------
      nphy = (F_jphy-1) * F_dim_ni
      ijp  = min(nphy+F_iphy, F_lcl_ni*F_lcl_nj)
      jmod = ijp/F_lcl_ni + min(1,mod(ijp, F_lcl_ni))
      mapphy2mod(1) = ijp - (jmod-1)*F_lcl_ni + F_offi
      mapphy2mod(2) = jmod + F_offj
      return
   end function mapphy2mod


   !/@*
   function mapmod2phy(F_imod, F_jmod, F_lcl_ni, F_dim_ni)
      implicit none
      !@Arguments
      integer :: F_imod, F_jmod, F_lcl_ni, F_dim_ni
      !@Return
      integer :: mapmod2phy(2)
      !*@/
      integer :: nmod,ijm,iphy,jphy
      !---------------------------------------------------------------
      nmod = (F_jmod-1) * F_lcl_ni
      ijm  = nmod+F_imod
      jphy = ijm/F_dim_ni + min(1, mod(ijm, F_dim_ni))
      iphy = ijm - (jphy-1)*F_dim_ni
      mapmod2phy(1) = iphy
      mapmod2phy(2) = jphy
      return
   end function mapmod2phy


end module phygridmap
