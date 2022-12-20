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

   subroutine itf_phy_geom (F_istat)
   use dcst
   use lam_options
   use gmm_table
   use gmm_geof
   use geomh
   use glb_ld
   use mem_nest
   use tdpack
   use, intrinsic :: iso_fortran_env
   implicit none
#include <arch_specific.hf>

   integer F_istat

#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   logical :: nest_it
   integer :: i,j,istat
   real, pointer :: wrk1(:,:)
   real(kind=REAL64), parameter :: deg2rad_8 = pi_8 / 180.d0
   type(gmm_metadata) :: mymeta
!
!-------------------------------------------------------------------
!
   F_istat = RMN_OK

   mymeta = GMM_NULL_METADATA
   mymeta%l(1) = gmm_layout(1,l_ni,0,0,l_ni)
   mymeta%l(2) = gmm_layout(1,l_nj,0,0,l_nj)

   nullify(wrk1)
   istat = gmm_create('DLAT',wrk1,mymeta,GMM_FLAG_RSTR)
   gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)='DLAT'; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='SF' ; GMM_tbl%fst(gmm_cnt)='LAT'
   if (RMN_IS_OK(istat)) then
      wrk1 = deg2rad_8*geomh_latrx
   else
      F_istat = RMN_ERR
      call msg(MSG_ERROR,'(itf_phy_geom) Problem creating DLAT')
   end if

   nullify(wrk1)
   istat = gmm_create('DLON',wrk1,mymeta,GMM_FLAG_RSTR)
   gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)='DLON'; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='SF' ; GMM_tbl%fst(gmm_cnt)='LON'
   if (RMN_IS_OK(istat)) then
      where(geomh_lonrx >= 0)
         wrk1 = deg2rad_8*geomh_lonrx
      elsewhere
         wrk1 = deg2rad_8*(geomh_lonrx+360.)
      endwhere
   else
      F_istat = RMN_ERR
      call msg(MSG_ERROR,'(itf_phy_geom) Problem creating DLON')
   end if

   nullify(wrk1)
   istat = gmm_create('DXDY',wrk1,mymeta,GMM_FLAG_RSTR)
   gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)='DXDY'; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='SF' ; GMM_tbl%fst(gmm_cnt)='DXDY'
   if (RMN_IS_OK(istat)) then
      do j = 1,l_nj
         do i = 1,l_ni
             wrk1(i,j) = geomh_hx_8 * geomh_hy_8 * Dcst_rayt_8**2 *geomh_cy_8(j)
         end do
      end do
   else
      F_istat = RMN_ERR
      call msg(MSG_ERROR,'(itf_phy_geom) Problem creating DXDY')
   end if

   nullify(wrk1)
   istat = gmm_create('TDMASK',wrk1,mymeta,GMM_FLAG_RSTR)
   gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)='TDMASK'; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='SF' ; GMM_tbl%fst(gmm_cnt)='PMSK'
   if (RMN_IS_OK(istat)) then
      wrk1(1:l_ni,1:l_nj) = 1.
      nest_it = ( Lam_0ptend_L .and. &
                ((Lam_blend_Hx > 0).or.(Lam_blend_Hy > 0)) )
      if ( nest_it ) then
         do j=1,l_nj
            do i=1,l_ni
               wrk1(i,j)= 1.-nest_weightm(i,j,G_nk+1)
            enddo
         enddo
      end if
   else
      F_istat= RMN_ERR
      call msg(MSG_ERROR,'(itf_phy_geom) Problem creating TDMASK')
   end if
!
!-------------------------------------------------------------------
!
   return
   end subroutine itf_phy_geom
