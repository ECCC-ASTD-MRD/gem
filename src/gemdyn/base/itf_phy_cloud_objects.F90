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

module itf_phy_cloud_objects
   private
#include <rmnlib_basics.hf>

   ! API parameters
   integer, public :: CLDOBJ_OK=0,CLDOBJ_ERROR=1

   ! API subprograms
   public :: cldobj_displace
   public :: cldobj_expand

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !**s/r cldobj_displace - Displace objects using forward trajectories: CMC/RPN package
   integer function cldobj_displace () result (F_status)
      use glb_ld, only: l_minx,l_maxx,l_miny,l_maxy
      use HORgrid_options, only: Grd_lphy_i0,Grd_lphy_in,Grd_lphy_j0,Grd_lphy_jn
      use dcst, only: Dcst_inv_rayt_8
      use geomh, only: geomh_invcy_8,geomh_y_8,geomh_hy_8,geomh_x_8,geomh_hx_8,geomh_yv_8,geomh_xu_8
      use cstv, only: Cstv_dt_8
      use phy_itf, only: phy_get,phy_put
      implicit none
#include <arch_specific.hf>

      ! Local variable declarations
      integer :: i,j,newi,newj,i0,in,j0,jn
      real :: newrlat,newrlon
      integer, dimension(l_minx:l_maxx,l_miny:l_maxy) :: srci,srcj
      real, dimension(l_minx:l_maxx,l_miny:l_maxy), target :: merged,rlat,rlon,xrnd,yrnd
      real, dimension(l_minx:l_maxx,l_miny:l_maxy), target :: uadv,vadv,objrlat,objrlon,objage,objmerge,objzlcl
      logical :: newsrc

      ! Set error return status
      F_status = CLDOBJ_ERROR


      ! Retrieve advecting winds
      if (get_obj_property(uadv,'coadvu') /= CLDOBJ_OK) return
      if (get_obj_property(vadv,'coadvv') /= CLDOBJ_OK) return

      ! Retrieve current object positions
      if (get_obj_property(objrlat,'colat',F_exchange=.false.) /= CLDOBJ_OK) return
      if (get_obj_property(objrlon,'colon',F_exchange=.false.) /= CLDOBJ_OK) return

      ! Compute positions of newly formed objects in rotated coordinates
      if (get_obj_property(objage,'coage') /= CLDOBJ_OK) return
      ! Random number generator is topology-dependent... use grid cell centres
!!$      if (get_obj_property(xrnd,'rand') /= CLDOBJ_OK
!!$      call random_number(yrnd(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn))
      xrnd = 0.5
      yrnd = 0.5
      do j=Grd_lphy_j0,Grd_lphy_jn
         do i=Grd_lphy_i0,Grd_lphy_in
            if (abs(objage(i,j)-Cstv_dt_8) < epsilon(1.)) then
               objrlat(i,j) = geomh_y_8(j) + (yrnd(i,j)-0.5)*geomh_hy_8
               objrlon(i,j) = geomh_x_8(i) + (xrnd(i,j)-0.5)*geomh_hx_8
            end if
         end do
      end do

      ! Distribute object status and position information to neighbours
      if (exchange_field(objrlat) /= CLDOBJ_OK) return
      if (exchange_field(objrlon) /= CLDOBJ_OK) return

      ! Retrieve primary object property field
      if (get_obj_property(objmerge,'cowlcl') /= CLDOBJ_OK) return

      ! Compute object displacements
      if (domain_bounds(i0,in,j0,jn) /= CLDOBJ_OK) return
      srci = i0-1
      do j=j0,jn
         do i=i0,in
            if (objage(i,j) < Cstv_dt_8) cycle
            newrlat = objrlat(i,j) + vadv(i,j)*Cstv_dt_8*Dcst_inv_rayt_8
            newrlon = objrlon(i,j) + uadv(i,j)*Cstv_dt_8*Dcst_inv_rayt_8*geomh_invcy_8(j)
            newj = j0
            do while (newj < jn .and. newrlat > geomh_yv_8(newj))
               newj = newj + 1
            end do
            newi = i0
            do while (newi < in .and. newrlon > geomh_xu_8(newi))
               newi = newi + 1
            end do
            newsrc = .false.
            if (srci(newi,newj) < i0) then
               newsrc = .true.
            else if (objmerge(i,j) > merged(newi,newj)) then
               newsrc = .true.
            end if
            if (newsrc) then
               srci(newi,newj) = i
               srcj(newi,newj) = j
               merged(newi,newj) = objmerge(i,j)
               rlat(newi,newj) = newrlat
               rlon(newi,newj) = newrlon
            end if
         end do
      end do

      ! Retrieve object property fields
      if (get_obj_property(objzlcl,'cozlcl') /= CLDOBJ_OK) return

      ! Update object properties
      objrlat = 0.
      objrlon = 0.
      do j=Grd_lphy_j0,Grd_lphy_jn
         do i=Grd_lphy_i0,Grd_lphy_in
            if (srci(i,j) >= i0) then
               objrlat(i,j) = rlat(i,j)
               objrlon(i,j) = rlon(i,j)
            end if
         end do
      end do
      if (morph_field_to_physics(objmerge,'cowlcl',srci,srcj,i0) /= CLDOBJ_OK) return
      if (morph_field_to_physics(objage,'coage',srci,srcj,i0) /= CLDOBJ_OK) return
      if (morph_field_to_physics(objzlcl,'cozlcl',srci,srcj,i0) /= CLDOBJ_OK) return

      ! Ship object information back to the physics
      if (push_to_physics(objrlat,'colat') /= CLDOBJ_OK) return
      if (push_to_physics(objrlon,'colon') /= CLDOBJ_OK) return

      ! Successful completion
      F_status = CLDOBJ_OK

   end function cldobj_displace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !**s/r cldobj_expand - Expand each objects' radius of influence: CMC/RPN package
   integer function cldobj_expand (F_rmax, F_rtau) result (F_status)
      use glb_ld, only: l_minx,l_maxx,l_miny,l_maxy,l_nk
      use HORgrid_options, only: Grd_lphy_i0,Grd_lphy_in,Grd_lphy_j0,Grd_lphy_jn
      use dcst, only: Dcst_rayt_8,Dcst_inv_rayt_8
      use geomh, only: geomh_invcy_8,geomh_cy_8,geomh_y_8,geomh_hy_8,geomh_x_8,geomh_hx_8
      use cstv, only: Cstv_dt_8
      use phy_itf, only: phy_get,phy_put
      implicit none
#include <arch_specific.hf>

      ! Input arguments
      real, intent(in) :: F_rmax                                !Maximum radius of influence (m)
      real, intent(in) :: F_rtau                                !Time scale for radius increase (s)

      ! Local parameters
      real, parameter :: NRAD=2.                                !Multiple of radius to which to extend calculations
      real, parameter :: WLCL_SCALE=0.005                         !Multiplicative scaling of wlcl -> wklcl

      ! Local variable declarations
      integer :: i,j,k,is,js,i0,in,j0,jn,hw_i,hw_j
      real :: rad,di,dj,dist
      real, dimension(l_minx:l_maxx,l_miny:l_maxy), target :: objrlat,objrlon,objage,objwlcl,wklcl_adj
      real, dimension(:,:,:), pointer :: wklcl

      ! Set error return status
      F_status = CLDOBJ_ERROR

      ! Set computational domain bounds
      if (domain_bounds(i0,in,j0,jn) /= CLDOBJ_OK) return

      ! Retrieve current object positions
      if (get_obj_property(objrlat,'colat') /= CLDOBJ_OK) return
      if (get_obj_property(objrlon,'colon') /= CLDOBJ_OK) return

      ! Retrieve object property fields
      if (get_obj_property(objwlcl,'cowlcl') /= CLDOBJ_OK) return
      if (get_obj_property(objage,'coage') /= CLDOBJ_OK) return

      ! Compute trigger field adjustment near existing convection
      wklcl_adj = 0.
      do j=j0,jn
         do i=i0,in
            if (objage(i,j) < Cstv_dt_8) cycle
            ! Compute radius of influence
            rad = F_rmax * (1.-exp(-objage(i,j)/F_rtau))
            ! Find half-width of subdomain on the grid
            hw_i = ceiling(NRAD*rad*Dcst_inv_rayt_8*geomh_invcy_8(j) / geomh_hx_8)
            hw_j = ceiling(NRAD*rad*Dcst_inv_rayt_8 / geomh_hy_8)
            do js=max(j-hw_j,j0),min(j+hw_j,jn)
               do is=max(i-hw_i,i0),min(i+hw_i,in)
                  ! Compute distance of each point on subdomain from cloud centre
                  di = (objrlon(i,j)-geomh_x_8(is))*Dcst_rayt_8*geomh_cy_8(js)
                  dj = (objrlat(i,j)-geomh_y_8(js))*Dcst_rayt_8
                  dist = sqrt(di**2 + dj**2)
                  ! Adjust trigger velocity using a Gaussian radial decay
                  wklcl_adj(is,js) = min(wklcl_adj(is,js), -WLCL_SCALE*objwlcl(i,j)*exp(-0.5*(dist/rad)**2))
               end do
            end do
         end do
      end do

      ! Adjust convective trigger field
      nullify(wklcl)
      if (.not.RMN_IS_OK(phy_get(wklcl,'tr/wt:P',F_npath='V',F_bpath='D',F_quiet=.false.))) return
      do k=1,l_nk
         wklcl(:,:,k) = wklcl(:,:,k) + wklcl_adj(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn)
      end do
      if (.not.RMN_IS_OK(phy_put(wklcl,'tr/wt:P',F_npath='V',F_bpath='D',F_quiet=.false.))) return

      ! Successful completion
      F_status = CLDOBJ_OK
   end function cldobj_expand

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Update field by "morphing" from source points
   integer function morph_field_to_physics(F_fld,F_name,F_srci,F_srcj,F_i0) result(F_status)
      use HORgrid_options, only: Grd_lphy_j0,Grd_lphy_jn,Grd_lphy_i0,Grd_lphy_in
      use glb_ld, only: l_minx,l_maxx,l_miny,l_maxy
      implicit none
      real, dimension(l_minx:l_maxx,l_miny:l_maxy), intent(in) :: F_fld
      character(len=*), intent(in) :: F_name
      integer, intent(in) :: F_i0
      integer, dimension(l_minx:l_maxx,l_miny:l_maxy), intent(in) :: F_srci,F_srcj
      integer :: i,j
      real, dimension(l_minx:l_maxx,l_miny:l_maxy) :: newfld

      ! Set error return status
      F_status = CLDOBJ_ERROR

      ! Transform field following trajectory
      newfld = 0.
      do j=Grd_lphy_j0,Grd_lphy_jn
         do i=Grd_lphy_i0,Grd_lphy_in
            if (F_srci(i,j) > F_i0) then
               newfld(i,j) = F_fld(F_srci(i,j),F_srcj(i,j))
            end if
         end do
      end do

      ! Push revised object field to the physics
      if (push_to_physics(newfld,F_name) /= CLDOBJ_OK) return

      ! Successful completion
      F_status = CLDOBJ_OK

   end function morph_field_to_physics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Retrieve a property of the physics object
   integer function get_obj_property(F_fld,F_name,F_exchange) result(F_status)
      use HORgrid_options, only: Grd_lphy_j0,Grd_lphy_jn,Grd_lphy_i0,Grd_lphy_in
      use glb_ld, only: l_minx,l_maxx,l_miny,l_maxy
      use phy_itf, only: phy_get
      implicit none
      real, dimension(l_minx:l_maxx,l_miny:l_maxy), target, intent(out) :: F_fld
      character(len=*), intent(in) :: F_name
      logical, intent(in), optional :: F_exchange
      real, dimension(:,:), pointer :: fld
      logical :: exchange

      ! Set error return status
      F_status = CLDOBJ_ERROR

      ! Process optional arguments
      exchange = .true.
      if (present(F_exchange)) exchange = F_exchange

      ! Retrieve object from the permanent bus
      F_fld = 0.
      fld => F_fld(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn)
      if (.not.RMN_IS_OK(phy_get(fld,F_name,F_npath='V',F_bpath='P',F_quiet=.false.))) return

      ! Distribute primary object property field to neighbours
      if (exchange) then
         if (exchange_field(F_fld) /= CLDOBJ_OK) return
      end if

      ! Successful completion
      F_status = CLDOBJ_OK
   end function get_obj_property

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Push a field to the permanent bus
   integer function push_to_physics(F_fld,F_name) result(F_status)
      use HORgrid_options, only: Grd_lphy_j0,Grd_lphy_jn,Grd_lphy_i0,Grd_lphy_in
      use glb_ld, only: l_minx,l_maxx,l_miny,l_maxy
      use phy_itf, only: phy_put
      implicit none
      real, dimension(l_minx:l_maxx,l_miny:l_maxy), target, intent(in) :: F_fld
      character(len=*), intent(in) :: F_name
      real, dimension(:,:), pointer :: fld

      ! Set error return status
      F_status = CLDOBJ_ERROR

      ! Push field section to the permanent bus
      fld => F_fld(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn)
      if (.not.RMN_IS_OK(phy_put(fld,F_name,F_npath='V',F_bpath='P',F_quiet=.false.))) return

      ! Successful completion
      F_status = CLDOBJ_OK
   end function push_to_physics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Exchange a field with neighbours
   integer function exchange_field(F_fld) result(F_status)
      use gem_options, only: G_halox,G_haloy
      use glb_ld, only: l_ni,l_nj,G_periodx,G_periody,l_minx,l_maxx,l_miny,l_maxy
      use HORgrid_options, only : Grd_yinyang_L
      implicit none
      real, dimension(l_minx:l_maxx,l_miny:l_maxy), target, intent(in) :: F_fld

      ! Set error return status
      F_status = CLDOBJ_ERROR

      ! Distribute object property field to neighbours
      if (Grd_yinyang_L) then
         call yyg_xchng(F_fld,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,1,.false., &
              'NEAREST',.true.)
      else
         call rpn_comm_xch_halo(F_fld,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,1, &
           G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
      end if
      ! Successful completion
      F_status = CLDOBJ_OK
   end function exchange_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute maximum region for valid calculations
   integer function domain_bounds(F_i0,F_in,F_j0,F_jn) result(F_status)
      use gem_options, only: G_halox,G_haloy
      use glb_ld, only: l_ni,l_nj,pil_w,pil_e,pil_n,pil_s
      integer, intent(out) :: F_i0,F_in,F_j0,F_jn

      ! Set error return status
      F_status = CLDOBJ_ERROR

      ! Compute maximum region for calculations
      F_i0 = 1 - G_halox + pil_w
      F_in = l_ni + G_halox - pil_e
      F_j0 = 1 - G_haloy + pil_s
      F_jn = l_nj + G_haloy - pil_n

      ! Successful completion
      F_status = CLDOBJ_OK
   end function domain_bounds

end module itf_phy_cloud_objects
