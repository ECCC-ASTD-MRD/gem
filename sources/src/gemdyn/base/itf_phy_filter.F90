!---------------------------------- LICENCE BEGIN -----------------------------
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
!---------------------------------- LICENCE END -------------------------------
!/@*
module itf_phy_filter
   use HORgrid_options, only: Grd_lphy_i0,Grd_lphy_in,Grd_lphy_j0,Grd_lphy_jn
   ! Provide spatial filtering options for physics variables
   private

   ! External decalarations
#include <msg.h>

   ! Module variables
   integer, save :: hw = 0                           !Filter half-width (gridpoints)
   real, pointer, save :: wt_gauss1(:,:,:) => NULL() !Filter weights for each point
   real, pointer, save :: wt_gauss2(:,:,:) => NULL() !Filter weights for each point
   logical, save :: initialized = .false.            !Status indicators for package

   ! Public API
   public :: ipf_init
   public :: ipf_smooth_fld
   public :: ipf_smooth_tend

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function ipf_init(F_sig, F_sig2) result(F_status)
      use gem_options, only: G_halox,G_haloy
      ! Prepare filter weights for smoothing
      implicit none

      ! Common block variables
#include <rmnlib_basics.hf>

      ! Argument declaration
      real, intent(in), optional :: F_sig               !Standard deviation of Gaussian filter (increase for more smoothing) [1 gridpoint]
      real, intent(in), optional :: F_sig2              !Standard deviation of Gaussian filter (increase for more smoothing) [1 gridpoint]

      ! Local variables
      real :: mysig1, mysig2

      ! Initialize return status
      F_status = RMN_ERR

      ! Handle optional arguments
      mysig1 = 1.
      if (present(F_sig)) mysig1 = F_sig
      mysig2 = -1.
      if (present(F_sig2)) mysig2 = F_sig2

      ! Confirm package initialization status
      if (initialized) then
         call msg(MSG_WARNING,'(itf_phy_filter::ipf_init) Attempted re-initialization of package')
         return
      end if
      initialized = .true.

      ! No weighting if no smoothing is to be done
      if (mysig1 <= 0. .and. mysig2 <= 0.) then
         F_status = RMN_OK
         return
      end if

      ! Set filter half-width
      hw = min(G_halox,G_haloy)

      ! Initialize all filters
      call ipf_init_stencil(mysig1, wt_gauss1)
      call ipf_init_stencil(mysig2, wt_gauss2)

      ! End of subprogram
      F_status = RMN_OK
      return
   end function ipf_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine ipf_init_stencil(F_sig, F_weight)
      use geomh
      use glb_ld
      use cstv
      use dcst !#TODO: use TDPACK?
      implicit none

      ! Argument declaration
      real, intent(in) :: F_sig               !Standard deviation of Gaussian filter (increase for more smoothing) [1 gridpoint]
      real, pointer :: F_weight(:,:,:)        !Filter weight stencil

      ! Local variables
      integer :: i,j,ii,jj,cnt
      real :: d_euclid,d_sphere,sigm,wt_sum

      nullify(F_weight)
      if (F_sig <= 0.) return

      ! Create space for filter weights
      allocate(F_weight(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,(2*hw+1)**2))

      ! Apply smoothing
      sigm = F_sig*geomh_hy_8*Dcst_rayt_8       !standard deviation converted to meters
      do j=Grd_lphy_j0,Grd_lphy_jn
         do i=Grd_lphy_i0,Grd_lphy_in
            wt_sum = 0.; cnt = 0
            do jj=-hw,hw
               do ii=-hw,hw
                  d_euclid = Dcst_rayt_8*sqrt((geomh_x_8(i)-geomh_x_8(i+ii))**2+(geomh_y_8(j)-geomh_y_8(j+jj))**2)
                  d_sphere = Dcst_rayt_8*asin(d_euclid/(2.*Dcst_rayt_8**2)*sqrt(4.*Dcst_rayt_8**2-d_euclid**2))
                  cnt = cnt+1
                  F_weight(i,j,cnt) = exp(-d_sphere**2/(2.*sigm**2))
                  wt_sum = wt_sum + F_weight(i,j,cnt)
               end do
            end do
            F_weight(i,j,:) = F_weight(i,j,:) / wt_sum
         end do
      end do

      return
   end subroutine ipf_init_stencil

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function ipf_smooth_tend(F_data,F_name,F_fid) result(F_status)
      use phy_itf, only: phy_get,phy_put
      use HORgrid_options
      use glb_ld
      use cstv
      use gmm_itf_mod
      ! Smooth selected scheme-based tendency from full physics tendency.
      implicit none

      ! Common block variables
#include <rmnlib_basics.hf>

      ! Argument declaration
      real, dimension(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,G_nk), intent(inout) :: F_data    !Full physics tendency on physics grid
      character(len=*), intent(in) :: F_name                                                            !Name of V-bus tendency variable to smooth
      integer, intent(in), optional :: F_fid                                                            !Filter ID (currently 1 or 2)

      ! Local variables
      real, dimension(:,:), pointer :: tdmask,tdmask1
      real, dimension(:,:,:), pointer :: tend
      real, dimension(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,G_nk), target :: tends
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,G_nk), target :: data_tend
      integer :: iend(3), fid, k
      logical :: apply_filter

      ! Initialize return status
      F_status = RMN_ERR

      ! Check for a valid stencil
      fid = 1
      if (present(F_fid)) fid = F_fid
      apply_filter = .false.
      if (fid == 1 .and. associated(wt_gauss1)) apply_filter = .true.
      if (fid == 2 .and. associated(wt_gauss2)) apply_filter = .true.

      ! Only apply smoothing if necessary
      if (.not.apply_filter) then
         F_status = RMN_OK
         return
      end if

      ! Retrieve tendency mask
      nullify(tdmask, tdmask1)
      F_status = gmm_get('TDMASK',tdmask)
      if (F_status < 0 .or. .not.associated(tdmask)) then
         call msg(MSG_WARNING,'(itf_phy_filter::ipf_smooth_tend) Problem getting TDMASK - Call to gaussian_filter failed')
         return
      end if
      tdmask1 => tdmask(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn)

      ! Retrieve selected physics tendency
      data_tend = 0.
      tend => data_tend(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:l_nk)
      iend = (/-1,-1,l_nk/)
      F_status = phy_get(tend,F_name,F_npath='V',F_bpath='V',F_end=iend,F_quiet=.false.)
      if (F_status < 0) return

      ! Remove selected tendency from total physics tendency
      do k=1,G_nk
         F_data(:,:,k) = F_data(:,:,k) - tend(:,:,k)*Cstv_dt_8*tdmask1(:,:)
      end do

      ! Apply smoothing operator
      if (gaussian_filter(tends,data_tend,fid) /= RMN_OK) then
         call msg(MSG_WARNING,'(itf_phy_filter::ipf_smooth_tend) Call to gaussian_filter failed')
         do k=1,G_nk
            F_data(:,:,k) = F_data(:,:,k) + tend(:,:,k)*Cstv_dt_8*tdmask1(:,:)
         end do
         return
      end if

      ! Add smoothed tendency to total physics tendency
      do k=1,G_nk
         F_data(:,:,k) = F_data(:,:,k) + tends(:,:,k)*Cstv_dt_8*tdmask1(:,:)
      end do
      tend => tends
      F_status = phy_put(tend,F_name,F_npath='V',F_bpath='VP',F_end=(/-1,-1,l_nk/))

      ! End of subprogram
      F_status = RMN_OK
      return
   end function ipf_smooth_tend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function ipf_smooth_fld(F_iname,F_oname,F_fid) result(F_status)
      use phy_itf, only: phy_get,phy_put
      use HORgrid_options
      use glb_ld
      ! Smooth selected field and place it in a separate bus entry
      implicit none

      ! Common block variables
#include <rmnlib_basics.hf>

      ! Argument declaration
      character(len=*), intent(in) :: F_iname                   !Name of the field to smooth
      character(len=*), intent(in) :: F_oname                   !Name of the field in which to store smoothed result
      integer, intent(in), optional :: F_fid                    !Filter ID (currently 1 or 2)

      ! Local variables
      real, dimension(:,:,:), pointer :: fld
      real, dimension(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,G_nk), target :: flds
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,G_nk), target :: ifld
      integer :: iend(3), fid
      logical :: apply_filter

      ! Initialize return status
      F_status = RMN_ERR

      ! Check for a valid stencil
      fid = 1
      if (present(F_fid)) fid = F_fid
      apply_filter = .false.
      if (fid == 1 .and. associated(wt_gauss1)) apply_filter = .true.
      if (fid == 2 .and. associated(wt_gauss2)) apply_filter = .true.

      ! Only apply smoothing if necessary
      if (.not.apply_filter) then
         F_status = RMN_OK
         return
      end if

      ! Retrieve selected input field
      ifld = 0.
      fld => ifld(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:l_nk)
      iend = (/-1,-1,l_nk/)
      F_status = phy_get(fld,F_iname,F_npath='V',F_bpath='PVD',F_end=iend,F_quiet=.false.)
      if (F_status < 0) return

      ! Apply smoothing operator
      if (gaussian_filter(flds,ifld,fid) /= RMN_OK) then
         call msg(MSG_WARNING,'(itf_phy_filter::ipf_smooth_fld) Call to gaussian_filter failed')
         return
      end if

      ! Post result to output field
      fld => flds
      F_status = phy_put(fld,F_oname,F_npath='V',F_bpath='PDV',F_end=iend)

      ! End of subprogram
      F_status = RMN_OK
   end function ipf_smooth_fld

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function gaussian_filter(F_ofld,F_ifld,F_fid) result(F_status)
      use gem_options
      use HORgrid_options
      use geomh
      use glb_ld
      use cstv
      use dcst

      ! PRIVATE - Apply a Gaussian filter to the selected field
      implicit none

      ! Common block variables
#include <rmnlib_basics.hf>

      ! Argument declaration
      real, dimension(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,G_nk), intent(out) :: F_ofld !Smoothed output field
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,G_nk), intent(in) :: F_ifld                      !Input field to smooth
      integer, intent(in) :: F_fid                                                                 !Filter ID (currently 1 or 2)

      ! Local variables
      integer :: i,j,k,ii,jj,cnt
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,G_nk) :: ifld
      real, dimension(:,:,:), pointer :: wt_gauss

      ! Initialize return status
      F_status = RMN_ERR

      ! Confirm package initialization status
      if (.not.initialized) then
         call msg(MSG_WARNING,'(itf_phy_filter::gaussian_filter) Filtering attempted before package initialization')
         F_ofld = F_ifld(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,:)
         return
      end if

      ! Check for a valid stencil
      nullify(wt_gauss)
      if (F_fid == 1 .and. associated(wt_gauss1)) wt_gauss => wt_gauss1
      if (F_fid == 2 .and. associated(wt_gauss2)) wt_gauss => wt_gauss2

      ! Only apply smoothing if necessary
      if (.not.associated(wt_gauss)) then
         F_ofld = F_ifld(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,:)
         F_status = RMN_OK
         return
      end if

      ! Exchange halos for smoothing operation
      ifld = F_ifld
      if (Grd_yinyang_L) then
         call yyg_xchng(ifld,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk,&
                        .false.,'PHYSI', .true.)
      else
      call rpn_comm_xch_halo(ifld,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk, &
           G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
      end if

      ! Apply smoothing
      F_ofld = 0.
!$omp parallel private(k,j,i,cnt,ii,jj) shared(hw,wt_gauss,ifld)
!$omp do
      do k=1,l_nk
         do j=Grd_lphy_j0,Grd_lphy_jn
            do i=Grd_lphy_i0,Grd_lphy_in
               cnt = 0
               do jj=-hw,hw
                  do ii=-hw,hw
                     cnt = cnt+1
                     F_ofld(i,j,k) = F_ofld(i,j,k) + wt_gauss(i,j,cnt)*ifld(i+ii,j+jj,k)
                 end do
               end do
            end do
         end do
      end do
!$omp enddo
!$omp end parallel

      ! End of subprogram
      F_status = RMN_OK
      return
   end function gaussian_filter

end module itf_phy_filter
