
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
       subroutine dfilter(F_iname,F_oname,F_nk)
      use HORgrid_options, only: Grd_lphy_i0,Grd_lphy_in,Grd_lphy_j0,Grd_lphy_jn
      use gem_options
      use itf_phy_filter, only: sfcflxfilt_o, sfcflxfilt_i
      use phy_itf, only: phy_get,phy_put
      use glb_ld
      use lun
      implicit none
#include <rmnlib_basics.hf>


      ! Local variables
      character(len=*), intent(in) :: F_iname                   !Name of the field to smooth
      character(len=*), intent(in) :: F_oname                   !Name of the field in which to store smoothed result
      integer , intent(in) :: F_nk

      ! Local variables
      real, dimension(:,:,:), pointer :: fld3d
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,F_nk), target :: fld3ds
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,F_nk), target :: ifld3d
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,F_nk) :: ifld
      integer :: istat, it, np

!
      if (Lun_debug_L) then
         if (Lun_out > 0) write( Lun_out,1000) F_iname
      endif

      if(sfcflxfilt_o /=2 .and. sfcflxfilt_o /=4) then
         write(Lun_out, 2000)
         call gem_error(-1,'dfilter','the order of the digital filter should be 2 or 4')
     endif

      np=sfcflxfilt_o/2
         ifld3d = 0.
         fld3d => ifld3d(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:F_nk)
         istat = phy_get(fld3d,F_iname,F_npath='V',F_bpath='V',F_quiet=.false.)
         call gem_error(istat,'dfilter','problem in phy_get: '//trim(F_iname))

        ! Apply sfcflxfilt_i passess of smoothing operator
         ifld=0.
         ifld=ifld3d
         do it=1,sfcflxfilt_i
            if (it >1)then
               ifld=0.
               ifld(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,:)=  &
                                fld3ds(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,:)
            endif
            ! Shapiro Filter
            call  shapiro_filter(fld3ds,ifld,F_nk,np)
         enddo
            fld3d =>fld3ds(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:F_nk)
            istat = phy_put(fld3d,F_oname,F_npath='V',F_bpath='P',F_quiet=.false.)
            call gem_error(istat,'dfilter','problem in phy_put: '//trim(F_oname))

1000 format(/' Apply digital filter to surface flux ',a)
2000 format(/' S/R  dfilter : The order of the digital filter should be 2 or 4')

         return

      contains

!! Shapiro filter
      subroutine  shapiro_filter(F_ofld,F_ifld,F_nk,F_np)
      use gem_options
      use HORgrid_options
      use geomh
      use glb_ld

      ! - Apply a Shapiro filter to the selected field
      implicit none

      ! Common block variables
#include <rmnlib_basics.hf>

      ! Argument declaration
      integer, intent(in) :: F_nk, F_np
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,F_nk), intent(out) :: F_ofld !Smoothed output field
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,F_nk), intent(in) :: F_ifld                      !Input field to smooth

      ! Local variables
      integer :: i,j,k
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,F_nk) :: ifld
      real(kind=REAL64), parameter ::a10=1.d0/2.d0 , a11=1.d0/4.d0
      real(kind=REAL64), parameter ::a20=10.d0/16.d0 , a21=4.d0/16.d0, a22=-1.d0/16.d0

      ! Exchange halos for smoothing operation
      ifld = F_ifld

      call rpn_comm_xch_halo(ifld,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,F_nk, &
           G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

      if(F_np==2) then
         ! Apply smoothing  X direction
         F_ofld=0.
         do k=1,F_nk
            do j=1,l_nj
               do i=1,l_ni
                        F_ofld(i,j,k) = a20*ifld(i,j,k) + a21*(ifld(i-1,j,k) + ifld(i+1,j,k)) &
                                                        + a22*(ifld(i-2,j,k) + ifld(i+2,j,k))
               end do
            end do
         end do


         ! Exchange halos
         call rpn_comm_xch_halo(F_ofld,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,F_nk, &
                                G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

         ! Apply smoothing  Y direction
         do k=1,F_nk
            do j=1,l_nj
               do i=1,l_ni
                        F_ofld(i,j,k) = a20*F_ofld(i,j,k) + a21*(F_ofld(i,j-1,k) + F_ofld(i,j+1,k)) &
                                                          + a22*(F_ofld(i,j-2,k) + F_ofld(i,j+2,k))
               end do
            end do
         end do
      elseif(F_np==1) then
         ! Apply smoothing  X direction
         F_ofld=0.
         do k=1,F_nk
            do j=1,l_nj
               do i=1,l_ni
                        F_ofld(i,j,k) = a10*ifld(i,j,k) + a11*(ifld(i-1,j,k) + ifld(i+1,j,k))
               end do
            end do
         end do

      ! Exchange halos
         call rpn_comm_xch_halo(F_ofld,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,F_nk, &
                                G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

         ! Apply smoothing  Y direction
         do k=1,F_nk
            do j=1,l_nj
               do i=1,l_ni
                        F_ofld(i,j,k) = a10*F_ofld(i,j,k) + a11*(F_ofld(i,j-1,k) + F_ofld(i,j+1,k))
               end do
            end do
         end do
      endif

      return

    end subroutine  shapiro_filter

  end subroutine dfilter
