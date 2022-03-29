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

!**s/r heights_uv - Compute heights on u and v points
      
      subroutine heights_uv ()
      use gmm_geof
      use tdpack
      use glb_ld
      use metric
      use cstv
      use ver
      use gem_options
      implicit none

      integer :: i,j,k
      real(kind=REAL64) :: aa,bb
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,2) :: me_uv,sls_uv
!
!     ---------------------------------------------------------------
!
!!$      ! cubic interpolation of fis0
!!$      aa = -0.0625d0
!!$      bb = +0.5625d0
!!$      do j= 1-(G_haloy-1)*south , l_nj+(G_haloy-2)*north
!!$         do i= 1-(G_halox-1)*west , l_ni+(G_halox-2)*east
!!$            me_uv(i,j,1) =  aa * (fis0(i-1,j) + fis0(i+2,j)) &
!!$                          + bb * (fis0(i  ,j) + fis0(i+1,j))
!!$            me_uv(i,j,2) =  aa * (fis0(i,j-1) + fis0(i,j+2)) &
!!$                          + bb * (fis0(i,j  ) + fis0(i,j+1))
!!$         end do
!!$      end do
!!$      call gem_xch_halo (me_uv, l_minx, l_maxx, l_miny, l_maxy, 2)
!!$      if (l_south) then
!!$         me_uv(:,1-G_haloy,1) = fis0(:,1-G_haloy)
!!$         me_uv(:,1-G_haloy,2) = fis0(:,1-G_haloy)
!!$      endif
!!$      if (l_west) then
!!$         me_uv(1-G_halox,:,1) = fis0(1-G_halox,:)
!!$         me_uv(1-G_halox,:,2) = fis0(1-G_halox,:)
!!$      endif
!!$      if (l_north) then
!!$         me_uv(:,l_nj+G_haloy-1:l_nj+G_haloy,1) = fis0(:,l_nj+G_haloy-1:l_nj+G_haloy)
!!$         me_uv(:,l_nj+G_haloy-1:l_nj+G_haloy,2) = fis0(:,l_nj+G_haloy-1:l_nj+G_haloy)
!!$      endif
!!$      if (l_east) then
!!$         me_uv(l_ni+G_halox-1:l_ni+G_halox,:,1) = fis0(l_ni+G_halox-1:l_ni+G_halox,:)
!!$         me_uv(l_ni+G_halox-1:l_ni+G_halox,:,2) = fis0(l_ni+G_halox-1:l_ni+G_halox,:)
!!$      endif
      ! linear interpolation of fis0 and sls
      aa = 0.5d0
      do j= 1-G_haloy*south , l_nj+(G_haloy-1)*north
         do i= 1-G_halox*west , l_ni+(G_halox-1)*east
            me_uv (i,j,1) =  aa * (fis0(i,j) + fis0(i+1,j))
            me_uv (i,j,2) =  aa * (fis0(i,j) + fis0(i,j+1))
            sls_uv(i,j,1) =  aa * (sls (i,j) + sls (i+1,j))
            sls_uv(i,j,2) =  aa * (sls (i,j) + sls (i,j+1))
         end do
      end do
      call gem_xch_halo (me_uv , l_minx, l_maxx, l_miny, l_maxy, 2)
      call gem_xch_halo (sls_uv, l_minx, l_maxx, l_miny, l_maxy, 2)
      if (l_north) then
         me_uv (:,l_nj+G_haloy,1) = fis0(:,l_nj+G_haloy)
         me_uv (:,l_nj+G_haloy,2) = fis0(:,l_nj+G_haloy)
         sls_uv(:,l_nj+G_haloy,1) = sls (:,l_nj+G_haloy)
         sls_uv(:,l_nj+G_haloy,2) = sls (:,l_nj+G_haloy)
      endif
      if (l_east) then
         me_uv (l_ni+G_halox,:,1) = fis0(l_ni+G_halox,:)
         me_uv (l_ni+G_halox,:,2) = fis0(l_ni+G_halox,:)
         sls_uv(l_ni+G_halox,:,1) = sls (l_ni+G_halox,:)
         sls_uv(l_ni+G_halox,:,2) = sls (l_ni+G_halox,:)
      endif
      
      GVM%ztht_u(:,:,0)=ver_z_8%m(0) ; GVM%ztht_v(:,:,0)=ver_z_8%m(0)
      GVM%zmom_u(:,:,0)=ver_z_8%m(0) ; GVM%zmom_v(:,:,0)=ver_z_8%m(0)

      do k=1,G_nk
         do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox
            GVM%zmom_u(i,j,k)=ver_z_8%m(k)+Cstv_bar1_8*(Ver_b_8%m(k)*me_uv(i,j,1)+Ver_c_8%m(k)*sls_uv(i,j,1))/grav_8
            GVM%ztht_u(i,j,k)=ver_z_8%t(k)+Cstv_bar1_8*(Ver_b_8%t(k)*me_uv(i,j,1)+Ver_c_8%t(k)*sls_uv(i,j,1))/grav_8
            GVM%zmom_v(i,j,k)=ver_z_8%m(k)+Cstv_bar1_8*(Ver_b_8%m(k)*me_uv(i,j,2)+Ver_c_8%m(k)*sls_uv(i,j,2))/grav_8
            GVM%ztht_v(i,j,k)=ver_z_8%t(k)+Cstv_bar1_8*(Ver_b_8%t(k)*me_uv(i,j,2)+Ver_c_8%t(k)*sls_uv(i,j,2))/grav_8
         end do
         end do
      end do
      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox      
            GVM%zmom_u(i,j,G_nk+1)=Cstv_bar1_8*me_uv(i,j,1)/grav_8
            GVM%ztht_u(i,j,G_nk+1)=Cstv_bar1_8*me_uv(i,j,1)/grav_8
            GVM%zmom_v(i,j,G_nk+1)=Cstv_bar1_8*me_uv(i,j,2)/grav_8
            GVM%ztht_v(i,j,G_nk+1)=Cstv_bar1_8*me_uv(i,j,2)/grav_8
         end do
      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine heights_uv
