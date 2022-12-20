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

      subroutine itf_phy_UVupdate
      use gmm_vt1
      use gmm_pw
      use glb_ld
      use rmn_gmm
      implicit none
#include <arch_specific.hf>

      integer istat,k
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,G_nk), target :: tdu,tdv
!
!-----------------------------------------------------------------
!
      istat = gmm_get (gmmk_ut1_s, ut1)
      istat = gmm_get (gmmk_vt1_s, vt1)
      istat = gmm_get (gmmk_pw_uu_copy_s,pw_uu_copy)
      istat = gmm_get (gmmk_pw_vv_copy_s,pw_vv_copy)
      istat = gmm_get (gmmk_pw_uu_plus_s,pw_uu_plus)
      istat = gmm_get (gmmk_pw_vv_plus_s,pw_vv_plus)

!$omp parallel private(k)
!$omp do
      do k= 1, G_nk
         tdu(l_minx:l_maxx,l_miny:0,k) = 0. ; tdu(l_minx:l_maxx,l_nj+1:l_maxy,k) = 0.
         tdv(l_minx:l_maxx,l_miny:0,k) = 0. ; tdv(l_minx:l_maxx,l_nj+1:l_maxy,k) = 0.
         tdu(l_minx:0,l_miny:l_maxy,k) = 0. ; tdu(l_ni+1:l_maxx,l_miny:l_maxy,k) = 0.
         tdv(l_minx:0,l_miny:l_maxy,k) = 0. ; tdv(l_ni+1:l_maxx,l_miny:l_maxy,k) = 0.
         tdu(1:l_ni,1:l_nj,k) = (pw_uu_plus(1:l_ni,1:l_nj,k)- &
                                 pw_uu_copy(1:l_ni,1:l_nj,k))
         tdv(1:l_ni,1:l_nj,k) = (pw_VV_plus(1:l_ni,1:l_nj,k)- &
                                 pw_VV_copy(1:l_ni,1:l_nj,k))
      end do
!$omp end do
!$omp end parallel

      call hwnd_stag ( pw_uu_copy,pw_vv_copy, tdu,tdv, &
                       l_minx,l_maxx,l_miny,l_maxy,l_nk,.true. )

!$omp parallel do  private(k)
      do k= 1, G_nk
         ut1(1:l_niu,1:l_nj,k) = ut1(1:l_niu,1:l_nj,k) + pw_uu_copy(1:l_niu,1:l_nj,k)
         vt1(1:l_ni,1:l_njv,k) = vt1(1:l_ni,1:l_njv,k) + pw_vv_copy(1:l_ni,1:l_njv,k)
      end do
!$omp end parallel do

!
!-----------------------------------------------------------------
!
   return
   end subroutine itf_phy_UVupdate
