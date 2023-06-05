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

!**s/r hzd_theta - applies horizontal diffusion on theta
!
      subroutine hzd_theta
      use hzd_exp
      use gmm_pw
      use gmm_vt1
      use tdpack
      use glb_ld
      use lun
      implicit none

      integer k
      real, parameter :: p_naught=100000., eps=1.0e-5
      real :: pres_t(l_ni,l_nj,G_nk),th(l_minx:l_maxx,l_miny:l_maxy,G_nk)
!
!-------------------------------------------------------------------
!

      do k=1,G_nk
         pres_t(1:l_ni,1:l_nj,k) = (p_naught/pw_pt_plus(1:l_ni,1:l_nj,k))**cappa_8
         th(1:l_ni,1:l_nj,k) = tt1(1:l_ni,1:l_nj,k) * pres_t(1:l_ni,1:l_nj,k)
!NOTE: painting the halo regions in order to avoid float error when dble(X)
!      under yyg_xchng: IF EVER we remove dble in yyg_xchng, we do not need this.
         th(l_minx:0     ,:     ,k) = tcdk_8
         th(l_ni+1:l_maxx,:     ,k) = tcdk_8
         th(1:l_ni,l_miny:0     ,k) = tcdk_8
         th(1:l_ni,l_nj+1:l_maxy,k) = tcdk_8
      end do

      call hzd_exp_visco ( th, 'S_THETA', l_minx,l_maxx,l_miny,l_maxy, G_nk )

      do k=1,G_nk
         tt1(1:l_ni,1:l_nj,k) = th    (1:l_ni,1:l_nj,k) / &
                                pres_t(1:l_ni,1:l_nj,k)
      end do

!
!-------------------------------------------------------------------
!
      return
      end
