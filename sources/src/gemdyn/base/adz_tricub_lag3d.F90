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

!**s/r adz_tricub_lag3d - SL tricubic interpolation except linear
!                         vertical iterpolation near vertical boundaries

      subroutine adz_tricub_lag3d (F_cub, F_in, F_xyz, F_num, F_lev_S)
      use ISO_C_BINDING
      use glb_ld
      use ver
      use adz_mem
      implicit none
#include <arch_specific.hf>

      include "tricublin_beta.inc"
      character(len=*), intent(in) :: F_lev_S
      integer, intent(in) :: F_num
      real,dimension(*), intent(in)  :: F_in
      real,dimension(*), intent(out) :: F_cub
      real,dimension(3,F_num), intent(in) :: F_xyz

      logical(C_BOOL) :: zcubic_L
      integer n, o1, indxij
      integer :: ii1, jj1, kk1, kp1
      real*8, parameter :: ov6 = 1.d0/6.d0, ov2 = 1.d0/2.d0
      real*8 , dimension(:),pointer, contiguous :: &
                         p_zabcd_8, p_zbacd_8, &
                         p_zcabd_8, p_zdabc_8, posz, dz, odz
      real*8 :: px(4),py(4),pz(4), del,del12,del10
      real*8 :: triprd, zb, zc, zd, ra, rb, rc, rd
      real za, zt
      triprd(za,zb,zc,zd)=(za-zb)*(za-zc)*(za-zd)
!
!---------------------------------------------------------------------
!
      if (F_lev_S == 'm') then
         posz => Ver_z_8%m
         dz   => Adz_delz_m
         odz  => Adz_odelz_m
         p_zabcd_8 => Adz_zabcd_8%m
         p_zbacd_8 => Adz_zbacd_8%m
         p_zcabd_8 => Adz_zcabd_8%m
         p_zdabc_8 => Adz_zdabc_8%m
      elseif  (F_lev_S == 't') then
         posz => Ver_z_8%t
         dz   => Adz_delz_t
         odz  => Adz_odelz_t
         p_zabcd_8 => Adz_zabcd_8%t
         p_zbacd_8 => Adz_zbacd_8%t
         p_zcabd_8 => Adz_zcabd_8%t
         p_zdabc_8 => Adz_zdabc_8%t
      endif

      do n=1, F_num
         ii1   = F_xyz(1,n)
         del   = F_xyz(1,n) - dble(ii1)
         del12 = (del-1)*(del-2)
         del10 = (del+1)* del
         px(1) = -del   *del12 *ov6
         px(2) = (del+1)*del12 *ov2
         px(3) = -del10*(del-2)*ov2
         px(4) =  del10*(del-1)*ov6

         jj1   = F_xyz(2,n)
         del   = F_xyz(2,n) - dble(jj1)
         del12 = (del-1)*(del-2)
         del10 = (del+1)* del
         py(1) = -del   *del12 *ov6
         py(2) = (del+1)*del12 *ov2
         py(3) = -del10*(del-2)*ov2
         py(4) =  del10*(del-1)*ov6

         ii1 = ii1 - Adz_lminx - 1
         jj1 = jj1 - Adz_lminy - 1
         indxij = (jj1-Adz_joff)*Adz_nit + (ii1-Adz_ioff) - Adz_nit

         kk1 = F_xyz(3,n)
         zt = posz(kk1) + (F_xyz(3,n)-dble(kk1))*dz(kk1)
         kk1 = min(Adz_kkmax,max(0,kk1 - 1))

         zcubic_L = ( kk1 > 0) .and. (kk1 < Adz_kkmax)

         kp1 = kk1 + 1
         if (zcubic_L) then
            ra = posz(kk1  )
            rb = posz(kp1  )
            rc = posz(kk1+2)
            rd = posz(kk1+3)
            pz(1) = triprd(zt,rb,rc,rd)*p_zabcd_8(kp1)
            pz(2) = triprd(zt,ra,rc,rd)*p_zbacd_8(kp1)
            pz(3) = triprd(zt,ra,rb,rd)*p_zcabd_8(kp1)
            pz(4) = triprd(zt,ra,rb,rc)*p_zdabc_8(kp1)
         else
            kk1   = min(l_nk-1,kk1+1)
            pz(1) = 0.d0 ; pz(4) = 0.d0
            pz(3) = (zt-posz(kk1))*odz(kk1)
            pz(2) = 1.d0 - pz(3)
         endif

         o1 = (kk1-1)*Adz_nij + indxij
         call tricublin_zyxf_beta (F_cub(n),F_in(o1), px,py,pz,&
                                   Adz_nit,Adz_nij,zcubic_L)
      enddo
!
!---------------------------------------------------------------------
!
      end subroutine adz_tricub_lag3d
