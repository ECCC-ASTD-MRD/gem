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
!
      subroutine adv_param ()
      use adv_grid
      use adv_interp
      use glb_ld
      use ver
      implicit none
#include <arch_specific.hf>

   ! objective set 1-D interpolation
   !
   ! Rabah Aider   august 2017

      real*8, parameter :: LARGE_8 = 1.D20 , TWO_8 = 2.D0  , SIX_8 = 6.D0
      integer :: i, j, k, i0, j0, k0, pnx, pny
      real*8 :: ra,rb,rc,rd,rx,re
      real*8 :: prhxmn, prhymn, prhzmn, pdfi
      real*8 :: whx(G_ni+2*adv_halox)
      real*8 :: why(G_nj+2*adv_haloy)
      real*8, dimension(0:l_nk+1) :: verz_m,whzt,whzm,whzx
      real   :: alpha
      real*8, dimension(l_nk) :: verz_t,verz_x
!
!     ---------------------------------------------------------------
!
      sig=(Ver_z_8%m(l_nk)-Ver_z_8%m(1))/(abs(  Ver_z_8%m(l_nk)-Ver_z_8%m(1) ))
      verz_m(0)=sig*Ver_z_8%m(0)
      verz_m(l_nk+1)=sig*Ver_z_8%m(l_nk+1)

      do k=1,l_nk
         verz_m(k)=sig*Ver_z_8%m(k)
         verz_t(k)=sig*Ver_z_8%t(k)
         verz_x(k)=sig*Ver_z_8%x(k)
      end do

      adv_xbc_8 = 1.D0/(adv_xg_8(adv_gminx+1)-adv_xg_8(adv_gminx))
      adv_xabcd_8=-adv_xbc_8**3/SIX_8
      adv_xdabc_8=-adv_xabcd_8
      adv_xbacd_8= adv_xbc_8**3/TWO_8
      adv_xcabd_8=-adv_xbacd_8

      adv_ybc_8 = 1.D0/(adv_yg_8(adv_gminx+1)-adv_yg_8(adv_gminx))
      adv_yabcd_8=-adv_ybc_8**3/SIX_8
      adv_ydabc_8=-adv_yabcd_8
      adv_ybacd_8=adv_ybc_8**3/TWO_8
      adv_ycabd_8=-adv_ybacd_8

      adv_x00_8 = adv_xg_8(adv_gminx)
      adv_y00_8 = adv_yg_8(adv_gminy)

      prhxmn = LARGE_8 ; prhymn = LARGE_8 ; prhzmn = LARGE_8

      do i = adv_gminx,adv_gmaxx-1
         whx(adv_halox+i) = adv_xg_8(i+1) - adv_xg_8(i)
         prhxmn = min(whx(adv_halox+i), prhxmn)
      end do

      do j = adv_gminy,adv_gmaxy-1
         why(adv_haloy+j) = adv_yg_8(j+1) - adv_yg_8(j)
         prhymn = min(why(adv_haloy+j), prhymn)
      end do

! Prepare zeta on super vertical grid

      alpha=1.0
      if(sig < 0.) alpha=verz_m(0)

      whzt(0    ) = alpha; whzt(l_nk  ) = alpha; whzt(l_nk+1) = alpha
      do k = 1,l_nk-1
         whzt(k) = verz_t(k+1) - verz_t(k)
         prhzmn = min(whzt(k), prhzmn)
      end do

      whzm(0    ) = alpha ; whzm(l_nk  ) = alpha ; whzm(l_nk+1) = alpha
      do k = 1,l_nk-1
         whzm(k) = verz_m(k+1) - verz_m(k)
         prhzmn = min(whzm(k), prhzmn)
      end do

      whzx(0    ) = alpha ; whzx(l_nk  ) = alpha ; whzx(l_nk+1) = alpha
      do k = 1,l_nk-1
         whzx(k) = verz_x(k+1) - verz_x(k)
         prhzmn = min(whzx(k), prhzmn)
      end do

      adv_ovdx_8 = 1.0d0/prhxmn
      adv_ovdy_8 = 1.0d0/prhymn
      adv_ovdz_8 = 1.0d0/prhzmn

      pnx = int(1.0+(adv_xg_8(adv_gmaxx)-adv_x00_8)   *adv_ovdx_8)
      pny = int(1.0+(adv_yg_8(adv_gmaxy)-adv_y00_8)   *adv_ovdy_8)
      pnz = nint(1.0+(verz_m(l_nk+1)-verz_m(0))*adv_ovdz_8)

      allocate( adv_lcx(pnx), adv_lcy(pny),  adv_bsx_8(G_ni+2*adv_halox),&
              adv_dlx_8(G_ni+2*adv_halox), adv_bsy_8(G_nj+2*adv_haloy),  &
              adv_dly_8(G_nj+2*adv_haloy), adv_lcz%t(pnz),               &
              adv_lcz%m(pnz), adv_lcz%x(pnz),                            &
              adv_bsz_8%t(0:l_nk-1), adv_bsz_8%m(0:l_nk-1),              &
              adv_bsz_8%x(0:l_nk-1),                                     &
              adv_dlz_8%t(-1:l_nk) ,adv_dlz_8%m(-1:l_nk) ,               &
              adv_dlz_8%x(-1:l_nk) ,adv_diz_8(-1:l_nk))

      i0 = 1
      do i=1,pnx
         pdfi = adv_xg_8(adv_gminx) + (i-1) * prhxmn
         if (pdfi > adv_xg_8(i0+1-adv_halox)) i0 = min(G_ni+2*adv_halox-1,i0+1)
         adv_lcx(i) = i0
      end do
      do i = adv_gminx,adv_gmaxx-1
         adv_dlx_8(adv_halox+i) = whx(adv_halox+i)
      end do
      do i = adv_gminx,adv_gmaxx
         adv_bsx_8(adv_halox+i) = adv_xg_8(i)
      end do

      j0 = 1
      do j = 1,pny
         pdfi = adv_yg_8(adv_gminy) + (j-1) * prhymn
         if (pdfi > adv_yg_8(j0+1-adv_haloy)) j0 = min(G_nj+2*adv_haloy-1,j0+1)
         adv_lcy(j) = j0
      end do
      do j = adv_gminy,adv_gmaxy-1
         adv_dly_8(adv_haloy+j) = why(adv_haloy+j)
      end do

      do j = adv_gminy,adv_gmaxy
         adv_bsy_8(adv_haloy+j) = adv_yg_8(j)
      end do

      k0 = 1
      do k = 1,pnz
         pdfi = verz_m(0) + (k-1) * prhzmn
         if (pdfi > verz_t(k0+1)) k0 = min(l_nk-2, k0+1)
         adv_lcz%t(k) = k0
      end do

      do k = 0,l_nk+1            !! warning note the shift in k !!
         adv_dlz_8%t(k-1) =       whzt(k)
      end do

      do k = 1,l_nk
         adv_bsz_8%t(k-1) = ver_z_8%t(k)
      end do

      k0 = 1
      do k = 1,pnz
         pdfi = verz_m(0) + (k-1) * prhzmn
         if (pdfi > verz_m(k0+1)) k0 = min(l_nk-2, k0+1)
         adv_lcz%m(k) = k0
      end do

      do k = 0,l_nk+1            !! warning note the shift in k !!
         adv_dlz_8%m(k-1) =       whzm(k)
      end do

      do k = 1,l_nk
         adv_bsz_8%m(k-1) = ver_z_8%m(k)
      end do

      k0 = 1
      do k = 1,pnz
         pdfi = verz_m(0) + (k-1) * prhzmn
         if (pdfi > verz_x(k0+1)) k0 = min(l_nk-2, k0+1)
         adv_lcz%x(k) = k0
      end do

      do k = 0,l_nk+1            !! warning note the shift in k !!
         adv_dlz_8%x(k-1) =       whzx(k)
      end do

      do k = 1,l_nk
         adv_bsz_8%x(k-1) = ver_z_8%x(k)
      end do

      do k = 0,l_nk+1                  !! warning note the shift in k !
         adv_diz_8(k-1) = sig/whzm(k)
      end do

      allocate( &
      adv_zbc_8%t(l_nk),   adv_zbc_8%m(l_nk),    adv_zbc_8%x(l_nk), &
      adv_zabcd_8%t(l_nk), adv_zabcd_8%m(l_nk),  adv_zabcd_8%x(l_nk), &
      adv_zbacd_8%t(l_nk), adv_zbacd_8%m(l_nk),  adv_zbacd_8%x(l_nk), &
      adv_zcabd_8%t(l_nk), adv_zcabd_8%m(l_nk),  adv_zcabd_8%x(l_nk), &
      adv_zdabc_8%t(l_nk), adv_zdabc_8%m(l_nk),  adv_zdabc_8%x(l_nk) )

      allocate( &
      adv_zxabcde_8%t(l_nk), adv_zaxbcde_8%t(l_nk), adv_zbxacde_8%t(l_nk), &
      adv_zcxabde_8%t(l_nk), adv_zdxabce_8%t(l_nk), adv_zexabcd_8%t(l_nk))

      do k = 2,l_nk-2
         ra = Ver_z_8%m(k-1)
         rb = Ver_z_8%m(k)
         rc = Ver_z_8%m(k+1)
         rd = Ver_z_8%m(k+2)
         adv_zabcd_8%m(k) = 1.0/triprod(ra,rb,rc,rd)
         adv_zbacd_8%m(k) = 1.0/triprod(rb,ra,rc,rd)
         adv_zcabd_8%m(k) = 1.0/triprod(rc,ra,rb,rd)
         adv_zdabc_8%m(k) = 1.0/triprod(rd,ra,rb,rc)

         ra = Ver_z_8%t(k-1)
         rb = Ver_z_8%t(k)
         rc = Ver_z_8%t(k+1)
         rd = Ver_z_8%t(k+2)
         adv_zabcd_8%t(k) = 1.0/triprod(ra,rb,rc,rd)
         adv_zbacd_8%t(k) = 1.0/triprod(rb,ra,rc,rd)
         adv_zcabd_8%t(k) = 1.0/triprod(rc,ra,rb,rd)
         adv_zdabc_8%t(k) = 1.0/triprod(rd,ra,rb,rc)

         ra = Ver_z_8%x(k-1)
         rb = Ver_z_8%x(k)
         rc = Ver_z_8%x(k+1)
         rd = Ver_z_8%x(k+2)
         adv_zabcd_8%x(k) = 1.0/triprod(ra,rb,rc,rd)
         adv_zbacd_8%x(k) = 1.0/triprod(rb,ra,rc,rd)
         adv_zcabd_8%x(k) = 1.0/triprod(rc,ra,rb,rd)
         adv_zdabc_8%x(k) = 1.0/triprod(rd,ra,rb,rc)
      end do

      do k = 3,l_nk-3
         rx = Ver_z_8%t(k-2)
         ra = Ver_z_8%t(k-1)
         rb = Ver_z_8%t(k)
         rc = Ver_z_8%t(k+1)
         rd = Ver_z_8%t(k+2)
         re = Ver_z_8%t(k+3)
         adv_zxabcde_8%t(k) = 1.0/quiprod(rx,ra,rb,rc,rd,re)
         adv_zaxbcde_8%t(k) = 1.0/quiprod(ra,rx,rb,rc,rd,re)
         adv_zbxacde_8%t(k) = 1.0/quiprod(rb,rx,ra,rc,rd,re)
         adv_zcxabde_8%t(k) = 1.0/quiprod(rc,rx,ra,rb,rd,re)
         adv_zdxabce_8%t(k) = 1.0/quiprod(rd,rx,ra,rb,rc,re)
         adv_zexabcd_8%t(k) = 1.0/quiprod(re,rx,ra,rb,rc,rd)
      end do

      do k = 1,l_nk-1
         rb = Ver_z_8%m(k)
         rc = Ver_z_8%m(k+1)
         adv_zbc_8%m(k) = 1.0/(rc-rb)

         rb = Ver_z_8%t(k)
         rc = Ver_z_8%t(k+1)
         adv_zbc_8%t(k) = 1.0/(rc-rb)

         rb = Ver_z_8%x(k)
         rc = Ver_z_8%x(k+1)
         adv_zbc_8%x(k) = 1.0/(rc-rb)
      end do

!
!     ---------------------------------------------------------------
!
      return

contains

      real*8 function triprod(za,zb,zc,zd)
         real*8, intent(in) :: za,zb,zc,zd
         triprod = ((za-zb)*(za-zc)*(za-zd))
      end function triprod

      real*8 function quiprod(zx,za,zb,zc,zd,ze)
         real*8, intent(in) :: zx,za,zb,zc,zd,ze
         quiprod = ((zx-za)*(zx-zb)*(zx-zc)*(zx-zd)*(zx-ze))
      end function quiprod

      end subroutine adv_param
