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

!**s/r bicubHQV_lin_min_max - Horizontal:Cubic Vertical:Quintic + Store Lin/Min/Max

      subroutine bicubHQV_lin_min_max (F_qut,F_lin,F_min,F_max,F_in, &
                                       F_x,F_y,F_z,F_num,F_nind,ii,F_nk,F_lev_S)

      use adv_grid
      use adv_interp
      use dynkernel_options

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      character(len=*), intent(IN) :: F_lev_S
      integer, intent(IN) :: F_num,F_nind,F_nk
      real,dimension(F_num), intent(OUT) :: F_qut ! High-order SL solution
      real,dimension(F_num), intent(OUT) :: F_lin ! Low-order SL solution
      real,dimension(F_num), intent(OUT) :: F_min ! MIN over cell
      real,dimension(F_num), intent(OUT) :: F_max ! MAX over cell
      real,dimension(*),     intent(IN)  :: F_in  ! Field to interpolate
      real,dimension(F_num), intent(IN)  :: F_x,F_y,F_z
      integer,dimension(F_nind*4),intent(IN) :: ii

      real(kind=REAL64) , dimension(:), pointer, contiguous :: p_bsz_8,p_zbc_8,    &
                                                               p_zabcd_8,p_zbacd_8,&
                                                               p_zcabd_8,p_zdabc_8,&
                                                               p_zxabcde_8,p_zaxbcde_8,p_zbxacde_8,&
                                                               p_zcxabde_8,p_zdxabce_8,p_zexabcd_8

      integer :: n0, nx, ny, nz, m1, o1, o2, o3, o4, kkmax, n, id
      real(kind=REAL64)  :: a1, a2, a3, a4, b1, b2, b3, b4, &
                 c1, c2, c3, c4, d1, d2, d3, d4, &
                 p1, p2, p3, p4, ra, rb, rc, rd, &
                 x1, x2, x3, x4, e1, e2, e3, e4, &
                 p0, p5, rx, re, &
                 prf1, prf2, capx, capy, capz
      real    :: prmin, prmax
      logical :: zcubic_L, zqutic_L
!
!---------------------------------------------------------------------
!
      if  (F_lev_S /= 't') call gem_error (-1,'bicubHQV_lin_min_max','NOT AVAILABLE')

      p_bsz_8   => adv_bsz_8%t
      p_zabcd_8 => adv_zabcd_8%t
      p_zbacd_8 => adv_zbacd_8%t
      p_zcabd_8 => adv_zcabd_8%t
      p_zdabc_8 => adv_zdabc_8%t
      p_zbc_8   => adv_zbc_8%t

      p_zxabcde_8 => adv_zxabcde_8%t
      p_zaxbcde_8 => adv_zaxbcde_8%t
      p_zbxacde_8 => adv_zbxacde_8%t
      p_zcxabde_8 => adv_zcxabde_8%t
      p_zdxabce_8 => adv_zdxabce_8%t
      p_zexabcd_8 => adv_zexabcd_8%t

      kkmax = F_nk - 1

      do id=1,F_nind

            m1=id*4
            n0=m1-3
            nx=m1-2
            ny=m1-1
            nz=m1

            n=ii(n0)

            zqutic_L = (ii(nz) > 1) .and. (ii(nz) < kkmax-2)
            zcubic_L =((ii(nz)== 1) .or.   ii(nz)== kkmax-2) .and. .not.Schm_autobar_L

            !- x interpolation
            ra = adv_bsx_8(ii(nx)-1)
            rb = adv_bsx_8(ii(nx)  )
            rc = adv_bsx_8(ii(nx)+1)
            rd = adv_bsx_8(ii(nx)+2)
            p1 = triprd(F_x(n),rb,rc,rd)*adv_xabcd_8
            p2 = triprd(F_x(n),ra,rc,rd)*adv_xbacd_8
            p3 = triprd(F_x(n),ra,rb,rd)*adv_xcabd_8
            p4 = triprd(F_x(n),ra,rb,rc)*adv_xdabc_8

            o2 = (ii(nz)-1)*adv_nijag + (ii(ny)-adv_int_j_off-1)*adv_nit + (ii(nx)-adv_int_i_off)
            o1 = o2-adv_nit
            o3 = o2+adv_nit
            o4 = o3+adv_nit

            if (zqutic_L) then
               o1 = o1 - adv_nijag
               o2 = o2 - adv_nijag
               o3 = o3 - adv_nijag
               o4 = o4 - adv_nijag

               x1 = p1 * F_in(o1-1) + p2 * F_in(o1) + p3 * F_in(o1+1) + p4 * F_in(o1+2)
               x2 = p1 * F_in(o2-1) + p2 * F_in(o2) + p3 * F_in(o2+1) + p4 * F_in(o2+2)
               x3 = p1 * F_in(o3-1) + p2 * F_in(o3) + p3 * F_in(o3+1) + p4 * F_in(o3+2)
               x4 = p1 * F_in(o4-1) + p2 * F_in(o4) + p3 * F_in(o4+1) + p4 * F_in(o4+2)

               o1 = o1 + adv_nijag
               o2 = o2 + adv_nijag
               o3 = o3 + adv_nijag
               o4 = o4 + adv_nijag
            endif

            if (zqutic_L.or.zcubic_L) then
               a1 = p1 * F_in(o1-1) + p2 * F_in(o1) + p3 * F_in(o1+1) + p4 * F_in(o1+2)
               a2 = p1 * F_in(o2-1) + p2 * F_in(o2) + p3 * F_in(o2+1) + p4 * F_in(o2+2)
               a3 = p1 * F_in(o3-1) + p2 * F_in(o3) + p3 * F_in(o3+1) + p4 * F_in(o3+2)
               a4 = p1 * F_in(o4-1) + p2 * F_in(o4) + p3 * F_in(o4+1) + p4 * F_in(o4+2)
            endif

            o1 = o1 + adv_nijag
            o2 = o2 + adv_nijag
            o3 = o3 + adv_nijag
            o4 = o4 + adv_nijag

            prmax = max(F_in(o2),F_in(o2+1),F_in(o3),F_in(o3+1))
            prmin = min(F_in(o2),F_in(o2+1),F_in(o3),F_in(o3+1))

            b1 = p1 * F_in(o1-1) + p2 * F_in(o1) + p3 * F_in(o1+1) + p4 * F_in(o1+2)
            b2 = p1 * F_in(o2-1) + p2 * F_in(o2) + p3 * F_in(o2+1) + p4 * F_in(o2+2)
            b3 = p1 * F_in(o3-1) + p2 * F_in(o3) + p3 * F_in(o3+1) + p4 * F_in(o3+2)
            b4 = p1 * F_in(o4-1) + p2 * F_in(o4) + p3 * F_in(o4+1) + p4 * F_in(o4+2)

            o1 = o1 + adv_nijag
            o2 = o2 + adv_nijag
            o3 = o3 + adv_nijag
            o4 = o4 + adv_nijag

            prmax = max(prmax,F_in(o2),F_in(o2+1),F_in(o3),F_in(o3+1))
            prmin = min(prmin,F_in(o2),F_in(o2+1),F_in(o3),F_in(o3+1))

            c1 = p1 * F_in(o1-1) + p2 * F_in(o1) + p3 * F_in(o1+1) + p4 * F_in(o1+2)
            c2 = p1 * F_in(o2-1) + p2 * F_in(o2) + p3 * F_in(o2+1) + p4 * F_in(o2+2)
            c3 = p1 * F_in(o3-1) + p2 * F_in(o3) + p3 * F_in(o3+1) + p4 * F_in(o3+2)
            c4 = p1 * F_in(o4-1) + p2 * F_in(o4) + p3 * F_in(o4+1) + p4 * F_in(o4+2)

            if (zqutic_L.or.zcubic_L) then
               o1 = o1 + adv_nijag
               o2 = o2 + adv_nijag
               o3 = o3 + adv_nijag
               o4 = o4 + adv_nijag

               d1 = p1 * F_in(o1-1) + p2 * F_in(o1) + p3 * F_in(o1+1) + p4 * F_in(o1+2)
               d2 = p1 * F_in(o2-1) + p2 * F_in(o2) + p3 * F_in(o2+1) + p4 * F_in(o2+2)
               d3 = p1 * F_in(o3-1) + p2 * F_in(o3) + p3 * F_in(o3+1) + p4 * F_in(o3+2)
               d4 = p1 * F_in(o4-1) + p2 * F_in(o4) + p3 * F_in(o4+1) + p4 * F_in(o4+2)
            endif

            if (zqutic_L) then
               o1 = o1 + adv_nijag
               o2 = o2 + adv_nijag
               o3 = o3 + adv_nijag
               o4 = o4 + adv_nijag

               e1 = p1 * F_in(o1-1) + p2 * F_in(o1) + p3 * F_in(o1+1) + p4 * F_in(o1+2)
               e2 = p1 * F_in(o2-1) + p2 * F_in(o2) + p3 * F_in(o2+1) + p4 * F_in(o2+2)
               e3 = p1 * F_in(o3-1) + p2 * F_in(o3) + p3 * F_in(o3+1) + p4 * F_in(o3+2)
               e4 = p1 * F_in(o4-1) + p2 * F_in(o4) + p3 * F_in(o4+1) + p4 * F_in(o4+2)
            endif

            !- y interpolation
            ra = adv_bsy_8(ii(ny)-1)
            rb = adv_bsy_8(ii(ny)  )
            rc = adv_bsy_8(ii(ny)+1)
            rd = adv_bsy_8(ii(ny)+2)
            p1 = triprd(F_y(n),rb,rc,rd)*adv_yabcd_8
            p2 = triprd(F_y(n),ra,rc,rd)*adv_ybacd_8
            p3 = triprd(F_y(n),ra,rb,rd)*adv_ycabd_8
            p4 = triprd(F_y(n),ra,rb,rc)*adv_ydabc_8

            b1 = p1 * b1 + p2 * b2 + p3 * b3 + p4 * b4
            c1 = p1 * c1 + p2 * c2 + p3 * c3 + p4 * c4

            !- z interpolation
            if (zqutic_L) then
               x1 = p1 * x1 + p2 * x2 + p3 * x3 + p4 * x4
               a1 = p1 * a1 + p2 * a2 + p3 * a3 + p4 * a4
               d1 = p1 * d1 + p2 * d2 + p3 * d3 + p4 * d4
               e1 = p1 * e1 + p2 * e2 + p3 * e3 + p4 * e4

               rx = p_bsz_8(ii(nz)-2)
               ra = p_bsz_8(ii(nz)-1)
               rb = p_bsz_8(ii(nz)  )
               rc = p_bsz_8(ii(nz)+1)
               rd = p_bsz_8(ii(nz)+2)
               re = p_bsz_8(ii(nz)+3)

               p0 = quiprd(F_z(n),ra,rb,rc,rd,re)*p_zxabcde_8(ii(nz)+1)
               p1 = quiprd(F_z(n),rx,rb,rc,rd,re)*p_zaxbcde_8(ii(nz)+1)
               p2 = quiprd(F_z(n),rx,ra,rc,rd,re)*p_zbxacde_8(ii(nz)+1)
               p3 = quiprd(F_z(n),rx,ra,rb,rd,re)*p_zcxabde_8(ii(nz)+1)
               p4 = quiprd(F_z(n),rx,ra,rb,rc,re)*p_zdxabce_8(ii(nz)+1)
               p5 = quiprd(F_z(n),rx,ra,rb,rc,rd)*p_zexabcd_8(ii(nz)+1)

               F_qut(n) = p0 * x1 + p1 * a1 + p2 * b1 + p3 * c1 + p4 * d1 + p5 * e1

            elseif (zcubic_L) then
               a1 = p1 * a1 + p2 * a2 + p3 * a3 + p4 * a4
               d1 = p1 * d1 + p2 * d2 + p3 * d3 + p4 * d4
               ra = p_bsz_8(ii(nz)-1)
               rb = p_bsz_8(ii(nz)  )
               rc = p_bsz_8(ii(nz)+1)
               rd = p_bsz_8(ii(nz)+2)
               p1 = triprd(F_z(n),rb,rc,rd)*p_zabcd_8(ii(nz)+1)
               p2 = triprd(F_z(n),ra,rc,rd)*p_zbacd_8(ii(nz)+1)
               p3 = triprd(F_z(n),ra,rb,rd)*p_zcabd_8(ii(nz)+1)
               p4 = triprd(F_z(n),ra,rb,rc)*p_zdabc_8(ii(nz)+1)

               F_qut(n) = p1 * a1 + p2 * b1 + p3 * c1 + p4 * d1
            else
               p3 = (F_z(n)-p_bsz_8(ii(nz)))*p_zbc_8(ii(nz)+1)
               p2 = 1. - p3

               F_qut(n) = p2 * b1 + p3 * c1
            endif

            capx = (F_x(n)-adv_bsx_8(ii(nx))) * adv_xbc_8
            capy = (F_y(n)-adv_bsy_8(ii(ny))) * adv_ybc_8
            capz = (F_z(n)-  p_bsz_8(ii(nz))) * p_zbc_8(ii(nz)+1)

            o2 = (ii(nz))*adv_nijag + (ii(ny)-adv_int_j_off-1)*adv_nit + (ii(nx)-adv_int_i_off)
            o3 = o2 + adv_nit

            prf1 = (1. - capy) * ( (1. - capx) * F_in(o2) + capx * F_in(o2+1) ) + &
                         capy  * ( (1. - capx) * F_in(o3) + capx * F_in(o3+1) )

            o2 = o2 + adv_nijag
            o3 = o3 + adv_nijag

            prf2 = (1. - capy) * ( (1. - capx) * F_in(o2) + capx * F_in(o2+1) ) + &
                         capy  * ( (1. - capx) * F_in(o3) + capx * F_in(o3+1) )

            F_lin (n) = max((1. - capz) * prf1 + capz * prf2,0.0d0)

            F_min (n) = prmin
            F_max (n) = prmax

      end do
!
!---------------------------------------------------------------------
!
      return

contains

      real(kind=REAL64) function triprd(za,zb,zc,zd)
         real, intent(in) :: za
         real(kind=REAL64), intent(in) :: zb,zc,zd
         triprd = (za-zb)*(za-zc)*(za-zd)
      end function triprd

      real(kind=REAL64) function quiprd(zx,za,zb,zc,zd,ze)
         real, intent(in) :: zx
         real(kind=REAL64), intent(in) :: za,zb,zc,zd,ze
         quiprd = (zx-za)*(zx-zb)*(zx-zc)*(zx-zd)*(zx-ze)
      end function quiprd

      end subroutine bicubHQV_lin_min_max
