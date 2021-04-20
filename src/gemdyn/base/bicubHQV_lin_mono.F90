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

!**s/r bicubHQV_lin_min_max - SL Bi-Cubic H interpolation and Quintic V interpolation except
!                             Cubic or Linear V interpolation near vertical boundaries
!                             (Based on tricublin_beta.c): Positions in Index space

      subroutine bicubHQV_lin_min_max (F_qut, F_lin, F_min, F_max, F_in, F_xyz, F_num, &
                                       F_minx, F_maxx, F_miny, F_maxy)

      use adz_mem
      use dynkernel_options
      use ver

      use, intrinsic :: iso_fortran_env

      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_num,F_minx,F_maxx,F_miny,F_maxy
      real,dimension(*), intent(in)  :: F_in
      real,dimension(*), intent(out) :: F_qut,F_lin,F_min,F_max
      real,dimension(3,F_num), intent(in) :: F_xyz

      logical :: zcubic_L,zqutic_L
      integer :: n,o1,indxij,ii1,jj1,kk1,nit,nij
      real(kind=REAL64), parameter :: ov6 = 1.d0/6.d0, ov2 = 1.d0/2.d0
      real(kind=REAL64), dimension(:), pointer ::&
                         p_zabcd_8, p_zbacd_8,& 
                         p_zcabd_8, p_zdabc_8, posz, dz, odz,&
                         p_zxabcde_8, p_zaxbcde_8, p_zbxacde_8,&
                         p_zcxabde_8, p_zdxabce_8, p_zexabcd_8
      real(kind=REAL64) :: px(4),py(4),qz(6),del,del12,del10,&
                           ra,rb,rc,rd,rx,re,lx(2),ly(2),lz(2),zt
!
!---------------------------------------------------------------------
!
      posz => Ver_z_8%t
      dz   => Adz_delz_t
      odz  => Adz_odelz_t

      p_zabcd_8 => Adz_zabcd_8%t
      p_zbacd_8 => Adz_zbacd_8%t
      p_zcabd_8 => Adz_zcabd_8%t
      p_zdabc_8 => Adz_zdabc_8%t

      p_zxabcde_8 => Adz_zxabcde_8%t
      p_zaxbcde_8 => Adz_zaxbcde_8%t
      p_zbxacde_8 => Adz_zbxacde_8%t
      p_zcxabde_8 => Adz_zcxabde_8%t
      p_zdxabce_8 => Adz_zdxabce_8%t
      p_zexabcd_8 => Adz_zexabcd_8%t

      nit = (F_maxx - F_minx + 1)
      nij = (F_maxx - F_minx + 1)*(F_maxy - F_miny + 1)

      do n=1,F_num

         ii1   = F_xyz(1,n)
         del   = F_xyz(1,n) - dble(ii1)
         del12 = (del-1)*(del-2)
         del10 = (del+1)* del
         px(1) = -del   *del12 *ov6
         px(2) = (del+1)*del12 *ov2
         px(3) = -del10*(del-2)*ov2
         px(4) =  del10*(del-1)*ov6

         lx(1) = 1.0 - del ; lx(2) = del 

         jj1   = F_xyz(2,n)
         del   = F_xyz(2,n) - dble(jj1)
         del12 = (del-1)*(del-2)
         del10 = (del+1)* del
         py(1) = -del   *del12 *ov6
         py(2) = (del+1)*del12 *ov2
         py(3) = -del10*(del-2)*ov2
         py(4) =  del10*(del-1)*ov6

         ly(1) = 1.0 - del ; ly(2) = del 

         ii1 = ii1 - F_minx - 1
         jj1 = jj1 - F_miny - 1
         indxij = (jj1-Adz_joff)*nit + (ii1-Adz_ioff) - nit

         kk1 = F_xyz(3,n)

         zqutic_L = (kk1 >= 3) .and.  (kk1 <= l_nk-3)  .and. .not.Schm_autobar_L
         zcubic_L =((kk1 == 2) .or.   (kk1 == l_nk-2)) .and. .not.Schm_autobar_L

         if (.not.zqutic_L.and..not.zcubic_L) kk1 = min(l_nk-1,max(1,kk1)) 

         zt = posz(kk1) + (F_xyz(3,n)-dble(kk1))*dz(kk1)

         lz(2) = (zt-posz(kk1))*odz(kk1); lz(1) = 1.0 - lz(2)

         if (zqutic_L) then

            o1 = (kk1-3)*nij + indxij

            rx = posz(kk1-2)
            ra = posz(kk1-1)
            rb = posz(kk1  )
            rc = posz(kk1+1)
            rd = posz(kk1+2)
            re = posz(kk1+3)

            qz(1) = quiprd(zt,ra,rb,rc,rd,re)*p_zxabcde_8(kk1)
            qz(2) = quiprd(zt,rx,rb,rc,rd,re)*p_zaxbcde_8(kk1)
            qz(3) = quiprd(zt,rx,ra,rc,rd,re)*p_zbxacde_8(kk1)
            qz(4) = quiprd(zt,rx,ra,rb,rd,re)*p_zcxabde_8(kk1)
            qz(5) = quiprd(zt,rx,ra,rb,rc,re)*p_zdxabce_8(kk1)
            qz(6) = quiprd(zt,rx,ra,rb,rc,rd)*p_zexabcd_8(kk1)

         elseif (zcubic_L) then

            o1 = (kk1-2)*nij + indxij

            ra = posz(kk1-1)
            rb = posz(kk1  )
            rc = posz(kk1+1)
            rd = posz(kk1+2)

            qz(1) = 0.
            qz(2) = triprd(zt,rb,rc,rd)*p_zabcd_8(kk1)
            qz(3) = triprd(zt,ra,rc,rd)*p_zbacd_8(kk1)
            qz(4) = triprd(zt,ra,rb,rd)*p_zcabd_8(kk1)
            qz(5) = triprd(zt,ra,rb,rc)*p_zdabc_8(kk1)
            qz(6) = 0.

         else

            o1 = (kk1-1)*nij + indxij

            qz(1) = 0.
            qz(2) = 0.
            qz(4) = (zt-posz(kk1))*odz(kk1)
            qz(3) = 1.d0 - qz(4)
            qz(5) = 0. 
            qz(6) = 0.

         end if
 
         call bicubHQV_lin_min_max_zyx (F_qut(n),F_lin(n),F_min(n),F_max(n),F_in(o1),px,py,qz,lx,ly,lz,&
                                        nit,nij,zqutic_L,zcubic_L)

      end do
!
!---------------------------------------------------------------------
!
      return

contains

      real(kind=REAL64) function triprd(za,zb,zc,zd)
         real(kind=REAL64), intent(in) :: za,zb,zc,zd
         triprd = ((za-zb)*(za-zc)*(za-zd))
      end function triprd

      real(kind=REAL64) function quiprd(zx,za,zb,zc,zd,ze)
         real(kind=REAL64), intent(in) :: zx,za,zb,zc,zd,ze
         quiprd = (zx-za)*(zx-zb)*(zx-zc)*(zx-zd)*(zx-ze)
      end function quiprd

      end subroutine bicubHQV_lin_min_max 
