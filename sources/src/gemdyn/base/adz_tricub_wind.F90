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

!**s/r adz_tricub_wind - SL tricubic interpolation of 3D winds

      subroutine adz_tricub_wind (F_u, F_v, F_w, F_uvw_in, F_xyz, F_num)
      use ISO_C_BINDING
      use adz_mem
      use ver
      implicit none
#include <arch_specific.hf>

      include "tricublin_beta.inc"

      integer, intent(in) :: F_num
      real,dimension(3,F_num), intent(in)  :: F_xyz
      real,dimension(3,*  ), intent(in)  :: F_uvw_in
      real,dimension(F_num), intent(out) :: F_u, F_v, F_w

      logical(C_BOOL) :: zcubic_L
      integer n, o1, indxij
      integer :: ii1, jj1, kk1, kp1
      real*8, parameter :: ov6 = 1.d0/6.d0, ov2 = 1.d0/2.d0
      real*8  :: px(4),py(4),pz(4), del,del12,del10
      real*8  :: triprd, zb, zc, zd, ra, rb, rc, rd
      real    :: xyz(3), za, zt
      triprd(za,zb,zc,zd)=(za-zb)*(za-zc)*(za-zd)
!
!--------------------------------------------------------
!
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
               zt = Ver_z_8%m(kk1) + (F_xyz(3,n)-dble(kk1))*Adz_delz_m(kk1)
               kk1 = min(Adz_kkmax,max(0,kk1 - 1))

               zcubic_L = ( kk1 > 0) .and. (kk1 < Adz_kkmax)

               kp1 = kk1 + 1
               if (zcubic_L) then
                  ra    = ver_z_8%m(kk1  )
                  rb    = ver_z_8%m(kp1  )
                  rc    = ver_z_8%m(kk1+2)
                  rd    = ver_z_8%m(kk1+3)
                  pz(1) = triprd(zt,rb,rc,rd)*Adz_zabcd_8%m(kp1)
                  pz(2) = triprd(zt,ra,rc,rd)*Adz_zbacd_8%m(kp1)
                  pz(3) = triprd(zt,ra,rb,rd)*Adz_zcabd_8%m(kp1)
                  pz(4) = triprd(zt,ra,rb,rc)*Adz_zdabc_8%m(kp1)
               else
                  kk1   = min(l_nk-1,kk1+1)
                  pz(1) = 0.d0 ; pz(4) = 0.d0
                  pz(3) = (zt-ver_z_8%m(kk1))*Adz_odelz_m(kk1)
                  pz(2) = 1.d0 - pz(3)
               endif

               o1 = (kk1-1)*Adz_nij + indxij
               call tricublin_zyx3f_beta (xyz,F_uvw_in(1,o1),&
                                          px,py,pz,Adz_nit,Adz_nij,zcubic_L)
               F_u(n) = xyz(1)
               F_v(n) = xyz(2)
               F_w(n) = xyz(3)

      enddo
!
!--------------------------------------------------------
!
      end subroutine adz_tricub_wind
