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

!**s/r vspng_abc -- Prepares matrices aix,bix,cix,dix,aiy,biy,ciy
!
      subroutine vspng_abc2(F_aix_8 , F_bix_8, F_cix_8, F_dix_8 , &
                            F_aiy_8 , F_biy_8, F_ciy_8, F_cy2_8 , &
                            F_xp0_8 , F_xp2_8, F_yp0_8, F_yp2_8 , &
                            F_coef_8, F_njpole, Gni, Gnj,F_gnjv, NK)
      use glb_ld
      use trp
      use ptopo
      implicit none
#include <arch_specific.hf>

      integer F_njpole,Gni,Gnj,F_gnjv,NK
      real*8 F_aix_8(*), F_bix_8(*), F_cix_8(*), F_dix_8(*), &
             F_aiy_8(*), F_biy_8(*), F_ciy_8(*), F_cy2_8(*), &
             F_xp0_8(Gni,3), F_xp2_8(Gni,3), &
             F_yp0_8(Gnj,3), F_yp2_8(Gnj,3), &
             F_coef_8(NK)

      integer i, j, k, jj, j2, in, jn
      real*8 ax_8(NK,trp_12emax,G_ni), bx_8(NK,trp_12emax,G_ni), &
             cx_8(NK,trp_12emax,G_ni), ay_8(NK,trp_22emax,G_nj), &
             by_8(NK,trp_22emax,G_nj), cy_8(NK,trp_22emax,G_nj), &
             diy_8,mdifc(Gnj)
      real*8 ZERO_8,ONE_8,HALF_8
      parameter ( ZERO_8 = 0.0 , ONE_8 = 1.0 , HALF_8 = 0.5 )
!
!     ---------------------------------------------------------------
!
      if (F_njpole < 1) then
         mdifc = ONE_8
      else
         do j = 1, F_njpole
            mdifc(j) = dble(j-1)/dble(F_njpole)
         end do
         do j = F_njpole+1, F_gnjv-F_njpole
            mdifc(j) = ONE_8
         end do
         do j = F_gnjv-F_njpole+1, Gnj
            mdifc(j) = max(ZERO_8,dble(F_gnjv-j)/dble(F_njpole))
         end do
      end if
!
!     Calcul le long de X
!     calculate the ending point JN of where to fill the data
!     as the tile is ldnh_maxy size (l_maxy size)
!     jn = trp_12en
!88   j2 = Ptopo_gindx(3,Ptopo_myproc+1) + Trp_12en0 + jn - 2
!     if (j2 > Gnj) then
!        jn = jn - 1
!        goto 88
!     end if
      jn = trp_12en
      j2 = Ptopo_gindx(3,Ptopo_myproc+1) + Trp_12en0 + jn - 2
      if (j2 > Gnj) jn = jn - (j2-Gnj)

!     Insure that any filling on the end of the tile is within the tile
!     in case JN is negative
      jn = max(0,jn)

!$omp parallel private(jj,j2,j,k)
!$omp do
      do i = 1, G_ni
         do j = 1, jn
            jj = Trp_12en0 + j - 1
            j2 = Ptopo_gindx(3,Ptopo_myproc+1) + jj - 1
            do k = 1, NK
            ax_8(k,j,i) = F_xp0_8(i,1) - F_xp2_8(i,1) * mdifc(j2) &
                         *F_coef_8(k) / cos(F_cy2_8(jj))**2
            bx_8(k,j,i) = F_xp0_8(i,2) - F_xp2_8(i,2) * mdifc(j2) &
                         *F_coef_8(k) / cos(F_cy2_8(jj))**2
            cx_8(k,j,i) = F_xp0_8(i,3) - F_xp2_8(i,3) * mdifc(j2) &
                         *F_coef_8(k) / cos(F_cy2_8(jj))**2
            end do
         end do
!
         do j = jn+1,trp_12emax
            do k = 1, NK
               bx_8(k,j,i)=  ONE_8
               cx_8(k,j,i)= ZERO_8
               ax_8(k,j,i)= ZERO_8
            end do
         end do
      end do
!$omp enddo
!$omp end parallel
!
      call set_trig21 (F_aix_8,F_bix_8,F_cix_8,F_dix_8, ax_8,bx_8,cx_8,  &
                       NK*trp_12emax, 1, G_ni,  &
                       NK*trp_12en, .true.)
!
!     Calcul le long de Y
!
!     calculate the ending point IN of where to fill the data
!     as the tile is ldnh_maxx size (l_maxx size)
!     in = trp_22en
!99   j2 = Ptopo_gindx(1,Ptopo_myproc+1) + trp_22en0 + in - 2
!     if (j2 > G_ni) then
!        in = in - 1
!        goto 99
!     end if
      in = trp_22en
      j2 = Ptopo_gindx(1,Ptopo_myproc+1) + trp_22en0 + in - 2
      if (j2 > G_ni) in = in - (j2-G_ni)

!     Insure that any filling on the end of the tile is within the tile
!     in case IN is negative
      in = max(0,in)

!$omp parallel private(i,k)
!$omp do
      do j= 1, F_gnjv
         do i= 1, in
         do k= 1, NK
            ay_8(k,i,j) =  F_yp0_8(j,1) - F_yp2_8(j,1)  &
                          *F_coef_8(k) * mdifc(j)
            by_8(k,i,j) =  F_yp0_8(j,2) - F_yp2_8(j,2)  &
                          *F_coef_8(k) * mdifc(j)
            cy_8(k,i,j) =  F_yp0_8(j,3) - F_yp2_8(j,3)  &
                          *F_coef_8(k) * mdifc(j)
         end do
         end do
      end do
!$omp enddo
!$omp do
      do i= 1, in
         do k= 1, NK
            ay_8(k,i,1) = ZERO_8
         end do
      end do
!$omp enddo

!$omp do
      do j = 1, F_gnjv
      do i = in+1,trp_22emax
      do k = 1, NK
         by_8(k,i,j)=  ONE_8
         cy_8(k,i,j)= ZERO_8
         ay_8(k,i,j)= ZERO_8
      end do
      end do
      end do
!$omp enddo
!$omp end parallel

      call set_trig21 (F_aiy_8,F_biy_8,F_ciy_8,diy_8, ay_8,by_8,cy_8,  &
                       NK*trp_22emax, 1, F_gnjv,  &
                       NK*trp_22en, .false.)
!
!     ---------------------------------------------------------------
!
      return
      end
