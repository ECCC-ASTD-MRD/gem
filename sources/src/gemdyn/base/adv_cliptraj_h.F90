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

!**s/r adv_cliptraj_h - Same as adv_cliptraj but with positions with HALO

      subroutine adv_cliptraj_h (F_x,F_y,F_minx,F_maxx,F_miny,F_maxy,F_ni,F_nj,F_nk,i0,in,j0,jn,k0,mesg)

      use adv_grid
      use glb_ld
      use HORgrid_options

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      character(len=*) :: mesg
      integer,intent(in) :: F_minx,F_maxx,F_miny,F_maxy
      integer,intent(in) :: i0,in,j0,jn,k0 !I, scope of the operator
      integer,intent(in) :: F_ni,F_nj,F_nk
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk) :: F_x, F_y !I/O, upstream pos

      !object
      !=======================================================
      !     Clip SL hor. trajectories to either fit inside the
      !     physical domain of the processor or to the
      !     actual maximum allowed COURANT number (LAM)
      !=======================================================

#include "stop_mpi.h"
#include "msg.h"

      !------------------------------------------------------------------------------------------------
      real*8,  parameter :: EPS_8 = 1.D-5
      integer :: BCS_BASE       ! BCS points for Yin-Yang, normal LAM

      character(len=MSG_MAXLEN) :: msg_S
      integer :: i,j,k,cnt,sum_cnt,cto,sum_cto,err
      real :: minposx,maxposx,minposy,maxposy

      !------------------------------------------------------------------------------------------------

      call timing_start2 (35, 'ADV_CLIP', 34)

      BCS_BASE= 4
      if (Grd_yinyang_L) BCS_BASE = 3
      minposx = adv_xx_8(adv_lminx+1) + EPS_8
      if (l_west)  minposx = adv_xx_8(1+BCS_BASE) + EPS_8
      maxposx = adv_xx_8(adv_lmaxx-1) - EPS_8
      if (l_east)  maxposx = adv_xx_8(F_ni-BCS_BASE) - EPS_8
      minposy = adv_yy_8(adv_lminy+1) + EPS_8
      if (l_south) minposy = adv_yy_8(1+BCS_BASE) + EPS_8
      maxposy = adv_yy_8(adv_lmaxy-1) - EPS_8
      if (l_north) maxposy = adv_yy_8(F_nj-BCS_BASE) - EPS_8

      cnt=0
      cto=0

      !Clipping to processor boundary
      !------------------------------
      do k=k0,F_nk
         do j=j0,jn
            do i=i0,in

            cto=cto+1

            if ((F_x(i,j,k)<minposx).or.(F_x(i,j,k)>maxposx).or. &
                (F_y(i,j,k)<minposy).or.(F_y(i,j,k)>maxposy) ) then
               cnt=cnt+1
            end if

            F_x(i,j,k) = min(max(F_x(i,j,k),minposx),maxposx)

            F_y(i,j,k) = min(max(F_y(i,j,k),minposy),maxposy)

         end do
      end do
      end do

      call rpn_comm_Allreduce(cnt,sum_cnt,1,"MPI_INTEGER", "MPI_SUM","grid",err)
      call rpn_comm_Allreduce(cto,sum_cto,1,"MPI_INTEGER", "MPI_SUM","grid",err)

      if (trim(mesg) /= "" .and. sum_cnt>0) then
         write(msg_S,'(a,i6,a,f6.2,2x,a)')  &
         ' ADW trajtrunc: npts=',sum_cnt, &
         ', %=',real(sum_cnt)/real(sum_cto)*100., &
         mesg
         call msg(MSG_INFO,msg_S)
      end if

      call timing_stop (35)

      return
      end subroutine adv_cliptraj_h

