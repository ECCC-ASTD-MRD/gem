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

!**s/r yyg_xchng_vec_uv2uv - Interpolate and exchange vectors
!                            fields from u,v to u,v points

      subroutine yyg_xchng_vec_uv2uv ( F_u, F_v, Minx,Maxx,Miny,Maxy,NK )
      use ISO_C_BINDING
      use gem_timing
      use yyg_param
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Minx,Maxx,Miny,Maxy,NK
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(inout) :: F_u, F_v

      integer ii,jj,k,kk,m,mm,adr
!
!----------------------------------------------------------------------
!
      call gemtime_start ( 6, 'YYG_XCHNG', 0)

      call yyg_SendRecv_v ( F_u, F_v, YYG_PILT_uv2u, YYG_PILT_uv2v, &
                            Minx,Maxx,Miny,Maxy,NK )

! Fill my results buffers if I have received something
      if (YYG_PILT_recv_uv > 0) then

         do kk=1, YYG_PILT_uv2u%recvmaxproc
            mm=0
            do m=1, YYG_PILT_uv2u%recv_len(kk)
               adr= YYG_PILT_uv2u%recv_adr(kk)+m
               ii = YYG_PILT_uv2u%recv_i(adr)
               jj = YYG_PILT_uv2u%recv_j(adr)
               do k=1,Nk
                  mm=mm+1
                  F_u(ii,jj,k)= YYG_urecvuv(mm,KK)
               end do
            end do
         end do

         do kk=1, YYG_PILT_uv2v%recvmaxproc
            mm=0
            do m= 1,YYG_PILT_uv2v%recv_len(kk)
               adr= YYG_PILT_uv2v%recv_adr(kk)+m
               ii = YYG_PILT_uv2v%recv_i(adr)
               jj = YYG_PILT_uv2v%recv_j(adr)
               do k=1,Nk
                  mm=mm+1
                  F_v(ii,jj,k)= YYG_vrecvuv(mm,KK)
               end do
            end do
         end do

      end if

      call gemtime_stop (6)
!
!----------------------------------------------------------------------
!
      return
      end subroutine yyg_xchng_vec_uv2uv

