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

!**s/r yyg_blend_sca - Blend a scalar variable on the overlap region

      subroutine yyg_blend_sca ( F_src, F_comm, Minx,Maxx,Miny,Maxy,Nk )
      use ISO_C_BINDING
      use yyg_param
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Minx,Maxx,Miny,Maxy, Nk
      real, intent(inout) :: F_src (Minx:Maxx,Miny:Maxy,Nk)
      type(YYG_comm_param), intent(inout) :: F_comm

      integer ii,jj,k,kk,m,mm,adr
!
!----------------------------------------------------------------------
!
      call yyg_SendRecv_s ( F_src, F_comm, Minx,Maxx,Miny,Maxy, &
                            NK, 'CUBIC', .false. )

! Fill my results buffers if I have received something
      if (F_comm%maxrecv > 0) then

         do kk=1, F_comm%recvmaxproc
            mm=0
            do m=1 ,F_comm%recv_len(kk)
               adr= F_comm%recv_adr(kk)+m
               ii = F_comm%recv_i(adr)
               jj = F_comm%recv_j(adr)
               do k=1,Nk
                  mm=mm+1
                  F_src(ii,jj,k)= F_src(ii,jj,k) * 0.5 +&
                               YYG_uvrecv(mm,KK) * 0.5
               end do
            end do
         end do

      end if
!
!----------------------------------------------------------------------
!
      return
      end subroutine yyg_blend_sca

