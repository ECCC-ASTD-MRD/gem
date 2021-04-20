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

!**s/r yyg_xchng_sca_q2q - Interpolate and exchange scalars fields
!                          from q to q points

      subroutine yyg_xchng_sca_q2q ( F_src, F_comm, Minx,Maxx,Miny,Maxy,&
                                     F_ni, F_nj, Nk, F_interpo_S       ,&
                                     mono_L, do_xch )
      use ISO_C_BINDING
      use gem_timing
      use yyg_param
      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in)  :: F_interpo_S
      logical, intent(in) :: mono_L, do_xch
      integer, intent(in) :: Minx,Maxx,Miny,Maxy, Nk, F_ni,F_nj
      real, intent(inout) :: F_src (Minx:Maxx,Miny:Maxy,Nk)
      type(YYG_comm_param), intent(inout) :: F_comm

      integer ii,jj,k,kk,m,mm,adr
!
!----------------------------------------------------------------------
!
      call gemtime_start ( 6, 'YYG_XCHNG', 0)

      call yyg_SendRecv_s ( F_src, F_comm, Minx,Maxx,Miny,Maxy, &
                            NK, F_interpo_S, mono_L )

! Fill my results buffers if I have received something
      if (F_comm%maxrecv > 0) then

         do kk=1, F_comm%recvmaxproc
            mm=0
            do m= 1,F_comm%recv_len(kk)
               adr= F_comm%recv_adr(kk)+m
               ii = F_comm%recv_i(adr)
               jj = F_comm%recv_j(adr)
               do k=1,Nk
                  mm=mm+1
                  F_src(ii,jj,k) = YYG_uvrecv(mm,KK)
               end do
            end do
         end do

      end if

      if (do_xch) then
          if (Glb_pilotcirc_L) then
              call rpn_comm_propagate_pilot_circular ( F_src,&
                       l_minx,l_maxx,l_miny,l_maxy, F_ni,F_nj,nk,&
                       Glb_pil_e,Glb_pil_s,G_halox,G_haloy )
          else
              call rpn_comm_xch_halo(F_src,l_minx,l_maxx,l_miny,l_maxy,&
               l_ni,l_nj,Nk,G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
          endif
      endif

      call gemtime_stop (6)
!
!----------------------------------------------------------------------
!
      return
      end subroutine yyg_xchng_sca_q2q

