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

!**s/r yyg_xchng_vec_q2q - Interpolate and exchange u,v vectors
!                          fields from q to q points

      subroutine yyg_xchng_vec_q2q ( F_u, F_v, Minx,Maxx,Miny,Maxy, NK )
      use ISO_C_BINDING
      use gem_options
      use glb_ld
      use ptopo
      use yyg_param
      implicit none
#include <arch_specific.hf>

      include "intrp_bicub_yx.inc"

      integer, intent(in) :: Minx,Maxx,Miny,Maxy,NK
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(inout) :: F_u, F_v

      integer k,kk,kk_proc,adr,m,mm,ierr
      integer tag1,ireq,ind1,ind2
      integer request(Ptopo_numproc*2)
!
!----------------------------------------------------------------------
!
      call rpn_comm_xch_halo(F_u, Minx,Maxx,Miny,Maxy,l_ni,l_nj,Nk, &
                             G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
      call rpn_comm_xch_halo(F_v, Minx,Maxx,Miny,Maxy,l_ni,l_nj,Nk, &
                             G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

      tag1=14 ; ireq=0

! Posting the receive requests
      do kk= 1, YYG_PILT_q2q%recvmaxproc

         kk_proc= yyg_proc(YYG_PILT_q2q%recvproc(kk))
         m= YYG_PILT_q2q%recv_len(kk)
         if (m > 0) then
            ireq = ireq+1
            call RPN_COMM_IRecv (YYG_uvrecv(1,KK), m*NK*2        ,&
                                 'MPI_REAL', kk_proc,tag1+kk_proc,&
                                 'MULTIGRID',request(ireq),ierr)
         end if

      end do

! Interpolation and shipping
      do kk= 1, YYG_PILT_q2q%sendmaxproc

         kk_proc= yyg_proc(YYG_PILT_q2q%sendproc(kk))
         mm= YYG_PILT_q2q%send_len(kk)
         if (mm > 0) then

             adr=YYG_PILT_q2q%send_adr(kk)

             do m=1,YYG_PILT_q2q%send_len(KK)
                ind1= (m-1)*2*nk+1
                ind2= adr+m
                call intrp_bicub_yx_uv ( F_u(1,1,1), F_v(1,1,1) ,&
                      YYG_uvsend(ind1,kk),YYG_uvsend(ind1+nk,kk),&
                  (Maxx-Minx+1), (Maxx-Minx+1)*(Maxy-Miny+1), nk,&
                                      YYG_PILT_q2q%send_xu(ind2),&
                                      YYG_PILT_q2q%send_yu(ind2),&
                               YYG_PILT_q2q%send_s((ind2-1)*4+1) )
             end do

             ireq = ireq+1
             call RPN_COMM_ISend ( YYG_uvsend(1,KK),mm*NK*2,'MPI_REAL',&
                                   kk_proc,tag1+Ptopo_world_myproc,&
                                   'MULTIGRID', request(ireq),ierr )
         end if

      end do

      call RPN_COMM_waitall_nostat(ireq,request,ierr)

! Fill my results buffers if I have received something
      if (YYG_PILT_q2q%maxrecv > 0) then

         do kk=1, YYG_PILT_q2q%recvmaxproc
            do m=1,YYG_PILT_q2q%recv_len(kk)
               adr=YYG_PILT_q2q%recv_adr(kk)+m
               do k=1,Nk
                  ind1= (m-1)*2*nk+k
                  ind2= ind1+nk
                  F_u(YYG_PILT_q2q%recv_i(adr),&
                      YYG_PILT_q2q%recv_j(adr),k) = YYG_uvrecv(ind1,KK)
                  F_v(YYG_PILT_q2q%recv_i(adr),&
                      YYG_PILT_q2q%recv_j(adr),k) = YYG_uvrecv(ind2,KK)
               end do
            end do
         end do

      end if
!
!----------------------------------------------------------------------
!
      return
      end subroutine yyg_xchng_vec_q2q

