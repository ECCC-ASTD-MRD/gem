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
!*s/r yyg_rhs_xchng - Interpolate and exchange RHS data

      subroutine yyg_rhs_xchng ( F_rhs_8, F_sol_8, minx,maxx,miny,maxy,&
                                 lminx,lmaxx,lminy,lmaxy              ,&
                                 NK,iter,Linfini )
      use geomh
      use glb_ld
      use sol
      use ptopo
      use yyg_rhs
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: minx,maxx,miny,maxy,lminx,lmaxx,lminy,lmaxy, NK,iter
      real(kind=REAL64) , dimension(lminx:lmaxx,lminy:lmaxy,Nk), intent(in) :: F_sol_8
      real(kind=REAL64) , dimension( minx:maxx , miny:maxy ,Nk), intent(inout) :: F_rhs_8
      real, intent(out) :: Linfini

      integer status,ierr,i,k,kk,kk_proc,m,mm,adr
      integer tag1,tag2,recvlen,sendlen,ireq
      integer request(Ptopo_numproc*2)
      real L1,L2,L3(2),L4(2),infini(Nk)
      real  , dimension (:,:), allocatable :: recv_pil,send_pil
      real(kind=REAL64), dimension (:  ), allocatable :: send_Rhsx_8
!
!----------------------------------------------------------------------
!
      tag2= 14 ; tag1= 13 ; infini= 0. ; linfini= 0.0

      sendlen=0 ; recvlen=0
      ireq=0
      do kk=1,Rhsx_sendmaxproc
         sendlen=max(sendlen,Rhsx_send_len(kk))
      end do
      do kk=1,Rhsx_recvmaxproc
         recvlen=max(recvlen,Rhsx_recv_len(kk))
      end do

      if (sendlen > 0) then
          allocate(send_pil(sendlen*NK,Rhsx_sendmaxproc))
          allocate(send_Rhsx_8(sendlen*NK))
      end if
      if (recvlen > 0) then
          allocate(recv_pil(recvlen*NK,Rhsx_recvmaxproc))
      end if

!     For each processor (in other colour)
      do kk= 1, Rhsx_recvmaxproc
         if (Ptopo_couleur == 0) then
             kk_proc = Rhsx_recvproc(kk)+Ptopo_numproc-1
         else
             kk_proc = Rhsx_recvproc(kk)-1
         end if
         if (Rhsx_recv_len(kk) > 0) then
             ireq = ireq+1
             call RPN_COMM_IRecv(recv_pil(1,KK),Rhsx_recv_len(kk)*NK,'MPI_REAL',&
                    kk_proc,tag2+kk_proc,'MULTIGRID',request(ireq),ierr)
         end if
      end do

!     For each processor (in other colour)
      do kk= 1, Rhsx_sendmaxproc
         if (Ptopo_couleur == 0) then
             kk_proc = Rhsx_sendproc(kk)+Ptopo_numproc-1
         else
             kk_proc = Rhsx_sendproc(kk)-1
         end if

         if (Rhsx_send_len(kk) > 0) then
             adr=Rhsx_send_adr(kk)+1
             call yyg_int_cub88 ( send_Rhsx_8,F_sol_8                ,&
                             Rhsx_send_imx(adr:adr+Rhsx_send_len(kk)-1),&
                             Rhsx_send_imy(adr:adr+Rhsx_send_len(kk)-1),&
                             geomh_x_8,geomh_y_8,l_minx              ,&
                             l_maxx,l_miny,l_maxy, NK                ,&
                             Rhsx_send_xxr(adr:adr+Rhsx_send_len(kk)-1),&
                             Rhsx_send_yyr(adr:adr+Rhsx_send_len(kk)-1),&
                             Rhsx_send_len(kk) )
             mm=0
             do m=1,Rhsx_send_len(kk)
                do k=1,NK
                   mm=mm+1
                   infini(k)=max(infini(k),abs(real(send_Rhsx_8(mm))-Sol_rhs(mm,1,kk)))
                   Sol_rhs(mm,1,kk)=real(send_Rhsx_8(mm))
                   send_pil(mm,KK)=real(send_Rhsx_8(mm)*Rhsx_send_sten(adr+m-1))
                end do
             end do

             linfini=infini(1)
             do k=2,NK
                linfini=max(linfini,infini(k))
             end do

             ireq = ireq+1
             call RPN_COMM_ISend (send_pil(1,KK),Rhsx_send_len(kk)*NK,&
                                  'MPI_REAL', kk_proc,tag2+Ptopo_world_myproc, &
                                  'MULTIGRID',request(ireq),ierr)
         end if

      end do

!Wait for all done sending and receiving
      call RPN_COMM_waitall_nostat(ireq,request,ierr)

      if (sendlen > 0) deallocate(send_pil,send_Rhsx_8)

      if (iter > 0) then
         l2=0
         do kk= 1, Rhsx_sendmaxproc
         do  i= 1, sendlen*NK
            l2 =max(l2,abs(Sol_rhs(i,1,kk)))
         end do
         end do
         call RPN_COMM_allreduce(L2,L1,1,"MPI_REAL","MPI_MAX","grid",ierr)
         L3(2)= L1 ; L2= L1

!     Check the precision
         call RPN_COMM_allreduce(Linfini,L1,1,"MPI_REAL","MPI_MAX","grid",ierr)
         L3(1)= L1 ; L4= 0.
         Linfini=L1

         if (Ptopo_myproc == 0.and.Ptopo_couleur == 0) then
            call RPN_COMM_Send(L3,2,'MPI_REAL',1,tag1,'GRIDPEERS',ierr)
            call RPN_COMM_Recv(L4,2,'MPI_REAL',1,tag2,'GRIDPEERS',status,ierr)
         end if
         if (Ptopo_myproc == 0.and.Ptopo_couleur == 1) then
            call RPN_COMM_Recv(L4,2,'MPI_REAL',0,tag1,'GRIDPEERS',status,ierr)
            call RPN_COMM_Send(L3,2,'MPI_REAL',0,tag2,'GRIDPEERS',ierr)
         end if
         Linfini=max(Linfini,L4(1))
         L2=max(L2,L4(2))
         if (iter > 1) linfini=linfini/L2
         call RPN_COMM_bcast(linfini,1,'MPI_REAL',0,'GRID',ierr)
      end if

! Now fill my results if I have received something

      if (recvlen > 0) then

         do  kk= 1, Rhsx_recvmaxproc
            mm=0
            do m= 1, Rhsx_recv_len(kk)
               adr=Rhsx_recv_adr(kk)+m
               do k=1,NK
                  mm=mm+1
                  F_rhs_8(Rhsx_recv_i(adr),Rhsx_recv_j(adr),k)= &
                  F_rhs_8(Rhsx_recv_i(adr),Rhsx_recv_j(adr),k)+recv_pil(mm,KK)
               end do
            end do
         end do

         deallocate(recv_pil)

      end if
!
!----------------------------------------------------------------------
!
      return
      end subroutine yyg_rhs_xchng

