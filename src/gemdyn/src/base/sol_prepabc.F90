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

!**s/r sol_prepabc - initialisation de l'operateur unidimensionnel.
!                    direction de travail: suivant la latitude.
!
      subroutine sol_prepabc (F_ai_8, F_bi_8, F_ci_8, F_hcon_8, NK)
      use glb_ld
      use glb_pil
      use numa
      use opr
      use sol
      use trp
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
   
      integer NK
      real(kind=REAL64), intent(in) :: F_hcon_8(NK)
      real(kind=REAL64), dimension(Sol_mink:Sol_maxk,G_ni,Sol_miny:Sol_maxy  ), intent(out) :: F_ai_8, F_bi_8
      real(kind=REAL64), dimension(Sol_mink:Sol_maxk,G_ni,Sol_miny-1:Sol_maxy), intent(out) :: F_ci_8

      integer i, j, k, j0, jn, dim, err, tag
      integer status(MPI_STATUS_SIZE)
      real(kind=REAL128) di_8
      real(kind=REAL64) abc(3,Sol_mink:Sol_maxk,G_ni,Sol_miny:Sol_maxy)
!     
!     ---------------------------------------------------------------
!
      dim= Sol_sock_nk*G_ni ; tag=95
      if (Numa_sockrank==0) then
         j0= max(Sol_miny,2   +Lam_pil_s)
         abc(1,:,:,:)= 0.d0 ; abc(2,:,:,:)= 1.d0 ; abc(3,:,:,:)= 0.d0
         do j = j0, min(Sol_maxy,G_nj-Lam_pil_n)
            di_8 = Opr_opsyp0_8(G_nj+j) / cos( G_yg_8(j) )**2
            do i = 1+Lam_pil_w, G_ni-Lam_pil_e
               do k = Sol_mink,Sol_maxk
                  abc(2,k,i,j)= Opr_xeval_8(i) * di_8 + &
                                F_hcon_8(k)*Opr_opsyp0_8(  G_nj+j  )+Opr_opsyp2_8(  G_nj+j  )
                  abc(3,k,i,j)= F_hcon_8(k)*Opr_opsyp0_8(2*G_nj+j  )+Opr_opsyp2_8(2*G_nj+j  )
                  abc(1,k,i,j)= F_hcon_8(k)*Opr_opsyp0_8(2*G_nj+j-1)+Opr_opsyp2_8(2*G_nj+j-1)
               end do
            end do
         end do
         if (Sol_miny == 1) then
            di_8 = Opr_opsyp0_8(G_nj+1+Lam_pil_s) / cos( G_yg_8(1+Lam_pil_s) )**2
            do i = 1+Lam_pil_w, G_ni-Lam_pil_e
               do k = Sol_mink,Sol_maxk
                  abc(1,k,i,1+Lam_pil_s)= 0.d0
                  abc(2,k,i,1+Lam_pil_s)= Opr_xeval_8(i) * di_8 + &
                                          F_hcon_8(k)*Opr_opsyp0_8(  G_nj+1+Lam_pil_s)+Opr_opsyp2_8(  G_nj+1+Lam_pil_s)
                  abc(3,k,i,1+Lam_pil_s)= F_hcon_8(k)*Opr_opsyp0_8(2*G_nj+1+Lam_pil_s)+Opr_opsyp2_8(2*G_nj+1+Lam_pil_s)
               end do
            end do
         endif
         if (Sol_maxy == G_nj) abc(3,:,:,G_nj-Lam_pil_n)= 0.d0
         
         if (Sol_miny == 1) then
            j= 1+Lam_pil_s
            do k= Sol_mink, Sol_maxk
               do i= 1+Lam_pil_w, G_ni-Lam_pil_e
                  F_bi_8(k,i,j) = 1.d0 / abc(2,k,i,j)
                  F_ci_8(k,i,j) = abc(3,k,i,j) * F_bi_8(k,i,j)
               end do
            end do
         else
            call MPI_RECV(F_ci_8(Sol_mink,1,Sol_miny-1), dim, MPI_DOUBLE_PRECISION, Numa_peerrank-1,&
                          tag+Numa_peerrank-1 , Numa_peercomm, status, err )
         endif

         jn= min(Sol_maxy,G_nj-Lam_pil_n-1)
         do j= j0, jn
            do k= Sol_mink, Sol_maxk
               do i= 1+Lam_pil_w, G_ni-Lam_pil_e
                  F_bi_8(k,i,j) = 1.d0 / ( abc(2,k,i,j)-abc(1,k,i,j) * F_ci_8(k,i,j-1) )
                  F_ci_8(k,i,j) = F_bi_8(k,i,j) * abc(3,k,i,j)
                  F_ai_8(k,i,j) = F_bi_8(k,i,j) * abc(1,k,i,j)
               end do
            end do
         end do
         if (Sol_maxy == G_nj) then
            j= G_nj-Lam_pil_n
            do k= Sol_mink, Sol_maxk
               do i= 1+Lam_pil_w, G_ni-Lam_pil_e
                  F_bi_8(k,i,j) = 1.d0 / ( abc(2,k,i,j)-abc(1,k,i,j) * F_ci_8(k,i,j-1) )
                  F_ai_8(k,i,j) = F_bi_8(k,i,j) * abc(1,k,i,j)
               end do
            end do
         else
            call MPI_SEND(F_ci_8(Sol_mink,1,jn), dim, MPI_DOUBLE_PRECISION, Numa_peerrank+1,&
                          tag+Numa_peerrank, Numa_peercomm, err )
         endif
      endif
      call MPI_Barrier ( Numa_sockcomm, err )
!
!     ---------------------------------------------------------------
!
      return
      end subroutine sol_prepabc 

