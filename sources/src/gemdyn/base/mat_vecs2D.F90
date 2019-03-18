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

!**s/r mat_vec2D - 2D_elliptic matrix_vector's computation
!

      subroutine  mat_vecs2D ( F_Sol, F_Rhs, Minx, Maxx, Miny, Maxy,nil,  &
                             njl,minx1, maxx1, minx2, maxx2,Nk,fdg1 )
      use gem_options
      use geomh
      use glb_ld
      use cstv
      use sol
      use opr
      implicit none
#include <arch_specific.hf>
!
      integer, intent(in):: Minx, Maxx, Miny, Maxy,nil, njl, &
                            minx1, maxx1, minx2, maxx2, Nk
      real*8, intent(out) ::  F_Rhs(Minx:Maxx,Miny:Maxy,Nk)
      real*8, intent(in) :: F_Sol(Minx:Maxx,Miny:Maxy,Nk)
      real*8, intent(out) :: fdg1(minx1:maxx1, minx2:maxx2, Nk)
!author
!       Abdessamad Qaddouri - December  2006
!
!revision
! v3_30 - Qaddouri A.       - initial version


      integer j,i,k,ii,jj,halox,haloy
      real*8  stencil1,stencil2,stencil3,stencil4,stencil5,cst,di_8
!
!     ---------------------------------------------------------------
!
      halox=1
      haloy=halox
!$omp parallel private (i,j,jj,ii,di_8,stencil1,&
!$omp        stencil2,stencil3,stencil4,stencil5,cst)
!$omp do
      do k = 1, nk
         fdg1(:,:,k) = .0d0
         do j=1+sol_pil_s, njl-sol_pil_n
         do i=1+sol_pil_w, nil-sol_pil_e
            fdg1(i,j,k)=F_Sol(i,j,k)
         end do
         end do
      end do
!$omp enddo

!$omp single
      call rpn_comm_xch_halon (fdg1,minx1,maxx1,minx2,maxx2,nil,njl, &
                         Nk,halox,haloy,G_periodx,G_periody,nil,0,2)
!$omp end single

!$omp do
      do k=1, NK
         cst= (Cstv_hco1_8+Cstv_hco0_8*Opr_zeval_8(k))
         do j=1+sol_pil_s, njl-sol_pil_n
            jj=j+l_j0-1
            di_8= Opr_opsyp0_8(G_nj+jj) * geomh_invcy2_8(j)
            do i=1+sol_pil_w, nil-sol_pil_e
               ii=i+l_i0-1
               stencil1=cst*Opr_opsxp0_8(G_ni+ii)* &
                 Opr_opsyp0_8(G_nj+jj) +Opr_opsxp2_8(G_ni+ii)*di_8+ &
                 Opr_opsxp0_8(G_ni+ii)*Opr_opsyp2_8(G_nj+jj)
               stencil2= Opr_opsxp2_8(ii)*di_8
               stencil3= Opr_opsxp2_8(2*G_ni+ii)*di_8
               stencil4= Opr_opsxp0_8(G_ni+ii)*Opr_opsyp2_8(jj)
               stencil5= Opr_opsxp0_8(G_ni+ii)*Opr_opsyp2_8(2*G_nj+jj)

               F_Rhs(i,j,k) =stencil1*fdg1(i  ,j  ,k) + &
                             stencil2*fdg1(i-1,j  ,k) + &
                             stencil5*fdg1(i  ,j+1,k) + &
                             stencil4*fdg1(i  ,j-1,k) + &
                             stencil3*fdg1(i+1,j  ,k)
            end do
         end do
      end do
!$omp enddo
!$omp end parallel
!
!     ---------------------------------------------------------------
!
      return
      end

