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

!** s/r tfilt - Asselin time filter

      subroutine tfilt ()
      use glb_ld
      use glb_pil
      use gmm_vt2
      use gmm_vt1
      use gmm_vt0
      implicit none

      integer i,j,k
      integer i0, in, j0, jn
      real c1, grtf
!     
!     ---------------------------------------------------------------
!
      grtf = 0.05

      if (grtf.le.0.0) then
         return
      else
         c1=(1.-2.*grtf)
      endif
      
      i0= 1   +pil_w
      in= l_ni-pil_e
      j0= 1   +pil_s
      jn= l_nj-pil_n
!$omp do
      do k= 1,l_nk
         do j= j0, jn
         do i= i0, in
            ut1(i,j,k) = c1*ut1(i,j,k) + grtf*(ut0(i,j,k)+ut2(i,j,k))
            vt1(i,j,k) = c1*vt1(i,j,k) + grtf*(vt0(i,j,k)+vt2(i,j,k))
            tt1(i,j,k) = c1*tt1(i,j,k) + grtf*(tt0(i,j,k)+tt2(i,j,k))
            wt1(i,j,k) = c1*wt1(i,j,k) + grtf*(wt0(i,j,k)+wt2(i,j,k))
            zdt1(i,j,k)= c1*zdt1(i,j,k)+ grtf*(zdt0(i,j,k)+zdt2(i,j,k))
            qt1(i,j,k) = c1*qt1(i,j,k) + grtf*(qt0(i,j,k)+qt2(i,j,k))
         end do
         end do
      end do
!$omp enddo
!
!     ---------------------------------------------------------------
!
      return
      end subroutine tfilt
