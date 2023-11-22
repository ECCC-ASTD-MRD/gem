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
!**s/r yyg_initstencils - Computes Yin-Yang BC stencils for the elliptic solver
!


      subroutine yyg_initstencils()
      use gem_options
      use glb_ld
      use glb_pil
      use sol_mem
      use opr
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
!
!author
!     A. Qaddouri - initial MPI Yin-Yang version
!
!object
!     see ID section above
!
      real(kind=REAL64), parameter :: one  = 1.d0, half = .5d0
      real(kind=REAL64) di_8,xagauche_8,xadroite_8,yasud_8,yanord_8
      integer Gni,Gnj,i,j,ii,jj

!
!To compute Boundary Condition stencils for Yin-Yang

      Gni = G_ni-Lam_pil_w-Lam_pil_e
      Gnj = G_nj-Lam_pil_s-Lam_pil_n

      allocate (Sol_stencil2_8 (Gnj))
      allocate (Sol_stencil3_8 (Gnj))
      allocate (Sol_stencil4_8 (Gni))
      allocate (Sol_stencil5_8 (Gni))


      xagauche_8=one/(G_xg_8(Lam_pil_w+1) - G_xg_8(Lam_pil_w))
      xadroite_8=one/(G_xg_8(G_ni-Lam_pil_w+1) - G_xg_8(G_ni-Lam_pil_w))

      do j=1,Gnj
         jj=j+Lam_pil_s
         di_8 = sin((G_yg_8(jj+1)+G_yg_8(jj  ))* HALF)- &
                sin((G_yg_8(jj  )+G_yg_8(jj-1))* HALF)
         di_8 = Opr_opsyp0_8(G_nj+jj) / cos( G_yg_8 (jj) )**2

         Sol_stencil2_8(j)=di_8*xagauche_8
         Sol_stencil3_8(j)=di_8*xadroite_8
      end do

      di_8= (sin (G_yg_8(Lam_pil_s+1))-sin(G_yg_8(Lam_pil_s)))/ &
            (cos ((G_yg_8(Lam_pil_s+1)+G_yg_8(Lam_pil_s))*HALF)**2)
      yasud_8 =one/di_8

      di_8=(sin (G_yg_8(G_nj-Lam_pil_n+1))-sin(G_yg_8(G_nj-Lam_pil_n)))/ &
           (cos((G_yg_8(G_nj-Lam_pil_n+1)+G_yg_8(G_nj-Lam_pil_n))*HALF)**2)
      yanord_8=one/di_8

      do i=1,Gni
         ii=i+lam_pil_w
         di_8=Opr_opsxp0_8(G_ni+ii)
         Sol_stencil4_8(i)=yasud_8*di_8
         Sol_stencil5_8(i)=yanord_8*di_8

      end do

      return

      end
