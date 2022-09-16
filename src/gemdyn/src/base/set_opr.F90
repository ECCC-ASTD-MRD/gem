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

!**s/r set_opr - initialize the commons containing model operators

      subroutine set_opr
      use dynkernel_options
      use HORgrid_options
      use gem_options
      use tdpack
      use glb_ld
      use lun
      use glb_pil
      use sol_mem
      use opr
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer i, j, dim, Gni,Gnj,jj
      real(kind=REAL64)  sc_8, gdx_8, aab_8
      real(kind=REAL64), dimension(:)  ,allocatable :: wk_8, wk2_8
!
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write(Lun_out,1000)

!     Initialize projection operators to 0.d0

      dim = 3*G_ni
      allocate (Opr_opsxp0_8(dim), Opr_opsxp2_8(dim))
      dim = 3*G_nj
      allocate (Opr_opsyp0_8(dim), Opr_opsyp2_8(dim))
      dim = 3*G_nk
      allocate (Opr_opszp0_8(dim), Opr_opszpm_8(dim), &
                Opr_opszpl_8(dim), Opr_opszp2_8(dim))

      Opr_opsxp0_8 (1:3*G_ni) = 0.d0
      Opr_opsxp2_8 (1:3*G_ni) = 0.d0
      Opr_opsyp0_8 (1:3*G_nj) = 0.d0
      Opr_opsyp2_8 (1:3*G_nj) = 0.d0
      Opr_opszp0_8 (1:3*G_nk) = 0.d0
      Opr_opszpm_8 (1:3*G_nk) = 0.d0
      Opr_opszpl_8 (1:3*G_nk) = 0.d0
      Opr_opszp2_8 (1:3*G_nk) = 0.d0

!     Allocate memory for eigenvectors

      allocate (Opr_xeval_8(G_ni), Opr_zevec_8 (G_nk*G_nk), &
                Opr_zeval_8(G_nk), Opr_lzevec_8(G_nk*G_nk))

!     Dimension without pilot region

      Gni = G_ni-Lam_pil_w-Lam_pil_e
      Gnj = G_nj-Lam_pil_s-Lam_pil_n

!     Prepare projection operators

      dim = max(Gni,Gnj)
      allocate ( wk_8(dim) )

      do i = 1+Lam_pil_w, G_ni-Lam_pil_e
         Opr_opsxp0_8(G_ni+i) = (G_xg_8(i+1) - G_xg_8(i-1)) * 0.5d0
         wk_8(i-Lam_pil_w)    =  G_xg_8(i+1) - G_xg_8(i)
      end do

      allocate ( wk2_8(Gni*3) )
      call set_ops8 (wk2_8,wk_8,1.d0,G_periodx,Gni,Gni,1)
      do i=1,Gni
         Opr_opsxp2_8(i+Lam_pil_w)=wk2_8(i)
         Opr_opsxp2_8(G_ni+i+Lam_pil_w)=wk2_8(Gni+i)
         Opr_opsxp2_8(G_ni*2+i+Lam_pil_w)=wk2_8(Gni*2+i)
      end do
      if (Grd_yinyang_L) then
         Opr_opsxp2_8(G_ni  +1+Lam_pil_w)=&
                   2*Opr_opsxp2_8(G_ni  +1+Lam_pil_w)
         Opr_opsxp2_8(G_ni+Gni+Lam_pil_w)=&
                   2*Opr_opsxp2_8(G_ni+Gni+Lam_pil_w)
      end if
      deallocate ( wk2_8 )

      do j = 1+Lam_pil_s, G_nj-Lam_pil_n
         Opr_opsyp0_8(G_nj+j) = &
              sin((G_yg_8(j+1)+G_yg_8(j  )) * 0.5d0)-  &
              sin((G_yg_8(j  )+G_yg_8(j-1)) * 0.5d0)
         wk_8(j-Lam_pil_s) = (sin  (G_yg_8(j+1)) -sin(G_yg_8(j))) /  &
                             (cos ((G_yg_8(j+1)+G_yg_8(j)) * 0.5d0)**2)
      end do

      allocate ( wk2_8(Gnj*3) )
      call set_ops8(wk2_8,wk_8,0.d0,G_periody,Gnj,Gnj,1)
      do j=1,Gnj
         Opr_opsyp2_8(j+Lam_pil_s)       = wk2_8(      j)
         Opr_opsyp2_8(G_nj+j+Lam_pil_s)  = wk2_8(Gnj  +j)
         Opr_opsyp2_8(G_nj*2+j+Lam_pil_s)= wk2_8(Gnj*2+j)
      end do

      if (Grd_yinyang_L) then
         jj=1
         j=Lam_pil_s
         aab_8=     (sin (G_yg_8(j+1))-sin(G_yg_8(j))) / &
                    (cos ((G_yg_8(j+1)+G_yg_8(j)) * 0.5d0)**2)
         Opr_opsyp2_8(G_nj+jj+Lam_pil_s)=Opr_opsyp2_8(G_nj+jj+Lam_pil_s)&
                                         -1.0D0/aab_8
         jj=Gnj
         j=G_nj-1-Lam_pil_n+1
         aab_8=     (sin (G_yg_8(j+1))-sin(G_yg_8(j))) / &
                    (cos ((G_yg_8(j+1)+G_yg_8(j)) * 0.5d0)**2)

         Opr_opsyp2_8(G_nj+jj+Lam_pil_s)=Opr_opsyp2_8(G_nj+jj+Lam_pil_s)&
                                         -1.d0/aab_8
      end if

      deallocate ( wk_8, wk2_8 )

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P' .or. trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then

   !     Compute eigenvalues and eigenvector for the generalized
   !     eigenvalue problem in East-West direction

         sc_8 = pi_8 / dble( Gni )
         gdx_8 = (G_xg_8(G_ni-Lam_pil_e)-G_xg_8(Lam_pil_w) )/dble(Gni)

         if (Lun_debug_L) print *,'gdx=',gdx_8

         do i=1,1+Lam_pil_w
            Opr_xeval_8(i) = 0.d0
         end do

         do i=G_ni-Lam_pil_e+1,G_ni
            Opr_xeval_8(i) = 0.d0
         end do

         do i = 2+Lam_pil_w, G_ni-Lam_pil_e
            Opr_xeval_8(i) = - (2*sin(float(i-Lam_pil_w-1)*sc_8/2)/&
                                gdx_8)**2
         end do

         if (Grd_yinyang_L) then
            call set_poic  (Opr_xeval_8, Opr_opsxp0_8, &
                            Opr_opsxp2_8, Gni, G_ni)
         end if

      end if

 1000 format(/,'INITIALIZATING MODEL OPERATORS    (S/R SET_OPR)', &
             /,'=============================================')
!
!     ---------------------------------------------------------------
!
      return
      end
