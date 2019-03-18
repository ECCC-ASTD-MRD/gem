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

      subroutine set_opr (F_eigen_filename_S)
      use dynkernel_options
      use HORgrid_options
      use gem_options
      use tdpack
      use glb_ld
      use lun
      use glb_pil
      use fft
      use sol
      use opr
      use trp
      implicit none
#include <arch_specific.hf>
      character(len=*) F_eigen_filename_S

      real*8 ZERO_8, ONE_8, HALF_8
      parameter( ZERO_8 = 0.0 )
      parameter( ONE_8  = 1.0 )
      parameter( HALF_8 = 0.5 )

      integer i, j, dim, Gni,Gnj,jj
      real*8  sc_8, gdx_8, aab_8
      real*8, dimension(:)  ,allocatable :: wk_8, wk2_8
!
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write(Lun_out,1000)

!     Initialize projection operators to ZERO_8

      dim = 3*G_ni
      allocate (Opr_opsxp0_8(dim), Opr_opsxp2_8(dim))
      dim = 3*G_nj
      allocate (Opr_opsyp0_8(dim), Opr_opsyp2_8(dim))
      dim = 3*G_nk
      allocate (Opr_opszp0_8(dim), Opr_opszpm_8(dim), &
                Opr_opszpl_8(dim), Opr_opszp2_8(dim))

      Opr_opsxp0_8 (1:3*G_ni) = ZERO_8
      Opr_opsxp2_8 (1:3*G_ni) = ZERO_8
      Opr_opsyp0_8 (1:3*G_nj) = ZERO_8
      Opr_opsyp2_8 (1:3*G_nj) = ZERO_8
      Opr_opszp0_8 (1:3*G_nk) = ZERO_8
      Opr_opszpm_8 (1:3*G_nk) = ZERO_8
      Opr_opszpl_8 (1:3*G_nk) = ZERO_8
      Opr_opszp2_8 (1:3*G_nk) = ZERO_8

!     Allocate memory for eigenvectors

      if ( .not. Fft_fast_L ) allocate (Opr_xevec_8(G_ni*G_ni))
      allocate (Opr_xeval_8(G_ni), Opr_zevec_8 (G_nk*G_nk), &
                Opr_zeval_8(G_nk), Opr_lzevec_8(G_nk*G_nk))

!     Dimension without pilot region

      Gni = G_ni-Lam_pil_w-Lam_pil_e
      Gnj = G_nj-Lam_pil_s-Lam_pil_n

!     Prepare projection operators

      dim = max(Gni,Gnj)
      allocate ( wk_8(dim) )

      do i = 1+Lam_pil_w, G_ni-Lam_pil_e
         Opr_opsxp0_8(G_ni+i) = (G_xg_8(i+1) - G_xg_8(i-1))*HALF_8
         wk_8(i-Lam_pil_w)    =  G_xg_8(i+1) - G_xg_8(i)
      end do

      allocate ( wk2_8(Gni*3) )
      call set_ops8 (wk2_8,wk_8,ONE_8,G_periodx,Gni,Gni,1)
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
              sin((G_yg_8(j+1)+G_yg_8(j  ))* HALF_8)-  &
              sin((G_yg_8(j  )+G_yg_8(j-1))* HALF_8)
         wk_8(j-Lam_pil_s) = (sin  (G_yg_8(j+1)) -sin(G_yg_8(j))) /  &
                             (cos ((G_yg_8(j+1)+G_yg_8(j))*HALF_8)**2)
      end do

      allocate ( wk2_8(Gnj*3) )
      call set_ops8(wk2_8,wk_8,ZERO_8,G_periody,Gnj,Gnj,1)
      do j=1,Gnj
         Opr_opsyp2_8(j+Lam_pil_s)       = wk2_8(      j)
         Opr_opsyp2_8(G_nj+j+Lam_pil_s)  = wk2_8(Gnj  +j)
         Opr_opsyp2_8(G_nj*2+j+Lam_pil_s)= wk2_8(Gnj*2+j)
      end do

      if (Grd_yinyang_L) then
         jj=1
         j=Lam_pil_s
         aab_8=     (sin (G_yg_8(j+1))-sin(G_yg_8(j))) / &
                    (cos ((G_yg_8(j+1)+G_yg_8(j))*HALF_8)**2)
         Opr_opsyp2_8(G_nj+jj+Lam_pil_s)=Opr_opsyp2_8(G_nj+jj+Lam_pil_s)&
                                         -1.0D0/aab_8
         jj=Gnj
         j=G_nj-1-Lam_pil_n+1
         aab_8=     (sin (G_yg_8(j+1))-sin(G_yg_8(j))) / &
                    (cos ((G_yg_8(j+1)+G_yg_8(j))*HALF_8)**2)

         Opr_opsyp2_8(G_nj+jj+Lam_pil_s)=Opr_opsyp2_8(G_nj+jj+Lam_pil_s)&
                                         -1.d0/aab_8
      end if

      deallocate ( wk_8, wk2_8 )

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P' .or. trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then

         dim = (trp_12smax-trp_12smin+1)*(trp_22max-trp_22min+1)*G_nj
         allocate (Sol_ai_8(dim),Sol_bi_8(dim),Sol_ci_8(dim))

   !     Compute eigenvalues and eigenvector for the generalized
   !     eigenvalue problem in East-West direction

         if ( .not. Fft_fast_L ) then

            call set_poic2 (Opr_xeval_8, Opr_xevec_8 , Opr_opsxp0_8, &
                            Opr_opsxp2_8, Gni, G_ni, F_eigen_filename_S)

         else

            sc_8 = pi_8 / dble( Gni )
            gdx_8 = (G_xg_8(G_ni-Lam_pil_e)-G_xg_8(Lam_pil_w) )/dble(Gni)
            if (Lun_debug_L) print *,'gdx=',gdx_8
            do i=1,1+Lam_pil_w
               Opr_xeval_8(i)    = ZERO_8
            end do
            do i=G_ni-Lam_pil_e+1,G_ni
               Opr_xeval_8(i)    = ZERO_8
            end do
            do i = 2+Lam_pil_w, G_ni-Lam_pil_e
               Opr_xeval_8(i) = - (2*sin(float(i-Lam_pil_w-1)*sc_8/2)/&
                                   gdx_8)**2
            end do
            if (Grd_yinyang_L) then
               allocate (Opr_xevec_8(G_ni*G_ni))
               call set_poic2  (Opr_xeval_8, Opr_xevec_8 , Opr_opsxp0_8,&
                                Opr_opsxp2_8, Gni, G_ni, F_eigen_filename_S)
            end if

         end if

      end if

 1000 format(/,'INITIALIZATING MODEL OPERATORS    (S/R SET_OPR)', &
             /,'=============================================')
!
!     ---------------------------------------------------------------
!
      return
      end
