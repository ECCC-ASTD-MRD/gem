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

!** matvec_hlt - 3D Matrix-vector product subroutines (H coordinates)

      subroutine matvec_hlt ( F_vector, F_minx,F_maxx,F_miny,F_maxy,&
                              F_prod  , F_i0,F_in,F_j0,F_jn, F_nk )
      use dyn_fisl_options
      use geomh
      use gem_options
      use HORgrid_options
      use dynkernel_options
      use lam_options
      use glb_ld
      use ldnh
      use sol_mem
      use metric
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: F_minx,F_maxx,F_miny,F_maxy,F_i0,F_in,F_j0,F_jn,F_nk
      real(kind=REAL64), dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in) :: F_vector
      real(kind=REAL64), dimension(F_i0:F_in,F_j0:F_jn,F_nk), intent(out) :: F_prod

      integer :: i, j, k, k0, k0t, km, kp
!
!     ---------------------------------------------------------------
!
      k0=1+Lam_gbpil_T
      k0t=k0 ; if (Schm_opentop_L) k0t=k0-1

!$omp do collapse(2)
      do k = k0, l_nk
         do j=1+sol_pil_s, l_nj-sol_pil_n
            do i=1+sol_pil_w, l_ni-sol_pil_e
               fdg2(i,j,k)= F_vector(i,j,k)
            end do
         end do
      end do
!$omp enddo

!$omp do
      do j=1+sol_pil_s, l_nj-sol_pil_n
         do i=1+sol_pil_w, l_ni-sol_pil_e
            fdg2(i,j,l_nk+1)=GVM%mc_alfas_H_8(i,j) * F_vector(i,j,l_nk) &
                            -GVM%mc_betas_H_8(i,j) * F_vector(i,j,l_nk-1)
         end do
      end do
!$omp enddo

      if (Schm_opentop_L) then
!$omp do
         do j=1+sol_pil_s, l_nj-sol_pil_n
            do i=1+sol_pil_w, l_ni-sol_pil_e
               fdg2(i,j,k0t) = GVM%mc_alfat_8(i,j)* F_vector(i,j,k0)
            end do
         end do
!$omp enddo
      endif

!$omp single
      if ( Grd_yinyang_L) then
         call yyg_xchng (fdg2, l_minx,l_maxx,l_miny,l_maxy, &
                         l_ni,l_nj, l_nk+1, .false., 'CUBIC', .true.)
      else
         call rpn_comm_xch_halo(fdg2,l_minx,l_maxx,l_miny,l_maxy,&
                l_ni,l_nj,l_nk+1, 1,1,G_periodx,G_periody,l_ni,0 )
      endif
!$omp end single

!$omp do collapse(2)
      do k=k0,l_nk
         do j=1+sol_pil_s, l_nj-sol_pil_n
            km=max(k-1,1)
            kp=k+1
            do i=1+sol_pil_w, l_ni-sol_pil_e
               F_prod(i,j,k)= &
                  Sol_stencilh_8 (i,j,k, 1)*fdg2(i,  j,  k  ) &
               +  Sol_stencilh_8 (i,j,k, 2)*fdg2(i-1,j,  k  ) &
               +  Sol_stencilh_8 (i,j,k, 3)*fdg2(i+1,j,  k  ) &
               +  Sol_stencilh_8 (i,j,k, 4)*fdg2(i,  j,  km ) &
               +  Sol_stencilh_8 (i,j,k, 5)*fdg2(i,  j,  kp ) &
               +  Sol_stencilh_8 (i,j,k, 6)*fdg2(i-1,j,  km ) &
               +  Sol_stencilh_8 (i,j,k, 7)*fdg2(i-1,j,  kp ) &
               +  Sol_stencilh_8 (i,j,k, 8)*fdg2(i+1,j,  km ) &
               +  Sol_stencilh_8 (i,j,k, 9)*fdg2(i+1,j,  kp ) &
               +  Sol_stencilh_8 (i,j,k,10)*fdg2(i,  j-1,k  ) &
               +  Sol_stencilh_8 (i,j,k,11)*fdg2(i,  j+1,k  ) &
               +  Sol_stencilh_8 (i,j,k,12)*fdg2(i,  j-1,km ) &
               +  Sol_stencilh_8 (i,j,k,13)*fdg2(i,  j-1,kp ) &
               +  Sol_stencilh_8 (i,j,k,14)*fdg2(i,  j+1,km ) &
               +  Sol_stencilh_8 (i,j,k,15)*fdg2(i,  j+1,kp )
            end do
         end do
      end do
!$omp enddo
!     
!     ---------------------------------------------------------------
!     
      return
      end subroutine matvec_hlt
