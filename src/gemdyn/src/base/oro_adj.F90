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

!**s/r oro_adj - Adjust orography

      subroutine oro_adj ()
      use gem_options
      use glb_ld
      use HORgrid_options
      use lam_options
      use step_options
      use mem_nest
      use gmm_geof
      use, intrinsic :: iso_fortran_env
      implicit none

      integer j
!
!     ---------------------------------------------------------------
!
      if ( .not. Grd_yinyang_L .and. .not. Lam_ctebcs_L) then
         if (Vtopo_mustadj_L) then
!$omp do
            do j=1, l_nj
               fis0 (1:l_ni,j)= nest_fullme(1:l_ni,j,1)
               orols(1:l_ni,j)= nest_fullme(1:l_ni,j,2)
            end do
!$omp enddo
            call oro_xchg_adj ()
         endif
      else
         if (Vtopo_L .and. (Lctl_step >= Vtopo_start)) then
!$omp single
            ! must create var_topo_omp
            call var_topo (fis0, orols, real(Lctl_step),&
                           l_minx,l_maxx,l_miny,l_maxy)
!$omp end single
            call oro_xchg_adj ()
         end if
      end if
!
!     ---------------------------------------------------------------
!
      return
      end subroutine oro_adj

      subroutine oro_xchg_adj ()
      use dyn_fisl_options
      use gem_options
      use glb_ld
      use glb_pil
      use dynkernel_options
      use HORgrid_options
      use lam_options
      use sol_options
      use step_options
      use mem_nest
      use gmm_geof
      use tdpack
      use, intrinsic :: iso_fortran_env
      implicit none

      integer j
!
!     ---------------------------------------------------------------
!
!$omp single
      if (Grd_yinyang_L) then
         call yyg_xchng (fis0, l_minx,l_maxx,l_miny,l_maxy, &
                         l_ni,l_nj, 1, .false., 'CUBIC', .true.)
         call yyg_xchng (orols, l_minx,l_maxx,l_miny,l_maxy, &
                         l_ni,l_nj, 1, .false., 'CUBIC', .true.)
      else
         call rpn_comm_xch_halo (fis0,l_minx,l_maxx,l_miny,l_maxy,&
                                 l_ni,l_nj,1,G_halox,G_haloy,&
                                 G_periodx,G_periody,l_ni,0)
         call rpn_comm_xch_halo (orols,l_minx,l_maxx,l_miny,l_maxy,&
                                 l_ni,l_nj,1,G_halox,G_haloy,&
                                 G_periodx,G_periody,l_ni,0)
      end if
!$omp end single
      if (Schm_sleve_L) then
         call update_sls_omp (orols,sls,l_minx,l_maxx,l_miny,l_maxy)
!$omp single
         if (Grd_yinyang_L) then
            call yyg_xchng (sls, l_minx,l_maxx,l_miny,l_maxy, &
                            l_ni,l_nj, 1, .false., 'CUBIC', .true.)
         else
            call rpn_comm_xch_halo (sls,l_minx,l_maxx,l_miny,l_maxy,&
                                    l_ni,l_nj,1,G_halox,G_haloy,&
                                     G_periodx,G_periody,l_ni,0)
         endif
!$omp end single
      endif

      if (Dynamics_hauteur_L) then
         call vertical_metric_omp (GVM, fis0, sls, l_minx,l_maxx,l_miny,l_maxy)

         if (trim(Sol_type_S) == 'ITERATIVE_3D') call matvec3d_init_hlt()
!$omp do
         do j= 1-G_haloy+1, l_nj+G_haloy-1
            me_full (1-G_halox+1:l_ni+G_halox-1,j) = fis0(1-G_halox+1:l_ni+G_halox-1,j) / grav_8
            me_large(1-G_halox+1:l_ni+G_halox-1,j) =  sls(1-G_halox+1:l_ni+G_halox-1,j) / grav_8
         end do
!$omp enddo
      endif
!
!     ---------------------------------------------------------------
!
      return
      endsubroutine oro_xchg_adj
