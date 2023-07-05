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

!**s/r adz_post_tr_mono - Apply Clipping or ILMC Shape-Preserving schemes or Finalization

      subroutine adz_post_tr_mono_hlt (F_wp)

      use adz_mem
      use adz_options
      use omp_timing
      use HORgrid_options
      use lun
      use mem_tracers
      use masshlt
      use tr3d

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !object
      !==========================================================================================================
      !     Apply a posteriori Bermejo-Conde Mass-fixer/Clipping or ILMC Shape-Preserving schemes or Finalization
      !     -----------------------------------------------------------------------------------------------------
      !     F_wp=1: Apply Clipping/ILMC for Tr_3CWP and Store Localization in Adz_bc for all tracers using B-C
      !     F_wp=2: Apply Clipping/ILMC for Tr_BQWP and Store Localization in Adz_bc for all tracers using B-C
      !     F_wp=0: Apply Bermejo-Conde mass-fixer for all tracers in Adz_bc
      !==========================================================================================================

      integer, intent(in) :: F_wp

      character(len=8) :: name_tr_S
      logical :: Bermejo_Conde_L,Clip_L,ILMC_L,LAM_L
      integer :: n,i,j,k,i0,in,j0,jn,k0,F_nptr,wp,deb,num,i0_c,in_c,j0_c,jn_c,n_bc
      type(meta_tracers), dimension(:), pointer :: F_post
!!$omp threadprivate()
!
!---------------------------------------------------------------------
!
      LAM_L = .not.Grd_yinyang_L

!      call gtmg_start (38, 'C_TR_POST', 33)

      !--------------------------------------------------------------------------------------------------
      !Apply Clipping or ILMC Shape-preserving and Store Localization in Adz_bc for all tracers using B-C
      !--------------------------------------------------------------------------------------------------

      if (F_wp==1) F_post => Tr_3CWP
      if (F_wp==2) F_post => Tr_BQWP

      if (F_wp==1) F_nptr = Tr3d_ntrTRICUB_WP
      if (F_wp==2) F_nptr = Tr3d_ntrBICHQV_WP

      n_bc = 0
      do n=1,F_nptr

            !Initialization
            !--------------
            Bermejo_Conde_L = F_post(n)%mass==1 .or. (F_post(n)%mass>=111.and.F_post(n)%mass<=139)

            Clip_L = F_post(n)%mono == 1
            ILMC_L = F_post(n)%mono == 2

            !Apply Clipping Shape-preserving
            !-------------------------------
            if (Clip_L) then
!$omp single
               do k= Adz_k0t, l_nk
                  do j= Adz_j0,Adz_jn
                     do i= Adz_i0,Adz_in
                        Adz_stack(n)%dst(i,j,k) = &
                        max(Adz_post(n)%min(i,j,k), &
                        min(Adz_post(n)%max(i,j,k), &
                        Adz_stack(n)%dst(i,j,k)) )
                     end do
                  end do
               end do
!$omp end single
            end if

            name_tr_S = 'TR/'//trim(F_post(n)%name)//':M'

!$omp single
            if (Adz_verbose>0) call adz_post_tr_write_hlt (1)
!$omp end single

            !Apply ILMC Shape-preserving: Reset Monotonicity without changing Mass: Sorensen et al,ILMC, 2013,GMD
            !----------------------------------------------------------------------------------------------------
            if (ILMC_L) call ILMC_LAM_hlt ( n, name_tr_S, Adz_i0, Adz_in, Adz_j0, Adz_jn )

            !Store Localization in Adz_bc for all tracers using Bermejo-Conde
            !----------------------------------------------------------------
            if (Bermejo_Conde_L) then

              n_bc = n_bc + 1

!$omp single
              Adz_bc(n_bc)%n  = n
              Adz_bc(n_bc)%wp = F_wp
!$omp end single

            endif

      end do

      n_bc_TOTAL = n_bc

!     print *,  'n_bc_TOTAL =', n_bc_TOTAL

!$omp single
      if (Adz_verbose>0) call adz_post_tr_write_hlt (2)
!$omp end single

!     call gtmg_stop (38)
!
!---------------------------------------------------------------------
!
      return

contains

!**s/r adz_post_tr_write - Write adz_post_tr diagnostics based on F_numero if verbose is activated

      subroutine adz_post_tr_write_hlt (F_numero)

      implicit none

      !arguments
      !---------
      integer, intent(in) :: F_numero

      if (F_numero==1.and.Lun_out>0) then

         write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'
         write(Lun_out,*) 'TRACERS: High-order SL advection: ',name_tr_S(4:6)

         if (.not.ILMC_L) then

            if (.not.Clip_L) then
               write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'
               write(Lun_out,*) 'TRACERS: MONO (CLIPPING) is NOT activated: ',name_tr_S(4:6)
            else
               write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'
               write(Lun_out,*) 'TRACERS: MONO (CLIPPING) is activated: ',name_tr_S(4:6)
            end if

         end if

      elseif (F_numero==2.and.Lun_out>0) then

         write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'
         write(Lun_out,*) ''

      end if

      return

      end subroutine adz_post_tr_write_hlt

      end subroutine adz_post_tr_mono_hlt
