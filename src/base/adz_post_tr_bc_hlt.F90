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

!**s/r adz_post_tr - Apply a posteriori Bermejo-Conde Mass-fixer

      subroutine adz_post_tr_bc_hlt (F_wp)

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
      integer :: n,i,j,k,i0,in,j0,jn,k0,F_nptr,wp,deb,num,i0_c,in_c,j0_c,jn_c
      type(meta_tracers), dimension(:), pointer :: F_post
!
!---------------------------------------------------------------------
!
      LAM_L = .not.Grd_yinyang_L

      if (n_bc_TOTAL==0) return

      !--------------------------------------------------------------------
      !Prepare and Apply Bermejo-Conde mass-fixer for all tracers in Adz_bc
      !--------------------------------------------------------------------

      !Set CORE limits
      !---------------
      i0_c = 1 + pil_w ; in_c = l_ni - pil_e ; j0_c = 1 + pil_s ; jn_c = l_nj - pil_n

      i0 = i0_c ; in = in_c ; j0 = j0_c ; jn = jn_c ; k0 = Adz_k0t

      !Bermejo-Conde LAM Flux ZLF
      !--------------------------
      if (Adz_BC_LAM_flux==2) then

         i0 = Adz_i0b ; in = Adz_inb ; j0 = Adz_j0b ; jn = Adz_jnb ; k0 = 1

      end if

      do num=1,n_bc_TOTAL

         !Obtain Localization in Adz_bc for given tracer
         !----------------------------------------------
         n  = Adz_bc(num)%n
         wp = Adz_bc(num)%wp
!        print *, 'ADZ_POST_TR,num, n, wp:',num, n, wp

         !Set pointers when Tr_3CWP
         !-------------------------
         if (wp==1) then

              F_post => Tr_3CWP
            Adz_post => Adz_post_3CWP
            Adz_flux => Adz_flux_3CWP

            deb = Tr3d_debTRICUB_WP

         !Set pointers when Tr_BQWP
         !-------------------------
         else if (wp==2) then

              F_post => Tr_BQWP
            Adz_post => Adz_post_BQWP
            Adz_flux => Adz_flux_BQWP

            deb = Tr3d_debBICHQV_WP

         end if

!$omp single
         Adz_bc(num)%p  => tracers_P(deb+n-1)%pntr !As Pointer stack%src
         Adz_bc(num)%m  => tracers_M(deb+n-1)%pntr !As Pointer stack%dst

         Adz_bc(num)%lin => Adz_post(n)%lin
         Adz_bc(num)%min => Adz_post(n)%min
         Adz_bc(num)%max => Adz_post(n)%max

         if (LAM_L.and.Adz_BC_LAM_flux==1) then
            Adz_bc(num)%fo => Adz_flux(n)%fo
            Adz_bc(num)%fi => Adz_flux(n)%fi
         end if
!$omp end single

         !Obtain Bermejo-Conde's settings for given tracer
         !------------------------------------------------
         Adz_bc(num)%name_S = 'TR/'//trim(F_post(n)%name)//':M'

         Adz_bc(num)%mass_deficit = F_post(n)%BC_mass_deficit

         Adz_bc(num)%LEGACY_L = F_post(n)%mass == 1

         Adz_bc(num)%weight = 1
         Adz_bc(num)%pexp_n = 0

         if (.not.Adz_bc(num)%LEGACY_L) then
            Adz_bc(num)%weight = (F_post(n)%mass - 100)/10
            Adz_bc(num)%pexp_n = (F_post(n)%mass - (100 + Adz_bc(num)%weight*10))
         end if

      end do

      !Apply Bermejo-Conde mass-fixer for all tracers in Adz_bc
      !--------------------------------------------------------
      call Bermejo_Conde_hlt ( Adz_bc, n_bc_TOTAL, i0, in, j0, jn, k0 )

      !Store Mass deficit (EPSILON) induced when applying Bermejo-Conde
      !----------------------------------------------------------------
      do num=1,n_bc_TOTAL
!$omp single

         n  = Adz_bc(num)%n
         wp = Adz_bc(num)%wp

         if (wp==1) F_post => Tr_3CWP
         if (wp==2) F_post => Tr_BQWP

         F_post(n)%BC_mass_deficit = Adz_bc(num)%mass_deficit

!$omp end single
      end do

!      call gtmg_stop (38)
!
!---------------------------------------------------------------------
!
      return

      end subroutine adz_post_tr_bc_hlt
