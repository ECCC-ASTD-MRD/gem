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

!**s/r adz_post_tr - Apply a posteriori Bermejo-Conde Mass-fixer/Clipping or ILMC Shape-Preserving schemes or Finalization

      subroutine adz_post_tr (F_wp)

      use adz_mem
      use adz_options
      use gem_timing
      use HORgrid_options
      use lun
      use mem_tracers
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
      integer :: n,i,j,k,i0,in,j0,jn,F_nptr,wp,deb,num,i0_c,in_c,j0_c,jn_c
      integer, save :: n_bc = 0
      integer, save :: n_bc_TOTAL
      logical, save :: done_store_BC_L = .false.
      type(meta_tracers), dimension(:), pointer :: F_post
!
!---------------------------------------------------------------------
!
      LAM_L = .not.Grd_yinyang_L

      if (F_wp /=0) then

         call gemtime_start (38, 'C_TR_POST', 33)

         !--------------------------------------------------------------------------------------------------
         !Apply Clipping or ILMC Shape-preserving and Store Localization in Adz_bc for all tracers using B-C
         !--------------------------------------------------------------------------------------------------

         if (F_wp==1) F_post => Tr_3CWP
         if (F_wp==2) F_post => Tr_BQWP

         if (F_wp==1) F_nptr = Tr3d_ntrTRICUB_WP
         if (F_wp==2) F_nptr = Tr3d_ntrBICHQV_WP

         do n=1,F_nptr

            !Initialization
            !--------------
            Bermejo_Conde_L = F_post(n)%mass==1 .or. (F_post(n)%mass>=111.and.F_post(n)%mass<=139)

            Clip_L = F_post(n)%mono == 1
            ILMC_L = F_post(n)%mono == 2

            !Apply Clipping Shape-preserving
            !-------------------------------
            if (Clip_L) then

               do k= Adz_k0, l_nk
                  do j= Adz_j0,Adz_jn
                     do i= Adz_i0,Adz_in
                        Adz_stack(n)%dst(i,j,k) = &
                        max(Adz_post(n)%min(i,j,k), &
                        min(Adz_post(n)%max(i,j,k), &
                        Adz_stack(n)%dst(i,j,k)) )
                     end do
                  end do
               end do

            end if

            name_tr_S = 'TR/'//trim(F_post(n)%name)//':M'

            if (Adz_verbose>0) call adz_post_tr_write (1)

            if (.not.Bermejo_Conde_L.and..not.ILMC_L) cycle

            !Apply ILMC Shape-preserving: Reset Monotonicity without changing Mass: Sorensen et al,ILMC, 2013,GMD
            !----------------------------------------------------------------------------------------------------
            if (ILMC_L) call ILMC_LAM ( n, name_tr_S, Adz_i0, Adz_in, Adz_j0, Adz_jn )

            if (.not.Bermejo_Conde_L) cycle

            if (done_store_BC_L) cycle

            !Store Localization in Adz_bc for all tracers using Bermejo-Conde
            !----------------------------------------------------------------
            n_bc = n_bc + 1

            Adz_bc(n_bc)%n  = n
            Adz_bc(n_bc)%wp = F_wp

         end do

         if (.not.done_store_BC_L) n_bc_TOTAL = n_bc

         if (Adz_verbose>0) call adz_post_tr_write (2)

         call gemtime_stop (38)

         return

      endif

      done_store_BC_L = .true.

      if (n_bc_TOTAL==0) return

      call gemtime_start (38, 'C_TR_POST', 33)

      !--------------------------------------------------------------------
      !Prepare and Apply Bermejo-Conde mass-fixer for all tracers in Adz_bc
      !--------------------------------------------------------------------

      !Set CORE limits
      !---------------
      i0_c = 1 + pil_w ; in_c = l_ni - pil_e ; j0_c = 1 + pil_s ; jn_c = l_nj - pil_n

      i0 = i0_c ; in = in_c ; j0 = j0_c ; jn = jn_c

      !Bermejo-Conde LAM Flux ZLF
      !--------------------------
      if (Adz_BC_LAM_flux==2) then

         i0 = Adz_i0b ; in = Adz_inb ; j0 = Adz_j0b ; jn = Adz_jnb

      end if

      do num=1,n_bc_TOTAL

         !Obtain Localization in Adz_bc for given tracer
         !----------------------------------------------
         n  = Adz_bc(num)%n
         wp = Adz_bc(num)%wp

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

         Adz_bc(num)%p  => tracers_P(deb+n-1)%pntr !As Pointer stack%src
         Adz_bc(num)%m  => tracers_M(deb+n-1)%pntr !As Pointer stack%dst

         Adz_bc(num)%lin => Adz_post(n)%lin
         Adz_bc(num)%min => Adz_post(n)%min
         Adz_bc(num)%max => Adz_post(n)%max

         if (LAM_L.and.Adz_BC_LAM_flux==1) then
            Adz_bc(num)%fo => Adz_flux(n)%fo
            Adz_bc(num)%fi => Adz_flux(n)%fi
         end if

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
      call Bermejo_Conde ( Adz_bc, n_bc_TOTAL, i0, in, j0, jn )

      !Store Mass deficit (EPSILON) induced when applying Bermejo-Conde
      !----------------------------------------------------------------
      do num=1,n_bc_TOTAL

         n  = Adz_bc(num)%n
         wp = Adz_bc(num)%wp

         if (wp==1) F_post => Tr_3CWP
         if (wp==2) F_post => Tr_BQWP

         F_post(n)%BC_mass_deficit = Adz_bc(num)%mass_deficit

      end do

      call gemtime_stop (38)
!
!---------------------------------------------------------------------
!
      return

contains

!**s/r adz_post_tr_write - Write adz_post_tr diagnostics based on F_numero if verbose is activated

      subroutine adz_post_tr_write (F_numero)

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

      end subroutine adz_post_tr_write

      end subroutine adz_post_tr
