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

!**s/r adz_post_tr - Apply a posteriori Mass-fixer/Shape-Preserving schemes or Finalization

      subroutine adz_post_tr ( F_post_S, F_post_itr, F_nptr )

      use adz_options
      use adz_mem
      use HORgrid_options
      use lun
      use tr3d

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      integer, intent(IN) :: F_nptr
      integer, dimension(F_nptr), intent(IN) :: F_post_itr
      character(len=4), dimension(F_nptr), intent(IN) :: F_post_S

      !object
      !===========================================================================
      !     Apply a posteriori Mass-fixer/Shape-Preserving schemes or Finalization
      !===========================================================================

      character(len=8) :: name_S
      logical :: Bermejo_Conde_L,Clip_L,ILMC_L
      integer :: n,itr,i,j,k,i0,in,j0,jn
!
!     ---------------------------------------------------------------
!
      do n=1,F_nptr

         i0 = Adz_i0 ; in = Adz_in ; j0 = Adz_j0 ; jn = Adz_jn

         itr = F_post_itr(n)

         !Bermejo-Conde LAM Flux ZLF
         !--------------------------
         if (Tr3d_mass(itr)==1.and.Adz_BC_LAM_flux==2) then

            i0 = Adz_i0b ; in = Adz_inb ; j0 = Adz_j0b ; jn = Adz_jnb

         end if

         if (Tr3d_mono(itr) == 1) then
            do k= Adz_k0, l_nk
               do j= j0,jn
                  do i= i0,in
                     Adz_stack(n)%dst(i,j,k) = &
                     max(Adz_post(i,j,k,(n-1)*3+2), &
                     min(Adz_post(i,j,k,(n-1)*3+3), &
                     Adz_stack(n)%dst(i,j,k)) )
                  end do
               end do
            end do
         end if

         name_S = 'TR/'//trim(F_post_S(n))//':M'

         !Initialization
         !--------------
         Bermejo_Conde_L = Tr3d_mass(itr)==1 .or. (Tr3d_mass(itr)>=111.and.Tr3d_mass(itr)<=139)

         Clip_L = Tr3d_mono(itr) == 1
         ILMC_L = Tr3d_mono(itr) == 2

         if (Bermejo_Conde_L) then

            Adz_BC_LEGACY_L = Tr3d_mass(itr) == 1

            Adz_BC_weight = 1
            Adz_BC_pexp_n = 0

            if (.not.Adz_BC_LEGACY_L) then

               Adz_BC_weight = (Tr3d_mass(itr) - 100)/10
               Adz_BC_pexp_n = (Tr3d_mass(itr) - (100 + Adz_BC_weight*10))

            end if

         end if

         !Printing
         !--------
         if (Adz_verbose>0.and.Lun_out>0) then

            write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'
            write(Lun_out,*) 'TRACERS: High-order SL advection: ',name_S(4:6)

            if (.not.ILMC_L) then

               if (.not.Clip_L) then
                  write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'
                  write(Lun_out,*) 'TRACERS: MONO (CLIPPING) is NOT activated: ',name_S(4:6)
               else
                  write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'
                  write(Lun_out,*) 'TRACERS: MONO (CLIPPING) is activated: ',name_S(4:6)
               end if

            end if

         end if

         if (.not.Bermejo_Conde_L.and..not.ILMC_L) cycle

         !Apply ILMC shape-preserving: Reset Monotonicity without changing Mass: Sorensen et al,ILMC, 2013,GMD
         !----------------------------------------------------------------------------------------------------
         if (ILMC_L) call ILMC_LAM ( n, itr, i0, in, j0, jn )

         !Apply Bermejo-Conde mass-fixer: Bermejo and Conde,2002,MWR and return
         !---------------------------------------------------------------------
         if (Bermejo_Conde_L) call Bermejo_Conde ( n, itr, i0, in, j0, jn, BC_mass_deficit(itr) )

      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_post_tr
