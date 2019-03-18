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

!**s/r adv_a_posteriori - Apply a posteriori Mass-fixer/Shape-Preserving schemes or Finalization

      subroutine adv_a_posteriori ( F_name_S, F_output, F_high, F_mono, F_lin, F_min, F_max, F_input, &
                                    F_for_flux_o, F_for_flux_i, F_minx, F_maxx, F_miny, F_maxy, F_nk, &
                                    F_i0, F_in, F_j0, F_jn, F_k0, F_mono_kind, F_mass_kind )

      use adv_options
      use HORgrid_options
      use lun

      implicit none

#include <arch_specific.hf>

      character(len=*), intent(in) :: F_name_S                                       !Name of the interpolated field
      integer,          intent(in) :: F_nk                                           !Number of vertical levels
      integer,          intent(in) :: F_i0,F_in,F_j0,F_jn,F_k0                       !Scope of operator
      integer,          intent(in) :: F_mono_kind                                    !Kind of Shape-preserving
      integer,          intent(in) :: F_mass_kind                                    !Kind of Mass conservation
      integer,          intent(in) :: F_minx,F_maxx,F_miny,F_maxy                    !Dimension H
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(out)  :: F_output    !A posteriori solution
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in)   :: F_high      !High-order SL solution
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(inout):: F_mono      !High-order SL solution Monotonic
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in)   :: F_lin       !Linear SL solution
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in)   :: F_min       !MIN over cell
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in)   :: F_max       !MAX over cell
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in)   :: F_input     !Field at previous time step
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in)   :: F_for_flux_o!Advected mixing ratio with 0 in NEST
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in)   :: F_for_flux_i!Advected mixing ratio with 0 in CORE

      !object
      !===========================================================================
      !     Apply a posteriori Mass-fixer/Shape-Preserving schemes or Finalization
      !===========================================================================

      !---------------------------------------------------------------------------

      logical :: verbose_L, Lun_verbose_L
      logical :: Cubic_L, Bermejo_Conde_L, Clip_L, ILMC_L, BC_LAM_Aranami_L, LCSL_L
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk) :: high_w

      !---------------------------------------------------------------------------

      !Initialization
      !--------------
      verbose_L       = Adv_verbose /= 0
      Lun_verbose_L   = Lun_out>0.and.verbose_L

      Bermejo_Conde_L = F_mass_kind == 1
      LCSL_L          = F_mass_kind >  1
      BC_LAM_Aranami_L= Bermejo_Conde_L.and..not.Grd_yinyang_L.and.Adv_BC_LAM_flux==1
      Cubic_L         = .not.LCSL_L

      Clip_L          = F_mono_kind == 1
      ILMC_L          = F_mono_kind == 2

      !Printing
      !--------
      if (Cubic_L.and.Lun_verbose_L) then

         write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'
         write(Lun_out,*) 'TRACERS: Cubic SL advection: ',F_name_S(4:6)

         if (.not. ILMC_L) then

            if (.not.Clip_L) then
               write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'
               write(Lun_out,*) 'TRACERS: MONO (CLIPPING) is NOT activated: ',F_name_S(4:6)
            else
               write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'
               write(Lun_out,*) 'TRACERS: MONO (CLIPPING) is activated: ',F_name_S(4:6)
            end if

         end if

      else if (LCSL_L.and.Lun_verbose_L) then

         write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'
         write(Lun_out,*) 'TRACERS: Local Conservation SL advection (LCSL): ',F_name_S(4:6)

      end if

      !Finalize SL advection and return
      !--------------------------------
      if (.not.Bermejo_Conde_L.and..not.ILMC_L) then

         if (Clip_L) then
            F_output(F_i0:F_in,F_j0:F_jn,F_k0:F_nk) = F_mono(F_i0:F_in,F_j0:F_jn,F_k0:F_nk)
         else
            F_output(F_i0:F_in,F_j0:F_jn,F_k0:F_nk) = F_high(F_i0:F_in,F_j0:F_jn,F_k0:F_nk)
         end if

         return

      end if

      !Apply ILMC shape-preserving: Reset Monotonicity without changing Mass: Sorensen et al,ILMC, 2013,GMD
      !----------------------------------------------------------------------------------------------------
      if (ILMC_L) call ILMC_LAM (F_name_S,F_mono,F_high,F_min,F_max,F_minx,F_maxx,F_miny,F_maxy,F_nk, &
                                 F_i0,F_in,F_j0,F_jn,F_k0,Adv_ILMC_min_max_L,Adv_ILMC_sweep_max,verbose_L)

      !Apply Bermejo-Conde mass-fixer: Bermejo and Conde,2002,MWR and return
      !---------------------------------------------------------------------
      if (Bermejo_Conde_L) then

         high_w = F_high
         if (Clip_L.or.ILMC_L) high_w = F_mono

         call Bermejo_Conde (F_name_S,F_output,high_w,F_lin,F_min,F_max,F_input,F_for_flux_o,F_for_flux_i, &
                             F_minx,F_maxx,F_miny,F_maxy,F_nk,F_i0,F_in,F_j0,F_jn,F_k0, &
                             BC_LAM_Aranami_L,Adv_BC_min_max_L,verbose_L)

         return

      !Finalize ILMC without Bermejo-Conde mass-fixer and return
      !---------------------------------------------------------
      else

         F_output(F_i0:F_in,F_j0:F_jn,F_k0:F_nk) = F_mono(F_i0:F_in,F_j0:F_jn,F_k0:F_nk)

         return

      end if

      return

      end subroutine adv_a_posteriori
