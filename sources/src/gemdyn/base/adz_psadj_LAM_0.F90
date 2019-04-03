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

!**s/r adz_psadj_LAM_0 - Estimate FLUX_out/FLUX_in based on Aranami et al. (2015) and
!                        Call psadj_LAM to adjust surface pressure

      subroutine adz_psadj_LAM_0 ()

      use adz_mem
      use adz_options
      use gmm_itf_mod
      use gmm_tracers

      implicit none

#include <arch_specific.hf>

      !object
      !=================================================================
      !     Estimate FLUX_out/FLUX_in based on Aranami et al. (2015) and
      !     Call psadj_LAM to adjust surface pressure
      !     ------------------------------------------------------------
      !     Reference: Aranami et al.,2015,QJRMS,141,1795-1803
      !=================================================================

      !-----------------------------------------------------------------

      integer :: empty_i,err
      real, dimension(1,1,1) :: empty

      !-----------------------------------------------------------------

      !Estimate FLUX_out/FLUX_in using Tracer=1
      !----------------------------------------
      Adz_Mass_Cons_tr_L = .true.
      Adz_BC_LAM_flux_n = 2
      Adz_pos_reset = 1

      call adz_cubic (empty,empty,Adz_pxyzt                     ,&
                      l_ni,l_nj,l_nk,l_minx,l_maxx,l_miny,l_maxy,&
                      empty_i,empty_i,empty_i,empty_i,Adz_k0,'t',.false.)

      !Reset to DEFAULT
      !----------------
      Adz_Mass_Cons_tr_L = .false.
      Adz_BC_LAM_flux_n = 0
      Adz_pos_reset = 0

      !Adjust surface pressure
      !-----------------------
      err = gmm_get(gmmk_cub_o_s,fld_cub_o)
      err = gmm_get(gmmk_cub_i_s,fld_cub_i)

      call psadj_LAM (fld_cub_o,fld_cub_i,l_minx,l_maxx,l_miny,l_maxy,l_nk,Adz_k0)

      return

      end subroutine adz_psadj_LAM_0
