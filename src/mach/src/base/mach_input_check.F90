!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2013 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------

!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet/Project : GEM-MACH
! Fichier/File   : mach_input_check.ftn90
! Creation       : S. Menard, Feb 2008.
! Description    : This subroutine is intend to do quality check of meteorological
!                  INPUT field. Here we modified the values of vertical diffusion
!                  obtain from the GEM model.
!
! Extra info     :
!
! Arguments      :
!                 IN
!                     metvar3d(i, k, MV3D_ZPLUS)       => Thermodynamic vertical levels height (m)
!                     metvar3d(i, k, MV3D_KT)           => Vertical diffusion from GEM
!                     metvar2d(:, MV2D_H)               => Boundary layer height
!                     landuse(:, lucprm)               => landuse
!                 OUT
!                     busvol(sm(sp_KTN)), or KT_NEW     => vertical diffusion for GEM-MACH
!
!============================================================================
!
!!if_on
subroutine mach_input_check(busvol, metvar2d, metvar3d, landuse)
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use mach_drydep_mod,      only: lucprm
!!if_off
   use chm_metvar_mod,       only: MV3D_KT, MV2D_H, MV3D_ZPLUS
   use chm_utils_mod,        only: ik
   use chm_nml_mod,          only: chm_vert_diff_s, chm_kt_minmax, &
                                   chm_urban_abl_min, chm_debug_2d_i
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_KTN, dbg_2d
   implicit none
!!if_on
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent   (in) :: metvar2d(chm_ni,SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d(chm_ni,chm_nk,SIZE_MV3D)
   real(kind=4),    intent   (in) :: landuse (chm_ni, lucprm)
!!if_off
!
!  Local variables.
!
   real(kind=4), dimension(chm_ni, chm_nk) :: kt_new
   real(kind=4)                            :: ktmin, urbfc, pblmin
   integer(kind=4)                         :: i, k
!
!
!  The following limiting values are set in chm_nml.ftn90:
!   chm_kt_minmax(1): minimum value (used for either all land use types or non-urban land use types)
!   chm_kt_minmax(2): maximum value (used for all land use types
!   chm_kt_minmax(3): minimum value for bilinear interpolation of minimum KT from urban
!                     land use fraction, chm_vert_diff_s = "RPNPHY_U"
   select case (chm_vert_diff_s)

!   Use urban land use fraction to determine the minimum KT value:
!   A bilinear interpolation dependent on the land-use fraction is used.
   case ('RPNPHY_U')
      do k = 1, chm_nk
         do i = 1, chm_ni
            if (chm_debug_2d_i >= 1 .and. k == chm_nk) &
               busvol(sm(dbg_2d(1)) % out_offset + i-1) = landuse(i,15)
            kt_new(i, k) = metvar3d(i, k, MV3D_KT)
            urbfc = max(min(landuse(i,15), 0.9999), 1.0e-5)
!  First, apply the lower value of KTMIN everywhere, up to the model's predicted PBL height:
            if (metvar3d(i, k, MV3D_ZPLUS) <= metvar2d(i, MV2D_H)) &
                kt_new(i, k) = max(chm_kt_minmax(1), metvar3d(i, k, MV3D_KT))
!  Next, when the urban land-use is more than 0.005:
!  (1) Determine the urban PBL height as the larger of the model's predicted PBL height sans TEB,
!      and the increment in the ABL from the average of the ABL height increments in Ren et al.
!  (2) Determine the new ktmin based on the bilinear interpolation between the non-ruban land-use
!      based on urban land use fraction
!  (3) Apply the new minimum up to the estimated urban PBL height based on Ren et al
            if (urbfc > 0.005) then
               pblmin = max(metvar2d(i, MV2D_H), chm_urban_abl_min)
               ktmin = chm_kt_minmax(1) * (1.0 - urbfc) + &
                       chm_kt_minmax(3) * urbfc
               if (metvar3d(i, k, MV3D_ZPLUS) <= pblmin) &
                  kt_new(i, k) = max(ktmin, kt_new(i,k))
            endif
         end do
      end do
!!
! Same minimum value of KT used everywhere, and only below the PBL height:
!
   case ('RPNPHY', 'RPNPHY_I')
      do k = 1, chm_nk
         do i = 1, chm_ni
            kt_new(i, k) = metvar3d(i, k, MV3D_KT)
            if (metvar3d(i, k, MV3D_ZPLUS) <= metvar2d(i, MV2D_H)) &
               kt_new(i, k) = max(chm_kt_minmax(1), metvar3d(i, k, MV3D_KT))
         end do
      end do
!
!  Set minimum and maximum values for vertical diffusion (KT) (default case)
!
   case ('FLUX', 'BOUNDARY')
   case DEFAULT
      do k = 1, chm_nk
         do i = 1, chm_ni
            if (metvar3d(i, k, MV3D_ZPLUS) <= metvar2d(i, MV2D_H)) then
               kt_new(i, k) = max(chm_kt_minmax(1), metvar3d(i, k, MV3D_KT))
            else
               kt_new(i, k) = min(metvar3d(i, k, MV3D_KT), chm_kt_minmax(2))
            end if
         end do
      end do
   end select

   do k = 1, chm_nk
      do i = 1, chm_ni
         busvol(sm(sp_KTN) % out_offset + ik(i, k, chm_ni)) = kt_new(i, k)
      end do
   end do

  return

end subroutine mach_input_check
