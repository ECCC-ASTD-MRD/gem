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

!**s/r stat_mass_tracers - Calculate and print the mass of each tracer scaled by area

      subroutine stat_mass_tracers (F_time,F_comment_S)

      use adz_options
      use dyn_fisl_options
      use glb_ld
      use gmm_itf_mod
      use HORgrid_options
      use lam_options
      use lun
      use ptopo
      use tdpack
      use tr3d

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer,          intent(in) :: F_time       !Time 0 or Time 1
      character(len=*), intent(in) :: F_comment_S  !Comment

      !object
      !===============================================================
      !     Calculate and print the mass of each tracer scaled by area
      !===============================================================

      !------------------------------------------------------------------------------

      integer :: err,n,k0,count,i0,in,j0,jn,i0_sb,in_sb,j0_sb,jn_sb
      real*8  :: tracer_8
      logical :: do_subset_GY_L
      real, pointer, dimension (:,:,:) :: fld_tr
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk) :: air_mass,bidon,fld_ONE,w_tr
      character(len=21) :: type_S
      character(len= 7) :: time_S
      character(len=GMM_MAXNAMELENGTH) :: in_S

      !------------------------------------------------------------------------------

      i0 = 1+pil_w ; j0 = 1+pil_s ; in = l_ni-pil_e ; jn = l_nj-pil_n

      k0 = Lam_gbpil_T+1

      do_subset_GY_L = .false.

      !Evaluation of Mass on SUBSET of YIN available only if ONE processor on Yin
      !--------------------------------------------------------------------------
      if (Grd_yinyang_L.and.Adz_pil_sub_s>0.and.Ptopo_numproc==1) do_subset_GY_L = .true.

      i0_sb = 1+Adz_pil_sub_w ; j0_sb = 1+Adz_pil_sub_s ; in_sb = l_ni-Adz_pil_sub_e; jn_sb = l_nj-Adz_pil_sub_n

      !Terminator: Initialization
      !--------------------------
      call canonical_terminator_0 (count)

      call get_density (bidon,air_mass,F_time,l_minx,l_maxx,l_miny,l_maxy,l_nk,k0)

      if (F_time==1) time_S = "TIME T1"
      if (F_time==0) time_S = "TIME T0"

      type_S = "Mass of Mixing  (WET)"
      if (Schm_dry_mixing_ratio_L) type_S = "Mass of Mixing  (DRY)"

      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,*) ''

      !-------------
      !Treat tracers
      !-------------
      do n=1,Tr3d_ntr

         if (F_time==1) in_S = 'TR/'//trim(Tr3d_name_S(n))//':P'
         if (F_time==0) in_S = 'TR/'//trim(Tr3d_name_S(n))//':M'

         err = gmm_get(in_S, fld_tr)

         !Terminator: Prepare CLY
         !-----------------------
         call canonical_terminator_1 (fld_tr,in_S,count,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk)

         !Print the mass of each tracer scaled by area (CORE)
         !---------------------------------------------------
         call mass_tr (tracer_8,fld_tr,air_mass,l_minx,l_maxx,l_miny,l_maxy,l_nk,i0,in,j0,jn,k0)

         if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1002) 'TRACERS: ',type_S,time_S,'  C= ', &
                                                                  tracer_8/Adz_gc_area_8,Tr3d_name_S(n)(1:4),F_comment_S

         !Print the mass of each tracer scaled by area (SUBSET of YIN)
         !------------------------------------------------------------
         if (do_subset_GY_L) then

            w_tr(:,:,:) = fld_tr(:,:,:)

            if (Ptopo_couleur/=0) w_tr = 0.

            call mass_tr (tracer_8,w_tr,air_mass,l_minx,l_maxx,l_miny,l_maxy,l_nk,i0_sb,in_sb,j0_sb,jn_sb,k0)

            if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1002) 'TRACERS: ',type_S,time_S,' SB= ', &
                                                                     tracer_8/Adz_gs_area_8,Tr3d_name_S(n)(1:4),F_comment_S

         end if

         !Terminator: Print the mass of CLY scaled by area (CORE)
         !-------------------------------------------------------
         call canonical_terminator_2 (air_mass,tracer_8,count,l_minx,l_maxx,l_miny,l_maxy,l_nk,k0,Lun_out, &
                                      Ptopo_couleur,type_S,time_S,F_comment_S)

      end do

      !-----------------
      !Treat Air density
      !-----------------
      fld_ONE = 1.

      type_S = "Mass of Density (WET)"
      if (Schm_dry_mixing_ratio_L) type_S = "Mass of Density (DRY)"

      !Print Total Air mass scaled by area (CORE)
      !------------------------------------------
      call mass_tr (tracer_8,fld_ONE,air_mass,l_minx,l_maxx,l_miny,l_maxy,l_nk,i0,in,j0,jn,k0)

      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1002) 'TRACERS: ',type_S,time_S,'  C= ', &
                                                               tracer_8/Adz_gc_area_8,'RHO ',F_comment_S

      !Same scaling as in STAT_PSADJ
      !-----------------------------
      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1002) 'TRACERS: ',type_S,time_S,'  C= ', &
                                                               tracer_8/Adz_gc_area_8*grav_8,'RHOg',F_comment_S

      !Print Total Air mass scaled by area (SUBSET of YIN)
      !---------------------------------------------------
      if (do_subset_GY_L) then

         w_tr(:,:,:) = fld_ONE(:,:,:)

         if (Ptopo_couleur/=0) w_tr = 0.

         call mass_tr (tracer_8,w_tr,air_mass,l_minx,l_maxx,l_miny,l_maxy,l_nk,i0_sb,in_sb,j0_sb,jn_sb,k0)

         if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1002) 'TRACERS: ',type_S,time_S,' SB= ', &
                                                                  tracer_8/Adz_gs_area_8,'RHO ',F_comment_S

      end if

      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,*) ''

      return

1002  format(1X,A9,A21,1X,A7,A5,E19.12,1X,A4,1X,A16)

      end subroutine stat_mass_tracers
