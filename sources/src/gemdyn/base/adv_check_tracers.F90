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

!**s/r adv_check_tracers - Initialize Variables for Mass Conservation/Shape-Preserving of Chemical Tracers (not in NAMELIST)

      subroutine adv_check_tracers ()
      use adv_grid
      use adv_options
      use cstv
      use dyn_fisl_options
      use dynkernel_options
      use geomh
      use glb_ld
      use HORgrid_options
      use lun
      use ptopo
      use tdpack
      use tr3d

      implicit none

#include <arch_specific.hf>

      !object
      !====================================================================================
      !     Initialize Variables for Mass Conservation/Shape-Preserving of Chemical Tracers
      !     (not in NAMELIST)
      !====================================================================================

      !----------------------------------------------------

      integer :: n,err,i,j

      real*8  :: c_area_8,s_area_8,gc_area_8,gs_area_8

      logical :: BC_LAM_L,psadj_LAM_flux_L,do_subset_GY_L,BC_activated_L

      real*8, parameter :: QUATRO_8 = 4.0

      !----------------------------------------------------

      Adv_Mass_Cons_L   = .false.
      BC_LAM_L          = .false.
      Adv_LCSL_option_L = .false.

      do n=1,Tr3d_ntr

         BC_activated_L = Tr3d_mass(n)==1.or.(Tr3d_mass(n)>=111.and.Tr3d_mass(n)<=139)

         if (Tr3d_intp(n)=='QUINTIC'.and.Schm_autobar_L) &
            call handle_error(-1,'ADV_CHECK_TRACERS','INTP (QUINTIC_V) not valid when AUTOBAR')

         if (Tr3d_intp(n)/='NONE'.and.Tr3d_intp(n)/='CUBIC'.and.Tr3d_intp(n)/='QUINTIC') &
            call handle_error(-1,'ADV_CHECK_TRACERS','INTP not valid')

         if (Tr3d_mono(n)>1.or.Tr3d_mass(n)>0) Adv_Mass_Cons_L = .true.

         if (BC_activated_L.and..not.Grd_yinyang_L) BC_LAM_L = .true.

         if (.NOT.BC_activated_L.and.Tr3d_mass(n)>1) Adv_LCSL_option_L = .true.

      end do

      if (Adv_Mass_Cons_L.and..not.Schm_cub_traj_L) &
         call handle_error(-1,'ADV_CHECK_TRACERS','Schm_cub_traj_L should be .T. if Mass Conservation/Shape-preserving')

      psadj_LAM_flux_L = .false.
      if (Schm_psadj>0.and..not.Grd_yinyang_L) psadj_LAM_flux_L = .true.

      Adv_extension_L = Adv_LCSL_option_L.or.BC_LAM_L.or.psadj_LAM_flux_L

      if (Adv_extension_L.and.Lun_out>0) write(Lun_out,1000)

      !Bermejo-Conde LAM Flux=1/PSADJ LAM: Set mask_o/mask_i for Flux calculations based on Aranami et al. (2015)
      !----------------------------------------------------------------------------------------------------------
      if ((BC_LAM_L.and.Adv_BC_LAM_flux==1).or.psadj_LAM_flux_L) then

         allocate (adv_BC_LAM_mask_o(adv_lminx:adv_lmaxx,adv_lminy:adv_lmaxy,l_nk), &
                   adv_BC_LAM_mask_i(adv_lminx:adv_lmaxx,adv_lminy:adv_lmaxy,l_nk))

         call adv_BC_LAM_setup (adv_BC_LAM_mask_o,adv_BC_LAM_mask_i,     &
                                adv_lminx,adv_lmaxx,adv_lminy,adv_lmaxy, &
                                l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk)

      end if

      !If PIL_SUB is defined in namelist adv_cfgs, adjust according to topology
      !------------------------------------------------------------------------
      if (Adv_pil_sub_s /= -1) then
         if (.not.l_north) Adv_pil_sub_n = 0
         if (.not.l_south) Adv_pil_sub_s = 0
         if (.not.l_west ) Adv_pil_sub_w = 0
         if (.not.l_east ) Adv_pil_sub_e = 0
      end if

      !Consider SUBSET of YIN only if ONE processor on Yin
      !---------------------------------------------------
      do_subset_GY_L = .false.

      if (Grd_yinyang_L.and.Adv_pil_sub_s>0.and.Ptopo_numproc==1) do_subset_GY_L = .true.

      !Evaluate CORE/SUBSET areas
      !--------------------------
      if (Grd_yinyang_L) then

         gc_area_8 = QUATRO_8 * pi_8

         if (do_subset_GY_L) then

            s_area_8 = 0.0d0

            do j=1+Adv_pil_sub_s,l_nj-Adv_pil_sub_n
            do i=1+Adv_pil_sub_w,l_ni-Adv_pil_sub_e
               s_area_8 = s_area_8 + geomh_area_8(i,j)
            end do
            end do

            call rpn_comm_ALLREDUCE (s_area_8,gs_area_8,1,"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",err)

            Adv_gs_area_8 = gs_area_8

            if (Schm_autobar_L)  Adv_gs_area_8 = Adv_gs_area_8 * (Cstv_pref_8-Cstv_ptop_8) / grav_8

         end if

      else

         c_area_8 = 0.0d0

         do j=1+pil_s,l_nj-pil_n
         do i=1+pil_w,l_ni-pil_e
            c_area_8 = c_area_8 + geomh_area_8(i,j)
         end do
         end do

         call rpn_comm_ALLREDUCE(c_area_8,gc_area_8,1,"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",err )

      end if

      Adv_gc_area_8 = gc_area_8

      if (Schm_autobar_L)  Adv_gc_area_8 = Adv_gc_area_8 * (Cstv_pref_8-Cstv_ptop_8) / grav_8

      return

 1000 format(/,'EXTENDED ADVECTION OPERATIONS REQUIRED BEYOND CORE ZONE (S/R ADV_CHECK_TRACERS)',/,79('='))

      end subroutine adv_check_tracers
