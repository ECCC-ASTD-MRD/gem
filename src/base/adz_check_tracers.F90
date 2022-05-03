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

!**s/r adz_check_tracers - Initialize Variables for Mass Conservation/Shape-Preservation
!                          (not in NAMELIST)

      subroutine adz_check_tracers ()

      use adz_mem
      use adz_options
      use cstv
      use dyn_fisl_options
      use dynkernel_options
      use gem_options
      use geomh
      use HORgrid_options
      use ptopo
      use tdpack
      use tr3d

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !object
      !===========================================================================
      !     Initialize Variables for Mass Conservation/Shape-Preserving of Tracers
      !     (not in NAMELIST)
      !===========================================================================

      include 'mpif.h'

      include 'rpn_comm.inc'

      integer :: n,err,i,j,g_i0,g_in,g_j0,g_jn,i0_sb,in_sb,j0_sb,jn_sb,comm,flux_keep

      real(kind=REAL64)  :: c_area_8,s_area_8,gc_area_8

      logical :: BC_LAM_L,psadj_LAM_flux_L,do_subset_GY_L,BC_activated_L

      real(kind=REAL64), parameter :: QUATRO_8 = 4.0

      real(kind=REAL64) :: gathS(Ptopo_numproc)
!
!---------------------------------------------------------------------
!
      if (Adz_BC_LAM_flux<0.or.Adz_BC_LAM_flux>2) &
         call handle_error(-1,'ADZ_CHECK_TRACERS','BC_LAM_FLUX not valid-1')

      flux_keep = Adz_BC_LAM_flux

      BC_LAM_L = .false.

      do n=1,Tr3d_ntr

         if (Tr3d_intp(n)/='NONE'.and.Tr3d_intp(n)/='TRICUB'.and.Tr3d_intp(n)/='BICUBH_QV') &
            call handle_error(-1,'ADZ_CHECK_TRACERS','INTP not valid')

         BC_activated_L = Tr3d_mass(n)==1.or.(Tr3d_mass(n)>=111.and.Tr3d_mass(n)<=139)

         if (BC_activated_L.and..not.Grd_yinyang_L) BC_LAM_L = .true.

      end do

      if (BC_LAM_L.and.Adz_BC_LAM_flux==0) & 
         call handle_error(-1,'ADZ_CHECK_TRACERS','BC_LAM_FLUX not valid-2')

      if (Grd_yinyang_L) Adz_BC_LAM_flux = 0

      psadj_LAM_flux_L = .false.
      if (Schm_psadj>0.and..not.Grd_yinyang_L) psadj_LAM_flux_L = .true.

      if (.not.BC_LAM_L.and.psadj_LAM_flux_L.and.Adz_BC_LAM_flux==0) Adz_BC_LAM_flux = 1

      if (.not.BC_LAM_L.and..not.psadj_LAM_flux_L) Adz_BC_LAM_flux = 0

      if (Adz_verbose/=0.and.flux_keep/=Adz_BC_LAM_flux.and.Lun_out>0) &
         write(Lun_out,*) 'Revised Adz_BC_LAM_flux = ',Adz_BC_LAM_flux,'(INTERNAL PURPOSE)'

      !Bermejo-Conde LAM Flux=1/PSADJ LAM: Set mask_o/mask_i for Flux calculations based on Aranami et al. (2015)
      !----------------------------------------------------------------------------------------------------------
!!$      if (((BC_LAM_L.and.Adz_BC_LAM_flux==1).or.psadj_LAM_flux_L).and..not.ADZ_OD_L) then

         allocate (Adz_BC_LAM_mask_o(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,l_nk), &
                   Adz_BC_LAM_mask_i(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,l_nk))

         call adz_BC_LAM_setup (Adz_BC_LAM_mask_o,Adz_BC_LAM_mask_i,     &
                                Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy, &
                                l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk)

!!$      else if (((BC_LAM_L.and.Adz_BC_LAM_flux==1).or.psadj_LAM_flux_L).and.ADZ_OD_L) then
!!$
!!$         allocate (Adz_BC_LAM_mask_o(l_minx:l_maxx,l_miny:l_maxy,l_nk), &
!!$                   Adz_BC_LAM_mask_i(l_minx:l_maxx,l_miny:l_maxy,l_nk))
!!$
!!$         call adz_od_BC_LAM_setup (Adz_BC_LAM_mask_o,Adz_BC_LAM_mask_i, &
!!$                                   l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk)
!!$
!!$      end if

      !Check if Bermejo-Conde LAM Flux=2 is required for at least one tracer
      !---------------------------------------------------------------------
      if (BC_LAM_L.and.Adz_BC_LAM_flux==2) Adz_BC_LAM_zlf_L = .true.

      !Prepare SUBSET of Yin when PIL_SUB is prescribed in namelist adz_cfgs
      !---------------------------------------------------------------------
      do_subset_GY_L = Grd_yinyang_L.and.adz_pil_sub_s_g /= -1

      if (do_subset_GY_L) then

         g_i0 =    1 + adz_pil_sub_w_g
         g_in = G_ni - adz_pil_sub_e_g
         g_j0 =    1 + adz_pil_sub_s_g
         g_jn = G_nj - adz_pil_sub_n_g

         i0_sb = max(g_i0 - Ptopo_gindx(1,Ptopo_myproc+1) + 1, 1   -west *G_halox)
         in_sb = min(g_in - Ptopo_gindx(1,Ptopo_myproc+1) + 1, l_ni+east *G_halox)
         j0_sb = max(g_j0 - Ptopo_gindx(3,Ptopo_myproc+1) + 1, 1   -south*G_haloy)
         jn_sb = min(g_jn - Ptopo_gindx(3,Ptopo_myproc+1) + 1, l_nj+north*G_haloy)

         Adz_pil_sub_w = i0_sb - 1
         Adz_pil_sub_e = l_ni - in_sb
         Adz_pil_sub_s = j0_sb - 1
         Adz_pil_sub_n = l_nj - jn_sb

      end if

      comm = RPN_COMM_comm ('GRID')

      !Evaluate CORE/SUBSET areas
      !--------------------------
      if (Grd_yinyang_L) then

         gc_area_8 = QUATRO_8 * pi_8

         if (do_subset_GY_L) then

            s_area_8 = 0.0d0

            do j=j0_sb,jn_sb
            do i=i0_sb,in_sb
               s_area_8 = s_area_8 + geomh_area_8(i,j)
            end do
            end do

            call MPI_Allgather(s_area_8,1,MPI_DOUBLE_PRECISION,gathS,1,MPI_DOUBLE_PRECISION,comm,err)

            Adz_gs_area_8 = sum(gathS) 

         end if

      else

         c_area_8 = 0.0d0

         do j=1+pil_s,l_nj-pil_n
         do i=1+pil_w,l_ni-pil_e
            c_area_8 = c_area_8 + geomh_area_8(i,j)
         end do
         end do

         call MPI_Allgather(c_area_8,1,MPI_DOUBLE_PRECISION,gathS,1,MPI_DOUBLE_PRECISION,comm,err)

         gc_area_8 = sum(gathS)

      end if

      Adz_gc_area_8 = gc_area_8
!
!---------------------------------------------------------------------
!
      return
      end subroutine adz_check_tracers
