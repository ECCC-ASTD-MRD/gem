!---------------------------------- LICENCE BEGIN ------------------------------- ! GEM - Library of kernel routines for the GEM numerical atmospheric model
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

!**s/r adv_BC_LAM_Aranami: Estimate FLUX_out/FLUX_in based on Aranami et al. (2015)

      subroutine adv_BC_LAM_Aranami (F_cub_o,F_in_o,F_cub_i,F_in_i,F_x,F_y,F_z,F_num,F_k0,F_nk,F_lev)

      use adv_grid
      use adv_interp
      use adz_options
      use gem_options
      use glb_ld
      use HORgrid_options
      use lun
      use ptopo
      use ver

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      character(len=*),      intent(in)  :: F_lev         ! m/t : Momemtum/thermo level
      integer,               intent(in)  :: F_num         ! Number of points
      integer,               intent(in)  :: F_nk          ! Number of vertical levels
      integer,               intent(in)  :: F_k0          ! Scope of operator
      real,dimension(F_num), intent(in)  :: F_x, F_y, F_z ! Interpolation target x,y,z coordinates
      real,dimension(*),     intent(in)  :: F_in_o        ! Field to interpolate   FLUX_out
      real,dimension(*),     intent(in)  :: F_in_i        ! Field to interpolate   FLUX_in
      real,dimension(F_num), intent(out) :: F_cub_o       ! High-order SL solution FLUX_out
      real,dimension(F_num), intent(out) :: F_cub_i       ! High-order SL solution FLUX_in

      !object
      !=================================================================================
      !     Estimate FLUX_out/FLUX_in based on Aranami et al.,QJRMS,141,1795-1803 (2015)
      !=================================================================================

      !-----------------------------------------------------------------------------------------------------------

      integer, save :: nind_w_o,nind_s_o,nind_e_o,nind_n_o,nind_w_i,nind_s_i,nind_e_i,nind_n_i, &
                       i0_w_i,in_w_i,j0_w_i,jn_w_i,i0_w_o,in_w_o,j0_w_o,jn_w_o,                 &
                       i0_s_i,in_s_i,j0_s_i,jn_s_i,i0_s_o,in_s_o,j0_s_o,jn_s_o,                 &
                       i0_e_i,in_e_i,j0_e_i,jn_e_i,i0_e_o,in_e_o,j0_e_o,jn_e_o,                 &
                       i0_n_i,in_n_i,j0_n_i,jn_n_i,i0_n_o,in_n_o,j0_n_o,jn_n_o,                 &
                       num_PE_near_NE,num_PE_near_NW, &
                       num_PE_near_SE,num_PE_near_SW,num_e,i_near_E,i_near_W

      integer, save, allocatable, dimension(:), target :: ii_w_o, ii_s_o, ii_e_o, ii_n_o, &
                                                          ii_w_i, ii_s_i, ii_e_i, ii_n_i

      real,    save, allocatable, dimension(:) :: px,py,pz

      logical, save :: done_BC_LAM_alloc_L = .false.

      logical, save :: narrow_i_w_L,narrow_j_s_L,narrow_i_e_L,narrow_j_n_L,narrow_L

      !-----------------------------------------------------------------------------------------------------------

      integer :: nind_o,nind_i,n0,nx,ny,nz,m1,id,n,kkmax,o1,o2,o3,o4, &
                 i0_e,in_e,j0_e,jn_e,jext,err,i0_e_out,in_e_out,j0_e_out,jn_e_out, &
                 g_narrow_i_e,g_narrow_j_n,narrow_i_e,narrow_j_n, &
                 g_narrow_i_w,g_narrow_j_s,narrow_i_w,narrow_j_s, &
                 minx_e,maxx_e,miny_e,maxy_e,i_near,ne,i,j,k,i0_i,in_i,j0_i,jn_i

      integer, pointer, dimension(:) :: ii_o, ii_i

      real(kind=REAL64) :: a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4, &
                d1,d2,d3,d4,p1,p2,p3,p4,ra,rb,rc,rd

      real(kind=REAL64), dimension(:), pointer :: p_bsz_8,p_zbc_8,p_zabcd_8,p_zbacd_8,p_zcabd_8,p_zdabc_8
      integer,dimension(:), pointer :: p_lcz

      logical :: zcubic_L

      real :: i_h(adv_lminx:adv_lmaxx,adv_lminy:adv_lmaxy,F_nk)

      !---------------------------------------------------------------------------------------------------------

      nullify(ii_o, ii_i)

      kkmax = F_nk - 1

      if (F_lev == 'm') then
        p_lcz     => adv_lcz%m
        p_bsz_8   => adv_bsz_8%m
        p_zabcd_8 => adv_zabcd_8%m
        p_zbacd_8 => adv_zbacd_8%m
        p_zcabd_8 => adv_zcabd_8%m
        p_zdabc_8 => adv_zdabc_8%m
        p_zbc_8   => adv_zbc_8%m
      else if (F_lev  == 't') then
        p_lcz     => adv_lcz%t
        p_bsz_8   => adv_bsz_8%t
        p_zabcd_8 => adv_zabcd_8%t
        p_zbacd_8 => adv_zbacd_8%t
        p_zcabd_8 => adv_zcabd_8%t
        p_zdabc_8 => adv_zdabc_8%t
        p_zbc_8   => adv_zbc_8%t
      else if (F_lev == 'x') then
        p_lcz     => adv_lcz%x
        p_bsz_8   => adv_bsz_8%x
        p_zabcd_8 => adv_zabcd_8%x
        p_zbacd_8 => adv_zbacd_8%x
        p_zcabd_8 => adv_zcabd_8%x
        p_zdabc_8 => adv_zdabc_8%x
        p_zbc_8   => adv_zbc_8%x
      end if

      jext = Grd_maxcfl + 1

      !E grid: HALO_E and extended advection limits (in nesting zone)
      !--------------------------------------------------------------
      minx_e = adv_lminx
      maxx_e = adv_lmaxx
      miny_e = adv_lminy
      maxy_e = adv_lmaxy

      call adz_get_ij0n_ext (i0_e,in_e,j0_e,jn_e,1) !EXTENSION (CFL)

      !We do not modify how FLUX_out is calculated
      !-------------------------------------------
      i0_e_out = i0_e
      in_e_out = in_e
      j0_e_out = j0_e
      jn_e_out = jn_e

      if (.not.l_west)  i0_e = minx_e
      if (.not.l_east)  in_e = maxx_e
      if (.not.l_south) j0_e = miny_e
      if (.not.l_north) jn_e = maxy_e

      !Computations related to upstream positions done at each timestep
      !----------------------------------------------------------------
      if (Adz_pos_reset==1) then

         !Allocation
         !----------
         if (.not.done_BC_LAM_alloc_L) then

            !Check if PEs at WEST/SOUTH/EAST/NORTH boundaries are too narrow to calculate FLUX_in
            !------------------------------------------------------------------------------------
            narrow_i_w = 0
            narrow_j_s = 0
            narrow_i_e = 0
            narrow_j_n = 0

            if (l_west) then

               in_w_i = 1+pil_w+jext

               if (in_w_i > l_ni) narrow_i_w = 1

            end if

            if (l_south) then

               jn_s_i = 1+pil_s+jext

               if (jn_s_i > l_nj) narrow_j_s = 1

            end if

            if (l_east) then

               i0_e_i = l_ni-pil_e-jext

               if (i0_e_i < 1) narrow_i_e = 1

            end if

            if (l_north) then

               j0_n_i = l_nj-pil_n-jext

               if (j0_n_i < 1) narrow_j_n = 1

            end if

            call rpn_comm_ALLREDUCE(narrow_i_w,g_narrow_i_w,1,"MPI_INTEGER","MPI_MAX","GRID",err)
            call rpn_comm_ALLREDUCE(narrow_j_s,g_narrow_j_s,1,"MPI_INTEGER","MPI_MAX","GRID",err)
            call rpn_comm_ALLREDUCE(narrow_i_e,g_narrow_i_e,1,"MPI_INTEGER","MPI_MAX","GRID",err)
            call rpn_comm_ALLREDUCE(narrow_j_n,g_narrow_j_n,1,"MPI_INTEGER","MPI_MAX","GRID",err)

            narrow_i_w_L = g_narrow_i_w/=0
            narrow_j_s_L = g_narrow_j_s/=0
            narrow_i_e_L = g_narrow_i_e/=0
            narrow_j_n_L = g_narrow_j_n/=0

            narrow_L = narrow_i_w_L.or.narrow_j_s_L.or.narrow_i_e_L.or.narrow_j_n_L

            if (Lun_out>0) then
                                 write(Lun_out,*) ''
               if (narrow_i_w_L) write(Lun_out,*) 'ADV_TRICUB_LAG3D_FLUX: PEs WEST  are narrow: We communicate with neighbors'
               if (narrow_j_s_L) write(Lun_out,*) 'ADV_TRICUB_LAG3D_FLUX: PEs SOUTH are narrow: We communicate with neighbors'
               if (narrow_i_e_L) write(Lun_out,*) 'ADV_TRICUB_LAG3D_FLUX: PEs EAST  are narrow: We communicate with neighbors'
               if (narrow_j_n_L) write(Lun_out,*) 'ADV_TRICUB_LAG3D_FLUX: PEs NORTH are narrow: We communicate with neighbors'
                                 write(Lun_out,*) ''
            end if

            if (narrow_L) then

               num_e = (maxx_e-minx_e+1)*(maxy_e-miny_e+1)*F_nk

               allocate (px(num_e),py(num_e),pz(num_e))

            end if

            if (l_west) then

               i0_w_o = i0_e_out
               in_w_o = pil_w
               j0_w_o = j0_e_out
               jn_w_o = jn_e_out

               nind_w_o = (in_w_o-i0_w_o+1)*(jn_w_o-j0_w_o+1)*(F_nk-F_k0+1)

               allocate (ii_w_o(4*nind_w_o))

               i0_w_i = 1+pil_w
               in_w_i = 1+pil_w+jext
               j0_w_i = 1+pil_s
               jn_w_i = l_nj-pil_n

               if (narrow_i_w_L) in_w_i = max(in_w_i,maxx_e)

               nind_w_i = (in_w_i-i0_w_i+1)*(jn_w_i-j0_w_i+1)*(F_nk-F_k0+1)

               allocate (ii_w_i(4*nind_w_i))

            end if

            if (l_south) then

               i0_s_o = i0_e_out
               in_s_o = in_e_out
               j0_s_o = j0_e_out
               jn_s_o = pil_s

               nind_s_o = (in_s_o-i0_s_o+1)*(jn_s_o-j0_s_o+1)*(F_nk-F_k0+1)

               allocate (ii_s_o(4*nind_s_o))

               i0_s_i = 1+pil_w
               in_s_i = l_ni-pil_e
               j0_s_i = 1+pil_s
               jn_s_i = 1+pil_s+jext

               if (narrow_j_s_L) jn_s_i = max(jn_s_i,maxy_e)

               nind_s_i = (in_s_i-i0_s_i+1)*(jn_s_i-j0_s_i+1)*(F_nk-F_k0+1)

               allocate (ii_s_i(4*nind_s_i))

            end if

            if (l_east) then

               i0_e_o = l_ni-pil_e+1
               in_e_o = in_e_out
               j0_e_o = j0_e_out
               jn_e_o = jn_e_out

               nind_e_o = (in_e_o-i0_e_o+1)*(jn_e_o-j0_e_o+1)*(F_nk-F_k0+1)

               allocate (ii_e_o(4*nind_e_o))

               i0_e_i = l_ni-pil_e-jext
               in_e_i = l_ni-pil_e
               j0_e_i = 1+pil_s
               jn_e_i = l_nj-pil_n

               if (narrow_i_e_L) i0_e_i = max(i0_e_i,minx_e)

               nind_e_i = (in_e_i-i0_e_i+1)*(jn_e_i-j0_e_i+1)*(F_nk-F_k0+1)

               allocate (ii_e_i(4*nind_e_i))

            end if

            if (l_north) then

               i0_n_o = i0_e_out
               in_n_o = in_e_out
               j0_n_o = l_nj-pil_n+1
               jn_n_o = jn_e_out

               nind_n_o = (in_n_o-i0_n_o+1)*(jn_n_o-j0_n_o+1)*(F_nk-F_k0+1)

               allocate (ii_n_o(4*nind_n_o))

               i0_n_i = 1+pil_w
               in_n_i = l_ni-pil_e
               j0_n_i = l_nj-pil_n-jext
               jn_n_i = l_nj-pil_n

               if (narrow_j_n_L) j0_n_i = max(j0_n_i,miny_e)

               nind_n_i = (in_n_i-i0_n_i+1)*(jn_n_i-j0_n_i+1)*(F_nk-F_k0+1)

               allocate (ii_n_i(4*nind_n_i))

            end if

            num_PE_near_SW =  1 + 1
            num_PE_near_SE =  Ptopo_npex - 1 - 1
            num_PE_near_NW = (Ptopo_npex * Ptopo_npey -1) - Ptopo_npex + 1
            num_PE_near_NE = (Ptopo_npex * Ptopo_npey -1) - 1

            !Initialize i_near_E = i0_e_i of PE at E
            !---------------------------------------
            if (narrow_i_e_L.and.(narrow_j_s_L.or.narrow_j_n_L)) then
               i_near = i0_e_i
               call RPN_COMM_allreduce(i_near,i_near_E,1,"mpi_integer","mpi_min","GRID",err)
            end if

            !Initialize i_near_W = in_w_i of PE at W
            !---------------------------------------
            if (narrow_i_w_L.and.(narrow_j_s_L.or.narrow_j_n_L)) then
               i_near = in_w_i
               call RPN_COMM_allreduce(i_near,i_near_W,1,"mpi_integer","mpi_min","GRID",err)
            end if

            done_BC_LAM_alloc_L = .true.

         end if !END Allocation
         !---------------------

         !Exchange HALO_E for positions
         !-----------------------------
         if (narrow_L) then

!$omp parallel private(i,j,k,n,ne)&
!$omp shared(px,py,pz,minx_e,maxx_e,miny_e,maxy_e)
!$omp do
            do k=F_k0,F_nk
               do j=1,l_nj
                  do i=1,l_ni

                     n  = (k-1)*l_ni*l_nj + (j-1)*l_ni + i
                     ne = (k-1)*(maxy_e-miny_e+1)*(maxx_e-minx_e+1) + (j-miny_e)*(maxx_e-minx_e+1) + i-minx_e + 1

                     px(ne) = F_x(n)
                     py(ne) = F_y(n)
                     pz(ne) = F_z(n)

                  end do
               end do
            end do
!$omp end do
!$omp end parallel

            call rpn_comm_xch_halo(px,minx_e,maxx_e,miny_e,maxy_e,l_ni,l_nj,F_nk,1-minx_e,1-miny_e,.false.,.false.,l_ni,0)
            call rpn_comm_xch_halo(py,minx_e,maxx_e,miny_e,maxy_e,l_ni,l_nj,F_nk,1-minx_e,1-miny_e,.false.,.false.,l_ni,0)
            call rpn_comm_xch_halo(pz,minx_e,maxx_e,miny_e,maxy_e,l_ni,l_nj,F_nk,1-minx_e,1-miny_e,.false.,.false.,l_ni,0)

            !Clipping to limit the upstream positions
            !----------------------------------------
            call adv_cliptraj_h (px,py,minx_e,maxx_e,miny_e,maxy_e,l_ni,l_nj,F_nk,i0_e,in_e,j0_e,jn_e,F_k0,'FLUX')

         end if

         if (l_west) then

            call adv_get_indices (ii_w_o, F_x, F_y, F_z, F_num , nind_w_o, &
                                  i0_w_o, in_w_o, j0_w_o, jn_w_o, F_k0 , F_nk, F_lev)

            if (.not.narrow_i_w_L) &
            call adv_get_indices (ii_w_i, F_x, F_y, F_z, F_num , nind_w_i, &
                                  i0_w_i, in_w_i, j0_w_i, jn_w_i, F_k0 , F_nk, F_lev)

         end if

         if (l_south) then

            call adv_get_indices (ii_s_o, F_x, F_y, F_z, F_num , nind_s_o, &
                                  i0_s_o, in_s_o, j0_s_o, jn_s_o, F_k0 , F_nk, F_lev)

            if (.not.narrow_j_s_L) &
            call adv_get_indices (ii_s_i, F_x, F_y, F_z, F_num , nind_s_i, &
                                  i0_s_i, in_s_i, j0_s_i, jn_s_i, F_k0 , F_nk, F_lev)

         end if

         if (l_east) then

            call adv_get_indices (ii_e_o, F_x, F_y, F_z, F_num , nind_e_o, &
                                  i0_e_o, in_e_o, j0_e_o, jn_e_o, F_k0 , F_nk, F_lev)

            if (.not.narrow_i_e_L) &
            call adv_get_indices (ii_e_i, F_x, F_y, F_z, F_num , nind_e_i, &
                                  i0_e_i, in_e_i, j0_e_i, jn_e_i, F_k0 , F_nk, F_lev)

         end if

         if (l_north) then

            call adv_get_indices (ii_n_o, F_x, F_y, F_z, F_num , nind_n_o, &
                                  i0_n_o, in_n_o, j0_n_o, jn_n_o, F_k0 , F_nk, F_lev)

            if (.not.narrow_j_n_L) &
            call adv_get_indices (ii_n_i, F_x, F_y, F_z, F_num , nind_n_i, &
                                  i0_n_i, in_n_i, j0_n_i, jn_n_i, F_k0 , F_nk, F_lev)

         end if

         Adz_pos_reset = 0

      end if  !END pos_reset
      !---------------------

      !--------------------------------------------
      !Bermejo-Conde LAM: Estimate FLUX_out/FLUX_in
      !--------------------------------------------
      F_cub_o = 0.0
      F_cub_i = 0.0

      i_h = 0.0

      nind_o = -1
      nind_i = -1

      if (l_west) then

         nind_o = nind_w_o
         nind_i = nind_w_i

         ii_o => ii_w_o
         ii_i => ii_w_i

         i0_i = i0_w_i
         in_i = in_w_i
         j0_i = j0_w_i
         jn_i = jn_w_i

         if (narrow_i_w_L) &
         call adv_BC_LAM_flux_loop_iH ()

      end if

      call adv_BC_LAM_flux_loop_o ()
      if (.not.narrow_i_w_L) &
      call adv_BC_LAM_flux_loop_i ()

      nind_o = -1
      nind_i = -1

      if (l_south) then

         nind_o = nind_s_o
         nind_i = nind_s_i

         ii_o => ii_s_o
         ii_i => ii_s_i

         i0_i = i0_s_i
         in_i = in_s_i
         j0_i = j0_s_i
         jn_i = jn_s_i

         if (narrow_j_s_L) &
         call adv_BC_LAM_flux_loop_iH ()

      end if

      call adv_BC_LAM_flux_loop_o ()
      if (.NOT.narrow_j_s_L) &
      call adv_BC_LAM_flux_loop_i ()

      nind_o = -1
      nind_i = -1

      if (l_east) then

         nind_o = nind_e_o
         nind_i = nind_e_i

         ii_o => ii_e_o
         ii_i => ii_e_i

         i0_i = i0_e_i
         in_i = in_e_i
         j0_i = j0_e_i
         jn_i = jn_e_i

         if (narrow_i_e_L) &
         call adv_BC_LAM_flux_loop_iH ()

      end if

      call adv_BC_LAM_flux_loop_o ()
      if (.NOT.narrow_i_e_L) &
      call adv_BC_LAM_flux_loop_i ()

      nind_o = -1
      nind_i = -1

      if (l_north) then

         nind_o = nind_n_o
         nind_i = nind_n_i

         ii_o => ii_n_o
         ii_i => ii_n_i

         i0_i = i0_n_i
         in_i = in_n_i
         j0_i = j0_n_i
         jn_i = jn_n_i

         if (narrow_j_n_L) &
         call adv_BC_LAM_flux_loop_iH ()

      end if

      call adv_BC_LAM_flux_loop_o ()
      if (.NOT.narrow_j_n_L) &
      call adv_BC_LAM_flux_loop_i ()

      if (l_west.and.     narrow_i_w_L.and.l_north.and..not.narrow_j_n_L) i_h(i0_w_i:in_w_i,j0_n_i:jn_n_i,F_k0:F_nk) = 0.
      if (l_east.and.     narrow_i_e_L.and.l_north.and..not.narrow_j_n_L) i_h(i0_e_i:in_e_i,j0_n_i:jn_n_i,F_k0:F_nk) = 0.
      if (l_west.and.     narrow_i_w_L.and.l_south.and..not.narrow_j_s_L) i_h(i0_w_i:in_w_i,j0_s_i:jn_s_i,F_k0:F_nk) = 0.
      if (l_east.and.     narrow_i_e_L.and.l_south.and..not.narrow_j_s_L) i_h(i0_e_i:in_e_i,j0_s_i:jn_s_i,F_k0:F_nk) = 0.

      if (l_west.and..not.narrow_i_w_L.and.l_north.and.     narrow_j_n_L) i_h(i0_w_i:in_w_i,j0_n_i:jn_n_i,F_k0:F_nk) = 0.
      if (l_east.and..not.narrow_i_e_L.and.l_north.and.     narrow_j_n_L) i_h(i0_e_i:in_e_i,j0_n_i:jn_n_i,F_k0:F_nk) = 0.
      if (l_west.and..not.narrow_i_w_L.and.l_south.and.     narrow_j_s_L) i_h(i0_w_i:in_w_i,j0_s_i:jn_s_i,F_k0:F_nk) = 0.
      if (l_east.and..not.narrow_i_e_L.and.l_south.and.     narrow_j_s_L) i_h(i0_e_i:in_e_i,j0_s_i:jn_s_i,F_k0:F_nk) = 0.

      if (l_west.and.     narrow_i_w_L.and.l_north.and.     narrow_j_n_L) i_h(l_ni+1:in_w_i,     1:jn_n_i,F_k0:F_nk) = 0.
      if (l_west.and.     narrow_i_w_L.and.l_north.and.     narrow_j_n_L) i_h(i0_w_i:  l_ni,j0_n_i:     0,F_k0:F_nk) = 0.

      if (l_east.and.     narrow_i_e_L.and.l_north.and.     narrow_j_n_L) i_h(i0_e_i:     0,     1:jn_n_i,F_k0:F_nk) = 0.
      if (l_east.and.     narrow_i_e_L.and.l_north.and.     narrow_j_n_L) i_h(     1:in_e_i,j0_n_i:     0,F_k0:F_nk) = 0.

      if (l_west.and.     narrow_i_w_L.and.l_south.and.     narrow_j_s_L) i_h(l_ni+1:in_w_i,j0_s_i:  l_nj,F_k0:F_nk) = 0.
      if (l_west.and.     narrow_i_w_L.and.l_south.and.     narrow_j_s_L) i_h(i0_w_i:  l_ni,l_nj+1:jn_s_i,F_k0:F_nk) = 0.

      if (l_east.and.     narrow_i_e_L.and.l_south.and.     narrow_j_s_L) i_h(i0_e_i:     0,j0_s_i:  l_nj,F_k0:F_nk) = 0.
      if (l_east.and.     narrow_i_e_L.and.l_south.and.     narrow_j_s_L) i_h(     1:in_e_i,l_nj+1:jn_s_i,F_k0:F_nk) = 0.

      if (num_PE_near_NW==Ptopo_myproc.and.narrow_i_w_L.and.narrow_j_n_L) i_h(            1:i_near_W-l_ni,j0_n_i:0,     F_k0:F_nk) = 0.

      if (num_PE_near_NE==Ptopo_myproc.and.narrow_i_e_L.and.narrow_j_n_L) i_h(i_near_E+l_ni:l_ni,         j0_n_i:0,     F_k0:F_nk) = 0.

      if (num_PE_near_SW==Ptopo_myproc.and.narrow_i_w_L.and.narrow_j_s_L) i_h(            1:i_near_W-l_ni,l_nj+1:jn_s_i,F_k0:F_nk) = 0.

      if (num_PE_near_SE==Ptopo_myproc.and.narrow_i_e_L.and.narrow_j_s_L) i_h(i_near_E+l_ni:l_ni,         l_nj+1:jn_s_i,F_k0:F_nk) = 0.

      if (narrow_L) then

         !Adjoint of Fill East-West/North-South Halo
         !------------------------------------------
         call rpn_comm_adj_halo (i_h,minx_e,maxx_e,miny_e,maxy_e,l_ni,l_nj,F_nk,1-minx_e,1-miny_e,.false.,.false.,l_ni,0)

         do k=1,F_nk
         do j=1,l_nj
         do i=1,l_ni
            n = (k-1)*l_ni*l_nj + (j-1)*l_ni + i
            F_cub_i(n) = i_h(i,j,k) + F_cub_i(n)
         end do
         end do
         end do

      end if

      return

contains

      real(kind=REAL64) function triprd(za,zb,zc,zd)
         real, intent(in) :: za
         real(kind=REAL64), intent(in) :: zb,zc,zd
         triprd = ((dble(za)-zb)*(za-zc)*(za-zd))
      end function triprd

      include 'adv_BC_LAM_flux_loop_o.inc'
      include 'adv_BC_LAM_flux_loop_i.inc'

      include 'adv_BC_LAM_flux_loop_iH.inc'

      end subroutine adv_BC_LAM_Aranami
