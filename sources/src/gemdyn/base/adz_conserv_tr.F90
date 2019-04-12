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

!**s/r adz_conserv_tr - Advection of Tracers + Mass Conservation/Shape-Preservation

      subroutine adz_conserv_tr ()

      use adz_mem
      use adz_options
      use adv_pos
      use gem_options
      use gem_timing
      use gmm_itf_mod
      use gmm_tracers
      use gmm_vt0
      use HORgrid_options
      use tr3d

      implicit none

#include <arch_specific.hf>

      !-----------------------------------------------------------------------------

      logical BC_activated_L
      integer n,err,mass_KEEP,i,j,k
      integer num,nind_0,nind_2
      integer i0_0,j0_0,in_0,jn_0,i0_2,j0_2,in_2,jn_2,i0_w,j0_w,in_w,jn_w
      integer, pointer, contiguous, dimension(:) :: ii_0,ii_2
      real, pointer, contiguous, dimension (:,:,:) :: src, dst
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk)   :: bidon,dst_w,store_pilot
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk+1) :: pr_m,pr_t
      real, dimension(l_minx:l_maxx,l_miny:l_maxy)        :: pr_p0
      logical :: verbose_L,BC_LAM_flux_1_L,BC_LAM_flux_2_L
      real, dimension(1,1,1), target :: empty

      !-----------------------------------------------------------------------------

      if (.not.Adz_Mass_Cons_L) return

      verbose_L = Adz_verbose/=0

      if (verbose_L) call stat_mass_tracers (1,"BEFORE ADVECTION")

      call gemtime_start (35, 'ADZ_INTTR_B', 10)

      err = gmm_get(gmmk_mono_s ,fld_mono)
      err = gmm_get(gmmk_lin_s  ,fld_lin)
      err = gmm_get(gmmk_min_s  ,fld_min)
      err = gmm_get(gmmk_max_s  ,fld_max)
      err = gmm_get(gmmk_cub_o_s,fld_cub_o)
      err = gmm_get(gmmk_cub_i_s,fld_cub_i)

      !-----------------------------------------------
      !Pre-Compute indices to be used in interpolation
      !-----------------------------------------------
      Adz_pos_reset = 1

      call adz_get_ij0n_ext (i0_0,in_0,j0_0,jn_0,0) !CORE

      nind_0 = (in_0-i0_0+1)*(jn_0-j0_0+1)*(l_nk-Adz_k0+1)

      allocate (ii_0(4*nind_0))

      call adv_get_indices (ii_0,pxt,pyt,pzt,num,nind_0, &
                            i0_0,in_0,j0_0,jn_0,Adz_k0,l_nk,'t')

      if (Adz_BC_LAM_flux==2) then

         call adz_get_ij0n_ext (i0_2,in_2,j0_2,jn_2,2) !EXTENSION (BCS_BASE)

         nind_2 = (in_2-i0_2+1)*(jn_2-j0_2+1)*(l_nk-Adz_k0+1)

         allocate (ii_2(4*nind_2))

         call adv_get_indices(ii_2,pxt,pyt,pzt,num,nind_2, &
                              i0_2,in_2,j0_2,jn_2,Adz_k0,l_nk,'t')

      end if

      !-----------------------------------
      !Reset Air Mass at TIME 1 and TIME 0
      !-----------------------------------
      err = gmm_get(gmmk_airm1_s,airm1)
      err = gmm_get(gmmk_airm0_s,airm0)

      airm1 = 0. ; airm0 = 0.

      call get_density (bidon,airm1,1,l_minx,l_maxx,l_miny,l_maxy,l_nk,Adz_k0)
      call get_density (bidon,airm0,0,l_minx,l_maxx,l_miny,l_maxy,l_nk,Adz_k0)

      !Fill Halo of Air Mass at TIME 0 (required in ILMC)
      !--------------------------------------------------
      call rpn_comm_xch_halo (airm0,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

      !----------------------------------------------------
      !Reset pr_t(k)/pr_s at TIME 0 (used in Bermejo-Conde)
      !----------------------------------------------------
      err = gmm_get(gmmk_pkps_s,pkps)
      err = gmm_get(gmmk_st0_s, st0)

      call calc_pressure (pr_m,pr_t,pr_p0,st0,l_minx,l_maxx,l_miny,l_maxy,l_nk)

      pkps = 0.

      do k=1,l_nk
         do j=1,l_nj
            do i=1,l_ni
               pkps(i,j,k) = pr_t(i,j,k) / pr_p0(i,j)
            end do
         end do
      end do

      !----------------------------------------------------------------
      !Advection of Tracers + Global Mass Conservation/Shape-Preserving
      !----------------------------------------------------------------
      do n=1,Tr3d_ntr

         if ( (Tr3d_mass(n) == 0) .and. (Tr3d_mono(n) < 2) .and. (Tr3d_intp(n) == 'CUBIC') ) cycle

         Adz_Mass_Cons_tr_L = .not. ( (Tr3d_mass(n) == 0) .and. (Tr3d_mono(n) < 2) )

         Adz_intp_S = Tr3d_intp(n)

         BC_LAM_flux_1_L = Tr3d_mass(n)==1.and..not.Grd_yinyang_L.and.Adz_BC_LAM_flux==1
         BC_LAM_flux_2_L = Tr3d_mass(n)==1.and..not.Grd_yinyang_L.and.Adz_BC_LAM_flux==2

         Adz_BC_LAM_flux_n = 0
         if (BC_LAM_flux_1_L) Adz_BC_LAM_flux_n = 1

         !Initialize Bermejo-Conde parameters based on Tr3d_mass
         !------------------------------------------------------
         BC_activated_L = Tr3d_mass(n)==1.or.(Tr3d_mass(n)>=111.and.Tr3d_mass(n)<=139)

         if (BC_activated_L) then

            Adz_BC_LEGACY_L = Tr3d_mass(n)==1

            mass_KEEP = Tr3d_mass(n)

            if (.not.Adz_BC_LEGACY_L) then
               Adz_BC_weight  = (Tr3d_mass(n) - 100)/10
               Adz_BC_pexp_n  = (Tr3d_mass(n) - (100 + Adz_BC_weight*10))
               Tr3d_mass(n)   = 1
            else
               Adz_BC_weight  = 1
               Adz_BC_pexp_n  = 0
            end if

         end if

         err= gmm_get('TR/'//trim(Tr3d_name_S(n))//':P' ,src)
         err= gmm_get('TR/'//trim(Tr3d_name_S(n))//':M' ,dst)

         !Bermejo-Conde LAM Flux ZLF: Set ZERO piloting conditions outside EXTENSION (CFL)
         !--------------------------------------------------------------------------------
         if (BC_LAM_flux_2_L) call adz_BC_LAM_zlf (src,empty,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk,1)

         if (Tr3d_mass(n)==1.and.Adz_BC_LAM_flux==2) then
            ii_w => ii_2 ; i0_w = i0_2 ; in_w = in_2 ; j0_w = j0_2 ; jn_w = jn_2
         else
            ii_w => ii_0 ; i0_w = i0_0 ; in_w = in_0 ; j0_w = j0_0 ; jn_w = jn_0
         end if

         dst_w(:,:,:) = dst(:,:,:) !Copy piloting conditions

         !High-order interpolation + storage of MONO/LIN/MAX/MIN +
         !estimate FLUX_out/FLUX_in if Bermejo-Conde LAM Flux Aranami
         !-----------------------------------------------------------
         call adz_cubic (dst_w, src,                                &
                         l_ni,l_nj,l_nk,l_minx,l_maxx,l_miny,l_maxy,&
                         i0_w, in_w, j0_w, jn_w, Adz_k0,'t',.false.)

         !Bermejo-Conde LAM Flux ZLF
         !--------------------------
         if (BC_LAM_flux_2_L) then

             !Keep piloting conditions outside CORE
             !-------------------------------------
             call adz_BC_LAM_zlf (dst,store_pilot,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk,2)

             !ZERO piloting conditions outside EXTENSION (BCS_BASE)
             !-----------------------------------------------------
             call adz_BC_LAM_zlf (dst_w,empty,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk,3)

         end if

         !Apply a posteriori Mass-fixer/Shape-Preserving schemes
         !------------------------------------------------------
         call adz_a_posteriori ('TR/'//trim(Tr3d_name_S(n))//':M',dst,dst_w,fld_mono,fld_lin,fld_min,fld_max, &
                                 src, fld_cub_o, fld_cub_i, l_minx, l_maxx, l_miny, l_maxy, l_nk, &
                                 i0_w, in_w, j0_w, jn_w, Adz_k0, Tr3d_mono(n), Tr3d_mass(n))

         !Bermejo-Conde LAM Flux ZLF: Reset piloting conditions outside CORE
         !------------------------------------------------------------------
         if (BC_LAM_flux_2_L) call adz_BC_LAM_zlf (dst,store_pilot,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk,4)

         if (BC_activated_L) Tr3d_mass(n) = mass_KEEP

         !Reset to DEFAULT
         !----------------
         Adz_Mass_Cons_tr_L = .false.
         Adz_BC_LAM_flux_n = 0
         Adz_pos_reset = 0
         Adz_intp_S = 'CUBIC'

      end do

      deallocate(ii_0)

      if (Adz_BC_LAM_flux==2) deallocate(ii_2)

      call gemtime_stop (35)

      if (verbose_L) call stat_mass_tracers (0,"AFTER ADVECTION")

      return

      end subroutine adz_conserv_tr
