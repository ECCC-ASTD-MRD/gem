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
!
      subroutine adv_tracers (F_before_psadj_L)
      use adv
      use adv_options
      use adv_pos
      use gem_options
      use dyn_fisl_options
      use lam_options
      use glb_ld
      use gmm_itf_mod
      use gmm_tracers
      use gmm_vt0
      use HORgrid_options
      use tr3d
      use gem_timing
      implicit none
#include <arch_specific.hf>

      logical, intent(in) :: F_before_psadj_L
      logical qw_L,tr_not_before_psadj_L,tr_not_after_psadj_L,BC_activated_L
      integer n,jext,err,istat,mass_KEEP,i,j,k
      integer nind, nind_0, nind_3, nind_1, num
      integer i0,j0,in,jn,k0        ! scope of advection operations
      integer i0_0,j0_0,in_0,jn_0,i0_3,j0_3,in_3,jn_3,i0_1,j0_1,in_1,jn_1
      integer, pointer, contiguous, dimension(:) :: ii,ii_0,ii_3,ii_w !pre-computed index used in the tricubic lagrangian interp loop
      type(gmm_metadata) :: mymeta
      real, pointer, contiguous, dimension (:,:,:) :: fld_in, fld_out
      real, pointer, contiguous, dimension (:,:,:) :: px_w,py_w,pz_w
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk)   :: bidon
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk+1) :: pr_m,pr_t
      real, dimension(l_minx:l_maxx,l_miny:l_maxy)        :: pr_p0
!
!---------------------------------------------------------------------
!
      i0 = 1
      in = l_ni
      j0 = 1
      jn = l_nj

      jext=1
      if (Grd_yinyang_L) jext=2
      if (l_west)  i0 =        pil_w - jext
      if (l_south) j0 =        pil_s - jext
      if (l_east)  in = l_ni - pil_e + jext
      if (l_north) jn = l_nj - pil_n + jext

      i0_1 = i0 ; in_1 = in ; j0_1 = j0 ; jn_1 = jn

      k0=Lam_gbpil_t+1

      num    =  l_ni*l_nj*l_nk
      nind   = (in-i0+1)*(jn-j0+1)*(l_nk-k0+1)

      allocate (ii(4*nind))

      nind_1 = nind

      !Pre-compute indices ii used in: adv_tricub_lag3d_loop
      !-----------------------------------------------------
      call adv_get_indices(ii, pxt, pyt, pzt, num ,nind, &
                           i0, in, j0, jn, k0 , l_nk, 't')

      !Preparations Mass Conservation/Shape-Preserving
      !-----------------------------------------------
      if (.not.F_before_psadj_L.and.Adv_Mass_Cons_L) then

         !Pre-compute indices ii
         !----------------------
         call adv_get_ij0n_ext (i0_0,in_0,j0_0,jn_0,0) !CORE

         nind_0 = (in_0-i0_0+1)*(jn_0-j0_0+1)*(l_nk-k0+1)

         allocate (ii_0(4*nind_0))

         call adv_get_indices(ii_0, pxt, pyt, pzt, num ,nind_0, &
                              i0_0, in_0, j0_0, jn_0, k0 , l_nk, 't')

         !Bermejo-Conde LAM Flux ZLF
         !--------------------------
         if (Adv_BC_LAM_flux==2) then

            call adv_get_ij0n_ext (i0_3,in_3,j0_3,jn_3,3) !EXTENSION (BCS_BASE)

            nind_3 = (in_3-i0_3+1)*(jn_3-j0_3+1)*(l_nk-k0+1)

            allocate (ii_3(4*nind_3))

            call adv_get_indices(ii_3, pxt, pyt, pzt, num ,nind_3, &
                                 i0_3, in_3, j0_3, jn_3, k0 , l_nk, 't')

         end if

         !Reset AIR MASS at TIME 1 and TIME 0
         !-----------------------------------
         istat = gmm_get(gmmk_airm1_s,airm1)
         istat = gmm_get(gmmk_airm0_s,airm0)

         airm1 = 0. ; airm0 = 0.

         call get_density (bidon,airm1,1,l_minx,l_maxx,l_miny,l_maxy,l_nk,k0)
         call get_density (bidon,airm0,0,l_minx,l_maxx,l_miny,l_maxy,l_nk,k0)

         !Fill Halo of AIR MASS at TIME 0 (required in ILMC)
         !--------------------------------------------------
         call rpn_comm_xch_halo (airm0,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk, &
                                 G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

         !Reset pr_t(k)/pr_s at TIME 0 (used in Bermejo-Conde)
         !----------------------------------------------------
         istat = gmm_get(gmmk_pkps_s,pkps)
         istat = gmm_get(gmmk_st0_s, st0)

         !Obtain pressure levels at TIME 0
         !--------------------------------
         call calc_pressure (pr_m,pr_t,pr_p0,st0,l_minx,l_maxx,l_miny,l_maxy,l_nk)

         pkps = 0.

         do k=1,l_nk
         do j=1,l_nj
         do i=1,l_ni
            pkps(i,j,k) = pr_t(i,j,k)/pr_p0(i,j)
         end do
         end do
         end do

      end if

      Adv_component_S = 'INTP_TR'
      call gemtime_start (27, 'ADV_INTP_TR', 10)

      do n=1,Tr3d_ntr

         if (Tr3d_mono(n)>1.or.Tr3d_mass(n)/=0) then

         adv_intp_S = Tr3d_intp(n)

         qw_L= Tr3d_wload(n) .or. Tr3d_name_S(n)(1:2) == 'HU'

         if (qw_L .and. Schm_psadj==2 .and. (Tr3d_mono(n)>1 .or. Tr3d_mass(n)>0)) &
            call handle_error (-1,'ADV_TRACERS','CAUTION: Conservation TRACER (Water) + PSADJ DRY')

         if (F_before_psadj_L) then

            tr_not_before_psadj_L= .not. (qw_L .and. Schm_psadj==2)

            if (tr_not_before_psadj_L) cycle

         else

            tr_not_after_psadj_L= qw_L .and. Schm_psadj==2

            if (tr_not_after_psadj_L) cycle

         end if

         BC_activated_L = Tr3d_mass(n)==1.or.(Tr3d_mass(n)>=111.and.Tr3d_mass(n)<=139)

         if (BC_activated_L) then

            Adv_BC_LEGACY_L = Tr3d_mass(n)==1

            mass_KEEP = Tr3d_mass(n)

            if (.not.Adv_BC_LEGACY_L) then
               Adv_BC_weight  = (Tr3d_mass(n) - 100)/10
               Adv_BC_pexp_n  = (Tr3d_mass(n) - (100 + Adv_BC_weight*10))
               Tr3d_mass(n)   = 1
            else
               Adv_BC_weight  = 1
               Adv_BC_pexp_n  = 0
            end if

         end if

         err= gmm_get('TR/'//trim(Tr3d_name_S(n))//':P' ,fld_in ,mymeta)
         err= gmm_get('TR/'//trim(Tr3d_name_S(n))//':M' ,fld_out,mymeta)

         !No interpolation
         !----------------
         if (Tr3d_intp(n)=='NONE') then

            fld_out(i0:in,j0:jn,1:l_nk) = fld_in(i0:in,j0:jn,1:l_nk)

            cycle

         end if

         ii_w => ii ; i0 = i0_1 ; in = in_1 ; j0 = j0_1 ; jn = jn_1 ; nind = nind_1

         px_w => pxt ; py_w => pyt ; pz_w => pzt

         !Select pre-compute indices ii and positions for Mass Conservation/Shape-Preserving
         !----------------------------------------------------------------------------------
         if (Tr3d_mono(n)>1.or.Tr3d_mass(n)/=0) then

            ii_w => ii_0 ; i0 = i0_0 ; in = in_0 ; j0 = j0_0 ; jn = jn_0 ; nind = nind_0

            if (Tr3d_mass(n)==2.or.Tr3d_mass(n)==4) then !LCSL
               px_w => pxm_s ; py_w => pym_s ; pz_w => pzm_s
            end if

            !Bermejo-Conde LAM Flux ZLF
            !--------------------------
            if (Tr3d_mass(n)==1.and.Adv_BC_LAM_flux==2) then
               ii_w => ii_3 ; i0 = i0_3 ; in = in_3 ; j0 = j0_3 ; jn = jn_3 ; nind = nind_3
            end if

         end if

         call adv_cubic ('TR/'//trim(Tr3d_name_S(n))//':M', fld_out , fld_in, px_w, py_w, pz_w, &
                         l_ni, l_nj, l_nk, l_minx, l_maxx, l_miny, l_maxy, &
                         nind, ii_w, i0, in, j0, jn, k0,'t', Tr3d_mono(n), Tr3d_mass(n) )

         if (BC_activated_L) Tr3d_mass(n) = mass_KEEP

         end if

      end do
      call gemtime_stop (27)

      deallocate(ii)

      if (.not.F_before_psadj_L.and.Adv_Mass_Cons_L) deallocate(ii_0)
      if (.not.F_before_psadj_L.and.Adv_Mass_Cons_L.and.Adv_BC_LAM_flux==2) deallocate(ii_3)

!
!---------------------------------------------------------------------
!
      return
      end subroutine adv_tracers
