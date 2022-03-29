!i---------------------------------- LICENCE BEGIN -------------------------------
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

!**s/r hzd_smago_in_split - Applies horizontal Smagorinsky-type nonlinear diffusion
!                           in the split mode
!
      subroutine hzd_smago_in_split (F_u, F_v, F_w, F_t, F_zd, &
                               lminx, lmaxx, lminy, lmaxy, &
                               nk, smago_momentum_L)
      use cstv
      use dcst
      use gem_options
      use geomh
      use gmm_itf_mod
      use gmm_smag
      use gmm_pw
      use glb_ld
      use glb_pil
      use hvdif_options
      use HORgrid_options
      use hzd_mod
      use tdpack
      use tr3d
      use ver
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: lminx,lmaxx,lminy,lmaxy, nk
      real, dimension(lminx:lmaxx,lminy:lmaxy,nk), intent(inout) :: F_u, F_v, F_w, F_t,F_zd
!
!Author:  Claude Girard and Syed Husain
!
      integer :: i, j, k, istat, i0, in, j0, jn
      real, dimension(lminx:lmaxx,lminy:lmaxy) :: tension, shear_z, kt, kz
      real, dimension(lminx:lmaxx,lminy:lmaxy) :: smagcoef_z, smagcoef_u, smagcoef_v
      real, dimension(lminx:lmaxx,lminy:lmaxy) :: smagcoef_uo,smagcoef_vo, F_du,F_dv
      real, dimension(lminx:lmaxx,lminy:lmaxy) :: tension_u, shear_u, tension_v, shear_v
      real, dimension(lminx:lmaxx,lminy:lmaxy,G_nk) :: th
      real, dimension(lminx:lmaxx,lminy:lmaxy) :: hutmp, pres_t
      real, pointer, dimension (:,:,:) :: hu
      real, dimension(nk) :: base_coefM, base_coefT
      real :: cdelta2, tension_z, shear
      real :: fact, smagparam, ismagprandtl, ismagprandtl_hu
      real :: crit_coef
      logical :: switch_on_THETA, switch_on_hu, switch_on_fric_heat
      logical :: switch_on_W, smago_momentum_L, switch_on_wzd

      real, parameter :: p_naught=100000., eps=1.0e-5

      if( (hzd_smago_param <= 0.) .and. (hzd_smago_lnr(2) <=0.) ) return
      switch_on_wzd   = (Hzd_lnr <= 0.)

      smagparam= hzd_smago_param
      switch_on_THETA = Hzd_smago_prandtl > 0. .and. &
                        (.not. smago_momentum_L) .and. (Hzd_lnr_theta <= 0.)
      switch_on_hu    = (Hzd_smago_prandtl_hu > 0.) .and. &
                        (.not. smago_momentum_L) .and. (Hzd_lnr_tr <= 0.)
      switch_on_fric_heat = (Hzd_smago_fric_heat > 0.) .and. (.not. smago_momentum_L)
      switch_on_W = (.not. smago_momentum_L) .and. (switch_on_wzd)


      if (Grd_yinyang_L) then
         call yyg_xchng_vec_uv2uv (F_u,F_v,l_minx,l_maxx,l_miny,l_maxy,G_nk)
         call yyg_xchng (F_zd, l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj, G_nk, &
                         .false., 'CUBIC', .true.)
         call rpn_comm_xch_halo( F_u, l_minx,l_maxx,l_miny,l_maxy,l_niu,l_nj, G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
         call rpn_comm_xch_halo( F_v, l_minx,l_maxx,l_miny,l_maxy,l_ni ,l_njv, G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      else
         call rpn_comm_xch_halo( F_u, l_minx,l_maxx,l_miny,l_maxy,l_niu,l_nj, G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
         call rpn_comm_xch_halo( F_v, l_minx,l_maxx,l_miny,l_maxy,l_ni ,l_njv, G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
         call rpn_comm_xch_halo( F_zd, l_minx,l_maxx,l_miny,l_maxy,l_ni ,l_nj, G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      end if
         if (switch_on_theta) then
            istat = gmm_get(gmmk_pw_pt_plus_s, pw_pt_plus)
            if(Grd_yinyang_L)then
               call yyg_xchng (F_t, l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj, G_nk, &
                               .false., 'CUBIC', .true.)
               call yyg_xchng (pw_pt_plus, l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj, G_nk+1, &
                               .false., 'CUBIC', .true.)

            else
               call rpn_comm_xch_halo( F_t, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj, G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
               call rpn_comm_xch_halo( pw_pt_plus, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj, G_nk+1, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )

            end if
            ismagprandtl    = 1./Hzd_smago_prandtl
         end if

         if (switch_on_W) then
            if(Grd_yinyang_L) then
               call yyg_xchng (F_w, l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj, G_nk, &
                               .false., 'CUBIC', .true.)
            else
               call rpn_comm_xch_halo( F_w, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj, G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
            end if
         end if

         if (switch_on_hu) then
            istat = gmm_get('TR/'//trim(Tr3d_name_S(1))//':P' ,hu)
            if (Grd_yinyang_L)then

               call yyg_xchng (hu, l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj, G_nk, &
                               .false., 'CUBIC', .true.)
            else
               call rpn_comm_xch_halo( hu, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj, G_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
            end if
            ismagprandtl_hu = 1./Hzd_smago_prandtl_hu
         end if

      cdelta2 = (smagparam * Dcst_rayt_8 * geomh_hy_8)**2
      crit_coef = 0.25*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8

      i0  = 1    + pil_w
      in  = l_ni - pil_e

      j0  = 1    + pil_s
      jn  = l_nj - pil_n

      istat = gmm_get (gmmk_smag_s, smag)
      do k=1,nk
         base_coefM(k)=Hzd_smago_lnrM_8(k)*crit_coef
         base_coefT(k)=Hzd_smago_lnrT_8(k)*crit_coef
      end do

      do k=1,nk

         do j=j0-2, jn+2
            do i=i0-2, in+2

               tension(i,j) = ((F_u(i,j,k) - F_u(i-1,j,k)) * geomh_invDX_8(j)) &
                            - ((F_v(i,j,k) * geomh_invcyv_8(j) - F_v(i,j-1,k) * geomh_invcyv_8(j-1)) &
                                 * geomh_invDY_8 * geomh_cy_8(j))

               shear_z(i,j) = ((F_v(i+1,j,k) - F_v(i,j,k)) * geomh_invDXv_8(j)) &
                            + ((F_u(i,j+1,k) * geomh_invcy_8(j+1) - F_u(i,j,k) * geomh_invcy_8(j)) &
                                 * geomh_invDY_8 * geomh_cyv_8(j))
            end do
         end do

         if(hzd_smago_param <= 0. ) then
            do j=j0-1, jn+1
               do i=i0-1, in+1

                  smagcoef_u(i,j) = base_coefT(k) * geomh_invDX_8(j)
                  smagcoef_v(i,j) = base_coefT(k) * geomh_cyv_8(j) * geomh_invDY_8

                  kt(i,j) = geomh_cy2_8(j)  * base_coefM(k) * tension(i,j)
                  kz(i,j) = geomh_cyv2_8(j) * base_coefM(k) * shear_z(i,j)

                  smagcoef_uo(i,j)= base_coefM(k)
                  smagcoef_vo(i,j)= base_coefM(k)

               end do
            end do
         else
            do j=j0-1, jn+1
               do i=i0-1, in+1

                  tension_z = 0.25d0 * (tension(i,j) + tension(i+1,j) + tension(i,j+1) + tension(i+1,j+1))
                  shear     = 0.25d0 * (shear_z(i,j) + shear_z(i-1,j) + shear_z(i,j-1) + shear_z(i-1,j-1))
                  tension_u(i,j) = 0.5d0 * (tension(i,j) + tension(i+1,j))
                  shear_u(i,j)   = 0.5d0 * (shear_z(i,j) + shear_z(i,j-1))
                  tension_v(i,j) = 0.5d0 * (tension(i,j+1) + tension(i,j))
                  shear_v(i,j)   = 0.5d0 * (shear_z(i,j) + shear_z(i-1,j))

                  smag(i,j,k) = min((cdelta2 * sqrt(tension(i,j)**2 + shear**2)) + base_coefM(k), crit_coef)
                  smagcoef_z(i,j) = min((cdelta2 * sqrt(tension_z**2 + shear_z(i,j)**2)) + base_coefM(k), crit_coef)
                  smagcoef_u(i,j) = min((cdelta2 * sqrt(tension_u(i,j)**2 + shear_u(i,j)**2)) + base_coefT(k), crit_coef)
                  smagcoef_v(i,j) = min((cdelta2 * sqrt(tension_v(i,j)**2 + shear_v(i,j)**2)) + base_coefT(k), crit_coef)

                  smagcoef_uo(i,j)= smagcoef_u(i,j)
                  smagcoef_vo(i,j)= smagcoef_v(i,j)

                  smagcoef_u(i,j) = smagcoef_u(i,j) * geomh_invDX_8(j)
                  smagcoef_v(i,j) = smagcoef_v(i,j) * geomh_cyv_8(j) * geomh_invDY_8

                  kt(i,j) = geomh_cy2_8(j)  * smag(i,j,k) * tension(i,j)
                  kz(i,j) = geomh_cyv2_8(j) * smagcoef_z(i,j) * shear_z(i,j)

               end do
            end do
         end if

         fact=Cstv_dt_8
         do j=j0, jn
            do i=i0, in

               F_du(i,j) = F_u(i,j,k) + fact * smagcoef_uo(i,j) *      (       &
                              ( (F_u(i+1,j,k)-F_u(i,j,k))*geomh_invDX_8(j)  -  &
                                (F_u(i,j,k)-F_u(i-1,j,k))*geomh_invDX_8(j)) *  &
                                geomh_invDXv_8(j)       &
                            + ( (F_u(i,j+1,k)-F_u(i,j,k))*geomh_invDY_8     -  &
                                (F_u(i,j,k)-F_u(i,j-1,k))*geomh_invDY_8)    *  &
                                geomh_invDY_8           &
                            - ( (F_v(i+1,j,k)-F_v(i,j,k))*geomh_invDXv_8(j) +  &
                                (F_v(i,j,k)-F_v(i-1,j,k))*geomh_invDXv_8(j))*  &
                                geomh_tyoa_8(j)         &
                            - ( (F_u(i,j+1,k)-F_u(i,j,k))*geomh_invDY_8     +  &
                                (F_u(i,j,k)-F_u(i,j-1,k))*geomh_invDY_8)    *  &
                                geomh_tyoa_8(j)*0.5d0   &
                            + ( (2.0d0 - geomh_cy2_8(j))*F_u(i,j,k)*           &
                                Dcst_inv_rayt_8**2) ) + fact *         (       &
                            + ( (smagcoef_uo(i+1,j)-smagcoef_uo(i,j))*         &
                                 geomh_invDX_8(j) +                            &
                                (smagcoef_uo(i,j)-smagcoef_uo(i-1,j))*         &
                                 geomh_invDX_8(j))*0.05d0*tension_u(i,j)       &
                            + ( (smagcoef_uo(i,j+1)-smagcoef_uo(i,j))*         &
                                 geomh_invDY_8 +                               &
                                (smagcoef_uo(i,j)-smagcoef_uo(i,j-1))*         &
                                 geomh_invDY_8)*0.05d0*shear_u(i,j)    )


               F_dv(i,j) = F_v(i,j,k) + fact * smagcoef_vo(i,j) *      (          &
                              ( (F_v(i+1,j,k)-F_v(i,j,k))*geomh_invDXv_8(j)     - &
                                (F_v(i,j,k)-F_v(i-1,j,k))*geomh_invDXv_8(j))    * &
                                geomh_invDX_8(j)        &
                            + ( (F_v(i,j+1,k)-F_v(i,j,k))*geomh_invDY_8  -        &
                                (F_v(i,j,k)-F_v(i,j-1,k))*geomh_invDY_8)   *      &
                                geomh_invDY_8           &
                            + ( (F_u(i,j+1,k)-F_u(i-1,j+1,k))*geomh_invDX_8(j+1)+ &
                                (F_u(i,j,k)-F_u(i-1,j,k))*geomh_invDX_8(j))     * &
                                geomh_tyoav_8(j)         &
                            - ( (F_v(i,j+1,k)-F_v(i,j,k))*geomh_invDY_8  +        &
                                (F_v(i,j,k)-F_v(i,j-1,k))*geomh_invDY_8)   *      &
                                geomh_tyoav_8(j)*0.5d0   &
                            + ( (2.0d0 - geomh_cyv2_8(j))*F_v(i,j,k)*             &
                                Dcst_inv_rayt_8**2))  + fact *         (          &
                            + ( (smagcoef_vo(i+1,j)-smagcoef_vo(i,j))*            &
                                 geomh_invDXv_8(j) +                              &
                                (smagcoef_vo(i,j)-smagcoef_vo(i-1,j))*            &
                                 geomh_invDXv_8(j))*0.05d0*shear_v(i,j)           &
                            - ( (smagcoef_vo(i,j+1)-smagcoef_vo(i,j))*            &
                                 geomh_invDY_8 +                                  &
                                (smagcoef_vo(i,j)-smagcoef_vo(i,j-1))*            &
                                 geomh_invDY_8)*0.05d0*tension_v(i,j)  )

               F_u(i,j,k)=F_du(i,j)
               F_v(i,j,k)=F_dv(i,j)

            end do
         end do

         if (switch_on_wzd) then
            do j=j0, jn
               do i=i0, in
                  F_zd(i,j,k) = F_zd(i,j,k) + fact * ( &
                       geomh_invDXMu_8(j)*((smagcoef_u(i,j)   * (F_zd(i+1,j,k) - F_zd(i,j,k))) -   &
                                           (smagcoef_u(i-1,j) * (F_zd(i,j,k  ) - F_zd(i-1,j,k))) ) &
                     + geomh_invcy_8(j)*geomh_invDYMv_8(j)                                         &
                                         *((smagcoef_v(i,j)   * (F_zd(i,j+1,k) - F_zd(i,j,k))) -   &
                                           (smagcoef_v(i,j-1) * (F_zd(i,j,k  ) - F_zd(i,j-1,k))) ) )
               end do
            end do
         end if

         if (switch_on_W) then
            do j=j0, jn
               do i=i0, in
                  F_w(i,j,k) = F_w(i,j,k) + fact * ( &
                       geomh_invDXMu_8(j)*((smagcoef_u(i,j)   * (F_w(i+1,j,k) - F_w(i,j,k))) -   &
                                           (smagcoef_u(i-1,j) * (F_w(i  ,j,k) - F_w(i-1,j,k))) ) &
                     + geomh_invcy_8(j)*geomh_invDYMv_8(j)                                       &
                                         *((smagcoef_v(i,j)   * (F_w(i,j+1,k) - F_w(i,j,k))) -   &
                                           (smagcoef_v(i,j-1) * (F_w(i,j  ,k) - F_w(i,j-1,k))) ) )
              end do
           end do
         end if

         if (switch_on_THETA) then
            do j=j0-1, jn+1
               do i=i0-1, in+1
                  pres_t(i,j) = (p_naught/pw_pt_plus(i,j,k))**cappa_8
                  th(i,j,k)=F_t(i,j,k)*pres_t(i,j)
               end do
            end do

            fact=Cstv_dt_8*ismagprandtl
            do j=j0, jn
               do i=i0, in
                  th(i,j,k) = th(i,j,k) + fact *  ( &
                       geomh_invDXMu_8(j)*((smagcoef_u(i,j)   * (th(i+1,j,k) - th(i,j,k))) - &
                                           (smagcoef_u(i-1,j) * (th(i,j,k) - th(i-1,j,k))) ) &
                     + geomh_invcy_8(j)*geomh_invDYMv_8(j) &
                                         *((smagcoef_v(i,j)   * (th(i,j+1,k) - th(i,j,k))) - &
                                           (smagcoef_v(i,j-1) * (th(i,j,k) - th(i,j-1,k))) ) )

                  F_t(i,j,k) = th(i,j,k)/pres_t(i,j)
               end do
            end do
         end if

         if (switch_on_hu) then
            fact=Cstv_dt_8*ismagprandtl_hu
            do j=j0, jn
               do i=i0, in
                  hutmp(i,j)=fact * ( &
                          geomh_invDXMu_8(j)*((smagcoef_u(i,j)   * (hu(i+1,j,k) - hu(i,j,k))) - &
                                        (smagcoef_u(i-1,j) * (hu(i,j,k) - hu(i-1,j,k))) ) &
                        + geomh_invcy_8(j)*geomh_invDYMv_8(j) &
                                      *((smagcoef_v(i,j)   * (hu(i,j+1,k) - hu(i,j,k))) - &
                                        (smagcoef_v(i,j-1) * (hu(i,j,k) - hu(i,j-1,k))) ) )
                  end do
               end do
            do j=j0, jn
               do i=i0, in
                  hu(i,j,k)=hu(i,j,k)+hutmp(i,j)
               end do
            end do
         end if

      end do

      if (switch_on_fric_heat) then
         fact=Cstv_dt_8/cdelta2**2/cpd_8

         do k=1, nk
            if(k /= nk) then
               do j=j0, jn
                  do i=i0, in
                     F_t(i,j,k)=F_t(i,j,k) &
                               +fact*((0.5d0*(smag(i,j,k+1)+smag(i,j,k)))**3)
                  end do
               end do
            else
               do j=j0, jn
                  do i=i0, in
                     F_t(i,j,k)=F_t(i,j,k) + fact*smag(i,j,k)**3
                  end do
               end do
            end if
         end do

      end if

      return
      end subroutine hzd_smago_in_split
