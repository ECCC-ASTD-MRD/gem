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

      subroutine adz_cubic ( F_dest, F_src, F_ni,F_nj,F_nk,&
                             F_minx,F_maxx,F_miny,F_maxy          ,&
                             F_i0, F_in, F_j0, F_jn               ,&
                             F_k0, F_lev_S, F_mono_L )
      use adz_mem
      use adz_options
      use adv_pos
      use gmm_itf_mod
      use gmm_tracers

      implicit none

      character(len=1), intent(in) :: F_lev_S
      logical, intent(in) :: F_mono_L
      integer, intent(in) :: F_ni,F_nj,F_nk
      integer, intent(in) :: F_minx,F_maxx,F_miny, F_maxy
      integer, intent(in) :: F_i0, F_j0, F_in, F_jn, F_k0
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in) :: F_src
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(out) :: F_dest

      integer :: num, nind, err, i0_e, in_e, j0_e, jn_e, k
      real, dimension(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,F_nk) :: extended,adv_o,adv_i
      real, dimension(F_ni,F_nj,F_nk) :: wrkc,w_mono,w_lin,w_min,w_max,w_cub_o,w_cub_i
!
!---------------------------------------------------------------------
!
      if (.not.Adz_Mass_Cons_tr_L.and.Adz_intp_S=="CUBIC") then
         call handle_error(-1,'ADZ_CUBIC','ADZ_CUBIC OBSOLETE')
      end if

      if (.not.Adz_Mass_Cons_tr_L.and.Adz_intp_S=="QUINTIC") then

         extended = 0.

         call rpn_comm_xch_halox( F_src, l_minx, l_maxx, l_miny, l_maxy,&
               l_ni, l_nj, l_nk, Adz_halox, Adz_haloy, .false., .false.,&
               extended, Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy, l_ni, 0)

         num = F_ni*F_nj*F_nk

         nind = (F_in-F_i0+1)*(F_jn-F_j0+1)*(F_nk-F_k0+1)

         if (F_mono_L) then
            call adv_qutver_lag3d_mono (wrkc,extended, &
                                        pxt,pyt,pzt,num,nind,ii_w,F_nk,F_lev_S)
         else
            call adv_qutver_lag3d      (wrkc,extended, &
                                        pxt,pyt,pzt,num,nind,ii_w,F_nk,F_lev_S)
         end if

         F_dest(F_i0:F_in,F_j0:F_jn,F_k0:F_nk) = &
           wrkc(F_i0:F_in,F_j0:F_jn,F_k0:F_nk)

      else if (Adz_Mass_Cons_tr_L) then

         extended = 0.

         if (Adz_BC_LAM_flux_n/=2) then
            call rpn_comm_xch_halox( F_src, l_minx, l_maxx, l_miny, l_maxy,&
                  l_ni, l_nj, l_nk, Adz_halox, Adz_haloy, .false., .false.,&
                  extended, Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy, l_ni, 0)
         end if

         num = F_ni*F_nj*F_nk

         nind = (F_in-F_i0+1)*(F_jn-F_j0+1)*(F_nk-F_k0+1)

         !Calculate FLUX_out/FLUX_in
         !--------------------------
         if (Adz_BC_LAM_flux_n>0) then

            !Bermejo-Conde LAM: Apply mask_o/mask_i to Tracer
            !------------------------------------------------
            if (Adz_BC_LAM_flux_n==1) then

               call adz_BC_LAM_mask (adv_o,adv_i,extended, &
                                     Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy,F_nk)

            !PSADJ LAM: Use mask_o/mask_i based on Tracer=1
            !----------------------------------------------
            else if (Adz_BC_LAM_flux_n==2) then

               adv_o = Adz_BC_LAM_mask_o
               adv_i = Adz_BC_LAM_mask_i

            end if

            !Estimate FLUX_out/FLUX_in based on Aranami et al. (2015)
            !--------------------------------------------------------
            call adv_BC_LAM_Aranami (w_cub_o,adv_o,w_cub_i,adv_i, &
                                     pxt,pyt,pzt,num,F_k0,F_nk,F_lev_S)

            err = gmm_get(gmmk_cub_o_s,fld_cub_o)
            err = gmm_get(gmmk_cub_i_s,fld_cub_i)

            call adz_get_ij0n_ext (i0_e,in_e,j0_e,jn_e,1) !EXTENSION (CFL)

            fld_cub_o = 0.
            fld_cub_i = 0.

            do k = F_k0, F_nk
               fld_cub_o(i0_e:in_e,j0_e:jn_e,k)= w_cub_o(i0_e:in_e,j0_e:jn_e,k)
               fld_cub_i(i0_e:in_e,j0_e:jn_e,k)= w_cub_i(i0_e:in_e,j0_e:jn_e,k)
            end do

            if (Adz_BC_LAM_flux_n == 2) return

         end if

         !High-order SL advection of Tracer + Storage Mono/Lin/Min/Max
         !------------------------------------------------------------
         if (Adz_intp_S=="CUBIC") then
            call adv_tricub_lag3d_conserv (wrkc,w_mono,w_lin,w_min,w_max,extended, &
                                           pxt,pyt,pzt,num,nind,ii_w,F_nk,F_lev_S)
         else if (Adz_intp_S=="QUINTIC") then
            call adv_qutver_lag3d_conserv (wrkc,w_mono,w_lin,w_min,w_max,extended, &
                                           pxt,pyt,pzt,num,nind,ii_w,F_nk,F_lev_S)
         end if

         F_dest(F_i0:F_in,F_j0:F_jn,F_k0:F_nk) = &
           wrkc(F_i0:F_in,F_j0:F_jn,F_k0:F_nk)

         err = gmm_get(gmmk_mono_s,fld_mono)
         err = gmm_get(gmmk_lin_s ,fld_lin )
         err = gmm_get(gmmk_min_s ,fld_min )
         err = gmm_get(gmmk_max_s ,fld_max )

         fld_mono(F_i0:F_in,F_j0:F_jn,F_k0:F_nk) = &
           w_mono(F_i0:F_in,F_j0:F_jn,F_k0:F_nk)
         fld_lin (F_i0:F_in,F_j0:F_jn,F_k0:F_nk) = &
           w_lin (F_i0:F_in,F_j0:F_jn,F_k0:F_nk)
         fld_min (F_i0:F_in,F_j0:F_jn,F_k0:F_nk) = &
           w_min (F_i0:F_in,F_j0:F_jn,F_k0:F_nk)
         fld_max (F_i0:F_in,F_j0:F_jn,F_k0:F_nk) = &
           w_max (F_i0:F_in,F_j0:F_jn,F_k0:F_nk)

      end if
!
!---------------------------------------------------------------------
!
      return
      end subroutine adz_cubic
