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

!**s/r hzd_main_stag - Main controler for horizontal diffusion

      subroutine hzd_main_stag
      use hzd_ctrl
      use hzd_exp
      use gmm_vt1
      use HORgrid_options
      use hvdif_options
      use ens_options
      use glb_ld
      use lun
      use tr3d
      use gmm_itf_mod
      use gem_timing
      implicit none
#include <arch_specific.hf>

      character(len=GMM_MAXNAMELENGTH) :: tr_name
      real, pointer, dimension(:,:,:) :: tr1
      logical switch_on_UVW, switch_on_TR, switch_on_vrtspng_UVT    , &
              switch_on_vrtspng_W, switch_on_eqspng, switch_on_THETA
      logical xch_UV,xch_TT,xch_TR,xch_WZD
      integer i,j,k,istat,n
      integer, save :: depth
      real*8 pis2,fract
      real, dimension(:)    , pointer, save :: weight=> null()
      real, dimension(:,:,:), pointer       :: tr
!
!-------------------------------------------------------------------
!
      if (Lun_debug_L) write (Lun_out,1000)

      call gemtime_start ( 60, 'HZD_main', 1 )
      xch_UV = .false.
      xch_TT = .false.
      xch_TR = .false.
      xch_WZD= .false.
      if (ens_conf) then
          xch_UV = .true.
          xch_TT = .true.
      end if
      switch_on_UVW         = Hzd_lnr       > 0.
      switch_on_TR          =(Hzd_lnr_tr    > 0.).and.any(Tr3d_hzd)
      switch_on_THETA       = Hzd_lnr_theta > 0.
      switch_on_vrtspng_UVT = Vspng_nk      >=1
      switch_on_vrtspng_W   = Vspng_nk      >=1
      switch_on_eqspng      = Eq_nlev       > 1

      istat = gmm_get(gmmk_ut1_s,ut1)
      istat = gmm_get(gmmk_vt1_s,vt1)
      istat = gmm_get(gmmk_zdt1_s,zdt1)
      istat = gmm_get(gmmk_tt1_s,tt1)
      istat = gmm_get(gmmk_wt1_s,wt1)

      call itf_ens_hzd (ut1,vt1,tt1, l_minx,l_maxx,l_miny,l_maxy, G_nk)

!**********************************
!  Horizontal diffusion on theta  *
!**********************************

      if ( switch_on_THETA ) then
         call gemtime_start ( 61, 'HZD_theta', 60 )
         xch_TT = .true.
         call pw_update_GPW
         call hzd_theta
         call gemtime_stop  ( 61 )
      end if

!**********************************
!  Horizontal diffusion on tracers*
!**********************************

      if ( switch_on_TR ) then
         call gemtime_start ( 62, 'HZD_tracers', 60 )
         xch_TR = .true.
         do i=1, Tr3d_ntr
            if (Tr3d_hzd(i)) then
               nullify (tr)
               istat = gmm_get('TR/'//trim(Tr3d_name_S(i))//':P',tr)
               if (istat == 0) call hzd_ctrl4 &
                     (tr, 'S_TR', l_minx,l_maxx,l_miny,l_maxy,G_nk)
            end if
         end do
         call gemtime_stop  ( 62 )
      end if

!************************
!  Horizontal diffusion *
!************************

      if ( switch_on_UVW ) then
         call gemtime_start ( 64, 'HZD_bkgrnd', 60 )
         xch_UV = .true.
         xch_TT = .true.
         xch_WZD= .true.
         call hzd_ctrl4 ( ut1, vt1, l_minx,l_maxx,l_miny,l_maxy,G_nk)
         call hzd_ctrl4 (zdt1, 'S', l_minx,l_maxx,l_miny,l_maxy,G_nk)
         call hzd_ctrl4 ( wt1, 'S', l_minx,l_maxx,l_miny,l_maxy,G_nk)
         call gemtime_stop ( 64 )
      end if

!********************
!  Vertical sponge  *
!********************

      if ( switch_on_vrtspng_UVT ) then
         call gemtime_start ( 65, 'V_SPNG', 60 )
         xch_UV = .true.
         xch_TT = .true.
         call hzd_exp_deln ( ut1,  'U', l_minx,l_maxx,l_miny,l_maxy,&
                             Vspng_nk, F_VV=vt1, F_type_S='VSPNG' )
         call hzd_exp_deln ( tt1, 'M', l_minx,l_maxx,l_miny,l_maxy,&
                             Vspng_nk, F_type_S='VSPNG' )
         call gemtime_stop ( 65 )
      end if
      if ( switch_on_vrtspng_W ) then
         call gemtime_start ( 65, 'V_SPNG', 60 )
         xch_WZD= .true.
         if (Vspng_riley_L) then
            if (.not. associated(weight)) then
               pis2 = acos(0.0d0)
               depth= Vspng_nk+1
               allocate (weight(depth-1))
               weight(1)= 0.
               do k=2, depth-1
                  fract= (dble(k) - dble(depth))/ dble(depth-1)
                  weight(k)= cos(pis2*fract)**2
               end do
            end if
            do k=1, depth-1
               do j=1, l_nj
                  do i=1, l_ni
                     zdt1(i,j,k) = zdt1(i,j,k)*weight(k)
                      wt1(i,j,k) =  wt1(i,j,k)*weight(k)
                  end do
               end do
            end do
         else
            call hzd_exp_deln ( zdt1, 'M', l_minx,l_maxx,l_miny,l_maxy,&
                                Vspng_nk, F_type_S='VSPNG' )
            call hzd_exp_deln ( wt1, 'M', l_minx,l_maxx,l_miny,l_maxy,&
                                Vspng_nk, F_type_S='VSPNG' )
         end if
         call gemtime_stop ( 65 )
      end if

!**********************
!  Equatorial sponge  *
!**********************

      if ( switch_on_eqspng ) then
         call gemtime_start ( 67, 'EQUA_SPNG', 60)
         xch_UV= .true.
         call eqspng_drv (ut1,vt1,l_minx,l_maxx,l_miny,l_maxy,G_nk)
         call gemtime_stop ( 67 )
      end if

      call itf_ens_hzd ( ut1,vt1,tt1, l_minx,l_maxx,l_miny,l_maxy, G_nk )

!*********************************************************
!  Yin-Yang exchange pilot zones, blend wind overlap zones before physics*
!*********************************************************
      if (Grd_yinyang_L) then
         if (xch_UV) call yyg_xchng_vec_uv2uv (ut1,vt1,l_minx,l_maxx,l_miny,l_maxy,G_nk)
         if (xch_TT) call yyg_xchng (tt1, l_minx,l_maxx,l_miny,l_maxy, l_ni, l_nj,G_nk,&
                                     .false., 'CUBIC', .false.)
         if (xch_WZD) then
            call yyg_xchng (zdt1, l_minx,l_maxx,l_miny,l_maxy, l_ni, l_nj, G_nk,&
                            .false., 'CUBIC', .false.)
            call yyg_xchng (wt1 , l_minx,l_maxx,l_miny,l_maxy, l_ni, l_nj, G_nk,&
                            .false., 'CUBIC', .false. )
         end if
         if (xch_TR) then
            do n= 1, Tr3d_ntr
               tr_name = 'TR/'//trim(Tr3d_name_S(n))//':P'
               istat = gmm_get(tr_name, tr1)
               call yyg_xchng (tr1 , l_minx,l_maxx,l_miny,l_maxy,l_ni, l_nj, G_nk,&
                               .true., 'CUBIC', .false.)
            end do
          end if

      end if

      call gemtime_stop ( 60 )

 1000 format(3X,'MAIN HORIZONTAL DIFFUSION : (S/R HZD_MAIN)')
!
!-------------------------------------------------------------------
!
      return
      end
