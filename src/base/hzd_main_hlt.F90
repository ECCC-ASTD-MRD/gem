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

!**s/r hzd_main - Main controler for horizontal diffusion

      subroutine hzd_main_hlt
      use hzd_exp_hlt
      use glb_ld
      use gmm_vt1
      use gmm_pw
      use gmm_hzd
      use HORgrid_options
      use hvdif_options
      use ens_options
      use lun
      use tr3d
      use mem_tstp
      use mem_tracers
      use omp_timing
      use ver
      use lam_options
      use hzd_mod
      use ptopo 
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      logical switch_on_UVW,switch_on_UVW_alh,  switch_on_TR, switch_on_vrtspng_UVT    , &
              switch_on_vrtspng_W, switch_on_eqspng, switch_on_THETA, switch_on_THETA_alh
      logical xch_UV,xch_TT,xch_TR,xch_WZD
      integer i,ik,j,dim,k,k0,km,n
      integer kd0 
      !real, dimension(:,:,:), pointer :: wk2
      !real, dimension (l_minx:l_maxx,l_miny:l_maxy,1:hzd_hyb_bot) ::  u_wrk, v_wrk, w_wrk, zdt_wrk
      !real, dimension (l_minx:l_maxx,l_miny:l_maxy,1:hzd_hyb_top) ::  u_wrk1, v_wrk1, w_wrk1, zdt_wrk1
 


!-------------------------------------------------------------------
!
      call gtmg_start (60, 'HZD_main', 1 )
!  kd0 given by user      
      kd0=hzd_hyb_top 

      if (hzd_conserv_th) then 
         call get_air_density_hlt (1)
      endif

      if (Lun_debug_L) write (Lun_out,1000)

      xch_UV = .false.
      xch_TT = .false.
      xch_TR = .false.
      xch_WZD= .false.
      if (ens_conf) then
          xch_UV = .true.
          xch_TT = .true.
      end if
      switch_on_UVW         = Hzd_lnr         > 0. .and. hzd_hyb_bot < 1
      switch_on_UVW_alh     = Hzd_lnr_z       > 0.
      switch_on_TR          =(Hzd_lnr_tr      > 0.) .and. any(Tr3d_hzd)
      switch_on_THETA       = Hzd_lnr_theta   > 0. .and. hzd_hyb_bot < 1
      switch_on_THETA_alh   = Hzd_lnr_theta_z > 0.
      switch_on_vrtspng_UVT =(Vspng_nk      >=1 ) .and. (Vspng_niter>0)
      switch_on_vrtspng_W   = switch_on_vrtspng_UVT
      switch_on_eqspng      = Eq_nlev       > 1

!$omp single
      call itf_ens_hzd (ut1,vt1,tt1, l_minx,l_maxx,l_miny,l_maxy, G_nk)
!$omp end single

!**********************************
!  Horizontal diffusion on theta  *
!**********************************

      if ( switch_on_THETA ) then
         call gtmg_start (61, 'HZD_theta', 60)
         xch_TT = .true.
         call hzd_theta_hlt ()
         call gtmg_stop(61)
      end if
      if ( switch_on_THETA_alh ) then
         call gtmg_start (62, 'HZD_theta_alh', 60)
         xch_TT = .true.
         call hzd_theta_z ()
         call gtmg_stop(62)
      end if

!**********************************
!  Horizontal diffusion on tracers*
!**********************************

      if ( switch_on_TR ) then
         call gtmg_start (63, 'HZD_tracers', 60)
         xch_TR = .true.
         do i=1, Tr3d_ntr
            if (Tr3d_hzd(i)) then
               if (Hzd_tr_ALH_L) then
                  if ( hzd_conserv_tr )then
                     call hzd_uvwzd_alh(tracers_P(i)%pntr,Hzd_lnR_tr,Hzd_pwr_tr, &
                                        l_minx,l_maxx,l_miny,l_maxy,G_nk,4)
                  else
                     call hzd_uvwzd_alh(tracers_P(i)%pntr,Hzd_lnR_tr,Hzd_pwr_tr,&
                                        l_minx,l_maxx,l_miny,l_maxy,G_nk,3)
                  endif
               else
                  if ( hzd_conserv_tr )then
                    call hzd_expc_deln (tracers_P(i)%pntr,air_dens, Hzd_pwr_tr,Hzd_lnR_tr, &
                                        l_minx,l_maxx,l_miny,l_maxy,G_nk)
                  else
                    call hzd_exp_deln (tracers_P(i)%pntr, Hzd_pwr_tr,&
                          Hzd_lnR_tr, l_minx,l_maxx,l_miny,l_maxy,G_nk)
                  endif
               endif
            end if
         end do
         call gtmg_stop(63)
      end if

!************************
!  Horizontal diffusion *
!************************

      if ( switch_on_UVW ) then
         call gtmg_start (64, 'HZD_bkgrnd', 60)
         xch_UV = .true.
         xch_TT = .true.
         xch_WZD= .true.

          call hzd_exp_deln ( ut1, Hzd_pwr, Hzd_lnR, &
                          l_minx,l_maxx,l_miny,l_maxy,G_nk)
          call hzd_exp_deln ( vt1, Hzd_pwr, Hzd_lnR, &
                          l_minx,l_maxx,l_miny,l_maxy,G_nk)
          call hzd_exp_deln (zdt1, Hzd_pwr, Hzd_lnR, &
                            l_minx,l_maxx,l_miny,l_maxy,G_nk)
          call hzd_exp_deln ( wt1, Hzd_pwr, Hzd_lnR, &
                             l_minx,l_maxx,l_miny,l_maxy,G_nk)
      endif
      if ( switch_on_UVW_alh ) then
        call gtmg_start (65, 'HZD_alh', 60)
         xch_UV = .true.
         xch_TT = .true.
         xch_WZD= .true.
         ! Hybrid diffusion for hzd_hyb_bot > 0 
         if(hzd_hyb_bot > 0) then 
!$omp do collapse(2)
            do ik=1,hzd_hyb_bot
               do j=1-G_haloy,l_nj+G_haloy
                  do i=1-G_halox,l_ni+G_halox
                     u_wrk   (i,j,ik) = ut1 (i,j,l_nk+1-ik) 
                     v_wrk   (i,j,ik) = vt1 (i,j,l_nk+1-ik) 
                     w_wrk   (i,j,ik) = wt1 (i,j,l_nk+1-ik) 
                     zdt_wrk (i,j,ik) = zdt1(i,j,l_nk+1-ik) 
                  end do
               end do
            end do
!$omp end do 
        if(kd0.ne.1) then
!$omp do collapse(2)
            do ik=1,kd0
               do j=1-G_haloy,l_nj+G_haloy
                  do i=1-G_halox,l_ni+G_halox
                     u_wrk1   (i,j,ik) = ut1 (i,j,ik)
                     v_wrk1   (i,j,ik) = vt1 (i,j,ik)
                     w_wrk1   (i,j,ik) = wt1 (i,j,ik)
                     zdt_wrk1 (i,j,ik) = zdt1(i,j,ik)
                  end do
               end do
            end do
!$omp end do 
        endif
            call hzd_uvwzd_alh(ut1 ,Hzd_lnR_z,Hzd_pwr_z,l_minx,l_maxx,l_miny,l_maxy,G_nk,1)
            call hzd_uvwzd_alh(vt1 ,Hzd_lnR_z,Hzd_pwr_z,l_minx,l_maxx,l_miny,l_maxy,G_nk,2)
            call hzd_uvwzd_alh(wt1 ,Hzd_lnR_z,Hzd_pwr_z,l_minx,l_maxx,l_miny,l_maxy,G_nk,3)
            call hzd_uvwzd_alh(zdt1,Hzd_lnR_z,Hzd_pwr_z,l_minx,l_maxx,l_miny,l_maxy,G_nk,3)

            call hzd_exp_deln ( u_wrk, Hzd_pwr, Hzd_lnR,&
                          l_minx,l_maxx,l_miny,l_maxy,hzd_hyb_bot)
            call hzd_exp_deln ( v_wrk, Hzd_pwr, Hzd_lnR,&
                          l_minx,l_maxx,l_miny,l_maxy,hzd_hyb_bot)
            call hzd_exp_deln ( w_wrk, Hzd_pwr, Hzd_lnR,&
                          l_minx,l_maxx,l_miny,l_maxy,hzd_hyb_bot)
            call hzd_exp_deln ( zdt_wrk, Hzd_pwr, Hzd_lnR,&
                          l_minx,l_maxx,l_miny,l_maxy,hzd_hyb_bot)
!$omp do 
            do ik=1,hzd_hyb_bot
               do j=1-G_haloy,l_nj+G_haloy
                  do i=1-G_halox,l_ni+G_halox
                     ut1(i,j,l_nk+1-ik)  = u_wrk (i,j,ik)
                     vt1(i,j,l_nk+1-ik)  = v_wrk (i,j,ik)
                     wt1(i,j,l_nk+1-ik)  = w_wrk (i,j,ik)
                     zdt1(i,j,l_nk+1-ik) = zdt_wrk (i,j,ik)
                  end do
               end do
            end do
!$omp end do 
        if(kd0.ne.1) then
            call hzd_exp_deln ( u_wrk1, Hzd_pwr, Hzd_lnR, &
                          l_minx,l_maxx,l_miny,l_maxy,kd0)
            call hzd_exp_deln ( v_wrk1, Hzd_pwr, Hzd_lnR,&
                          l_minx,l_maxx,l_miny,l_maxy,kd0)
            call hzd_exp_deln ( w_wrk1, Hzd_pwr, Hzd_lnR, &
                          l_minx,l_maxx,l_miny,l_maxy,kd0)
            call hzd_exp_deln ( zdt_wrk1, Hzd_pwr, Hzd_lnR, &
                          l_minx,l_maxx,l_miny,l_maxy,kd0)
!$omp do collapse(2)
            do ik=1,kd0
               do j=1-G_haloy,l_nj+G_haloy
                  do i=1-G_halox,l_ni+G_halox
                     ut1(i,j,ik)  = u_wrk1 (i,j,ik)
                     vt1(i,j,ik)  = v_wrk1 (i,j,ik)
                     wt1(i,j,ik)  = w_wrk1 (i,j,ik)
                     zdt1(i,j,ik) = zdt_wrk1 (i,j,ik)
                  end do
               end do
            end do
!$omp end do 
         endif
         else
            call hzd_uvwzd_alh(ut1, Hzd_lnR_z,Hzd_pwr_z,l_minx,l_maxx,l_miny,l_maxy,G_nk,1)
            call hzd_uvwzd_alh(vt1, Hzd_lnR_z,Hzd_pwr_z,l_minx,l_maxx,l_miny,l_maxy,G_nk,2)
            call hzd_uvwzd_alh(wt1, Hzd_lnR_z,Hzd_pwr_z,l_minx,l_maxx,l_miny,l_maxy,G_nk,3)
            call hzd_uvwzd_alh(zdt1,Hzd_lnR_z,Hzd_pwr_z,l_minx,l_maxx,l_miny,l_maxy,G_nk,3)

         endif

         call gtmg_stop(65)
      end if

!  Vertical sponge  *
!********************

      if ( switch_on_vrtspng_UVT ) then
         call gtmg_start (66, 'V_SPNG', 60)
         xch_UV = .true.
         xch_TT = .true.
         call hzd_exp_del2 ( ut1,  'U', l_minx,l_maxx,l_miny,l_maxy,&
                             Vspng_nk, Hzd_geom_u, F_VV=vt1)
         call hzd_exp_del2 ( tt1, 'M', l_minx,l_maxx,l_miny,l_maxy,&
                             Vspng_nk, Hzd_geom_q)
         call gtmg_stop(66)
      end if

      if ( switch_on_vrtspng_W ) then
         call gtmg_start (66, 'V_SPNG', 60)
         xch_WZD= .true.
         call hzd_exp_del2 ( zdt1, 'M', l_minx,l_maxx,l_miny,l_maxy,&
                                Vspng_nk, Hzd_geom_q)
         call hzd_exp_del2 ( wt1, 'M', l_minx,l_maxx,l_miny,l_maxy,&
                                Vspng_nk, Hzd_geom_q)
         call gtmg_stop(66)
      end if

!**********************
!  Equatorial sponge  *
!**********************

      if ( switch_on_eqspng ) then
         call gtmg_start (67, 'EQUA_SPNG', 60)
         xch_UV= .true.
         call eqspng (ut1,vt1,l_minx,l_maxx,l_miny,l_maxy,G_nk)
         call gtmg_stop(67)
      end if

!$omp single
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
               call yyg_xchng (tracers_P(n)%pntr, l_minx,l_maxx,l_miny,l_maxy,&
                               l_ni, l_nj, G_nk, .true., 'CUBIC', .false.)
            end do
          end if

      end if

      call hzd_smago_main()
!$omp end single

      call gtmg_stop(60)

 1000 format(3X,'MAIN HORIZONTAL DIFFUSION : (S/R HZD_MAIN)')
!
!-------------------------------------------------------------------
!
      return
      end subroutine hzd_main_hlt
