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
!/@*
      subroutine itf_phy_update3 (F_apply_L)
      use phy_itf, only: phy_get, phymeta, phy_getmeta
      use itf_phy_filter, only: ipf_smooth_tend
      use gmm_vt1
      use gmm_pw
      use gmm_phy
      use HORgrid_options
      use dyn_fisl_options
      use dynkernel_options
      use adz_options, only: Adz_slt_winds
      use glb_ld
      use cstv
      use lun
      use tr3d
      use gmm_itf_mod
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      logical, intent(in) :: F_apply_L

      logical, parameter :: SMOOTH_EXPLICIT=.false.
      integer, parameter :: SMOOTH_GWD=2

      character(len=GMM_MAXNAMELENGTH) :: trname_S
      integer istat, i,j,k,n, cnt, iteration,iend(3)
      real, dimension(:,:), pointer :: ptr2d
      real, dimension(:,:,:), pointer :: data3d,minus,ptr3d
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,G_nk), target :: tdu,tdv,tv,pw_uu_plus0,pw_vv_plus0,pw_tt_plus0
      real,  dimension(l_ni,l_nj,G_nk) :: qw_phy,qw_dyn
      real(kind=REAL64),dimension(l_minx:l_maxx,l_miny:l_maxy)        :: pr_p0_8
      real(kind=REAL64),dimension(l_minx:l_maxx,l_miny:l_maxy,G_nk+1) :: pr_m_dyn_8,pr_m_phy_8,pr_t_8
      logical :: source_ps_L
!
!-----------------------------------------------------------------
!

   ! The correction due to sources and sinks of specific humidity
   ! is applied only for the case of dry air conservation
   source_ps_L = (Schm_psadj == 2)

   iend = (/-1,-1,l_nk/)

   ! Make diagnosed winds at the lowest thermodynamic level available for advection
   if (Adz_slt_winds) then
      istat = gmm_get (gmmk_pw_uslt_s,pw_uslt)
      ptr2d => pw_uslt(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn)
      istat = phy_get(ptr2d,gmmk_pw_uslt_S,F_npath='V',F_bpath='V')
      istat = gmm_get (gmmk_pw_vslt_s,pw_vslt)
      ptr2d => pw_vslt(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn)
      istat = phy_get(ptr2d,gmmk_pw_vslt_S,F_npath='V',F_bpath='V')
   end if

   ! Retrieve a copy of the PW state before the physics
   nullify(ptr3d)
   istat = gmm_get(gmmk_pw_uu_plus_s,ptr3d); pw_uu_plus0 = ptr3d
   istat = gmm_get(gmmk_pw_vv_plus_s,ptr3d); pw_vv_plus0 = ptr3d
   istat = gmm_get(gmmk_pw_tt_plus_s,ptr3d); pw_tt_plus0 = ptr3d
   if (F_apply_L) then

      if (source_ps_L) then
         qw_phy = 0. ; qw_dyn = 0.
         do n= 1, Tr3d_ntr
            trname_S = 'TR/'//trim(Tr3d_name_S(n))//':P'
            istat = gmm_get(trim(trname_S),data3d)
            if ( (Tr3d_name_S(n)(1:2) == 'HU') .or. &
                 (Schm_wload_L.and.Tr3d_wload(n)) )then
               ptr3d => tdu(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:l_nk)
               if ( phy_get (ptr3d, trim(trname_S), F_npath='V', F_bpath='D', &
                             F_end=iend, F_quiet=.true.) < 0 ) cycle
!$omp parallel private(i,j,k)
!$omp do
               do k=1, l_nk
                  do j=1+pil_s,l_nj-pil_n
                     do i=1+pil_w,l_ni-pil_e
                        qw_phy(i,j,k)= qw_phy(i,j,k) +    tdu(i,j,k)
                        qw_dyn(i,j,k)= qw_dyn(i,j,k) + data3d(i,j,k)
                        data3d(i,j,k)= tdu   (i,j,k)
                     end do
                  end do
               end do
!$omp enddo
!$omp end parallel
            else
               ptr3d => data3d(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:l_nk)
               istat = phy_get (ptr3d, trim(trname_S), F_npath='V', F_bpath='D',&
                                           F_end=iend, F_quiet=.true. )
            end if
         end do
      else
         do k= 1, Tr3d_ntr
            trname_S = 'TR/'//trim(Tr3d_name_S(k))//':P'
            istat = gmm_get(trim(trname_S),data3d)
            ptr3d => data3d(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:l_nk)
            istat = phy_get ( ptr3d, trim(trname_S), F_npath='V', F_bpath='D',&
                              F_end=iend, F_quiet=.true. )
            if (Tr3d_name_S(k)(1:2) == 'HU' .and. SMOOTH_EXPLICIT) istat = ipf_smooth_tend(ptr3d,'SQE')
         end do
      end if

      ! Apply horizontal filtering on tendencies if requeted
      istat = gmm_get (gmmk_pw_uu_plus_s,pw_uu_plus)
      istat = gmm_get (gmmk_pw_vv_plus_s,pw_vv_plus)
      istat = gmm_get (gmmk_pw_tt_plus_s,pw_tt_plus)

      ptr3d => pw_uu_plus(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:l_nk)
      istat = phy_get(ptr3d,gmmk_pw_uu_plus_s,F_npath='V',F_bpath='D',F_end=iend)
      istat = ipf_smooth_tend(ptr3d,'ugwd_td1',SMOOTH_GWD)

      ptr3d => pw_vv_plus(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:l_nk)
      istat = phy_get(ptr3d,gmmk_pw_vv_plus_s,F_npath='V',F_bpath='D',F_end=iend)
      istat = ipf_smooth_tend(ptr3d,'vgwd_td1',SMOOTH_GWD)

      ptr3d => pw_tt_plus(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:l_nk)
      istat = phy_get(ptr3d,gmmk_pw_tt_plus_s,F_npath='V',F_bpath='D',F_end=iend)
      if (SMOOTH_EXPLICIT) istat = ipf_smooth_tend(ptr3d,'STE')

      ! Compute moisture sources for dry air conservation
      if (source_ps_L) then

         iteration = 1

         istat = gmm_get(gmmk_st1_s,st1)

         !Obtain pressure levels
         !----------------------
         call calc_pressure_8 (pr_m_dyn_8,pr_t_8,pr_p0_8,st1,l_minx,l_maxx,l_miny,l_maxy,l_nk)

         pr_m_dyn_8(:,:,l_nk+1) = pr_p0_8(:,:)

         do while (iteration<4)

            !Obtain pressure levels
            !----------------------
            call calc_pressure_8 (pr_m_phy_8,pr_t_8,pr_p0_8,st1,l_minx,l_maxx,l_miny,l_maxy,l_nk)

            pr_m_phy_8(:,:,l_nk+1) = pr_p0_8(:,:)

            pr_p0_8(:,:) = pr_m_phy_8(:,:,1)

            !Estimate source of surface pressure due to fluxes of water:
            !-----------------------------------------------------------------------------------------------------
            !Vertical_Integral [d(p_phy)_k+1] = Vertical_Integral [ d(p_phy)_k q_phy + d(p_dyn) (1-q_dyn) based on
            !-----------------------------------------------------------------------------------------------------
            !d(ps) = Vertical_Integral [ d(qw)/(1-qw_phy)] d(pi) (Claude Girard)
            !-----------------------------------------------------------------------------------------------------
!$omp parallel private(i,j,k) &
!$omp shared(qw_dyn,qw_phy,pr_m_dyn_8,pr_m_phy_8,pr_p0_8)
!$omp do
            do j=1+pil_s,l_nj-pil_n
               do k=1,l_nk
                  do i=1+pil_w,l_ni-pil_e
                     pr_p0_8(i,j)= pr_p0_8(i,j) + (1.0-qw_dyn(i,j,k)) * (pr_m_dyn_8(i,j,k+1)-pr_m_dyn_8(i,j,k)) + &
                                                       qw_phy(i,j,k)  * (pr_m_phy_8(i,j,k+1)-pr_m_phy_8(i,j,k))
                  end do
               end do
            end do
!$omp enddo
!$omp end parallel

!$omp parallel private(i,j) shared(pr_p0_8,st1)
!$omp do
            do j=1+pil_s,l_nj-pil_n
               do i=1+pil_w,l_ni-pil_e
                  st1(i,j)= log(pr_p0_8(i,j)/Cstv_pref_8)
               end do
            end do
!$omp enddo
!$omp end parallel

            iteration = iteration + 1

         end do

         if (Lun_out>0) write(Lun_out,*) ''
         if (Lun_out>0) write(Lun_out,*) '--------------------------------------'
         if (Lun_out>0) write(Lun_out,*) 'SOURCE_PS is done for DRY AIR (REAL64)'
         if (Lun_out>0) write(Lun_out,*) '--------------------------------------'
         if (Lun_out>0) write(Lun_out,*) ''

         call pw_update_GPW()

      end if

      ! Compute tendencies and reset physical world if requested
      if (Schm_phycpl_S == 'RHS' .or. Schm_phycpl_S == 'AVG') then

         istat = gmm_get(gmmk_phy_uu_tend_s,phy_uu_tend)
         istat = gmm_get(gmmk_phy_vv_tend_s,phy_vv_tend)
         istat = gmm_get(gmmk_phy_tv_tend_s,phy_tv_tend)

         tdu(1:l_ni,1:l_nj,1:l_nk) = pw_uu_plus(1:l_ni,1:l_nj,1:l_nk) - pw_uu_plus0(1:l_ni,1:l_nj,1:l_nk)
         tdv(1:l_ni,1:l_nj,1:l_nk) = pw_vv_plus(1:l_ni,1:l_nj,1:l_nk) - pw_vv_plus0(1:l_ni,1:l_nj,1:l_nk)
         call hwnd_stag(phy_uu_tend,phy_vv_tend,tdu,tdv,l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

         istat = gmm_get(gmmk_tt1_s, tt1)
         call tt2virt (tv,.true.,l_minx,l_maxx,l_miny,l_maxy,l_nk)
         phy_tv_tend(1:l_ni,1:l_nj,1:l_nk) = tv(1:l_ni,1:l_nj,1:l_nk) - tt1(1:l_ni,1:l_nj,1:l_nk)

         phy_uu_tend = phy_uu_tend/Cstv_dt_8
         phy_vv_tend = phy_vv_tend/Cstv_dt_8
         phy_tv_tend = phy_tv_tend/Cstv_dt_8

         RESET_PW: if (Schm_phycpl_S == 'RHS') then
            pw_uu_plus = pw_uu_plus0
            pw_vv_plus = pw_vv_plus0
            pw_tt_plus = pw_tt_plus0
         else
            pw_uu_plus = pw_uu_plus0 + Cstv_bA_m_8*(pw_uu_plus-pw_uu_plus0)
            pw_vv_plus = pw_vv_plus0 + Cstv_bA_m_8*(pw_vv_plus-pw_vv_plus0)
            pw_tt_plus = pw_tt_plus0 + Cstv_bA_8*(pw_tt_plus-pw_tt_plus0)
         end if RESET_PW

      else
         if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then !For SPLIT in GEM-H
            istat = gmm_get(gmmk_phy_tv_tend_s,phy_tv_tend)
            istat = gmm_get(gmmk_tt1_s, tt1)
            call tt2virt (tv,.true.,l_minx,l_maxx,l_miny,l_maxy,l_nk)
            phy_tv_tend(1:l_ni,1:l_nj,1:l_nk) = tv(1:l_ni,1:l_nj,1:l_nk) - tt1(1:l_ni,1:l_nj,1:l_nk)
            phy_tv_tend = phy_tv_tend/Cstv_dt_8
         end if         
      end if

   else

      cnt = 0
      do k= 1, Tr3d_ntr
         if (trim(Tr3d_name_S(k)) == 'HU' .or.                          &
             any(NTR_Tr3d_name_S(1:NTR_Tr3d_ntr)==trim(Tr3d_name_S(k))))&
             cycle
         trname_S = 'TR/'//trim(Tr3d_name_S(k))//':P'
         istat = gmm_get(trim(trname_S),data3d)
         ptr3d => data3d(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:l_nk)
         if ( phy_get ( ptr3d, trim(trname_S), F_npath='V', F_bpath='D',&
                        F_end=iend, F_quiet=.true. ) < 0 ) cycle
         trname_S = 'TR/'//trim(Tr3d_name_S(k))//':M'
         if (Grd_yinyang_L) &
         call yyg_xchng (data3d, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj, &
                         G_nk,.true., 'CUBIC', .false.)

         istat = gmm_get(trim(trname_S),minus)
         minus = data3d
         cnt   = cnt + 1
      end do

      if (cnt > 0) then
         istat = gmm_get(gmmk_tt1_s, tt1)
         call tt2virt (tt1, .true., l_minx,l_maxx,l_miny,l_maxy,l_nk)
         if (Grd_yinyang_L) then
            call yyg_xchng (tt1, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj, &
                            G_nk, .false., 'CUBIC', .false.)
            call pw_update_T
         end if
         call pw_update_GPW
      end if

   end if
!
!-----------------------------------------------------------------
!
   return
   end subroutine itf_phy_update3
