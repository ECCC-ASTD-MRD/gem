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
      subroutine set_vt()
      use adz_options
      use dynkernel_options
      use gmm_contiguous
      use gmm_vt2
      use gmm_vt1
      use gmm_vt0
      use gmm_vth
      use mem_tracers
      use gmm_tracers
      use gmm_geof
      use gmm_pw
      use gmm_smag
      use gmm_phy
      use gem_options
      use glb_ld
      use lun
      use tr3d
      use gmm_itf_mod
      use var_gmm
      use, intrinsic :: iso_fortran_env
      implicit none

#include <msg.h>
#include "gmm_gem_flags.hf"

      type(gmm_metadata) :: mymeta, meta1d

      integer :: i,istat,dim,dimH
      integer :: flag_n, flag_r_n
!
!     ---------------------------------------------------------------
!
      flag_n   = GMM_FLAG_IZER
      flag_r_n = GMM_FLAG_RSTR+GMM_FLAG_IZER

      istat = GMM_OK !hopefully 0

      dimh= (l_maxx-l_minx+1) * (l_maxy-l_miny+1)
      dim = dimH * (6*l_nk + 2)
      gmm_nbplans = dim

      call gmm_build_meta1D ( mymeta, 1,dim,0,0,dim, &
                              0,GMM_NULL_FLAGS )

      nullify(dynt2,dynt1,dynt0)
      istat = min(gmm_create('DYNT2',dynt2,mymeta,flag_r_n),istat)
      istat = min(gmm_create('DYNT1',dynt1,mymeta,flag_r_n),istat)
      istat = min(gmm_create('DYNT0',dynt0,mymeta,flag_r_n),istat)

      allocate (timlvl2(7),timlvl1(7),timlvl0(7))

      dim = dimH * l_nk

      timlvl2(1)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => dynt2(      1:)
      timlvl2(2)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => dynt2(  dim+1:)
      timlvl2(3)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => dynt2(2*dim+1:)
      timlvl2(4)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk+1)=> dynt2(3*dim+1:)
      timlvl2(5)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => dynt2(4*dim+dimH+1:)
      timlvl2(6)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => dynt2(5*dim+dimH+1:)
      timlvl2(1)%pntr_2d(l_minx:l_maxx,l_miny:l_maxy) => dynt2(6*dim+dimH+1:)

      istat= min(gmm_create(gmmk_ut2_s ,timlvl2(1)%pntr_3d,meta3d_nk ,0),istat)
      istat= min(gmm_create(gmmk_vt2_s ,timlvl2(2)%pntr_3d,meta3d_nk ,0),istat)
      istat= min(gmm_create(gmmk_tt2_s ,timlvl2(3)%pntr_3d,meta3d_nk ,0),istat)
      istat= min(gmm_create(gmmk_qt2_s ,timlvl2(4)%pntr_3d,meta3d_nk1 ,0),istat)
      istat= min(gmm_create(gmmk_wt2_s ,timlvl2(5)%pntr_3d,meta3d_nk ,0),istat)
      istat= min(gmm_create(gmmk_zdt2_s ,timlvl2(6)%pntr_3d,meta3d_nk ,0),istat)
      istat= min(gmm_create(gmmk_st2_s ,timlvl2(1)%pntr_2d,meta2d ,0),istat)

      timlvl1(1)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => dynt1(      1:)
      timlvl1(2)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => dynt1(  dim+1:)
      timlvl1(3)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => dynt1(2*dim+1:)
      timlvl1(4)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk+1)=> dynt1(3*dim+1:)
      timlvl1(5)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => dynt1(4*dim+dimH+1:)
      timlvl1(6)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => dynt1(5*dim+dimH+1:)
      timlvl1(1)%pntr_2d(l_minx:l_maxx,l_miny:l_maxy) => dynt1(6*dim+dimH+1:)

      istat= min(gmm_create(gmmk_ut1_s ,timlvl1(1)%pntr_3d,meta3d_nk ,0),istat)
      istat= min(gmm_create(gmmk_vt1_s ,timlvl1(2)%pntr_3d,meta3d_nk ,0),istat)
      istat= min(gmm_create(gmmk_tt1_s ,timlvl1(3)%pntr_3d,meta3d_nk ,0),istat)
      istat= min(gmm_create(gmmk_qt1_s ,timlvl1(4)%pntr_3d,meta3d_nk1 ,0),istat)
      istat= min(gmm_create(gmmk_wt1_s ,timlvl1(5)%pntr_3d,meta3d_nk ,0),istat)
      istat= min(gmm_create(gmmk_zdt1_s ,timlvl1(6)%pntr_3d,meta3d_nk ,0),istat)
      istat= min(gmm_create(gmmk_st1_s ,timlvl1(1)%pntr_2d,meta2d ,0),istat)

      timlvl0(1)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => dynt0(      1:)
      timlvl0(2)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => dynt0(  dim+1:)
      timlvl0(3)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => dynt0(2*dim+1:)
      timlvl0(4)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk+1)=> dynt0(3*dim+1:)
      timlvl0(5)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => dynt0(4*dim+dimH+1:)
      timlvl0(6)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => dynt0(5*dim+dimH+1:)
      timlvl0(1)%pntr_2d(l_minx:l_maxx,l_miny:l_maxy) => dynt0(6*dim+dimH+1:)

      istat= min(gmm_create(gmmk_ut0_s ,timlvl0(1)%pntr_3d,meta3d_nk ,0),istat)
      istat= min(gmm_create(gmmk_vt0_s ,timlvl0(2)%pntr_3d,meta3d_nk ,0),istat)
      istat= min(gmm_create(gmmk_tt0_s ,timlvl0(3)%pntr_3d,meta3d_nk ,0),istat)
      istat= min(gmm_create(gmmk_qt0_s ,timlvl0(4)%pntr_3d,meta3d_nk1 ,0),istat)
      istat= min(gmm_create(gmmk_wt0_s ,timlvl0(5)%pntr_3d,meta3d_nk ,0),istat)
      istat= min(gmm_create(gmmk_zdt0_s ,timlvl0(6)%pntr_3d,meta3d_nk ,0),istat)
      istat= min(gmm_create(gmmk_st0_s ,timlvl0(1)%pntr_2d,meta2d ,0),istat)

      call gem_error(istat,'set_vt','ERROR at gmm_create')

      istat = gmm_get (gmmk_ut0_s , ut0)
      istat = gmm_get (gmmk_vt0_s , vt0)
      istat = gmm_get (gmmk_tt0_s , tt0)
      istat = gmm_get (gmmk_st0_s , st0)
      istat = gmm_get (gmmk_wt0_s , wt0)
      istat = gmm_get (gmmk_qt0_s , qt0)
      istat = gmm_get (gmmk_zdt0_s, zdt0)
      istat = gmm_get (gmmk_ut1_s , ut1)
      istat = gmm_get (gmmk_vt1_s , vt1)
      istat = gmm_get (gmmk_tt1_s , tt1)
      istat = gmm_get (gmmk_st1_s , st1)
      istat = gmm_get (gmmk_wt1_s , wt1)
      istat = gmm_get (gmmk_qt1_s , qt1)
      istat = gmm_get (gmmk_zdt1_s, zdt1)
      istat = gmm_get (gmmk_ut2_s , ut2)
      istat = gmm_get (gmmk_vt2_s , vt2)
      istat = gmm_get (gmmk_tt2_s , tt2)
      istat = gmm_get (gmmk_st2_s , st2)
      istat = gmm_get (gmmk_wt2_s , wt2)
      istat = gmm_get (gmmk_qt2_s , qt2)
      istat = gmm_get (gmmk_zdt2_s, zdt2)

      istat = GMM_OK

      call gmm_build_meta2D(meta1d, &
                            1,l_ni*l_nj*l_nk,0,0,l_ni*l_nj*l_nk, &
                            0,0,0,0,0,0,GMM_NULL_FLAGS)

      istat = min(gmm_create(gmmk_xth_s ,xth, meta1d, flag_r_n),istat)
      istat = min(gmm_create(gmmk_yth_s ,yth, meta1d, flag_r_n),istat)
      istat = min(gmm_create(gmmk_zth_s ,zth, meta1d, flag_r_n),istat)

      if (GMM_IS_ERROR(istat)) then
         call msg(MSG_ERROR,'set_vt ERROR at gmm_create(*th)')
      end if

      dim = (l_maxx-l_minx+1) * (l_maxy-l_miny+1) * l_nk * Tr3d_ntr
      call gmm_build_meta1D ( mymeta, 1,dim,0,0,dim, &
                              0,GMM_NULL_FLAGS )

      istat = GMM_OK
      nullify(trt0,trt1,trt2,trdf)
      istat = min(gmm_create('TRACERS:t0',trt0,mymeta,flag_r_n),istat)
      istat = min(gmm_create('TRACERS:t1',trt1,mymeta,flag_r_n),istat)
      istat = min(gmm_create('TRACERS:t2',trt2,mymeta,flag_r_n),istat)
      istat = min(gmm_create('TRACERS:df',trdf,mymeta,flag_r_n),istat)

      if (GMM_IS_ERROR(istat)) then
         call msg(MSG_ERROR,'set_vt ERROR at gmm_create(TR/*)')
      end if

      allocate (tracers_P(Tr3d_ntr), tracers_M(Tr3d_ntr), tracers_t2(Tr3d_ntr))

      dim = (l_maxx-l_minx+1) * (l_maxy-l_miny+1) * l_nk
      do i=1,Tr3d_ntr
         nullify(tracers_P(i)%pntr)
         tracers_P(i)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trt1((i-1)*dim+1:)
         tracers_M(i)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trt0((i-1)*dim+1:)
         tracers_t2(i)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trt2((i-1)*dim+1:)
! necessary to create those GMM tracers variables for phy_input
         istat= min(gmm_create('TR/'//trim(Tr3d_name_S(i))//':M',tracers_M(i)%pntr,meta3d_nk,0),istat)
         istat= min(gmm_create('TR/'//trim(Tr3d_name_S(i))//':P',tracers_P(i)%pntr,meta3d_nk,0),istat)
         istat= min(gmm_create('TR/'//trim(Tr3d_name_S(i))//':t2',tracers_t2(i)%pntr,meta3d_nk,0),istat)
      end do

      !Allocation if Bermejo-Conde LAM ZFL
      !-----------------------------------
      if (adz_BC_LAM_flux==2) then

         istat = GMM_OK
         nullify(trtb)
         istat = min(gmm_create('TRACERS:B',trtb,mymeta,flag_r_n),istat)

         if (GMM_IS_ERROR(istat)) then
            call msg(MSG_ERROR,'set_vt ERROR-B at gmm_create(TR/*)')
         end if

         allocate (tracers_B(Tr3d_ntr))

         dim = (l_maxx-l_minx+1) * (l_maxy-l_miny+1) * l_nk
         do i=1,Tr3d_ntr
            nullify(tracers_B(i)%pntr)
            tracers_B(i)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trtb((i-1)*dim+1:)
            istat= min(gmm_create('TR/'//trim(Tr3d_name_S(i))//':B',tracers_B(i)%pntr,meta3d_nk,0),istat)
         end do

      end if

      istat = GMM_OK

      istat = min(gmm_create(gmmk_airm1_s,airm1,meta3d_nk,flag_r_n),istat)
      istat = min(gmm_create(gmmk_airm0_s,airm0,meta3d_nk,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pkps_s ,pkps ,meta3d_nk,flag_r_n),istat)

      istat = min(gmm_create(gmmk_dgzm_s   ,dgzm  ,meta3d_nk  ,0),istat)
      istat = min(gmm_create(gmmk_dgzt_s   ,dgzt  ,meta3d_nk  ,0),istat)

      if (GMM_IS_ERROR(istat)) then
         call msg(MSG_ERROR,'set_vt ERROR at gmm_create(TR_CONS/*)')
      end if

      nullify(pw_uu_plus ,pw_vv_plus ,pw_wz_plus ,pw_tt_plus ,pw_pm_plus,pw_pt_plus,pw_gz_plus)
      nullify(pw_uu_moins,pw_vv_moins,            pw_tt_moins,pw_pm_moins,pw_gz_moins)
      nullify(pw_uu_copy ,pw_vv_copy, pw_log_pm, pw_log_pt)
      nullify(pw_pm_plus_8,pw_p0_plus_8,pw_pm_moins_8,pw_p0_moins_8)
      istat = GMM_OK

      istat = min(gmm_create(gmmk_pw_uu_plus_s   ,pw_uu_plus  ,meta3d_nk  ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_vv_plus_s   ,pw_vv_plus  ,meta3d_nk  ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_wz_plus_s   ,pw_wz_plus  ,meta3d_nk  ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_tt_plus_s   ,pw_tt_plus  ,meta3d_nk  ,flag_r_n),istat)

      istat = min(gmm_create(gmmk_pw_pt_plus_s   ,pw_pt_plus  ,meta3d_nk1 ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_gz_plus_s   ,pw_gz_plus  ,meta3d_nk  ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_pm_plus_s   ,pw_pm_plus  ,meta3d_nk1 ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_pm_plus_8_s ,pw_pm_plus_8,meta3d_nk1 ,flag_r_n),istat)

      istat = min(gmm_create(gmmk_pw_me_plus_s   ,pw_me_plus  ,meta2d     ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_p0_plus_s   ,pw_p0_plus  ,meta2d     ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_p0_plus_8_s ,pw_p0_plus_8,meta2d     ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_log_pm_s    ,pw_log_pm   ,meta3d_nk1 ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_log_pt_s    ,pw_log_pt   ,meta3d_nk1 ,flag_r_n),istat)

      istat = min(gmm_create(gmmk_pw_uu_moins_s  ,pw_uu_moins ,meta3d_nk  ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_vv_moins_s  ,pw_vv_moins ,meta3d_nk  ,flag_r_n),istat)

      istat = min(gmm_create(gmmk_pw_pt_moins_s  ,pw_pt_moins ,meta3d_nk1 ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_gz_moins_s  ,pw_gz_moins ,meta3d_nk  ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_pm_moins_s  ,pw_pm_moins ,meta3d_nk1 ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_pm_moins_8_s,pw_pm_moins_8,meta3d_nk1,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_tt_moins_s  ,pw_tt_moins ,meta3d_nk  ,flag_r_n),istat)

      istat = min(gmm_create(gmmk_pw_me_moins_s  ,pw_me_moins ,meta2d     ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_p0_moins_s  ,pw_p0_moins ,meta2d     ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_p0_moins_8_s,pw_p0_moins_8,meta2d    ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_p0_ls_s     ,pw_p0_ls    ,meta2d     ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_uslt_s      ,pw_uslt     ,meta2d     ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_vslt_s      ,pw_vslt     ,meta2d     ,flag_r_n),istat)

      istat = min(gmm_create(gmmk_pw_uu_copy_s   ,pw_uu_copy  ,meta3d_nk  ,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pw_vv_copy_s   ,pw_vv_copy  ,meta3d_nk  ,flag_r_n),istat)

      nullify(smag)
      istat = min(gmm_create(gmmk_smag_s   ,smag  ,meta3d_nk ,flag_n),istat)

      if (GMM_IS_ERROR(istat)) then
         call msg(MSG_ERROR,'set_vt ERROR at gmm_create(PW_*)')
      end if

      istat = min(gmm_create(gmmk_phy_cplm_s, phy_cplm, meta2d, flag_r_n),istat)
      istat = min(gmm_create(gmmk_phy_cplt_s, phy_cplt, meta2d, flag_r_n),istat)
      istat = min(gmm_create(gmmk_phy_uu_tend_s, phy_uu_tend, meta3d_nk, flag_r_n),istat)
      istat = min(gmm_create(gmmk_phy_vv_tend_s, phy_vv_tend, meta3d_nk, flag_r_n),istat)
      istat = min(gmm_create(gmmk_phy_tv_tend_s, phy_tv_tend, meta3d_nk, flag_r_n),istat)
      if (GMM_IS_ERROR(istat)) then
         call msg(MSG_ERROR,'set_vt ERROR at gmm_create(PHY)')
      end if

      call canonical_cases ("SET_VT")

!
!     ---------------------------------------------------------------
      return
      end subroutine set_vt