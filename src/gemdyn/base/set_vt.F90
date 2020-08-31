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
      use gmm_vt1
      use gmm_vt0
      use gmm_vt2
      use gmm_vth
      use mem_tracers
      use gmm_tracers
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
#include <arch_specific.hf>

!@objective Initialization of the commons for time-dependent variable.
!           Virtual Memory manager initialization
!@description
!  Notes:
!	The level at time t0 is not created explicitly in the
!	VMM manager tables.
!	It exists in the timestep by renaming the variables
!	created at time level t1 when they are not needed anymore
!	and will bear the same attributes.
!	However, for clarity, a complete comdeck is created
!	for all of the variables at time t0 and the keys
!	will be shuffled at run time according to the renaming
!	performed.
!	The user will therefore use two separate sets
!	of variables for clarity, but will only be using one in
!	memory or disk.
!       vt0 (VMM variables at time t0)
!       vt1 (VMM variables at time t1)
!       vt2 (VMM variables at time t2)
!       vth (VMM variables at time th [t0-dt/2])
!*@/

#include <msg.h>
#include "gmm_gem_flags.hf"

#define SET_GMMUSR_FLAG(MYMETA,MYFLAG) gmm_metadata(MYMETA%l,gmm_attributes(MYMETA%a%key,ior(MYMETA%a%uuid1,MYFLAG),MYMETA%a%uuid2,MYMETA%a%initmode,MYMETA%a%flags))

      type(gmm_metadata) :: mymeta, mymeta3d_nk_u, mymeta3d_nk_v, mymeta3d_nk_t, &
                            mymeta3d_nk_q, mymeta2d_s, mymeta_tracers, meta1d
      integer :: i,istat,dim
      integer :: flag_n, flag_r_n
      integer(kind=INT64) :: flag_m_t,flag_m_u,flag_m_v,flag_m_f,flag_s_f
!
!     ---------------------------------------------------------------
!
      !- Note: gmm_create does NOT keep flags other than GMM_FLAG_IZER+GMM_FLAG_INAN+GMM_FLAG_RSTR
      !  GMM_FLAG_STAG_* need to be added w/ gmm_updatemeta() or in uuid1/2

      flag_n   = GMM_FLAG_IZER
      flag_r_n = GMM_FLAG_RSTR+GMM_FLAG_IZER

      flag_m_t = FLAG_LVL_T                  !thermo   lvl, T
      flag_m_f = FLAG_LVL_M                  !momentum lvl, phi-pt
      flag_m_u = FLAG_LVL_M+GMM_FLAG_STAG_X  !momentum lvl, u-pt
      flag_m_v = FLAG_LVL_M+GMM_FLAG_STAG_Y  !momentum lvl, v-pt
      flag_s_f = 0                           !2d-surf lvl,  phi-pt

      mymeta3d_nk_u  = SET_GMMUSR_FLAG(meta3d_nk  ,flag_m_u)
      mymeta3d_nk_v  = SET_GMMUSR_FLAG(meta3d_nk  ,flag_m_v)
      mymeta3d_nk_t  = SET_GMMUSR_FLAG(meta3d_nk  ,flag_m_t)
      mymeta3d_nk_q  = SET_GMMUSR_FLAG(meta3d_nk1 ,flag_m_f)
      mymeta2d_s     = SET_GMMUSR_FLAG(meta2d     ,flag_s_f)

      gmmk_ut0_s   =  'URT0'
      gmmk_vt0_s   =  'VRT0'
      gmmk_tt0_s   =  'TT0'
      gmmk_st0_s   =  'ST0'
      gmmk_wt0_s   =  'WT0'
      gmmk_qt0_s   =  'QT0'
      gmmk_zdt0_s  = 'ZDT0'

      gmmk_ut1_s   =  'URT1'
      gmmk_vt1_s   =  'VRT1'
      gmmk_tt1_s   =  'TT1'
      gmmk_st1_s   =  'ST1'
      gmmk_wt1_s   =  'WT1'
      gmmk_qt1_s   =  'QT1'
      gmmk_zdt1_s  = 'ZDT1'

      istat = GMM_OK

      istat = min(gmm_create(gmmk_ut0_s ,  ut0, mymeta3d_nk_u, flag_r_n),istat)
      istat = min(gmm_create(gmmk_vt0_s ,  vt0, mymeta3d_nk_v, flag_r_n),istat)
      istat = min(gmm_create(gmmk_tt0_s ,  tt0, mymeta3d_nk_t, flag_r_n),istat)
      istat = min(gmm_create(gmmk_st0_s ,  st0, mymeta2d_s   , flag_r_n),istat)
      istat = min(gmm_create(gmmk_wt0_s ,  wt0, mymeta3d_nk_t, flag_r_n),istat)
      istat = min(gmm_create(gmmk_qt0_s ,  qt0, mymeta3d_nk_q, flag_r_n),istat)
      istat = min(gmm_create(gmmk_zdt0_s, zdt0, mymeta3d_nk_t, flag_r_n),istat)
      if (GMM_IS_ERROR(istat)) then
         call msg(MSG_ERROR,'set_vt ERROR at gmm_create(*t0)')
      end if

      istat = min(gmm_create(gmmk_ut1_s ,  ut1, mymeta3d_nk_u, flag_r_n),istat)
      istat = min(gmm_create(gmmk_vt1_s ,  vt1, mymeta3d_nk_v, flag_r_n),istat)
      istat = min(gmm_create(gmmk_tt1_s ,  tt1, mymeta3d_nk_t, flag_r_n),istat)
      istat = min(gmm_create(gmmk_st1_s ,  st1, mymeta2d_s   , flag_r_n),istat)
      istat = min(gmm_create(gmmk_wt1_s ,  wt1, mymeta3d_nk_t, flag_r_n),istat)
      istat = min(gmm_create(gmmk_qt1_s ,  qt1, mymeta3d_nk_q, flag_r_n),istat)
      istat = min(gmm_create(gmmk_zdt1_s, zdt1, mymeta3d_nk_t, flag_r_n),istat)

      if (GMM_IS_ERROR(istat)) then
         call msg(MSG_ERROR,'set_vt ERROR at gmm_create(*t0)')
      end if

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

      gmmk_xth_s = 'XTH'
      gmmk_yth_s = 'YTH'
      gmmk_zth_s = 'ZTH'

      istat = GMM_OK

      call gmm_build_meta2D(meta1d, &
                            1,l_ni*l_nj*l_nk,0,0,l_ni*l_nj*l_nk, &
                            0,0,0,0,0,0,GMM_NULL_FLAGS)

      mymeta  = SET_GMMUSR_FLAG(meta1d, flag_m_f)
      istat = min(gmm_create(gmmk_xth_s ,xth, mymeta, flag_r_n),istat)
      istat = min(gmm_create(gmmk_yth_s ,yth, mymeta, flag_r_n),istat)
      istat = min(gmm_create(gmmk_zth_s ,zth, mymeta, flag_r_n),istat)

      if (GMM_IS_ERROR(istat)) then
         call msg(MSG_ERROR,'set_vt ERROR at gmm_create(*th)')
      end if
      istat = gmm_get(gmmk_xth_s , xth)
      istat = gmm_get(gmmk_yth_s , yth)
      istat = gmm_get(gmmk_zth_s , zth)

      dim = (l_maxx-l_minx+1) * (l_maxy-l_miny+1) * l_nk * Tr3d_ntr
      call gmm_build_meta1D ( mymeta, 1,dim,0,0,dim, &
                              0,GMM_NULL_FLAGS )
      mymeta_tracers= SET_GMMUSR_FLAG(mymeta,flag_m_f)

      istat = GMM_OK
      nullify(trt0,trt1,trdf)
      istat = min(gmm_create('TRACERS:M',trt0,mymeta_tracers,flag_r_n),istat)
      istat = min(gmm_create('TRACERS:P',trt1,mymeta_tracers,flag_r_n),istat)
      istat = min(gmm_create('TRACERS:D',trdf,mymeta_tracers,flag_r_n),istat)

      if (GMM_IS_ERROR(istat)) then
         call msg(MSG_ERROR,'set_vt ERROR at gmm_create(TR/*)')
      end if

      istat = gmm_get('TRACERS:P',trt1)
      istat = gmm_get('TRACERS:M',trt0)
      istat = gmm_get('TRACERS:D',trdf)

      allocate (tracers_P(Tr3d_ntr), tracers_M(Tr3d_ntr))

      dim = (l_maxx-l_minx+1) * (l_maxy-l_miny+1) * l_nk
      do i=1,Tr3d_ntr
         nullify(tracers_P(i)%pntr)
         tracers_P(i)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trt1((i-1)*dim+1:)
         tracers_M(i)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trt0((i-1)*dim+1:)
!??? a venir peut-etre         tracers_D(i)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trdf((i-1)*dim+1:)
! necessary to create those GMM tracers variables for phy_input (but restart will NOT work)
         istat= min(gmm_create('TR/'//trim(Tr3d_name_S(i))//':M',tracers_M(i)%pntr,mymeta3d_nk_t,0),istat)
         istat= min(gmm_create('TR/'//trim(Tr3d_name_S(i))//':P',tracers_P(i)%pntr,mymeta3d_nk_t,0),istat)
      end do

      !Allocation if Bermejo-Conde LAM ZFL
      !-----------------------------------
      if (adz_BC_LAM_flux==2) then

         istat = GMM_OK
         nullify(trtb)
         istat = min(gmm_create('TRACERS:B',trtb,mymeta_tracers,flag_r_n),istat)

         if (GMM_IS_ERROR(istat)) then
            call msg(MSG_ERROR,'set_vt ERROR-B at gmm_create(TR/*)')
         end if

         istat = gmm_get('TRACERS:B',trtb)

         allocate (tracers_B(Tr3d_ntr))

         dim = (l_maxx-l_minx+1) * (l_maxy-l_miny+1) * l_nk
         do i=1,Tr3d_ntr
            nullify(tracers_B(i)%pntr)
            tracers_B(i)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trtb((i-1)*dim+1:)
            istat= min(gmm_create('TR/'//trim(Tr3d_name_S(i))//':B',tracers_B(i)%pntr,mymeta3d_nk_t,0),istat)
         end do

      end if

      gmmk_airm1_s= 'AIR1'
      gmmk_airm0_s= 'AIR0'
      gmmk_pkps_s = 'PKPS'

      istat = GMM_OK

      istat = min(gmm_create(gmmk_airm1_s,airm1,mymeta3d_nk_t,flag_r_n),istat)
      istat = min(gmm_create(gmmk_airm0_s,airm0,mymeta3d_nk_t,flag_r_n),istat)
      istat = min(gmm_create(gmmk_pkps_s ,pkps ,mymeta3d_nk_t,flag_r_n),istat)

      if (GMM_IS_ERROR(istat)) then
         call msg(MSG_ERROR,'set_vt ERROR at gmm_create(TR_CONS/*)')
      end if

      istat = gmm_get(gmmk_airm1_s,airm1)
      istat = gmm_get(gmmk_airm0_s,airm0)
      istat = gmm_get(gmmk_pkps_s, pkps )

      istat = gmm_get('TR/HU:M' ,hut1)
      istat = gmm_get('TR/HU:P' ,hut0)

      gmmk_pw_uu_plus_s  = 'PW_UU:P'
      gmmk_pw_vv_plus_s  = 'PW_VV:P'
      gmmk_pw_wz_plus_s  = 'PW_WZ:P'
      gmmk_pw_tt_plus_s  = 'PW_TT:P'
      gmmk_pw_pm_plus_s  = 'PW_PM:P'
      gmmk_pw_pm_plus_8_s= 'PW_PM8:P'
      gmmk_pw_pt_plus_s  = 'PW_PT:P'
      gmmk_pw_gz_plus_s  = 'PW_GZ:P'
      gmmk_pw_me_plus_s  = 'PW_ME:P'
      gmmk_pw_p0_plus_s  = 'PW_P0:P'
      gmmk_pw_p0_plus_8_s= 'PW_P08:P'
      gmmk_pw_log_pm_s   = 'PW_LNPM'
      gmmk_pw_log_pt_s   = 'PW_LNPT'

      gmmk_pw_uu_moins_s = 'PW_UU:M'
      gmmk_pw_vv_moins_s = 'PW_VV:M'
      gmmk_pw_tt_moins_s = 'PW_TT:M'
      gmmk_pw_pm_moins_s = 'PW_PM:M'
      gmmk_pw_pm_moins_8_s='PW_PM8:M'
      gmmk_pw_pt_moins_s = 'PW_PT:M'
      gmmk_pw_gz_moins_s = 'PW_GZ:M'
      gmmk_pw_me_moins_s = 'PW_ME:M'
      gmmk_pw_p0_moins_s = 'PW_P0:M'
      gmmk_pw_p0_moins_8_s='PW_P08:M'

      gmmk_pw_uu_copy_s  = 'PW_UU_COPY'
      gmmk_pw_vv_copy_s  = 'PW_VV_COPY'
      gmmk_pw_uslt_s     = 'PW_USLT'
      gmmk_pw_vslt_s     = 'PW_VSLT'

      gmmk_pw_p0_ls_s    = 'PW_P0_LS'

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

      gmmk_smag_s = 'SMAG'
      nullify(smag)
      istat = min(gmm_create(gmmk_smag_s   ,smag  ,meta3d_nk ,flag_n),istat)

      if (GMM_IS_ERROR(istat)) then
         call msg(MSG_ERROR,'set_vt ERROR at gmm_create(PW_*)')
      end if

      istat = gmm_get(gmmk_pw_uu_plus_s   ,pw_uu_plus)
      istat = gmm_get(gmmk_pw_vv_plus_s   ,pw_vv_plus)
      istat = gmm_get(gmmk_pw_wz_plus_s   ,pw_wz_plus)
      istat = gmm_get(gmmk_pw_tt_plus_s   ,pw_tt_plus)
      istat = gmm_get(gmmk_pw_pt_plus_s   ,pw_pt_plus)
      istat = gmm_get(gmmk_pw_gz_plus_s   ,pw_gz_plus)
      istat = gmm_get(gmmk_pw_pm_plus_s   ,pw_pm_plus)
      istat = gmm_get(gmmk_pw_pm_plus_8_s ,pw_pm_plus_8)
      istat = gmm_get(gmmk_pw_me_plus_s   ,pw_me_plus)
      istat = gmm_get(gmmk_pw_p0_plus_s   ,pw_p0_plus)
      istat = gmm_get(gmmk_pw_p0_plus_8_s ,pw_p0_plus_8)
      istat = gmm_get(gmmk_pw_log_pm_s    ,pw_log_pm)
      istat = gmm_get(gmmk_pw_log_pt_s    ,pw_log_pt)
      istat = gmm_get(gmmk_pw_uu_moins_s  ,pw_uu_moins)
      istat = gmm_get(gmmk_pw_vv_moins_s  ,pw_vv_moins)
      istat = gmm_get(gmmk_pw_tt_moins_s  ,pw_tt_moins)
      istat = gmm_get(gmmk_pw_pt_moins_s  ,pw_pt_moins)
      istat = gmm_get(gmmk_pw_gz_moins_s  ,pw_gz_moins)
      istat = gmm_get(gmmk_pw_pm_moins_s  ,pw_pm_moins)
      istat = gmm_get(gmmk_pw_pm_moins_8_s,pw_pm_moins_8)
      istat = gmm_get(gmmk_pw_me_moins_s  ,pw_me_moins)
      istat = gmm_get(gmmk_pw_p0_moins_s  ,pw_p0_moins)
      istat = gmm_get(gmmk_pw_p0_moins_8_s,pw_p0_moins_8)
      istat = gmm_get(gmmk_pw_p0_ls_s     ,pw_p0_ls)
      istat = gmm_get(gmmk_pw_uslt_s      ,pw_uslt)
      istat = gmm_get(gmmk_pw_vslt_s      ,pw_vslt)
      istat = gmm_get(gmmk_pw_uu_copy_s   ,pw_uu_copy)
      istat = gmm_get(gmmk_pw_vv_copy_s   ,pw_vv_copy)
      istat = gmm_get(gmmk_smag_s   ,smag)

      call canonical_cases ("SET_VT")

      gmmk_phy_cplm_s = 'PHYCPLM'
      gmmk_phy_cplt_s = 'PHYCPLT'
      istat = min(gmm_create(gmmk_phy_cplm_s, phy_cplm, meta2d, flag_r_n),istat)
      istat = min(gmm_create(gmmk_phy_cplt_s, phy_cplt, meta2d, flag_r_n),istat)
      gmmk_phy_uu_tend_s  = 'UPT'
      gmmk_phy_vv_tend_s  = 'VPT'
      gmmk_phy_tv_tend_s  = 'TVPT'
      istat = min(gmm_create(gmmk_phy_uu_tend_s, phy_uu_tend, mymeta3d_nk_u, flag_r_n),istat)
      istat = min(gmm_create(gmmk_phy_vv_tend_s, phy_vv_tend, mymeta3d_nk_v, flag_r_n),istat)
      istat = min(gmm_create(gmmk_phy_tv_tend_s, phy_tv_tend, mymeta3d_nk_t, flag_r_n),istat)
      if (GMM_IS_ERROR(istat)) then
         call msg(MSG_ERROR,'set_vt ERROR at gmm_create(PHY)')
      end if
      istat = gmm_get(gmmk_phy_cplm_s, phy_cplm)
      istat = gmm_get(gmmk_phy_cplt_s, phy_cplt)

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_EXPO_H') then
         call exp_set_vt()
      end if

!
!     ---------------------------------------------------------------
      return
      end subroutine set_vt
