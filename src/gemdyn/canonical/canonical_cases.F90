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

!**s/r canonical_cases - Various actions related to canonical cases (Williamson/DCMIP)

      subroutine canonical_cases (F_action_S)

      use canonical
      use cstv
      use ctrl
      use dcmip_options
      use dyn_fisl_options
      use dynkernel_options
      use glb_ld
      use gmm_itf_mod
      use gmm_vt0
      use gmm_vt1
      use gmm_pw
      use lun
      use mem_tracers
      use step_options
      use tdpack, only : rgasd_8, cpd_8
      use tr3d
      use var_gmm
      use ver
      use VERgrid_options
      use wil_options

      use, intrinsic :: iso_fortran_env
      implicit none

      !arguments
      character(len=*), intent(in) :: F_action_S

      !object
      !========================================================================
      !     F_action_S ='SET_ZETA': Print dcmip_HEIGHTS
      !     F_action_S ='SET_VT'  : Initialize gmm variables
      !     F_action_S ='BAC'     : Back subtitution (WINDS)
      !     F_action_S ='PHY'     : Physics
      !     F_action_S ='VRD'     : Vertical Diffusion
      !     F_action_S ='ERR'     : Evaluate Errors/Diagnostics
      !     F_action_S ='OUT'     : Output dependent variables on standard file
      !========================================================================

#include "gmm_gem_flags.hf"
#include <msg.h>

#define SET_GMMUSR_FLAG(MYMETA,MYFLAG) gmm_metadata(MYMETA%l,gmm_attributes(MYMETA%a%key,ior(MYMETA%a%uuid1,MYFLAG),MYMETA%a%uuid2,MYMETA%a%initmode,MYMETA%a%flags))

      type(gmm_metadata) :: mymeta3d_nk_u,mymeta3d_nk_v,mymeta3d_nk_t,mymeta2d_s
      integer ::  istat,flag_r_n,i,j,k,n,pnip1
      integer(kind=INT64) :: flag_m_t,flag_m_u,flag_m_v,flag_s_f
      real, pointer, dimension(:,:,:) :: hu,cl,cl2,tr,tr_r,tr_e
      character(len=8) :: dumc
      real :: dcmip_height,dcmip_heightp1
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,G_nk) :: bidon,th
      real(kind=REAL64) :: pr_8
      logical :: Terminator_L
      real, parameter :: CLY_REF = 4.*10.**(-6)
!
!     ---------------------------------------------------------------
!
      if (.not. Ctrl_testcases_L) return

      Terminator_L = Dcmip_Terminator_L.or.Williamson_Terminator_L

      !-------------------
      !Print dcmip_HEIGHTS
      !-------------------
      if (F_action_S=="SET_ZETA") then

         if (Dcmip_case==0.or.trim(Dynamics_Kernel_S)/='DYNAMICS_FISL_P') return

         if (Lun_out > 0) then
            write (Lun_out,1005) G_nk,Hyb_rcoef

            do k=1,G_nk
               dcmip_height  =-8780.2*alog(Ver_hyb%m(k))
               if (k < G_nk)&
               dcmip_heightp1 =-8780.2*alog(Ver_hyb%m(k+1))
               if (k == G_nk) dcmip_heightp1 = 0.
               call convip(pnip1,Ver_hyb%m(k),5,1,dumc,.false.)
               write (Lun_out,1006) k,Ver_hyb%m(k),dcmip_height, &
                                    dcmip_height-dcmip_heightp1,pnip1
            end do

         end if

      !------------------------
      !Initialize gmm variables
      !------------------------
      else if (F_action_S=="SET_VT") then

         flag_r_n = GMM_FLAG_RSTR+GMM_FLAG_IZER

         flag_m_t = FLAG_LVL_T                  !thermo   lvl, T
         flag_m_u = FLAG_LVL_M+GMM_FLAG_STAG_X  !momentum lvl, u-pt
         flag_m_v = FLAG_LVL_M+GMM_FLAG_STAG_Y  !momentum lvl, v-pt
         flag_s_f = 0                           !2d-surf  lvl, phi-pt

         mymeta3d_nk_u  = SET_GMMUSR_FLAG(meta3d_nk ,flag_m_u)
         mymeta3d_nk_v  = SET_GMMUSR_FLAG(meta3d_nk ,flag_m_v)
         mymeta3d_nk_t  = SET_GMMUSR_FLAG(meta3d_nk ,flag_m_t)
         mymeta2d_s     = SET_GMMUSR_FLAG(meta2d    ,flag_s_f)

         gmmk_pth_s   = 'PTH'
         gmmk_thbase_s= 'THBA'
         gmmk_thfull_s= 'THFU'
         gmmk_dtv_s   = 'DTV'

         gmmk_cly_s  = 'CLY'
         gmmk_acl_s  = 'ACL'
         gmmk_acl2_s = 'ACL2'
         gmmk_acly_s = 'ACLY'

         gmmk_irt_s  = 'IRT'
         gmmk_art_s  = 'ART'
         gmmk_wrt_s  = 'WRT'

         gmmk_uref_s = 'UREF'
         gmmk_vref_s = 'VREF'
         gmmk_wref_s = 'WREF'
         gmmk_zdref_s= 'ZDRF'
         gmmk_qvref_s= 'QVRF'
         gmmk_qcref_s= 'QCRF'
         gmmk_qrref_s= 'RWRF'
         gmmk_thref_s= 'THRF'

         gmmk_q1ref_s= 'QR1'
         gmmk_q2ref_s= 'QR2'
         gmmk_q3ref_s= 'QR3'
         gmmk_q4ref_s= 'QR4'

         gmmk_q1err_s= 'QE1'
         gmmk_q2err_s= 'QE2'
         gmmk_q3err_s= 'QE3'
         gmmk_q4err_s= 'QE4'

         gmmk_clyref_s= 'CLYR'
         gmmk_clyerr_s= 'CLYE'

         istat = GMM_OK

         istat = min(gmm_create(gmmk_pth_s,   pth,   mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_thbase_s,thbase,mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_thfull_s,thfull,mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_dtv_s,   dtv,   mymeta3d_nk_t, flag_r_n),istat)

         istat = min(gmm_create(gmmk_cly_s,   cly,   mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_acl_s,   acl,   mymeta2d_s   , flag_r_n),istat)
         istat = min(gmm_create(gmmk_acl2_s,  acl2,  mymeta2d_s   , flag_r_n),istat)
         istat = min(gmm_create(gmmk_acly_s,  acly,  mymeta2d_s   , flag_r_n),istat)

         istat = min(gmm_create(gmmk_irt_s,   irt,   mymeta2d_s   , flag_r_n),istat)
         istat = min(gmm_create(gmmk_art_s,   art,   mymeta2d_s   , flag_r_n),istat)
         istat = min(gmm_create(gmmk_wrt_s,   wrt,   mymeta2d_s   , flag_r_n),istat)

         istat = min(gmm_create(gmmk_uref_s,  uref,  mymeta3d_nk_u, flag_r_n),istat)
         istat = min(gmm_create(gmmk_vref_s,  vref,  mymeta3d_nk_v, flag_r_n),istat)
         istat = min(gmm_create(gmmk_wref_s,  wref,  mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_zdref_s, zdref, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_qvref_s, qvref, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_qcref_s, qcref, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_qrref_s, qrref, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_thref_s, thref, mymeta3d_nk_t, flag_r_n),istat)

         istat = min(gmm_create(gmmk_q1ref_s, q1ref, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_q2ref_s, q2ref, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_q3ref_s, q3ref, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_q4ref_s, q4ref, mymeta3d_nk_t, flag_r_n),istat)

         istat = min(gmm_create(gmmk_q1err_s, q1err, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_q2err_s, q2err, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_q3err_s, q3err, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_q4err_s, q4err, mymeta3d_nk_t, flag_r_n),istat)

         istat = min(gmm_create(gmmk_clyref_s,clyref,mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_clyerr_s,clyerr,mymeta3d_nk_t, flag_r_n),istat)

         if (GMM_IS_ERROR(istat)) &
             call msg(MSG_ERROR,'set_vt ERROR at gmm_create(CANO)')

         istat = gmm_get(gmmk_pth_s   , pth)
         istat = gmm_get(gmmk_thbase_s, thbase)
         istat = gmm_get(gmmk_thfull_s, thfull)
         istat = gmm_get(gmmk_dtv_s   , dtv)
         istat = gmm_get(gmmk_cly_s   , cly)
         istat = gmm_get(gmmk_acl_s   , acl)
         istat = gmm_get(gmmk_acl2_s  , acl2)
         istat = gmm_get(gmmk_acly_s  , acly)
         istat = gmm_get(gmmk_irt_s   , irt)
         istat = gmm_get(gmmk_art_s   , art)
         istat = gmm_get(gmmk_wrt_s   , wrt)
         istat = gmm_get(gmmk_uref_s  , uref)
         istat = gmm_get(gmmk_vref_s  , vref)
         istat = gmm_get(gmmk_wref_s  , wref)
         istat = gmm_get(gmmk_zdref_s , zdref)
         istat = gmm_get(gmmk_qvref_s , qvref)
         istat = gmm_get(gmmk_qcref_s , qcref)
         istat = gmm_get(gmmk_qrref_s , qrref)
         istat = gmm_get(gmmk_thref_s , thref)
         istat = gmm_get(gmmk_q1ref_s , q1ref)
         istat = gmm_get(gmmk_q2ref_s , q2ref)
         istat = gmm_get(gmmk_q3ref_s , q3ref)
         istat = gmm_get(gmmk_q4ref_s , q4ref)
         istat = gmm_get(gmmk_q1err_s , q1err)
         istat = gmm_get(gmmk_q2err_s , q2err)
         istat = gmm_get(gmmk_q3err_s , q3err)
         istat = gmm_get(gmmk_q4err_s , q4err)
         istat = gmm_get(gmmk_clyref_s, clyref)
         istat = gmm_get(gmmk_clyerr_s, clyerr)

      !------------------------
      !Back subtitution (WINDS)
      !------------------------
      else if (F_action_S=="BAC") then

         if (Williamson_case==1) then

            call wil_uvcase1 (ut0,vt0,l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.,Lctl_step)

            return

         end if

         if (Dcmip_case>=11.and.Dcmip_case<=13) then

             if (Dcmip_case==11) call dcmip_tracers11_transport (ut0,vt0,wt0,zdt0,bidon,bidon,bidon,bidon, &
                                                                 bidon,bidon,bidon,bidon,bidon,            &
                                                                 l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

             if (Dcmip_case==12) call dcmip_tracers12_transport (ut0,vt0,wt0,zdt0,bidon,bidon,bidon,bidon, &
                                                                 bidon,bidon,                              &
                                                                 l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

             if (Dcmip_case==13) call dcmip_tracers13_transport (ut0,vt0,wt0,zdt0,bidon,bidon,bidon,bidon, &
                                                                 bidon,bidon,bidon,bidon,bidon,            &
                                                                 l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)
             return

         end if

      !-------
      !Physics
      !-------
      else if (F_action_S=="PHY") then

         if (Dcmip_prec_type/=-1.or.Dcmip_pbl_type/=-1) call dcmip_2016_physics()

         if (Terminator_L) call canonical_Terminator()

      !------------------
      !Vertical diffusion
      !------------------
      else if (F_action_S=="VRD") then

         if (Dcmip_case>0.and.Dcmip_vrd_L) call dcmip_vrd_main()

      !---------------------------
      !Evaluate Errors/Diagnostics
      !---------------------------
      else if (F_action_S=="ERR") then

         if (Williamson_case==1) call wil_diagnostics (Lctl_step)
         if (Williamson_case==2) call wil_diagnostics (Lctl_step)

         if (Dcmip_case>0) call dcmip_diagnostics (Lctl_step)

      !-------------------------------------------
      !Output dependent variables on standard file
      !-------------------------------------------
      else if (F_action_S=="OUT") then

         if (Williamson_case==1) then

            do n=1,Tr3d_ntr

               if (.NOT.((Tr3d_name_S(n)(1:2)=='Q1').or. &
                         (Tr3d_name_S(n)(1:2)=='Q2').or. &
                         (Tr3d_name_S(n)(1:2)=='Q3').or. &
                         (Tr3d_name_S(n)(1:2)=='Q4'))) cycle

               if (Lctl_step/=Step_total    .and.(Williamson_NAIR==1.or.Williamson_NAIR==2)) cycle

               if (Tr3d_name_S(n)(1:2)/='Q1'.and.(Williamson_NAIR==0.or.Williamson_NAIR==3)) cycle

               tr => tracers_P(n)%pntr

               if (Tr3d_name_S(n)(1:2)=='Q1') tr_r => q1ref
               if (Tr3d_name_S(n)(1:2)=='Q2') tr_r => q2ref
               if (Tr3d_name_S(n)(1:2)=='Q3') tr_r => q3ref
               if (Tr3d_name_S(n)(1:2)=='Q4') tr_r => q4ref

               if (Tr3d_name_S(n)(1:2)=='Q1') tr_e => q1err
               if (Tr3d_name_S(n)(1:2)=='Q2') tr_e => q2err
               if (Tr3d_name_S(n)(1:2)=='Q3') tr_e => q3err
               if (Tr3d_name_S(n)(1:2)=='Q4') tr_e => q4err

               !Initialize REFERENCE at TIME>0
               !------------------------------
               if (Williamson_Nair==0) call wil_case1(tr_r,l_minx,l_maxx,l_miny,l_maxy,G_nk,0,Lctl_step)
               if (Williamson_Nair==3) call wil_case1(tr_r,l_minx,l_maxx,l_miny,l_maxy,G_nk,5,Lctl_step)

               !Initialize ERROR
               !----------------
               tr_e(1:l_ni,1:l_nj,1:G_nk) = tr(1:l_ni,1:l_nj,1:G_nk) - tr_r(1:l_ni,1:l_nj,1:G_nk)

            end do

         end if

         if (Terminator_L) then

            istat = tr_get('CL:P',cl)
            istat = tr_get('CL2:P',cl2)

            !Initialize CLY
            !--------------
            cly(1:l_ni,1:l_nj,1:G_nk) = cl(1:l_ni,1:l_nj,1:G_nk) + 2.0d0 * cl2(1:l_ni,1:l_nj,1:G_nk)

            !Initialize CLY REFERENCE
            !------------------------
            clyref(1:l_ni,1:l_nj,1:G_nk) =  CLY_REF

            !Initialize CLY ERROR
            !--------------------
            clyerr(1:l_ni,1:l_nj,1:G_nk) = cly(1:l_ni,1:l_nj,1:G_nk) - clyref(1:l_ni,1:l_nj,1:G_nk)

            !Initialize Average Column Integrated of CL/CL2/CLY
            !--------------------------------------------------
            if (Schm_autobar_L) then
               acl (:,:) = cl (:,:,1)
               acl2(:,:) = cl2(:,:,1)
               acly(:,:) = cly(:,:,1)
            else
               call dcmip_avg_column_integrated (acl ,cl, l_minx,l_maxx,l_miny,l_maxy,G_nk)
               call dcmip_avg_column_integrated (acl2,cl2,l_minx,l_maxx,l_miny,l_maxy,G_nk)
               call dcmip_avg_column_integrated (acly,cly,l_minx,l_maxx,l_miny,l_maxy,G_nk)
            end if

         end if

         if (Dcmip_case>=11.and.Dcmip_case<=13) then

            do n=1,Tr3d_ntr

               if (.NOT.((Tr3d_name_S(n)(1:2)=='Q1').or. &
                         (Tr3d_name_S(n)(1:2)=='Q2').or. &
                         (Tr3d_name_S(n)(1:2)=='Q3').or. &
                         (Tr3d_name_S(n)(1:2)=='Q4'))) cycle

               if (Tr3d_name_S(n)(1:2)/='Q1'.and.Dcmip_case==12) cycle

               if (Lctl_step/=Step_total) cycle

               tr => tracers_P(n)%pntr

               if (Tr3d_name_S(n)(1:2)=='Q1') tr_r => q1ref
               if (Tr3d_name_S(n)(1:2)=='Q2') tr_r => q2ref
               if (Tr3d_name_S(n)(1:2)=='Q3') tr_r => q3ref
               if (Tr3d_name_S(n)(1:2)=='Q4') tr_r => q4ref

               if (Tr3d_name_S(n)(1:2)=='Q1') tr_e => q1err
               if (Tr3d_name_S(n)(1:2)=='Q2') tr_e => q2err
               if (Tr3d_name_S(n)(1:2)=='Q3') tr_e => q3err
               if (Tr3d_name_S(n)(1:2)=='Q4') tr_e => q4err

               !Initialize ERROR
               !----------------
               tr_e(1:l_ni,1:l_nj,1:G_nk) = tr(1:l_ni,1:l_nj,1:G_nk) - tr_r(1:l_ni,1:l_nj,1:G_nk)

            end do

         end if

         !Prepare Perturbation of real potential Temperature
         !--------------------------------------------------
         if (Dcmip_case==31.or.Dcmip_case==163) then

            hu => tracers_P(Tr3d_hu)%pntr

            do k=1,G_nk
               do j=1,l_nj
                  do i=1,l_ni

                     !Real potential temperature
                     !--------------------------
                     pr_8 = pw_pt_plus(i,j,k)

                     th(i,j,k) = (tt1(i,j,k) / (1.d0 + 0.608d0 * hu(i,j,k))) * (Cstv_pref_8/pr_8) ** (rgasd_8/cpd_8)

                     pth(i,j,k) = th(i,j,k) - thbase(i,j,k)

                  end do
               end do
            end do

         end if

         !Prepare Perturbation of Virtual Temperature
         !-------------------------------------------
         dtv(1:l_ni,1:l_nj,1:G_nk) = tt1(1:l_ni,1:l_nj,1:G_nk) - Cstv_Tstr_8

      else

         call handle_error(-1,'CANONICAL_CASES','F_action_S unknown')

      end if
!
!---------------------------------------------------------------------
!
      return

 1005 format (/'STAGGERED VERTICAL LAYERING ON',I4,' MOMENTUM HYBRID LEVELS WITH ', &
               'Hyb_rcoef= ',4f7.2,':'/ &
               2x,'level',10x,'HYB',2x,'~dcmip_HEIGHTS',2x,'~dcmip_DELTA_Z',2x,'IP1')
 1006 format (1x,i4,3x,es16.4,2(6x,f6.0),4x,i10)

      end
