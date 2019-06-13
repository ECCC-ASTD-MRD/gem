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

!**s/r gemdm_config - Establish final model configuration
!
#include <msg.h>

      integer function gemdm_config()
      use dynkernel_options
      use timestr_mod
      use mu_jdate_mod
      use step_options
      use HORgrid_options
      use fislh_sol
      use hvdif_options
      use init_options
      use gem_options
      use lam_options
      use out_options
      use VERgrid_options
      use dyn_fisl_options
      use ctrl
      use grdc_options
      use tdpack
      use glb_ld
      use cstv
      use lun
      use inp_mod
      use out_mod
      use out3
      use levels
      use rstr
      use wb_itf_mod
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>
#include <rmnlib_basics.hf>

      character(len=16) :: dumc_S, datev
      integer :: i, ipcode, ipkind, err
      real :: pcode,deg_2_rad,sec
      real(kind=REAL64) :: dayfrac
      real(kind=REAL64), parameter :: sec_in_day = 86400.0d0
!
!-------------------------------------------------------------------
!
      gemdm_config = -1

      if (Step_runstrt_S == 'NIL') then
         if (lun_out>0) then
            write (Lun_out, 6005)
            write (Lun_out, 8000)
         end if
         return
      end if

      if ( (Step_nesdt <= 0) .and. (Grd_typ_S(1:1) == 'L') &
                             .and. (.not. Lam_ctebcs_L ) ) then
         if (Lun_out > 0) then
            write(Lun_out,*)  ' Fcst_nesdt_S must be specified in namelist &step'
         end if
         return
      end if
      if (.not.Step_leapyears_L) then
         call Ignore_LeapYear ()
         call mu_set_leap_year (MU_JDATE_LEAP_IGNORED)
         if (Lun_out>0) write(Lun_out,6010)
      end if

      dayfrac= - (Step_initial * Step_dt / sec_in_day)
      call incdatsd (datev, Step_runstrt_S, dayfrac)
      call datp2f ( Out3_date, datev )
      if (lun_out>0) write (Lun_out,6007) Step_runstrt_S, datev
      call datp2f ( Step_CMCdate0, Step_runstrt_S)

      Lun_debug_L = (Lctl_debug_L.and.Ptopo_myproc == 0)

!     Use newcode style:
      call convip ( ipcode, pcode, ipkind, 0, ' ', .false. )

      if( trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then
         if ( trim(sol_type_S) == 'DIRECT') then
            isol_d=1.0d0
            isol_i=0.0d0
         else
            isol_i=1.0d0
            isol_d=0.0d0
         end if
      end if

      if (Grd_yinyang_L) then
         Lam_blend_H  = 0
         Lam_ctebcs_L = .true.
         Lam_blendoro_L = .false.
      else
         Lam_blend_H = max(0,Lam_blend_H)
      end if
      Lam_blend_Hx = Lam_blend_H ; Lam_blend_Hy = Lam_blend_H
      if (Ctrl_theoc_L) Lam_blend_Hy = 0
      Lam_tdeb = 999999.0

      if (.not.Grd_yinyang_L) then
         Lam_gbpil_T = max(0,Lam_gbpil_T)
         Lam_blend_T = max(0,Lam_blend_T)
      else
         Lam_gbpil_T = 0
         Lam_blend_T = 0
      end if

      deg_2_rad = pi_8/180.

      P_lmvd_high_lat = min(90.,abs(P_lmvd_high_lat))
      P_lmvd_low_lat  = min(P_lmvd_high_lat,abs(P_lmvd_low_lat))
      P_lmvd_weigh_low_lat  = max(0.,min(1.,P_lmvd_weigh_low_lat ))
      P_lmvd_weigh_high_lat = max(0.,min(1.,P_lmvd_weigh_high_lat))

      Out_rewrit_L   = .false.

      if (Out3_close_interval_S == " ") Out3_close_interval_S= '1H'
      err = timestr_parse(Out3_close_interval,Out3_unit_S,Out3_close_interval_S)
      if (Out3_close_interval <= 0.) Out3_close_interval_S= '1H'
      err = timestr_parse(Out3_close_interval,Out3_unit_S,Out3_close_interval_S)
      if (.not.RMN_IS_OK(err)) then
         if (lun_out>0) write (Lun_out, *) "(gemdm_config) Invalid Out3_close_interval_S="//trim(Out3_close_interval_S)
         return
      end if
      Out3_postproc_fact = max(0,Out3_postproc_fact)

      if(Out3_nbitg < 0) then
         if (lun_out>0) write (Lun_out, 9154)
         Out3_nbitg=16
      end if
      Out3_lieb_nk = 0
      do i = 1, MAXELEM
         if (Out3_lieb_levels(i) <= 0.) EXIT
         Out3_lieb_nk = Out3_lieb_nk + 1
      end do

      Cstv_dt_8 = dble(Step_dt)

      if (Schm_autobar_L) Ctrl_phyms_L = .false.

      Dynamics_hauteur_L = Dynamics_Kernel_S(14:15) == '_H'

      if (Dynamics_hauteur_L) then
         Dynamics_hydro_L = .false.
         call fislh_hybrid ( hyb_H, G_nk)
      else
         call set_zeta ( hyb, G_nk )
      endif

      Schm_nith = G_nk

      err = 0
      err = min( timestr2step (Grdc_ndt,   Grdc_nfe    , Step_dt), err)
      err = min( timestr2step (Grdc_start, Grdc_start_S, Step_dt), err)
      err = min( timestr2step (Grdc_end,   Grdc_end_S  , Step_dt), err)
      if (err < 0) return

      if (Grdc_ndt   > -1) Grdc_ndt  = max( 1, Grdc_ndt)
      if (Grdc_start <  0) Grdc_start= Lctl_step
      if (Grdc_end   <  0) Grdc_end  = Step_total

      Grdc_maxcfl = max(1,Grdc_maxcfl)

      call low2up  (Lam_hint_S ,dumc_S)
      Lam_hint_S = dumc_S

      call low2up  (sol_type_S ,dumc_S)
      sol_type_S = trim(dumc_S)
      call low2up  (sol2D_precond_S ,dumc_S)
      sol2D_precond_S = trim(dumc_S)
      call low2up  (sol3D_precond_S ,dumc_S)
      sol3D_precond_S = trim(dumc_S)

      if ( (Sol_type_S /= 'ITERATIVE_2D') .and. &
           (Sol_type_S /= 'ITERATIVE_3D') .and. &
           (Sol_type_S /= 'DIRECT'   ) ) then
         if (lun_out>0) write (Lun_out, 9200) Sol_type_S
         return
      end if

      if (Sol_type_S(1:9) == 'ITERATIVE') then
         if (Sol2D_precond_S /= 'JACOBI' .and. Sol2D_precond_S /= 'REDBLACK') then
            if (lun_out>0) write (Lun_out, 9201) Sol2D_precond_S
            return
         end if
      end if

      if (Sol_type_S == 'ITERATIVE_3D') then
         if (Sol3D_krylov_S /= 'FGMRES' .and. Sol3D_krylov_S /= 'FBICGSTAB') then
            if (Lun_out > 0) then
               write(Lun_out, *) 'ABORT: WRONG CHOICE OF KRYLOV METHOD FOR 3D ITERATIVE SOLVER: Sol3D_krylov_S =', Sol3D_krylov_S
            end if
            return
         end if

         if (Sol3D_precond_S /= 'JACOBI'          &
             .and. Sol3D_precond_S /= 'REDBLACK'  &
             .and. Sol3D_precond_S /= 'MULTIGRID' &
             .and. Sol3D_precond_S /= 'GAUSS'     &
            ) then
            if (lun_out>0) then
               write (Lun_out, *) 'ABORT: WRONG CHOICE OF PRE-CONDITIONER FOR 3D ITERATIVE SOLVER: Sol3D_precond_S =', Sol3D_precond_S
            end if
            return
         end if
      end if

      if (Hzd_smago_param > 0. .or. Hzd_smago_lnr(2) > 0.) then
         if (Hzd_smago_lnr(1) > Hzd_smago_lnr(2) .or. Hzd_smago_lnr(1)<0.) Hzd_smago_lnr(1)=Hzd_smago_lnr(2)
         if (Hzd_smago_lnr(3) < Hzd_smago_lnr(2) .or. Hzd_smago_lnr(3)<0.) Hzd_smago_lnr(3)=Hzd_smago_lnr(2)
      end if

      G_ni  = Grd_ni
      G_nj  = Grd_nj

      G_niu = G_ni
      G_njv = G_nj - 1

      G_niu = G_ni - 1

      if ( Schm_psadj<0 .or. Schm_psadj>2 ) then
         if (lun_out>0) write (Lun_out, 9700)
         return
      end if

!     Check for open top (piloting and blending)
      Schm_opentop_L  = .false.
      if (Lam_gbpil_T > 0) then
         Schm_opentop_L = .true.
         if(lun_out>0.and.Vspng_nk/=0)write (Lun_out, 9570)
         Vspng_nk       = 0
      end if
      if (Lam_gbpil_T <= 0 .and. Lam_blend_T > 0) then
         if (lun_out>0)write (Lun_out, 9580)
         return
      end if


!     Check for modified epsilon/ super epsilon
      if ( Cstv_rE_8 /= 1.0d0 ) then
         if (Cstv_Tstr_8 < 0. ) then
             if (lun_out>0) write (Lun_out, 9680) 'Cstv_rE_8'
            return
         end if

         if (Schm_opentop_L) then
            if (lun_out>0) write (Lun_out, 9681) 'Cstv_rE_8'
            return
         end if
      end if

      if ( Cstv_bA_8 < Cstv_bA_m_8)  Cstv_bA_m_8=Cstv_bA_8
      if ( Cstv_bA_nh_8 < Cstv_bA_8) Cstv_bA_nh_8=Cstv_bA_8

!     Check for incompatible use of IAU and DF
      if (Iau_period > 0. .and. Init_balgm_L) then
         if (Lun_out>0) then
            write (Lun_out, 7040)
            write (Lun_out, 8000)
         end if
         return
      end if

      Ctrl_testcases_L = Ctrl_canonical_dcmip_L .or. &
                         Ctrl_canonical_williamson_L

!     Some common setups for AUTOBAROTROPIC runs
      if (Schm_autobar_L) then
         Dynamics_hydro_L = .true.
         call wil_set (Schm_topo_L,Ctrl_testcases_adv_L,Lun_out,err)
         if (err < 0) return
         if (Lun_out>0) write (Lun_out, 6100) Schm_topo_L
      end if

      err= 0
      err= min( timestr2step (Init_dfnp, Init_dflength_S, Step_dt), err)
      err= min( timestr2sec  (sec,       Init_dfpl_S,     Step_dt), err)
      if (err < 0) return

      Init_halfspan = -9999

      if (Init_dfnp <= 0) then
         if (Init_balgm_L .and. Lun_out > 0) write(Lun_out,6602)
         Init_balgm_L = .false.
      end if
      Init_dfpl_8 = sec
      if (Init_balgm_L) then
         Init_dfnp   = max(2,Init_dfnp)
         if ( mod(Init_dfnp,2) /= 1 ) Init_dfnp= Init_dfnp+1
         Init_halfspan = (Init_dfnp - 1) / 2
         if (Step_total<Init_halfspan) Init_balgm_L= .false.
      end if
      if (.not.Init_balgm_L) then
         Init_dfnp     = -9999
         Init_halfspan = -9999
      end if
      if (.not.Rstri_rstn_L) then
         Init_mode_L= .false.
         if (Init_balgm_L) Init_mode_L= .true.
      end if

      if (Lctl_debug_L) then
         call msg_set_minMessageLevel(MSG_INFOPLUS)
!        istat = gmm_verbosity(GMM_MSG_DEBUG)
!        istat = wb_verbosity(WB_MSG_DEBUG)
!        call msg_set_minMessageLevel(MSG_DEBUG)
         call handle_error_setdebug(Lctl_debug_L)
      end if

      if (Hzd_lnr_theta < 0.) then
         Hzd_lnr_theta = Hzd_lnr
         if (lun_out>0) write (Lun_out, 6200)
      end if
      if (Hzd_pwr_theta < 0 ) then
         if ((lun_out>0).and.(Hzd_lnr_theta > 0.)) write (Lun_out, 6201)
         Hzd_pwr_theta = Hzd_pwr
      end if

      gemdm_config = 1

      if ( Vtopo_start_S  == '' ) then
         Vtopo_start = Step_initial
      else
         err= min( timestr2step (Vtopo_start,Vtopo_start_S,Step_dt),err)
      end if
      if ( Vtopo_length_S  == '' ) then
         Vtopo_ndt = 0
      else
         err= min( timestr2step (Vtopo_ndt,Vtopo_length_S,Step_dt),err)
      end if
      Vtopo_L = (Vtopo_ndt > 0)

 6005 format (/' Starting time Step_runstrt_S not specified'/)
 6007 format (/X,48('#'),/,2X,'STARTING DATE OF THIS RUN IS : ',a16,/ &
                           2X,'00H FORECAST WILL BE VALID AT: ',a16,/ &
               X,48('#')/)
 6010 format (//'  =============================='/&
                '  = Leap Years will be ignored ='/&
                '  =============================='//)
 6100 format (//'  =================================='/&
                '  AUTO-BAROTROPIC RUN: Schm_topo = ',L1,/&
                '  =================================='//)
 6200 format ( ' Applying DEFAULT FOR THETA HOR. DIFFUSION Hzd_lnr_theta= Hzd_lnr')
 6201 format ( ' Applying DEFAULT FOR THETA HOR. DIFFUSION Hzd_pwr_theta= Hzd_pwr')
 6602 format (/' WARNING: Init_dfnp is <= 0; Settings Init_balgm_L=.F.'/)
 7040 format (/' OPTION Init_balgm_L=.true. NOT AVAILABLE if applying IAU (Iau_period>0)'/)
 8000 format (/,'========= ABORT IN S/R GEMDM_CONFIG ============='/)
 9154 format (/,' Out3_nbitg IS NEGATIVE, VALUE will be set to 16'/)
 9200 format (/'ABORT: WRONG CHOICE OF SOLVER for Helmholtz problem: Sol_type_S =',a/)
 9201 format (/'ABORT: WRONG CHOICE OF PRE-CONDITIONER FOR 2D ITERATIVE SOLVER: Sol2D_precond_S =',a/)
 9570 format (/,'WARNING: Vspng_nk set to zero since top piloting is used'/)
 9580 format (/,'ABORT: Non zero Lam_blend_T cannot be used without top piloting'/)
 9680 format (/,'ABORT: ',a,' cannot be less than 1.0 for T*<0'/)
 9681 format (/,'ABORT: ',a,' cannot be less than 1.0 for OPEN_TOP scheme'/)
 9700 format (/,'ABORT: Schm_psadj not valid'/)
!
!-------------------------------------------------------------------
!
      return
      end

