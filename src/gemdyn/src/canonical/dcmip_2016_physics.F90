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

!**s/r dcmip_2016_physics - DCMIP 2016 physics

      subroutine dcmip_2016_physics ()

      use canonical
      use cstv
      use DCMIP_2016_physics_module
      use dcmip_options
      use dyn_fisl_options
      use dynkernel_options
      use geomh
      use glb_ld
      use rmn_gmm
      use gmm_pw
      use gmm_vt1
      use gmm_phy
      use lun
      use mem_tracers
      use ptopo
      use step_options

      implicit none

#include <arch_specific.hf>

      !object
      !===============================================================|
      !  DCMIP 2016 physics                                           |
      !===============================================================|
      !  prec_type  | Type of precipitation/microphysics              |
      !             | ------------------------------------------------|
      !             |  0: Kessler Microphysics                        |
      !             |  1: Reed-Jablonowski Large-scale precipitation  |
      !             | -1: NONE                                        |
      !---------------------------------------------------------------|
      !  pbl_type   | Type of planetary boundary layer                |
      !             | ------------------------------------------------|
      !             |  0: Reed-Jablonowski Boundary layer             |
      !             |  1: Georges Bryan Boundary Layer                |
      !             | -1: NONE                                        |
      !---------------------------------------------------------------|

      !------------------------------------------------------------------------------!
      !Subroutines DCMIP 2016 physics (https://doi.org/10.5281/zenodo.1298671)       !
      !------------------------------------------------------------------------------!
      !GEM's adaptation to DCMIP 2016 physics:                                       !
      !1) Combines elements from 2 subroutines: simple_physics_v6/dcmip_physics_z_v1 !
      !2) PRESSURE not changed ((RHO is changed)                                     !
      !3) Vertical index increases from BOTTOM to TOP                                !
      !------------------------------------------------------------------------------!

      real(kind=REAL64) :: uu(l_nk),vv(l_nk),qsv(l_nk),qsc(l_nk),qsr(l_nk),precl,pm(0:l_nk),pt(l_nk+1),tt(l_nk)

      real(kind=REAL64) :: not_rotated_lat(l_ni,l_nj),lat,rlon_8,s_8(2,2),x_a_8,y_a_8

      integer :: i,j,k,istat,kk,step_reset,test

      real    :: period_4

      real, pointer, dimension (:,:,:) :: qsv_p,qsc_p,qsr_p,ptr3d

      real, dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk) :: tdu,tdv,tv_plus,pw_uu_plus0,pw_vv_plus0,pw_tt_plus0
!
!---------------------------------------------------------------------
!
      if (Lun_out>0) write (Lun_out,1000)

      !----------------------------------------------------!
      !INPUT/OUTPUT variables for DCMIP_2016 physics       !
      !----------------------------------------------------!
      !test      (IN) DCMIP2016 test id of SST (1,2,3)     !
      !uu     (INOUT) Zonal velocity (m/s)                 !
      !vv     (INOUT) Meridional velocity (m/s)            !
      !pt        (IN) Pressure (Pa) !THERMO                !
      !pm        (IN) Pressure (Pa) !MOMENTUM              !
      !qsv    (INOUT) Specific humidity (kg/kg)            !
      !qsc    (INOUT) Cloud water specific (kg/kg)         !
      !qsr    (INOUT) Rain water specific (kg/kg)          !
      !tt     (INOUT) Temperature (K)                      !
      !dt        (IN) Time step (s)                        !
      !lat       (IN) Latitude of column (radians)         !
      !l_nk      (IN) Number of levels in the column       !
      !precl    (OUT) Large-scale precip rate (m/s)        !
      !pbl_type  (IN) Type of planetary boundary layer     !
      !prec_type (IN) Type of precipitation/microphysics   !
      !----------------------------------------------------!

      !Copy Scalar winds (Required in itf_phy_UVupdate)
      !------------------------------------------------
      call itf_phy_copy ()

      !Averaging period for Averaged precipitation rate
      !------------------------------------------------
      period_4 = Step_total*Cstv_dt_8

      if (Dcmip_case==161) period_4 = 24.*3600.
      if (Dcmip_case==162) period_4 = 6.*3600.
      if (Dcmip_case==163) period_4 = 300.

      step_reset = period_4/Cstv_dt_8

      test = -1

      !Set Sea Surface Temperature (SST)
      !---------------------------------
      if (Dcmip_case== 43) test = 1
      if (Dcmip_case==161) test = 1
      if (Dcmip_case==162) test = 2
      if (Dcmip_case==163) test = 3

      if (test == -1) call handle_error(-1,'DCMIP_2016_PHYSICS','SST is not prescribed')

      !Get GMM variables
      !-----------------
      !pw_uu_plus !U wind component on Scalar grid
      !pw_vv_plus !V wind component on Scalar grid
      !pw_pm_plus !Pressure MOMENTUM
      !pw_pt_plus !Pressure THERMO
      !pw_p0_plus !Surface pressure
      !pw_tt_plus !Real temperature

      istat = tr_get('HU:P',qsv_p) !Specific humidity
      istat = tr_get('QC:P',qsc_p) !Cloud water specific
      istat = tr_get('RW:P',qsr_p) !Rain water specific

      !irt !Instantaneous precipitation rate
      !art !Averaged precipitation rate
      !wrt !Averaged precipitation rate (WORK FIELD)

      !Retrieve a copy of the PW state before the physics
      !--------------------------------------------------
      nullify(ptr3d)

      ptr3d => pw_uu_plus ; pw_uu_plus0 = ptr3d
      ptr3d => pw_vv_plus ; pw_vv_plus0 = ptr3d
      ptr3d => pw_tt_plus ; pw_tt_plus0 = ptr3d

      !--------------------------------
      !Evaluate latitudes (Not rotated)
      !--------------------------------
      do j = 1, l_nj

         lat   = geomh_y_8(j)
         y_a_8 = geomh_y_8(j)

         if (Ptopo_couleur == 0) then

             do i = 1, l_ni

                not_rotated_lat(i,j) = lat

             end do

         else

             do i = 1, l_ni

                x_a_8 = geomh_x_8(i) - acos(-1.D0)

                call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                not_rotated_lat(i,j) = lat

             end do

         end if

      end do

      !------------------------------------------------------------------------------!
      !DCMIP 2016 physics                                                            !
      !------------------------------------------------------------------------------!
      !Subroutines DCMIP 2016 physics (https://doi.org/10.5281/zenodo.1298671)       !
      !------------------------------------------------------------------------------!
      !GEM's adaptation to DCMIP 2016 physics:                                       !
      !1) Combines elements from 2 subroutines: simple_physics_v6/dcmip_physics_z_v1 !
      !2) PRESSURE not changed ((RHO is changed)                                     !
      !3) Vertical index increases from BOTTOM to TOP                                !
      !------------------------------------------------------------------------------!
      do j = 1,l_nj

         do i = 1,l_ni

            do k = 1,l_nk !Reverse TOP/BOTTOM

               kk = l_nk-k+1

               uu(kk) = pw_uu_plus(i,j,k)  !U wind component on Scalar grid
               vv(kk) = pw_vv_plus(i,j,k)  !V wind component on Scalar grid

               pt(kk) = pw_pt_plus(i,j,k)  !Pressure THERMO
               pm(kk) = pw_pm_plus(i,j,k)  !Pressure MOMENTUM

               tt(kk) = pw_tt_plus(i,j,k)  !Real temperature

              qsv(kk) = qsv_p(i,j,k)       !Specific Humidity
              qsc(kk) = qsc_p(i,j,k)       !Cloud water specific
              qsr(kk) = qsr_p(i,j,k)       !Rain water specific

            end do

            !Surface pressure
            !----------------
            pm(0) = pw_p0_plus(i,j)

            !Top (ESTIMATION)
            !----------------
            pt(l_nk+1) = Cstv_ptop_8

            call DCMIP2016_PHYSICS (test, uu, vv, pt, pm, qsv, qsc, qsr, tt,      &
                                    Cstv_dt_8, not_rotated_lat(i,j), l_nk, precl, &
                                    Dcmip_pbl_type, Dcmip_prec_type)

            do k = 1,l_nk !Reverse TOP/BOTTOM

               kk = l_nk-k+1

               pw_uu_plus(i,j,k) = uu(kk)  !U wind component on Scalar grid
               pw_vv_plus(i,j,k) = vv(kk)  !V wind component on Scalar grid

               pw_tt_plus(i,j,k) = tt(kk)  !Real temperature

                    qsv_p(i,j,k) = qsv(kk) !Specific Humidity
                    qsc_p(i,j,k) = qsc(kk) !Cloud water specific
                    qsr_p(i,j,k) = qsr(kk) !Rain water specific

            end do

            !Instantaneous precipitation rate (m/s):
            !Kessler Microphysics/Reed-Jablonowski Large-scale precipitation
            !---------------------------------------------------------------
            if (Dcmip_prec_type >= 0) then

               irt(i,j) = precl

            else

               irt(i,j) = 0.

            end if

            wrt(i,j) = wrt(i,j) + irt(i,j) !Accumulation (m/s) for Averaged precipitation rate

         end do

      end do

      !------------------------------------
      !Finalize Averaged precipitation rate
      !------------------------------------
      if (mod(Lctl_step,step_reset)==0) then

         art = wrt/float(step_reset) !Averaged precipitation rate (m/s)
         wrt = 0.

      end if

      call glbstat (irt,'IRT','LCPR', &
                     l_minx,l_maxx,l_miny,l_maxy,1,1,1,G_ni,1,G_nj,1,1)

      call glbstat (art,'ART','LCPR', &
                     l_minx,l_maxx,l_miny,l_maxy,1,1,1,G_ni,1,G_nj,1,1)

      istat = gmm_get(gmmk_phy_uu_tend_s,phy_uu_tend)
      istat = gmm_get(gmmk_phy_vv_tend_s,phy_vv_tend)
      istat = gmm_get(gmmk_phy_tv_tend_s,phy_tv_tend)

      !-------------------------------------------------------------------------------
      !Compute tendencies and reset physical world if requested (As in itf_phy_update)
      !-------------------------------------------------------------------------------
      if (Schm_phycpl_S == 'RHS' .or. Schm_phycpl_S == 'AVG') then

         tdu(1:l_ni,1:l_nj,1:l_nk) = pw_uu_plus(1:l_ni,1:l_nj,1:l_nk) - pw_uu_plus0(1:l_ni,1:l_nj,1:l_nk)
         tdv(1:l_ni,1:l_nj,1:l_nk) = pw_vv_plus(1:l_ni,1:l_nj,1:l_nk) - pw_vv_plus0(1:l_ni,1:l_nj,1:l_nk)
         call hwnd_stag(phy_uu_tend,phy_vv_tend,tdu,tdv,l_minx,l_maxx,l_miny,l_maxy,l_nk,.true.)

         call tt2virt (tv_plus,.true.,l_minx,l_maxx,l_miny,l_maxy,l_nk)
         phy_tv_tend(1:l_ni,1:l_nj,1:l_nk) = tv_plus(1:l_ni,1:l_nj,1:l_nk) - tt1(1:l_ni,1:l_nj,1:l_nk)

         phy_uu_tend = phy_uu_tend/Cstv_dt_8
         phy_vv_tend = phy_vv_tend/Cstv_dt_8
         phy_tv_tend = phy_tv_tend/Cstv_dt_8

         do k=1,l_nk
            pw_uu_plus(:,:,k) = pw_uu_plus0(:,:,k) + phy_cplm(:,:)*(pw_uu_plus(:,:,k)-pw_uu_plus0(:,:,k))
            pw_vv_plus(:,:,k) = pw_vv_plus0(:,:,k) + phy_cplm(:,:)*(pw_vv_plus(:,:,k)-pw_vv_plus0(:,:,k))
            pw_tt_plus(:,:,k) = pw_tt_plus0(:,:,k) + phy_cplt(:,:)*(pw_tt_plus(:,:,k)-pw_tt_plus0(:,:,k))
         end do

      else ! Pure split coupling

         phy_uu_tend = 0.
         phy_vv_tend = 0.
         phy_tv_tend = 0.

         if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then !For SPLIT in GEM-H
            call tt2virt (tv_plus,.true.,l_minx,l_maxx,l_miny,l_maxy,l_nk)
            phy_tv_tend(1:l_ni,1:l_nj,1:l_nk) = tv_plus(1:l_ni,1:l_nj,1:l_nk) - tt1(1:l_ni,1:l_nj,1:l_nk)
            phy_tv_tend = phy_tv_tend/Cstv_dt_8
         end if

      end if
!
!---------------------------------------------------------------------
!
      return

 1000 format( &
      /,'USE DCMIP 2016 PHYSICS : (S/R DCMIP_2016_physics)',/, &
        '=================================================')

      end subroutine dcmip_2016_physics
