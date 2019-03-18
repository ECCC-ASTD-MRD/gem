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

      use DCMIP_2016_physics_module
      use vertical_interpolation, only: vertint2
      use step_options
      use gmm_pw
      use canonical
      use dcmip_options
      use gem_options
      use geomh
      use tdpack, only : rgasd_8, rgasv_8, cpd_8, grav_8
      use glb_ld
      use cstv
      use lun
      use ver
      use gmm_itf_mod
      use ptopo

      implicit none

#include <arch_specific.hf>

      !object
      !======================================================================|
      !  DCMIP 2016 physics                                                  |
      !======================================================================|
      !  prec_type         | Type of precipitation/microphysics              |
      !                    | ------------------------------------------------|
      !                    |  0: Large-scale precipitation (Kessler)         |
      !                    |  1: Large-scale precipitation (Reed-Jablonowski)|
      !                    | -1: NONE                                        |
      !----------------------------------------------------------------------|
      !  pbl_type          | Type of planetary boundary layer                |
      !                    | ------------------------------------------------|
      !                    |  0: Reed-Jablonowski Boundary layer             |
      !                    |  1: Georges Bryan Planetary Boundary Layer      |
      !                    | -1: NONE                                        |
      !----------------------------------------------------------------------|

      !------------------------------------------------------------------------------------------------------------

      real(8)  :: uu(l_nk),vv(l_nk),qsv(l_nk),qv(l_nk),qsc(l_nk),qsr(l_nk),rho(l_nk),zm(0:l_nk),zt(0:l_nk),precl, &
                  pm(0:l_nk),pt(l_nk+1),tv(l_nk),tt(l_nk),wm(l_nk),wp(l_nk),theta(l_nk),exner(l_nk), &
                  dudt(l_nk),dvdt(l_nk),pdel_mid(l_nk)

      real(8)  :: dlnpint,zvir

      real(8) not_rotated_lat(l_ni,l_nj),lat,rlon_8,s_8(2,2),x_a_8,y_a_8

      integer i,j,k,istat,kk,step_reset,test

      real, pointer, dimension (:,:,:) :: qsv_p,qsc_p,qsr_p

      real  pt_p(l_minx:l_maxx,l_miny:l_maxy,l_nk+1),    pm_p(l_minx:l_maxx,l_miny:l_maxy,l_nk+1), &
            uu_p(l_minx:l_maxx,l_miny:l_maxy,l_nk)  ,    vv_p(l_minx:l_maxx,l_miny:l_maxy,l_nk)  , &
           gzt_p(l_minx:l_maxx,l_miny:l_maxy,l_nk+1),   gzm_p(l_minx:l_maxx,l_miny:l_maxy,l_nk+1), &
        log_pt_p(l_minx:l_maxx,l_miny:l_maxy,l_nk+1),log_pm_p(l_minx:l_maxx,l_miny:l_maxy,l_nk+1), &
            p0_p(l_minx:l_maxx,l_miny:l_maxy)

      real period_4

      real(8), parameter :: Rd   = Rgasd_8, & ! cte gaz - air sec   [J kg-1 K-1]
                            Rv   = Rgasv_8, & ! cte gaz - vap eau   [J kg-1 K-1]
                            cp_  = cpd_8,   & ! chal. spec. air sec [J kg-1 K-1]
                            grav = grav_8     ! acc. de gravite     [m s-2]

      !------------------------------------------------------------------------------------------------------------

      if (Lun_out>0) write (Lun_out,1000)

      !------------------------------------------------------------------------------
      !Preparation DCMIP physics
      !------------------------------------------------------------------------------
      !     u      - zonal wind at model levels (m/s)
      !     v      - meridional wind at model levels (m/s)
      !     qv     - water vapor mixing ratio (gm/gm) !WITH RESPECT TO DRY AIR
      !     qc     - cloud water mixing ratio (gm/gm) !WITH RESPECT TO DRY AIR
      !     qr     - rain  water mixing ratio (gm/gm) !WITH RESPECT TO DRY AIR
      !     rho    - dry air density (not mean state as in KW) (kg/m^3)
      !     dt     - time step (s)
      !     zm     - heights of MOMENTUM levels in the grid column (m)
      !     zt     - heights of THERMO levels in the grid column (m)
      !     pm     - pressure of MOMENTUM levels in the grid column (Pa)
      !     pt     - pressure of THERMO levels in the grid column (Pa)
      !     tv     - virtual temperature (K)
      !     wm     - weight Thermo to Momentum (above M)
      !     wp     - weight Thermo to Momentum (below M)
      !     l_nk   - number of levels in the column
      !     precl  - large-scale precip rate (m/s)
      !------------------------------------------------------------------------------

      zvir = (Rv/Rd) - 1 ! Constant for virtual temp. calc. is approx. 0.608

      !Copy Scalar winds
      !-----------------
      call itf_phy_copy ()

      !Averaging period for Averaged precipitation rate
      !------------------------------------------------
      period_4 = Step_total*Cstv_dt_8

      if (Dcmip_case==161) period_4 = 24.*3600.
      if (Dcmip_case==162) period_4 = 6.*3600.
      if (Dcmip_case==163) period_4 = 300.

      step_reset = period_4/Cstv_dt_8

      test = -1

      !Set Sea surface temperature
      !---------------------------
      if (Dcmip_case== 43) test = 1
      if (Dcmip_case==161) test = 1
      if (Dcmip_case==162) test = 2
      if (Dcmip_case==163) test = 3

      if (test == -1) call handle_error(-1,'DCMIP_2016_PHYSICS','SST is not prescribed')

      !Get VMM variables
      !-----------------
      istat = gmm_get(gmmk_pw_uu_plus_s , pw_uu_plus )
      istat = gmm_get(gmmk_pw_vv_plus_s , pw_vv_plus )
      istat = gmm_get(gmmk_pw_pm_plus_s , pw_pm_plus )
      istat = gmm_get(gmmk_pw_pt_plus_s , pw_pt_plus )
      istat = gmm_get(gmmk_pw_p0_plus_s , pw_p0_plus )
      istat = gmm_get(gmmk_pw_gz_plus_s , pw_gz_plus )
      istat = gmm_get(gmmk_pw_log_pm_s  , pw_log_pm  )
      istat = gmm_get(gmmk_pw_log_pt_s  , pw_log_pt  )
      istat = gmm_get(gmmk_pw_me_plus_s , pw_me_plus )
      istat = gmm_get(gmmk_pw_tt_plus_s , pw_tt_plus )

      istat = gmm_get('TR/'//'HU'//':P',qsv_p) !Specific humidity
      istat = gmm_get('TR/'//'QC'//':P',qsc_p) !Cloud water specific
      istat = gmm_get('TR/'//'RW'//':P',qsr_p) !Rain water specific

      istat = gmm_get(gmmk_irt_s,  irt)
      istat = gmm_get(gmmk_art_s,  art)
      istat = gmm_get(gmmk_wrt_s,  wrt)

      !Get Scalar winds
      !----------------
      uu_p(:,:,:) = pw_uu_plus(:,:,:)
      vv_p(:,:,:) = pw_vv_plus(:,:,:)

      !Get hydrostatic pressure on Momentum/Thermo levels
      !--------------------------------------------------
      pm_p(:,:,1:l_nk) = pw_pm_plus(:,:,1:l_nk)
      pt_p(:,:,1:l_nk) = pw_pt_plus(:,:,1:l_nk)
      p0_p(:,:)        = pw_p0_plus(:,:)

      pm_p(1:l_ni,1:l_nj,l_nk+1) = p0_p(1:l_ni,1:l_nj)
      pt_p(1:l_ni,1:l_nj,l_nk+1) = p0_p(1:l_ni,1:l_nj)

      !Calculate geopotential on Momentum/Thermo levels
      !------------------------------------------------
      gzm_p(:,:,1:l_nk)   = pw_gz_plus(:,:,1:l_nk)
      gzm_p(:,:,l_nk+1)   = pw_me_plus(:,:)

      log_pm_p(:,:,1:l_nk+1) = pw_log_pm (:,:,1:l_nk+1)
      log_pt_p(:,:,1:l_nk+1) = pw_log_pt (:,:,1:l_nk+1)

      gzt_p(1:l_ni,1:l_nj,l_nk+1) = gzm_p(1:l_ni,1:l_nj,l_nk+1)

      call vertint2 ( gzt_p, log_pt_p, l_nk, gzm_p, log_pm_p, l_nk+1, &
                      l_minx,l_maxx,l_miny,l_maxy,1,l_ni,1,l_nj )

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

      !------------------
      !DCMIP 2016 physics
      !------------------
      do j = 1,l_nj

         do i = 1,l_ni

            do k = 1,l_nk !Reverse TOP/BOTTOM

               kk = l_nk-k+1

               uu(kk) = uu_p(i,j,k)                                     !U wind component on Scalar grid
               vv(kk) = vv_p(i,j,k)                                     !V wind component on Scalar grid

              qsv(kk) = qsv_p(i,j,k)                                    !Specific Humidity (gm/gm)
              qsc(kk) = qsc_p(i,j,k)                                    !Cloud water specific (gm/gm)
              qsr(kk) = qsr_p(i,j,k)                                    !Rain water specific (gm/gm)

               qv(kk) = qsv(kk)/(1.0d0 - qsv(kk))                       !Conversion Specific Humidity to
                                                                        !Water vapor mixing ratio (gm/gm)

               pt(kk) = pt_p(i,j,k)                                     !Pressure on thermodynamic levels
               pm(kk) = pm_p(i,j,k)                                     !Pressure on momentum      levels

            end do

            !Surface pressure
            !----------------
            pm(0) = pm_p(i,j,l_nk+1)

            !Top (ESTIMATION)
            !----------------
            pt(l_nk+1) = Cstv_ptop_8

            do k = 1,l_nk !Reverse TOP/BOTTOM

               kk = l_nk-k+1

               !Coefficients for interpolating from Thermo to Momentum
               !------------------------------------------------------
               wp(kk) = (pm(kk) - pt(kk+1)) / (pt(kk) - pt(kk+1))        !Below M (LINEAR interpolation PRESSURE)
               wm(kk) = 1.0d0 - wp(kk)                                   !Above M (LINEAR interpolation PRESSURE)

               tt(kk) =  pw_tt_plus(i,j,k)                               !Real temperature
               tv(kk) =  (1.d0 + zvir * qsv(kk))*tt(kk)                  !Virtual temperature

               exner(kk) = (pt(kk)/Cstv_pref_8) ** (Rd/cp_)              !Exner pressure
                 rho(kk) =  pt(kk)/(Rd*tv(kk)*(1.d0 + qv(kk)))           !Dry air density (kg/m^3)
               theta(kk) =  tt(kk) / exner(kk)                           !Potential temperature

            pdel_mid(kk) =  pm(kk-1) - pm(kk)                            !For RJ large-scale prec.

            end do

            !------------------------------------------
            !Estimate heights of MOMENTUM/THERMO levels
            !------------------------------------------

               !Bottom
               !------
               zt(0) = 0. !ASSUME NO TOPOGRAPHY
               zm(0) = 0. !ASSUME NO TOPOGRAPHY

               !Estimate zm
               !-------------------------------------------------------------------
               !NOTE1: We did not use here DIAG_FI. It is another OPTION
               !NOTE2: The revised zm before PBL uses also the following estimation
               !-------------------------------------------------------------------
               do kk = 1,l_nk
                  dlnpint = log(pm(kk-1)) - log(pm(kk))
                  zm(kk)  = zm(kk-1) + Rd/grav*tt(kk)*(1.d0 + zvir * qsv(kk))*dlnpint
               end do

               !Estimate zt
               !-----------

                  do k = 1,l_nk !Reverse TOP/BOTTOM

                     kk = l_nk-k+1

                     zt(kk) = gzt_p(i,j,k)/grav

                  end do

            call DCMIP2016_PHYSICS (test, uu, vv, pt, pm, qsv, qsc, qsr, rho, theta, exner, tt, dudt, dvdt, pdel_mid, &
                                    Cstv_dt_8, zt, zm, not_rotated_lat(i,j), wm, wp, l_nk, precl, Dcmip_pbl_type, Dcmip_prec_type)

            !-------------------------------------------------------------------
            !                        UPDATE
            !-------------------------------------------------------------------
            !Note: PRESSURE is not changed (RHO is changed) in DCMIP2016_PHYSICS
            !-------------------------------------------------------------------

           do k = 1,l_nk !Reverse TOP/BOTTOM

              kk = l_nk-k+1

              pw_uu_plus(i,j,k) = uu(kk) !Scalar u
              pw_vv_plus(i,j,k) = vv(kk) !Scalar v

              qsv_p(i,j,k) = qsv(kk)
              qsc_p(i,j,k) = qsc(kk)
              qsr_p(i,j,k) = qsr(kk)

            end do

            do k = 1,l_nk !Reverse TOP/BOTTOM

               kk = l_nk-k+1

               pw_tt_plus(i,j,k) = tt(kk) !Real temperature

            end do

            !Instantaneous precipitation rate (m/s):
            !Kessler microphysics/Reed-Jablonowski large scale precipitation
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

      call glbstat2 (irt,'IRT','LCPR', &
                     l_minx,l_maxx,l_miny,l_maxy,1,1,1,G_ni,1,G_nj,1,1)

      call glbstat2 (art,'ART','LCPR', &
                     l_minx,l_maxx,l_miny,l_maxy,1,1,1,G_ni,1,G_nj,1,1)

      return

 1000 format( &
      /,'USE DCMIP 2016 PHYSICS : (S/R DCMIP_2016_physics)',/, &
        '=================================================')

      end subroutine dcmip_2016_physics
