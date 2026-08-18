!*****************************************************************************
!                Time-stamp: <2014-09-24 11:29:19 sander>
!*****************************************************************************

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program; if not, get it from:
! http://www.gnu.org/copyleft/gpl.html

!*****************************************************************************
! original source from Sander et al 2015 (doi:10.5194/gmd-7-2653-2014)
!-----------------
  ! 2016 Apr - Jack Chen implement for GEMMACH
  !   Modification to original code:
  !   - changed constants params. to be those from in GEMMACH
  !   ! USE messy_main_constants_mem, ONLY: N_A, R_gas, pi, g, k_B, M_air
  !   ! to  USE chm_consphychm_mod,   ONLY: avno, pi, grav
  !   - Combine/extract variables from original 3 modules into 1 'messy_module'
  !   - Convert all dynamic allocable array/pointers to static/autom. array
  !     (necessary for thread safty: OMP>1)
  !   - Expand "SUBROUTINE jvalues" to include "jval_cal_uv" and
  !     "#include "messy_jval_jvpp.inc" and remove all allocable arrays
  !     declaration from main module to static arrays in subroutine "jvalues"
  !   - add output arguments for SR:JVALUE (needed for OMP>1)
  !     INTENT(OUT) rh_o2_2d, rh_o3_2d, fhuv_2d, fhuvdna_2d, jval_2d
  !   - Change output pointer array to 3D static array in 'messy_jval_jvpp.inc'
  !     jval_2d(1:ip_max)%(:,:) => jval_2d(1:ip_out,:,:)
  !   - Add 'integer:ip_out', 'integer:jval_[xx]' and 'character:jvalname'
  !     in messy_module to index, label output array
  !   - Convert "pbllev" integer to 1D array input for all calculating grids
  !   - Add "dz2d" INTENT(IN) needed for part_col calculation.  This is now
  !     derived by grid height in memomentum levels (m)
  !   - shift proprocessor calls to be inside "jvalue"
  !        CALL messy_jval_init
  !        CALL aerosol_data( ... )
  !        CALL jval_solar_time_control( ... )
  !   - How To output additional Jx
  !      - edit "IP_OUT" integer for number of requested Jvalues
  !      - edit "lp(ip_out)" logical flag to turn on the required species
  !      - outputs are in "jval_2d(:,:,:)" array in sec-1
  !
  ! 2016 Jun - rename file/module
  !      "messy_jval.ftn90" -> "mach_messy_jval_mod.ftn90"
  !      "messy_jval_jvpp.inc" -> "mach_messy_jval_inc.ftn90"
  !      "messy_module.ftn90" -> "mach_pkg_messy_mod.ftn90"
  !      rename subroutine:
  !      "jvalues" -> 'mach_messy_jvalues"
!*****************************************************************************

MODULE mach_messy_jval_mod
USE chm_utils_mod,      only: chm_lun_out, global_debug
IMPLICIT NONE
CONTAINS
  ! **************************************************************************
  SUBROUTINE mach_messy_jvalues(v3_2d, cossza_1d, press_2d, rhum_2d, temp_2d, &
                                albedo_1d, aclc_2d, slf_1d, clp_2d, dz2d,     &
                                lmidatm, pbllev, iu0, kproma_day, jval_2d,    &
                                gni, gnk, relo3_2d, rh_o2_2d, rh_o3_2d,       &
                                fhuv_2d, fhuvdna_2d)

    ! Note that although relo3_2d has the dimension 1:klev+1, the value
    ! relo3_2d(1) is not used at all here. Also, note that relo3_2d is
    ! _only_ used for the heating rates. For the calculation of the
    ! J-values, only v3_2d is used.

    ! interface for calculation of J-values and heating rates
    ! (a) the input data are sorted in j = 1,kproma_day, so that
    !     for j = 1,kproma_day u0(j) > u0lim (daytime).
    ! (b) corrections is done for u0 because of
    !     spherical geometry of the earth.
    ! (c) for some interpolation temperature and pressure
    !     on a virtual level "0" is necessary
    !     (see e.g. press(j,0) and temp(j,0) in aero_2d).
    ! (d) aerosols are taken in to account very rudimentary
    !     in the lowest pbllev levels a mixture of
    !     rural and maritime aerosol is asumed
    ! (e) initialzation of J-values rates
    ! (f) resort the output of the band model to j = 1, chm_ni

    ! correction of the air mass factor
    ! F. Kasten and T. Young, Revised optical air mass tabels and
    ! approximation formula (1989) Applied Optics Vol 28, no. 22 p. 4735
    ! and J. Lenoble, Atmospheric Radiative Transfer (1993), p. 236

    USE mach_pkg_messy_mod
    USE chm_consphychm_mod  , ONLY: avno, pi, grav, kboltz

    IMPLICIT NONE

    !-------------------------------------------------------------------------
    ! I/O
    INTEGER,                                   INTENT(IN) :: gni, gnk
    REAL(SP), DIMENSION(gni, gnk+1),           INTENT(IN) :: v3_2d   ! ozone column
    REAL(SP), DIMENSION(gni),                  INTENT(IN) :: cossza_1d
    REAL(SP), DIMENSION(gni, gnk),             INTENT(IN) :: press_2d
    REAL(SP), DIMENSION(gni, gnk),             INTENT(IN) :: rhum_2d
    REAL(SP), DIMENSION(gni, gnk),             INTENT(IN) :: temp_2d
    REAL(SP), DIMENSION(gni),                  INTENT(IN) :: albedo_1d
    REAL(SP), DIMENSION(gni, gnk),             INTENT(IN) :: aclc_2d
    REAL(SP), DIMENSION(gni),                  INTENT(IN) :: slf_1d
    REAL(SP), DIMENSION(gni, gnk),             INTENT(IN) :: clp_2d
! JC new input layer thickness for particle concentration
    REAL(SP), DIMENSION(gni, gnk),             INTENT(IN) :: dz2d
    LOGICAL,                                   INTENT(IN) :: lmidatm
! JC change pbl from integer to array
!   INTEGER,              INTENT(IN) :: pbllev
    INTEGER, DIMENSION(gni),                   INTENT(IN) :: pbllev
    INTEGER, DIMENSION(gni),                   INTENT(IN) :: iu0
    INTEGER,                                   INTENT(IN) :: kproma_day
    REAL(DP), DIMENSION(MESSY_NJX, gni, gnk),  INTENT(OUT) :: jval_2d
    REAL(SP), DIMENSION(gni, gnk+1), OPTIONAL, INTENT (IN) :: relo3_2d
    REAL(DP), DIMENSION(gni, gnk), OPTIONAL,   INTENT(OUT) :: rh_o2_2d
    REAL(DP), DIMENSION(gni, gnk), OPTIONAL,   INTENT(OUT) :: rh_o3_2d
    REAL(DP), DIMENSION(gni, gnk), OPTIONAL,   INTENT(OUT) :: fhuv_2d
    REAL(DP), DIMENSION(gni, gnk), OPTIONAL,   INTENT(OUT) :: fhuvdna_2d
    !-------------------------------------------------------------------------
    ! LOCAL PARAMETERs
    ! <moved to messy_module.ftn90>
    !-------------------------------------------------------------------------

    ! LOCAL variables
    REAL, DIMENSION(kproma_day,0:gnk) :: &
      v2,    & ! O2 column density             [part./cm^2]
      v2s,   & ! diff. O2 column density       [part./cm^2]
      v3,    & ! O3 column density             [part./cm^2]
      v3s      ! diff. O3 column density       [part./cm^2]

    REAL, DIMENSION(kproma_day,gnk) :: &
      v2s2,     & ! v2s (interval 1, ma-corr)
      dh_O2,    & ! O2 heating rate                   [K/s]
      dh_o3,    & ! O3 heating rate                   [K/s]
      aclc,     & ! cloud fraction                    [1]
      clp,      & ! cloud liquid water path per layer [g/m^2]
      rhum,     & ! relative humudity                 [%]
      part_col    ! aerosol part. column per layer    [part/cm^2]

    REAL, DIMENSION(kproma_day) :: &
      slf, & ! fraction of sea surface
      u0     ! cosine of solar zenith angle

    REAL, DIMENSION(kproma_day,gnk,MAXWAV) :: &
      fact    ! actinic flux [photons/(cm^2 s)]
    REAL, DIMENSION(kproma_day,gnk,3:5) :: &
      facth   ! flux [photons/(cm^2 s)]

    REAL, DIMENSION(2) :: sig_h_o3, sig_h_o2 ! heating rates

    REAL, DIMENSION(kproma_day,MAXWAV)   :: albedo ! surface albedo

    INTEGER :: i, j, k

    REAL, PARAMETER :: spo2 = avno / (mwt_air_SI*grav) * 1.E-4 ! [part./cm^2 * 1/Pa]

    REAL :: zen, dz !, aa1
    REAL :: coszet, part_surf
    REAL :: dv2s1 ! dv2s (interval 0)
    REAL :: dv3s1 ! dv3s (interval 0)
    REAL :: v3s1  ! v3s (intervals 1-2)
    REAL :: v3s2  ! v3s (intervals 3,4,5,7)
    REAL :: dlv2_h, v3s_du_h
    REAL :: v3_du, sig0_o2k, sig0_o3k, dtau_0

    ! effective optical depths for the different intervals (tau_6 is not used)
    REAL :: tau_0(kproma_day), tau_1, tau_2, tau_3, tau_4, tau_5, tau_7

    !-------------------------------------------------------------------------
    ! indices for lookup table
    INTEGER, DIMENSION(kproma_day,gnk) :: i0, i1, i2, i3
    ! PBL layer index for simplified aerosol condition
    INTEGER, DIMENSION(kproma_day) :: pblind
    ! save values of the level above
    REAL, DIMENSION(kproma_day,0:gnk) :: sig_top_o2, sig_top_o3
    REAL, DIMENSION(kproma_day,0:gnk) :: temp,  & ! temperature           [K]
                                         press    ! pressure              [Pa]
    REAL, DIMENSION(kproma_day,gnk) :: v2s_m,    & ! slant O2 column [m/cm^2]
                                       v3_du1,   & ! slant O3 column [DU] (intervals 0-2)
                                       v3_du2,   & ! slant O3 column [DU] (intervals 3-7)
                                       tnorm_sr, & ! normalized T (T0 = 240 K)
                                       tnorm,    & ! normalized T (T0 = 250 K)
                                       dlv2,     & ! log(v2)
                                       dens,     & ! density [molec/cm3]
                                       r_m,      &
                                       r_o2
    REAL, DIMENSION(kproma_day,gnk,0:MAXWAV) :: fint ! integrated actinic flux
    REAL, DIMENSION(kproma_day,8) :: &
      fj_corr                     ! correction for zenith angles > sza_thr
    INTEGER :: klev        ! number of levels
  ! for GEMMACH debug
    logical local_dbg
    !-------------------------------------------------------------------------

    klev   = gnk

    local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))
    IF (local_dbg) write(chm_lun_out, *)'*** inside mach_messy_jvalues ***'

!! initialize output array
    jval_2d = 0.0_DP

    DO j = 1, kproma_day
      v3(j,:)     = v3_2d(iu0(j),:)
      press(j,1:) = press_2d(iu0(j),:)
      temp(j,1:)  = temp_2d(iu0(j),:)
      rhum(j,:)   = rhum_2d(iu0(j),:)
      aclc(j,:)   = aclc_2d(iu0(j),:)
      clp(j,:)    = clp_2d(iu0(j),:)
      albedo(j,:) = albedo_1d(iu0(j))
      slf(j)      = slf_1d(iu0(j))       ! land sea fraction
      pblind(j)   = pbllev(iu0(j))       ! land sea fraction

      ! Rodgers, The radiative heat budget of the troposphere and lower
      ! stratosphere, Res. Rep. n0. a2, planatary circualtions project (1967)
      coszet = MAX(cossza_1d(iu0(j)),cossza_thr)
      u0(j)  = SQRT(coszet**2*1224.+1.)*1./35.

      ! pressure and temperature at virtual level 0
      ! above level 0 an only absorbing astmosphere is assumed
      press(j,0) = 0.37 * press(j,1) ! 0.37 = 1/e
      temp(j,0)  = (temp(j,2)-temp(j,1))/(press(j,2)-press(j,1)) * &
                   (-0.63) * press(j,1) + temp(j,1)
      ! assume constant mixing ratio above toa
      ! fix ozone mixing ratio in the first nt_lev layers by haloe climatology
      ! assume constant mixing ratio above toa (level 1)

      ! aerosol particle column per layer [part./cm^2]
      ! Christoph's assumption about aerosol particle density:
      ! n(z) = par_surf * (press(z)/press(surface))^3
      ! part_surf: aerosol particle  at the surface [part./cm^3]
      part_surf = part_rural * slf(j) + part_sea * (1.-slf(j))
      DO k = 1,klev
        dz = dz2d(iu0(j),k)
        part_col(j,k) = part_surf * 0.5 * &
          ((press(j,k)/press(j,klev))**3 + (press(j,k-1)/press(j,klev))**3) &
          * dz * 100.
      ENDDO
    ENDDO

    !-------------------------------------------------------------------------
    ! calculate O2 and O3 column densities
    DO j = 1,kproma_day
      v3s(j,:) = v3(j,:) / u0(j)

    ! O2
      v2(j,:)  = spo2 * relo2 * press(j,:)
      ! slant columns
      v2s(j,:) = v2(j,:) / u0(j)
    ENDDO

    !-------------------------------------------------------------------------
    ! Practical improved flux method (Zdunkwoski) to calculate actinic fluxes
    CALL flux_cal(aclc, clp, temp, rhum, albedo, u0, slf, part_col, pblind, &
                  v2, v3, fact, facth, kproma_day, klev)
    IF (local_dbg) write(chm_lun_out, *)'mach_messy_jvalues done flux cal'


    ! ------------------------------------------------------------------------

    ! do first the calculation of the optical depth tau_0 above
    ! the model atmosphere

    ! initialization of tau_0
    k = 1
    DO j = 1,kproma_day

      dlv2_h   = MAX(44.5,MIN(56.,LOG(v2s(j,0))))
      v3s_du_h = MIN(300., v3s(j,0)* constant * 1.E+3)
      i0(j,k)       = MIN(58,INT(AINT((dlv2_h-44.5)/0.2) + 1.00001))
      tnorm_sr(j,k) = (temp(j,0)-240.) / 240.
      sig_top_O2(j,0) = p2(a0_1_1_O2(i0(j,k))*dlv2_h + a0_1_2_O2(i0(j,k)), &
        a0_2_1_O2(i0(j,k))*dlv2_h + a0_2_2_O2(i0(j,k)), &
        a0_3_1_O2(i0(j,k))*dlv2_h + a0_3_2_O2(i0(j,k)), &
        v3s_du_h) * &
        p2(1.,c0_1_1_O2(i0(j,k))*dlv2_h+c0_1_2_O2(i0(j,k)), &
        c0_2_1_O2(i0(j,k))*dlv2_h+c0_2_2_O2(i0(j,k)), &
        tnorm_sr(j,k))
      sig_top_o3(j,0) = p1(a0_1_1_o3(i0(j,k))*dlv2_h + a0_1_2_o3(i0(j,k)), &
        a0_2_1_o3(i0(j,k))*dlv2_h + a0_2_2_o3(i0(j,k)), &
        v3s_du_h) * &
        p2(1.,c0_1_1_o3(i0(j,k))*dlv2_h+c0_1_2_o3(i0(j,k)), &
        c0_2_1_o3(i0(j,k))*dlv2_h+c0_2_2_o3(i0(j,k)), &
        tnorm_sr(j,k))
      tau_0(j)     = p2(ct(1),ct(2),ct(3),(dlv2_h-47.)/47.)
    ENDDO

    ! basic species, minimal set, definition of arrays for tables and fits
    DO k = 1,klev ! altitude loop

      DO j  = 1,kproma_day ! longitude loop

        ! change of units [v2s_m] = meter, [v3_du] = DU
        v2s_m(j,k) = v2s(j,k) * constant * 1.E-2
        v3_du      = v3s(j,k) * constant * 1.E+3
        dlv2(j,k)  = LOG(v2s(j,k))

        ! allowed ranges for columns
        ! interval: 0
        IF (dlv2(j,k)>= 56.) THEN
          dlv2(j,k)  = 56.
          dv2s1 = 1.E+25
        ELSE
          dv2s1     = v2s(j,k)-v2s(j,k-1)
        ENDIF

        ! definition range of polynomial
        ! interval: 1
        IF (v2s_m(j,k)>= 500.) THEN
          v2s_m(j,k) = 500.
          v2s2(j,k)  = 1.3602E+24
        ELSE
          v2s2(j,k)   = v2s(j,k)
        ENDIF

        ! interval: 0 - 2
        IF (v3_du>= 300.) THEN
          dv3s1     = 1.E+19
          v3_du1(j,k) = 300.
          v3s1      = 8.161E+18
        ELSE
          dv3s1     = v3s(j,k)-v3s(j,k-1)
          v3_du1(j,k) = v3_du
          v3s1      = v3s(j,k)
        ENDIF

        ! interval: 3 - 7
        IF (v3_du>= 3000.) THEN
          v3_du2(j,k) = 3000.
          v3s2      = 8.161E+19
        ELSE
          v3_du2(j,k) = v3_du
          v3s2      = v3s(j,k)
        ENDIF

        ! scaling of temperature variable
        tnorm_sr(j,k) = (temp(j,k)-240.) / 240. ! for interval 0
        tnorm(j,k)    = (temp(j,k)-250.) / 250. ! for interval 1-7
        dens(j,k)  =  PRESS(J,K)/(kboltz*TEMP(J,K))*1.E-6

        ! indices for lookup table
        ! mz_rs_20140314+
        ! At very high altitudes, dlv2 can be smaller than 44.5 and the
        ! resulting index i0 will be out of the valid range. This bug
        ! fix sets i0=1 in that case.
        i0(j,k)  = MAX(1,MIN(58,INT(AINT((dlv2(j,k)-44.5)/0.2) + 1.00001)))
        !i0(j,k)  = MIN(58,INT(AINT((dlv2(j,k)-44.5)/0.2) + 1.00001))
        ! mz_rs_20140314-
        i1(j,k)  = MIN(INT((v3_du-0.5)/2.5) +1,115)
        i1(j,k)  = ifil(i1(j,k))
        i2(j,k)  = MIN(INT((v3_du-0.5)/2.5) +1,115)
        i2(j,k)  = ifil(i2(j,k))
        i3(j,k)  = MIN(INT((v3_du-5.)/25.) +1,115)
        i3(j,k)  = ifil(i3(j,k))

        ! calculation of optical depths and integrated actinic flux
        ! effective optical depths and integrated actinic fluxes
        sig0_O2k = &
          p2(a0_1_1_O2(i0(j,k))*dlv2(j,k) + a0_1_2_O2(i0(j,k)), &
          a0_2_1_O2(i0(j,k))*dlv2(j,k) + a0_2_2_O2(i0(j,k)), &
          a0_3_1_O2(i0(j,k))*dlv2(j,k) + a0_3_2_O2(i0(j,k)), &
          v3_du1(j,k)) * &
          p2(1.,c0_1_1_O2(i0(j,k))*dlv2(j,k)+c0_1_2_O2(i0(j,k)), &
          c0_2_1_O2(i0(j,k))*dlv2(j,k)+c0_2_2_O2(i0(j,k)), &
          tnorm_sr(j,k))

        sig0_o3k = &
          p1(a0_1_1_o3(i0(j,k))*dlv2(j,k) + a0_1_2_o3(i0(j,k)), &
          a0_2_1_o3(i0(j,k))*dlv2(j,k) + a0_2_2_o3(i0(j,k)), &
          v3_du1(j,k)) * &
          p2(1.,c0_1_1_o3(i0(j,k))*dlv2(j,k)+c0_1_2_o3(i0(j,k)), &
          c0_2_1_o3(i0(j,k))*dlv2(j,k)+c0_2_2_o3(i0(j,k)), &
          tnorm_sr(j,k))
        dtau_0     = &
          0.5*(sig_top_O2(j,k-1) + sig0_O2k) * dv2s1 + &
          0.5*(sig_top_o3(j,k-1) + sig0_o3k) * dv3s1
        tau_0(j)   = tau_0(j) + dtau_0

        sig_top_O2(j,k) = sig0_O2k
        sig_top_o3(j,k) = sig0_o3k

        ! calculate the effective optical depths tau_i
        ! intervals i = 1..3 from lookup tables
        tau_1 = &
          p1(taub1_1(i1(j,k)),taua1_1(i1(j,k)),v3_du1(j,k)) + &
          p1(taub1_2(i1(j,k)),taua1_2(i1(j,k)),v3_du1(j,k)) * &
          v2s_m(j,k) + &
          p1(taub1_3(i1(j,k)),taua1_3(i1(j,k)),v3_du1(j,k)) * &
          v2s_m(j,k)**2
        tau_2 = p1(taub2(i2(j,k)),taua2(i2(j,k)),v3_du1(j,k))
        tau_3 = p1(taub3(i3(j,k)),taua3(i3(j,k)),v3_du2(j,k))
        ! for intervals i = 4..7 from the polynomials with B_i_O3(j)
        tau_4 = p2(B_4_O3(1),B_4_O3(2),B_4_O3(3),v3_du2(j,k))
        tau_5 = p1(B_5_O3(1),B_5_O3(2),v3_du2(j,k))
        ! tau_6 is not used
        tau_7 = p1(B_7_O3(1),B_7_O3(2),v3_du2(j,k))
        fint(j,k,0) = EXP(-tau_0(j))* SR_toa_flux
        fint(j,k,1) = EXP(-tau_1 + crs_o3(1)*v3s1 + &
          crs_O2(1)*v2s2(j,k))*f0(1)*fact(j,k,1)
        fint(j,k,2) = EXP(-tau_2 + crs_o3(2)*v3s1)*f0(2)*fact(j,k,2)
        IF (v3_du2(j,k) <= 1000.) THEN
          fint(j,k,3) = EXP(-tau_3 + crs_o3(3)*v3s2)*f0(3)*fact(j,k,3)
        ELSE
          fint(j,k,3) = 0.
        ENDIF
        IF (v3_du2(j,k) <= 2500.) THEN
          fint(j,k,4) = EXP(-tau_4 + crs_o3(4)*v3s2)*f0(4)*fact(j,k,4)
        ELSE
          fint(j,k,4) = 0.
        ENDIF
        fint(j,k,5) = EXP(-tau_5 + crs_o3(5)*v3s2)*f0(5)*fact(j,k,5)
        fint(j,k,6) = f0(6)*fact(j,k,6)
        fint(j,k,7) = EXP(-tau_7 + crs_o3(7)*v3s2)*f0(7)*fact(j,k,7)
      ENDDO
    ENDDO

    IF (local_dbg) write(chm_lun_out, *)'mach_messy_jvaluesdone cal fint/tau'


    !-------------------------------------------------------------------------

    ! reduction factors r_O2 und r_m for lyman-alpha line
    ! S. Chabrillat and G. Kockarts grl 24, p. 2659-2662 [2633]
    ! with the corrected values from grl 25 p.79 [2634]

    ! zero for non-ma
    r_m(:,:)  = 0. ! for CO2, CH4, SF6, H2SO4, and H2O
    r_O2(:,:) = 0. ! only for O2
    ! lyman-alpha photolysis for ma (only relevant in the mesosphere)
    IF (lmidatm) THEN
      DO i = 1,3
        DO k = 1,klev
          DO j = 1,kproma_day
            r_m(j,k)  = r_m(j,k)  + b_la(i) * EXP(-c_la(i)*v2s2(j,k))
            r_O2(j,k) = r_O2(j,k) + d_la(i) * EXP(-e_la(i)*v2s2(j,k))
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    !-------------------------------------------------------------------------

    j_loop: DO j = 1,kproma_day
      ! correction for zenith angles > sza_thr
      ! ref2642 from Lamago et al. (ACP 3, 1981-1990, 2003)
      fj_corr(j,:) = 1.
      zen = ACOS(cossza_1d(iu0(j)))*180./pi
      IF (zen > sza_thr) THEN
        fj_corr(j,1) = EXP(    1.2*((82.-zen)/5.5+1.))
        fj_corr(j,2) = EXP(1.3*1.2*((82.-zen)/5.5+1.))
        fj_corr(j,3) = EXP(1.4*1.2*((82.-zen)/5.5+1.))
        fj_corr(j,4) = EXP(1.5*1.2*((82.-zen)/5.5+1.))
        fj_corr(j,5) = EXP(2.0*1.2*((82.-zen)/5.5+1.))
        fj_corr(j,6) = EXP(3.0*1.2*((82.-zen)/5.5+1.))
        fj_corr(j,7) = EXP(4.0*1.2*((82.-zen)/5.5+1.))
        fj_corr(j,8) = EXP(5.0*1.2*((82.-zen)/5.5+1.))
      ENDIF
    ENDDO j_loop

    !-------------------------------------------------------------------------

    ! UV heating rates (only far UV)

    IF (PRESENT(relo3_2d)) THEN
      rh_o3_2d      = 0.0_dp
      rh_o2_2d      = 0.0_dp

      DO k = 1,klev
        DO j = 1,kproma_day

          ! O2
          sig_h_o2(1) = p2(                                      &
            a0_1_1_h_O2(i0(j,k))*dlv2(j,k)+a0_1_2_h_O2(i0(j,k)), &
            a0_2_1_h_O2(i0(j,k))*dlv2(j,k)+a0_2_2_h_O2(i0(j,k)), &
            a0_3_1_h_O2(i0(j,k))*dlv2(j,k)+a0_3_2_h_O2(i0(j,k)), &
            v3_du1(j,k))                                         &
            * p2(1.,c0_1_1_h_O2(i0(j,k))*dlv2(j,k)               &
            + c0_1_2_h_O2(i0(j,k)),                              &
            c0_2_1_h_O2(i0(j,k))*dlv2(j,k)                       &
            + c0_2_2_h_O2(i0(j,k)),tnorm_sr(j,k) )
          sig_h_o2(2) = p1(                                      &
            b1_1_h_O2(i1(j,k)),a1_1_h_O2(i1(j,k)),v3_du1(j,k) )  &
            + p1( b1_2_h_O2(i1(j,k)),                            &
            a1_2_h_O2(i1(j,k)),v3_du1(j,k) )                     &
            * v2s_m(j,k)
          dh_O2(j,k) = (sig_h_o2(1) * fint(j,k,0) +              &
            sig_h_o2(2) * fint(j,k,1) ) * relo2
          rh_O2_2d(iu0(j),k) = REAL(MAX(0.0, dh_O2(j,k) *fj_corr(j,8)), dp) ! 5
          ! O3
          sig_h_o3(1) = p1(a0_1_1_h_o3(i0(j,k))*dlv2(j,k)        &
            + a0_1_2_h_o3(i0(j,k)),                              &
            a0_2_1_h_o3(i0(j,k))*dlv2(j,k)                       &
            + a0_2_2_h_o3(i0(j,k)),                              &
            v3_du1(j,k)) *                                       &
            p2(1.,c0_1_1_h_o3(i0(j,k))*dlv2(j,k)+                &
            c0_1_2_h_o3(i0(j,k)),                                &
            c0_2_1_h_o3(i0(j,k))*dlv2(j,k)+                      &
            c0_2_2_h_o3(i0(j,k)),tnorm_sr(j,k))
          sig_h_o3(2) = p1(b1_1_h_o3(i1(j,k))                    &
            ,a1_1_h_o3(i1(j,k)),v3_du1(j,k))+                    &
            p1(b1_2_h_o3(i1(j,k)),a1_2_h_o3(i1(j,k))             &
            ,v3_du1(j,k))*                                       &
            v2s_m(j,k)
          dh_o3(j,k)   =   (sig_h_o3(1)    * fint(j,k,0) +       &
            sig_h_o3(2)    * fint(j,k,1) ) * relo3_2d(iu0(j),k)*1.e6
          rh_o3_2d(iu0(j),k) = REAL(MAX(0.0, dh_o3(j,k) *fj_corr(j,6)), dp) ! 3
        ENDDO
      ENDDO

      IF (PRESENT(fhuv_2d) .AND. PRESENT(fhuvdna_2d)) then
         CALL jval_cal_uv(fint, fact, facth, v3_du2, fhuv_2d, fhuvdna_2d, &
                          iu0, kproma_day, gni, klev)
         IF (local_dbg) write(chm_lun_out, *)'mach_messy_jvalues done cal_uv'
      ENDIF
    ENDIF

    !-------------------------------------------------------------------------

    CALL jval_cal
    IF (local_dbg) write(chm_lun_out, *)'mach_messy_jvalues done jval_cal'

    !-------------------------------------------------------------------------

  CONTAINS

  ! **************************************************************************
  ! include the subroutines jval_cal and jval_cal_* which are
  ! dynamically generated by jvpp:
  include "mach_messy_jval_inc.cdk"

    !-------------------------------------------------------------------------
!! extend jvalues SR to here to avoid dynamic allocation of array for OMP>1
  END SUBROUTINE mach_messy_jvalues

! ****************************************************************************

  SUBROUTINE flux_cal(aclc, clp, temp, rhum, albedo, u0, slf, part_col, pblind, &
                      v2, v3, fact, facth, kproma_day, klev)
    USE mach_pkg_messy_mod,  only: MAXWAV, crray

    IMPLICIT NONE
    INTEGER ,                                 INTENT(IN)  :: kproma_day, klev
    REAL, DIMENSION(kproma_day, klev),        INTENT(IN)  :: aclc     ! cloud fraction                 [1]
    REAL, DIMENSION(kproma_day, klev),        INTENT(IN)  :: clp      ! cloud liquid water path per layer [g/m^2]
    REAL, DIMENSION(kproma_day, 0:klev),      INTENT(IN)  :: temp     ! Air temperatures               [K]
    REAL, DIMENSION(kproma_day, klev),        INTENT(IN)  :: rhum     ! relative humudity              [%]
    REAL, DIMENSION(kproma_day, MAXWAV),      INTENT(IN)  :: albedo   ! surface albedo
    REAL, DIMENSION(kproma_day),              INTENT(IN)  :: u0       ! cosine of solar zenith angle
    REAL, DIMENSION(kproma_day),              INTENT(IN)  :: slf      ! fraction of sea surface
    REAL, DIMENSION(kproma_day, klev),        INTENT(IN)  :: part_col ! aerosol part. column per layer [part/cm^2]
    REAL, DIMENSION(kproma_day, 0:klev),      INTENT(IN)  :: v2       ! O2 column density [part./cm^2]
    REAL, DIMENSION(kproma_day, 0:klev),      INTENT(IN)  :: v3       ! O3 column density [part./cm^2]
    INTEGER, DIMENSION(kproma_day),           INTENT(IN)  :: pblind   ! PBL layer index
    REAL, DIMENSION(kproma_day, klev,MAXWAV), INTENT(OUT) :: fact     ! actinic flux [photons/(cm^2 s)]
    REAL, DIMENSION(kproma_day, klev, 3:5),   INTENT(OUT) :: facth    ! flux [photons/(cm^2 s)]

    REAL, DIMENSION(kproma_day, klev,MAXWAV) :: taer_sca ! scattering optical depth of aerosol (Shettle, Fenn)
    REAL, DIMENSION(kproma_day, klev,MAXWAV) :: taer_abs ! absorption optical depth of aerosol (Shettle, Fenn)
    REAL, DIMENSION(kproma_day, klev,MAXWAV) :: gaer     ! asymmetry factor of aerosol (Shettle, Fenn)

    REAL, DIMENSION(kproma_day, klev,MAXWAV) :: taus_clr ! scattering optical depth (clear sky)
    REAL, DIMENSION(kproma_day, klev,MAXWAV) :: taua_clr ! absorption optical depth (clear sky)
    REAL, DIMENSION(kproma_day, klev,MAXWAV) :: g_clr    ! asymmetry factor of layer (clear sky)

    REAL, DIMENSION(kproma_day, klev) :: taus_cld ! scattering optical depth (slingo cloud parameter)
    REAL, DIMENSION(kproma_day, klev) :: taua_cld ! absorption optical depth (slingo cloud parameter)
    REAL, DIMENSION(kproma_day, klev) :: g_cld    ! asymmetry factor of layer(slingo cloud parameter)

    REAL, DIMENSION(kproma_day, MAXWAV) :: ftop  ! absorption correction above model

    REAL    :: ta_O2, ta_o3, temp_avg, dv2
    INTEGER :: j, k, l

      ! reference temperatures
    REAL, PARAMETER :: t1 = 226.
    REAL, PARAMETER :: t2 = 263.

    REAL, PARAMETER :: cr_O2(MAXWAV) = & ! O2 cross section
        (/ 7.6300E-24, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
        0.0000E+00, 0.0000E+00, 0.0000E+00/)
    REAL, PARAMETER :: cr_o3(MAXWAV) = & ! O3 cross section (T = 298 K)
        (/ 3.6500E-19, 1.7900E-18, 2.7700E-19, 1.0500E-19, &
        2.6000E-20, 0.0000E+00, 4.5500E-21/)
    REAL, PARAMETER :: x1_o3(MAXWAV) = & ! T-dependence of a0_o3
        (/-2.7027E-23, 5.4054E-22, 0.0000E+00, 0.0000E+00, &
        0.0000E+00, 0.0000E+00, 0.0000E+00/)
    REAL, PARAMETER :: x2_o3(MAXWAV) = & ! T-dependence of a0_o3
        (/-4.1828E-25, 8.3655E-24, 0.0000E+00, 0.0000E+00, &
        0.0000E+00, 0.0000E+00, 0.0000E+00/)

    ! calculate aerosol parameter
    CALL aero_2d(slf, rhum, part_col, pblind, taer_abs, taer_sca, gaer, &
                 kproma_day, klev)

      ! cloud parameter (slingo parameterization)
    CALL slingo(clp, taus_cld, taua_cld, g_cld, kproma_day, klev)

      ! correction for absorption above model toa
    DO l = 1,MAXWAV
      j_loop1: DO j = 1,kproma_day
        ta_O2 = cr_O2(l) * v2(j,0)
        ta_o3 = (((temp(j,0)-t2)*x2_o3(l) +x1_o3(l))*(temp(j,0)-t1) &
          +cr_o3(l)) * v3(j,0)
        ftop(j,l) = EXP(-1./u0(j)*(ta_O2 + ta_o3))
      ENDDO j_loop1
    ENDDO

    ! optical depths of each layer
    DO l = 1,MAXWAV
      taua_clr(:, 1, l) = 0.0
      taus_clr(:, 1, l) = 0.
      g_clr(:, 1, l)    = 0.

      k_loop1: DO k  = 2, klev
        j_loop2: DO j = 1, kproma_day

          temp_avg  = 0.5 * (temp(j,k-1) + temp(j,k))

          dv2 = v2(j,k) - v2(j,k-1)
          ta_O2     = cr_O2(l) * dv2
          ta_o3     = (((temp_avg-t2)*x2_o3(l) +x1_o3(l)) &
                      *(temp_avg-t1)+cr_o3(l)) * (v3(j,k) - v3(j,k-1))
          taua_clr(j,k,l) = ta_o3 + ta_O2 + taer_abs(j,k,l)

          ! for calculation of g_clr see eg.: Radiation and cloud processes
          ! in the atmosphere, K.N. Liou, p.155 eq. 3.8.7
          taus_clr(j,k,l) = crray(l)*1./0.21*dv2 + taer_sca(j,k,l)
          g_clr(j,k,l)    = gaer(j,k,l)*taer_sca(j,k,l) / taus_clr(j,k,l)
        ENDDO j_loop2
      ENDDO k_loop1
    ENDDO

    ! practical improved flux method (pifm) to calculate actinic fluxes
    CALL pifm(taus_cld, taua_cld, g_cld, taus_clr, taua_clr, g_clr, &
              aclc, u0, albedo, fact, facth, kproma_day, klev)

    DO k = 1,klev
      DO j = 1,kproma_day
        fact(j,k,:) = fact(j,k,:) * ftop(j,:)
        facth(j,k,3:5) = facth(j,k,3:5) * ftop(j,3:5)
      ENDDO
    ENDDO

  END SUBROUTINE flux_cal

  !-------------------------------------------------------------------------
  SUBROUTINE aero_2d(slf, rhum, part_col, pblind, taer_abs, taer_sca, gaer, &
                     kproma_day, klev)
    USE mach_pkg_messy_mod,  only: MAXWAV

    IMPLICIT NONE
    INTEGER ,                                 INTENT(IN)  :: kproma_day, klev
    REAL, DIMENSION(kproma_day),              INTENT(IN)  :: slf      ! fraction of sea surface
    REAL, DIMENSION(kproma_day, klev),        INTENT(IN)  :: rhum     ! relative humudity              [%]
    REAL, DIMENSION(kproma_day, klev),        INTENT(IN)  :: part_col ! aerosol part. column per layer [part/cm^2]
    INTEGER, DIMENSION(kproma_day),           INTENT(IN)  :: pblind   ! PBL layer index
    REAL, DIMENSION(kproma_day, klev,MAXWAV), INTENT(OUT) :: taer_sca ! scattering optical depth of aerosol (Shettle, Fenn)
    REAL, DIMENSION(kproma_day, klev,MAXWAV), INTENT(OUT) :: taer_abs ! absorption optical depth of aerosol (Shettle, Fenn)
    REAL, DIMENSION(kproma_day, klev,MAXWAV), INTENT(OUT) :: gaer     ! asymmetry factor of aerosol (Shettle, Fenn)

    REAL :: domin, a1_sc, a2_sc, a4_sc
    REAL :: b1_sc, b2_sc, b4_sc, a1_ab, a2_ab, a4_ab, b1_ab, b2_ab, b4_ab
    REAL :: a1_g, a2_g, a4_g, b1_g, b2_g, b4_g, bsa, baa, ga
    REAL :: &
        aext(MAXWAV,8,4),  & ! extinction coefficient [1/km]
        asca(MAXWAV,8,4),  & ! effective scattering cross section [cm^2/part.]
        aabs(MAXWAV,8,4),  & ! effective absorption cross section [cm^2/part.]
        ag(MAXWAV,8,4)     ! effective asymmetry factor
      ! aerosol data from: Models for aerosols of the lower atmosphere
      ! and the effects of humidity variations on their optical properties
      ! E. P. Shettle and R. W. Fenn (1979)
      ! Environmental Research Paper No. 676

      ! relative humidity of lookup table [%]
      REAL, PARAMETER :: rh_ref(8) = (/ 0., 50., 70., 80., 90., 95., 98., 99./)

      INTEGER :: i_ref ! interval index for relative humidity

      ! different aerosol types:
      ! 1 = rural aerosol
      ! 2 = maritime aerosol
      ! 3 = urban aerosol
      ! 4 = free troposphere aerosol

      INTEGER :: i, j, k, l
      REAL, PARAMETER :: pn_ref(4) = (/ 15000., 4000., 20000., 5000./)
      ! rural model
      ! relative humidity = 0.0 %
      aext(:,1,1) = (/ 3.0331E-01, 2.6173E-01, 2.5400E-01, 2.5011E-01, 2.4395E-01, 2.1722E-01, 1.3731E-01 /)
      aabs(:,1,1) = (/ 9.2798E-02, 2.2387E-02, 1.6832E-02, 1.5121E-02, 1.3489E-02, 1.1261E-02, 8.4159E-03 /)
      ag(:,1,1)   = (/ 7.5290E-01, 6.8383E-01, 6.7781E-01, 6.7575E-01, 6.7346E-01, 6.6742E-01, 6.4481E-01 /)
      ! relative humidity = 50.0 %
      aext(:,2,1) = (/ 3.1458E-01, 2.7099E-01, 2.6296E-01, 2.5893E-01, 2.5255E-01, 2.2496E-01, 1.4230E-01 /)
      aabs(:,2,1) = (/ 9.4048E-02, 2.2410E-02, 1.6798E-02, 1.5083E-02, 1.3467E-02, 1.1366E-02, 8.4195E-03 /)
      ag(:,2,1)   = (/ 7.5485E-01, 6.8880E-01, 6.8292E-01, 6.8088E-01, 6.7856E-01, 6.7231E-01, 6.5029E-01 /)
      ! relative humidity = 70.0 %
      aext(:,3,1) = (/ 3.3914E-01, 2.9169E-01, 2.8296E-01, 2.7858E-01, 2.7167E-01, 2.4184E-01, 1.5346E-01 /)
      aabs(:,3,1) = (/ 9.7106E-02, 2.2834E-02, 1.7015E-02, 1.5237E-02, 1.3562E-02, 1.1400E-02, 8.5197E-03 /)
      ag(:,3,1)   = (/ 7.5859E-01, 6.9748E-01, 6.9220E-01, 6.9041E-01, 6.8844E-01, 6.8315E-01, 6.6063E-01 /)
      ! relative humidity = 80.0 %
      aext(:,4,1) = (/ 4.6086E-01, 3.9492E-01, 3.8325E-01, 3.7747E-01, 3.6840E-01, 3.2934E-01, 2.1149E-01 /)
      aabs(:,4,1) = (/ 1.0993E-01, 2.4472E-02, 1.7877E-02, 1.5895E-02, 1.4081E-02, 1.2031E-02, 8.8703E-03 /)
      ag(:,4,1)   = (/ 7.6932E-01, 7.2720E-01, 7.2359E-01, 7.2237E-01, 7.2103E-01, 7.1720E-01, 6.9661E-01 /)
      ! relative humidity = 90.0 %
      aext(:,5,1) = (/ 6.7903E-01, 6.1749E-01, 5.5327E-01, 5.1363E-01, 4.4529E-01, 2.0029E-01, 2.6598E-01 /)
      aabs(:,5,1) = (/ 1.2497E-01, 2.5991E-02, 1.8534E-02, 1.6356E-02, 1.4461E-02, 1.2863E-02, 9.1785E-03 /)
      ag(:,5,1)   = (/ 7.7428E-01, 7.5237E-01, 7.5024E-01, 7.4945E-01, 7.4847E-01, 7.4511E-01, 7.2836E-01 /)
      ! relative humidity = 95.0 %
      aext(:,6,1) = (/ 8.0780E-01, 7.0467E-01, 6.8690E-01, 6.7817E-01, 6.6449E-01, 6.0441E-01, 4.0403E-01 /)
      aabs(:,6,1) = (/ 1.3235E-01, 2.6829E-02, 1.8930E-02, 1.6640E-02, 1.4677E-02, 1.3195E-02, 9.3323E-03 /)
      ag(:,6,1)   = (/ 7.7486E-01, 7.6184E-01, 7.6042E-01, 7.5986E-01, 7.5910E-01, 7.5619E-01, 7.4127E-01 /)
      ! relative humidity = 98.0 %
      aext(:,7,1) = (/ 1.0511E+00, 9.3767E-01, 9.1720E-01, 9.0699E-01, 8.9083E-01, 8.1871E-01, 5.6955E-01 /)
      aabs(:,7,1) = (/ 1.4878E-01, 2.9403E-02, 2.0470E-02, 1.7882E-02, 1.5667E-02, 1.4020E-02, 9.8336E-03 /)
      ag(:,7,1)   = (/ 7.7704E-01, 7.7311E-01, 7.7252E-01, 7.7225E-01, 7.7183E-01, 7.6996E-01, 7.5918E-01 /)
      ! relative humidity = 99.0 %
      aext(:,8,1) = (/ 1.2932E+00, 1.1771E+00, 1.1548E+00, 1.1435E+00, 1.1253E+00, 1.0428E+00, 7.4704E-01 /)
      aabs(:,8,1) = (/ 1.6170E-01, 3.0955E-02, 2.1241E-02, 1.8451E-02, 1.6103E-02, 1.4596E-02, 1.0031E-02 /)
      ag(:,8,1)   = (/ 7.7794E-01, 7.7937E-01, 7.7928E-01, 7.7918E-01, 7.7898E-01, 7.7782E-01, 7.7005E-01 /)
      ! maritime model
      ! relative humidity = 0.0 %
      aext(:,1,2) = (/ 1.1714E-01, 1.0744E-01, 1.0552E-01, 1.0455E-01, 1.0299E-01, 9.6330E-02, 7.7600E-02 /)
      aabs(:,1,2) = (/ 2.3293E-02, 4.7427E-03, 3.3183E-03, 2.8928E-03, 2.5072E-03, 2.0902E-03, 1.3660E-03 /)
      ag(:,1,2)   = (/ 7.4802E-01, 7.0052E-01, 6.9643E-01, 6.9505E-01, 6.9353E-01, 6.8918E-01, 6.7501E-01 /)
      ! relative humidity = 50.0 %
      aext(:,2,2) = (/ 1.2518E-01, 1.1512E-01, 1.1335E-01, 1.1248E-01, 1.1111E-01, 1.0502E-01, 8.4545E-02 /)
      aabs(:,2,2) = (/ 2.3540E-02, 4.7178E-03, 3.2864E-03, 2.8636E-03, 2.4880E-03, 2.1214E-03, 1.3557E-03 /)
      ag(:,2,2)   = (/ 7.5592E-01, 7.1030E-01, 7.0599E-01, 7.0441E-01, 7.0253E-01, 6.9753E-01, 6.9120E-01 /)
      ! relative humidity = 70.0 %
      aext(:,3,2) = (/ 1.4988E-01, 1.3954E-01, 1.3763E-01, 1.3667E-01, 1.3514E-01, 1.2846E-01, 1.0760E-01 /)
      aabs(:,3,2) = (/ 2.4193E-02, 4.7563E-03, 3.2808E-03, 2.8459E-03, 2.4610E-03, 2.0977E-03, 1.3526E-03 /)
      ag(:,3,2)   = (/ 7.6810E-01, 7.3178E-01, 7.2840E-01, 7.2718E-01, 7.2576E-01, 7.2243E-01, 7.2140E-01 /)
      ! relative humidity = 80.0 %
      aext(:,4,2) = (/ 2.6744E-01, 2.5365E-01, 2.5127E-01, 2.5009E-01, 2.4827E-01, 2.4046E-01, 2.1678E-01 /)
      aabs(:,4,2) = (/ 2.6889E-02, 4.8768E-03, 3.2419E-03, 2.7725E-03, 2.3777E-03, 2.1241E-03, 1.3276E-03 /)
      ag(:,4,2)   = (/ 7.9436E-01, 7.7966E-01, 7.7799E-01, 7.7730E-01, 7.7636E-01, 7.7352E-01, 7.7178E-01 /)
      ! relative humidity = 90.0 %
      aext(:,5,2) = (/ 3.8330E-01, 3.6488E-01, 3.6140E-01, 3.5965E-01, 3.5686E-01, 3.4492E-01, 3.1084E-01 /)
      aabs(:,5,2) = (/ 2.9979E-02, 4.9768E-03, 3.1745E-03, 2.6765E-03, 2.2897E-03, 2.2327E-03, 1.3043E-03 /)
      ag(:,5,2)   = (/ 8.0012E-01, 7.9255E-01, 7.9193E-01, 7.9173E-01, 7.9151E-01, 7.9080E-01, 7.8588E-01 /)
      ! relative humidity = 95.0 %
      aext(:,6,2) = (/ 5.1615E-01, 4.9529E-01, 4.9145E-01, 4.8953E-01, 4.8648E-01, 4.7318E-01, 4.3302E-01 /)
      aabs(:,6,2) = (/ 3.1363E-02, 5.0292E-03, 3.1472E-03, 2.6332E-03, 2.2444E-03, 2.2567E-03, 1.3038E-03 /)
      ag(:,6,2)   = (/ 8.0759E-01, 8.0372E-01, 8.0348E-01, 8.0342E-01, 8.0340E-01, 8.0318E-01, 7.9781E-01 /)
      ! relative humidity = 98.0 %
      aext(:,7,2) = (/ 7.8803E-01, 7.6471E-01, 7.6086E-01, 7.5899E-01, 7.5608E-01, 7.4288E-01, 6.9431E-01 /)
      aabs(:,7,2) = (/ 3.2908E-02, 5.1181E-03, 3.1457E-03, 2.6121E-03, 2.2174E-03, 2.2908E-03, 1.3037E-03 /)
      ag(:,7,2)   = (/ 8.1493E-01, 8.1833E-01, 8.1815E-01, 8.1795E-01, 8.1751E-01, 8.1522E-01, 8.0982E-01 /)
      ! relative humidity = 99.0 %
      aext(:,8,2) = (/ 1.1277E+00, 1.1062E+00, 1.1025E+00, 1.1006E+00, 1.0976E+00, 1.0838E+00, 1.0284E+00 /)
      aabs(:,8,2) = (/ 3.4070E-02, 5.1842E-03, 3.1481E-03, 2.6024E-03, 2.2079E-03, 2.3432E-03, 1.3072E-03 /)
      ag(:,8,2)   = (/ 8.2128E-01, 8.2690E-01, 8.2699E-01, 8.2690E-01, 8.2662E-01, 8.2476E-01, 8.1909E-01 /)
      ! urban model
      ! relative humidity = 0.0 %
      aext(:,1,3) = (/ 3.3148E-01, 2.9477E-01, 2.8744E-01, 2.8369E-01, 2.7766E-01, 2.5096E-01, 1.6723E-01 /)
      aabs(:,1,3) = (/ 1.3662E-01, 1.0755E-01, 1.0369E-01, 1.0197E-01, 9.9459E-02, 8.9558E-02, 6.0876E-02 /)
      ag(:,1,3)   = (/ 7.7487E-01, 7.2345E-01, 7.1743E-01, 7.1491E-01, 7.1141E-01, 6.9828E-01, 6.5684E-01 /)
      ! relative humidity = 50.0 %
      aext(:,2,3) = (/ 3.4926E-01, 3.0928E-01, 3.0148E-01, 2.9751E-01, 2.9115E-01, 2.6305E-01, 1.7487E-01 /)
      aabs(:,2,3) = (/ 1.4041E-01, 1.0979E-01, 1.0577E-01, 1.0398E-01, 1.0139E-01, 9.1213E-02, 6.1876E-02 /)
      ag(:,2,3)   = (/ 7.7764E-01, 7.2861E-01, 7.2286E-01, 7.2045E-01, 7.1711E-01, 7.0448E-01, 6.6367E-01 /)
      ! relative humidity = 70.0 %
      aext(:,3,3) = (/ 4.6116E-01, 4.0134E-01, 3.9057E-01, 3.8521E-01, 3.7675E-01, 3.3989E-01, 2.2481E-01 /)
      aabs(:,3,3) = (/ 1.6134E-01, 1.2160E-01, 1.1665E-01, 1.1452E-01, 1.1148E-01, 9.9897E-02, 6.7307E-02 /)
      ag(:,3,3)   = (/ 7.8806E-01, 7.5152E-01, 7.4701E-01, 7.4507E-01, 7.4232E-01, 7.3156E-01, 6.9530E-01 /)
      ! relative humidity = 80.0 %
      aext(:,4,3) = (/ 7.0148E-01, 6.0491E-01, 5.8850E-01, 5.8046E-01, 5.6793E-01, 5.1364E-01, 3.4072E-01 /)
      aabs(:,4,3) = (/ 1.9374E-01, 1.3916E-01, 1.3278E-01, 1.3012E-01, 1.2643E-01, 1.1285E-01, 7.5244E-02 /)
      ag(:,4,3)   = (/ 7.9357E-01, 7.7371E-01, 7.7093E-01, 7.6965E-01, 7.6776E-01, 7.5984E-01, 7.3021E-01 /)
      ! relative humidity = 90.0 %
      aext(:,5,3) = (/ 1.0358E+00, 9.0438E-01, 8.8195E-01, 8.7096E-01, 8.5379E-01, 7.7850E-01, 5.2733E-01 /)
      aabs(:,5,3) = (/ 2.2461E-01, 1.5611E-01, 1.4843E-01, 1.4529E-01, 1.4105E-01, 1.2581E-01, 8.3559E-02 /)
      ag(:,5,3)   = (/ 7.9206E-01, 7.8587E-01, 7.8450E-01, 7.8377E-01, 7.8259E-01, 7.7704E-01, 7.5405E-01 /)
      ! relative humidity = 95.0 %
      aext(:,6,3) = (/ 1.4636E+00, 1.3097E+00, 1.2821E+00, 1.2683E+00, 1.2466E+00, 1.1491E+00, 8.0463E-01 /)
      aabs(:,6,3) = (/ 2.5345E-01, 1.7323E-01, 1.6439E-01, 1.6081E-01, 1.5601E-01, 1.3902E-01, 9.2399E-02 /)
      ag(:,6,3)   = (/ 7.8870E-01, 7.9251E-01, 7.9212E-01, 7.9176E-01, 7.9103E-01, 7.8710E-01, 7.7032E-01 /)
      ! relative humidity = 98.0 %
      aext(:,7,3) = (/ 2.2504E+00, 2.0973E+00, 2.0654E+00, 2.0488E+00, 2.0219E+00, 1.8952E+00, 1.4029E+00 /)
      aabs(:,7,3) = (/ 2.9543E-01, 1.9976E-01, 1.8929E-01, 1.8507E-01, 1.7944E-01, 1.5968E-01, 1.0601E-01 /)
      ag(:,7,3)   = (/ 7.8442E-01, 7.9688E-01, 7.9744E-01, 7.9746E-01, 7.9725E-01, 7.9526E-01, 7.8555E-01 /)
      ! relative humidity = 99.0 %
      aext(:,8,3) = (/ 2.9810E+00, 2.8574E+00, 2.8252E+00, 2.8076E+00, 2.7782E+00, 2.6340E+00, 2.0363E+00 /)
      aabs(:,8,3) = (/ 3.3013E-01, 2.2188E-01, 2.1001E-01, 2.0522E-01, 1.9882E-01, 1.7643E-01, 1.1639E-01 /)
      ag(:,8,3)   = (/ 7.8264E-01, 7.9839E-01, 7.9940E-01, 7.9963E-01, 7.9972E-01, 7.9897E-01, 7.9345E-01 /)
      ! free troposphere model
      ! relative humidity = 0.0 %
      aext(:,1,4) = (/ 9.6917E-02, 8.2981E-02, 8.0379E-02, 7.9069E-02, 7.6992E-02, 6.7999E-02, 4.1242E-02 /)
      aabs(:,1,4) = (/ 2.9018E-02, 5.9614E-03, 4.1862E-03, 3.6543E-03, 3.1696E-03, 2.6301E-03, 1.7271E-03 /)
      ag(:,1,4)   = (/ 7.4651E-01, 6.7635E-01, 6.7031E-01, 6.6827E-01, 6.6602E-01, 6.5979E-01, 6.3029E-01 /)
      ! relative humidity = 50.0 %
      aext(:,2,4) = (/ 1.0052E-01, 8.5881E-02, 8.3187E-02, 8.1836E-02, 7.9699E-02, 7.0451E-02, 4.2745E-02 /)
      aabs(:,2,4) = (/ 2.9357E-02, 5.9318E-03, 4.1461E-03, 3.6171E-03, 3.1448E-03, 2.6693E-03, 1.7141E-03 /)
      ag(:,2,4)   = (/ 7.4861E-01, 6.8186E-01, 6.7579E-01, 6.7365E-01, 6.7116E-01, 6.6398E-01, 6.3632E-01 /)
      ! relative humidity = 70.0 %
      aext(:,3,4) = (/ 1.0825E-01, 9.2398E-02, 8.9474E-02, 8.8007E-02, 8.5687E-02, 7.5677E-02, 4.6051E-02 /)
      aabs(:,3,4) = (/ 3.0195E-02, 5.9803E-03, 4.1383E-03, 3.5940E-03, 3.1103E-03, 2.6408E-03, 1.7100E-03 /)
      ag(:,3,4)   = (/ 7.5245E-01, 6.9062E-01, 6.8517E-01, 6.8330E-01, 6.8118E-01, 6.7508E-01, 6.4685E-01 /)
      ! relative humidity = 80.0 %
      aext(:,4,4) = (/ 1.4685E-01, 1.2477E-01, 1.2085E-01, 1.1890E-01, 1.1585E-01, 1.0271E-01, 6.3278E-02 /)
      aabs(:,4,4) = (/ 3.3598E-02, 6.1340E-03, 4.0904E-03, 3.5023E-03, 3.0054E-03, 2.6727E-03, 1.6797E-03 /)
      ag(:,4,4)   = (/ 7.6350E-01, 7.2093E-01, 7.1716E-01, 7.1587E-01, 7.1438E-01, 7.0974E-01, 6.8392E-01 /)
      ! relative humidity = 90.0 %
      aext(:,5,4) = (/ 2.1349E-01, 1.8287E-01, 1.7767E-01, 1.7512E-01, 1.7115E-01, 1.5382E-01, 9.7287E-02 /)
      aabs(:,5,4) = (/ 3.7513E-02, 6.2539E-03, 3.9994E-03, 3.3760E-03, 2.8912E-03, 2.8144E-03, 1.6493E-03 /)
      ag(:,5,4)   = (/ 7.6918E-01, 7.4714E-01, 7.4493E-01, 7.4409E-01, 7.4303E-01, 7.3917E-01, 7.1877E-01 /)
      ! relative humidity = 95.0 %
      aext(:,6,4) = (/ 2.5805E-01, 2.2351E-01, 2.1757E-01, 2.1465E-01, 2.1008E-01, 1.9000E-01, 1.2284E-01 /)
      aabs(:,6,4) = (/ 3.9242E-02, 6.3251E-03, 3.9701E-03, 3.3259E-03, 2.8370E-03, 2.8409E-03, 1.6490E-03 /)
      ag(:,6,4)   = (/ 7.6997E-01, 7.5695E-01, 7.5552E-01, 7.5494E-01, 7.5416E-01, 7.5099E-01, 7.3275E-01 /)
      ! relative humidity = 98.0 %
      aext(:,7,4) = (/ 3.2929E-01, 2.9127E-01, 2.8443E-01, 2.8102E-01, 2.7563E-01, 2.5151E-01, 1.6782E-01 /)
      aabs(:,7,4) = (/ 4.1205E-02, 6.4373E-03, 3.9673E-03, 3.2981E-03, 2.8015E-03, 2.8823E-03, 1.6499E-03 /)
      ag(:,7,4)   = (/ 7.7007E-01, 7.6643E-01, 7.6581E-01, 7.6552E-01, 7.6504E-01, 7.6274E-01, 7.4798E-01 /)
      ! relative humidity = 99.0 %
      aext(:,8,4) = (/ 4.0049E-01, 3.6140E-01, 3.5393E-01, 3.5014E-01, 3.4407E-01, 3.1647E-01, 2.1729E-01 /)
      aabs(:,8,4) = (/ 4.2656E-02, 6.5175E-03, 3.9683E-03, 3.2845E-03, 2.7889E-03, 2.9497E-03, 1.6547E-03 /)
      ag(:,8,4)   = (/ 7.6981E-01, 7.7196E-01, 7.7188E-01, 7.7176E-01, 7.7151E-01, 7.6990E-01, 7.5824E-01 /)

      DO l = 1,MAXWAV
        DO i = 1,4
          asca(l,:,i) = aext(l,:,i) - aabs(l,:,i)   ! [1/km]
          ! scaling for particle density 1 part./cm^3 = effective cross sections
          asca(l,:,i) = asca(l,:,i)/pn_ref(i)*1.E-5
          aabs(l,:,i) = aabs(l,:,i)/pn_ref(i)*1.E-5
        ENDDO

        DO k = 1,klev
          DO j = 1,kproma_day

            ! DO i = 1,7
            ! IF (rh_ref(i)<= rhum(j,k)) i_ref = i
            ! ENDDO
            IF (rh_ref(1)<= rhum(j,k)) i_ref = 1
            IF (rh_ref(2)<= rhum(j,k)) i_ref = 2
            IF (rh_ref(3)<= rhum(j,k)) i_ref = 3
            IF (rh_ref(4)<= rhum(j,k)) i_ref = 4
            IF (rh_ref(5)<= rhum(j,k)) i_ref = 5
            IF (rh_ref(6)<= rhum(j,k)) i_ref = 6
            IF (rh_ref(7)<= rhum(j,k)) i_ref = 7

            domin = 1./(rh_ref(i_ref)-rh_ref(i_ref+1))

            a1_sc = (asca(l,i_ref,1)-asca(l,i_ref+1,1)) * domin
            a2_sc = (asca(l,i_ref,2)-asca(l,i_ref+1,2)) * domin
            a4_sc = (asca(l,i_ref,4)-asca(l,i_ref+1,4)) * domin

            b1_sc = asca(l,i_ref,1) - a1_sc*rh_ref(i_ref)
            b2_sc = asca(l,i_ref,2) - a2_sc*rh_ref(i_ref)
            b4_sc = asca(l,i_ref,4) - a4_sc*rh_ref(i_ref)

            a1_ab = (aabs(l,i_ref,1)-aabs(l,i_ref+1,1)) * domin
            a2_ab = (aabs(l,i_ref,2)-aabs(l,i_ref+1,2)) * domin
            a4_ab = (aabs(l,i_ref,4)-aabs(l,i_ref+1,4)) * domin

            b1_ab = aabs(l,i_ref,1) - a1_ab*rh_ref(i_ref)
            b2_ab = aabs(l,i_ref,2) - a2_ab*rh_ref(i_ref)
            b4_ab = aabs(l,i_ref,4) - a4_ab*rh_ref(i_ref)

            a1_g = (ag(l,i_ref,1)-ag(l,i_ref+1,1)) * domin
            a2_g = (ag(l,i_ref,2)-ag(l,i_ref+1,2)) * domin
            a4_g = (ag(l,i_ref,4)-ag(l,i_ref+1,4)) * domin

            b1_g = ag(l,i_ref,1) - a1_g*rh_ref(i_ref)
            b2_g = ag(l,i_ref,2) - a2_g*rh_ref(i_ref)
            b4_g = ag(l,i_ref,4) - a4_g*rh_ref(i_ref)

            bsa = 0.
            baa = 0.
            ga  = 0.

            ! aerosol: boundary layer assumed for the lowest pbllev layers
            IF (k > pblind(j)) THEN  ! free troposphere aerosol
              bsa = (a1_sc*rhum(j,k) + b1_sc) * slf(j) + &
                (a2_sc*rhum(j,k) + b2_sc) * (1.-slf(j))
              baa = (a1_ab*rhum(j,k) + b1_ab) * slf(j) + &
                (a2_ab*rhum(j,k) + b2_ab) * (1.-slf(j))
              ga  = (a1_g *rhum(j,k) + b1_g)  * slf(j) + &
                (a2_g *rhum(j,k) + b2_g)  * (1.-slf(j))
            ELSE  ! mix marine/rural aerosol
              bsa = (a4_sc*rhum(j,k) + b4_sc)
              baa = (a4_ab*rhum(j,k) + b4_ab)
              ga  = (a4_g *rhum(j,k) + b4_g)
            ENDIF

            taer_sca(j,k,l) = bsa * part_col(j,k)
            taer_abs(j,k,l) = baa * part_col(j,k)
            gaer(j,k,l) = ga
          ENDDO
        ENDDO

      ENDDO

  END SUBROUTINE aero_2d

    !-------------------------------------------------------------------------

  SUBROUTINE slingo(clp, tau_sca, tau_abs, gcld, kproma_day, klev)

      ! A. Slingo's data for cloud particle radiative properties (from 'a gcm
      ! parameterization for the shortwave properties of water clouds' jas
      ! vol. 46 may 1989 pp 1419-1427)
      ! here only for the spectral range 250nm - 690 nm
    IMPLICIT NONE
      INTEGER ,                          INTENT(IN)  :: kproma_day, klev
      REAL, DIMENSION(kproma_day, klev), INTENT(IN)  :: clp     ! cloud liquid water path per layer [g/m^2]
      REAL, DIMENSION(kproma_day, klev), INTENT(OUT) :: tau_sca ! scattering optical depth
      REAL, DIMENSION(kproma_day, klev), INTENT(OUT) :: tau_abs ! absorption optical depth
      REAL, DIMENSION(kproma_day, klev), INTENT(OUT) :: gcld    ! asymmetry factor

      REAL, PARAMETER :: &
        abar =  2.817E-02, & ! a coefficient for extinction optical depth
        bbar =  1.305,     & ! b coefficient for extinction optical depth
        cbar = -5.62E-08,  & ! c coefficient for single particle scat albedo
        dbar =  1.63E-07,  & ! d coefficient for single particle scat albedo
        ebar =  0.829,     & ! e coefficient for asymmetry parameter
        fbar =  2.482E-03, & ! f coefficient for asymmetry parameter
        cer  =  10.

      REAL :: &
        tau, & ! total optical depth
        wc,  & ! single scattering albedo
        cwc    ! co single scattering albedo 1-wc

      INTEGER :: j, k
      ! set cloud extinction optical depth, single scatter albedo,
      ! asymmetry parameter, and forward scattered fraction:
      ! do not let single scatter albedo be 1; delta-eddington solution
      ! for non-conservative case:
      wc   = MIN(1.-cbar-dbar*cer, 0.999999)
      cwc  = cbar + dbar * cer
      DO k = 1, klev
        DO j = 1, kproma_day
          tau           = clp(j,k) * (abar + bbar / cer)
          gcld(j,k)     = ebar + fbar * cer
          tau_abs(j,k)  = cwc * tau
          tau_sca(j,k)  = wc * tau
        ENDDO
      ENDDO

  END SUBROUTINE slingo

    !-------------------------------------------------------------------------

  SUBROUTINE pifm(taus_cld, taua_cld, g_cld, taus_clr, taua_clr, g_clr, &
                  aclc, u0_in, albedo, fact, facth, kproma_day, klev)

      ! practical improved flux method (pifm)
      ! to calculate actinic fluxes
      ! Zdunkowski,Welch,Korb: Beitr. Phys. Atmosph. vol. 53, p. 147 ff

      ! This version is not suitable for calculation for conserving
      ! scattering (w0 = 1). w0 is limited to w0 <= 1. - 1.E-15.
      ! for w0 = 1, al(4) and al(5) have to be calculated differently.
    USE mach_pkg_messy_mod,  only: MAXWAV, fraclim, flux

    IMPLICIT NONE
    INTEGER ,                                 INTENT(IN)  :: kproma_day, klev
    REAL, DIMENSION(kproma_day, klev,MAXWAV), INTENT(IN)  :: taus_clr ! scattering optical depth (clear sky)
    REAL, DIMENSION(kproma_day, klev,MAXWAV), INTENT(IN)  :: taua_clr ! absorption optical depth (clear sky)
    REAL, DIMENSION(kproma_day, klev,MAXWAV), INTENT(IN)  :: g_clr    ! asymmetry factor of layer (clear sky)
    REAL, DIMENSION(kproma_day, klev),        INTENT(IN)  :: taus_cld ! scattering optical depth (slingo cloud parameter)
    REAL, DIMENSION(kproma_day, klev),        INTENT(IN)  :: taua_cld ! absorption optical depth (slingo cloud parameter)
    REAL, DIMENSION(kproma_day, klev),        INTENT(IN)  :: g_cld    ! asymmetry factor of layer(slingo cloud parameter)
    REAL, DIMENSION(kproma_day, klev),        INTENT(IN)  :: aclc     ! cloud fraction [1]
    REAL, DIMENSION(kproma_day, MAXWAV),      INTENT(IN)  :: albedo   ! surface albedo
    REAL, DIMENSION(kproma_day),              INTENT(IN)  :: u0_in    ! cosine of solar zenith angle
    REAL, DIMENSION(kproma_day, klev,MAXWAV), INTENT(OUT) :: fact     ! actinic flux [photons/(cm^2 s)]
    REAL, DIMENSION(kproma_day, klev, 3:5),   INTENT(OUT) :: facth    ! flux [photons/(cm^2 s)]

    REAL, DIMENSION(kproma_day, klev, 4) :: bb  ! initialization coefficients of pifm
    REAL, DIMENSION(kproma_day)          :: u0  ! Internally modified u0
    REAL, DIMENSION(kproma_day,klev) :: &
        tu1, tu2, tu3, tu4, tu5, tu6, tu7, tu8, tu9

    REAL :: &
        al(kproma_day,2*klev,5),   & ! matrix coefficient
        rw(kproma_day,3*(klev+1)), & ! flux array for cloudy sky
        rf(kproma_day,3*(klev+1)), & ! flux array for clear sky
        sd,                    & ! parallel solar flux
        fd,                    & ! downward diffuse flux
        fu                       ! upward diffuse flux

      INTEGER :: nlev3p1,k3
      REAL :: f,tautot,tauscat,w0,p_1,smooth1,smooth2
      REAL :: b0,bu0,alph1,alph2,alph3,alph4
      REAL :: eps,factor,ueps2,gam1,gam2,e,rm,ga,gg,ha,hb,hc,hd,gb,gc,gd
      REAL :: td1,td2,td3,td4,td5,td6,td7,tds1,tds2,tds3,tus1

      REAL, PARAMETER :: u = 2.      ! diffusivity factor
      REAL, PARAMETER :: delu0  = 1.E-3
      REAL, PARAMETER :: resonc = 1.E-6
      REAL, PARAMETER :: w0min  = 1.E-7
    INTEGER :: j, k, l

      !-----------------------------------------------------------------------

      nlev3p1 = 3 * (klev + 1)

      ! initiatization of pifm

      CALL pifmini(aclc, bb, kproma_day, klev)

      ! option: maximum overlap

      ! calculation of the matrix coefficients a1,...,a5
      u0 = u0_in ! First create a copy of the cosine sza
      DO l  = 1,MAXWAV      ! wavel. loop

        ! first: clear sky

        DO k  = 1,klev      ! altitude loop
          DO j  = 1,kproma_day         ! longitude loop

            tautot = taus_clr(j,k,l)+taua_clr(j,k,l)

            ! single scattering albedo

            IF (tautot > 0.) THEN
              w0 = taus_clr(j,k,l) / tautot
            ELSE
              w0 = 0.
            ENDIF

            w0 = MIN(w0,1.-w0min)

            ! p_1: first expansion coefficient of the phase function
            p_1 = 3.*g_clr(j,k,l)

            ! f: fraction of radiation contained in diffraction peak
            f = g_clr(j,k,l)**2

            ! b0:  fractional mean backward scattering coefficient
            ! of diffuse light
            ! bu0: backward scattering coefficient of primary scattered
            ! parallel solar light
            ! for small p_1 smooth1,smooth2 manage the smooth change of
            ! b0 and bu0 to 0

            IF (p_1<= 0.1) THEN
              smooth1 = 1.33333333-p_1*3.3333333
              smooth2 = 10.*p_1
            ELSE
              smooth1 = 1.
              smooth2 = 1.
            ENDIF

            b0 = (3.-p_1)/8.  *smooth1

            ! alpha coefficient
            alph1 = u*(1.-(1.-b0)*w0)
            alph2 = u*b0*w0

            ! epsilon and gamma coefficient
            eps = SQRT(alph1**2-alph2**2)
            factor = 1.-w0*f

            ! check for resonance condition in gam1 and gam2, if fulfil then
            ! chance u0(j) and calculate ueps2, bu0, alph3, alph4 again.
            ueps2 = (u0(j)*eps)**2
            IF (ABS(ueps2-factor**2)<resonc) THEN
              IF (ueps2<factor**2) THEN
                u0(j) = u0(j)-delu0
              ELSE
                u0(j) = u0(j)+delu0
              ENDIF
              ueps2 = (u0(j)*eps)**2
            ENDIF

            bu0 = 0.5-u0(j)/4.*(p_1-3.*f)/(1.-f)  *smooth2
            alph3 = w0*bu0*(1.-f)
            alph4 = w0*(1.-bu0)*(1.-f)

            gam1 = ( factor*alph3-u0(j)*(alph1*alph3+alph2*alph4) ) / &
              (factor**2-ueps2)
            gam2 = (-factor*alph4-u0(j)*(alph1*alph4+alph2*alph3) ) / &
              (factor**2-ueps2)

            e = EXP(-eps*tautot)
            rm = alph2/(alph1+eps)

            al(j,k,4) = e*(1.-rm**2)/(1.-e**2 * rm**2)
            al(j,k,5) = rm*(1.-e**2)/(1.-e**2 * rm**2)
            al(j,k,1) = EXP(-factor*tautot/u0(j))
            al(j,k,2) = -al(j,k,4)*gam2-al(j,k,5)*gam1*al(j,k,1) +gam2*al(j,k,1)
            al(j,k,3) = -al(j,k,5)*gam2-al(j,k,4)*gam1*al(j,k,1)+gam1

          ENDDO
        ENDDO

        ! second: cloudy sky

        DO k  = klev+1,2*klev ! altitude loop

          DO j = 1,kproma_day   ! longitude loop

            al(j,k,1) = al(j,k-klev,1)
            al(j,k,2) = al(j,k-klev,2)
            al(j,k,3) = al(j,k-klev,3)
            al(j,k,4) = al(j,k-klev,4)
            al(j,k,5) = al(j,k-klev,5)

          ENDDO

          DO j = 1,kproma_day    ! longitude loop

            IF (aclc(j,k-klev) <= fraclim) CYCLE

            tauscat = taus_clr(j,k-klev,l) + taus_cld(j,k-klev)
            tautot  = taua_clr(j,k-klev,l) + taua_cld(j,k-klev) + tauscat
            gg      = g_cld(j,k-klev)*taus_cld(j,k-klev) / tauscat

            IF (tautot > 0.) THEN
              w0 = tauscat/tautot
            ELSE
              w0 = 0.
            ENDIF

            w0 = MIN(w0,1.-w0min)

            p_1 = 3.*gg
            f = gg**2

            IF (p_1<= 0.1) THEN
              smooth1 = 1.33333333-p_1*3.3333333
              smooth2 = 10.*p_1
            ELSE
              smooth1 = 1.
              smooth2 = 1.
            ENDIF

            b0 = (3.-p_1)/8.  *smooth1

            alph1 = u*(1.-(1.-b0)*w0)
            alph2 = u*b0*w0

            eps = SQRT(alph1**2-alph2**2)
            factor = 1.-w0*f

            ueps2 = (u0(j)*eps)**2
            IF (ABS(ueps2-factor**2)<resonc) THEN
              IF (ueps2<factor**2) THEN
                u0(j) = u0(j)-delu0
              ELSE
                u0(j) = u0(j)+delu0
              ENDIF
              ueps2 = (u0(j)*eps)**2
            ENDIF

            bu0 = 0.5-u0(j)/4.*(p_1-3.*f)/(1.-f)  *smooth2
            alph3 = w0*bu0*(1.-f)
            alph4 = w0*(1.-bu0)*(1.-f)

            gam1 = ( factor*alph3-u0(j)*(alph1*alph3+alph2*alph4)) / &
              (factor**2-ueps2)
            gam2 = (-factor*alph4-u0(j)*(alph1*alph4+alph2*alph3)) / &
              (factor**2-ueps2)

            e = EXP(-eps*tautot)
            rm = alph2/(alph1+eps)

            al(j,k,4) = e*(1.-rm**2)/(1.-e**2 * rm**2)
            al(j,k,5) = rm*(1.-e**2)/(1.-e**2 * rm**2)
            al(j,k,1) = EXP(-factor*tautot/u0(j))
            al(j,k,2) = -al(j,k,4)*gam2-al(j,k,5)*gam1*al(j,k,1) + gam2*al(j,k,1)
            al(j,k,3) = -al(j,k,5)*gam2-al(j,k,4)*gam1*al(j,k,1) + gam1

          ENDDO
        ENDDO

        ! matrix inversion

        DO j  = 1,kproma_day ! longitude loop

          ! direct solution of the first four equations

          rf(j,1) = u0(j)*flux(l)
          rw(j,1) = 0.
          rf(j,2) = 0.
          rw(j,2) = 0.

          ! 5th to 10th equation: bring matrix elements on the left of the main
          ! diagonal to the rhs:  save elements on the right of the main
          ! diagonal in array -tu(l,1)

          rf(j,3) = al(j,1,3) * bb(j,1,1) * rf(j,1)
          rf(j,4) = al(j,1,1) * bb(j,1,1) * rf(j,1)
          rf(j,5) = al(j,1,2) * bb(j,1,1) * rf(j,1)
          rw(j,3) = al(j,1+klev,3) * (1.-bb(j,1,1)) * rf(j,1)
          rw(j,4) = al(j,1+klev,1) * (1.-bb(j,1,1)) * rf(j,1)
          rw(j,5) = al(j,1+klev,2) * (1.-bb(j,1,1)) * rf(j,1)

          tu1(j,1) = 0.
          tu2(j,1) = al(j,1,4)      * bb(j,1,2)
          tu3(j,1) = al(j,1,4)      * (1.-bb(j,1,4))
          tu4(j,1) = al(j,1+klev,4) * (1.-bb(j,1,2))
          tu5(j,1) = al(j,1+klev,4) * bb(j,1,4)
          tu6(j,1) = al(j,1,5)      * bb(j,1,2)
          tu7(j,1) = al(j,1,5)      * (1.-bb(j,1,4))
          tu8(j,1) = al(j,1+klev,5) * (1.-bb(j,1,2))
          tu9(j,1) = al(j,1+klev,5) * bb(j,1,4)


        ENDDO

        ! blocks of 6 equations: eliminate left matrix elements, save right
        ! matrix elements in array -tu(l,i), calculate rhs.

        DO k = 2,klev
          k3 = 3*k
          DO j  = 1,kproma_day ! longitude loop

            ha  = bb(j,k,1)*tu6(j,k-1) + (1.-bb(j,k,3))*tu8(j,k-1)
            hb  = bb(j,k,1)*tu7(j,k-1) + (1.-bb(j,k,3))*tu9(j,k-1)
            hc  = (1.-bb(j,k,1))*tu6(j,k-1) + bb(j,k,3)*tu8(j,k-1)
            hd  = (1.-bb(j,k,1))*tu7(j,k-1) + bb(j,k,3)*tu9(j,k-1)
            ga  = bb(j,k,1)*rf(j,k3-2) + (1.-bb(j,k,3))*rw(j,k3-2)
            gb  = bb(j,k,1)*rf(j,k3-1) + (1.-bb(j,k,3))*rw(j,k3-1)
            gc  = (1.-bb(j,k,1))*rf(j,k3-2) + bb(j,k,3)*rw(j,k3-2)
            gd  = (1.-bb(j,k,1))*rf(j,k3-1) + bb(j,k,3)*rw(j,k3-1)

            td1        = 1./(1.-al(j,k,5)*ha)
            rf(j,k3)   = td1*(al(j,k,3)*ga + al(j,k,5)*gb)
            tu1(j,k)   = td1*al(j,k,5)*hb
            tu2(j,k)   = td1*al(j,k,4)*bb(j,k,2)
            tu3(j,k)   = td1*al(j,k,4)*(1.-bb(j,k,4))
            td2        = al(j,k+klev,5)*hc
            td3        = 1./(1.-al(j,k+klev,5)*hd-td2*tu1(j,k))
            rw(j,k3)   = td3*(al(j,k+klev,3)*gc + &
                         al(j,k+klev,5)*gd+td2*rf(j,k3))
            tu4(j,k)   = td3*(al(j,k+klev,4)*(1.-bb(j,k,2)) + &
                         td2*tu2(j,k))
            tu5(j,k)   = td3*(al(j,k+klev,4)*bb(j,k,4) + &
                         td2*tu3(j,k))
            rf(j,k3+1) = al(j,k,1)*ga
            rw(j,k3+1) = al(j,k+klev,1)*gc
            td4        = al(j,k,4)*ha
            td5        = al(j,k,4)*hb
            rf(j,k3+2) = al(j,k,2)*ga +al(j,k,4)*gb+td4* &
                         rf(j,k3) + td5*rw(j,k3)
            tu6(j,k)   = al(j,k,5)*bb(j,k,2) + td4*tu2(j,k) + &
                         td5* tu4(j,k)
            tu7(j,k)   = al(j,k,5)*(1.-bb(j,k,4))+td4*tu3(j,k)+ &
                         td5*tu5(j,k)
            td6        = al(j,k+klev,4)*hc
            td7        = al(j,k+klev,4)*hd
            tu8(j,k)   = al(j,k+klev,5)*(1.-bb(j,k,2)) + &
                         td6*tu2(j,k) + td7*tu4(j,k)
            tu9(j,k)   = al(j,k+klev,5)*bb(j,k,4) + td6*tu3(j,k)+ &
                         td7*tu5(j,k)
            rw(j,k3+2) = al(j,k+klev,2)*gc + al(j,k+klev,4)*gd + &
                         td6*rf(j,k3) + td7*rw(j,k3)

          ENDDO
        ENDDO

        ! last two equations: the same as before

        DO j = 1,kproma_day ! longitude loop
          tds1          = 1. / (1.-albedo(j,l)*tu6(j,klev))
          rf(j,nlev3p1) = tds1 * (albedo(j,l) * rf(j,nlev3p1-2)+ &
                          albedo(j,l) * rf(j,nlev3p1-1))
          tus1          = tds1*albedo(j,l)*tu7(j,klev)
          tds2          = albedo(j,l)*tu8(j,klev)
          tds3          = 1./(1.-albedo(j,l)*tu9(j,klev)-tds2*tus1)
          rw(j,nlev3p1) = tds3*(albedo(j,l)*rw(j,nlev3p1-2) + &
                          albedo(j,l)*rw(j,nlev3p1-1) + tds2*rf(j,nlev3p1))
          rf(j,nlev3p1) = rf(j,nlev3p1) + tus1*rw(j,nlev3p1)
        ENDDO

        ! now we have created an upper triangular matrix the elements of which
        ! are -tu(l,i), 0, or 1 (in the main diagonal). the 0 and 1 elements
        ! are not stored in an array. let us solve the system now and store the
        ! results in the arrays rf (fluxes clear sky) and rw (fluxes cloudy sky)

        DO k = klev,1,-1
          k3 = 3*k
          DO j  = 1,kproma_day ! longitude loop

            rw(j,k3+2) = rw(j,k3+2) + tu8(j,k)*rf(j,k3+3) + &
                         tu9(j,k)*rw(j,k3+3)
            rf(j,k3+2) = rf(j,k3+2) + tu6(j,k)*rf(j,k3+3) + &
                         tu7(j,k)*rw(j,k3+3)
            rw(j,k3)   = rw(j,k3) + tu4(j,k)*rf(j,k3+3)   + &
                         tu5(j,k)*rw(j,k3+3)
            rf(j,k3)   = rf(j,k3) + tu2(j,k)*rf(j,k3+3)   + &
                         tu3(j,k)*rw(j,k3+3) + tu1(j,k)*rw(j,k3)

            sd = rf(j,k3+1) + rw(j,k3+1)
            fd = rf(j,k3+2) + rw(j,k3+2)
            fu = rf(j,k3+3) + rw(j,k3+3)

            ! actinic flux shall not be caculated at toa
            fact(j,k,l) = MAX(0.e0 , sd/u0(j) + u * fd + u * fu)
            IF (l>2 .AND. l<6) facth(j,k,l) = MAX(0.e0 ,sd+fd)
          ENDDO
        ENDDO

      ENDDO

  END SUBROUTINE pifm

  !-------------------------------------------------------------------------

  SUBROUTINE pifmini(aclc, bb, kproma_day, klev)
      ! initialization of practical improved flux method
      ! Zdunkowski et al 1982
      ! and
      ! J.F. Geleyn and A. Hollingsworth Contrib. to Atm. Phys. 52 no.1 p.1
      !-----------------------------------------------------------------------
      ! optical parameters of clouds and aerosol
      ! calculate continuity matrices for fractional cloudiness which find the
      ! flux going through a level i; levels are numbered from top to bottom.
      ! the following notation is used
      !
      ! s(cld,b) = bb(1)*s(clr,t)      + bb(3)*s(cld,t)
      ! s(clr,b) = (1-bb(1))*s(clr,t)  + (1-bb(3))*s(cld,t)
      ! fu(cld,t) = bb(2)*fu(clr,b)     + bb(4)*fu(cld,b)
      ! fu(fr,t) = (1-bb(2))*fu(clr,b) + (1-bb(4))*fu(cld,b)
      !
      ! where
      !
      ! cld = cloudy
      ! clr = clear
      ! b   = bottom of level k
      ! t   = top of level k
      ! s   = direct or diffuse downward flux
      ! fu  = diffuse upward flux
    IMPLICIT NONE
    INTEGER ,                             INTENT(IN)  :: kproma_day, klev
    REAL, DIMENSION(kproma_day, klev),    INTENT(IN)  :: aclc ! cloud fraction [1]
    REAL, DIMENSION(kproma_day, klev, 4), INTENT(OUT) :: bb   ! initialization coeffs.

    INTEGER :: j, k, km, kp

      ! first for highest level k = 1
      DO j  = 1,kproma_day ! longitude loop
        bb(j,1,1) = 1.-aclc(j,1)
        bb(j,1,2) = 1.
        bb(j,1,3) = 1.
        bb(j,1,4) = 1.
      ENDDO

      ! now for levels 2 to klev-1
      DO k = 2,klev-1
        km = k-1
        kp = k+1

        DO j  = 1,kproma_day ! longitude loop

          IF (aclc(j,km)<1.) THEN
            bb(j,k,1) = (1. - MAX(aclc(j,k),aclc(j,km))) / (1. - aclc(j,km))
          ELSE
            bb(j,k,1) = 1.
          ENDIF

          IF (aclc(j,km)>0.) THEN
            bb(j,k,3) = MIN(aclc(j,k),aclc(j,km)) / aclc(j,km)
          ELSE
            bb(j,k,3) = 1.
          ENDIF

          IF (aclc(j,kp)<1.) THEN
            bb(j,k,2) = (1. - MAX(aclc(j,k),aclc(j,kp))) / (1. - aclc(j,kp))
          ELSE
            bb(j,k,2) = 1.
          ENDIF

          IF (aclc(j,kp)>0.) THEN
            bb(j,k,4) =  MIN(aclc(j,k),aclc(j,kp)) / aclc(j,kp)
          ELSE
            bb(j,k,4) = 1.
          ENDIF

        ENDDO
      ENDDO

      ! finally for lowest level k = klev
      k  = klev
      km = klev-1

      DO j  = 1,kproma_day ! longitude loop

        IF (aclc(j,km)<1.) THEN
          bb(j,k,1) = (1. - MAX(aclc(j,k),aclc(j,km))) / (1. - aclc(j,km))
        ELSE
          bb(j,k,1) = 1.
        ENDIF

        IF (aclc(j,km)>0.) THEN
          bb(j,k,3) = MIN(aclc(j,k),aclc(j,km)) / aclc(j,km)
        ELSE
          bb(j,k,3) = 1.
        ENDIF

        bb(j,klev,2) = 1.
        bb(j,klev,4) = 1.

      ENDDO

  END SUBROUTINE pifmini

  ! **************************************************************************

  SUBROUTINE jval_cal_uv(fint, fact, facth, v3_du2, fhuv_2d, fhuvdna_2d, &
                         iu0, kproma_day, gni, klev)
    USE mach_pkg_messy_mod,  only: MAXWAV, p1, p2, p3
    USE chm_utils_mod,       only: dp

    IMPLICIT NONE
    INTEGER ,                                    INTENT(IN)  :: kproma_day, klev, gni
    REAL, DIMENSION(kproma_day, klev, 0:MAXWAV), INTENT(IN)  :: fint  ! integrated actinic flux
    REAL, DIMENSION(kproma_day, klev, MAXWAV),   INTENT(IN)  :: fact  ! actinic flux [photons/(cm^2 s)]
    REAL, DIMENSION(kproma_day, klev, 3:5),      INTENT(IN)  :: facth ! flux [photons/(cm^2 s)]
    REAL, DIMENSION(kproma_day, klev),           INTENT(IN)  :: v3_du2! slant O3 column [DU] (intervals 3-7)
    INTEGER,  DIMENSION(gni),                    INTENT(IN)  :: iu0
    REAL(DP), DIMENSION(gni, klev),              INTENT(OUT) :: fhuv_2d
    REAL(DP), DIMENSION(gni, klev),              INTENT(OUT) :: fhuvdna_2d

    REAL, DIMENSION(kproma_day, klev, 3:5) :: finth
    INTEGER :: j, k
    REAL    :: dj,djb
    REAL, DIMENSION(3:5) :: coeff_uv, coeff_dna
    REAL, PARAMETER :: a3_uv(3) = (/6.5612E-19, -2.3893E-24,  2.6419E-28/)
    REAL, PARAMETER :: a4_uv(2) = (/6.3890E-19, -6.9346E-25/)
    REAL, PARAMETER :: a5_uv(4) = (/1.3218E-19,-6.9323E-23, &
         1.3036E-26,-8.3734E-31/)
    REAL, PARAMETER :: a3_dna(3) = (/8.3302E-20, -2.6503E-23, &
         3.0102E-27/)
    REAL, PARAMETER :: a4_dna(3) = (/7.8083E-21,-2.2999E-24, 2.3155E-28/)
    REAL, PARAMETER :: a5_dna(4) = (/1.3206E-22,-7.8328E-26, &
         1.6011E-29,-1.0860E-33/)

    fhuv_2d(:,:)    = 0.0_dp
    fhuvdna_2d(:,:) = 0.0_dp

    DO k = 1,klev
       DO j = 1,kproma_day
          finth(j,k,3) = fint(j,k,3) * facth(j,k,3) / fact(j,k,3)
          finth(j,k,4) = fint(j,k,4) * facth(j,k,4) / fact(j,k,4)
          finth(j,k,5) = fint(j,k,5) * facth(j,k,5) / fact(j,k,5)
          coeff_uv(3)=p2(a3_uv(1),a3_uv(2),a3_uv(3),v3_du2(j,k))
          coeff_dna(3)=p2(a3_dna(1),a3_dna(2),a3_dna(3),v3_du2(j,k))
          coeff_uv(4)=p1(a4_uv(1),a4_uv(2),v3_du2(j,k))
          coeff_dna(4)=p2(a4_dna(1),a4_dna(2),a4_dna(3),v3_du2(j,k))
          coeff_uv(5)=p3(a5_uv(1),a5_uv(2),a5_uv(3),a5_uv(4),v3_du2(j,k))
          coeff_dna(5)=p3(a5_dna(1),a5_dna(2),a5_dna(3),a5_dna(4),v3_du2(j,k))
          dj=1.e4*(coeff_uv(3)*finth(j,k,3)+coeff_uv(4)*finth(j,k,4)+ &
               coeff_uv(5)*finth(j,k,5))
          djb=1.e4*(coeff_dna(3)*finth(j,k,3)+coeff_dna(4)*finth(j,k,4)+ &
               coeff_dna(5)*finth(j,k,5))
          fhuv_2d(iu0(j),k)=REAL(MAX(0.0,dj),dp)
          fhuvdna_2d(iu0(j),k)=REAL(MAX(0.0,djb),dp)
       ENDDO
    ENDDO

  END SUBROUTINE jval_cal_uv


END MODULE mach_messy_jval_mod

! ****************************************************************************
