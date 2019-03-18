MODULE DCMIP_2016_physics_module
!-------------------------------------------------------------------------
!Subroutines DCMIP_2016 (https://github.com/ClimateGlobalChange/DCMIP2016)
!-------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------
!
!  Version:  1.1
!
!  Date:  June 1st, 2016
!
!  Change log:
!
!  SUBROUTINE DCMIP2016_PHYSICS(test, u, v, p, pint, qsv, qsc, qsr, rho, theta, exner, t, dudt, dvdt, pdel_mid, &
!                           dt, z, zi, lat, wm, wp, nz, precl, pbl_type, prec_type)
!
!  Input variables:
!     test      (IN) DCMIP2016 test id (1,2,3)
!     u      (INOUT) zonal velocity on model levels (m/s)
!     v      (INOUT) meridional velocity on model levels (m/s)
!     p      (INOUT) total pressure on model levels (Pa)      !THERMO
!     pint   (INOUT) total pressure on interfaces levels (Pa) !MOMENTUM
!     qsv    (INOUT) specific humidity on model levels (kg/kg)
!     qsc    (INOUT) cloud water specific on model levels (kg/kg)
!     qsr    (INOUT) rain water specific on model levels (kg/kg)
!     rho    (INOUT) dry air density on model levels (kg/m^3)
!     theta  (INOUT) Potential temperature (K)
!     exner  (INOUT) Exner function (p/p0)**(R/cp)
!     t      (INOUT) Temperature (K)
!     dudt   (INOUT) Zonal wind tendency      (For GEM staggering)
!     dvdt   (INOUT) Meridional wind tendency (For GEM staggering)
!     pdel_mid  (IN) Layer thickness (Pa) (Thermo) For RJ large-scale prec.
!     dt        (IN) time step (s)
!     z      (INOUT) heights of model levels in the grid column (m)            !THERMO
!     zi     (INOUT) heights of model interfaces levels in the grid column (m) !MOMENTUM
!     lat       (IN) latitude of column
!     wm        (IN) weight Thermo to Momentum (above M)
!     wp        (IN) weight Thermo to Momentum (below M)
!     nz        (IN) number of levels in the column
!     precl     (IN) large-scale precip rate (m/s)
!     pbl_type  (IN) type of planetary boundary layer to use (0,1)
!                    0 = Default Reed-Jablonowski boundary layer
!                    1 = Modified Bryan boundary layer
!                   -1 = NONE
!     prec_type (IN) type of precipitation/microphysics to use (0,1)
!                    0 = Default Kessler physics routines
!                    1 = Reed-Jablonowski microphysics
!                   -1 = NONE
!
!  Authors: Paul Ullrich
!           University of California, Davis
!           Email: paullrich@ucdavis.edu
!
!           Kevin Reed
!           Stony Brook University
!           Email: kevin.a.reed@stonybrook.edu
!
!           Kessler is based on a code by Joseph Klemp
!           (National Center for Atmospheric Research)
!
!  Reference:
!
!    Klemp, J. B., W. C. Skamarock, W. C., and S.-H. Park, 2015:
!    Idealized Global Nonhydrostatic Atmospheric Test Cases on a Reduced
!    Radius Sphere. Journal of Advances in Modeling Earth Systems.
!    doi:10.1002/2015MS000435
!
!    Note: PRESSURE is not changed (RHO is changed) in DCMIP2016_PHYSICS when GEM
!
!=======================================================================

SUBROUTINE DCMIP2016_PHYSICS(test, u, v, p, pint, qsv, qsc, qsr, rho, theta, exner, t, dudt, dvdt, pdel_mid, &
                             dt, z, zi, lat, wm , wp, nz, precl, pbl_type, prec_type)

  use Kessler_module

  use vertical_interpolation, only: vertint2
  use dcmip_options

  IMPLICIT NONE

  !------------------------------------------------
  !   Arguments
  !------------------------------------------------

  INTEGER, INTENT(IN) :: &
            test         ! DCMIP2016 test index

  REAL(8), DIMENSION(nz), INTENT(INOUT) :: &
            u       ,  & ! Zonal velocity on model levels (m/s)
            v       ,  & ! Meridional velocity on model levels (m/s)
            qsv     ,  & ! Specific humidity (kg/kg)
            qsc     ,  & ! Cloud water specific (kg/kg)
            qsr          ! Rain water specific (kg/kg)

  REAL(8), DIMENSION(nz+1), INTENT(INOUT) :: &
            p            ! Pressure on model levels (Pa)     !THERMO

  REAL(8), DIMENSION(0:nz), INTENT(INOUT) :: &
            pint         ! Pressure on interface levels (Pa) !MOMENTUM

  REAL(8), DIMENSION(nz), INTENT(INOUT) :: &
            rho          ! Dry air density on model levels (kg/m^3)

  REAL(8), DIMENSION(nz), INTENT(INOUT) :: &
            theta        ! Potential temperature (K)

  REAL(8), DIMENSION(nz), INTENT(INOUT) :: &
            exner        ! Exner function (p/p0)**(R/cp)

  REAL(8), DIMENSION(nz), INTENT(INOUT) :: &
            t            ! Temperature (K)

  REAL(8), DIMENSION(nz), INTENT(INOUT) :: &
            dudt         ! Zonal wind tendency       (For GEM staggering)

  REAL(8), DIMENSION(nz), INTENT(INOUT) :: &
            dvdt         ! Meridional  wind tendency (For GEM staggering)

  REAL(8), DIMENSION(nz), INTENT(IN) :: &
            pdel_mid     ! Layer thickness (Pa) (Thermo) For RJ large-scale prec.

  REAL(8), INTENT(IN) :: &
            dt           ! Time step (s)

  REAL(8), DIMENSION(0:nz), INTENT(INOUT) :: &
            z       ,  & ! Heights of model levels (m)     !THERMO
            zi           ! Heights of model interfaces (m) !MOMENTUM

  REAL(8), INTENT(IN) :: &
            lat          ! Latitude of column (radians)

  REAL(8), DIMENSION(nz), INTENT(IN) :: &
            wm      ,  & ! Weight Thermo to Momentum (above M)
            wp           ! Weight Thermo to Momentum (below M)

  INTEGER, INTENT(IN) :: &
            nz           ! Number of levels in grid column

  REAL(8), INTENT(INOUT) :: &
            precl        ! Large-scale precip beneath the grid column (mm)

  INTEGER, INTENT(IN) :: &
            pbl_type,    & ! Type of planetary boundary layer to use
            prec_type      ! Type of precipitation to use

  !------------------------------------------------
  ! Physical Constants - MAY BE MODEL DEPENDENT
  !------------------------------------------------
  REAL(8), PARAMETER ::      &
    one      = 1.d0,         & ! One
    gravit  = 9.80616d0,     & ! Gravity (m/s^2)
    rair    = 287.d0,        & ! Gas constant for dry air (J/kg/K)
    cpair   = 1004.5d0,      & ! Specific heat of dry air (J/kg/K)
    latvap  = 2.5d6,         & ! Latent heat of vaporization (J/kg)
    rh2o    = 461.5d0,       & ! Gas constant for water vapor (J/kg/K)
    epsilo  = rair/rh2o,     & ! Ratio of gas constants for dry air to vapor
    zvir    = (rh2o/rair) - 1.d0, & ! Constant for virtual temp. calc. (~0.608)
    a       = 6371220.0,     & ! Reference Earth's Radius (m)
    omega   = 7.29212d-5,    & ! Reference rotation rate of the Earth (s^-1)
    pi      = 4.d0*atan(1.d0)  ! Pi

  !------------------------------------------------
  ! Local Constants for Simple Physics
  !------------------------------------------------
  REAL(8), PARAMETER ::      &
    C        = 0.0011d0,     & ! From Simth and Vogl 2008
    SST_TC   = 302.15d0,     & ! Constant Value for SST
    T0       = 273.16d0,     & ! Control temp for calculation of qsat
    e0       = 610.78d0,     & ! Saturation vapor pressure at T0
    rhow     = 1000.0d0,     & ! Density of Liquid Water
    Cd0      = 0.0007d0,     & ! Constant for Cd calc. Simth and Vogl 2008
    Cd1      = 0.000065d0,   & ! Constant for Cd calc. Simth and Vogl 2008
    Cm       = 0.002d0,      & ! Constant for Cd calc. Simth and Vogl 2008
    v20      = 20.0d0,       & ! Thresh. wind spd. for Cd Smith and Vogl 2008
    p0       = 100000.0d0,   & ! Constant for potential temp calculation
    pbltop   = 85000.d0,     & ! Top of boundary layer in p
    zpbltop  = 1000.d0,      & ! Top of boundary layer in z
    pblconst = 10000.d0,     & ! Constant for the decay of diffusivity
    T00      = 288.0d0,      & ! Horizontal mean T at surf. for moist baro test
    u0       = 35.0d0,       & ! Zonal wind constant for moist baro test
    latw     = 2.0d0*pi/9.0d0, & ! Halfwidth for  for baro test
    eta0     = 0.252d0,      & ! Center of jets (hybrid) for baro test
    etav     = (1.d0-eta0)*0.5d0*pi, & ! Auxiliary variable for baro test
    q0       = 0.021d0,      & ! Maximum specific humidity for baro test
    kappa    = 0.4d0           ! von Karman constant

  !------------------------------------------------
  ! Variables used in calculation
  !------------------------------------------------
  INTEGER :: k,kk

  REAL(8) ::                  &
    zat,                      & ! Altitude of lowest model level (m) !THERMO
    zam,                      & ! Altitude of lowest model level (m) !MOMENTUM
    Tsurf,                    & ! Sea surface temperature (K)
    ps,                       & ! Surface pressure (Pa)
    wind,                     & ! Wind speed in the lowest model level (m/s)
    Cd,                       & ! Drag coefficient for momentum
    qsat,                     & ! Saturation specific humidity
    qsats                       ! Saturation specific humidity at sea surface

  ! Matrix coefficients for PBL scheme
  REAL(8) ::                  &
    CAm,                      & ! Matrix coefficients for PBL scheme
    CAE,                      & ! Matrix coefficients for PBL scheme
    CCm,                      & ! Matrix coefficients for PBL scheme
    CCE                         ! Matrix coefficients for PBL scheme

  ! Variables defined on model levels
  REAL(8), DIMENSION(nz) ::   &
    rhom,                     & ! Moist density on model levels
 !!!qsv,                      & ! Specific humidity (kg/kg)
 !!!t,                        & ! Temperature (K)
 !!!theta,                    & ! Potential temperature on model levels
 !!!exner,                    & ! Exner function (p/p0)**(R/cp)
    CEm,                      & ! Matrix coefficients for PBL scheme
    CEE,                      & ! Matrix coefficients for PBL scheme
    CFu,                      & ! Matrix coefficients for PBL scheme
    CFv,                      & ! Matrix coefficients for PBL scheme
    CFt,                      & ! Matrix coefficients for PBL scheme
    CFq                         ! Matrix coefficients for PBL scheme

  ! Variables defined at interfaces
  REAL(8), DIMENSION(1:nz+1) :: &
    Km                          ! Eddy diffusivity for boundary layer (on THERMO   levels)

  REAL(8), DIMENSION(0:nz)   :: &
    Ke                          ! Eddy diffusivity for boundary layer (on MOMENTUM levels)

  REAL(8), DIMENSION(nz)     :: &
    qv                          ! Water vapor mixing ratio

  REAL(8), DIMENSION(nz)     :: &
    qc                          ! Cloud water mixing ratio

  REAL(8), DIMENSION(nz)     :: &
    qr                          ! Rain water mixing ratio

  REAL(4), DIMENSION(nz)     :: &
    t4,q4                       ! Storage (TEMPORAIRE)

  REAL(8), DIMENSION(nz)     :: &
    dtdt,dqdt                   ! Storage (TEMPORAIRE)

  REAL(4), DIMENSION(nz+1)   :: &
    log_p_t_4,log_p_m_4         ! Storage (TEMPORAIRE)

  REAL(4), DIMENSION(nz+1)   :: &
    zt_4,zm_4                   ! Storage (TEMPORAIRE)

  REAL(8) rpdel_int, rpdel_mid, dlnpint, tmp, rho_, tv_

  !Initialize tendencies
  !---------------------
  dudt = 0.0d0
  dvdt = 0.0d0
  dtdt = 0.0d0
  dqdt = 0.0d0

  !Conversion from Specific to Mixing ratio
  !----------------------------------------
  do k = 1,nz
     qv(k) = qsv(k)/(1.0d0 - qsv(k))
     qc(k) = qsc(k)/(1.0d0 - qsv(k))
     qr(k) = qsr(k)/(1.0d0 - qsv(k))
  enddo

  !------------------------------------------------
  ! Store altitude of lowest model level
  !------------------------------------------------
  dlnpint = log(pint(0)) - log(pint(1))
  zat = rair/gravit*t(1)*(1.d0 + zvir * qsv(1))*0.5d0*dlnpint
  zam = 2.0d0*zat

  !------------------------------------------------
  ! Calculate sea surface temperature
  !------------------------------------------------

  ! Moist baroclinic wave test
  if (test == 1) then
    Tsurf = &
      (T00 + pi*u0/rair * 1.5d0 * sin(etav) * (cos(etav))**0.5d0 * &
      ((-2.d0*(sin(lat))**6 * ((cos(lat))**2 + 1.d0/3.d0) &
        + 10.d0/63.d0) * &
      u0 * (cos(etav))**1.5d0  + &
      (8.d0/5.d0*(cos(lat))**3 * ((sin(lat))**2 + 2.d0/3.d0) &
        - pi/4.d0)*a*omega*0.5d0 ))/ &
      (1.d0+zvir*q0*exp(-(lat/latw)**4))

  ! Tropical cyclone test
  elseif (test == 2) then
    Tsurf = SST_TC

  ! Supercell test
  elseif (test == 3) then
    Tsurf = 0.d0

  ! Invalid test
  else
    write(*,*) 'Invalid test specified in DCMIP2016_PHYSICS', test
    stop
  endif

  !------------------------------------------------
  ! Initialize precipitation rate to zero
  !------------------------------------------------
  precl = 0.d0

  !------------------------------------------------
  ! Calculate moist density
  !------------------------------------------------
  do k = 1, nz
    rhom(k) = rho(k) * (1.0 + qv(k))
  enddo

  !--------------------------------------------------------------------------
  ! Large-scale precipitation (Reed-Jablonowski): As in SIMPLIFIED_PHYSICS_V6
  !--------------------------------------------------------------------------
  if (prec_type == 1) then

    qsc = 0.0d0
    qsr = 0.0d0
    qc  = 0.0d0
    qr  = 0.0d0

    !Calculate Tendencies
    !---------------------
    do k=nz,1,-1 !TOP to BOTTOM
      qsat = epsilo*e0/p(k)*exp(-latvap/rh2o*((1.0d0/t(k))-1.0d0/T0))
      if (qsv(k) > qsat) then
        tmp  = 1.0d0/dt*(qsv(k)-qsat)/(1.0d0+(latvap/cpair)*(epsilo*latvap*qsat/(rair*t(k)**2)))
        dtdt(k) = dtdt(k)+latvap/cpair*tmp
        dqdt(k) = dqdt(k)-tmp
        precl = precl + tmp*pdel_mid(k)/(gravit*rhow)
      end if
    end do

    !Update moisture and temperature fields from Larger-Scale Precipitation Scheme
    !-----------------------------------------------------------------------------
    do k=1,nz
      t(k) = t(k) + dtdt(k)*dt
      qsv(k) = qsv(k) + dqdt(k)*dt
      qv(k) = qsv(k) / (1.0 - qsv(k))
    end do

  !------------------------------------------------
  ! Large-scale precipitation (Kessler)
  !------------------------------------------------
  elseif (prec_type == 0) then

    do k=1,nz
      qv(k) = max(qv(k),0.0d0)
      qc(k) = max(qc(k),0.0d0)
      qr(k) = max(qr(k),0.0d0)
    enddo

    CALL KESSLER(   &
      theta,        &
      qv,           &
      qc,           &
      qr,           &
      rho,          &
      exner,        &
      dt,           &
      z(1),         & !Above BOTTOM
      nz,           &
      precl)

    !Conversion from Mixing ratio to Specific
    !----------------------------------------

    do k = 1,nz
       qsv(k)= qv(k)/(1.0d0 + qv(k))
       qsc(k)= qc(k)/(1.0d0 + qv(k))
       qsr(k)= qr(k)/(1.0d0 + qv(k))
    enddo

    !Convert theta to temperature
    !----------------------------
    do k = 1,nz
      t(k)  = theta(k) * exner(k)
      rho(k)= p(k)/(rair*t(k)*(one + zvir * qsv(k))*(one + qv(k)))
      rhom(k) = rho(k) / (one - qsv(k))
    enddo

  elseif (prec_type /= -1) then
    write(*,*) 'Invalid prec_type specified in DCMIP2016_PHYSICS', prec_type
    stop
  endif

  !----------------------------------------------------
  ! Do not apply surface fluxes or PBL for supercell
  !----------------------------------------------------
  if (test == 3) then
    return
  endif

  !----------------------------------------------------
  ! Do not apply surface fluxes or PBL
  !----------------------------------------------------
  if (pbl_type == -1) return

  !----------------------------------------------------
  ! Reset Heights of model levels (m)
  !----------------------------------------------------

     !Bottom
     !------
     z (0) = 0. !ASSUME NO TOPOGRAPHY
     zi(0) = 0. !ASSUME NO TOPOGRAPHY

     !Estimate zi
     !-----------
     do kk = 1,nz
        dlnpint = log(pint(kk-1)) - log(pint(kk))
        zi(kk) = zi(kk-1) + rair/gravit*t(kk)*(1.d0 + zvir * qsv(kk))*dlnpint
     enddo

     !Estimate zt
     !-----------
     log_p_t_4(nz+1) = log(pint(0))
     log_p_m_4(nz+1) = log(pint(0))
     zm_4(nz+1) = zi(0)

     do k=nz,1,-1

        kk = nz-k+1

        log_p_t_4(k) = log(   p(kk))
        log_p_m_4(k) = log(pint(kk))

        zm_4(k) = zi(kk)

     end do

     call vertint2 ( zt_4, log_p_t_4, nz, zm_4, log_p_m_4, nz+1, &
                        1,1,1,1,1,1,1,1 )

     do k=nz,1,-1

        kk = nz-k+1

        z(kk) = zt_4(k)

     end do

  !------------------------------------------------
  ! Store altitude of lowest model level
  !------------------------------------------------
  dlnpint = log(pint(0)) - log(pint(1))
  zat = rair/gravit*t(1)*(1.d0 + zvir * qsv(1))*0.5d0*dlnpint
  zam = 2.0d0*zat

  !------------------------------------------------
  ! Turbulent mixing coefficients
  !------------------------------------------------
  wind = sqrt(u(1)**2 + v(1)**2)

  if (wind < v20) then
    Cd = Cd0 + Cd1 * wind
  else
    Cd = Cm
  endif

  Km = 0.
  Ke = 0.

  ! Reed-Jablonowski Boundary layer
  if (pbl_type == 0) then
    do k = 1, nz

      ! ----------------------------
      ! Km Thermo
      ! ----------------------------
      if (p(k) >= pbltop) then
        Km(k) = Cd * wind * zam
      else
        Km(k) = Cd * wind * zam * exp(-(pbltop-p(k))**2/pblconst**2)
      endif

      ! ----------------------------
      ! Ke Momentum
      ! ----------------------------
      if (pint(k) >= pbltop) then
        Ke(k) = C  * wind * zat
      else
        Ke(k) = C  * wind * zat * exp(-(pbltop-pint(k))**2/pblconst**2)
      endif

    enddo

  ! Bryan Planetary Boundary Layer
  elseif (pbl_type == 1) then
    do k = 1, nz

      ! ----------------------------
      ! Km Thermo
      ! ----------------------------
      if (z(k) <= zpbltop) then
        Km(k) = kappa * sqrt(Cd) * wind * z(k) &
              * (one - z(k)/zpbltop) * (one - z(k)/zpbltop)
      else
        Km(k) = 0.d0
      endif

      ! ----------------------------
      ! Ke Momentum
      ! ----------------------------
      if (zi(k) <= zpbltop) then
        Ke(k) = kappa * sqrt(C) * wind * zi(k) &
              * (one - zi(k)/zpbltop) * (one - zi(k)/zpbltop)
      else
        Ke(k) = 0.d0
      endif

    enddo

  ! Invalid PBL
  else
    write(*,*) 'Invalid pbl_type specified in DCMIP2016_PHYSICS', pbl_type
    stop
  endif

  !------------------------------------------------
  ! Surface fluxes
  !------------------------------------------------

  ! Hydrostatic surface pressure
  ps = pint(0)
  qsats = epsilo * e0 / ps * exp(-latvap / rh2o * ((one/Tsurf)-(one/T0)))

  dudt(1) = dudt(1) + (u(1) / (one + Cd * wind * dt / zam) - u(1))/dt
  dvdt(1) = dvdt(1) + (v(1) / (one + Cd * wind * dt / zam) - v(1))/dt
  u(1) = u(1) / (one + Cd * wind * dt / zam)
  v(1) = v(1) / (one + Cd * wind * dt / zam)
  qsv(1) = (qsv(1) + C * wind * qsats * dt / zat) / (one + C * wind * dt / zat)
  t(1) = (t(1) + C * wind * Tsurf * dt / zat) / (one + C * wind * dt / zat)

  qv(1) = qsv(1) / (one - qsv(1))

  rho(1)= p(1)/(rair*t(1)*(one + zvir * qsv(1))*(one + qv(1)))
  rhom(1) = rho(1) / (1.0 - qsv(1))

  !------------------------------------------------
  ! Boundary layer
  !------------------------------------------------

  !
  ! print *,'VIV: calc tstag,qstag'
  do k=2, nz
    t4(k-1)= wp(k-1)*t(k-1) + wm(k-1)*t(k)
    q4(k-1)= wp(k-1)*qsv(k-1) + wm(k-1)*qsv(k)
  enddo

  !-----------------------------------------------------
  !AS IN SIMPLIFIED_PHYSICS_V6
  !-----------------------------------------------------
  do k = 1, nz

    CAE = 0.d0
    CCE = 0.d0
    CAm = 0.d0
    CCm = 0.d0

    ! A,C for T,Q ! So Ke Momentum
    ! ----------------------------
    rpdel_mid = pint(k-1) - pint(k)
    rpdel_mid = 1.0d0/rpdel_mid

    if (k/=1) then
       tv_ = t4(k-1)*(1.0d0+zvir*q4(k-1))
       rho_= (pint(k-1)/(rair*tv_))
       CAE = rpdel_mid*dt*gravit*gravit*Ke(k-1)*rho_*rho_/(p(k-1)-p(k))
    endif

    if (k/=nz) then
       tv_ = t4(k)*(1.0d0+zvir*q4(k))
       rho_= (pint(k)/(rair*tv_))
       CCE = rpdel_mid*dt*gravit*gravit*Ke(k)*rho_*rho_/(p(k)-p(k+1))
    endif

    ! A,C for U,V ! So Km Thermo
    ! ----------------------------
    rpdel_int = p(k) - p(k+1)
    rpdel_int = 1.0d0/rpdel_int
    if (k/=1) then
    rho_= (p(k)/(rair*t(k)*(1.0d0+zvir*qsv(k))))
    CAm = rpdel_int*dt*gravit*gravit*Km(k)*rho_*rho_/(pint(k-1)-pint(k))
    endif
    if (k/=nz) then
    rho_= (p(k+1)/(rair*t(k+1)*(1.0d0+zvir*qsv(k+1))))
    CCm = rpdel_int*dt*gravit*gravit*Km(k+1)*rho_*rho_/(pint(k)-pint(k+1))
    endif

    if (k == 1) then
       CEE(k) = CCE/(1.0d0+CAE+CCE)
       CEm(k) = CCm/(1.0d0+CAm+CCm)
       CFu(k) = u(k)/(1.0d0+CAm+CCm)
       CFv(k) = v(k)/(1.0d0+CAm+CCm)
       CFt(k) = ((p0/p(k))**(rair/cpair)*t(k))/(1.0d0+CAE+CCE)
       CFq(k) = qsv(k)/(1.0d0+CAE+CCE)
    else
       CEE(k) = CCE/(1.0d0+CAE+CCE-CAE*CEE(k-1))
       CEm(k) = CCm/(1.0d0+CAm+CCm-CAm*CEm(k-1))
       CFu(k) = (u(k)+CAm*CFu(k-1))/(1.0d0+CAm+CCm-CAm*CEm(k-1))
       CFv(k) = (v(k)+CAm*CFv(k-1))/(1.0d0+CAm+CCm-CAm*CEm(k-1))
       CFt(k) = ((p0/p(k))**(rair/cpair)*t(k)+CAE*CFt(k-1))/(1.0d0+CAE+CCE-CAE*CEE(k-1))
       CFq(k) = (qsv(k)+CAE*CFq(k-1))/(1.0d0+CAE+CCE-CAE*CEE(k-1))
    endif

  enddo

  !Calculate the updated temperaure and specific humidity and wind tendencies
  !--------------------------------------------------------------------------

  !First we need to calculate the tendencies at the top model level
  !----------------------------------------------------------------

  dudt(nz) = dudt(nz)+(CFu(nz)-u(nz))/dt
  dvdt(nz) = dvdt(nz)+(CFv(nz)-v(nz))/dt
  u(nz)    = CFu(nz)
  v(nz)    = CFv(nz)
  t(nz)    = CFt(nz)*(p(nz)/p0)**(rair/cpair)
  qsv(nz)  = CFq(nz)

  do k=nz-1,1,-1
     dudt(k) = dudt(k)+(CEm(k)*u(k+1)+CFu(k)-u(k))/dt
     dvdt(k) = dvdt(k)+(CEm(k)*v(k+1)+CFv(k)-v(k))/dt
     u(k)    = CEm(k)*u(k+1)+CFu(k)
     v(k)    = CEm(k)*v(k+1)+CFv(k)
     t(k)    =(CEE(k)*t(k+1)*(p0/p(k+1))**(rair/cpair)+CFt(k))*(p(k)/p0)**(rair/cpair)
     qsv(k)  = CEE(k)*qsv(k+1)+CFq(k)
  end do

  do k = 1,nz
    qv(k) = qsv(k) / (one - qsv(k))
  enddo

  !Conversion from Mixing ratio to Specific
  !----------------------------------------

  do k = 1,nz
     qsc(k)= qc(k)/(1.0d0 + qv(k))
     qsr(k)= qr(k)/(1.0d0 + qv(k))
  enddo

  return

END SUBROUTINE DCMIP2016_PHYSICS

END MODULE DCMIP_2016_physics_module
