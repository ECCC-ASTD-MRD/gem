MODULE DCMIP_2016_physics_module

!-----------------------------------------------------------------------------!
!Subroutines DCMIP_2016 (https://doi.org/10.5281/zenodo.1298671)              !
!-----------------------------------------------------------------------------!
!GEM's adaptation to DCMIP_2016 physics:                                      !
!1) Combines elements from 2 subroutines: simple_physics_v6/dcmip_physics_z_v1!
!2) PRESSURE not changed ((RHO is changed)                                    !
!3) Vertical index increases from BOTTOM to TOP                               !
!-----------------------------------------------------------------------------!

CONTAINS

!------------------------------------------------------------------------------
!
!  Version:  1.2
!
!  Date:  July 28th, 2017
!
!  Change log:
!
!  SUBROUTINE DCMIP2016_PHYSICS(test, u, v, pt, pm, qsv, qsc, qsr, t, &
!                           dt, lat, nz, precl, pbl_type, prec_type)
!
!  Input variables:
!     test      (IN) DCMIP2016 test index of SST (1,2,3)
!     u      (INOUT) Zonal velocity (m/s)
!     v      (INOUT) Meridional velocity (m/s)
!     pt        (IN) Pressure (Pa) !THERMO
!     pm        (IN) Pressure (Pa) !MOMENTUM
!     qsv    (INOUT) Specific humidity (kg/kg)
!     qsc    (INOUT) Cloud water specific (kg/kg)
!     qsr    (INOUT) Rain water specific (kg/kg)
!     t      (INOUT) Temperature (K)
!     dt        (IN) Time step (s)
!     lat       (IN) Latitude of column (radians)
!     nz        (IN) Number of levels in the column
!     precl    (OUT) Large-scale precip rate (m/s)
!     pbl_type  (IN) Type of planetary boundary layer
!                    ------------------------------------------------
!                     0 = Reed-Jablonowski Boundary layer
!                     1 = Georges Bryan Boundary layer
!                    -1 = NONE
!                    ------------------------------------------------
!     prec_type (IN) Type of precipitation/microphysics
!                    ------------------------------------------------
!                     0 = Kessler Microphysics
!                     1 = Reed-Jablonowski Large-scale precipitation
!                    -1 = NONE
!                    ------------------------------------------------
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
!=======================================================================

SUBROUTINE DCMIP2016_PHYSICS(test, u, v, pt, pm, qsv, qsc, qsr, t, &
                             dt, lat, nz, precl, pbl_type, prec_type)

  use Kessler_module

  implicit none

  !------------------------------------------------
  !   Arguments
  !------------------------------------------------

  integer, intent(in) :: &
            test         ! DCMIP2016 test index of SST (1,2,3)

  integer, intent(in) :: &
            nz           ! Number of levels in grid column

  real(kind=REAL64), dimension(nz), intent(inout) :: &
            u       ,  & ! Zonal velocity (m/s)
            v       ,  & ! Meridional velocity (m/s)
            qsv     ,  & ! Specific humidity (kg/kg)
            qsc     ,  & ! Cloud water specific (kg/kg)
            qsr          ! Rain water specific (kg/kg)

  real(kind=REAL64), dimension(nz+1),intent(in) :: &
            pt           ! Pressure (Pa) !THERMO

  real(kind=REAL64), dimension(0:nz),intent(in) :: &
            pm           ! Pressure (Pa) !MOMENTUM

  real(kind=REAL64), dimension(nz), intent(inout) :: &
            t            ! Temperature (K)

  real(kind=REAL64), intent(in) :: &
            dt           ! Time step (s)

  real(kind=REAL64), intent(in) :: &
            lat          ! Latitude of column (radians)

  real(kind=REAL64), intent(inout) :: &
            precl        ! Large-scale precip beneath the grid column (mm)

  integer, intent(in) :: &
            pbl_type,  & ! Type of planetary boundary layer
            prec_type    ! Type of precipitation/microphysics

  !------------------------------------------------
  ! Physical Constants - MAY BE MODEL DEPENDENT
  !------------------------------------------------
  real(kind=REAL64), parameter ::      &
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
  real(kind=REAL64), parameter ::      &
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
  integer :: k

  real(kind=REAL64) ::                  &
    zat,                      & ! Altitude of lowest model level (m) !THERMO
    zam,                      & ! Altitude of lowest model level (m) !MOMENTUM
    Tsurf,                    & ! Sea surface temperature (K)
    ps,                       & ! Surface pressure (Pa)
    wind,                     & ! Wind speed in the lowest model level (m/s)
    Cd,                       & ! Drag coefficient for momentum
    qsat,                     & ! Saturation specific humidity
    qsats                       ! Saturation specific humidity at sea surface

  ! Matrix coefficients for PBL scheme
  real(kind=REAL64) ::                  &
    CAm,                      & ! Matrix coefficients for PBL scheme
    CAE,                      & ! Matrix coefficients for PBL scheme
    CCm,                      & ! Matrix coefficients for PBL scheme
    CCE                         ! Matrix coefficients for PBL scheme

  ! Matrix coefficients for PBL scheme (Array)
  real(kind=REAL64), dimension(nz) ::   &
    CEm,                      & ! Matrix coefficients for PBL scheme
    CEE,                      & ! Matrix coefficients for PBL scheme
    CFu,                      & ! Matrix coefficients for PBL scheme
    CFv,                      & ! Matrix coefficients for PBL scheme
    CFt,                      & ! Matrix coefficients for PBL scheme
    CFq                         ! Matrix coefficients for PBL scheme

  ! Eddy diffusivity for boundary layer
  real(kind=REAL64), dimension(nz)     :: &
    Km                          ! THERMO

  real(kind=REAL64), dimension(nz)     :: &
    Ke                          ! MOMENTUM

  real(kind=REAL64), dimension(nz)     :: &
    theta,                    & ! Potential temperature
    exner                       ! Exner function (p/p0)**(R/cp)

  real(kind=REAL64), dimension(nz)     :: &
    qv                          ! Water vapor mixing ratio

  real(kind=REAL64), dimension(nz)     :: &
    qc                          ! Cloud water mixing ratio

  real(kind=REAL64), dimension(nz)     :: &
    qr                          ! Rain water mixing ratio

  real(kind=REAL64), dimension(0:nz)   :: &
    zt, &                       ! Heights (m) !THERMO
    zm                          ! Heights (m) !MOMENTUM

  real(kind=REAL64), dimension(nz)     :: &
   wm, &                        ! Weight int. THERMO to MOMENTUM (above M)
   wp                           ! Weight int. THERMO to MOMENTUM (below M)

  real(kind=REAL64), dimension(nz)     :: &
    t_m,qsv_m                   ! Temperature/Specific Humidity MOMENTUM

  real(kind=REAL64), dimension(nz)     :: &
    dtdt,dqdt                   ! Tendencies used in RJ large-scale prec.

  real(kind=REAL64), dimension(nz)     :: &
    rho                         ! Dry air density used in Kessler large-scale prec.

  real(kind=REAL64), dimension(nz)     :: &
    dp_t                        ! Layer thickness (Pa) (Thermo) used in RJ large-scale prec.

  real(kind=REAL64) :: r_dp_m, r_dp_t, dlnp_m, dlnp_t, rho_m, rho_t, tv_m, tmp

  !Conversion from Specific to Mixing ratio
  !----------------------------------------
  do k=1,nz
     qv(k) = qsv(k)/(1.0d0 - qsv(k))
     qc(k) = qsc(k)/(1.0d0 - qsv(k))
     qr(k) = qsr(k)/(1.0d0 - qsv(k))
  end do

  !Coefficients for linear interpolating from THERMO to MOMENTUM
  !wm=Above M /wp=Below M
  !-------------------------------------------------------------
  do k=1,nz
     wp(k) = (pm(k) - pt(k+1)) / (pt(k) - pt(k+1))
     wm(k) = 1.0d0 - wp(k)
  end do

  !Estimate Temperature/Specific Humidity MOMENTUM
  !-----------------------------------------------
  do k=1,nz-1
       t_m(k)= wp(k)*  t(k) + wm(k)*  t(k+1)
     qsv_m(k)= wp(k)*qsv(k) + wm(k)*qsv(k+1)
  end do

  !Estimate zm = Heights MOMENTUM
  !------------------------------
  zm(0) = 0.

  do k=1,nz
     dlnp_t = log(pm(k-1)) - log(pm(k))
     zm(k) = zm(k-1) + rair/gravit*t(k)*(1.d0 + zvir * qsv(k))*dlnp_t
  end do

  !Estimate zt = Heights THERMO
  !----------------------------
  zt(0) = 0.
  zt(1) = 0.5 * zm(1)

  do k=2,nz
     dlnp_m = log(pt(k-1)) - log(pt(k))
     zt(k) = zt(k-1) + rair/gravit*t_m(k-1)*(1.d0 + zvir * qsv_m(k-1))*dlnp_m
  end do

  !------------------------------------------------
  ! Store altitude of lowest model level
  !------------------------------------------------
  zat = zt(1)
  zam = zm(1)

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
  end if

  !------------------------------------------------
  ! Initialize precipitation rate to zero
  !------------------------------------------------
  precl = 0.d0

  !-------------------------------------------------------------------------
  ! Reed-Jablonowski Large-scale precipitation (As in SIMPLIFIED_PHYSICS_V6)
  !-------------------------------------------------------------------------
  if (prec_type == 1) then

    qsc = 0.0d0
    qsr = 0.0d0
    qc  = 0.0d0
    qr  = 0.0d0

    dtdt = 0.0d0
    dqdt = 0.0d0

    do k=1,nz
       dp_t(k) =  pm(k-1) - pm(k)
    end do

    !Calculate Tendencies
    !---------------------
    do k=nz,1,-1 !TOP to BOTTOM
      qsat = epsilo*e0/pt(k)*exp(-latvap/rh2o*((1.0d0/t(k))-1.0d0/T0))
      if (qsv(k) > qsat) then
        tmp  = 1.0d0/dt*(qsv(k)-qsat)/(1.0d0+(latvap/cpair)*(epsilo*latvap*qsat/(rair*t(k)**2)))
        dtdt(k) = dtdt(k)+latvap/cpair*tmp
        dqdt(k) = dqdt(k)-tmp
        precl = precl + tmp*dp_t(k)/(gravit*rhow)
      end if
    end do

    !Update moisture and temperature fields from Larger-Scale Precipitation Scheme
    !-----------------------------------------------------------------------------
    do k=1,nz
      t(k) = t(k) + dtdt(k)*dt
      qsv(k) = qsv(k) + dqdt(k)*dt
      qv(k) = qsv(k) / (1.0 - qsv(k))
    end do

  !---------------------
  ! Kessler Microphysics
  !---------------------
  elseif (prec_type == 0) then

    do k=1,nz
      qv(k) = max(qv(k),0.0d0)
      qc(k) = max(qc(k),0.0d0)
      qr(k) = max(qr(k),0.0d0)
    end do

    !Estimate dry air density
    !------------------------
    do k=1,nz
       rho(k) = pt(k)/(rair*t(k)*(one + zvir * qsv(k))*(one + qv(k)))
    end do

    !Estimate Exner function and Potential temperature
    !-------------------------------------------------
    do k=1,nz
       exner(k) = (pt(k) / p0)**(rair/cpair)
       theta(k) = t(k) / exner(k)
    end do

    call kessler(   &
      theta,        &
      qv,           &
      qc,           &
      qr,           &
      rho,          &
      exner,        &
      dt,           &
      zt(1),        & !Above BOTTOM
      nz,           &
      precl)

    !Conversion from Mixing ratio to Specific
    !----------------------------------------
    do k=1,nz
       qsv(k)= qv(k)/(1.0d0 + qv(k))
       qsc(k)= qc(k)/(1.0d0 + qv(k))
       qsr(k)= qr(k)/(1.0d0 + qv(k))
    end do

    !Convert theta to temperature
    !----------------------------
    do k=1,nz
      t(k)  = theta(k) * exner(k)
    end do

  elseif (prec_type /= -1) then
    write(*,*) 'Invalid prec_type specified in DCMIP2016_PHYSICS', prec_type
    stop
  end if

  !----------------------------------------------------
  ! Do not apply surface fluxes or PBL for supercell
  !----------------------------------------------------
  if (test == 3) then
    return
  end if

  !----------------------------------------------------
  ! Do not apply surface fluxes or PBL
  !----------------------------------------------------
  if (pbl_type == -1) return

  !Estimate Temperature/Specific Humidity MOMENTUM
  !-----------------------------------------------
  do k=1,nz-1
       t_m(k)= wp(k)*  t(k) + wm(k)*  t(k+1)
     qsv_m(k)= wp(k)*qsv(k) + wm(k)*qsv(k+1)
  end do

  !Estimate zm = Heights MOMENTUM
  !------------------------------
  zm(0) = 0.

  do k=1,nz
     dlnp_t = log(pm(k-1)) - log(pm(k))
     zm(k) = zm(k-1) + rair/gravit*t(k)*(1.d0 + zvir * qsv(k))*dlnp_t
  end do

  !Estimate zt = Heights THERMO
  !----------------------------
  zt(0) = 0.
  zt(1) = 0.5 * zm(1)

  do k=2,nz
     dlnp_m = log(pt(k-1)) - log(pt(k))
     zt(k) = zt(k-1) + rair/gravit*t_m(k-1)*(1.d0 + zvir * qsv_m(k-1))*dlnp_m
  end do

  !------------------------------------------------
  ! Store altitude of lowest model level
  !------------------------------------------------
  zat = zt(1)
  zam = zm(1)

  !------------------------------------------------
  ! Turbulent mixing coefficients
  !------------------------------------------------
  wind = sqrt(u(1)**2 + v(1)**2)

  if (wind < v20) then
    Cd = Cd0 + Cd1 * wind
  else
    Cd = Cm
  end if

  Km = 0.
  Ke = 0.

  !Reed-Jablonowski Boundary layer
  !-------------------------------
  if (pbl_type == 0) then

    do k=1,nz

      !Km Thermo
      !---------
      if (pt(k) >= pbltop) then
        Km(k) = Cd * wind * zam
      else
        Km(k) = Cd * wind * zam * exp(-(pbltop-pt(k))**2/pblconst**2)
      end if

      !Ke Momentum
      !-----------
      if (pm(k) >= pbltop) then
        Ke(k) = C  * wind * zat
      else
        Ke(k) = C  * wind * zat * exp(-(pbltop-pm(k))**2/pblconst**2)
      end if

    end do

  !Georges Bryan Boundary Layer
  !----------------------------
  elseif (pbl_type == 1) then

    do k=1,nz

      !Km Thermo
      !---------
      if (zt(k) <= zpbltop) then
        Km(k) = kappa * sqrt(Cd) * wind * zt(k) &
              * (one - zt(k)/zpbltop) * (one - zt(k)/zpbltop)
      else
        Km(k) = 0.d0
      end if

      !Ke Momentum
      !-----------
      if (zm(k) <= zpbltop) then
        Ke(k) = kappa * sqrt(C) * wind * zm(k) &
              * (one - zm(k)/zpbltop) * (one - zm(k)/zpbltop)
      else
        Ke(k) = 0.d0
      end if

    end do

  ! Invalid PBL
  else
    write(*,*) 'Invalid pbl_type specified in DCMIP2016_PHYSICS', pbl_type
    stop
  end if

  !------------------------------------------------
  ! Surface fluxes
  !------------------------------------------------

  !Hydrostatic surface pressure
  !----------------------------
  ps = pm(0)

  qsats = epsilo * e0 / ps * exp(-latvap / rh2o * ((one/Tsurf)-(one/T0)))

  u(1) = u(1) / (one + Cd * wind * dt / zam)
  v(1) = v(1) / (one + Cd * wind * dt / zam)
  qsv(1) = (qsv(1) + C * wind * qsats * dt / zat) / (one + C * wind * dt / zat)
  t(1) = (t(1) + C * wind * Tsurf * dt / zat) / (one + C * wind * dt / zat)

  qv(1) = qsv(1) / (one - qsv(1))

  !Estimate Temperature/Specific Humidity MOMENTUM
  !-----------------------------------------------
    t_m(1)= wp(1)*  t(1) + wm(1)*  t(2)
  qsv_m(1)= wp(1)*qsv(1) + wm(1)*qsv(2)

  !------------------------------------------------
  ! Boundary layer (AS IN SIMPLIFIED_PHYSICS_V6)
  !------------------------------------------------

  do k=1,nz

    CAE = 0.d0
    CCE = 0.d0
    CAm = 0.d0
    CCm = 0.d0

    !A,C for T,Q ! So Ke Momentum
    !----------------------------
    r_dp_t = pm(k-1) - pm(k)
    r_dp_t = 1.0d0/r_dp_t

    if (k/=1) then
       tv_m = t_m(k-1)*(1.0d0+zvir*qsv_m(k-1))
       rho_m= (pm(k-1)/(rair*tv_m))
       CAE  = r_dp_t*dt*gravit*gravit*Ke(k-1)*rho_m*rho_m/(pt(k-1)-pt(k))
    end if

    if (k/=nz) then
       tv_m = t_m(k)*(1.0d0+zvir*qsv_m(k))
       rho_m= (pm(k)/(rair*tv_m))
       CCE  = r_dp_t*dt*gravit*gravit*Ke(k)*rho_m*rho_m/(pt(k)-pt(k+1))
    end if

    !A,C for U,V ! So Km Thermo
    !--------------------------
    r_dp_m = pt(k) - pt(k+1)
    r_dp_m = 1.0d0/r_dp_m

    if (k/=1) then
       rho_t= (pt(k)/(rair*t(k)*(1.0d0+zvir*qsv(k))))
       CAm  = r_dp_m*dt*gravit*gravit*Km(k)*rho_t*rho_t/(pm(k-1)-pm(k))
    end if

    if (k/=nz) then
       rho_t= (pt(k+1)/(rair*t(k+1)*(1.0d0+zvir*qsv(k+1))))
       CCm  = r_dp_m*dt*gravit*gravit*Km(k+1)*rho_t*rho_t/(pm(k)-pm(k+1))
    end if

    if (k==1) then
       CEE(k) = CCE/(1.0d0+CAE+CCE)
       CEm(k) = CCm/(1.0d0+CAm+CCm)
       CFu(k) = u(k)/(1.0d0+CAm+CCm)
       CFv(k) = v(k)/(1.0d0+CAm+CCm)
       CFt(k) = ((p0/pt(k))**(rair/cpair)*t(k))/(1.0d0+CAE+CCE)
       CFq(k) = qsv(k)/(1.0d0+CAE+CCE)
    else
       CEE(k) = CCE/(1.0d0+CAE+CCE-CAE*CEE(k-1))
       CEm(k) = CCm/(1.0d0+CAm+CCm-CAm*CEm(k-1))
       CFu(k) = (u(k)+CAm*CFu(k-1))/(1.0d0+CAm+CCm-CAm*CEm(k-1))
       CFv(k) = (v(k)+CAm*CFv(k-1))/(1.0d0+CAm+CCm-CAm*CEm(k-1))
       CFt(k) = ((p0/pt(k))**(rair/cpair)*t(k)+CAE*CFt(k-1))/(1.0d0+CAE+CCE-CAE*CEE(k-1))
       CFq(k) = (qsv(k)+CAE*CFq(k-1))/(1.0d0+CAE+CCE-CAE*CEE(k-1))
    end if

  end do

  !---------------------------------------------------------------------------
  !Calculate the updated temperature and specific humidity and wind tendencies
  !---------------------------------------------------------------------------

  !First we need to calculate the tendencies at the top model level
  !----------------------------------------------------------------
  u(nz)   = CFu(nz)
  v(nz)   = CFv(nz)
  t(nz)   = CFt(nz)*(pt(nz)/p0)**(rair/cpair)
  qsv(nz) = CFq(nz)

  do k=nz-1,1,-1
     u(k)   =  CEm(k)*u(k+1)+CFu(k)
     v(k)   =  CEm(k)*v(k+1)+CFv(k)
     t(k)   = (CEE(k)*t(k+1)*(p0/pt(k+1))**(rair/cpair)+CFt(k))*(pt(k)/p0)**(rair/cpair)
     qsv(k) =  CEE(k)*qsv(k+1)+CFq(k)
  end do

  !Convert qsv to qv
  !-----------------
  do k=1,nz
    qv(k) = qsv(k) / (one - qsv(k))
  end do

  !Conversion from Mixing ratio to Specific for qc/qr
  !--------------------------------------------------
  do k=1,nz
     qsc(k)= qc(k)/(1.0d0 + qv(k))
     qsr(k)= qr(k)/(1.0d0 + qv(k))
  end do

  return

END SUBROUTINE DCMIP2016_PHYSICS

END MODULE DCMIP_2016_physics_module
