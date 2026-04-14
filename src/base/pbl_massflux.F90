module pbl_massflux
  implicit none
  private

  ! Internal parameters
  real, parameter :: NO_DRY_THERMAL=-1.         !No thermal for dry mass-flux scheme

  ! External API
  public :: massflux                            !Compute mass-flux in thermals      

contains

  subroutine massflux(F_pvars, F_kount, F_ni, F_nkm1)
    use phymem, only: phyvar
    use phybusidx, except=>buoy
    use tdpack, only: DELTA, CAPPA, GRAV, RGASD, CPD, CHLC
    use pbl_utils, only: WSTAR_MIN
    use vintphy, only: vint_thermo2mom
    implicit none

    ! Compute mass-flux profiles for coherent turbulence based on de Rooy et al. (2022).  The
    ! fluxes are computed and output on thermodynamic levels.

    ! Arguments
    integer, intent(in) :: F_ni                         !horizontal dimension
    integer, intent(in) :: F_nkm1                       !vertical dimension
    integer, intent(in) :: F_kount                      !time step number
    type(phyvar), pointer, contiguous :: F_pvars(:)     !physics buses

    ! External definitions
#include "phymkptr.hf"
#include "nocld.cdk"
#include "surface.cdk"
    
    ! Local parameters
    integer, parameter :: NITER=2                       !number of interations for updraft
    real, parameter :: W_INIT=0.1                       !initial updraft speed [m/s]
    real, parameter :: AU_INIT=0.1                      !dry updraft area fraction
    real, parameter :: DGAUSS=1.754                     !percentile for AU_INIT=0.1 (Neggers09 Table A1)
    real, parameter :: ACOEF=10./7.                     !"a" coefficient (Table 2 of deRooy22)
    real, parameter :: BCOEF=5./7.                      !"b" coefficient (Table 2 of deRooy22)
    real, parameter :: TAUE=400.                        !test entrainment timescale [s] (Neggers09 Eq 19)
    real, parameter :: CDRY=0.4                         !entrainment coefficient (deRooy22 Eq 9)
    real, parameter :: A1=40.                           !surface entrainment limiter [m] (deRooy22 App B)
    real, parameter :: A2=1.                            !inversion entrainment limiter [m] (deRooy22 App B)
    real, parameter :: C_WCASC=0.5                      !scaling for TKE source term (deRooy22 Eq 18) 
    logical, parameter :: W_LIN_AVG=.true.              !use linear vertical averaging in w instead of w^2 
    
    ! Local variables
    integer :: i, k, it
    real :: inroot, bdzt, buoy, entr, flq, flt, edz, zlevu, zentr, thpu, &
         hupu, dz, wpu, the, hue, thve, bterm, cterm, imp
    real, dimension(:), pointer :: zwstar, zfc_ag, zfv_ag, zpmoins, zzid
    real, dimension(F_ni) :: zlcl
    real, dimension(F_ni,F_nkm1) :: dzt, th, thp, hup, wp, exner, &
         prest, rho, thm, hum
    real, dimension(:,:), pointer :: ztmoins, zhumoins, zumoins, zvmoins, &
         zgzmom, zgztherm, zsigt, zqcmoins, zfxp, zwtng, zwqng, zwcascd
    real, pointer, dimension(:,:,:), contiguous :: zvcoef
    logical, dimension(F_ni) :: cloud

    ! Extract needed data from the buses
    MKPTR1D(zpmoins, pmoins, F_pvars)
    MKPTR1D(zwstar, wstar, F_pvars)
    MKPTR1D(zzid, zid, F_pvars)

    MKPTR1DK(zfc_ag, fc, indx_agrege, F_pvars)
    MKPTR1DK(zfv_ag, fv, indx_agrege, F_pvars)
    
    MKPTR2D(zgzmom, gzmom, F_pvars)
    MKPTR2D(zgztherm, gztherm, F_pvars)
    MKPTR2D(zfxp, fxp, F_pvars)
    MKPTR2D(zhumoins, humoins, F_pvars)
    MKPTR2D(zqcmoins, qcmoins, F_pvars)
    MKPTR2D(zsigt, sigt, F_pvars)
    MKPTR2D(ztmoins, tmoins, F_pvars)
    MKPTR2D(zumoins, umoins, F_pvars)
    MKPTR2D(zvmoins, vmoins, F_pvars)
    MKPTR2D(zwcascd, wcascd, F_pvars)
    MKPTR2D(zwqng, wqng, F_pvars)
    MKPTR2D(zwtng, wtng, F_pvars)

    MKPTR3D(zvcoef, vcoef, F_pvars)

    ! Basic Initializations
    do k=2,F_nkm1
       do i=1,F_ni
          dzt(i,k) = zgztherm(i,k-1) - zgztherm(i,k)
          prest(i,k) = zpmoins(i) * zsigt(i,k)
          exner(i,k) = zsigt(i,k)**CAPPA
          th(i,k) = ztmoins(i,k) / exner(i,k)
          rho(i,k) = prest(i,k) / (RGASD * ztmoins(i,k))
       enddo
    enddo
    call vint_thermo2mom(thm, th, zvcoef, F_ni, F_nkm1)
    call vint_thermo2mom(hum, zhumoins, zvcoef, F_ni, F_nkm1)

    ! Initialize parcel properties (no k-dependence as per deRooy22 p1518)
    do i=1,F_ni
       if (zwstar(i) > 0) then
          flq = zfv_ag(i) / (rho(i,F_nkm1) * CHLC)
          flt = zfc_ag(i) / (rho(i,F_nkm1) * CPD)
          wp(i,F_nkm1) = W_INIT
       else
          flq = 0.
          flt = 0.
          wp(i,F_nkm1) = 0.
       endif
       hup(i,F_nkm1) = zhumoins(i,F_nkm1) + DGAUSS * max(flq, 0.) / max(zwstar(i), WSTAR_MIN)
       thp(i,F_nkm1) = th(i,F_nkm1) + DGAUSS * max(flt, 0.) / max(zwstar(i), WSTAR_MIN)       
       do k=1,F_nkm1-1
          wp(i,k) = 0.
       enddo
    enddo

    ! Dry test parcel ascent for first estimate of inversion height (zzid)
    zzid(:) = NO_DRY_THERMAL
    cloud(:) = .false.
    do i=1,F_ni
       k = F_nkm1
       TEST_PARCEL: do while (wp(i,k) > 0 .and. k > 2 .and. .not.cloud(i))
          
          ! Confirm no clouds in the environment
          if (associated(zqcmoins)) then
             if (zqcmoins(i,k) > epsilon(zqcmoins)) cloud(i) = .true.
          endif
          if (associated(zfxp)) then
             if (zfxp(i,k) > CLDFTH) cloud(i) = .true.
          endif
          if (cloud(i)) cycle
          
          ! Compute updraft speed using buoyancy at the current departure level
          thve = th(i,k) * (1. + DELTA * zhumoins(i,k))
          buoy = GRAV / thve * (thp(i,k) * (1. + DELTA*hup(i,k)) - thve)
          bdzt = BCOEF * dzt(i,k) / TAUE
          inroot = bdzt**2 + 4. * (2. * dzt(i,k) * ACOEF * buoy + wp(i,k) * (wp(i,k) - bdzt))
          if (inroot > 0) then
             wp(i,k-1) = 0.5 * (-bdzt + sqrt(inroot))
          else             
             wp(i,k-1) = 0.
          endif

          ! Compute parcel thermodynamic properties
          entr = 1. / (TAUE * max(0.5*(wp(i,k) + wp(i,k-1)), epsilon(wp)))
          edz = entr * dzt(i,k)
          thp(i,k-1) = (thp(i,k) * (1. - 0.5*edz) + edz * thm(i,k)) / (1. + 0.5*edz)
          hup(i,k-1) = (hup(i,k) * (1. - 0.5*edz) + edz * hum(i,k)) / (1. + 0.5*edz)

          ! Estimate LCL (cloud forms) or inversion level (dry)
          if (wp(i,k-1) > 0) then
             cloud(i) = check_lcl(zzid(i), thp(i,k-1), thp(i,k), hup(i,k-1), hup(i,k), &
                  prest(i,k-1), prest(i,k), exner(i,k-1), exner(i,k), zgztherm(i,k-1), &
                  zgztherm(i,k))
          elseif (k < F_nkm1) then             
             zzid(i) = zgztherm(i,k) + wp(i,k)**2 / ((BCOEF / TAUE) * wp(i,k) - 2. * ACOEF * buoy)
          endif
          
          ! Move to the next level
          k = k-1
       enddo TEST_PARCEL
    enddo
    
    ! Iterate on dry ascent for updraft properties
    DRY_UPDRAFT_ITER: do it=1,NITER
       cloud(:) = .false.
       do i=1,F_ni
          if (nint(zzid(i)) == NO_DRY_THERMAL) cycle
          k = F_nkm1
          DRY_UPDRAFT: do while (wp(i,k) > 0 .and. k > 2 .and. .not.cloud(i))
             
             ! Compute parcel thermodynamic properties
             zlevu = min(zgztherm(i,k-1), zzid(i))
             !dz = (zlevu - zgztherm(i,k))
             dz = max( 1.e-3, (zlevu - zgztherm(i,k)) )
             zentr = 0.5 * (zlevu + zgztherm(i,k))
             !entr = CDRY * (1. / (zentr + A1) + 1. / (zzid(i) - zentr + A2))
             entr = CDRY * (1. / (zentr + A1) + 1. / max( 1.e-3, (zzid(i) - zentr + A2) ))
             edz = entr * dz
             the = th(i,k) + (th(i,k-1) - th(i,k)) &
                  / (zgztherm(i,k-1) - zgztherm(i,k)) * dz
             hue = zhumoins(i,k) + (zhumoins(i,k-1) - zhumoins(i,k)) &
                  / (zgztherm(i,k-1) - zgztherm(i,k)) * dz
             thpu = (thp(i,k) * (1. - 0.5*edz) + edz * the) / (1. + 0.5*edz)
             hupu = (hup(i,k) * (1. - 0.5*edz) + edz * hue) / (1. + 0.5*edz)
 
             ! Compute updraft speed using layer-average buoyancy
             thve =  the * (1. + DELTA * hue)
             buoy = GRAV / thve * ( &
                  0.5 * (thp(i,k) * (1. + DELTA*hup(i,k)) + thpu * (1. + DELTA*hupu)) - thve )
             if (W_LIN_AVG) then
                ! Vertical average in w
                imp = 1. / (2.*dz) + BCOEF * entr / 4.
                bterm = (BCOEF * entr * wp(i,k) / 2.) / imp
                cterm = -((1. / (2.*dz) - BCOEF * entr / 4.) * wp(i,k)**2 + ACOEF * buoy) / imp
                inroot = bterm**2 - 4. * cterm
                if (inroot > 0) then
                   wpu = 0.5 * (-bterm + sqrt(inroot))
                else
                   wpu = 0.
                endif
             else
                ! Vertical average in w^2
                inroot = (ACOEF * buoy + (1./(2.*dz) - BCOEF * entr) * wp(i,k)**2) / &
                     (1./(2.*dz) + BCOEF * entr) 
                if (inroot > 0) then
                   wpu = sqrt(inroot)
                else
                   wpu = 0.
                endif
             endif
             
             ! Confirm no condensation in the updraft or revise inversion level
             if (wpu > 0) then
                thp(i,k-1) = thpu
                hup(i,k-1) = hupu
                wp(i,k-1) = wpu
                cloud(i) = check_lcl(zzid(i), thp(i,k-1), thp(i,k), hup(i,k-1), hup(i,k), &
                     prest(i,k-1), prest(i,k), exner(i,k-1), exner(i,k), zgztherm(i,k-1), &
                     zgztherm(i,k))                
                if (cloud(i) .or. zzid(i) < zgztherm(i,k-1)) wp(i,1:k-1) = 0.
             elseif (k < F_nkm1) then
                if (W_LIN_AVG) then
                   zzid(i) = zgztherm(i,k) + wp(i,k)**2 / &
                        (BCOEF * entr * wp(i,k)**2 / 2. - 2. * ACOEF * buoy)
                else
                   zzid(i) = zgztherm(i,k) + wp(i,k)**2 / &
                        (BCOEF * entr * wp(i,k)**2 - 2. * ACOEF * buoy)
                endif
                wp(i,1:k-1) = 0.
             endif
             
             ! Move to the next level
             k = k-1
          enddo DRY_UPDRAFT
       enddo
       
    enddo DRY_UPDRAFT_ITER

    ! Compute dry convective temperature and moisture fluxes (deRooy22 Eq. 3)
    zwtng(:,:) = 0.
    zwqng(:,:) = 0.
    do i=1,F_ni
       if (nint(zzid(i)) == NO_DRY_THERMAL) cycle
       k = F_nkm1
       do while (wp(i,k) > 0. .and. k > 1)
          zwtng(i,k) = AU_INIT * wp(i,k) * (thp(i,k) - th(i,k))
          zwqng(i,k) = AU_INIT * wp(i,k) * (hup(i,k) - zhumoins(i,k))
          k = k-1
       enddo
    enddo

    ! Compute dry energy cascade TKE source term (deRooy22 Eq. 18)
    zwcascd(:,:) = 0.
    do i=1,F_ni
       if (nint(zzid(i)) == NO_DRY_THERMAL) cycle
       k = F_nkm1
       do while (zgztherm(i,k) <= zzid(i) .and. k > 1)
          entr = CDRY * (1. / (zgztherm(i,k) + A1) + 1. &
               / (zzid(i) - zgztherm(i,k) + A2))
          zwcascd(i,k) = C_WCASC * entr * rho(i,k) * AU_INIT * wp(i,k)**3
          k = k-1
       enddo
    enddo
    
    ! End of subprogram
    return
  end subroutine massflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical function check_lcl(zi, thpu, thp, hupu, hup, prestu, prest, &
       exneru, exner, gzthermu, gztherm) result(cloud)
    use tdpack, only: foqst
    implicit none
    
    ! Check for condensation and return LCL height if present

    ! Arguments
    real, intent(in) :: thpu            !Parcel potential temperature at k-1 (K)
    real, intent(in) :: thp             !Parcel potential temperature at k (K)
    real, intent(in) :: hupu            !Parcel mixing ratio at k-1 (kg/kg)
    real, intent(in) :: hup             !Parcel mixing ratio at k (kg/kg)
    real, intent(in) :: prestu          !Thermo-level pressure at k-1 (Pa)
    real, intent(in) :: prest           !Thermo-level pressure at k (Pa)
    real, intent(in) :: exneru          !Exner function at k-1
    real, intent(in) :: exner           !Exner function at k
    real, intent(in) :: gzthermu        !Thermo-level height at k-1 (m)
    real, intent(in) :: gztherm         !Thermo-level height at k (m)
    real, intent(inout) :: zi           !Height of LCL if cloud forms (m)

    ! Internal variables
    real :: tp, qsatp, rhp, rhpu, lcl

    ! Check for condensation and compute LCL
    tp = thpu * exneru
    qsatp = foqst(tp, prestu)    
    if (hupu >= qsatp) then
       cloud = .true.
       rhp = hup / foqst(thp*exner, prest)
       rhpu = hupu / qsatp
       lcl = gztherm + (1. - rhp) * (gzthermu - gztherm) / (rhpu - rhp)
       if (nint(zi) == NO_DRY_THERMAL) then
          zi = lcl
       else
          zi = min(lcl, zi)
       endif
    else
       cloud = .false.
    endif

    ! End of subprogram
    return
  end function check_lcl
  
end module pbl_massflux
