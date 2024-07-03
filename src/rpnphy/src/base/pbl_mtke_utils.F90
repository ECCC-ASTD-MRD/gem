module pbl_mtke_utils
  use vintphy, only: vint_thermo2mom2, vint_thermo2mom3
  implicit none
  private

  ! API functions
  public :: baktotq                                     !Convert conserved variable tendencies to state tendencies
  public :: blcloud                                     !Buoyancy flux calculation
  public :: thermco                                     !Compute thermodynamic coefficients
  public :: tkealg                                      !Algebraic solution of the TKE equation

  ! Internal functions
  private :: clsgs                                      !Diagnose PBL clouds

  ! Internal parameters
  real, parameter :: CLEFC4 = 4.5                       !Coefficient from higher-order closure
  real, parameter :: CLEFC8 = 6.5                       !Coefficient from higher-order closure
  
  ! API parameters
  real, parameter, public :: BLCONST_CK = 0.516         !Coefficient for eddy diffusivity calculation
  real, parameter, public :: BLCONST_CE = 0.14          !Coefficient for dissipation term
  real, parameter, public :: BLCONST_CU = 3.75          !Coefficient for surface horizontal TKE source
  real, parameter, public :: BLCONST_CW = 0.2           !Coefficient for surface vertical TKE source
  real, parameter, public :: CLEFAE = 3.*CLEFC4/CLEFC8  !Coefficient for diffusion of TKE
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine baktotq(dt,dqv,dqc,thl,qw,dthl,dqw,qc,s,sw,ps,gztherm,tif,fice,&
       tve,hpbl,hflux,qcbl,fnn,fn,fngauss,fnnonloc,c1,zn,ze,mg,mrk2, &
       vcoef,pblsigs,pblq1,tau,n,nk)
    use tdpack_const, only: CPD, CAPPA, CHLC, CHLF
    use ens_perturb, only: ens_nc2d
    implicit none
!!!#include <arch_specific.hf>

    !@Arguments
    integer, intent(in) :: n                             !horizontal dimension
    integer, intent(in) :: nk                            !vertical dimension
    real, intent(in) :: tau                              !time step length (s)
    real, dimension(n,nk), intent(in) :: thl             !liquid water potential temperature (K; theta_l)
    real, dimension(n,nk), intent(in) :: qw              !total water mixing ratio (kg/kg; q_tot)
    real, dimension(n,nk), intent(in) :: qc              !PBL cloud water content (kg/kg)
    real, dimension(n,nk), intent(in) :: s               !sigma for full levels
    real, dimension(n,nk), intent(in) :: sw              !sigma for working levels
    real, dimension(n),    intent(in) :: ps              !surface pressure (Pa)
    real, dimension(n,nk), intent(in) :: gztherm         !height of thermodynamic levels (m)
    real, dimension(n,nk), intent(in) :: tif             !temperature used for ice fraction (K)
    real, dimension(n,nk), intent(in) :: fice            !ice fraction
    real, dimension(n,nk), intent(in) :: dthl            !tendency of theta_l (K/s)
    real, dimension(n,nk), intent(in) :: dqw             !tendency of q_tot (kg/kg/s)
    real, dimension(n,nk), intent(in) :: tve             !virtual temperature on e-levs (K)
    real, dimension(n), intent(in) :: hpbl               !boundary layer height (m)
    real, dimension(n), intent(in) :: hflux              !surface heat flux (W/m2)
    real, dimension(n,nk), intent(in) :: c1              !constant C1 in second-order moment closure
    real, dimension(n,nk), intent(in) :: zn              !mixing length (m)
    real, dimension(n,nk), intent(in) :: ze              !dissipation length (m)
    real, dimension(n),    intent(in) :: mg              !land-sea mask
    real, dimension(n,ens_nc2d), intent(in) :: mrk2      !Markov chains for stochastic parameters
    real, dimension(*), intent(in) :: vcoef              !coefficients for vertical interpolation   
    real, dimension(n,nk), intent(inout) :: fnn          !flux enhancement * cloud fraction
    real, dimension(n,nk), intent(inout) :: fnnonloc     !nonlocal cloud fraction
    real, dimension(n,nk), intent(out) :: fn             !cloud fraction
    real, dimension(n,nk), intent(out) :: qcbl           !water content of PBL clouds (kg/kg)
    real, dimension(n,nk), intent(out) :: dt             !tendency of dry air temperature (K/s)
    real, dimension(n,nk), intent(out) :: dqv            !tendency of specific humidity (kg/kg/s)
    real, dimension(n,nk), intent(out) :: dqc            !tendency of cloud water content (kg/kg/s)
    real, dimension(n,nk), intent(out) :: fngauss        !Gaussian cloud fraction
    real, dimension(n,nk), intent(out) :: pblsigs        !Subgrid moisture variance
    real, dimension(n,nk), intent(out) :: pblq1          !Normalized saturation deficit

    !@Author  J. Mailhot (Nov 2000)
    !@Revision
    ! 001      A.-M. Leduc (Oct 2001) Automatic arrays
    ! 002      B. Bilodeau and J. Mailhot (Dec 2001) Add a test to
    !                      check the presence of advected explicit cloud water.
    ! 003      J. Mailhot (Nov 2000) Cleanup of routine
    ! 004      J. Mailhot (Feb 2003) - MOISTKE option based on implicit clouds only
    ! 005      A-M. Leduc (Jun 2003) - pass ps to clsgs---> clsgs2
    ! 006      J. P. Toviessi ( Oct. 2003) - IBM conversion
    !               - calls to exponen4 (to calculate power function '**')
    !               - etc.
    ! 007      B. Bilodeau (Dec 2003)   More optimizations for IBM
    !                                   - Call to vspown1
    !                                   - Replace divisions by multiplications
    ! 008      L. Spacek (Dec 2007) - add "vertical staggering" option
    !                                 change the name to baktotq3
    ! 009      A. Zadra (Oct 2015) -- add land-water mask (MG) to input, which
    !                                 is then passed on to CLSGS4
    ! 010      A. Zadra, R. McT-C (Sep 2016) - implement non-local scaling
    !                      deveoped by J. Mailhot/A. Lock (Aug 2012)
    !@Object
    !          Transform conservative variables and their tendencies
    !          back to non-conservative variables and tendencies.
    !          Calculate the boundary layer cloud properties (cloud fraction, cloud
    !          water content, flux enhancement factor).

    ! Local parameter definitions
    integer, parameter :: IMPLICIT_CLOUD=0
    logical, parameter :: COMPUTE_FROM_STATE=.false.

    ! Local variables
    real :: cpdinv,tauinv
    real, dimension(n,nk) :: exner,thl_star,qw_star,acoef,bcoef,ccoef,alpha,beta,qcp,unused,qv

    ! Precompute inverses
    CPDINV = 1./CPD
    TAUINV = 1./TAU

    ! Update conserved variables following diffusion
    exner = sw**CAPPA
    thl_star = thl + tau*dthl
    qw_star = qw + tau*dqw

    ! Compute thermodynamic coefficients from conserved variables
    call thermco(unused,unused,unused,sw,ps,tif,fice,fnn,thl_star,qw_star,acoef,bcoef,ccoef,alpha,beta,&
         IMPLICIT_CLOUD,COMPUTE_FROM_STATE,n,nk)

    ! Retrive updated cloud water content (qcp) from conserved variables
    call clsgs(thl_star,tve,qw_star,qcp,fn,fnn,fngauss,fnnonloc,c1,zn,ze,hpbl,hflux,s,ps,gztherm,&
         mg,mrk2,acoef,bcoef,ccoef,vcoef,pblsigs,pblq1,n,nk)

    ! Convert back to state variables and tendencies
    qv = qw - max(0.,qc)
    dqc = (max(0.,qcp) - max(0.,qc)) * tauinv
    dqc = max(dqc,-max(0.0,qc)*tauinv) !prevent negative values of qc
    dqc = min(dqc,dqw + max(0.0,qv)*tauinv) !prevent over-depletion of water vapour by PBL cloud
    qcbl = max(0.,qc) + dqc*tau
    dt = exner*dthl + ((CHLC+fice*CHLF)*CPDINV)*dqc
    dqv = max((dqw-dqc),-max(0.0,qv)*tauinv)

    return
  end subroutine baktotq
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine blcloud(u,v,t,tve,qv,qc,fnn,frac,fngauss,fnnonloc,w_cld,&
       wb_ng,wthl_ng,wqw_ng,uw_ng,vw_ng,f_cs,dudz,dvdz,&
       hpar,frv,z0m,fb_surf,gzmom,ze,s,sw,ps,dudz2,ri,&
       dthv,tau,vcoef,n,nk,ncld)
    use tdpack_const, only: DELTA, GRAV, RGASD
    use phy_options
    use pbl_utils, only: ficemxp, dvrtdf

    implicit none
!!!#include <arch_specific.hf>
    !
    ! Arguments
    integer, intent(in) :: n                             !horizontal dimension
    integer, intent(in) :: nk                            !vertical dimension
    integer, intent(in) :: ncld                          !number of cloud profiles for length scales
    real, intent(in) :: tau                              !time step length (s)
    real, dimension(n,nk), intent(in) :: u               !u-component wind (m/s)
    real, dimension(n,nk), intent(in) :: v               !v-component wind (m/s)
    real, dimension(n,nk), intent(in) :: t               !dry air temperature (K)
    real, dimension(n,nk), intent(in) :: tve             !virt. temperature on e-lev (K)
    real, dimension(n,nk), intent(in) :: qv              !specific humidity (kg/kg)
    real, dimension(n,nk), intent(in) :: qc              !PBL cloud water content (kg/kg)
    real, dimension(n,nk), intent(in) :: fngauss         !Gaussian (local) cloud fraction
    real, dimension(n), intent(in) :: frv                !friction velocity (m/s)
    real, dimension(n), intent(in) :: z0m                !roughness length (m)
    real, dimension(n), intent(in) :: fb_surf            !surface buoyancy flux
    real, dimension(n,nk), intent(in) :: gzmom           !height of momentum levels (m)
    real, dimension(n,nk), intent(in) :: ze              !height of e-levs (m)
    real, dimension(n,nk), intent(in) :: s               !sigma for full levels
    real, dimension(n,nk), intent(in) :: sw              !sigma for working levels
    real, dimension(n), intent(in) :: ps                 !surface pressure (Pa)
    real, dimension(*), intent(in) :: vcoef              !coefficients for vertical interpolation
    real, dimension(n), intent(inout) :: hpar            !height of parcel ascent (m)
    real, dimension(n,nk), intent(inout) :: fnn          !flux enhancement * cloud fraction
    real, dimension(n,nk), intent(out) :: frac           !cloud fraction (nonlocal)
    real, dimension(n,nk,ncld), intent(out) :: w_cld     !cloud-layer velocity scales
    real, dimension(n,nk), intent(out) :: wb_ng          !non-gradient buoyancy flux (nonlocal)
    real, dimension(n,nk), intent(out) :: wthl_ng        !non-gradient theta_l flux (nonlocal)
    real, dimension(n,nk), intent(out) :: wqw_ng         !non-gradient total water flux (nonlocal)
    real, dimension(n,nk), intent(out) :: uw_ng          !non-gradient x-component stress (nonlocal)
    real, dimension(n,nk), intent(out) :: vw_ng          !non-gradient y-component stress (nonlocal)
    real, dimension(n,nk), intent(out) :: dudz           !x-compontent vertical wind shear (s^-1)
    real, dimension(n,nk), intent(out) :: dvdz           !y-compontent vertical wind shear (s^-1)
    real, dimension(n,nk), intent(out) :: dudz2          !square of vertical wind shear (s^-2)
    real, dimension(n,nk), intent(out) :: ri             !gradient Richardson number
    real, dimension(n,nk), intent(out) :: dthv           !buoyancy flux (m2/s2)
    real, dimension(n,nk), intent(out) :: fnnonloc       !nonlocal cloud fraction
    real, dimension(n,nk), intent(out) :: f_cs           !factor for l_s coeff c_s (nonlocal)

    !Author
    !          J. Mailhot (Nov 2000)

    !Revision
    ! 001      A.-M. Leduc (Oct 2001) Automatic arrays
    ! 002      J. Mailhot (Jun 2002) Change calling sequence and rename BLCLOUD1
    ! 003      J. Mailhot (Feb 2003) Change calling sequence and rename BLCLOUD2
    ! 004      L. Spacek (Dec 2007) - add "vertical staggering" option
    ! 005                             change the name to blcloud3
    ! 006      J. Mailhot/A. Lock (Aug 2012) Revisions for non-local scaling option
    !                      Change calling sequence and rename BLCLOUD4

    !Object
    !          Calculate the boundary layer buoyancy parameters (virtual potential
    !          temperature, buoyancy flux) and the vertical shear squared.

    !Notes
    !          Implicit (i.e. subgrid-scale) cloudiness scheme for unified
    !             description of stratiform and shallow, nonprecipitating
    !             cumulus convection appropriate for a low-order turbulence
    !             model based on Bechtold et al.:
    !            - Bechtold and Siebesma 1998, JAS 55, 888-895
    !            - Cuijpers and Bechtold 1995, JAS 52, 2486-2490
    !            - Bechtold et al. 1995, JAS 52, 455-463
    !            - Bechtold et al. 1992, JAS 49, 1723-1744
    !          The boundary layer cloud properties (cloud fraction, cloud water
    !            content) are computed in the companion S/R CLSGS.

    ! Local parameter definitions
    integer, parameter :: IMPLICIT_CLOUD=0
    logical, parameter :: COMPUTE_FROM_CONSERVED=.true.

    ! Local variable declarations
    integer :: j, k
    real, dimension(n,nk) :: thl,qw,alpha,beta,a,b,c,dz,dqwdz,dthldz,coefthl,coefqw,ficelocal, &
         unused,thv,tmom,qvmom,qcmom

    ! Compute ice fraction from temperature profile
    call ficemxp(ficelocal,unused,unused,t,n,nk)

    ! Compute thermodynamic coefficients following Bechtold and Siebsma (JAS 1998)
    call thermco(t,qv,qc,sw,ps,t,ficelocal,fnn,thl,qw,a,b,c,alpha,beta, &
         IMPLICIT_CLOUD,COMPUTE_FROM_CONSERVED,n,nk)

    ! Compute layer thicknesses
    do k=1,nk-1
       dz(:,k) = -RGASD*tve(:,k)*alog(s(:,k+1)/s(:,k))/GRAV
    end do
    dz(:,nk) = 0.

    ! Compute non-gradient fluxes and PBL cloud fraction in nonlocal formulation
    NONLOCAL_CLOUD: IF ( pbl_nonloc == 'LOCK06' ) THEN
       ! Convert to non-staggered state for nonlocal cloud estimates (FIXME)
       call vint_thermo2mom3(tmom, qvmom, qcmom, &
            &                t,    qv,    qc,    vcoef,n,nk)
       ! Compute non-local cloud properties
       call nlcalc(fnn,frac,fngauss,fnnonloc,wb_ng,wthl_ng,wqw_ng,w_cld,u,v,uw_ng,vw_ng, &
            f_cs,tmom,qvmom,qcmom,frv,z0m,fb_surf,hpar,tau,gzmom,ze,ps,s,n,n,nk,size(w_cld,dim=3))
    else
       w_cld = 0.
       f_cs = 1.
    endif NONLOCAL_CLOUD

    ! Compute terms of buoyancy flux equation (Eq. 4 of Bechtold and Siebsma (JAS 1998))
    thv = thl + alpha*qw + beta*qc
    coefthl = 1.0 + DELTA*qw  - beta*b*fnn
    coefqw = alpha + beta*a*fnn

    ! Add to WB_NG (non-gradient flux of buoyancy) the contributions from
    ! WTHL_NG and WQW_NG (non-gradient fluxes of theta_l and q_w)
    NONLOCAL_FLUX: IF (pbl_nonloc == 'LOCK06') THEN
       do k=1,nk-1
          do j=1,n
             wb_ng(j,k) = wb_ng(j,k) + ( coefthl(j,k)*wthl_ng(j,k) &
                  + coefqw(j,k)*wqw_ng(j,k) ) &
                  * ( grav / thv(j,k) )
          end do
       end do
       wb_ng(:,nk) = 0.
    endif NONLOCAL_FLUX

    ! Compute vertical derivatives of conserved variables
    call vint_thermo2mom2(thl, qw, &
         &                thl, qw, vcoef, n, nk)
    call dvrtdf(dthldz,thl,dz,n,n,n,nk)
    call dvrtdf(dqwdz,qw,dz,n,n,n,nk)

    ! Compute buoyancy flux
    do k=1,nk-1
       do j=1,n
          dthv(j,k) = (coefthl(j,k)*dthldz(j,k) + coefqw(j,k)*dqwdz(j,k)) * (GRAV/thv(j,k))
       enddo
    enddo
    dthv(:,nk) = 0.

    ! Compute vertical wind shear
    call dvrtdf(dudz,u,dz,n,n,n,nk)
    call dvrtdf(dvdz,v,dz,n,n,n,nk)
    dudz2 = dudz**2 + dvdz**2

    ! Compute gradient Richardson number
    ri = dthv / (dudz2+1e-6)  !FIXME: contains large background shear
    ri(:,nk) = 0.

    return
  end subroutine blcloud
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine thermco(t,qv,qc,sw,ps,tif,fice,fnn,thl,qw,acoef,bcoef,ccoef,alpha,beta, &
       type,inmode,n,nk)
    use tdpack
    use pbl_utils, only: ficemxp
    implicit none
!!!#include <arch_specific.hf>

    ! Arguments
    integer, intent(in) :: n                             !horizontal dimension
    integer, intent(in) :: nk                            !vertical dimension
    integer, intent(in) :: type                          !cloud type switch (0=implicit; 1=explicit; 2=both)
    logical, intent(in) :: inmode                        !input mode (T=compute from t,qv,qc; F=compute from thl,qw)
    real, dimension(n,nk), intent(in) :: t               !dry air temperature (K)
    real, dimension(n,nk), intent(in) :: qv              !specific humidity (kg/kg)
    real, dimension(n,nk), intent(in) :: qc              !PBL cloud water content (kg/kg)
    real, dimension(n,nk), intent(in) :: sw              !sigma for working levels
    real, dimension(n), intent(in) :: ps                 !surface pressure (Pa)
    real, dimension(n,nk), intent(in) :: tif             !temperature used for ice fraction (K)
    real, dimension(n,nk), intent(in) :: fice            !ice fraction
    real, dimension(n,nk), intent(in) :: fnn             !flux enhancement * cloud fraction
    real, dimension(n,nk), intent(inout) :: thl          !liquid water potential temperature (K; theta_l)
    real, dimension(n,nk), intent(inout) :: qw           !total water content (kg/kg; q_tot)
    real, dimension(n,nk), intent(out) :: acoef          !thermodynamic coefficient A
    real, dimension(n,nk), intent(out) :: bcoef          !thermodynamic coefficient B
    real, dimension(n,nk), intent(out) :: ccoef          !thermodynamic coefficient C
    real, dimension(n,nk), intent(out) :: alpha          !thermodynamic coefficient alpha
    real, dimension(n,nk), intent(out) :: beta           !thermodynamic coefficient beta

    !Author
    !          J. Mailhot (Nov 1999)

    !Revision
    ! 001      J. Mailhot  (Jan 2000) - Changes to add mixed-phase mode
    ! 002      A.-M. Leduc (Oct 2001) Automatic arrays
    ! 003      J. Mailhot  (Jun 2002) - Add cloud type and input mode
    !                       Change calling sequence and rename THERMCO2
    ! 004      A. Plante   (May 2003) - IBM conversion
    !                         - calls to exponen4 (to calculate power function '**')
    !                         - divisions replaced by reciprocals
    !                         - calls to optimized routine mfdlesmx
    ! 005      B. Bilodeau (Aug 2003) - exponen4 replaced by vspown1

    !Object
    !          Calculate the thermodynamic coefficients used in the presence of clouds
    !          and the conservative variables.

    !Notes
    !          See definitions in:
    !          - Bechtold and Siebesma 1998, JAS 55, 888-895

    ! Local variable declarations
    integer :: j,k
    real, dimension(n,nk) :: pres,exner,qsat,dqsat,th,tl,ffice,tfice,dfice,work

    ! Precompute pressure and Exner function
    do k=1,nk
       pres(:,k) = sw(:,k)*ps
       exner(:,k) = sw(:,k)**CAPPA
    enddo
    ffice = fice

    ! Compute conserved variables from state inputs
    COMPUTE_CONSERVED: if (inmode) then
       if ( type .eq. 0 ) call ficemxp(ffice,tfice,dfice,tif,n,nk)
       th = dble(t) / dble(exner)
       thl = th * (1.0-((CHLC+ffice*CHLF)/CPD) * (qc/t))
       qw = qv + qc
    endif COMPUTE_CONSERVED

    ! Compute liquid water temperature
    tl = exner*thl

    ! Compute saturation specific humidity for selected cloud type
    CLOUD_TYPE: select case (type)

    case(0) !Implicit clouds
       call ficemxp(work,tfice,dfice,tif,n,nk)
       do k=1,nk
          do j=1,n
             qsat(j,k) = fqsmx(tl(j,k),pres(j,k),tfice(j,k))
          enddo
       enddo
       call mfdlesmx(work,tl,tfice,dfice,n,nk)
       do k=1,nk
          do j=1,n
             dqsat(j,k) = fdqsmx(qsat(j,k),work(j,k) )
          enddo
       enddo

    case (1) !Explicit clouds
       do k=1,nk
          do j=1,n
             qsat(j,k) = foqsa(tl(j,k),pres(j,k))
             dqsat(j,k) = fodqa(qsat(j,k),tl(j,k))
          enddo
       enddo

    case (2) !Combined implicit and explicit clouds
       call ficemxp(work,tfice,dfice,tif,n,nk)
       call mfdlesmx(work,tl,tfice,dfice,n,nk)
       do k=1,nk
          do j=1,n
             if (fnn(j,k) < 1.0) then
                qsat(j,k) = fqsmx(tl(j,k),pres(j,k),tfice(j,k))
                dqsat(j,k) = fdqsmx(qsat(j,k),work(j,k) )
             else
                qsat(j,k) = foqsa(tl(j,k),pres(j,k))
                dqsat(j,k) = fodqa(qsat(j,k),tl(j,k))
             endif
          enddo
       enddo

    end select CLOUD_TYPE

    ! Compute thermodynamic coefficients following Bechtold and Siebesma (JAS 1998) Appendix A
    acoef = 1.0/( 1.0 + ((CHLC+ffice*CHLF)/CPD)*dqsat)
    bcoef = acoef*exner*dqsat
    ccoef = acoef*(qw-qsat)
    if (inmode) then
       alpha = DELTA*th
       beta = ((CHLC+ffice*CHLF)/CPD)/exner - (1.0+DELTA)*th
    else
       alpha = 0.
       beta = 0.
    endif

    return
  end subroutine thermco
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tkealg(estar,en,zn,ze,dvdz2,buoy_flux,diss,shr,buoy,tau,n,nk)
    use phy_options, only: etrmin2

    implicit none
!!!#include <arch_specific.hf>

    ! Argument declaration
    integer, intent(in) :: n                          !horizontal dimension
    integer, intent(in) :: nk                         !vertical dimension
    real, intent(in) :: tau                           !time step (s)
    real, dimension(n,nk), intent(in) :: en           !turbulent kinetic energy (m2/s2)
    real, dimension(n,nk), intent(in) :: zn           !mixing length scale (m)
    real, dimension(n,nk), intent(in) :: ze           !dissipation length scale (m)
    real, dimension(n,nk), intent(in) :: dvdz2        !square of vertical shear (/s2)
    real, dimension(n,nk), intent(in) :: buoy_flux    !buoyancy flux (/s2)
    real, dimension(n,nk), intent(out) :: estar       !updated turbulent kinetic energy (time *; m2/s2)
    real, dimension(n,nk), intent(out) :: diss        !viscous dissipation term (m2/s3)
    real, dimension(n,nk), intent(out) :: shr         !shear generation term (m2/s3)
    real, dimension(n,nk), intent(out) :: buoy        !buoyancy generation/suppression term (m2/s3)

    !Author
    !          J. Mailhot (August 1999)
    !Revision
    ! 001      A.-M. Leduc (Oct 2001) Automatic arrays
    ! 002      B. Bilodeau (Mar 2003) Include machcon.cdk
    ! 003      A. Plante   (May 2003) IBM conversion
    !             - calls to vsqrt routine (from massvp4 library)
    !             - calls to vstan routine (from massvp4 library)
    !             - calls to vslog routine (from massvp4 library)
    ! 004      a. Plante   (juillet 2003) correction of bugs in IBM conversion
    !             - optimization of tan removed
    !             - optimization of log removed
    !
    !Object
    !          Solve the algebraic part of the TKE equation

    ! Local variable declaration
    integer :: j,k
    real :: clamda,sqrt_tke,tanlim
    real, dimension(n,nk) :: b,c,b_over_c
    real(kind=8) :: r8tmp

    ! Compute coefficients (dble precision) and solve the algebraic TKE equation
    tanlim = exp(12. * log(2.))
    do k=1,nk-1
       do j=1,n
          b(j,k) = BLCONST_CK*zn(j,k)*(dvdz2(j,k) - buoy_flux(j,k))
          c(j,k) = BLCONST_CE/ze(j,k)
          r8tmp = (abs(b(j,k)/c(j,k)))
          b_over_c(j,k) = sqrt(r8tmp)
          r8tmp = en(j,k)
          sqrt_tke = sqrt(r8tmp)
          if( abs(b(j,k)) < epsilon(b) ) then
             estar(j,k) = sqrt_tke / (1.+0.5*sqrt_tke*c(j,k)*tau)
          elseif( b(j,k) > 0. ) then
             estar(j,k) = b_over_c(j,k)*( -1.0+2.0/(1.0+exp(-b_over_c(j,k)*c(j,k)*tau) * &
                  (-1.0+2.0/(1.0+sqrt_tke/b_over_c(j,k)))) )
          else
             estar(j,k)=b_over_c(j,k)*tan(min(tanlim,max(atan(sqrt_tke/b_over_c(j,k)) - &
                  0.5*b_over_c(j,k)*c(j,k)*tau , 0.0 ) ) )
          endif
          estar(j,k) = max(etrmin2, estar(j,k)**2)
       enddo
    enddo

    ! Compute TKE equation terms for output purposes (time series)
    do k=1,nk-1
       do j=1,n
          b_over_c(j,k) = b(j,k)/c(j,k)
          b(j,k)=0.0
          if( en(j,k) /= b_over_c(j,k) .and. estar(j,k) /= b_over_c(j,k) ) then
             b(j,k)=alog( abs( (estar(j,k)-b_over_c(j,k)) / (en(j,k)-b_over_c(j,k)) ) )
          endif
          clamda = blconst_ck*zn(j,k)*b(j,k)/(c(j,k)*tau)
          shr(j,k) = -clamda * dvdz2(j,k)
          buoy(j,k) = clamda * buoy_flux(j,k)
          diss(j,k) = (estar(j,k)-en(j,k))/tau-(shr(j,k)+buoy(j,k))
       enddo
    enddo

    ! Set surface values
    estar(:,nk) = 0.
    shr(:,nk) = 0.
    buoy(:,nk) = 0.
    diss(:,nk) = 0.

    ! End of subprogram
    return
  end subroutine tkealg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine clsgs(thl,tve,qw,qc,frac,fnn,fngauss,fnnonloc,c1,zn,ze, &
       hpbl,hflux,s,ps,gztherm,mg,mrk2,acoef,bcoef,ccoef,vcoef,pblsigs,pblq1,n,nk)
    use tdpack_const
    use phy_options
    use ens_perturb, only: ens_nc2d, ens_spp_get
    use pbl_utils, only: blweight, dvrtdf
    implicit none
!!!#include <arch_specific.hf>

    !@Arguments
    integer, intent(in) :: n                             !horizontal dimension
    integer, intent(in) :: nk                            !vertical dimension
    real, dimension(n,nk), intent(in) :: thl             !liquid water potential temperature (K; theta_l)
    real, dimension(n,nk), intent(in) :: tve             !virtual temperature on e-levs (K)
    real, dimension(n,nk), intent(in) :: qw              !total water mixing ratio (kg/kg; q_tot)
    real, dimension(n,nk), intent(in) :: c1              !coefficient C1 in second-order moment closure
    real, dimension(n,nk), intent(in) :: zn              !mixing length (m)
    real, dimension(n,nk), intent(in) :: ze              !dissipation length (m)
    real, dimension(n), intent(in) :: hpbl               !boundary layer height (m)
    real, dimension(n), intent(in) :: hflux              !surface heat flux (W/m2)
    real, dimension(n,nk), intent(in) :: s               !sigma for full levels
    real, dimension(n),    intent(in) :: ps              !surface pressure (Pa)
    real, dimension(n,nk), intent(in) :: gztherm         !height of thermodynamic levels (m)
    real, dimension(n),    intent(in) :: mg              !land-sea mask
    real, dimension(n,ens_nc2d), intent(in) :: mrk2      !Markov chains for stochastic parameters
    real, dimension(n,nk), intent(in) :: acoef           !thermodynamic coefficient A
    real, dimension(n,nk), intent(in) :: bcoef           !thermodynamic coefficient B
    real, dimension(n,nk), intent(in) :: ccoef           !thermodynamic coefficient C
    real, dimension(*), intent(in) :: vcoef              !coefficients for vertical interpolation
    real, dimension(n,nk), intent(inout) :: fnnonloc     !nonlocal cloud fraction
    real, dimension(n,nk), intent(out) :: qc             !PBL cloud water content (kg/kg)
    real, dimension(n,nk), intent(out) :: frac           !cloud fraction
    real, dimension(n,nk), intent(out) :: fnn            !flux enhancement * cloud fraction
    real, dimension(n,nk), intent(out) :: fngauss        !Gaussian cloud fraction
    real, dimension(n,nk), intent(out) :: pblsigs        !Subgrid moisture variance
    real, dimension(n,nk), intent(out) :: pblq1          !Normalized saturation deficit

    !@Author
    !          J. Mailhot (Jun 2002)
    !@Revision
    ! 001      J. Mailhot (Feb 2003) Clipping at upper levels
    ! 002      S. Belair  (Apr 2003) Minimum values of 50 m for ZE and ZN
    !                                in calculation of sigmase.
    ! 003      A-M. Leduc (Jun 2003) Pass ps to blweight ---> blweight2
    ! 004      J. P. Toviessi ( Oct. 2003) - IBM conversion
    !               - calls to vslog routine (from massvp4 library)
    !               - unnecessary calculations removed
    !               - etc.
    ! 005      L. Spacek (Dec 2007) - add "vertical staggering" option
    !                                 change the name to clsgs3
    ! 006      A. Zadra (Oct 2015) -- add land-water mask (MG) to input and
    !               add a new user-defined reduction parameter
    !               (which may or may not depend on the land-water fraction)
    !               to control the flux enhancement factor (fnn)
    ! 007      A. Zadra / R. McT-C (Sep 2016) - Revisions from Adrian Lock
    !               and Jocelyn Mailhot (2012) (non-local scalings,
    !               clips for sigmas and Q1, remove min of 50m for ZE and ZN)
    !@Object Calculate the boundary layer sub-grid-scale cloud properties
    !@Notes
    !          Implicit (i.e. subgrid-scale) cloudiness scheme for unified
    !             description of stratiform and shallow, nonprecipitating
    !             cumulus convection appropriate for a low-order turbulence
    !             model based on Bechtold et al.:
    !            - Bechtold and Siebesma 1998, JAS 55, 888-895 (BS98)
    !            - Cuijpers and Bechtold 1995, JAS 52, 2486-2490
    !            - Bechtold et al. 1995, JAS 52, 455-463 (BCMT95)
    !            - Bechtold et al. 1992, JAS 49, 1723-1744
    !          Revisions for non-local scalings based on:
    !            - Lock and Mailhot 2006, BLM 121, 313-338 (LM06)
    !
    !    The parameters fnn_mask and fnn_reduc are the two parameters
    !    that should be read from the gem_settings
    !       fnn_mask = T/F means "use (do not use) the land-water fraction
    !            to modulate the fnn reduction
    !       fnn_reduc = is the reduction parameter that should vary within
    !            the range 0. to 1.; 1 means "keep the original estimate";
    !            any value smaller than 1 means " multiply the original fnn
    !            by a factor fnn_reduc"
    !
    !            The default values should be
    !              fnn_mask = .false.
    !              fnn_reduc = 1
    !
    !            The values tested and chosen for the RDPS-10km are
    !              fnn_mask = .true.
    !              fnn_reduc = 0.8

    ! Local parameters
    real, parameter :: EPS=1e-10,QCMIN=1e-6,QCMAX=1e-3,WHMIN=0.5,WHMAX=1.2

    ! Local variables
    integer :: j,k
    real :: fnn_weight,sigmas_cu,hfc,lnsig,c1coef
    real(kind=8) :: gravinv
    real, dimension(n) :: wf,fnnreduc
    real, dimension(n,nk) :: dz,dqwdz,dthldz,sigmas,q1,weight,thlm,qwm,wh

    ! Retrieve stochastic information for parameter perturbations
    fnnreduc(:) = ens_spp_get('fnnreduc', mrk2, fnn_reduc)

    ! Pre-compute constants and vertical derivatives
    gravinv = 1.0/dble(GRAV)
    call vint_thermo2mom2(thlm, qwm, &
         &                thl,  qw,  vcoef, n, nk)
    do k=1,nk-1
       do j=1,n
          lnsig   = log(s(j,k+1)/s(j,k))
          dz(j,k) = -RGASD*tve(j,k)*lnsig*gravinv
       end do
    end do
    dz(:,nk) = 0.0
    call dvrtdf( dthldz, thlm, dz, n, n, n, nk)
    call dvrtdf( dqwdz , qwm , dz, n, n, n, nk)

    ! Compute subgrid standard deviation of supersaturation (s) on e-levels following Eq. 10 of BCMT95
    sigmas = 0.
    do k=1,nk-1
       do j=1,n          
          c1coef = sqrt(c1(j,k)*max(zn(j,k),50.)*max(ze(j,k),50.))
          sigmas(j,k) = c1coef * abs(acoef(j,k)*dqwdz(j,k) - bcoef(j,k)*dthldz(j,k))
       end do
    end do

    ! Compute the normalized saturation defecit (Q1)
    q1(:,1:nk-1) = ccoef(:,1:nk-1) / ( sigmas(:,1:nk-1) + EPS )
    q1(:,1) = 0. ; q1(:,nk) = 0. ;
    q1=max(-6.,min(4.,q1))
    sigmas(:,1) = 0.; sigmas(:,nk) = 0.
    pblq1 = q1
    pblsigs = sigmas

    ! Compute cloud properties for local or nonlocal scalings
    do k=2,nk-1
       do j=1,n

          NONLOCAL: if( pbl_nonloc == 'LOCK06' ) then

             ! Compute cloud fractions (LM06)
             if( q1(j,k) > 6.0 ) then
                fngauss(j,k) = 1.0
                fnn(j,k)  = 1.0
             elseif ( q1(j,k) .gt. -6.0 ) then
                ! Represent Bechtold FN function as 0.5(1+tanh(0.8*Q1))
                ! (gives much the same shape as atan but without need for max/min)
                fngauss(j,k) = 0.5*( 1.0 + tanh(0.8*q1(j,k)) )
                ! Use shifted approximate Gaussian for flux enhancement (Eq. 5)
                fnn(j,k)  = 0.5*( 1.0 + tanh(0.8*(q1(j,k)+0.5)) )
             else
                fngauss(j,k) = 0.0
                fnn(j,k)  = 0.0
             endif
             ! Combine Gaussian and non-local cumulus cloud fractions
             frac(j,k) = max(fngauss(j,k),fnnonloc(j,k))

             ! Now increase sigmas to include cumulus contribution in the low cloudiness regime (C<0).
             ! This then gives a less negative Q1 and so more realistic QC
             if ( ccoef(j,k) < 0. .and. fnnonloc(j,k) > 0.0001 .and. fnnonloc(j,k) < 0.4999) then
                !         ! Take inverse of tanh function to find sigmas_cu given FNNONLOC
                sigmas_cu = 1.6 * ccoef(j,k) / alog( fnnonloc(j,k)/(1.-fnnonloc(j,k)) )
                sigmas(j,k) = max( sigmas(j,k), sigmas_cu )
                q1(j,k)     = ccoef(j,k)/sigmas(j,k)
             endif

          else

             ! Compute cloud fraction (BS98, Eq. B1)
             if( q1(j,k) > -1.2) then
                frac(j,k) = max(0.,min(1.,0.5 + 0.36*atan(1.55*q1(j,k))))
             elseif( q1(j,k) >= -6.0) then
                frac(j,k) = exp(q1(j,k)-1.0)
             else
                frac(j,k) = 0.
             endif
             fnnonloc(j,k) = frac(j,k)

          endif NONLOCAL

          ! Compute liquid water specific humidity (BS98, Eq. B2)
          !JM start (modification to LM06 Eq. 6) - no impact
!!$        if( q1(j,k) >= 2.0 ) then
!!$          qc(j,k) = 2.032 + 0.9*(q1(j,k)-2.)
!!$        elseif( q1(j,k) >= 0.0 ) then
          if (q1(j,k) >= 0.0) then
             !JM end
             qc(j,k) = exp(-1.) + 0.66*q1(j,k) + 0.086*q1(j,k)**2
          elseif( q1(j,k) >= -6.0 ) then
             qc(j,k) = exp(1.2*q1(j,k)-1.)
          else
             qc(j,k) = 0.
          endif
          !JM start - no impact
!!$        qc(j,k) = min ( qc(j,k)*sigmas(j,k) , qcmax )
          qc(j,k) = min(qc(j,k)*(sigmas(j,k)+EPS),qcmax)
          !JM end

          ! Compute cloud-induced flux enhancement factor * cloud fraction
          if (pbl_nonloc == 'NIL') then
             TRADE_WIND_CU: if (pbl_cucloud) then
                ! Compute flux enhancement factor (BS98, approx. final Eq. of Appendix B, Fig. 4a)
                fnn(j,k) = 1.0
                if (q1(j,k) < 1.0 .and. q1(j,k) >= -1.68) then
                   fnn(j,k) = exp(-0.3*(q1(j,k)-1.0))
                elseif (q1(j,k) < -1.68 .and. q1(j,k) >= -2.5) then
                   fnn(j,k) = exp(-2.9*(q1(j,k)+1.4))
                elseif (q1(j,k) < -2.5) then
                   fnn(j,k) = 23.9 + exp(-1.6*(q1(j,k)+2.5))
                endif
                ! Adjust for cloud fraction (BS98, Fig. 4b)
                fnn(j,k) = fnn(j,k)*frac(j,k)
                if( q1(j,k).le.-2.39 .and. q1(j,k).ge.-4.0 ) then
                   fnn(j,k) = 0.60
                elseif( q1(j,k).lt.-4.0 .and. q1(j,k).ge.-6.0 ) then
                   fnn(j,k) = 0.30*( q1(j,k)+6.0 )
                elseif( q1(j,k).lt.-6. ) then
                   fnn(j,k) = 0.0
                endif
             else
                ! Without a trade wind cumulus regime, FnN = N as in BS98 Appendix B
                fnn(j,k) = frac(j,k)
             endif TRADE_WIND_CU
          endif

          ! Permit PBL moist processes only under specified conditions
          if (.not.pbl_moistke_legacy_cloud) then
             ! Linear ramp-down of cloud effects from 5 W/m2 to -5 W/m2
             hfc   = 5.
             wf(j) = 0.5*( 1. + hflux(j)/hfc)
             wf(j) = max(0.,min(wf(j),1.))
             frac(j,k)     = wf(j)*frac(j,k)
             fnnonloc(j,k) = wf(j)*fnnonloc(j,k)
             fnn(j,k)      = wf(j)*fnn(j,k)
             qc(j,k)       = wf(j)*qc(j,k)
          endif

          ! Adjust FNN by user-controlled scaling factor potentially land-sea masked
          fnn_weight = 1.
          if (fnn_mask) then
             fnn_weight = fnnreduc(j) + (1.-fnnreduc(j))*max(0.,min(1.,mg(j)))
          else
             fnn_weight = fnnreduc(j)
          endif
          fnn(j,k) = fnn_weight*fnn(j,k)

       end do
    end do

    ! Fill array extrema
    frac(1:n,1) = 0. ; frac(1:n,nk) = 0.
    fnn (1:n,1) = 0. ; fnn (1:n,nk) = 0.
    qc  (1:n,1) = 0. ; qc  (1:n,nk) = 0.
    fnnonloc(1:n,1) = 0.; fnnonloc(1:n,nk) = 0.
    fngauss(1:n,1) = 0.; fngauss(1:n,nk) = 0.

    ! Create a vertical function to eliminate PBL cloud effects at upper levels
    call blweight(weight,s,ps,n,nk)

    ! Do not allow PBL cloud effects to reach beyond 1.5 times the PBL depth
    if (.not.pbl_moistke_legacy_cloud) then
       do k=1,nk
          ! Linear ramp-down of cloud effects from 1.5*hpbl to 2*hpbl
          wh(:,k) = (2.*hpbl(:) - gztherm(:,k))/(0.5*hpbl(:))
          wh(:,k) = max(0., min(1., wh(:,k)))
          weight(:,k) = weight(:,k)*wh(:,k)
       enddo
    endif

    ! Scale all PBL cloud quantities
    frac = frac*weight
    fnnonloc = fnnonloc*weight
    fnn = fnn*weight
    qc = qc*weight
    fngauss = fngauss*weight

    return
  end subroutine clsgs
  
end module pbl_mtke_utils
