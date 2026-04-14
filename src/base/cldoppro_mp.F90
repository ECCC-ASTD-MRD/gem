
module cldoppro_MP
   implicit none
   private
   public :: cldoppro_MP3

contains

   !/@*
  subroutine cldoppro_MP3(pvars, &
       taucs, omcs, gcs, taucl, omcl, gcl, &
       liqwcin, icewcin, &
       liqwpin, icewpin, cldfrac, &
       tt, sig, ps, dontneedall, ni, nkm1, nk, kount)
    use, intrinsic :: iso_fortran_env, only: INT64
    use debug_mod, only: init2nan
    use tdpack_const, only: GRAV, RGASD, TCDK
    use phy_options
    use phy_status, only: physeterror
    use phybusidx
    use phymem, only: phyvar
    use cldwin_mod, only: cldwin1
    use ens_perturb, only: ens_nc2d, ens_spp_get
    use cldop_utils, only: cldop_rei, cldop_rew, cldop_proplw, cldop_propli, &
           cldop_propsw, cldop_propsi
    implicit none

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include "nbsnbl.cdk"

    !@Arguments
    integer, intent(in) :: ni, nkm1, nk, kount

    type(phyvar), pointer, contiguous :: pvars(:)
    real, intent(out), dimension(ni,nkm1,NBS) :: taucs, omcs, gcs
    real, intent(out), dimension(ni,nkm1,NBL) :: taucl, omcl, gcl
    real, intent(inout), dimension(ni,nkm1) :: liqwcin, icewcin
    real, intent(inout), dimension(ni,nkm1) :: liqwpin, icewpin
    real, intent(inout) :: cldfrac(ni,nkm1)
    real, intent(in)    :: tt(ni,nkm1), sig(ni,nkm1), ps(ni)
    logical, intent(in) :: dontneedall                   ! .true. = need only band1 in SW and band6 in LW

    !          - input/output -
    ! pvars    list of all phy vars (meta + slab data)
    !          - output -
    ! taucs    cloud solar optical thickness
    ! omcs     cloud solar scattering albedo
    ! gcs      cloud solar asymmetry factor
    ! taucl    cloud longwave optical thickness
    ! omcl     cloud longwave scattering albedo
    ! gcl      cloud longwave asymmetry factor
    ! topthw   total integrated optical thickness of water in the visible
    ! topthi   total integrated optical thickness of ice in the visible
    ! ctp      cloud top pressure
    ! ctt      cloud top temperature
    ! ecc      effective cloud cover (nt)
    !          - input / output -
    ! liqwpin  liquid water path in g/m2
    ! icewpin  solid water path in g/m2
    ! liqwcin  in-cloud liquid water content
    !          in kg water/kg air
    ! icewcin  in-cloud ice water content
    !          in kg water/kg air
    ! cldfrac  layer cloud amount (0. to 1.) (ni,nkm1)
    !          - input -      
    ! tt       layer temperature (k) (ni,nkm1)
    ! sig      sigma levels (0. to 1.) (ni,nkm1; local sigma)
    ! ps       surface pressure (n/m2) (ni)
    ! mg       ground cover (ocean=0.0,land <= 1.)  (ni)
    ! ml       fraction of lakes (0.-1.) (ni)
    ! ni       first dimension of temperature (usually ni)
    ! nkm1       number of layers
    ! kount    time step number
    !@Authors
    !        D. Paquin-Ricard june 2017 continuing adaptation
    !        p. vaillancourt, mai 2016 (merged simplified/adapted prep_cw_rad ,cldroppro and diagno_cw_rad)
    !
    !@Object
    !        calculate cloud optical properties for ccc radiative transfer scheme for MP schemes
    !        calculate vertical integrals of liquid/ice hydrometeors
    !        output merged fields of liquid water content,  ice water content and cloud fraction
    !*@/
#include <rmn/msg.h>
#include "phymkptr.hf"

    include "phyinput.inc"
    include "surface.cdk"
    include "nocld.cdk"

    logical, dimension(ni,nkm1) :: nocloud

    real, dimension(:), allocatable :: tausimp, omsimp, gsimp, taulimp, omlimp, glimp

    real, dimension(ni,nkm1) :: rew, rei, rewxp, reixp, vs1, dp
    real, dimension(ni,nkm1) :: lwpinmp, cldfmp, cldfxp, lwcinmp, iwpinmps, iwcinmps
    real, dimension(ni,nkm1) :: wexp,wimp
    real, dimension(ni) :: reifac
    real, dimension(ni) :: rewfac

    real, dimension(:,:,:), allocatable :: iwcinmp, iwpinmp, effradi

    logical :: nostrlwc, readfield_L
    integer :: i, j, k, l, mpcat, istat1, istat2, swi,swf,lwi,lwf
    real :: rec_grav,  press
    real :: rew1, dg, tausw, omsw, gsw, tausi, omsi, gsi, y1, y2, y3
    real :: taulw, omlw, glw, tauli, omli, gli
    real :: tauswmp, omswmp, gswmp
    real :: taulwmp, omlwmp, glwmp
    real :: mpcldth
    real :: epsilon, epsilon2, betan, betad

    real, pointer, dimension(:,:), contiguous :: zqcmoins, zqimoins, zqnmoins, zqgmoins
    real, pointer, dimension(:,:), contiguous :: zqti1m, zqti2m, zqti3m, zqti4m
    real, pointer, dimension(:,:), contiguous :: zeffradc, zeffradi1 , zeffradi2 , zeffradi3 , zeffradi4
    real, pointer, dimension(:), contiguous :: ztopthw,ztopthi
    real, pointer, dimension(:,:), contiguous :: ztose,ztosi,ztole,ztoli
    real, pointer, dimension(:), contiguous :: zmg,zml,zdlat
    real, pointer, dimension(:,:), contiguous :: ziwcimp,zlwcimp,zhumoins,ztmoins,zpmoins,zsigw,zfxp,zfmp
    real, pointer, dimension(:), contiguous   :: ztlwp, ztiwp,ztlwpin,ztiwpin
    real, pointer, dimension(:,:), contiguous   :: zrewx, zreix, zrewi, zreii
    real, pointer, dimension(:,:), contiguous :: zlwcrad, ziwcrad, zcldrad, zqcplus, zgraupel, &
         zqiplus, zsnow, zqi_cat1, zqi_cat2, zqi_cat3, zqi_cat4,zmrk2
    !----------------------------------------------------------------
    call msg_toall(MSG_DEBUG, 'cldoppro_MP [BEGIN]')

    MKPTR1D(zdlat, dlat, pvars)
    MKPTR1D(zmg, mg, pvars)
    MKPTR1D(zml, ml, pvars)
    MKPTR1D(ztiwp, tiwp, pvars)
    MKPTR1D(ztiwpin, tiwpin, pvars)
    MKPTR1D(ztlwp, tlwp, pvars)
    MKPTR1D(ztlwpin, tlwpin, pvars)
    MKPTR1D(ztopthi, topthi, pvars)
    MKPTR1D(ztopthw, topthw, pvars)
    MKPTR2Dm1(zrewx, rewx, pvars)
    MKPTR2Dm1(zreix, reix, pvars)
    MKPTR2Dm1(zrewi, rewi, pvars)
    MKPTR2Dm1(zreii, reii, pvars)
    MKPTR2Dm1(ztose, tose, pvars)
    MKPTR2Dm1(ztosi, tosi, pvars)
    MKPTR2Dm1(ztole, tole, pvars)
    MKPTR2Dm1(ztoli, toli, pvars)
    MKPTR2Dm1(zcldrad, cldrad, pvars)
    MKPTR2Dm1(zeffradc, effradc, pvars)
    MKPTR2Dm1(zeffradi1, effradi1, pvars)
    MKPTR2Dm1(zeffradi2, effradi2, pvars)
    MKPTR2Dm1(zeffradi3, effradi3, pvars)
    MKPTR2Dm1(zeffradi4, effradi4, pvars)
    MKPTR2Dm1(zfmp, fmp, pvars)
    MKPTR2Dm1(zfxp, fxp, pvars)
    MKPTR2Dm1(zgraupel, qgplus, pvars)
    MKPTR2Dm1(zhumoins, humoins, pvars)
    MKPTR2Dm1(zqti1m, qti1moins, pvars)
    MKPTR2Dm1(zqti2m, qti2moins, pvars)
    MKPTR2Dm1(zqti3m, qti3moins, pvars)
    MKPTR2Dm1(zqti4m, qti4moins, pvars)
    MKPTR2Dm1(ziwcimp, iwcimp, pvars)
    MKPTR2Dm1(ziwcrad, iwcrad, pvars)
    MKPTR2Dm1(zlwcimp, lwcimp, pvars)
    MKPTR2Dm1(zlwcrad, lwcrad, pvars)
    MKPTR2Dm1(zpmoins, pmoins, pvars)
    MKPTR2Dm1(zqcmoins, qcmoins, pvars)
    MKPTR2Dm1(zqcplus, qcplus, pvars)
    MKPTR2Dm1(zqgmoins, qgmoins, pvars)
    MKPTR2Dm1(zqi_cat1, qti1plus, pvars)
    MKPTR2Dm1(zqi_cat2, qti2plus, pvars)
    MKPTR2Dm1(zqi_cat3, qti3plus, pvars)
    MKPTR2Dm1(zqi_cat4, qti4plus, pvars)
    MKPTR2Dm1(zqimoins, qimoins, pvars)
    MKPTR2Dm1(zqiplus, qiplus, pvars)
    MKPTR2Dm1(zqnmoins, qnmoins, pvars)
    MKPTR2Dm1(zsigw, sigw, pvars)
    MKPTR2Dm1(zsnow, qnplus, pvars)
    MKPTR2Dm1(ztmoins, tmoins, pvars)
    MKPTR2D(zmrk2, mrk2, pvars)

    call init2nan(rew, rei, rewxp, reixp, vs1, dp)
    call init2nan(lwpinmp, cldfmp, cldfxp, lwcinmp, iwpinmps, iwcinmps)
    call init2nan(reifac, rewfac)

    ! Allocate cateory-dependent space
    if (stcond(1:5)=='MP_P3')  mpcat = p3_ncat
    if (stcond(1:6)=='MP_MY2') mpcat = 3
    allocate(tausimp(mpcat), omsimp(mpcat), gsimp(mpcat), taulimp(mpcat), &
         omlimp(mpcat), glimp(mpcat), stat=istat1)
    allocate(iwcinmp(ni,nkm1,mpcat), iwpinmp(ni,nkm1,mpcat), &
         effradi(ni,nkm1,mpcat), stat=istat2)
    if (istat1 /= 0 .or. istat2 /= 0) then
       call physeterror('cldoppro_mp', &
            'Cannot allocate category-dependent tables')
       return
    endif

    call init2nan(iwcinmp, iwpinmp, effradi)    
    call init2nan(tausimp, omsimp, gsimp, taulimp, omlimp, glimp)    

    rec_grav = 1./GRAV
    nostrlwc = (climat.or.stratos)
    zrewx    = 0.0
    zreix    = 0.0
    zrewi    = 0.0
    zreii    = 0.0
    ztlwp    = 0.0
    ztiwp    = 0.0
    ztlwpin  = 0.0
    ztiwpin  = 0.0
    zlwcrad  = 0.0
    ziwcrad  = 0.0
    zcldrad  = 0.0
    iwcinmps = 0.0
    iwcinmp  = 0.0
    lwcinmp  = 0.0

    ! Set cloud detection threshold
    if (rad_mpagg == 'BINARY') then
       mpcldth = 0.01 
    else
       mpcldth = cldfth ! use standard parameter
    endif

    ! Assume implicit sources have been agregated into zlwc and ziwc in subroutine prep_cw_MP
    do k=1,nkm1
       do i=1,ni
          liqwcin(i,k) = max(zlwcimp(i,k), 0.)
          icewcin(i,k) = max(ziwcimp(i,k), 0.)
          lwcinmp(i,k) = max(zqcmoins(i,k), 0.)
          cldfmp(i,k) = zfmp(i,k)  !implicit
          cldfxp(i,k) = zfxp(i,k)  !explicit
       enddo
    enddo

    ! Reshape ice mass and radius tables
    if (associated(zqti1m)) then
       iwcinmp(:,:,1) = max(zqti1m(:,:), 0.)
       effradi(:,:,1) = max(zeffradi1(:,:), 0.)
    endif
    if (associated(zqti2m)) then
       iwcinmp(:,:,2) = max(zqti2m(:,:), 0.)
       effradi(:,:,2) = max(zeffradi2(:,:), 0.)
    endif
    if (associated(zqti3m)) then
       iwcinmp(:,:,3) = max(zqti3m(:,:), 0.)
       effradi(:,:,3) = max(zeffradi3(:,:), 0.)
    endif
    if (associated(zqti4m)) then
       iwcinmp(:,:,4) = max(zqti4m(:,:), 0.)
       effradi(:,:,4) = max(zeffradi4(:,:), 0.)
    endif
    if (associated(zqimoins)) then
       iwcinmp(:,:,1) = max(zqimoins(:,:), 0.)
       effradi(:,:,1) = max(zeffradi1(:,:), 0.)
    endif
    if (associated(zqnmoins)) then
       iwcinmp(:,:,2) = max(zqnmoins(:,:), 0.)
       effradi(:,:,2) = max(zeffradi2(:,:), 0.)
    endif
    if (associated(zqgmoins)) then
       iwcinmp(:,:,3) = max(zqgmoins(:,:), 0.)
       effradi(:,:,3) = max(zeffradi3(:,:), 0.)
    endif

    if (kount == 0) then
       if (inilwc) then
          if (stcond /= 'NIL') then
             ! initialiser le champ d'eau nuageuse ainsi que la fraction
             ! nuageuse pour l'appel a la radiation a kount=0seulement.
             ! ces valeurs seront remplacees par celles calculees dans
             ! les modules de condensation.
             call cldwin1(zfmp,zlwcimp,ztmoins,zhumoins,zpmoins, &
                  zsigw,ni,nkm1,satuco)
             do k=1,nkm1
                do i=1,ni
                   cldfmp(i,k) = zfmp(i,k)
                   liqwcin(i,k) = zlwcimp(i,k)
                enddo
             enddo
          endif
       endif
    endif

    readfield_L = .false.
    if (ISDYNIN('qc') .or. ISPHYIN('tr/mpqc:p')) then
       readfield_L = .true.
       do k=1,nkm1
          do i=1,ni
             lwcinmp(i,k) = max(zqcplus(i,k), 0.)
          enddo
       enddo
    endif

    IF_MY2: if (stcond(1:6) == 'MP_MY2' .and. &
         (ISDYNIN('mpqi') .or. ISPHYIN('tr/mpqi:p')) .and. &
         (ISDYNIN('mpqs') .or. ISPHYIN('tr/mpqs:p')) .and. &
         (ISDYNIN('mpqg') .or. ISPHYIN('tr/mpqg:p'))) then

       readfield_L = .true.
       !ziwc = zqiplus + zsnow
       do k=1,nkm1
          do i=1,ni
             iwcinmp(i,k,1) = max(zqiplus(i,k), 0.)
             iwcinmp(i,k,2) = max(zsnow(i,k), 0.)
             iwcinmp(i,k,3) = max(zgraupel(i,k), 0.)
          enddo
       enddo

    elseif (stcond(1:5) == 'MP_P3' .and. &
         (ISDYNIN('qti1') .or. ISPHYIN('tr/qti1:p'))) then

       readfield_L = .true.
       !ziwc = zqi_cat1
       do k=1,nkm1
          do i=1,ni
             iwcinmp(i,k,1) = max(zqi_cat1(i,k), 0.)
          enddo
       enddo

       IF_NCAT2: if (p3_ncat >= 2 .and. &
            (ISDYNIN('qti2') .or. ISPHYIN('tr/qti2:p'))) then

          !ziwc = ziwc + zqi_cat2
          do k=1,nkm1
             do i=1,ni
                iwcinmp(i,k,2) = max(zqi_cat2(i,k), 0.)
             enddo
          enddo

          IF_NCAT3: if (p3_ncat >= 3 .and. &
               (ISDYNIN('qti3') .or. ISPHYIN('tr/qti3:p'))) then

             !ziwc = ziwc + zqi_cat3
             do k=1,nkm1
                do i=1,ni
                   iwcinmp(i,k,3) = max(zqi_cat3(i,k), 0.)
                enddo
             enddo

             IF_NCAT4: if (p3_ncat >= 4 .and. &
                  (ISDYNIN('qti4') .or. ISPHYIN('tr/qti4:p'))) then

                !ziwc = ziwc + zqi_cat4
                do k=1,nkm1
                   do i=1,ni
                      iwcinmp(i,k,4) = max(zqi_cat4(i,k), 0.)
                   enddo
                enddo

             endif IF_NCAT4
          endif IF_NCAT3
       endif IF_NCAT2
    endif IF_MY2

    if (readfield_L) then
       do k=1,nkm1
          do i=1,ni
             do l=1,mpcat
                iwcinmps(i,k) = iwcinmps(i,k) + iwcinmp(i,k,l)
             enddo
             if (lwcinmp(i,k) + iwcinmps(i,k) > 1.e-6) then   ! not coherent with the rest
                cldfxp(i,k) = 1
             else
                cldfxp(i,k) = 0
             endif
          enddo
       enddo
    endif

    !..."no stratospheric lwc" mode when CLIMAT or STRATOS = true
    !...no clouds above TOPC (see nocld.cdk)
    if (nostrlwc) then
       do k=1,nkm1
          do i=1,ni
             press = sig(i,k)*ps(i)
             if (topc > press) then
                cldfxp(i,k) = 0.0   ! explicit
                cldfmp(i,k) = 0.0   ! implicit
                liqwcin (i,k) = 0.0 ! implicit
                icewcin (i,k) = 0.0 ! implicit
                lwcinmp (i,k) = 0.0 ! explicit
                do l=1,mpcat
                   iwcinmp (i,k,l) = 0.0 !explicit
                enddo
             endif
          enddo
       enddo
    endif

    ! Filter exp and imp clouds fracs below mpcldth
    ! Normalize water contents to get in-cloud values
    DO_K: do k=1,nkm1
       do i=1,ni
          if (cldfmp(i,k) < mpcldth) then 
             liqwcin(i,k) = 0.
             icewcin(i,k) = 0.
             cldfmp(i,k) = 0.
          else
             liqwcin(i,k) = liqwcin(i,k)/cldfmp(i,k)
             icewcin(i,k) = icewcin(i,k)/cldfmp(i,k)
          endif

          ! remove ice water content if effective radius >= 1e-4,
          ! since it will not be seen by radiation 
          !best choice? how frequent? rei is limited to 70microns lower down
          ! test showed practically no impact of this filtering
          do l=1,mpcat
             if (kount == 0) effradi(i,k,l) = 5.e-5
             if (effradi(i,k,l) >= 1.e-4) then  
                iwcinmp(i,k,l) = 0.0
             endif
          enddo
          if (cldfxp(i,k) < mpcldth) then 
             lwcinmp(i,k) = 0.
             do l=1,mpcat
                iwcinmp(i,k,l) = 0.0
             enddo
             cldfxp(i,k) = 0.
          else
             lwcinmp(i,k) = lwcinmp(i,k)/cldfxp(i,k)
             do l=1,mpcat
                iwcinmp(i,k,l) = iwcinmp(i,k,l)/cldfxp(i,k)
             enddo
          endif
       end do
    end do DO_K

    ! Calculate in-cloud liquid and ice water paths in each layer (converted to g/m^2)

    do i=1,ni
       dp(i,1) = 0.5*(sig(i,1) + sig(i,2))*rec_grav*1000.
       dp(i,nkm1) = (1. - 0.5*(sig(i,nkm1) + sig(i,nkm1-1)))*rec_grav*1000.
       dp(i,1) = max(dp(i,1)*ps(i), 0.)
       dp(i,nkm1) = max(dp(i,nkm1)*ps(i), 0.)
    end do
    do k=2,nkm1-1
       do i=1,ni
          dp(i,k) = 0.5*(sig(i,k+1) - sig(i,k-1))*rec_grav*1000.
          dp(i,k) = max(dp(i,k)*ps(i), 0.)
       end do
    end do

    do k=1,nkm1
       do i=1,ni
          icewpin(i,k)  = icewcin(i,k)*dp(i,k) !imp
          liqwpin(i,k)  = liqwcin(i,k)*dp(i,k) !imp
          lwpinmp(i,k)  = lwcinmp(i,k)*dp(i,k) !exp
          iwpinmps(i,k) = 0.
          iwcinmps(i,k) = 0.
       enddo
    enddo

    do k=1,nkm1
       do i=1,ni
          do l=1,mpcat
             iwpinmp(i,k,l) = iwcinmp(i,k,l)*dp(i,k)
             iwpinmps(i,k) = iwpinmps(i,k) + iwpinmp(i,k,l) !exp
             iwcinmps(i,k) = iwcinmps(i,k) + iwcinmp(i,k,l) !exp
          enddo
       end do
    end do

    do k=1,nkm1
       do i=1,ni
          ztlwpin(i) = ztlwpin(i) + liqwpin(i,k) + lwpinmp(i,k)
          ztiwpin(i) = ztiwpin(i) + icewpin(i,k) + iwpinmps(i,k)
       enddo
    enddo

    do k=1,nkm1
       do i=1,ni
          if (lwpinmp(i,k) <= wpth .and. iwpinmps(i,k) <= wpth) then
             cldfxp(i,k) = 0.
          endif
          if (liqwpin(i,k) <= wpth .and. icewpin(i,k) <= wpth) then
             cldfmp(i,k) = 0.
          endif
!          nocloud(i,k) = (liqwpin(i,k)+lwpinmp(i,k)+icewpin(i,k)+iwpinmps(i,k) <= wpth)    ! would be more coherent if based on cldfxp and cldfmp
          nocloud(i,k) = (cldfxp(i,k) < mpcldth .and. cldfmp(i,k) < mpcldth)    ! would be more coherent if based on cldfxp and cldfmp
          ztlwp(i)     = ztlwp(i)   + liqwpin(i,k)*cldfmp(i,k) + lwpinmp(i,k)*cldfxp(i,k)  !diag output
          ztiwp(i)     = ztiwp(i)   + icewpin(i,k)*cldfmp(i,k) + iwpinmps(i,k)*cldfxp(i,k) !diag output
          zlwcrad(i,k) = liqwcin(i,k)*cldfmp(i,k) + lwcinmp(i,k)*cldfxp(i,k)               !diag output
          ziwcrad(i,k) = icewcin(i,k)*cldfmp(i,k) + iwcinmps(i,k)*cldfxp(i,k)              !diag output
          cldfrac(i,k) = max(min(cldfmp(i,k)+cldfxp(i,k), 1.0), 0.0)       !no overlap=ops !diag output
!          cldfrac(i,k) = max(min(max(cldfmp(i,k),cldfxp(i,k)), 1.0), 0.0)                ! max overlap
!          cldfrac(i,k) = min(1., max(0.,  1. - (1.-cldfmp(i,k))*(1.-cldfxp(i,k))))       ! random overlap
          if (nocloud(i,k)) then
             cldfrac(i,k) = 0.0
             zlwcrad(i,k) = 0.0
             ziwcrad(i,k) = 0.0
          endif
          zcldrad(i,k) = cldfrac(i,k)
       end do
    end do

    ! Set weights for aggregation of clouds from different sources
    wexp(:,:) = 1.
    wimp(:,:) = 1.
    select case (rad_mpagg)
    case ('RELATIVE') ! use relative weighting to combine optical thicknesses of imp and exp clouds
       do k=1,nkm1
          do i=1,ni
             if((cldfmp(i,k)+cldfxp(i,k)) >= mpcldth) then
                wexp(i,k) = cldfxp(i,k)/(cldfmp(i,k)+cldfxp(i,k))
                wimp(i,k) = cldfmp(i,k)/(cldfmp(i,k)+cldfxp(i,k))
             endif
          end do
       end do
    case ('MERGED') ! use scheme-generated cloud fractions
       do k=1,nkm1
          do i=1,ni
             if((cldfmp(i,k)+cldfxp(i,k)) >= mpcldth) then
                wexp(i,k) = cldfxp(i,k)/min(cldfmp(i,k)+cldfxp(i,k), 1.)
                wimp(i,k) = cldfmp(i,k)/min(cldfmp(i,k)+cldfxp(i,k), 1.)
             endif
          end do
       end do
    end select

    !     conversion d'unites : tlwp et tiwp en kg/m2
    do i=1,ni
       ztlwp(i) = ztlwp(i) * 0.001
       ztiwp(i) = ztiwp(i) * 0.001
    enddo

    ! Initialize output fields
    do i = 1, ni
       ztopthw(i) = 0.0
       ztopthi(i) = 0.0
    end do

    ! Effective radius for IMPLICIT WATER clouds
    rew(:,:) = cldop_rew(tt, liqwcin, sig, ps, zmg, zml, rew_const, rad_cond_rew, ni, nkm1)
    zrewi(:,:) = rew(:,:) * 1e-6

    ! Adjust the effective radius using stochastic perturbations
    rewfac(:) = ens_spp_get('rew_mult', zmrk2, default=1.)
    do k=1, nkm1
       rew(:,k) = rewfac(:) * rew(:,k)
    enddo

    ! Effective radius for EXPLICIT WATER clouds
    rewxp(:,:) = cldop_rew(tt, lwcinmp, sig, ps, zmg, zml, rewx_const, rad_exp_rew, ni, nkm1, remp=zeffradc)
    zrewx(:,:) = rewxp(:,:) * 1e-6

    ! Effective radius for IMPLICIT ICE clouds
    rei(:,:) = cldop_rei(tt, icewcin, sig, ps, zdlat, rei_const, rad_cond_rei, ni, nkm1)
    zreii(:,:) = rei(:,:) * 1e-6
    
    ! Adjust the effective radius using stochastic perturbations
    reifac(:) = ens_spp_get('rei_mult', zmrk2, default=1.)
    do k=1, nkm1
       rei(:,k) = reifac(:) * rei(:,k)
    enddo
    
    ! Effective radius for EXPLICIT ICE clouds
    reixp(:,:) = cldop_rei(tt, iwcinmp, sig, ps, zdlat, reix_const, rad_exp_rei, ni, nkm1, remp=effradi(:,:,1))
    zreix(:,:) = reixp(:,:) * 1e-6
   
    ! mask output effective radii where there is no cloud
    where (cldfxp(:,:) < mpcldth)
       zreix(:,:)= 0.0
       zrewx(:,:)= 0.0
    endwhere
    where (cldfmp(:,:) < mpcldth)
       zreii(:,:)= 0.0
       zrewi(:,:)= 0.0
    endwhere
    
    ! end-code : FOR EFFECTIVE RADII OF IMPLICIT CLOUDS (NON-mp SOURCES)
    !
    !----------------------------------------------------------------------
    !     cloud radiative properties for radiation.
    !     taucs, omcs, gcs (taucl, omcl, gcl): optical depth, single
    !     scattering albedo, asymmetry factor for solar (infrared).
    !     rew: effective radius (in micrometer) for water cloud
    !     rei: effective radius (in micrometer) for ice cloud
    !     dg: geometry length for ice cloud
    !     liqwcin  (icewcin): liquid water (ice) content (in gram / m^3)
    !     liqwpin (icewpin): liquid water (ice) path length (in gram / m^2)
    !     cloud: cloud fraction
    !     parameterization for water cloud:
    !     dobbie, etc. 1999, jgr, 104, 2067-2079
    !     lindner, t. h. and j. li., 2000, j. clim., 13, 1797-1805.
    !     parameterization for ice cloud:
    !     fu 1996, j. clim., 9, 2223-2337.
    !     fu et al. 1998 j. clim., 11, 2223-2337.
    !     NOTE:  In two Fu papers, opt prop are tested for DG from 15-130microns or REI from 10-90 microns(using dg=1.5395xrei)
    !            but in practice, rei must be below 70microns or the model crashes
    !----------------------------------------------------------------------

    if(dontneedall) then
     swi=1;swf=1;lwi=6;lwf=6
    else
     swi=1;swf=nbs;lwi=1;lwf=nbl
    endif

    ! Shortwave cloud radiative properties
    DO_NBS: do j = swi, swf
       do k = 1, nkm1
          do i = 1, ni
             IF_NOCLOUD: if (nocloud(i,k)) then
                taucs(i,k,j) = 0.
                omcs(i,k,j)  = 0.
                gcs(i,k,j)   = 0.
             else

                ! Shortwave properties of implicit liquid water clouds
                call cldop_propsw(tausw, omsw, gsw, liqwpin(i,k), rew(i,k), j, liqwpin(i,k) > wpth)

                ! Shortwave properties of implicit ice clouds
                call cldop_propsi(tausi, omsi, gsi, icewpin(i,k), rei(i,k), j, icewpin(i,k) > wpth)

                ! Shortwave properties of explicit liquid water clouds
                rew1 = rewxp(i,k)
                if (kount == 0) rew1 = 10.    ![microns] assign value at step zero (in case QC is non-zero at initial conditions)
                rew1 = min(max(4., rew1), 40.0) !where does this 40 come from?
                call cldop_propsw(tauswmp, omswmp, gswmp, lwpinmp(i,k), rew1, j, lwpinmp(i,k) > wpth)

                ! Shortwave properties of explicit ice clouds
                do l=1,mpcat
                   call cldop_propsi(tausimp(l), omsimp(l), gsimp(l), iwpinmp(i,k,l), reixp(i,k), j, &
                        iwpinmp(i,k,l) > wpth .and. effradi(i,k,l) < 1.e-4, rlim=(/15.,110./))
                enddo

                !PV agregate SW optical properties for liq-imp + ice-imp + liq-exp + ice-exp

                ! cccma recipe  
                !              taucs(i,k,j)  = tausw + tausi
                !                y1          = omsw * tausw
                !                y2          = omsi * tausi
                !                omcs(i,k,j) = (y1 + y2) / taucs(i,k,j)
                !                gcs (i,k,j) = (y1 * gsw + y2 * gsi) / (y1 + y2)

                tausw=tausw*wimp(i,k)
                tausi=tausi*wimp(i,k)
                tauswmp=tauswmp*wexp(i,k)

                y1 = tausw + tausi + tauswmp
                y2 = tausw*omsw + tausi*omsi + tauswmp*omswmp
                y3 = tausw*omsw*gsw + tausi*omsi*gsi + tauswmp*omswmp*gswmp

                do l=1,mpcat
                   tausimp(l)=tausimp(l)*wexp(i,k)
                   y1 = y1 +  tausimp(l)
                   y2 = y2 +  tausimp(l)*omsimp(l)
                   y3 = y3 +  tausimp(l)*omsimp(l)*gsimp(l)
                enddo

                taucs(i,k,j) = y1

                if (y1 > 0.0) then
                   omcs(i,k,j) = y2/y1
                   gcs (i,k,j) = y3/y2
                else
                   omcs(i,k,j) = 0.
                   gcs (i,k,j) = 0.
                endif

                !  optical depth for water and ice cloud in visible for output
                if (j == 1) then
                   ztopthw(i) = ztopthw(i) + tausw + tauswmp
                   ztopthi(i) = ztopthi(i) + tausi
                   do l=1,mpcat
                      ztopthi(i) = ztopthi(i) + tausimp(l)
                   enddo
                   if (associated(ztole)) ztole(i,k) = tauswmp
                   if (associated(ztose)) ztose(i,k) = tausimp(1)
                   if (associated(ztoli)) ztoli(i,k) = tausw
                   if (associated(ztosi)) ztosi(i,k) = tausi
                endif

             endif IF_NOCLOUD

          enddo
       enddo
    enddo DO_NBS

    ! Longwave cloud radiative properties
    DO_NBL: do j = lwi, lwf
       do k = 1, nkm1
          do i = 1, ni
             IF_NOCLOUD2: if (nocloud(i,k)) then
                taucl(i,k,j) = 0.
                omcl(i,k,j)  = 0.
                gcl(i,k,j)   = 0.
             else

                ! Longwave properties of implicit liquid water clouds
                call cldop_proplw(taulw, omlw, glw, liqwpin(i,k), rew(i,k), j, liqwpin(i,k) > wpth)

                ! Longwave properties of implicit ice clouds
                call cldop_propli(tauli, omli, gli, icewpin(i,k), rei(i,k), j, icewpin(i,k) > wpth)

                ! Longwave properties of explicit liquid water clouds
                rew1 = rewxp(i,k)
                if (kount == 0) rew1 = 10.    ![microns] aslign value at step zero (in case QC is non-zero at initial conditions)
                rew1 = min(max(4., rew1), 40.0) !where does this 40 come from?
                call cldop_proplw(taulwmp, omlwmp, glwmp, lwpinmp(i,k), rew1, j, lwpinmp(i,k) > wpth)

                ! Longwave properties of explicit ice clouds
                dg = reixp(i,k)
                if (kount == 0) dg = 50.
                do l=1,mpcat
                   call cldop_propli(taulimp(l), omlimp(l), glimp(l), iwpinmp(i,k,l), dg, j, &
                        iwpinmp(i,k,l) > wpth .and. effradi(i,k,l) < 1.e-4, rlim=(/15.,110./))
                enddo

                ! Agregate LW optical properties for liq-imp + ice-imp + liq-exp + ice-exp
                taulw=taulw*wimp(i,k)
                tauli=tauli*wimp(i,k)
                taulwmp=taulwmp*wexp(i,k)
                y1 = taulw + tauli + taulwmp
                y2 = taulw*omlw + tauli*omli + taulwmp*omlwmp
                y3 = taulw*omlw*glw + tauli*omli*gli + taulwmp*omlwmp*glwmp

                do l=1,mpcat
                   taulimp(l)=taulimp(l)*wexp(i,k)
                   y1 = y1 +  taulimp(l)
                   y2 = y2 +  taulimp(l)*omlimp(l)
                   y3 = y3 +  taulimp(l)*omlimp(l)*glimp(l)
                enddo

                taucl(i,k,j) = y1

                if (taucl(i,k,j) < 0.) then
                   print*, '*** taucl: ',taucl(i,k,j),taulw,tauli,taulwmp,taulimp(1)  !#TODO: avoid printing in MPI/OMP regions, render the listing unusable
                   call physeterror('cldoppro_MP', 'Found negative taucl values.')
                   return
                endif

                if (y1 > 0.0) then
                   omcl(i,k,j) = y2/y1
                   gcl (i,k,j) = y3/y2
                else
                   omcl(i,k,j) = 0.
                   gcl (i,k,j) = 0.
                endif

             endif IF_NOCLOUD2
          enddo
       enddo
    enddo DO_NBL

    ! Garbage collection
    deallocate(tausimp, omsimp, gsimp, taulimp, omlimp, glimp)
    deallocate(iwcinmp, iwpinmp, effradi)
    
    call msg_toall(MSG_DEBUG, 'cldoppro_MP [END]')
    !----------------------------------------------------------------
    return
  end subroutine cldoppro_MP3

end module cldoppro_MP

