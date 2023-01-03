!-------------------------------------a LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

module cldoppro_MP
   implicit none
   private
   public :: cldoppro_MP3

contains

   !/@*
  subroutine cldoppro_MP3(dbus, fbus, vbus, &
       taucs, omcs, gcs, taucl, omcl, gcl, &
       liqwcin, icewcin, &
       liqwpin, icewpin, cldfrac, &
       tt, sig, ps, ni, nkm1, nk, kount)
    use, intrinsic :: iso_fortran_env, only: INT64
    use debug_mod, only: init2nan
    use tdpack_const, only: GRAV, RGASD
    use phy_options
    use phybus
    implicit none

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
    include "nbsnbl.cdk"

    !@Arguments
    integer, intent(in) :: ni, nkm1, nk, kount

    real, pointer, contiguous :: dbus(:), fbus(:), vbus(:)
    real, intent(out), dimension(ni,nkm1,nbs) :: taucs, omcs, gcs
    real, intent(out), dimension(ni,nkm1,nbl) :: taucl, omcl, gcl
    real, intent(inout), dimension(ni,nkm1) :: liqwcin, icewcin
    real, intent(inout), dimension(ni,nkm1) :: liqwpin, icewpin
    real, intent(inout) :: cldfrac(ni,nkm1)
    real, intent(in)    :: tt(ni,nkm1), sig(ni,nkm1), ps(ni)

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

    include "cldop.cdk"
    include "phyinput.inc"
    include "surface.cdk"
    include "nocld.cdk"

    external :: cldoppro_data

    integer, parameter :: IOPTREW = 1
    integer, parameter :: IOPTREI = 2

    real, parameter :: THIRD = 1./3.

    logical, dimension(ni,nkm1) :: nocloud


    real, dimension(:), allocatable :: tausimp, omsimp, gsimp, taulimp, omlimp, glimp

    real, dimension(ni,nkm1) :: aird, rew, rei, rec_cdd, vs1, dp
    real, dimension(ni,nkm1) :: lwpinmp, cldfmp, cldfxp, lwcinmp, iwpinmps, iwcinmps

    real, dimension(:,:,:), allocatable :: iwcinmp, iwpinmp, effradi

    logical :: nostrlwc, readfield_L
    integer :: i, j, k, l, mpcat, istat1, istat2
    real :: rec_grav,  press
    real :: rew1, rew2, rew3, dg, dg2, dg3, tausw, omsw, gsw, tausi, omsi, gsi, y1, y2, y3
    real :: taulw, omlw, glw, tauli, omli, gli
    real :: tauswmp, omswmp, gswmp
    real :: taulwmp, omlwmp, glwmp

    real, pointer, dimension(:,:), contiguous :: zqcmoins, zqimoins, zqnmoins, zqgmoins
    real, pointer, dimension(:,:), contiguous :: zqti1m, zqti2m, zqti3m, zqti4m
    real, pointer, dimension(:,:), contiguous :: zeffradc, zeffradi1 , zeffradi2 , zeffradi3 , zeffradi4
    real, pointer, dimension(:), contiguous :: ztopthw,ztopthi
    real, pointer, dimension(:), contiguous :: zmg,zml
    real, pointer, dimension(:,:), contiguous :: ziwcimp,zlwcimp,zhumoins,ztmoins,zpmoins,zsigw,zfxp,zfmp
    real, pointer, dimension(:), contiguous   :: ztlwp, ztiwp,ztlwpin,ztiwpin
    real, pointer, dimension(:,:), contiguous :: zlwcrad, ziwcrad, zcldrad, zqcplus, zgraupel, &
         zqiplus, zsnow, zqi_cat1, zqi_cat2, zqi_cat3, zqi_cat4

    !----------------------------------------------------------------
    call msg_toall(MSG_DEBUG, 'cldoppro_MP [BEGIN]')

    MKPTR1D(zmg, mg, fbus)
    MKPTR1D(zml, ml, fbus)
    MKPTR1D(ztiwp, tiwp, fbus)
    MKPTR1D(ztiwpin, tiwpin, fbus)
    MKPTR1D(ztlwp, tlwp, fbus)
    MKPTR1D(ztlwpin, tlwpin, fbus)
    MKPTR1D(ztopthi, topthi, fbus)
    MKPTR1D(ztopthw, topthw, fbus)
    MKPTR2Dm1(zcldrad, cldrad, vbus)
    MKPTR2Dm1(zeffradc, effradc, fbus)
    MKPTR2Dm1(zeffradi1, effradi1, fbus)
    MKPTR2Dm1(zeffradi2, effradi2, fbus)
    MKPTR2Dm1(zeffradi3, effradi3, fbus)
    MKPTR2Dm1(zeffradi4, effradi4, fbus)
    MKPTR2Dm1(zfmp, fmp, fbus)
    MKPTR2Dm1(zfxp, fxp, fbus)
    MKPTR2Dm1(zgraupel, qgplus, dbus)
    MKPTR2Dm1(zhumoins, humoins, dbus)
    MKPTR2Dm1(zqti1m, qti1moins, dbus)
    MKPTR2Dm1(zqti2m, qti2moins, dbus)
    MKPTR2Dm1(zqti3m, qti3moins, dbus)
    MKPTR2Dm1(zqti4m, qti4moins, dbus)
    MKPTR2Dm1(ziwcimp, iwcimp, fbus)
    MKPTR2Dm1(ziwcrad, iwcrad, vbus)
    MKPTR2Dm1(zlwcimp, lwcimp, fbus)
    MKPTR2Dm1(zlwcrad, lwcrad, vbus)
    MKPTR2Dm1(zpmoins, pmoins, fbus)
    MKPTR2Dm1(zqcmoins, qcmoins, dbus)
    MKPTR2Dm1(zqcplus, qcplus, dbus)
    MKPTR2Dm1(zqgmoins, qgmoins, dbus)
    MKPTR2Dm1(zqi_cat1, qti1plus, dbus)
    MKPTR2Dm1(zqi_cat2, qti2plus, dbus)
    MKPTR2Dm1(zqi_cat3, qti3plus, dbus)
    MKPTR2Dm1(zqi_cat4, qti4plus, dbus)
    MKPTR2Dm1(zqimoins, qimoins, dbus)
    MKPTR2Dm1(zqiplus, qiplus, dbus)
    MKPTR2Dm1(zqnmoins, qnmoins, dbus)
    MKPTR2Dm1(zsigw, sigw, dbus)
    MKPTR2Dm1(zsnow, qnplus, dbus)
    MKPTR2Dm1(ztmoins, tmoins, dbus)

    call init2nan(aird, rew, rei, rec_cdd, vs1, dp)
    call init2nan(lwpinmp, cldfmp, cldfxp, lwcinmp, iwpinmps, iwcinmps)

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

    ! assume implicit sources have been agregated into zlwc and ziwc in subroutine prep_cw_MP
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
    if (any(dyninread_list_s == 'qc') .or. &
         any(phyinread_list_s(1:phyinread_n) == 'tr/mpqc:p')) then
       readfield_L = .true.
       do k=1,nkm1
          do i=1,ni
             lwcinmp(i,k) = max(zqcplus(i,k), 0.)
          enddo
       enddo
    endif

    IF_MY2: if (stcond(1:6) == 'MP_MY2' .and. &
         (any(dyninread_list_s == 'mpqi') .or. any(phyinread_list_s(1:phyinread_n) == 'tr/mpqi:p')) .and. &
         (any(dyninread_list_s == 'mpqs') .or. any(phyinread_list_s(1:phyinread_n) == 'tr/mpqs:p')) .and. &
         (any(dyninread_list_s == 'mpqg') .or. any(phyinread_list_s(1:phyinread_n) == 'tr/mpqg:p'))) then

       readfield_L = .true.
       !ziwc = zqiplus + zsnow
       do k=1,nkm1
          do i=1,ni
             iwcinmp(i,k,1) = max(zqiplus(i,k), 0.)
             iwcinmp(i,k,2) = max(zsnow(i,k), 0.)
             iwcinmp(i,k,3) = max(zgraupel(i,k), 0.)
          enddo
       enddo

    elseif (stcond=='MP_P3' .and. &
         (any(dyninread_list_s == 'qti1') .or. any(phyinread_list_s(1:phyinread_n) == 'tr/qti1:p'))) then

       readfield_L = .true.
       !ziwc = zqi_cat1
       do k=1,nkm1
          do i=1,ni
             iwcinmp(i,k,1) = max(zqi_cat1(i,k), 0.)
          enddo
       enddo

       IF_NCAT2: if (p3_ncat >= 2 .and. &
            (any(dyninread_list_s == 'qti2') .or. any(phyinread_list_s(1:phyinread_n) == 'tr/qti2:p'))) then

          !ziwc = ziwc + zqi_cat2
          do k=1,nkm1
             do i=1,ni
                iwcinmp(i,k,2) = max(zqi_cat2(i,k), 0.)
             enddo
          enddo

          IF_NCAT3: if (p3_ncat >= 3 .and. &
               (any(dyninread_list_s == 'qti3') .or. any(phyinread_list_s(1:phyinread_n) == 'tr/qti3:p')) ) then

             !ziwc = ziwc + zqi_cat3
             do k=1,nkm1
                do i=1,ni
                   iwcinmp(i,k,3) = max(zqi_cat3(i,k), 0.)
                enddo
             enddo

             IF_NCAT4: if (p3_ncat >= 4 .and. &
                  (any(dyninread_list_s == 'qti4') .or. any(phyinread_list_s(1:phyinread_n) == 'tr/qti4:p'))) then

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

    DO_K: do k=1,nkm1
       do i=1,ni

          ! Normalize water contents to get in-cloud values
!          if (cldfmp(i,k) < cldfth) then ! must revert to a18 value until agregation of implicit amd explicit op properties problem is fixed
          if (cldfmp(i,k) < .01) then
             liqwcin(i,k) = 0.
             icewcin(i,k) = 0.
             cldfmp(i,k) = 0.
          else
             liqwcin(i,k) = liqwcin(i,k)/cldfmp(i,k)
             icewcin(i,k) = icewcin(i,k)/cldfmp(i,k)
          endif
          ! remove ice water content if effective radius >= 1e-4,
          ! since it will not be seen by radiation
          do l=1,mpcat
             if (kount == 0) effradi(i,k,l) = 5.e-5
             if (effradi(i,k,l) >= 1.e-4) then
                iwcinmp(i,k,l) = 0.0
             endif
          enddo

!          if (cldfxp(i,k) < cldfth) then ! must revert to a18 value until agregation of implicit amd explicit op properties problem is fixed
          if (cldfxp(i,k) < .01) then
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
       dp(i,nkm1) = 1. - 0.5*(sig(i,nkm1) + sig(i,nkm1-1))*rec_grav*1000.
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
          icewpin(i,k)  = icewcin(i,k)*dp(i,k)
          liqwpin(i,k)  = liqwcin(i,k)*dp(i,k)
          lwpinmp(i,k)  = lwcinmp(i,k)*dp(i,k)
          iwpinmps(i,k) = 0.
          iwcinmps(i,k) = 0.
       enddo
    enddo

    do k=1,nkm1
       do i=1,ni
          do l=1,mpcat
             iwpinmp(i,k,l) = iwcinmp(i,k,l)*dp(i,k)
             iwpinmps(i,k) = iwpinmps(i,k) + iwpinmp(i,k,l)
             iwcinmps(i,k) = iwcinmps(i,k) + iwcinmp(i,k,l)
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
          nocloud(i,k) = (liqwpin(i,k)+lwpinmp(i,k)+icewpin(i,k)+iwpinmps(i,k) <= wpth)
          ztlwp(i)     = ztlwp(i)   + liqwpin(i,k)*cldfmp(i,k) + lwpinmp(i,k)*cldfxp(i,k)
          ztiwp(i)     = ztiwp(i)   + icewpin(i,k)*cldfmp(i,k) + iwpinmps(i,k)*cldfxp(i,k)
          zlwcrad(i,k) = liqwcin(i,k)*cldfmp(i,k) + lwcinmp(i,k)*cldfxp(i,k)
          ziwcrad(i,k) = icewcin(i,k)*cldfmp(i,k) + iwcinmps(i,k)*cldfxp(i,k)
          cldfrac(i,k) = max(min(cldfmp(i,k)+cldfxp(i,k), 1.0), 0.0)                    ! no overlap=ops
!          cldfrac(i,k) = max(min(max(cldfmp(i,k),cldfxp(i,k)), 1.0), 0.0)                ! max overlap
!          cldfrac(i,k) = min(1., max(0.,  1. - (1.-cldfmp(i,k))*(1.-cldfxp(i,k))))       ! random overlap
          if (nocloud(i,k)) then
             cldfrac(i,k) = 0.0
          endif
          zcldrad(i,k) = cldfrac(i,k)
       end do
    end do

    !     conversion d'unites : tlwp et tiwp en kg/m2
    do i=1,ni
       ztlwp(i) = ztlwp(i) * 0.001
       ztiwp(i) = ztiwp(i) * 0.001
    enddo

    !     initialize output fields
    !
    do i = 1, ni
       ztopthw(i) = 0.0
       ztopthi(i) = 0.0
    end do

    ! begin-code : FOR EFFECTIVE RADII OF IMPLICIT CLOUDS (NON-mp SOURCES)
    !           For calculation of effective radius of water clouds,
    !           set cdd=cloud droplet concentration per cm^3 to 100 over oceans and 500 over land
    !...   ioptrew=1 : as in newrad - from h. barker, based on aircraft data, range 4-17 micron from slingo

    do k = 1, nkm1
       do i = 1, ni
          if (zmg(i) <= 0.5 .and. zml(i) <= 0.5) then
             rec_cdd(i,k) = 0.01
          else
             rec_cdd(i,k) = 0.002
          endif
          aird(i,k) = sig(i,k) * ps(i) / ( tt(i,k) * RGASD )  !aird is air density in kg/m3
          rew(i,k) = (liqwcin(i,k) * aird(i,k) * rec_cdd(i,k))**THIRD
          rew(i,k) = min(max(4., 754.6 * rew(i,k)), 17.0)
       end do
    end do

    !...   choice of effective radius for implicit ice clouds
    !...   ioptrei=1 : from cccma (icewcin must be in g/m3)
    !...   ioptrei=2 : constant in microns

    if (ioptrei == 1) then
       do k = 1, nkm1
          do i = 1, ni
             rei(i,k) = (1000. * icewcin(i,k) * aird(i,k))**0.216
             rei(i,k) = max(min(83.8*rei(i,k), 50.0), 20.0)
          end do
       end do
    elseif (ioptrei == 2) then
       rei = 30.
    endif
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


    DO_NBS: do j = 1, nbs
       do k = 1, nkm1
          do i = 1, ni
             IF_NOCLOUD: if (nocloud(i,k)) then
                taucs(i,k,j) = 0.
                omcs(i,k,j)  = 0.
                gcs(i,k,j)   = 0.
             else
                rew2 = rew(i,k) * rew(i,k)
                rew3 = rew2 * rew(i,k)
                if (liqwpin(i,k) > wpth) then  !implicit
                   tausw = liqwpin(i,k) * &
                        (aws(1,j) + aws(2,j) / rew(i,k) + &
                        aws(3,j) / rew2 + aws(4,j) / rew3)
                   omsw  = 1.0 - (bws(1,j) + bws(2,j) * rew(i,k) + &
                        bws(3,j) * rew2 + bws(4,j) * rew3)
                   gsw   = cws(1,j) + cws(2,j) * rew(i,k) + &
                        cws(3,j) * rew2 + cws(4,j) * rew3
                else
                   tausw = 0.
                   omsw  = 0.
                   gsw   = 0.
                endif


                if (icewpin(i,k) > wpth) then  !implicit
                   dg   = 1.5396 * rei(i,k)
                   dg2  = dg  * dg
                   dg3  = dg2 * dg
                   tausi = icewpin(i,k) * ( ais(1,j) + ais(2,j) / dg )
                   omsi  = 1.0 - (bis(1,j) + bis(2,j) * dg + &
                        bis(3,j) * dg2 + bis(4,j) * dg3)
                   gsi   = cis(1,j) + cis(2,j) * dg + cis(3,j) * dg2 + &
                        cis(4,j) * dg3
                else
                   tausi = 0.
                   omsi  = 0.
                   gsi   = 0.
                endif


                if (lwpinmp(i,k) > wpth) then  !explicit
                   rew1 = zeffradc(i,k) * 1.e+6
                   if (kount == 0) rew1 = 10.    ![microns] assign value at step zero (in case QC is non-zero at initial conditions)
                   rew1 = min(max(4., rew1), 40.0)
                   rew2 = rew1*rew1
                   rew3 = rew1*rew1*rew1
                   tauswmp = lwpinmp(i,k) * &
                        (aws(1,j) + aws(2,j) / rew1 + &
                        aws(3,j) / rew2 + aws(4,j) / rew3)
                   omswmp  = 1.0 - (bws(1,j) + bws(2,j) * rew1 + &
                        bws(3,j) * rew2 + bws(4,j) * rew3)
                   gswmp   = cws(1,j) + cws(2,j) * rew1 + &
                        cws(3,j) * rew2 + cws(4,j) * rew3
                else
                   tauswmp = 0.
                   omswmp  = 0.
                   gswmp   = 0.
                endif

                do l=1,mpcat  !explicit
                   if (iwpinmp(i,k,l) > wpth .and. effradi(i,k,l) < 1.e-4) then
                      dg   = 1.5396 * effradi(i,k,l)*1.e6  ![microns]
                      !if (kount == 0) dg = 1.5396*50.   ![microns] assign value at step zero (in case "QI" is non-zero at initial conditions)
                      dg   = min(max(dg, 15.), 110.)  ! max value is lower otherwise, model is crashing
                      dg2  = dg  * dg
                      dg3  = dg * dg *dg
                      tausimp(l) = iwpinmp(i,k,l) * ( ais(1,j) + ais(2,j) / dg )
                      omsimp(l)  = 1.0 - (bis(1,j) + bis(2,j) * dg + &
                           bis(3,j) * dg2 + bis(4,j) * dg3)
                      gsimp(l)   = cis(1,j) + cis(2,j) * dg + cis(3,j) * dg2 + &
                           cis(4,j) * dg3
                   else
                      tausimp(l) = 0.
                      omsimp(l)  = 0.
                      gsimp(l)   = 0.
                   endif
                enddo

                !PV agregate SW optical properties for liq-imp + ice-imp + liq-exp + ice-exp

                ! old recipe - see cldoppro3
                !              taucs(i,k,j)  = tausw + tausi
                !                y1          = omsw * tausw
                !                y2          = omsi * tausi
                !                omcs(i,k,j) = (y1 + y2) / taucs(i,k,j)
                !                gcs (i,k,j) = (y1 * gsw + y2 * gsi) / (y1 + y2)

                y1 = tausw + tausi
                y2 = tausw*omsw + tausi*omsi
                y3 = tausw*omsw*gsw + tausi*omsi*gsi

                y1 = y1 + tauswmp
                y2 = y2 + tauswmp*omswmp
                y3 = y3 + tauswmp*omswmp*gswmp

                do l=1,mpcat
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
                !
                !     calculate the optical depth for water and ice cloud in visible
                !
                if (j == 1) then
                   ztopthw(i) = ztopthw(i) + tausw + tauswmp
                   ztopthi(i) = ztopthi(i) + tausi
                   do l=1,mpcat
                      ztopthi(i) = ztopthi(i) + tausimp(l)
                   enddo
                endif

             endif IF_NOCLOUD

          enddo
       enddo
    enddo DO_NBS

    DO_NBL: do j = 1, nbl
       do k = 1, nkm1
          do i = 1, ni
             IF_NOCLOUD2: if (nocloud(i,k)) then
                taucl(i,k,j) = 0.
                omcl(i,k,j)  = 0.
                gcl(i,k,j)   = 0.
             else
                rew2 = rew(i,k) * rew(i,k)
                rew3 = rew2 * rew(i,k)

                if (liqwpin(i,k) > wpth) then    !implicit
                   taulw = liqwpin(i,k) * (awl(1,j) + awl(2,j) * rew(i,k)+ &
                        awl(3,j) / rew(i,k) + awl(4,j) / rew2 + &
                        awl(5,j) / rew3)
                   omlw  = 1.0 - (bwl(1,j) + bwl(2,j) / rew(i,k) + &
                        bwl(3,j) * rew(i,k) + bwl(4,j) * rew2)
                   glw   = cwl(1,j) + cwl(2,j) / rew(i,k) + &
                        cwl(3,j) * rew(i,k) + cwl(4,j) * rew2
                else
                   taulw = 0.
                   omlw  = 0.
                   glw   = 0.
                endif

                !----------------------------------------------------------------------
                !     since in fu etc. the param. is for absorptance, so need a factor
                !     icewpin(i,k) / tauli for single scattering albedo
                !----------------------------------------------------------------------

                if (icewpin(i,k) > wpth) then    !implicit
                   dg    = 1.5396 * rei(i,k)
                   dg2   = dg  * dg
                   dg3   = dg2 * dg
                   tauli = icewpin(i,k) * (ail(1,j) + ail(2,j) / dg + &
                        ail(3,j) / dg2)
                   omli  = 1.0 - (bil(1,j) / dg + bil(2,j) + &
                        bil(3,j) * dg + bil(4,j) * dg2) * &
                        icewpin(i,k) / tauli
                   gli   = cil(1,j) + cil(2,j) * dg + cil(3,j) * dg2 + &
                        cil(4,j) * dg3
                else
                   tauli = 0.
                   omli  = 0.
                   gli   = 0.
                endif

                if (lwpinmp(i,k) > wpth) then    !explicit
                   rew1 = zeffradc(i,k) * 1.e6
                   rew1 = min(max(4., rew1), 40.0)  !TEST, JM
                   if (kount == 0) rew1 = 10.    !assign value at step zero (in case QC is non-zero at initial conditions)
                   rew2 = rew1*rew1
                   rew3 = rew1*rew1*rew1
                   taulwmp = lwpinmp(i,k) * (awl(1,j) + awl(2,j) * rew1+ &
                        awl(3,j) / rew1 + awl(4,j) / rew2 + &
                        awl(5,j) / rew3)
                   if (taulwmp < 0.) then
                      print *, '** calc, taulwp: ', taulwmp, rew1, zeffradc(i,k)  !#TODO: avoid printing in MPI/OMP regions, render the listing unusable
                      call physeterror('cldoppro_MP', 'Found negative taulwp values.')
                      return
                   endif

                   omlwmp  = 1.0 - (bwl(1,j) + bwl(2,j) / rew1 + &
                        bwl(3,j) * rew1 + bwl(4,j) * rew2)
                   glwmp   = cwl(1,j) + cwl(2,j) / rew1 + &
                        cwl(3,j) * rew1 + cwl(4,j) * rew2
                else
                   taulwmp = 0.
                   omlwmp  = 0.
                   glwmp   = 0.
                endif

                do l=1,mpcat
                   if (iwpinmp(i,k,l) > wpth .and. effradi(i,k,l) < 1.e-4) then    !explicit
                      dg = 1.5396 * effradi(i,k,l)*1.e6
                      if (kount == 0) dg = 1.5396*50.     !assign value at step zero (in case "QI" is non-zero at initial conditions)
                      dg  = min(max(dg, 15.), 110.)  ! max value is lower to avoid model crashing!
                      dg2 = dg  * dg
                      dg3 = dg2 * dg
                      taulimp(l) = iwpinmp(i,k,l) * (ail(1,j) + ail(2,j) / dg + &
                           ail(3,j) / dg2)
                      omlimp(l) = 1.0 - (bil(1,j) / dg + bil(2,j) + &
                           bil(3,j) * dg + bil(4,j) * dg2) * &
                           iwpinmp(i,k,l) / taulimp(l)
                      glimp(l)  = cil(1,j) + cil(2,j) * dg + cil(3,j) * dg2 + &
                           cil(4,j) * dg3
                   else
                      taulimp(l) = 0.
                      omlimp(l)  = 0.
                      glimp(l)   = 0.
                   endif
                enddo

                !PV agregate LW optical properties for liq-imp + ice-imp + liq-exp + ice-exp

                y1 = taulw + tauli
                y2 = taulw*omlw + tauli*omli
                y3 = taulw*omlw*glw + tauli*omli*gli

                y1 = y1 + taulwmp
                y2 = y2 + taulwmp*omlwmp
                y3 = y3 + taulwmp*omlwmp*glwmp

                do l=1,mpcat
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

