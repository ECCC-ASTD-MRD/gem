!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------

module ccc2_uv_raddriv
   private
   public :: ccc2_uv_raddriv1

contains
   
subroutine ccc2_uv_raddriv1(fatb, fadb, fafb, fctb, fcdb, fcfb, &
     fslo, fsamoon, ps, shtj, sig, &
     tt, o3, o3top, &
     qq, co2, ch4, &
     o2, rmu, r0r, salb, taucs, &
     omcs, gcs, &
     cldfrac, tauae, exta, exoma, exomga, &
     fa, mrk2, &
     il1, il2, ilg, lay, lev)
   use tdpack_const
   use phy_options, only: RAD_NUVBRANDS, rad_atmpath, rad_sun_angle_fix_l
   use ens_perturb, only: ens_nc2d
   implicit none
!!!#include <arch_specific.hf>
#include "nbsnbl.cdk"

   integer ilg,lay,lev,il1,il2

   real, dimension(ilg,RAD_NUVBRANDS)  ::   fatb,fadb,fafb,fctb,fcdb,fcfb

   real ps(ilg), shtj(ilg,lev), sig(ilg,lay), &
        tt(ilg,lay), o3(ilg,lay), &
        o3top(ilg), qq(ilg,lay), rmu(ilg), r0r, salb(ilg,nbs), &
        co2(ilg,lay),ch4(ilg,lay), o2(ilg,lay)

   real taucs(ilg,lay,nbs), omcs(ilg,lay,nbs), gcs(ilg,lay,nbs), &
        cldfrac(ilg,lay), fslo(ilg), fsamoon(ilg)

   real, dimension(ilg, ens_nc2d) :: mrk2


   !@Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !        JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)
   !
   !@Revisions
   !  001    P.Vaillancourt, M.Lazare (sep 2006) : displace a1(i,5)
   !  002    P.Vaillancourt           (Apr 08) : use integer variables(ilg1,ilg2) instead of actual integers
   !  003    P.Vaillancourt           (Feb 12) : assume temperature is isothermal above model top
   !  004    P.Vaillancourt           (Feb 12) : impose min on humidity mixing ratio of 1.5e-6 kg/kg for sw and lw
   !  005    P.Vaillancourt           (Sep 14) : Update to gcm17
   !  006    V.Lee                    (Dec 17) : mcont, mcontg are vectors for MPI bit reproducibility
   !
   !@Object
   !        Main subroutine that executes ccc radiative transfer
   !        for infrared and solar radiation
   !
   !@Arguments
   !              - Output -
   !              - Input -
   ! ps           pressure at ground in unit pa
   ! shtj         sigma at model levels
   ! sig          sigma at model layer center
   ! tfull/tt     temperature at model level / layer center
   ! gt           ground temperature
   ! o3           ozone mass mixing ratio in (g/g)
   ! o3top        accumulated ozone mass above the model top
   ! qq           water vapor specific humidity (mass mixing ratio in
   !              some versions)
   ! rmu          cosine of solar zenith angle
   ! r0r          calculate the variation of solar constant
   ! salb         surface albedo
   ! taucs/taucl  cloud optical depth for solar/longwave
   ! omcs/omcl    cloud single scattering albedo for solar / longwave
   ! gcs/gcl      cloud asymmetry factor for solar/longwave
   ! cldfrac      cloud fraction
   ! tauae        background aerosol optical depth for solar and
   !              longwave
   ! exta         extinction coefficient for solar
   ! exoma        extinction coefficient times single scattering
   !              albedo for solar
   ! exomga       exoma times asymmetry factor for solar
   ! fa           square of asymmetry factor for solar
   ! absa         absorption coefficient for longwave
   ! fslo         solar incoming flux at infrared range (0-2500cm-1)
   ! fsamoon      the energy absorbed between toa and model top level
   ! lcsw         logical key to control call to sw radiative transfer
   ! lclw         logical key to control call to lw radiative transfer
   ! mrk2         Markov chains for stochastic parameter perturbations
   ! il1          1
   ! il2          horizontal dimension (ni)
   ! ilg          horizontal dimension (il2-il1+1)
   ! lay          number of model levels
   ! lev          number of flux levels (lay+1)

#include "ccc_tracegases.cdk"
#include "ccc_aeros.cdk"
#include "tables.cdk"
include "nocld.cdk"

   integer, dimension(ilg) :: mtop
   real, dimension(ilg) :: c1
   real, dimension(ilg) :: c2
   real, dimension(ilg) :: bs
   real, dimension(ilg,lay) :: pg
   real, dimension(ilg,lay) :: qg
   real, dimension(ilg,lay) :: qgs
   real, dimension(ilg,lay) :: pp
   real, dimension(ilg,lay) :: dp
   real, dimension(ilg,lay) :: dps
   real, dimension(ilg,lay) :: taur
   real, dimension(ilg,lay) :: taug
   real, dimension(ilg,lay) :: taua
   real, dimension(ilg,lev) :: pfull
   real, dimension(ilg,lay) :: f1
   real, dimension(ilg,lay) :: f2
   real, dimension(ilg,lay) :: anu
   real, dimension(ilg,lay) :: tauoma
   real, dimension(ilg,lay) :: tauomga
   real, dimension(ilg,lay) :: dip
   real, dimension(ilg,lay) :: dt
   real, dimension(ilg,lay) :: dts
   real, dimension(ilg,2,lev) :: refl
   real, dimension(ilg,2,lev) :: tran
   real, dimension(ilg,lay,5) :: tauae

   !     gathered and other work arrays used generally by solar.

   real, dimension(ilg,12) :: a1
   real, dimension(ilg,12) :: a1g
   real, dimension(ilg,4,lev) :: cumdtr
   real, dimension(ilg,lay,nbs) :: exta
   real, dimension(ilg,lay,nbs) :: exoma
   real, dimension(ilg,lay,nbs) :: exomga
   real, dimension(ilg,lay,nbs) :: fa
   real, dimension(ilg,lay) :: taucsg
   real, dimension(ilg,lay) :: tauomc
   real, dimension(ilg,lay) :: tauomgc
   real, dimension(ilg,lev) :: pfullg
   real, dimension(ilg,lay) :: o3g
   real, dimension(ilg,lay) :: co2g
   real, dimension(ilg,lay) :: ch4g
   real, dimension(ilg,lay) :: o2g
   real, dimension(ilg,lay) :: cldg
   real, dimension(ilg,lay) :: cldmg
   real, dimension(ilg) :: o3topg
   real, dimension(ilg) :: albsur
   real, dimension(ilg) :: rmug, rmu0
   real, dimension(ilg) :: dmix
   integer, dimension(ilg,lay) :: inptg
   integer, dimension(ilg,lay) :: inptmg
   integer, dimension(ilg,lay) :: nblk
   integer, dimension(ilg) :: isun
   integer, dimension(ilg) :: mcont

   !     work arrays used generally by longwave.

   real, dimension(ilg,lay) :: tauci
   real, dimension(ilg,lay) :: omci
   real, dimension(ilg,lay) :: cldm
   real, dimension(ilg,lev) :: bf
   integer, dimension(ilg,lay) :: inpt
   integer, dimension(ilg,lay) :: inptm
   integer, dimension(ilg,lay) :: ncd
   integer, dimension(ilg,lay) :: ncu
   integer, dimension(ilg) :: mcontg !(size is lengath)
   integer, dimension(ilg) :: nct
   integer, dimension(ilg) :: nctg
   integer, dimension(lay) :: ncum
   integer, dimension(lay) :: ncdm

   !     band information.


   real, dimension(nbl) :: sfinptl
   integer, dimension(nbs) :: kgs
   integer, dimension(nbs) :: kgsgh
   integer, dimension(nbl) :: kglgh

   real a11, a12, a13, a21, a22, a23, a31, a32, a33, c20, c30
   real solarc, fracs, x, gw, rgw
   real cut, seuil,specirr,qmr,qmin
   integer i, k, lev1, maxc, jyes, lengath, j, kp1, ig
   logical gh
   integer ilg1,ilg2 !subsize of il1,il2

   integer, parameter :: IB = 1

   parameter (seuil=1.e-3)
   parameter (qmin=1.5e-6)

   !----------------------------------------------------------------------
   !     for hrcoef, 9.80665 / 1004.64 / 100 = 9.761357e-05, in (k / sec),
   !     since we use dp (diff in pressure) instead of diff in meter,
   !     there is a factor 1.02. thus 9.761357e-05 * 1.02 = 9.9565841e-05
   !     uu3 = 3 * u * u, u = 1 / e^0.5
   !----------------------------------------------------------------------

   data specirr /1367.9396/

   !----------------------------------------------------------------------
   !     this code can be extended to about 100 km, if the model top level
   !     is lower than the maximum height, the calculation can be
   !     simplified with less numbers of kgsgh and kglgh accounted
   !----------------------------------------------------------------------

   data kgs   / 6, 4, 6, 4 /
   cut = cldfth ! specified in nocld.cdk, .001 in CCCMA code


   if (std_p_prof(1).lt.1000.0) then
      !   for maximum height about 0.005 hPa
      !        data kgsgh / 3, 4, 4, 9 /
      !        data kglgh / 5, 1, 3, 5, 4, 0, 7, 3, 6 /
      kgsgh(1)=3
      kgsgh(2)=4
      kgsgh(3)=4
      kgsgh(4)=9
      kglgh(1)=5
      kglgh(2)=1
      kglgh(3)=3
      kglgh(4)=5
      kglgh(5)=4
      kglgh(6)=0
      kglgh(7)=7
      kglgh(8)=3
      kglgh(9)=6
   else
      !   if model top level is close to 1 mb
      !        data kgsgh / 3, 3, 3, 6 /
      !        data kglgh / 2, 1, 2, 4, 3, 0, 6, 2, 3 /
      kgsgh(1)=3
      kgsgh(2)=3
      kgsgh(3)=3
      kgsgh(4)=6
      kglgh(1)=2
      kglgh(2)=1
      kglgh(3)=2
      kglgh(4)=4
      kglgh(5)=3
      kglgh(6)=0
      kglgh(7)=6
      kglgh(8)=2
      kglgh(9)=3
   endif

   !----------------------------------------------------------------------
   !     scale mean (annual) value of solar constant by r0r accounting
   !     for eccentricity (passed through common block "eccent" - see
   !     routine sdet2). the spectral irradiance for model is 1366.2035
   !     w/m^2  which is the solar energy contained in the spectral
   !     region 0.2 - 10 um (50000 - 1000 cm).
   !     for longwave, from band1 to band4, the solar and infrared
   !     interaction is considered. the total solar energy considered in
   !     the infrared region is 11.9096 w / m^2. sfinptl is the input
   !     solar flux in each longwave band
   !     the solar input in shortwave region is 1366.2035 - 11.9096 =
   !     1354.3029, the solar fractions for each band are set in gasopts
   !----------------------------------------------------------------------

   solarc                    =  CONSOL2
   fracs                     =  r0r * solarc / specirr
   x                         =  fracs / pi

   sfinptl(1)                =  3.67839 * x
   sfinptl(2)                =  2.79694 * x
   sfinptl(3)                =  3.20284 * x
   sfinptl(4)                =  1.13984 * x
   sfinptl(5)                =  0.31893 * x
   sfinptl(6)                =  0.35404 * x
   sfinptl(7)                =  0.29578 * x
   sfinptl(8)                =  0.99624e-01 * x
   sfinptl(9)                =  0.23220e-01 * x


   !----------------------------------------------------------------------
   !     initialization
   !----------------------------------------------------------------------

   do i = il1, il2
      fsamoon(i)              =  0.0
      fslo(i)                 =  11.9096 * rmu(i) * fracs
      !       shtj(i,lev) = 1. ci-dessous
      pfull(i,lev)            =  0.01 * ps(i) * shtj(i,lev)
      mtop(i)                 =  0
      isun(i)                 =  1
   enddo

   fatb(IL1:IL2,1:RAD_NUVBRANDS) = 0.0
   fadb(IL1:IL2,1:RAD_NUVBRANDS) = 0.0
   fafb(IL1:IL2,1:RAD_NUVBRANDS) = 0.0
   fctb(IL1:IL2,1:RAD_NUVBRANDS) = 0.0
   fcdb(IL1:IL2,1:RAD_NUVBRANDS) = 0.0
   fcfb(IL1:IL2,1:RAD_NUVBRANDS) = 0.0

   do k = 1, lay
      kp1 = k + 1
      do i = il1, il2
         taug(i,k)               =  0.0
         tran(i,1,k)             =  0.0
         tran(i,2,k)             =  0.0
         x                       =  0.01 * ps(i)
         pp(i,k)                 =  sig (i,k) * x
         pfull(i,k)              =  shtj(i,k) * x

         !----------------------------------------------------------------------
         !     specific humidity to mixing ratio.
         !----------------------------------------------------------------------

         qmr                     =  qq(i,k) / (1.0 - qq(i,k))
         qg(i,k)                 =  max(qmr,qmin)
         dp(i,k)                 =  0.0102 * ps(i) * &
              (shtj(i,kp1) - shtj(i,k))
         dt(i,k)                 =  tt(i,k) - 250.0
      enddo
   enddo

   !----------------------------------------------------------------------C
   !     DETERMINATION OF THE HIGHEST PRESSURE LEVEL FOR CONTINUUM        C
   !     CALCULATIONS (> 200 MB). REUSING SPACES OF MTOP AND ISUN.        C
   !----------------------------------------------------------------------C

   mcont = lev
   do k = 1, lev
      do i = il1, il2
         if (pfull(i,k) .ge. 200.) THEN
            mtop(i)               =  mtop(i) + 1
            if (mtop(i) .eq. 1) mcont(i) =  max(k-1,1)
         endif
      enddo
   enddo


   !----------------------------------------------------------------------c
   ! define the spectral sampling for the shortwave and longwave.         c
   ! note that this applies only to mcica, so otherwise we set            c
   ! the output to unity and only one pass done.
   !----------------------------------------------------------------------c


   !----------------------------------------------------------------------
   !     initialize the band-dependant optical property arrays
   !----------------------------------------------------------------------

   !    now done in subroutine aerooppro called by cccmarad

   !----------------------------------------------------------------------
   !     calculate the cloud parameters for swtran and lwtran
   !     reusing inptg, inptmg, tauomgc space
   !----------------------------------------------------------------------

   call ccc_cldifm1 (cldm, tauomgc, anu, a1, ncd, &
        ncu, inptg, nct, ncum, ncdm, &
        cldfrac, pfull, mrk2, lev1, cut, maxc, &
        il1, il2, ilg, lay, lev)


   !----------------------------------------------------------------------
   !     determination of the interpolation points in pressure. inpt for
   !     28 reference levels and inptm for 18 levels
   !     note : remove commented lines at the end of preintp if top is less
   !            than .0005
   !----------------------------------------------------------------------

   call ccc2_preintp(inpt, inptm, dip, a1(1,12), pp, il1, il2, ilg, lay)


      !----------------------------------------------------------------------
      !     determine whether grid points are in daylight. gather the
      !     required field for daylight region
      !----------------------------------------------------------------------

      jyes = 0
      do i = il1, il2
         if (rmu(i) .gt. seuil) then
            jyes                  =  jyes + 1
            isun(jyes)            =  i
         endif
      enddo
      lengath = jyes

      !----------------------------------------------------------------------
      !     skip unnecessary solar
      !----------------------------------------------------------------------

      if (lengath .eq. 0) go to 499

      !     use integer variables instead of actual integers
      ilg1=1
      ilg2=lengath

      ! Set the effecitve solar path length
      select case (rad_atmpath)
      case ('RODGERS67')
         do i=ilg1,ilg2
            j = isun(i)
            rmug(i) =  sqrt (1224.0 * rmu(j) * rmu(j) + 1.0) / 35.0
         enddo
      case ('LI06')
         do i=ilg1,ilg2
            j = isun(i)
            rmug(i) = (2.0 * rmu(j) + sqrt(498.5225 * rmu(j) * rmu(j) + 1.0)) / 24.35
         enddo
      end select
      if (.not.rad_sun_angle_fix_l) then
         rmu0(ilg1:ilg2) = rmug(ilg1:ilg2)
      else
         rmu0(ilg1:ilg2) = rmu(ilg1:ilg2)
      endif
      
      DO230: do i = ilg1, ilg2
         j = isun(i)
         mcontg(i)               =  mcont(j) !mcontg is subset of mcont
         o3topg(i)               =  o3top(j)

         !----------------------------------------------------------------------
         !     c1 and c2 are coefficients for swtran
         !     reusing bf for a factor of anu
         !     reusing dmix for a factor of rmu
         !----------------------------------------------------------------------

         c1(i)                   =  0.75 * rmug(i)
         c2(i)                   =  2.0 * c1(i) * rmug(i)

         a1g(i,1)                =  a1(j,1)
         a1g(i,2)                =  a1(j,2)
         a1g(i,3)                =  a1(j,3)
         a1g(i,4)                =  a1(j,4)
         a1g(i,5)                =  a1(j,5)
         a1g(i,6)                =  a1(j,6)
         a1g(i,7)                =  1.0 - a1g(i,1) - a1g(i,2) - a1g(i,3)
         if (a1g(i,2) .ge. cut) then
            a1g(i,8)              =  a1g(i,4) / a1g(i,2)
         else
            a1g(i,8)              =  0.0
         endif

         a1g(i,9)                =  0.0
         a1g(i,10)               =  0.0
         a1g(i,11)               =  0.0
         x                       =  a1g(i,3) + a1g(i,5) + a1g(i,6)
         if (x .ge. cut) then
            if (a1g(i,1) .ge. cut) then
               a1g(i,9)            =  a1g(i,6) / (x * a1g(i,1))
            endif
            if (a1g(i,2) .ge. cut) then
               a1g(i,10)           =  a1g(i,5) / (x * a1g(i,2))
            endif
            a1g(i,11)             =  a1g(i,3) / x
         endif

         a1g(i,12)               =  a1(j,12)
         nctg(i)                 =  nct(j)
         pfullg(i,lev)           =  pfull(j,lev)
         bf(i,lev)               =  0.0
         dmix(i)                 = (2.0 - rmug(i)) ** 0.40

         !----------------------------------------------------------------------
         !     using a1(i,3) for rmu3
         !----------------------------------------------------------------------

         x                       =  1.0 - rmug(i)
         a1(i,3)                 =  x * x * x
         a1(i,4)                 =  0.0
         !----------------------------------------------------------------------
         !     reusing a1(i,5) for dt0
         !----------------------------------------------------------------------
         ! The following line extrapolates the temperature above model top for moon layer temperature
         !        a1(i,5)                 =  2.0 * tt(j,1) - tt(j,2) - 250.0
         ! The following line assumes an isothermal temperature above model top for moon layer temperature
         a1(i,5)                 =  tt(j,1) - 250.0

      enddo DO230

      DO255: do k = 1, lay
         kp1 = k + 1
         DO250: do i = ilg1, ilg2
            j = isun(i)
            pfullg(i,k)           =  pfull(j,k)

            !----------------------------------------------------------------------
            !     convert from specific humidity to mixing ratio.
            !     reusing omci for dipg
            !----------------------------------------------------------------------

            qgs(i,k)              =  qg(j,k)
            cldmg(i,k)            =  tauomgc(j,k)
            cldg(i,k)             =  cldfrac(j,k)
            nblk(i,k)             =  inptg(j,k)

            o2g(i,k)              =  o2(j,k)
            o3g(i,k)              =  o3(j,k)
            co2g(i,k)             =  co2(j,k)
            ch4g(i,k)             =  ch4(j,k)
            dts(i,k)              =  dt(j,k)
            pg(i,k)               =  pp(j,k)
            omci(i,k)             =  dip(j,k)

            inptg(i,k)            =  inpt(j,k)
            inptmg(i,k)           =  inptm(j,k)

            !----------------------------------------------------------------------
            !     here dp = difp / g = rho * dz, where difp is the layer pressure
            !     difference (in mb), g is the gravity constant, rho is air
            !     density, and dz is layer thickness (in cm). therefore gas mixing
            !     ratio * dp = gas mass * dz. or we can call dp as the air mass
            !     path for a model layer.
            !     0.0102 = 1.02 * 0.01
            !     1mb = 100 pascal = 1000 dynes / cm^2,
            !     1.02 = (1000 dynes / cm^2) / (980 cm / (second^2)).
            !     ps, surface pressure in unit pascal, so with 0.01 factor

            !     reusing bf as a factor for cloud subgrid variability in solar
            !----------------------------------------------------------------------

            dps(i,k)              =  dp(j,k)
            if (cldg(i,k) .lt. cut) then
               bf(i,k)             =  0.0
            else
               bf(i,k)             =  1.0 / (1.0 + 5.68 * anu(j,k) ** 1.4)
            endif
         enddo DO250
      enddo DO255

      !----------------------------------------------------------------------
      !     solar: 4 band for cloud, aerosol, and rayleigh,
      !     20 + 15 (20) monochromatic calculations for gas and radiative
      !     transfer
      !     fatb:  ALL SKY ,DOWNWARD AT THE SURFACE DIR+DIF FLUX, for 6 VIS-UV bands
      !     fadb:  ALL SKY ,DOWNWARD AT THE SURFACE DIRECT FLUX, for 6 VIS-UV bands
      !     fafb:  ALL SKY ,DOWNWARD AT THE SURFACE DIFFUSE FLUX, for 6 VIS-UV bands
      !     fctb:  CLEAR SKY ,DOWNWARD AT THE SURFACE DIR+DIF FLUX, for 6 VIS-UV bands
      !     fcdb:  CLEAR SKY ,DOWNWARD AT THE SURFACE DIRECT FLUX, for 6 VIS-UV bands
      !     fcfb:  CLEAR SKY ,DOWNWARD AT THE SURFACE DIFFUSE FLUX, for 6 VIS-UV bands
      !----------------------------------------------------------------------
         
         do i = ilg1, ilg2
            j = isun(i)
            albsur(i)               =  salb(j,1)
         enddo

         !----------------------------------------------------------------------
         !     scaling aerosol optical properties. taua is aerosol optical depth
         !----------------------------------------------------------------------

         DO310: do k = 1, lay
            do i = ilg1, ilg2
               j = isun(i)
               a11                   =  tauae(j,k,1) * extab(1,1)
               a12                   =  tauae(j,k,2) * extab(1,2)
               a13                   =  tauae(j,k,3) * extab(1,3)
               taua(i,k)             =  a11 + a12 + a13 + &
                    exta(j,k,1) * dps(i,k)

               a21                   =  a11 * omab(1,1)
               a22                   =  a12 * omab(1,2)
               a23                   =  a13 * omab(1,3)
               tauoma(i,k)           =  a21 + a22 + a23 + &
                    exoma(j,k,1) * dps(i,k)

               a31                   =  a21 * gab(1,1)
               a32                   =  a22 * gab(1,2)
               a33                   =  a23 * gab(1,3)
               tauomga(i,k)          =  a31 + a32 + a33 + &
                    exomga(j,k,1) * dps(i,k)

               f1(i,k)               =  a31 * gab(1,1) + a32 * gab(1,2) + &
                    a33 * gab(1,3) + fa(j,k,1)

               !----------------------------------------------------------------------
               !     scaling the cloud optical properties due to subgrid variability
               !     and standard scaling for radiative transfer
               !----------------------------------------------------------------------

               if (cldg(i,k) .ge. cut) then
                  if (k .eq. 1) then
                     tauci(i,k)        =  taucs(j,k,1)
                     x                 =  taucs(j,k,1) + &
                          9.2 * sqrt(taucs(j,k,1))
                  else
                     tauci(i,k)        =  tauci(i,k-1) + taucs(j,k,1)
                     x                 =  taucs(j,k,1) + &
                          9.2 * sqrt(tauci(i,k-1))
                  endif

                  taucsg(i,k)         =  taucs(j,k,1) / (1.0 + 0.185 * &
                       x * dmix(i) * bf(i,k))

                  c20                 =  taucsg(i,k) * omcs(j,k,1)
                  tauomc(i,k)         =  tauoma(i,k) + c20

                  c30                 =  c20 * gcs(j,k,1)
                  tauomgc(i,k)        =  tauomga(i,k) + c30
                  f2(i,k)             =  f1(i,k) + c30 * gcs(j,k,1)
               else
                  tauci(i,k)          =  0.0
                  taucsg(i,k)         =  0.0
                  tauomc(i,k)         =  0.0
                  tauomgc(i,k)        =  0.0
                  f2(i,k)             =  0.0
               endif
            enddo
         enddo DO310

         !----------------------------------------------------------------------
         !     raylei, near-ir rayleigh scattering, it is independent of ig.
         !     reusing a1(i,1) for moon layer attenuation
         !----------------------------------------------------------------------

         gh = .false.

          ! simplified version of do 440 loop for uvindex flux calculations
           DO400b: do ig = 3, kgs(1)

               !----------------------------------------------------------------------
               !     raylev, visible rayleigh scattering, it is dependant on ig.
               !----------------------------------------------------------------------

               call ccc2_raylev2(taur, ig, dps, a1(1,3), ilg1, ilg2, ilg, lay)

               !----------------------------------------------------------------------
               !     solar attenuation above the model top lay. only apply to band
               !     one for o3 and o2. this is true only for model top level above
               !     about 1 mb, water vapor contribution is small.
               !----------------------------------------------------------------------

               call ccc2_sattenu4(a1, IB, ig, rmug, o3topg, co2g, ch4g, o2g, &
                    pfullg, a1g(1,12), dts, a1(1,5), inptg, gh, &
                    ilg1, ilg2, ilg)

            !----------------------------------------------------------------------
            !     downward flux above 1 mb, further flux attenuation factor for
            !     the lower region
            !----------------------------------------------------------------------

            if (lev1 .gt. 1) then
               call ccc2_strandn3(tran, bs, a1, rmug, dps, dts, o3g, o2g, a1(1,3), &
                    IB, ig, lev1, ilg1, ilg2, ilg, lay, lev)

            else
               do i = ilg1, ilg2
                  bs(i)             =  a1(i,1)
               enddo
            endif

            call ccc2_gasopts5(taug, gw, dps, IB, ig, o3g, qgs, co2g, ch4g, o2g, &
                 inptmg, mcontg, omci, dts, a1(1,3), lev1, gh, &
                 ilg1, ilg2, ilg, lay)


            call ccc_swtran (refl, tran, cumdtr, bs, taua, &
                 taur, taug, tauoma, tauomga, f1, &
                 f2, taucsg, tauomc, tauomgc, cldg, &
                 cldmg, a1g, rmug, c1, c2, &
                 albsur, nblk, nctg, cut, lev1, &
                 ilg1, ilg2, ilg, lay, lev)

            if (lev1 .gt. 1) then
               call ccc2_stranup3(refl, dps, dts, o3g, o2g, IB, ig, lev1, &
                    ilg1, ilg2, ilg, lay, lev)
            endif


            !----------------------------------------------------------------------
            !     gather back the required fields
            !----------------------------------------------------------------------

            rgw = gw * fracs
            DO350b: do i = ilg1, ilg2
               j = isun(i)
               x                   =  a1g(i,7) * cumdtr(i,1,lev) + &
                    a1g(i,1) * cumdtr(i,2,lev) + &
                    a1g(i,2) * cumdtr(i,3,lev) + &
                    a1g(i,3) * cumdtr(i,4,lev)
               a1(i,2)             =  rgw * rmu0(i)
               !PV fluxes in VIS_UV sub-bands
                  fatb(J,IG)          =   TRAN(I,2,LEV)* A1(I,2)
                  fadb(J,IG)          =   X * BS(I) * A1(I,2)
                  fafb(J,IG)          =   fatb(J,IG) - fadb(J,IG)
                  fctb(J,IG)          =   TRAN(I,1,LEV)* A1(I,2)
                  fcdb(J,IG)          =   CUMDTR(I,1,LEV) * BS(I) * A1(I,2)
                  fcfb(J,IG)          =   fctb(J,IG) - fcdb(J,IG)
            enddo DO350b

           enddo DO400b

   499   continue
   return
end subroutine ccc2_uv_raddriv1

end module ccc2_uv_raddriv
