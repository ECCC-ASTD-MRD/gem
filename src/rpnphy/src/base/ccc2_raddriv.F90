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


subroutine ccc2_raddriv3(fsg, fsd, fsf, fsv, fsi, &
     fatb,fadb,fafb,fctb,fcdb,fcfb, &
     albpla, fdl, ful, hrs, hrl, &
     cst, csb, clt, clb, par, &
     flxds,flxus,flxdl,flxul, &
     fslo, fsamoon, ps, shtj, sig, &
     tfull, tt, gt, o3, o3top, &
     qq, co2, ch4, an2o, f11, &
     f12, f113, f114,o2,rmu, r0r, salb, em0, taucs, &
     omcs, gcs, taucl, omcl, gcl, &
     cldfrac, tauae, exta, exoma, exomga, &
     fa, absa, lcsw, lclw, mrk2, luvonly, &
     il1, il2, ilg, lay, lev)
   use tdpack_const
   use phy_options, only: RAD_NUVBRANDS, rad_atmpath, rad_sun_angle_fix_l
   use ens_perturb, only: ens_nc2d
   implicit none
!!!#include <arch_specific.hf>
#include "nbsnbl.cdk"

   integer ilg,lay,lev,il1,il2
   real fsg(ilg), fsd(ilg), fsf(ilg), fsv(ilg), fsi(ilg), &
        albpla(ilg), fdl(ilg), ful(ilg), hrs(ilg,lay), hrl(ilg,lay), &
        cst(ilg), csb(ilg), clt(ilg), clb(ilg), par(ilg), em0(ilg)

   real, dimension(ilg,RAD_NUVBRANDS)  ::   fatb,fadb,fafb,fctb,fcdb,fcfb

   real ps(ilg), shtj(ilg,lev), sig(ilg,lay), &
        tfull(ilg,lev), tt(ilg,lay), gt(ilg), o3(ilg,lay), &
        o3top(ilg), qq(ilg,lay), rmu(ilg), r0r, salb(ilg,nbs), &
        co2(ilg,lay),ch4(ilg,lay), an2o(ilg,lay), f11(ilg,lay),f12(ilg,lay), &
        f113(ilg,lay), f114(ilg,lay),o2(ilg,lay)

   real taucs(ilg,lay,nbs), omcs(ilg,lay,nbs), gcs(ilg,lay,nbs), &
        taucl(ilg,lay,nbl), omcl(ilg,lay,nbl), gcl(ilg,lay,nbl), &
        cldfrac(ilg,lay), fslo(ilg), fsamoon(ilg)

   real, dimension(ilg, ens_nc2d) :: mrk2
   logical lcsw, lclw, luvonly

   real flxds(ilg,lev),flxus(ilg,lev),flxdl(ilg,lev),flxul(ilg,lev)

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
   ! fsg          downward flux absorbed by ground.
   ! fsd          direct downward flux at the surface.
   ! fsf          diffuse downward flux at the surface.
   ! fsv          visible downward flux at the surface.
   ! fsi          near infrared downward flux at the surface.
   ! albpla       planetary albedo.
   ! ful/fdl      upward lw flux at the top / surface
   ! hrs/hrl      solar heating rate / longwave cooling rate
   ! cst/csb      net clear sky solar flux at top / surface
   ! clt/clb      net clear sky longwave flux at top/surface
   ! par          photosynthetic active radiation.
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
   real, dimension(ilg,lev) :: flxu
   real, dimension(ilg,lev) :: flxd
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
   real, dimension(ilg,lay) :: urbf
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

   real, dimension(ilg,lay,nbl) :: absa
   real, dimension(ilg,lay) :: tauci
   real, dimension(ilg,lay) :: omci
   real, dimension(ilg,lay) :: gci
   real, dimension(ilg,lay) :: cldm
   real, dimension(ilg,lev) :: bf
   integer, dimension(ilg,lay) :: inpt
   integer, dimension(ilg,lay) :: inptm
   integer, dimension(ilg,lay) :: inpr
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
   integer, dimension(nbl) :: kgl
   integer, dimension(nbl) :: kglgh

   real a11, a12, a13, a21, a22, a23, a31, a32, a33, c20, c30, tran0
   real solarc, fracs, x, gw, rgw, dfnet, gwgh, rsolarc, pgw
   real ubeta0, epsd0, hrcoef, uu3, cut, seuil,specirr,qmr,qmin
   integer i, k, ib, lev1, maxc, jyes, lengath, j, kp1, ig
   logical gh
   integer ilg1,ilg2 !subsize of il1,il2

   parameter (seuil=1.e-3)
   parameter (qmin=1.5e-6)

   !----------------------------------------------------------------------
   !     for hrcoef, 9.80665 / 1004.64 / 100 = 9.761357e-05, in (k / sec),
   !     since we use dp (diff in pressure) instead of diff in meter,
   !     there is a factor 1.02. thus 9.761357e-05 * 1.02 = 9.9565841e-05
   !     uu3 = 3 * u * u, u = 1 / e^0.5
   !----------------------------------------------------------------------

   data hrcoef, uu3 / 9.9565841e-05, 1.1036383 /
   !old      data specirr /1366.2035/
   data specirr /1367.9396/

   !----------------------------------------------------------------------
   !     this code can be extended to about 100 km, if the model top level
   !     is lower than the maximum height, the calculation can be
   !     simplified with less numbers of kgsgh and kglgh accounted
   !----------------------------------------------------------------------

   data kgs   / 6, 4, 6, 4 /
   data kgl   / 1, 1, 2, 5, 2, 3, 3, 6, 4 /

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


   !print*,'ps',ps
   !print*,'shtj',shtj
   !print*,'sig',sig
   !print*,'tfull',tfull
   !print*,'tt',tt
   !print*,'gt',gt
   !print*,'o3',o3
   !print*,'o3top',o3top
   !print*,'qq',qq
   !print*,'rmu',rmu
   !print*,'r0r',r0r
   !print*,'salb',salb
   !!print*,'taucs',taucs
   !!print*,'taucl',taucl
   !!print*,'omcs',omcs
   !!print*,'omcl',omcl
   !!print*,'gcs',gcs
   !!print*,'gcl',gcl
   !print*,'cldfrac',cldfrac
   !!print*,'tauae',tauae
   !!print*,'exta',exta
   !!print*,'exomga',exomga
   !!print*,'absa',absa


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
      fsg(i)                  =  0.0
      fsd(i)                  =  0.0
      fsf(i)                  =  0.0
      fsi(i)                  =  0.0
      fsv(i)                  =  0.0
      cst(i)                  =  0.0
      csb(i)                  =  0.0
      par(i)                  =  0.0
      fsamoon(i)              =  0.0
      fslo(i)                 =  11.9096 * rmu(i) * fracs
      albpla(i)               =  0.0
      !       shtj(i,lev) = 1. ci-dessous
      pfull(i,lev)            =  0.01 * ps(i) * shtj(i,lev)
      flxds(i,lev)            =  0.0
      flxus(i,lev)            =  0.0
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
         hrs(i,k)                =  0.0
         hrl(i,k)                =  0.0
         x                       =  0.01 * ps(i)
         pp(i,k)                 =  sig (i,k) * x
         pfull(i,k)              =  shtj(i,k) * x
         flxds(i,k)              =  0.0
         flxus(i,k)              =  0.0

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

   !new      if(mcica.ne.0) then
   !new        call initspecsampl(nsample_sw, nsample_lw, nbs, nbl,
   !new     1                     maxng, ivers)
   !new      else
   !new        do igh = 1, 2
   !new        do ig = 1, maxng
   !new          do ib = 1, nbs
   !new            nsample_sw(ib,ig,igh) = 1
   !new          end do ! ib
   !new          do ib = 1, nbl
   !new            nsample_lw(ib,ig,igh) = 1
   !new          end do ! ib
   !new        end do ! ig
   !new        end do ! igh
   !new      endif


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

   if (lcsw) then

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
         flxu(i,lev)             =  0.0
         flxd(i,lev)             =  0.0
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
            flxu(i,k)             =  0.0
            flxd(i,k)             =  0.0
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

      !     flxu:   all sky sw upward flux.
      !     flxd:   all sky sw downward flux.
      !     fsg:    downward flux absorbed by ground.
      !     fsd:    direct downward flux at the surface.
      !     fsf:    diffuse downward flux at the surface.
      !     fsv:    visible downward flux at the surface.
      !     fsi:    near infrared downward flux at the surface.
      !     par:    photosynthetic active radiation.
      !     albpla: planetary albedo.
      !     cst:    net clear sky flux at top.
      !     csb:    net clear sky flux at surface.
      !     fatb:  ALL SKY ,DOWNWARD AT THE SURFACE DIR+DIF FLUX, for 6 VIS-UV bands
      !     fadb:  ALL SKY ,DOWNWARD AT THE SURFACE DIRECT FLUX, for 6 VIS-UV bands
      !     fafb:  ALL SKY ,DOWNWARD AT THE SURFACE DIFFUSE FLUX, for 6 VIS-UV bands
      !     fctb:  CLEAR SKY ,DOWNWARD AT THE SURFACE DIR+DIF FLUX, for 6 VIS-UV bands
      !     fcdb:  CLEAR SKY ,DOWNWARD AT THE SURFACE DIRECT FLUX, for 6 VIS-UV bands
      !     fcfb:  CLEAR SKY ,DOWNWARD AT THE SURFACE DIFFUSE FLUX, for 6 VIS-UV bands
      !----------------------------------------------------------------------

      DONBS: do ib = 1, nbs

         do i = ilg1, ilg2
            j = isun(i)
            albsur(i)               =  salb(j,ib)
         enddo

         !----------------------------------------------------------------------
         !     scaling aerosol optical properties. taua is aerosol optical depth
         !----------------------------------------------------------------------

         DO310: do k = 1, lay
            do i = ilg1, ilg2
               j = isun(i)
               a11                   =  tauae(j,k,1) * extab(ib,1)
               a12                   =  tauae(j,k,2) * extab(ib,2)
               a13                   =  tauae(j,k,3) * extab(ib,3)
               taua(i,k)             =  a11 + a12 + a13 + &
                    exta(j,k,ib) * dps(i,k)

               a21                   =  a11 * omab(ib,1)
               a22                   =  a12 * omab(ib,2)
               a23                   =  a13 * omab(ib,3)
               tauoma(i,k)           =  a21 + a22 + a23 + &
                    exoma(j,k,ib) * dps(i,k)

               a31                   =  a21 * gab(ib,1)
               a32                   =  a22 * gab(ib,2)
               a33                   =  a23 * gab(ib,3)
               tauomga(i,k)          =  a31 + a32 + a33 + &
                    exomga(j,k,ib) * dps(i,k)

               f1(i,k)               =  a31 * gab(ib,1) + a32 * gab(ib,2) + &
                    a33 * gab(ib,3) + fa(j,k,ib)

               !----------------------------------------------------------------------
               !     scaling the cloud optical properties due to subgrid variability
               !     and standard scaling for radiative transfer
               !----------------------------------------------------------------------

               if (cldg(i,k) .ge. cut) then
                  if (k .eq. 1) then
                     tauci(i,k)        =  taucs(j,k,ib)
                     x                 =  taucs(j,k,ib) + &
                          9.2 * sqrt(taucs(j,k,ib))
                  else
                     tauci(i,k)        =  tauci(i,k-1) + taucs(j,k,ib)
                     x                 =  taucs(j,k,ib) + &
                          9.2 * sqrt(tauci(i,k-1))
                  endif

                  taucsg(i,k)         =  taucs(j,k,ib) / (1.0 + 0.185 * &
                       x * dmix(i) * bf(i,k))

                  c20                 =  taucsg(i,k) * omcs(j,k,ib)
                  tauomc(i,k)         =  tauoma(i,k) + c20

                  c30                 =  c20 * gcs(j,k,ib)
                  tauomgc(i,k)        =  tauomga(i,k) + c30
                  f2(i,k)             =  f1(i,k) + c30 * gcs(j,k,ib)
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

         if (ib .ne. 1) then
            call ccc_raylei (taur, ib, dps, ilg1, ilg2, ilg, lay)
         endif

         gh = .false.

         DOMAJSW: do ig = 1, kgs(ib)

            if (ib .eq. 1) then

               !----------------------------------------------------------------------
               !     raylev, visible rayleigh scattering, it is dependant on ig.
               !----------------------------------------------------------------------

               call ccc2_raylev2(taur, ig, dps, a1(1,3), ilg1, ilg2, ilg, lay)

               !----------------------------------------------------------------------
               !     solar attenuation above the model top lay. only apply to band
               !     one for o3 and o2. this is true only for model top level above
               !     about 1 mb, water vapor contribution is small.
               !----------------------------------------------------------------------

               call ccc2_sattenu4(a1, ib, ig, rmug, o3topg, co2g, ch4g, o2g, &
                    pfullg, a1g(1,12), dts, a1(1,5), inptg, gh, &
                    ilg1, ilg2, ilg)
            else
               do i = ilg1, ilg2
                  a1(i,1)           =  1.0
               enddo
            endif

            !----------------------------------------------------------------------
            !     downward flux above 1 mb, further flux attenuation factor for
            !     the lower region
            !----------------------------------------------------------------------

            if (lev1 .gt. 1) then
               call ccc2_strandn3(tran, bs, a1, rmug, dps, dts, o3g, o2g, a1(1,3), &
                    ib, ig, lev1, ilg1, ilg2, ilg, lay, lev)

            else
               do i = ilg1, ilg2
                  bs(i)             =  a1(i,1)
               enddo
            endif

            call ccc2_gasopts5(taug, gw, dps, ib, ig, o3g, qgs, co2g, ch4g, o2g, &
                 inptmg, mcontg, omci, dts, a1(1,3), lev1, gh, &
                 ilg1, ilg2, ilg, lay)


            call ccc_swtran (refl, tran, cumdtr, bs, taua, &
                 taur, taug, tauoma, tauomga, f1, &
                 f2, taucsg, tauomc, tauomgc, cldg, &
                 cldmg, a1g, rmug, c1, c2, &
                 albsur, nblk, nctg, cut, lev1, &
                 ilg1, ilg2, ilg, lay, lev)

            !new          if(mcica.eq.0) then
            !new            call swtran2(refl, tran, cumdtr, bs, taua, taur, taug,
            !new     1                   tauoma, tauomga, f1, f2, taucsg,
            !new     2                   tauomc, tauomgc, cldg, cldmg, a1g,
            !new     3                   rmug, c1, c2, albsur, csalg, nblkg, nctg,
            !new     4                   cut, lev1, ilg1, ilg2, ilg, lay, lev)
            !new
            if (lev1 .gt. 1) then
               call ccc2_stranup3(refl, dps, dts, o3g, o2g, ib, ig, lev1, &
                    ilg1, ilg2, ilg, lay, lev)
            endif
            !new          else
            !new            call swtran_mcica(refl, tran, cumdtr, bs, taua, taur, taug,
            !new     1                        tauoma, tauomga, f1, f2, taucsg,
            !new     2                        tauomc, tauomgc, cldg,
            !new     3                        rmug, c1, c2, albsur, csalg, nctg,
            !new     4                        cut, lev1, ilg1, ilg2, ilg, lay, lev)
            !c
            !new            if (lev1 .gt. 1) then
            !new             call stranup3(refl, dps, dts, o3g, o2g, ib, ig, lev1,
            !new     1                     ilg1, ilg2, ilg, lay, lev)
            !new            endif

            !           * for the total sky fluxes weight the clear and cloudy sky
            !           * fluxes by the total vertically projected cloud fraction.

            !new            do k = 1,lev
            !new            do i = ilg1, ilg2
            !new              j = isun(i)
            !new              if (cldt(j).lt.1.0) then
            !new                 refl(i,2,k) = (1.0 - cldt(j)) * refl(i,1,k) +
            !new     1                         cldt(j)  * refl(i,2,k)
            !new                 tran(i,2,k) = (1.0 - cldt(j)) * tran(i,1,k) +
            !new     1                         cldt(j)  * tran(i,2,k)
            !new              end if
            !new            end do ! i
            !new            end do ! k
            !new          endif


            !----------------------------------------------------------------------
            !     gather back the required fields
            !----------------------------------------------------------------------

            rgw = gw * fracs
            DO350: do i = ilg1, ilg2
               j = isun(i)
               x                   =  a1g(i,7) * cumdtr(i,1,lev) + &
                    a1g(i,1) * cumdtr(i,2,lev) + &
                    a1g(i,2) * cumdtr(i,3,lev) + &
                    a1g(i,3) * cumdtr(i,4,lev)
               a1(i,2)             =  rgw * rmu0(i)
               fsd(j)              =  fsd(j) + x * bs(i) * a1(i,2)
               cst(j)              =  cst(j) + (1.0 - refl(i,1,1) * &
                    a1(i,1)) * a1(i,2)
               csb(j)              =  csb(j) + (tran(i,1,lev) - &
                    refl(i,1,lev)) * a1(i,2)

               flxu(i,1)           =  flxu(i,1) + refl(i,2,1) * a1(i,2)
               flxd(i,1)           =  flxd(i,1) + tran(i,2,1) * a1(i,2)

               !PV fluxes in VIS_UV sub-bands
               if (ib .eq. 1) then
                  fatb(J,IG)          =   TRAN(I,2,LEV)* A1(I,2)
                  fadb(J,IG)          =   X * BS(I) * A1(I,2)
                  fafb(J,IG)          =   fatb(J,IG) - fadb(J,IG)
                  fctb(J,IG)          =   TRAN(I,1,LEV)* A1(I,2)
                  fcdb(J,IG)          =   CUMDTR(I,1,LEV) * BS(I) * A1(I,2)
                  fcfb(J,IG)          =   fctb(J,IG) - fcdb(J,IG)
               endif

            enddo DO350

            !new          rgw = gw * www * (fracs)
            !new          do 350 i = ilg1, ilg2
            !new            j = isun(i)
            !new            if(mcica.eq.0) then
            !new              x                 =  a1g(i,7) * cumdtr(i,1,lev) +
            !new     1                             a1g(i,1) * cumdtr(i,2,lev) +
            !new     2                             a1g(i,2) * cumdtr(i,3,lev) +
            !new     3                             a1g(i,3) * cumdtr(i,4,lev)
            !new            else
            !new              x                 =  (1.0 - cldt(j)) * cumdtr(i,1,lev) +
            !new     1                             cldt(j) * cumdtr(i,2,lev)
            !new            endif
            !new            a1(i,2)             =  rgw * rmug(i)
            !new            fsd(j)              =  fsd(j) + x * bs(i) * a1(i,2)
            !new            cst(j)              =  cst(j) + (1.0 - refl(i,1,1) *
            !new     1                             a1(i,1)) * a1(i,2)
            !new            csb(j)              =  csb(j) + (tran(i,1,lev) -
            !new     1                             refl(i,1,lev)) * a1(i,2)
            !new
            !new            csd(j)              =  csd(j) +
            !new     1                             cumdtr(i,1,lev) * bs(i) * a1(i,2)
            !new            csf(j)              =  csf(j) + tran(i,1,lev) * a1(i,2)
            !new            flxu(i,1)           =  flxu(i,1) + refl(i,2,1) * a1(i,2)
            !new            flxd(i,1)           =  flxd(i,1) + tran(i,2,1) * a1(i,2)

            !new  350     continue


            !----------------------------------------------------------------------
            !     heating rate calculation, for stability in calculation, each ig
            !     is done separately. heating rate in (k / sec),
            !----------------------------------------------------------------------

            do k = 1, lay
               kp1 = k + 1
               do i = ilg1, ilg2
                  j = isun(i)
                  dfnet             = (tran(i,2,k) - tran(i,2,kp1) - &
                       refl(i,2,k) + refl(i,2,kp1)) * &
                       a1(i,2)
                  hrs(j,k)        =  hrs(j,k) + hrcoef * max(dfnet,0.) / dps(i,k)


                  flxu(i,kp1)       =  flxu(i,kp1) + refl(i,2,kp1) * a1(i,2)
                  flxd(i,kp1)       =  flxd(i,kp1) + tran(i,2,kp1) * a1(i,2)
               enddo
            enddo

            !----------------------------------------------------------------------
            !     fsamoon is the energy absorbed between toa and model top level.
            !     a1(i,4) is the adjustment for upward flux from model top level
            !     to toa used for planetary albedo
            !----------------------------------------------------------------------

            if (ib .eq. 1) then
               do i = ilg1, ilg2
                  j = isun(i)
                  x                 = (1.0 - a1(i,1)) * a1(i,2)
                  fsamoon(j)        =  fsamoon(j) + x * (1.0 + refl(i,2,1))
                  a1(i,4)           =  a1(i,4) - x * refl(i,2,1)
               enddo
            endif

            if (ib .eq. 1 .and. ig .eq. 2) then
               do i = ilg1, ilg2
                  par(isun(i))      =   flxd(i,lev)
               enddo
            endif

         enddo DOMAJSW

        if (ib == 1 .and. luvonly) return

         !----------------------------------------------------------------------
         !     in accumulated space with interval close to 1, the extinction
         !     coefficients is extremely large, the calculation process can be
         !     simplified by ignoring scattering, reflection, cloud and aerosol.
         !----------------------------------------------------------------------

         gh = .true.

         DOMINSW: do ig = 1, kgsgh(ib)

            call ccc2_sattenu4(a1, ib, ig, rmug, o3topg, co2g, ch4g, o2g, &
                 pfullg, a1g(1,12), dts, a1(1,5), inptg, gh, &
                 ilg1, ilg2, ilg)

            call ccc2_strandngh4(tran, gwgh, a1, taua, tauoma, taucsg, tauomc, &
                 cldg, rmug, dps, o3g, qgs, co2g, ch4g, o2g, ib, &
                 ig, inptg, omci, dts, lev1, gh, cut, &
                 ilg1, ilg2, ilg, lay, lev)

            !new          if(mcica.ne.0) then
            !new
            !new           * for the total sky fluxes weight the clear and cloudy sky
            !new           * fluxes by the total vertically projected cloud fraction.
            !new
            !new            do k = 1, lev
            !new            do i = ilg1, ilg2
            !new              j = isun(i)
            !new              if (cldt(j).lt.1.0) then
            !new                tran(i,2,k) = (1.0 - cldt(j)) * tran(i,1,k) +
            !new     1                        cldt(j)  * tran(i,2,k)
            !new              end if
            !new            end do ! i
            !new            end do ! k
            !new          endif


            rgw = gwgh * fracs
            !new          rgw = gwgh * www* fracs

            do i = ilg1, ilg2
               j = isun(i)
               a1(i,2)             =  rgw * rmu0(i)
               cst(j)              =  cst(j) + a1(i,2)
               csb(j)              =  csb(j) + tran(i,1,lev) * a1(i,2)

               fsamoon(j)          =  fsamoon(j) + &
                    a1(i,2) * (1.0 - tran(i,2,1))
               flxd(i,1)           =  flxd(i,1) + a1(i,2) * tran(i,2,1)
            enddo
            do k = 1, lay
               kp1 = k + 1
               do i = ilg1, ilg2
                  j = isun(i)
                  flxd(i,kp1)       =  flxd(i,kp1) + a1(i,2) * tran(i,2,kp1)
                  dfnet             =  tran(i,2,k) - tran(i,2,kp1)
                  hrs(j,k)          =  hrs(j,k) + hrcoef * a1(i,2) * &
                       max(dfnet,0.) / dps(i,k)

               enddo
            enddo
         enddo DOMINSW

         if (ib .eq. 1) then
            do i = ilg1, ilg2
               fsv(isun(i))        =  flxd(i,lev)
            enddo
         endif

      enddo DONBS

      !----------------------------------------------------------------------
      !     gather back required field. for planetary albedo the incoming
      !     energy of 11.9096 * fracs is totally absorbed in longwave part
      !----------------------------------------------------------------------

      rsolarc = r0r * solarc
      do i = ilg1, ilg2
         j = isun(i)
         fsg(j)                  =  flxd(i,lev) - flxu(i,lev)
         fsi(j)                  =  flxd(i,lev) - fsv(j)
         fsf(j)                  =  flxd(i,lev) - fsd(j)

         !        cst(j)                  =  cst(j) + fslo(j)
         albpla(j)               = (flxu(i,1) + a1(i,4)) / &
              (rsolarc * rmu0(i))
      enddo
      !     on veut les flux en sortie
      !     make sure that sw heating rate is never negative
      do k = 1,lev
         do i = ilg1, ilg2
            j = isun(i)
            flxds(j,k)=flxd(i,k)
            flxus(j,k)=flxu(i,k)
         enddo
      enddo
      do k = 1,lay
         do i = ilg1, ilg2
            j = isun(i)
            hrs(j,k)=max(hrs(j,k),0.)
         enddo
      enddo

499   continue

      if (luvonly) return

   endif
   !     (lcsw)

   !----------------------------------------------------------------------
   !     longwave: 9 band for cloud, aerosol, continuum, and planck.
   !     24+22 monochromatic calculations for gas and radiative transfer

   !     flxu: all sky lw upward flux.
   !     flxd: all sky lw downward flux.
   !     ful:  upward lw flux at the top.
   !     fdl:  down lw flux received at the ground.
   !     clt:  net clear sky upward flux at the top.
   !     clb:  net clear sky downward flux at the surface.
   !----------------------------------------------------------------------

   if (lclw) then

      !----------------------------------------------------------------------
      !     convert from specific humidity to mixing ratio (bounded) and
      !     bound temperature for planck calculation.
      !----------------------------------------------------------------------

      do i = il1, il2
         ! The following line extrapolates the temperature above model top for moon layer temperature
         !        a1(i,5)                 =  2.0 * tt(i,1) - tt(i,2) - 250.0
         ! The following line assumes an isothermal temperature above model top for moon layer temperature
         a1(i,5)                 =  tt(i,1) - 250.0
         clt(i)                  =  0.0
         clb(i)                  =  0.0
      enddo
      do k = 1, lev
         do i = il1, il2
            flxu(i,k)               =  0.0
            flxd(i,k)               =  0.0
         enddo
      enddo
      !----------------------------------------------------------------------
      !     determination of the interpolation points in the ratio of co2
      !     to water vapor for tlinehc. reuse the space of pg for dir
      !     and reuse tauomc as a work array
      !----------------------------------------------------------------------

      call ccc2_preintr3 (inpr, pg, qg, co2, tauomc, il1, il2, ilg, lay)

      DONBL: do ib = 1, nbl

         !----------------------------------------------------------------------
         !     using c1 space for slwf which is the input solar energy in the
         !     infrared region. total 11.9096 w / m^2 from standard
         !     calculation
         !     scaling cloud optical properties for ir scattering calculation
         !----------------------------------------------------------------------

         do i = il1, il2
            if (rmu(i) .gt. 0.0) then
               c1(i)               =  rmu(i) * sfinptl(ib)
            else
               c1(i)               =  0.0
            endif
         enddo

         DO610: do k = 1, lay
            do i = il1, il2
               taua(i,k)             =  absa(i,k,ib) * dp(i,k) + &
                    tauae(i,k,1) * absab(ib,1) + &
                    tauae(i,k,2) * absab(ib,2) + &
                    tauae(i,k,3) * absab(ib,3)
               tauci(i,k)            =  0.0
               omci(i,k)             =  0.0
               gci(i,k)              =  0.0
               f2(i,k)               =  0.0

               if (cldfrac(i,k) .ge. cut) then
                  tauci(i,k)          =  taucl(i,k,ib)
                  omci(i,k)           =  omcl(i,k,ib) * tauci(i,k)
                  f2(i,k)             =  gcl(i,k,ib) * gcl(i,k,ib)
                  gci(i,k)            = (gcl(i,k,ib) - f2(i,k)) / &
                       (1.0 - f2(i,k))
                  gci(i,k)            =  - 0.5 * (1.0 - uu3 * gci(i,k))
               endif
            enddo
         enddo DO610

         !----------------------------------------------------------------------
         !    reusing space o3g for dbf
         !----------------------------------------------------------------------

         call ccc2_planck2(bf, bs, urbf, a1(1,2), a1(1,3), o3g, tfull, gt, ib, &
              il1, il2, ilg, lay, lev)

         gh = .false.

         DOMAJLW: do ig = 1, kgl(ib)

            call ccc2_gasoptl7(taug, gw, dp, ib, ig, o3, qg, co2, ch4, an2o, &
                 f11, f12, f113, f114, inpr, inptm, mcont, pg, &
                 dip, dt, lev1, gh, il1, il2, ilg, lay)

            pgw = pi * gw
            !new           if(mcica.eq.0) then
            call ccc2_lwtran4(refl, tran, c1, tauci, omci, &
                 gci,  f2, taua, taug, bf, &
                 bs, urbf, o3g, em0,cldfrac, &
                 cldm, anu, nct, ncd, ncu, &
                 ncum, ncdm,lev1, cut, maxc, &
                 il1, il2, ilg, lay, lev)
            !new           else
            !new            call lwtran_mcica2(refl, tran, c1, tauci, omci, gci,
            !new     1                         f2, taua, taug, bf, bs, urbf, o3g, em0,
            !new     2                         cldg, nct, lev1, cut, maxc,
            !new     3                         il1, il2, ilg, lay, lev)

            !           * for the total sky fluxes weight the clear and cloudy sky
            !           * fluxes by the total vertically projected cloud fraction.

            !new            do k = lev1,lev
            !new            do i = il1,il2
            !new              if (cldt(i).lt.1.0) then
            !new                refl(i,2,k) = (1.0 - cldt(i)) * refl(i,1,k) +
            !new     1                         cldt(i)  * refl(i,2,k)
            !new                tran(i,2,k) = (1.0 - cldt(i)) * tran(i,1,k) +
            !new     1                         cldt(i)  * tran(i,2,k)
            !new              end if
            !new            end do ! i
            !new            end do ! k
            !new           endif
            !new           pgw = pi * gw * www

            do k = lev1, lay
               kp1 = k + 1
               do i = il1, il2
                  flxu(i,k)         =  flxu(i,k) + refl(i,2,k) * pgw
                  flxd(i,k)         =  flxd(i,k) + tran(i,2,k) * pgw

                  dfnet             =  tran(i,2,k) - tran(i,2,kp1) - &
                       refl(i,2,k) + refl(i,2,kp1)
                  hrl(i,k)          =  hrl(i,k) + &
                       hrcoef * dfnet / dp(i,k) * pgw
               enddo
            enddo

            do i = il1, il2
               flxu(i,lev)         =  flxu(i,lev) + refl(i,2,lev) * pgw
               flxd(i,lev)         =  flxd(i,lev) + tran(i,2,lev) * pgw

               clt(i)              =  clt(i) - refl(i,1,lev1) * pgw
               clb(i)              =  clb(i) - &
                    (refl(i,1,lev) - tran(i,1,lev)) * pgw
            enddo

            if (lev1 .gt. 1) then
               do k = lev1 - 1, 1, - 1
                  kp1 =  k + 1
                  do i = il1, il2
                     flxu(i,k)       =  flxu(i,k) + refl(i,2,lev1) * pgw
                     flxd(i,k)       =  flxd(i,k) + c1(i) * pgw
                  enddo
               enddo
            endif

         enddo DOMAJLW

         if (ib .ne. 6) then

            gh = .true.

            DOMINLW: do ig = 1, kglgh(ib)

               call ccc2_gasoptlgh7(taug, gwgh, dp, ib, ig, o3, qg, co2, ch4, &
                    an2o, inpt, mcont, dip, dt, lev1, gh, &
                    il1, il2, ilg, lay)


               !----------------------------------------------------------------------
               !     consider the attenuation for the downward flux above the model
               !     top level. this is important to get the correct cooling rate. if
               !     the model top level pressure is lower than 0.01. this is not
               !     necessary
               !----------------------------------------------------------------------

               call ccc2_lattenu4(a1, ib, ig, o3top, qg, co2, pfull, a1(1,12), &
                    dt, a1(1,5), inpt, il1, il2, ilg)

               do i = il1, il2
                  tran0 = exp(- a1(i,1))
                  if (pfull(i,1) .gt. 0.001) then
                     x               =  max(a1(i,1), 1.e-10)
                     ubeta0          = 1.6487213 * a1(i,3) / x
                     epsd0           = ubeta0 + 1.0
                     if (abs(epsd0) .gt. 0.001) then
                        c2(i)         =  c1(i) * tran0 + &
                             (bf(i,1) - a1(i,2) * tran0) / epsd0
                     else
                        c2(i)         =  c1(i)*tran0+x*a1(i,2)*tran0
                     endif
                  else
                     c2(i)           =  c1(i) * tran0
                  endif
               enddo

               call ccc2_lwtragh4(refl, tran, c2, tauci, omci, taua, taug, &
                    bf, urbf, cldfrac, em0, bs, cut, &
                    il1, il2, ilg, lay, lev)

               pgw = pi * gwgh
               !new          if(mcica.ne.0) then

               !           * for the total sky fluxes weight the clear and cloudy sky
               !           * fluxes by the total vertically projected cloud fraction.

               !new            do k = 1, lev
               !new            do i = il1,il2
               !new              if (cldt(i).lt.1.0) then
               !new                  refl(i,2,k) = (1.0 - cldt(i))*refl(i,1,k)
               !new     1                           + cldt(i)  * refl(i,2,k)
               !new                  tran(i,2,k) = (1.0 - cldt(i))*tran(i,1,k)
               !new     1                           + cldt(i)  * tran(i,2,k)
               !new              end if
               !new            end do ! i
               !new            end do ! k
               !new           endif

               !new           pgw = pi * gwgh * www

               do k = 1, lay
                  kp1 = k + 1
                  do i = il1, il2
                     flxu(i,k)       =  flxu(i,k) + refl(i,2,k) * pgw
                     flxd(i,k)       =  flxd(i,k) + tran(i,2,k) * pgw
                     dfnet           =  tran(i,2,k) - tran(i,2,kp1) - &
                          refl(i,2,k) + refl(i,2,kp1)
                     hrl(i,k)        =  hrl(i,k) + &
                          hrcoef * dfnet / dp(i,k) * pgw
                  enddo
               enddo

               !----------------------------------------------------------------------
               !     the attenuation for the upward flux above the model top is not
               !     considered, since the impact on upward flux is very small if the
               !     model top is about 1 mb or higher
               !----------------------------------------------------------------------

               do i = il1, il2
                  flxu(i,lev)       =  flxu(i,lev) + refl(i,2,lev) * pgw
                  flxd(i,lev)       =  flxd(i,lev) + tran(i,2,lev) * pgw
                  clt(i)            =  clt(i) -  refl(i,1,1) * pgw
                  clb(i)            =  clb(i) - &
                       (refl(i,1,lev) - tran(i,1,lev)) * pgw
               enddo

            enddo DOMINLW

         endif
      enddo DONBL

      do i = il1, il2
         fdl(i)                  =  flxd(i,lev)
         ful(i)                  =  flxu(i,1)
      enddo

      !     on veut les flux en sortie
      do k = 1,lev
         do i = il1,il2
            flxdl(i,k)=flxd(i,k)
            flxul(i,k)=flxu(i,k)
         enddo
      enddo

      !     decommenter cette partie si on fait lclw = false sinon ca plante
            else
               do i = il1, il2
                 fdl(i)       = 0.0
                 ful(i)       = 0.0
                 clt(i)       = 0.0
                 clb(i)       = 0.0
                 flxdl(i,lev) = 0.0
                 flxul(i,lev) = 0.0
               enddo
               do k = 1, lay
               do i = il1, il2
                 flxdl(i,k)   = 0.0
                 flxul(i,k)   = 0.0
                 hrl(i,k)     = 0.0
               enddo
               enddo
   endif
   !     (lclw)

   return
end subroutine ccc2_raddriv3
