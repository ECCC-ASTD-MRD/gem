!-------------------------------------- LICENCE BEGIN -------------------------
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

subroutine ccc_swtran(refl, tran, cumdtr, tran0, taua, &
     taur, taug, tauoma, tauomga, f1, &
     f2, taucs, tauomc, tauomgc, cldfrac, &
     cldm, a1, rmu, c1, c2, &
     albsur, nblk, nct, cut, lev1, &
     il1, il2, ilg, lay, lev)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer ilg, lay, lev, lev1, il1, il2, k, i, km1, l, lp1
   real cut
   real(REAL64) :: extopt, omars, ssalb, sf1, sf, tau1, om1, cow
   real(REAL64) :: ssgas1, cowg, alamd, u3, u2, uu, efun, efun2, rn
   real(REAL64) :: x1,x2,x3,x4,x5,x6,x7,x8,x9,yy
   real(REAL64) :: dm, gscw, appgm, apmgm, omarcs, sf2
   real(REAL64) :: tau2, om2, ssgas2, sdtr, srdf, stdf, srdr, stdr
   real xx,y1
   real atran0, dmm, fmm, fpp, umm, upp, tranpp, reflpp, dpp
   real refl(ilg,2,lev), tran(ilg,2,lev), cumdtr(ilg,4,lev), &
        tran0(ilg)
   real taua(ilg,lay), taur(ilg,lay), taug(ilg,lay), tauoma(ilg,lay), &
        tauomga(ilg,lay), f1(ilg,lay), f2(ilg,lay), taucs(ilg,lay), &
        tauomc(ilg,lay), tauomgc(ilg,lay), cldfrac(ilg,lay), &
        cldm(ilg,lay), a1(ilg,11),  rmu(ilg), c1(ilg), c2(ilg), &
        albsur(ilg)
   real(REAL64) :: rdf(ilg,4,lay), tdf(ilg,4,lay), rdr(ilg,4,lay), &
        tdr(ilg,4,lay), dtr(ilg,4,lay), rmdf(ilg,4,lev), &
        tmdr(ilg,4,lev), rmur(ilg,4,lev), rmuf(ilg,4,lev)
   integer nblk(ilg, lay), nct(ilg)

   !@Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !        JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)
   !@Revisions
   ! 001    P.Vaillancourt, M.Lazarre (Oct 2006) : singularity problem
   !        (loop 200) in clear and all sky solved with real(REAL64) promotion
   ! 002    P.Vaillancourt (Feb 2012) : singularity problem again; increase min of yy to 1.e-10
   !@Object DELTA-EDDINGTON APPROXIMATION
   !        Delta-eddington approximation and adding process for clear and
   !        all sky, the adding method by coakley et al (1983). This code
   !        can deal with solar radiative transfer through atmosphere with
   !        proper treatment of cloud overlap (random + maximum or random +
   !        slantwise) and cloud sub-grid variability. the theory for adding,
   !        Cloud overlap Li and Dobbie (2003). cloud sub-grid variability
   !        similar to with adjustment of cloud optical depth
   !@Arguments
   !          - Output -
   ! refl     reflectivity (1) clear sky; (2) all sky
   ! tran     transmitivity
   ! cumdtr   direct transmission for mult-layers
   ! rdf      layer diffuse reflection
   ! tdf      layer diffuse transmission
   ! rdr      layer direct reflection
   ! tdr      layer direct transmission
   ! dtr      direct transmission
   ! rmdf     block diffuse reflection from model top level
   ! tmdr     block direct transmission from model top level
   ! rmur     block direct reflection from model bottom level
   ! rmuf     block diffuse reflection from model bottom level
   !          - Input -
   ! tran0    attenuation factor for reducing solar flux from toa to
   !          model top level
   ! taua     aerosol optical depth
   ! taur     rayleigh optical depth
   ! taug     gaseous optical depth
   ! tauoma   aerosol optical depth times aerosol single scattering
   !          albedo
   ! tauomga  tauoma times aerosol asymmetry factor
   ! f1       square of aerosol asymmetry factor
   ! f2       square of cloud asymmetry factor
   ! taucs    cloud optical depth
   ! tauomc   cloud optical depth times cloud single scattering albedo
   ! tauomgc  tauomc times cloud asymmetry factor
   ! cldfrac  cloud fraction
   ! cldm     maximum portion in each cloud block, in which the exact
   !          solution for subgrid variability is applied
   ! a1       various relation for cloud overlap
   ! rmu      cos of solar zenith angle
   ! c1 , c2  two factors not dependent on ib and ig calculated
   !          outside for efficiency
   ! albsur   surface albedo
   ! nblk     number of cloud blocks accounted from surface
   ! nct      the highest cloud top level for the longitude and
   !          latitude loop (ilg)
   ! cut      cloud fraction limit below which no cloud is considered
   ! lev1     a level close to 1 mb, below it the swtran start to work
   ! il1      1
   ! il2      horizontal dimension
   ! ilg      horizontal dimension
   ! lay      number of model levels
   ! lev      number of flux levels (lay+1)
   
   !----------------------------------------------------------------------
   !     combine the optical properties for solar,
   !     1, aerosol + rayleigh + gas; 2, cloud + aerosol + rayleigh + gas
   !     calculate the direct and diffuse reflection and transmission in
   !     the scattering layers using the delta-eddington method.
   !----------------------------------------------------------------------
   
   DO200: do k = lev1, lay
      do i = il1, il2
         extopt                    =  taua(i,k) + taur(i,k) + taug(i,k) + &
              1.0e-20
         omars                     =  tauoma(i,k) + taur(i,k)
         ssalb                     =  omars / extopt
         sf1                       =  f1(i,k) / omars
         sf                        =  ssalb * sf1
         tau1                      =  extopt * (1.0 - sf)
         om1                       = (ssalb - sf) / (1.0 - sf)
         cow                       =  1.0 - om1 + 1.e-10
         ssgas1                    = (tauomga(i,k) / omars - sf1) / &
              (1.0 - sf1)
         cowg                      =  1.0 - om1 * ssgas1

         dtr(i,1,k)                =  exp( - tau1 / rmu(i))
         alamd                     =  sqrt(3.0 * cow * cowg)
         u3                        =  1.50 * cowg / alamd
         u2                        =  u3 + u3
         uu                        =  u3 * u3
         efun                      =  exp(- alamd * tau1)
         efun2                     =  efun * efun
         x1                        = (uu - u2 + 1.0) * efun2
         rn                        =  1.0 / (uu + u2 + 1.0 - x1)
         rdf(i,1,k)                = (uu - 1.0) * (1.0  - efun2) * rn
         tdf(i,1,k)                = (u2 + u2) * efun * rn
         x2                        =  alamd * rmu(i)
         yy                        =  1.0 - x2 * x2
         yy                        = dsign (max (dabs(yy),1.d-10),yy)
         dm                        =  om1 / yy
         gscw                      =  ssgas1 * cow
         appgm                     = (c1(i) + 0.50 + &
              gscw * (c1(i) + c2(i))) * dm
         apmgm                     = (c1(i) - 0.50 + &
              gscw * (c1(i) - c2(i))) * dm
         rdr(i,1,k)                =  appgm * rdf(i,1,k) + apmgm * &
              (tdf(i,1,k) * dtr(i,1,k) - 1.0)
         tdr(i,1,k)                =  appgm * tdf(i,1,k) + &
              (apmgm * rdf(i,1,k) - appgm + 1.0) * &
              dtr(i,1,k)

         if (cldfrac(i,k) .lt. cut)                                  then
            rdf(i,2,k)              =  rdf(i,1,k)
            tdf(i,2,k)              =  tdf(i,1,k)
            rdr(i,2,k)              =  rdr(i,1,k)
            tdr(i,2,k)              =  tdr(i,1,k)
            dtr(i,2,k)              =  dtr(i,1,k)
            rdf(i,3,k)              =  rdf(i,1,k)
            tdf(i,3,k)              =  tdf(i,1,k)
            rdr(i,3,k)              =  rdr(i,1,k)
            tdr(i,3,k)              =  tdr(i,1,k)
            dtr(i,3,k)              =  dtr(i,1,k)
            rdf(i,4,k)              =  rdf(i,1,k)
            tdf(i,4,k)              =  tdf(i,1,k)
            rdr(i,4,k)              =  rdr(i,1,k)
            tdr(i,4,k)              =  tdr(i,1,k)
            dtr(i,4,k)              =  dtr(i,1,k)
         else
            extopt                  =  taucs(i,k) + extopt
            omarcs                  =  tauomc(i,k) + taur(i,k)
            ssalb                   =  omarcs / extopt
            sf2                     =  f2(i,k) / omarcs
            sf                      =  ssalb * sf2
            tau2                    =  extopt * (1.0 - sf)
            om2                     = (ssalb - sf) / (1.0 - sf)
            cow                     =  1.0 - om2
            ssgas2                  = (tauomgc(i,k) / omarcs - sf2) / &
                 (1.0 - sf2)
            cowg                    =  1.0 - om2 * ssgas2
            alamd                   =  sqrt(3.0 * cow * cowg)
            u3                      =  1.50 * cowg / alamd
            u2                      =  u3 + u3
            uu                      =  u3 * u3
            sdtr                    =  exp(- tau2 / rmu(i))
            efun                    =  exp(- alamd * tau2)
            efun2                   =  efun * efun
            x3                      = (uu - u2 + 1.0) * efun2
            rn                      =  1.0 / (uu + u2 + 1.0 - x3)
            x4                      =  alamd * rmu(i)
            yy                      =  1.0 - x4 * x4
            yy                      = dsign (max (dabs(yy),1.d-10),yy)
            dm                      =  om2 / yy
            gscw                    =  ssgas2 * cow
            appgm                   = (c1(i) + 0.50 + &
                 gscw * (c1(i) + c2(i))) * dm
            apmgm                   = (c1(i) - 0.50 + &
                 gscw * (c1(i) - c2(i))) * dm
            srdf                    = (uu - 1.0) * (1.0 - efun2) * rn
            stdf                    = (u2 + u2) * efun * rn
            srdr                    =  appgm * srdf + apmgm * &
                 (stdf * sdtr - 1.0)
            stdr                    =  appgm * stdf + (apmgm * srdf - &
                 appgm + 1.0) * sdtr
            if (nblk(i,k) .eq. 3)                                     then
               x5                    =  a1(i,9) * cldfrac(i,k)
               rdf(i,2,k)            =  rdf(i,1,k) + x5 * &
                    (srdf - rdf(i,1,k))
               tdf(i,2,k)            =  tdf(i,1,k) + x5 * &
                    (stdf - tdf(i,1,k))
               rdr(i,2,k)            =  rdr(i,1,k) + x5 * &
                    (srdr - rdr(i,1,k))
               tdr(i,2,k)            =  tdr(i,1,k) + x5 * &
                    (stdr - tdr(i,1,k))
               dtr(i,2,k)            =  dtr(i,1,k) + x5 * &
                    (sdtr - dtr(i,1,k))
               x6                    =  a1(i,10) * cldfrac(i,k)
               rdf(i,3,k)            =  rdf(i,1,k) + x6 * &
                    (srdf - rdf(i,1,k))
               tdf(i,3,k)            =  tdf(i,1,k) + x6 * &
                    (stdf - tdf(i,1,k))
               rdr(i,3,k)            =  rdr(i,1,k) + x6 * &
                    (srdr - rdr(i,1,k))
               tdr(i,3,k)            =  tdr(i,1,k) + x6 * &
                    (stdr - tdr(i,1,k))
               dtr(i,3,k)            =  dtr(i,1,k) + x6 * &
                    (sdtr - dtr(i,1,k))
               x7                    =  a1(i,11) * cldfrac(i,k)
               rdf(i,4,k)            =  rdf(i,1,k) + x7 * &
                    (srdf - rdf(i,1,k))
               tdf(i,4,k)            =  tdf(i,1,k) + x7 * &
                    (stdf - tdf(i,1,k))
               rdr(i,4,k)            =  rdr(i,1,k) + x7 * &
                    (srdr - rdr(i,1,k))
               tdr(i,4,k)            =  tdr(i,1,k) + x7 * &
                    (stdr - tdr(i,1,k))
               dtr(i,4,k)            =  sdtr
            else if (nblk(i,k) .eq. 1)                                then
               rdf(i,4,k)            =  rdf(i,1,k)
               tdf(i,4,k)            =  tdf(i,1,k)
               rdr(i,4,k)            =  rdr(i,1,k)
               tdr(i,4,k)            =  tdr(i,1,k)
               dtr(i,4,k)            =  dtr(i,1,k)
               x8                    =  cldfrac(i,k) / cldm(i,k)
               rdf(i,2,k)            =  rdf(i,1,k) + x8 * &
                    (srdf - rdf(i,1,k))
               tdf(i,2,k)            =  tdf(i,1,k) + x8 * &
                    (stdf - tdf(i,1,k))
               rdr(i,2,k)            =  rdr(i,1,k) + x8 * &
                    (srdr - rdr(i,1,k))
               tdr(i,2,k)            =  tdr(i,1,k) + x8 * &
                    (stdr - tdr(i,1,k))
               dtr(i,2,k)            =  sdtr
               if (a1(i,2) .ge. cut)                                   then
                  yy                  =  x8 * a1(i,8)
                  rdf(i,3,k)          =  rdf(i,1,k) + yy * &
                       (srdf - rdf(i,1,k))
                  tdf(i,3,k)          =  tdf(i,1,k) + yy * &
                       (stdf - tdf(i,1,k))
                  rdr(i,3,k)          =  rdr(i,1,k) + yy * &
                       (srdr - rdr(i,1,k))
                  tdr(i,3,k)          =  tdr(i,1,k) + yy * &
                       (stdr - tdr(i,1,k))
                  dtr(i,3,k)          =  dtr(i,1,k) + yy * &
                       (sdtr - dtr(i,1,k))
               else
                  rdf(i,3,k)          =  rdf(i,1,k)
                  tdf(i,3,k)          =  tdf(i,1,k)
                  rdr(i,3,k)          =  rdr(i,1,k)
                  tdr(i,3,k)          =  tdr(i,1,k)
                  dtr(i,3,k)          =  dtr(i,1,k)
               endif
            else if (nblk(i,k) .eq. 2)                                then
               rdf(i,2,k)            =  rdf(i,1,k)
               tdf(i,2,k)            =  tdf(i,1,k)
               rdr(i,2,k)            =  rdr(i,1,k)
               tdr(i,2,k)            =  tdr(i,1,k)
               dtr(i,2,k)            =  dtr(i,1,k)
               rdf(i,4,k)            =  rdf(i,1,k)
               tdf(i,4,k)            =  tdf(i,1,k)
               rdr(i,4,k)            =  rdr(i,1,k)
               tdr(i,4,k)            =  tdr(i,1,k)
               dtr(i,4,k)            =  dtr(i,1,k)
               x9                    =  cldfrac(i,k) / cldm(i,k)
               rdf(i,3,k)            =  rdf(i,1,k) + x9 * &
                    (srdf - rdf(i,1,k))
               tdf(i,3,k)            =  tdf(i,1,k) + x9 * &
                    (stdf - tdf(i,1,k))
               rdr(i,3,k)            =  rdr(i,1,k) + x9 * &
                    (srdr - rdr(i,1,k))
               tdr(i,3,k)            =  tdr(i,1,k) + x9 * &
                    (stdr - tdr(i,1,k))
               dtr(i,3,k)            =  sdtr
            endif
         endif
      enddo
   enddo DO200
 
   DO300: do i = il1, il2

      !----------------------------------------------------------------------
      !     initialization for the first level (lev1).
      !----------------------------------------------------------------------

      atran0                    =  1.0 - tran0(i)
      tmdr(i,1,lev1)            =  tran0(i)
      rmdf(i,1,lev1)            =  atran0
      cumdtr(i,1,lev1)          =  tran0(i)
      tmdr(i,2,lev1)            =  tran0(i)
      rmdf(i,2,lev1)            =  atran0
      cumdtr(i,2,lev1)          =  tran0(i)
      tmdr(i,3,lev1)            =  tran0(i)
      rmdf(i,3,lev1)            =  atran0
      cumdtr(i,3,lev1)          =  tran0(i)
      tmdr(i,4,lev1)            =  tran0(i)
      rmdf(i,4,lev1)            =  atran0
      cumdtr(i,4,lev1)          =  tran0(i)

      !----------------------------------------------------------------------
      !     initialization for the ground layer.
      !----------------------------------------------------------------------

      rmur(i,1,lev)             =  albsur(i)
      rmuf(i,1,lev)             =  albsur(i)
      rmur(i,2,lev)             =  albsur(i)
      rmuf(i,2,lev)             =  albsur(i)
      rmur(i,3,lev)             =  albsur(i)
      rmuf(i,3,lev)             =  albsur(i)
      rmur(i,4,lev)             =  albsur(i)
      rmuf(i,4,lev)             =  albsur(i)
   enddo DO300
   
   !----------------------------------------------------------------------
   !     add the layers downward from the second layer to the surface.
   !----------------------------------------------------------------------
   
   DO450: do k = lev1 + 1, lev
      km1 = k - 1
      l = lev - k + lev1
      lp1 = l + 1
      DO400: do i = il1, il2
         dmm                     =  tdf(i,1,km1) / &
              (1.0 - rdf(i,1,km1) * rmdf(i,1,km1))
         fmm                     =  rmdf(i,1,km1) * dmm
         tmdr(i,1,k)             =  cumdtr(i,1,km1) * (tdr(i,1,km1) + &
              rdr(i,1,km1) * fmm) + &
              (tmdr(i,1,km1) - cumdtr(i,1,km1)) * &
              dmm
         rmdf(i,1,k)             =  rdf(i,1,km1) + tdf(i,1,km1) * fmm
         cumdtr(i,1,k)           =  cumdtr(i,1,km1) * dtr(i,1,km1)

         if (a1(i,1) .ge. cut)                                     then
            if (k .le. nct(i))                                      then
               tmdr(i,2,k)         =  tmdr(i,1,k)
               rmdf(i,2,k)         =  rmdf(i,1,k)
               cumdtr(i,2,k)       =  cumdtr(i,1,k)
            else
               dpp                 =  tdf(i,2,km1) / &
                    (1.0 - rmdf(i,2,km1) * rdf(i,2,km1))
               fpp                 =  rmdf(i,2,km1) * dpp
               tmdr(i,2,k)         =  cumdtr(i,2,km1) * (tdr(i,2,km1) + &
                    rdr(i,2,km1) * fpp) + &
                    (tmdr(i,2,km1) - cumdtr(i,2,km1)) * &
                    dpp
               rmdf(i,2,k)         =  rdf(i,2,km1) + tdf(i,2,km1) * fpp
               cumdtr(i,2,k)       =  cumdtr(i,2,km1) * dtr(i,2,km1)
            endif
         else
            tmdr(i,2,k)           =  1.0
            rmdf(i,2,k)           =  0.0
            cumdtr(i,2,k)         =  0.0
         endif

         if (a1(i,2) .ge. cut)                                     then
            if (k .le. nct(i))                                      then
               tmdr(i,3,k)         =  tmdr(i,1,k)
               rmdf(i,3,k)         =  rmdf(i,1,k)
               cumdtr(i,3,k)       =  cumdtr(i,1,k)
            else
               dpp                 =  tdf(i,3,km1) / &
                    (1.0 - rmdf(i,3,km1) * rdf(i,3,km1))
               fpp                 =  rmdf(i,3,km1) * dpp
               tmdr(i,3,k)         =  cumdtr(i,3,km1) * (tdr(i,3,km1) + &
                    rdr(i,3,km1) * fpp) + &
                    (tmdr(i,3,km1) - cumdtr(i,3,km1)) * &
                    dpp
               rmdf(i,3,k)         =  rdf(i,3,km1) + tdf(i,3,km1) * fpp
               cumdtr(i,3,k)       =  cumdtr(i,3,km1) * dtr(i,3,km1)
            endif

            if (a1(i,3) .ge. cut)                                   then
               if (k .le. nct(i))                                    then
                  tmdr(i,4,k)       =  tmdr(i,1,k)
                  rmdf(i,4,k)       =  rmdf(i,1,k)
                  cumdtr(i,4,k)     =  cumdtr(i,1,k)
               else
                  dpp               =  tdf(i,4,km1) / &
                       (1.0 - rmdf(i,4,km1) * rdf(i,4,km1))
                  fpp               =  rmdf(i,4,km1) * dpp
                  tmdr(i,4,k)       =  cumdtr(i,4,km1) * (tdr(i,4,km1) + &
                       rdr(i,4,km1) * fpp) + &
                       (tmdr(i,4,km1) - cumdtr(i,4,km1)) * &
                       dpp
                  rmdf(i,4,k)       =  rdf(i,4,km1) + tdf(i,4,km1) * fpp
                  cumdtr(i,4,k)     =  cumdtr(i,4,km1) * dtr(i,4,km1)
               endif
            else
               tmdr(i,4,k)         =  1.0
               rmdf(i,4,k)         =  0.0
               cumdtr(i,4,k)       =  0.0
            endif
         else
            tmdr(i,3,k)           =  1.0
            rmdf(i,3,k)           =  0.0
            cumdtr(i,3,k)         =  0.0
            tmdr(i,4,k)           =  1.0
            rmdf(i,4,k)           =  0.0
            cumdtr(i,4,k)         =  0.0
         endif

         !----------------------------------------------------------------------
         !     add the layers upward from one layer above surface to the lev1.
         !----------------------------------------------------------------------

         umm                     =  tdf(i,1,l) / &
              (1.0 - rdf(i,1,l) * rmuf(i,1,lp1))
         fmm                     =  rmuf(i,1,lp1) * umm
         rmur(i,1,l)             =  rdr(i,1,l) + dtr(i,1,l) * &
              rmur(i,1,lp1) * umm + (tdr(i,1,l) - &
              dtr(i,1,l)) * fmm
         rmuf(i,1,l)             =  rdf(i,1,l) + tdf(i,1,l) * fmm

         if (a1(i,1) .ge. cut)                                     then
            upp                   =  tdf(i,2,l) / &
                 (1.0 - rmuf(i,2,lp1) * rdf(i,2,l))
            fpp                   =  rmuf(i,2,lp1) * upp
            rmur(i,2,l)           =  rdr(i,2,l) + dtr(i,2,l) * &
                 rmur(i,2,lp1) * upp + (tdr(i,2,l) - &
                 dtr(i,2,l)) * fpp
            rmuf(i,2,l)           =  rdf(i,2,l) + tdf(i,2,l) * fpp
         else
            rmur(i,2,l)           =  0.0
            rmuf(i,2,l)           =  0.0
         endif

         if (a1(i,2) .ge. cut)                                     then
            upp                   =  tdf(i,3,l) / &
                 (1.0 - rmuf(i,3,lp1) * rdf(i,3,l))
            fpp                   =  rmuf(i,3,lp1) * upp
            rmur(i,3,l)           =  rdr(i,3,l) + dtr(i,3,l) * &
                 rmur(i,3,lp1) * upp + (tdr(i,3,l) - &
                 dtr(i,3,l)) * fpp
            rmuf(i,3,l)           =  rdf(i,3,l) + tdf(i,3,l) * fpp

            if (a1(i,3) .ge. cut)                                   then
               upp                 =  tdf(i,4,l) / &
                    (1.0 - rmuf(i,4,lp1) * rdf(i,4,l))
               fpp                 =  rmuf(i,4,lp1) * upp
               rmur(i,4,l)         =  rdr(i,4,l) + dtr(i,4,l) * &
                    rmur(i,4,lp1) * upp + (tdr(i,4,l) - &
                    dtr(i,4,l)) * fpp
               rmuf(i,4,l)         =  rdf(i,4,l) + tdf(i,4,l) * fpp
            else
               rmur(i,4,l)         =  0.0
               rmuf(i,4,l)         =  0.0
            endif
         else
            rmur(i,3,l)           =  0.0
            rmuf(i,3,l)           =  0.0
            rmur(i,4,l)           =  0.0
            rmuf(i,4,l)           =  0.0
         endif
      enddo DO400
   enddo DO450
 
   !----------------------------------------------------------------------
   !     add downward to calculate the resultant reflectance and
   !     transmittance at flux levels.
   !----------------------------------------------------------------------
 
   DO550: do k = lev1, lev
      do i = il1, il2
         dmm                     =  1.0 / &
              (1.0 - rmuf(i,1,k) * rmdf(i,1,k))
         xx                      =  cumdtr(i,1,k) * rmur(i,1,k)
         y1                      =  tmdr(i,1,k) - cumdtr(i,1,k)
         tran(i,1,k)             =  cumdtr(i,1,k) + &
              (xx * rmdf(i,1,k) + y1) * dmm
         refl(i,1,k)             = (xx + y1 * rmuf(i,1,k)) * dmm

         if (a1(i,1) .ge. cut)                                     then
            dpp                   =  1.0 / &
                 (1.0 - rmuf(i,2,k) * rmdf(i,2,k))
            xx                    =  cumdtr(i,2,k) * rmur(i,2,k)
            y1                    =  tmdr(i,2,k) - cumdtr(i,2,k)
            tran(i,2,k)           =  a1(i,1) * (cumdtr(i,2,k) + &
                 (xx * rmdf(i,2,k) + y1) * dpp) + &
                 a1(i,7) * tran(i,1,k)
            refl(i,2,k)           =  a1(i,1) * (xx + y1 * rmuf(i,2,k)) * &
                 dpp + a1(i,7) * refl(i,1,k)
         else
            tran(i,2,k)           =  a1(i,7) * tran(i,1,k)
            refl(i,2,k)           =  a1(i,7) * refl(i,1,k)
         endif

         if (a1(i,2) .ge. cut)                                     then
            dpp                   =  1.0 / &
                 (1.0 - rmuf(i,3,k) * rmdf(i,3,k))
            xx                    =  cumdtr(i,3,k) * rmur(i,3,k)
            y1                    =  tmdr(i,3,k) - cumdtr(i,3,k)
            tranpp                =  cumdtr(i,3,k) + &
                 (xx * rmdf(i,3,k) + y1) * dpp
            reflpp                = (xx + y1 * rmuf(i,3,k)) * dpp
            tran(i,2,k)           =  a1(i,2) * tranpp + tran(i,2,k)
            refl(i,2,k)           =  a1(i,2) * reflpp + refl(i,2,k)

            if (a1(i,3) .ge. cut)                                   then
               dpp                 =  1.0 / &
                    (1.0 - rmuf(i,4,k) * rmdf(i,4,k))
               xx                  =  cumdtr(i,4,k) * rmur(i,4,k)
               y1                  =  tmdr(i,4,k) - cumdtr(i,4,k)
               tranpp              =  cumdtr(i,4,k) + &
                    (xx * rmdf(i,4,k) + y1) * dpp
               reflpp              = (xx + y1 * rmuf(i,4,k)) * dpp
               tran(i,2,k)         =  a1(i,3) * tranpp + tran(i,2,k)
               refl(i,2,k)         =  a1(i,3) * reflpp + refl(i,2,k)
            endif
         endif
      enddo
   enddo DO550

   return
end subroutine ccc_swtran
