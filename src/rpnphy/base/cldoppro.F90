!-------------------------------------- LICENCE BEGIN --------------------------
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
!-------------------------------------- LICENCE END ----------------------------


!/@*
subroutine cldoppro4(taucs, omcs, gcs, taucl, omcl, gcl, &
                           topthw, topthi,ecc,tcc, &
                           nlow,nmid,nhigh, &
                           ctp, ctt, liqwcin, icewcin, &
                           liqwpin, icewpin, cldfrac, &
                           tt, sig, ps, mg, ml, mrk2, m, &
                           lmx, nk, nkp)
   use tdpack_const
   use phy_options
   use ens_perturb, only: ens_nc2d, ens_spp_get
   implicit none
!!!#include <arch_specific.hf>
#include "nbsnbl.cdk"

   integer, intent(in) :: lmx, m, nk, nkp
   real taucs(lmx,nk,nbs), omcs(lmx,nk,nbs), gcs(lmx,nk,nbs), &
        taucl(lmx,nk,nbl), omcl(lmx,nk,nbl), gcl(lmx,nk,nbl)
   real topthw(lmx), topthi(lmx)
   real ctp(lmx), ctt(lmx)
   real liqwcin(lmx,nk), icewcin(lmx,nk)
   real liqwpin(lmx,nk), icewpin(lmx,nk)
   real cldfrac(lmx,nk), tt(m,nk),sig(lmx,nk),ps(lmx),mg(lmx),ml(lmx)
   real, dimension(lmx,ens_nc2d), intent(in) :: mrk2

!Authors
!        p. vaillancourt, d. talbot, j. li, rpn, cmc, cccma; (may 2006)
!
!Revisions
! 001    p.vaillancourt - code 3 options (ioptrew) for liquid effective radius
!                         and  2 options (ioptrei) for ice effective radius
! 002    p.vaillancourt - (april 2008) correct bug for diagnostic variables ctp (BP)and ctt (BE)
! 003    p.vaillancourt - (june 2008) add calculation of effective cloud cover (ecc=NT)
!
!Object
!        calculate cloud optical properties for ccc radiative transfer scheme
!        (adapted from cldoptx4)
!
!Arguments
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
!          - input -
! liqwpin  liquid water path in g/m2
! icewpin  solid water path in g/m2
! liqwcin  in-cloud liquid water content
!          in kg water/kg air
! icewcin  in-cloud ice water content
!          in kg water/kg air
! cldfrac  layer cloud amount (0. to 1.) (lmx,nk)
! tt       layer temperature (k) (m,nk)
! sig      sigma levels (0. to 1.) (lmx,nk; local sigma)
! ps       surface pressure (n/m2) (lmx)
! mg       ground cover (ocean=0.0,land <= 1.)  (lmx)
! ml       fraction of lakes (0.-1.) (lmx)
! mrk2     Markov chains for stochastic perturbations   
! lmx      number of profiles to compute
! m        first dimension of temperature (usually lmx)
! nk       number of layers
!
!***********************************************************************
!*@/

      external cldoppro_data

#include "cldop.cdk"

!*********************************************************
!     AUTOMATIC ARRAYS
!*********************************************************
      real, dimension(lmx,nk) :: transmissint
      real, dimension(lmx,nk) :: trans_exp
      logical, dimension(lmx) :: top
      real, dimension(m,nk) :: zrieff
      real, dimension(m,nk) :: aird
      real, dimension(m,nk) :: rew
      real, dimension(m,nk) :: rei
      real, dimension(m,nk) :: rec_cdd
      real, dimension(m,nk) :: vs1
      real, dimension(lmx) :: trmin
      real, dimension(lmx) :: tmem
      real, dimension(lmx) :: ecc
      real, dimension(lmx) :: tcc
      real, dimension(lmx) :: nlow
      real, dimension(lmx) :: nmid
      real, dimension(lmx) :: nhigh
      real, dimension(lmx) :: reifac
      real, dimension(lmx) :: rewfac
      real, dimension(lmx,nkp,nkp) :: ff
      integer, dimension(lmx     ) :: ih
      integer, dimension(lmx     ) :: ib


!*********************************************************

      real, parameter :: THIRD = 0.3333333 !#TODO test with 1./3. (bit pattern change)
      
      integer i, j, k, kind,ip,l
      real rec_grav, cut
      real epsilon,epsilon2,betan,betad
      real rew2,rew3,dg,dg2,dg3,tausw,omsw,gsw,tausi,omsi,gsi,y1,y2,taulw
      real omlw,glw,tauli,omli,gli
      real xnu

      rec_grav=1./grav
!
!
!     initialize output fields
!
      do i =1, lmx
        topthw(i) =  0.0
        topthi(i) =  0.0
      end do
!
      do k = 1, nk
          do i =1, lmx
!
!           effective radius for water clouds, set number of drops per
!           cm^3, 100 for water and 500 for land
!
            if (mg(i) .le. 0.5 .and. ml(i) .le. 0.5)                  then
!              cdd=50.
              rec_cdd(i,k) =  0.01
            else
!              cdd=250.
              rec_cdd(i,k) =  0.002
            endif
!
!           aird is air density in kg/m3
!
            aird(i,k) =  sig(i,k) * ps(i) / ( tt(i,k) * rgasd )
!
          end do
      end do

      ! Effective radius of droplets in liquid clouds
      select case (rad_cond_rew)
      case ('BARKER')
         ! Radius as in newrad: from H. Barker based on aircraft data (range 4-17um from Slingo)
         rew(:,:) = min(max(4., 754.6 * (liqwcin*aird*rec_cdd)**THIRD), 17.0)
      case ('NEWRAD')
         ! Radius as in newrad: corresponds to so called new optical properties
         vs1(:,:) = (1.0 + liqwcin(:,:) * 1.e4) &
              * liqwcin(:,:) * aird(:,:) * rec_cdd(:,:)
         rew(:,:) =  min(max(2.5, 3000. * vs1**THIRD), 50.0)
      case ('ROTSTAYN03')
         ! Radius according to Rotstayn and Liu (2003)
         do k = 1, nk
            do i = 1, lmx
               epsilon =  1.0 - 0.7 * exp(- 0.001 / rec_cdd(i,k))
               epsilon2 =  epsilon * epsilon
               betad =  1.0 + epsilon2
               betan =  betad + epsilon2
               rew(i,k) = 620.3504944*((betan*betan*liqwcin(i,k)*aird(i,k)) &
                    / (betad / rec_cdd(i,k)) )**third
               rew(i,k) =  min (max (2.5, rew(i,k)), 17.0)
            end do
         end do
      case DEFAULT
         ! Radius is a user-specified constant (in microns)
         rew = rew_const
      end select

      ! Adjust the effective radius using stochastic perturbations
      rewfac(:) = ens_spp_get('rew_mult', mrk2, default=1.)
      do k=1, nk
         rew(:,k) = rewfac(:) * rew(:,k)
      enddo


      ! Effective radius of crystals in ice clouds
      select case (rad_cond_rei)
      case ('CCCMA') 
         ! Units of icewcin must be in g/m3 for this parameterization of rei (in microns)
         zrieff(:,:) = (1000. * icewcin * aird)**0.216
         where (icewcin(:,:) >= 1.e-9)
            zrieff(:,:) = 83.8 * zrieff(:,:)
         elsewhere
            zrieff(:,:) = 20.
         endwhere
         rei(:,:) =  max(min(zrieff(:,:), 50.0), 20.0)
      case ('SIGMA')
         ! Radius varies from 60um (near-surface) to 15um (upper-troposphere)
         rei(:,:) = max(sig(:,:)-0.25, 0.0)*60. + 15.
      case DEFAULT
         ! Radius is a user-specified constant (in microns)
         rei(:,:) = rei_const
      end select

      ! Adjust the effective radius using stochastic perturbations
      reifac(:) = ens_spp_get('rei_mult', mrk2, default=1.)
      do k=1, nk
         rei(:,k) = reifac(:) * rei(:,k)
      enddo

!----------------------------------------------------------------------
!     cloud radiative properties for radiation.
!     taucs, omcs, gcs (taucl, omcl, gcl): optical depth, single
!     scattering albedo, asymmetry factor for solar (infrared).
!     rew: effective radiu (in micrometer) for water cloud
!     rei: effective radiu (in micrometer) for ice cloud
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
!----------------------------------------------------------------------
!
      cut = 0.001
      do 290 j = 1, nbs
        do 280 k = 1, nk
          do 270 i = 1, lmx
            if (cldfrac(i,k) .le. cut)                              then
              taucs(i,k,j) =  0.
              omcs(i,k,j)  =  0.
              gcs(i,k,j)   =  0.
            else
              rew2 =  rew(i,k) * rew(i,k)
              rew3 =  rew2 * rew(i,k)
              dg   =  1.5396 * rei(i,k)
              dg2  =  dg  * dg
              dg3  =  dg2 * dg
!
              if (liqwpin(i,k) .gt. 0.001)                          then
                tausw =  liqwpin(i,k) * &
                        (aws(1,j) + aws(2,j) / rew(i,k) + &
                         aws(3,j) / rew2 + aws(4,j) / rew3)
                omsw  =  1.0 - (bws(1,j) + bws(2,j) * rew(i,k) + &
                         bws(3,j) * rew2 + bws(4,j) * rew3)
                gsw   =  cws(1,j) + cws(2,j) * rew(i,k) + &
                         cws(3,j) * rew2 + cws(4,j) * rew3
              else
                tausw =  0.
                omsw  =  0.
                gsw   =  0.
              endif
!
              if (icewpin(i,k) .gt. 0.001)                          then
                tausi =  icewpin(i,k) * ( ais(1,j) + ais(2,j) / dg )
                omsi  =  1.0 - (bis(1,j) + bis(2,j) * dg + &
                         bis(3,j) * dg2 + bis(4,j) * dg3)
                gsi   =  cis(1,j) + cis(2,j) * dg + cis(3,j) * dg2 + &
                         cis(4,j) * dg3
              else
                tausi =  0.
                omsi  =  0.
                gsi   =  0.
              endif
!
              taucs(i,k,j)  =  tausw + tausi
              if (taucs(i,k,j) .gt. 0.0)                            then
                y1          =  omsw * tausw
                y2          =  omsi * tausi
                omcs(i,k,j) = (y1 + y2) / taucs(i,k,j)
                gcs (i,k,j) = (y1 * gsw + y2 * gsi) / (y1 + y2)
              else
                omcs(i,k,j) =  0.
                gcs (i,k,j) =  0.
              endif
!
!     calculate the optical depth for water and ice cloud in visible
!
              if (j .eq. 1)                                         then
                topthw(i) =  topthw(i) + tausw
                topthi(i) =  topthi(i) + tausi
              endif
            endif
270       continue
280     continue
290   continue
!
      do 390 j = 1, nbl
        do 380 k = 1, nk
          do 370 i = 1, lmx
            if (cldfrac(i,k) .le. cut)                              then
              taucl(i,k,j) =  0.
              omcl(i,k,j)  =  0.
              gcl(i,k,j)   =  0.
            else
              rew2 =  rew(i,k) * rew(i,k)
              rew3 =  rew2 * rew(i,k)
              dg   =  1.5396 * rei(i,k)
              dg2  =  dg  * dg
              dg3  =  dg2 * dg
!
              if (liqwpin(i,k) .gt. 0.001)                          then
                taulw =  liqwpin(i,k) * (awl(1,j) + awl(2,j) * rew(i,k)+ &
                         awl(3,j) / rew(i,k) + awl(4,j) / rew2 + &
                         awl(5,j) / rew3)
                omlw  =  1.0 - (bwl(1,j) + bwl(2,j) / rew(i,k) + &
                         bwl(3,j) * rew(i,k) + bwl(4,j) * rew2)
                glw   =  cwl(1,j) + cwl(2,j) / rew(i,k) + &
                         cwl(3,j) * rew(i,k) + cwl(4,j) * rew2
              else
                taulw =  0.
                omlw  =  0.
                glw   =  0.
              endif
!
!----------------------------------------------------------------------
!     since in fu etc. the param. is for absorptance, so need a factor
!     icewpin(i,k) / tauli for single scattering albedo
!----------------------------------------------------------------------
!
              if (icewpin(i,k) .gt. 0.001)                          then
                tauli =  icewpin(i,k) * (ail(1,j) + ail(2,j) / dg + &
                         ail(3,j) / dg2)
                omli  =  1.0 - (bil(1,j) / dg + bil(2,j) + &
                         bil(3,j) * dg + bil(4,j) * dg2) * &
                         icewpin(i,k) / tauli
                gli   =  cil(1,j) + cil(2,j) * dg + cil(3,j) * dg2 + &
                         cil(4,j) * dg3
              else
                tauli =  0.
                omli  =  0.
                gli   =  0.
              endif
!
              taucl(i,k,j)   =  taulw + tauli
              if (taucl(i,k,j) .gt. 0.0)                            then
                y1           =  omlw * taulw
                y2           =  omli * tauli
                omcl(i,k,j)  = (y1 + y2) / taucl(i,k,j)
                gcl (i,k,j)  = (glw * y1 + gli * y2) / (y1 + y2)
              else
                omcl(i,k,j)  =  0.
                gcl (i,k,j)  =  0.
              endif
            endif
370       continue
380     continue
390   continue
!
!     diagnostics: cloud top pressure (ctp) and temperature (ctt)
!     using the cloud optical depth at window region (band 6) to
!     calculate the emissivity
!
!     calcul des indices IH et IB pour nuages 2-D
!     IH = niveau le plus pres de sigma=0.4
!     IB = niveau le plus pres de sigma=0.7
!
        do k = 1, nk
          do i = 1, lmx
            trans_exp(i,k) = exp(- 1.64872 * taucl(i,k,6))
            if (sig(i,k).le.0.4) ih(i) = k
            if (sig(i,k).le.0.7) ib(i) = k

          enddo
        enddo

        do i = 1, lmx
          ctp (i)   = 110000.
          ctt (i)   = 310.
          top(i) = .true.
          transmissint(i,1) = 1. - cldfrac(i,1) * (1.-trans_exp(i,1) )
          if ( (1. - transmissint(i,1)) .gt. 0.99 .and. top(i) )    then
            ctp(i) = sig(i,1)*ps(i)
            ctt(i) = tt(i,1)
            top(i) = .false.
          end if
        end do
!
        do k = 2, nk
          do i = 1, lmx
!           transmissint(i,k)=transmissint(i,k-1) * (1. - cldfrac(i,k) *
!    1                          exp (- 1.64872 * taucl(i,k,6)))
            transmissint(i,k)=transmissint(i,k-1) * (1. - cldfrac(i,k) * &
                                (1.-trans_exp(i,k) ) )
            if ( (1. - transmissint(i,k)) .gt. 0.99 .and. top(i) )  then
              ctp(i) = sig(i,k)*ps(i)
              ctt(i) = tt(i,k)
              top(i) = .false.
            end if
          end do
        end do
!
!...  compute total, high cloud, middle cloud and low cloud effective cloud cover (nt) as in radir7
!     using the cloud optical depth at window region (band 6) to
!     calculate the emissivity

      do l=1,nkp-1
       do i=1,lmx
            ff(i,l,l)=1.
            tmem(i)=1.
            trmin(i)=1.
       enddo
       ip=l+1
       do k=ip,nkp
            kind=k-2
            kind=max0(kind,1)
            do i=1,lmx
               xnu=1.-cldfrac(i,k-1)*(1.-trans_exp(i,k-1) )
               if(cldfrac(i,kind).lt.0.01) then
                 tmem(i)= ff(i,l,k-1)
                 trmin(i)= xnu
               else
                 trmin(i)=min(trmin(i),xnu)
               endif
               ff(i,l,k)= tmem(i) * trmin(i)
            enddo
       enddo
      enddo
!
      do  i=1,lmx
         ecc(i)=1.-ff(i,1,nkp)
         nhigh(i) = 1. - ff(i, 1   ,IH(i))
         nmid (i) = 1. - ff(i,IH(i),IB(i))
         nlow (i) = 1. - ff(i,IB(i),nkp   )
      enddo

!...  compute total true cloud cover using maximum-random cloud overlap assumption

      do i=1,lmx
            ff(i,1,1)=1.
            tmem(i)=1.
            trmin(i)=1.
      enddo
      do k=2,nkp
            kind=k-2
            kind=max0(kind,1)
            do i=1,lmx
               xnu=1.-cldfrac(i,k-1)
               if(cldfrac(i,kind).lt.0.01) then
                 tmem(i)= ff(i,1,k-1)
                 trmin(i)= xnu
               else
                 trmin(i)=min(trmin(i),xnu)
               endif
               ff(i,1,k)= tmem(i) * trmin(i)
            enddo
      enddo
!
      do  i=1,lmx
         tcc(i)=1.-ff(i,1,nkp)
      enddo
!

      return
      end
!
      block data cldoppro_data
#include "nbsnbl.cdk"
#include "cldop.cdk"
!
!     --------------------------------------------------------
!     new water properties for sw (dobbie, li, and chylek 1999, jgr)
!     --------------------------------------------------------
!
        data ((aws(i,j), i = 1, 4), j = 1, nbs)           / &
        4.554e-04,   1.500e+00,   7.190e-01,  -9.419e-01, &
        3.859e-04,   1.508e+00,   9.512e-01,  -1.053e+00, &
       -3.946e-05,   1.538e+00,   1.035e+00,   2.638e-01, &
        2.936e-04,   1.541e+00,   1.698e-00,   1.521e+00  /
        data ((bws(i,j), i = 1, 4), j = 1, nbs)           / &
        6.481e-08,   1.553e-07,  -7.755e-10,   7.616e-12, &
        1.072e-05,   1.345e-05,  -1.799e-08,  -3.146e-11, &
        4.078e-04,   2.169e-03,  -2.177e-05,   1.506e-07, &
        2.013e-01,   1.109e-02,  -2.897e-04,   3.055e-06  /
        data ((cws(i,j), i = 1, 4), j = 1, nbs)           / &
        8.069e-01,   6.188e-03,  -2.065e-04,   2.352e-06, &
        7.685e-01,   9.337e-03,  -3.101e-04,   3.527e-06, &
        7.471e-01,   9.440e-03,  -2.616e-04,   2.614e-06, &
        7.956e-01,   8.138e-03,  -1.861e-04,   1.611e-06  /
!
!     --------------------------------------------------------
!     new water properties for lw (lindner and li, 2000 jcl)
!     --------------------------------------------------------
!
        data ((awl(i,j), i = 1, 5), j = 1, nbl)                      / &
         -.21671e-01, .79578e-03, .14899e+01, .62606e+01,-.12705e+02, &
         -.14126e+00, .28208e-02, .35125e+01,-.34541e+01,-.22679e+01, &
         -.18829e+00, .34065e-02, .46731e+01,-.11664e+02, .87105e+01, &
         -.16383e+00, .26574e-02, .48670e+01,-.16442e+02, .16128e+02, &
         -.20294e-01, .85110e-04, .28650e+01,-.11202e+02, .12047e+02, &
          .28752e-01,-.37315e-03, .14591e+01,-.48788e+01, .49725e+01, &
         -.40386e-01, .80822e-03, .25318e+01,-.64641e+01, .55609e+01, &
         -.48716e-01, .81275e-03, .30390e+01,-.97845e+01, .95101e+01, &
          .64794e-01,-.98530e-03, .12797e+01,-.55272e+01, .62599e+01 /
        data ((bwl(i,j), i = 1, 4), j = 1, nbl)          / &
          .36899e-02,-.54184e-03, .14561e-01,-.18451e-03, &
          .62141e-02, .61190e-01, .21127e-01,-.29731e-03, &
          .87326e-01, .29908e+00, .22928e-01,-.35569e-03, &
         -.37551e-01, .70237e+00, .26945e-01,-.37999e-03, &
          .51671e-01, .10199e+01, .18296e-01,-.21209e-03, &
          .52184e+00, .72352e+00,-.48090e-02, .10414e-03, &
          .57688e+00, .63008e+00,-.56325e-02, .87852e-04, &
          .50346e+00, .79407e+00,-.13179e-02, .25467e-04, &
          .67792e+00, .68259e+00,-.12136e-01, .20941e-03 /
        data ((cwl(i,j), i = 1, 4), j = 1, nbl)          / &
          .73147e+00, .11761e+00, .86402e-02,-.10761e-03, &
          .81284e+00,-.60287e-01, .45367e-02,-.33372e-04, &
          .92468e+00,-.39653e+00, .30494e-03, .20980e-04, &
          .10006e+01,-.71422e+00,-.46784e-02, .10114e-03, &
          .10635e+01,-.10097e+01,-.58726e-02, .97485e-04, &
          .10762e+01,-.12482e+01,-.40343e-02, .54330e-04, &
          .97445e+00,-.13875e+01, .79204e-03,-.27995e-04, &
          .79053e+00,-.13566e+01, .10452e-01,-.18111e-03, &
          .35512e+00,-.80671e+00, .30384e-01,-.47204e-03 /
!
!----------------------------------------------------------------------
!    ice fu 1997 jcl, fu et al. 1998 jcl
!----------------------------------------------------------------------
!
        data ((ais(i,j), i = 1, 2), j = 1, nbs) / &
       -0.24276e-04, 2.51884e+00, &
       -0.48500e-04, 2.52275e+00, &
       -0.98503e-05, 2.52048e+00, &
        0.24435e-03, 2.49116e+00                /
!
        data ((bis(i,j), i = 1, 4), j = 1, nbs)            / &
        0.13031e-06, 0.94102e-07,-0.75971e-10, 0.33977e-12, &
       -0.77603e-06, 0.73420e-05, 0.11514e-09,-0.90818e-12, &
        0.10007e-02, 0.10992e-02,-0.45043e-05, 0.12637e-07, &
        0.21201e+00, 0.25713e-02,-0.19228e-04, 0.62183e-07 /
!
        data ((cis(i,j), i = 1, 4), j = 1, nbs)            / &
        0.74821e+00, 0.92318e-03,-0.72862e-06,-0.95642e-08, &
        0.75227e+00, 0.10653e-02,-0.24930e-05,-0.29114e-08, &
        0.75553e+00, 0.17297e-02,-0.87585e-05, 0.19201e-07, &
        0.84323e+00, 0.20925e-02,-0.18302e-04, 0.60381e-07 /
!
        data ((ail(i,j), i = 1, 3), j = 1, nbl)    / &
       -.8839455e-03,  .2662598e+01,  .2196338e+01, &
       -.2066995e-02,  .2787904e+01,  .1397838e+01, &
       -.3085730e-02,  .2906257e+01, -.1911363e+01, &
       -.6968920e-02,  .3284275e+01, -.6973825e+01, &
       -.8372696e-02,  .3455018e+01, -.1516692e+02, &
       -.1691632e-02,  .2765756e+01, -.8331033e+01, &
       -.7098616e-02,  .3343404e+01, -.8144649e+01, &
       -.1041746e-01,  .3824226e+01, -.2255945e+02, &
        .5689700e-02,  .2285636e+01, -.1430752e+02 /
!
        data ((bil(i,j), i = 1, 4), j = 1, nbl)                   / &
        .5723611e+00,  .1627863e-01, -.1684272e-03,  .6061332e-06, &
        .4402328e+00,  .1736939e-01, -.1656608e-03,  .5709622e-06, &
        .8802908e+00,  .1249744e-01, -.1550609e-03,  .6105065e-06, &
        .6351165e+00,  .1781519e-01, -.1979682e-03,  .6000892e-06, &
        .5409536e-00,  .1949649e-01, -.2050908e-03,  .7364680e-06, &
        .1195515e+01,  .3350616e-02, -.5266996e-04,  .2233377e-06, &
        .1186334e+01,  .6213290e-02, -.1044277e-03,  .2233377e-06, &
        .2279562e+00,  .2017007e-01, -.1756872e-03,  .5703918e-06, &
        .7718967e+00,  .2120626e-01, -.2587649e-03,  .9878070e-06 /
!
        data ((cil(i,j), i = 1, 4), j = 1, nbl)                   / &
        .7975757e+00,  .3843973e-02, -.3540463e-04,  .1179791e-06, &
        .7947997e+00,  .3190423e-02, -.2386042e-04,  .6691811e-07, &
        .8737279e+00,  .2465886e-02, -.2468764e-04,  .8686448e-07, &
        .8577221e+00,  .2321034e-02, -.1897764e-04,  .8641223e-07, &
        .8906280e-00,  .1903269e-02, -.1733552e-04,  .5855071e-07, &
        .8663385e-00,  .2797934e-02, -.3187011e-04,  .1217209e-06, &
        .7644037e+00,  .4427001e-02, -.4494615e-04,  .1217209e-06, &
        .7200100e+00,  .3220301e-02, -.2195542e-04,  .6604318e-07, &
        .5355918e+00,  .1127081e-01, -.1234705e-03,  .4567953e-06 /
!
      end block data cldoppro_data
