!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!-------------------------------------- LICENCE END --------------------------------------
!** S/P AEROOPPRO - OPTICAL PROPERTIES OF AEROSOLS
!
      subroutine ccc_aerooppro (tauae,exta,exoma,exomga,fa,absa, &
                            tt,shtj,sig,ps,lat,mg,ml,pbl, &
                            aerosolback,il1,il2,ilg,lay,lev)
!
      implicit none
#include <arch_specific.hf>
#include "nbsnbl.cdk"
#include "tdpack_const.hf"
#include "ccc_aeros.cdk"
!
      integer ilg,lay,lev,il1,il2
      real tt(ilg,lay),sig(ilg,lay),shtj(ilg,lev)
      real ps(ilg),lat(ilg)
      real mg(ilg),ml(ilg),pbl(ilg)
      real tauae(ilg,lay,5)
      logical aerosolback
!
!     gathered and other work arrays used generally by solar.
      real exta(ilg,lay,nbs), exoma(ilg,lay,nbs), exomga(ilg,lay,nbs), &
           fa(ilg,lay,nbs)
!
!     work arrays used generally by longwave.
      real absa(ilg,lay,nbl)
!
!Author
!
!        P. Vaillancourt, RPN (May 2006)
!
!Revisions
!
! 001    P. Vaillancourt (Jan 07) - correct bug in calculation of dp closest to surface
!
!Object
!
!        Centralizes the calculation of the optical properties of
!        aerosols. as of May 2004, only climatological background
!        aerosols are considered.
!
!Arguments
!
!              - Output -
! tauae        background aerosol optical depth for solar and longwave
! exta         extinction coefficient for solar
! exoma        extinction coefficient times single scattering albedo
!              for solar
! exomga       exoma times asymmetry factor for solar
! fa           square of asymmetry factor for solar
! absa         absorption coefficient for longwave
!              - Input -
! tt           temperature
! sig          sigma levels
! shtj         sigma levels
! ps           surface pressure
! lat          latitude(radians)
! mg           land sea mask (0.0 - 1.0)
! ml           lake mask
! pbl          planetary boundary layer height
! aerosolback  logical key to control presence of background
!              aerosols
! il1          1
! il2          horizontal dimension
! ilg          horizontal dimension
! lay          number of model levels
! lev          number of flux levels (lay+1)
!
!*
!
!*********************************************************
!     AUTOMATIC ARRAYS
!*********************************************************
      real, dimension(ilg,lay) :: dz
      real, dimension(ilg) :: sdz
      real, dimension(ilg) :: ipbl
!*********************************************************

      real aird,dp
      real aero,ct
      real rec_grav,rec_rgasd
      integer i, k, ib, l, j
!
!
!     for background aerosol based on rpn rad, broad band results
!     for solar, the effect for longwave is small and neglected
!
!---------------------------------------------------------------------
!     initialize the band-dependant optical property arrays
!---------------------------------------------------------------------
!
      do 10 ib = 1, nbs
      do 10 k  = 1, lay
      do 10 i  = il1, il2
        exta  (i,k,ib)          =  1.0e-20
        exoma (i,k,ib)          =  1.0e-20
        exomga(i,k,ib)          =  1.0e-20
        fa    (i,k,ib)          =  0.0
  10  continue
!
      do 20 ib = 1, nbl
      do 20 k  = 1, lay
      do 20 i  = il1, il2
        absa  (i,k,ib)          =  1.0e-20
  20  continue
!
      do 30 l = 1, 5
      do 30 k = 1, lay
      do 30 i = il1, il2
        tauae(i,k,l)            =  1.e-10
  30  continue
!
!
!     background aerosol
!
      if (aerosolback)                                              then
!
      rec_grav=1./grav
      rec_rgasd=1./rgasd
!
      do 110 i  = il1, il2
         sdz(i)   = 0.
         ipbl(i)  = 0.
  110 continue
!
      do 120 k  = 1, lay
      do 120 i  = il1, il2
         aird   = sig(i,k)*ps(i)/tt(i,k)*rec_rgasd
         dp     = shtj(i,k+1) - shtj(i,k)
         dp     = max(dp*ps(i),0.)
         dz(i,k)= dp/aird*rec_grav
  120 continue
!
      do 130 k=lay,1,-1
      do 130 i=il1,il2
           if ( int(ipbl(i)).eq.0 ) then
!                   pbl heiht recomputed as sum of layer thicknesses
             sdz(i)=sdz(i)+dz(i,k)
!                   level closest to pbl
             if (sdz(i) .gt. pbl(i)) ipbl(i) = float(k)
           endif
 130  continue
!
      ct = 2.0 / pi
      do 150 k = 1, lay
      do 150 i = il1, il2
             if ( k.ge.int(ipbl(i)) ) then
               if ( mg(i).ge.0.5.or.ml(i).ge.0.5) then
!
!                  over land
!
                      aero = 0.25 - 0.2*ct * abs(lat(i))
                      tauae(i,k,1)= max(aero * dz(i,k)/sdz(i), 1.e-10)
                else
!
!                  over ocean
!
                      aero = 0.13 - 0.1*ct * abs(lat(i))
                      tauae(i,k,2)= max(aero  *dz(i,k)/sdz(i), 1.e-10)
                endif
             endif
  150 continue
!
!     else
!
!---------------------------------------------------------------------
!     the following subroutines are commented out because there are
!     no sources defined for these aerosols at the moment (may 2004)
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!     calculate the sulfate aerosol optical properties, the aerosol is
!     hygroscopic growth dependent
!---------------------------------------------------------------------
!
!      call sulfaero (exta, exoma, exomga, fa, absa, rh, so4load,
!     1               il1, il2, ilg, lay,
!     2               f1, f2, urbf, taucsg, tauoma, tauomga, tauomc)
!
!---------------------------------------------------------------------
!     calculate the sea salt (ss) aerosol optical properties,
!     the aerosol is hygroscopic growth dependent.
!     2 mode  1 (fine) and 2 (coarse)
!---------------------------------------------------------------------
!
!      call ssaltaero (exta, exoma, exomga, fa, absa, rh, ssload1,
!     1                ssload2, il1, il2, ilg, lay,
!     2                urbf, taucsg, tauoma, tauomga)
!
!---------------------------------------------------------------------
!     calculate the dust (ds) aerosol optical properties,
!     2 mode: 1 (fine) and 2 (coarse)
!---------------------------------------------------------------------
!
      endif
!
!
!---------------------------------------------------------------------
!     scaling aerosol optical properties. taua is aerosol optical depth
!---------------------------------------------------------------------
!
      return
      end
