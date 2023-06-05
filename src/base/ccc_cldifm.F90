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
!**S/P CLDIFM - GASEOUS CALCULATION
!
      subroutine ccc_cldifm1 (cldmin, cldmax, anu, a1, ncd, &
                         ncu, nblk, nct, ncum, ncdm, &
                         cldfrac, pfull, mrk2, lev1, cut, maxc, &
                         il1, il2, ilg, lay, lev)
        use ens_perturb, only: ens_nc2d, ens_spp_get
!
      implicit none
!!!#include <arch_specific.hf>
!
      integer ilg, lay, lev, lev1, maxc, il1, il2, i,k, km1, l, lp1
      real cut, x, y, z
      real, dimension(ilg) :: anumult
      real cldmin(ilg,lay), cldmax(ilg,lay), anu(ilg,lay), a1(ilg,12), &
           cldfrac(ilg,lay), pfull(ilg,lev), c1(ilg), mrk2(ilg,ens_nc2d)
!
      integer ncd(ilg,lay), ncu(ilg,lay), nblk(ilg,lay), nct(ilg), &
              ncum(lay), ncdm(lay), levc(ilg,lay), intg1(ilg), &
              intg2(ilg)
!
!Authors
!
!        J. Li, M. Lazare, CCCMA, rt code for gcm4
!        (Ref: J. Li, H. W. Barker, 2005:
!        JAS Vol. 62, no. 2, pp. 286\226309)
!        P. Vaillancourt, D. Talbot, RPN/CMC;
!        adapted for CMC/RPN physics (May 2006)
!
!Revisions
!
! 001    P.Vaillancourt, A.Plante (Feb 12) : maxc=max(maxc,lev1) ; necessary if cloud may occur close to 1mb
!                                            this problem can also be avoided through pressure threshold in nocld.cdk 
!
!Object
!
!        This subroutine determines the info for cloud and level info
!        for gaseous calculation
!
!Arguments
!
!         - Input -
! cldfrac cloud fraction
! pfull   pressure at model levels
! mrk2    Markov chains for parameter perturbations      
! cut     cloud fraction limit below which no cloud is considered
! il1     1
! il2     horizontal dimension
! ilg     horizontal dimension
! lay     number of model levels
!
!         - Output -
! cldmax  maximum cloud fraction for each cloud block.
! cldmin  minimum cloud fraction for each cloud block
! anu     nu factor for cloud subgrid variability
! a1      cloud fractions
! ncd     # of adjacent layers inside cloud block counted
!         from cloud top
! ncu     # of adjacent layers inside cloud block counted
!         from cloud base
! nblk    number of cloud blocks counted from surface
! nct     the level of the highest cloud at a model grid
! ncum    maximum loop number cloud vertical correlation accounted
!         from lower level to higher level (max of ncu)
! ncdm    maximum loop number cloud vertical correlation accounted
!         from higher level to lower level (max of ncd)
! lev1    a level close to 1 mb, below it the swtran start to work
! maxc    the level for the highest cloud in latitude-longitude
!         chain
! levc    a level close to 1 mb(2d)
! intg1   cloud base counter
! intg2   number of cloud layers in column
!
!*
!
      do 10 i = il1, il2
        intg1(i)                =  0
        intg2(i)                =  0
        nct(i)                  =  lev
        c1(i)                   =  0.0
        levc(i,1)               =  1
   10 continue

! Retrieve stochastic parameter information on request
      anumult = ens_spp_get('hetero_mult', mrk2, 1.)        
!
!----------------------------------------------------------------------
!     determine the highest cloud location. nct is the upper level of
!     the highest cloud,
!     determine the nu (anu) factor for cloud sub-grid variability
!     based on cloud fraction.
!----------------------------------------------------------------------
!
      do 25 k = 1, lay
        km1 = k - 1
        do 20 i = il1, il2
          if (cldfrac(i,k) .le. 0.9)                                then
            anu(i,k)            =  1.0
          elseif (cldfrac(i,k) .gt. 0.9 .and. cldfrac(i,k) .lt. 1.0)then
            anu(i,k)            =  2.0
          else
            anu(i,k)            =  4.0
          endif
          anu(i,k) = anu(i,k) * anumult(i)
!
!----------------------------------------------------------------------
!     minimum anu, it is extremely important to ensure consistency
!     between the definitions of anu here and their subsequent use in
!     lwtran,
!     cldmax the maximum cloud fraction for each cloud block.
!     cldmin the minimum cloud fraction for each cloud block.
!----------------------------------------------------------------------
!
          if (cldfrac(i,k) .lt. cut)                                then
            anu(i,k)            =  1000.0
            cldmax(i,k)         =  0.0
          else
            if (k .eq. 1)                                           then
              cldmax(i,k)       =  cldfrac(i,k)
            else
              anu(i,k)          =  min (anu(i,km1), anu(i,k))
              cldmax(i,k)       =  max (cldmax(i,km1), cldfrac(i,k))
            endif
!
            intg2(i)            =  intg2(i) + 1
            if (intg2(i) .eq. 1) nct(i) = k
          endif
!
!----------------------------------------------------------------------
!     determine lev1 for solar radiation
!----------------------------------------------------------------------
!
          if (pfull(i,k) .ge. 0.99)                                 then
            c1(i)               =  c1(i) + 1.0
            if (c1(i) .eq. 1.0)  levc(i,1) =  k
          endif
!
   20   continue
   25 continue
!
      lev1                      =  lev
      maxc                      =  lev
!
      a1(:,:)                   =  0.

      do 40 i = il1, il2
        maxc                    =  min (nct(i), maxc)
        lev1                    =  min (lev1, levc(i,1))
!
!----------------------------------------------------------------------
!     determine the layer order for each cloud block through down and
!     up paths, ncd and ncu.
!     determine the total cloud fractions looking from top and surface
!     for one cloud block (a cloud occupy several layers, choose the
!     minimum value of nu for the block.
!     nct is the top level number for the highest cloud
!     determine the minimum anu
!----------------------------------------------------------------------
!
        levc(i,1)               =  0
        levc(i,2)               =  0
   40 continue
        maxc=max(maxc,lev1)
!
      do 65 k = 2, lev
        km1 = k - 1
        l = lev - k + 1
        lp1 = l + 1
        do 60 i = il1, il2
          if (cldfrac(i,km1) .lt. cut)                              then
            levc(i,1)           =  0
            ncd(i,km1)          =  0
          else
            levc(i,1)           =  levc(i,1) + 1
            ncd(i,km1)          =  levc(i,1)
          endif
!
          if (cldfrac(i,l) .ge. cut .and. l .lt. lay)               then
            anu(i,l)            =  min (anu(i,lp1), anu (i,l))
            cldmax(i,l)         =  max (cldmax(i,lp1), cldmax(i,l))
          endif
   60   continue
   65 continue
!
      do 75 l = lay, 1, -1
        lp1 = l + 1
        do 70 i = il1, il2
          if (cldfrac(i,l) .lt. cut)                                then
            levc(i,2)           =  0
            ncu(i,l)            =  0
            nblk(i,l)           =  0
            cldmin(i,l)         =  1.
          else
            levc(i,2)           =  levc(i,2) + 1
            ncu(i,l)            =  levc(i,2)
            if (ncu(i,l) .eq. 1)                                    then
              intg1(i)          =  intg1(i) + 1
              nblk(i,l)         =  intg1(i)
              if (nblk(i,l) .gt. 3)  nblk(i,l) =  3
              if (nblk(i,l) .eq. 1)  a1(i,1)   =  cldmax(i,l)
              if (nblk(i,l) .eq. 2)  a1(i,2)   =  cldmax(i,l)
              if (nblk(i,l) .eq. 3)  a1(i,3)   = &
                                     max (a1(i,3), cldmax(i,l))
            else
              nblk(i,l)         =  nblk(i,lp1)
            endif
!
            if (ncu(i,l) .eq. 1)                                    then
              cldmin(i,l)       =  cldfrac(i,l)
            else
              cldmin(i,l)       =  min (cldmin(i,lp1), cldfrac(i,l))
            endif
          endif
   70   continue
   75 continue
!
      do 80 i = il1, il2
        x                       =  a1(i,3) * (1.0 - a1(i,1)) * &
                                  (1.0 - a1(i,2))
        a1(i,4)                 =  a1(i,1) * a1(i,2)
        a1(i,1)                 =  a1(i,1) * (1.0 - a1(i,2))
        if (a1(i,3) .ge. x + a1(i,2))                               then
          y                     =  a1(i,2)
          z                     =  a1(i,3) - x - a1(i,2)
        else
          y                     =  a1(i,3) - x
          z                     =  0.
        endif
!
        if (a1(i,3) .ge. x + a1(i,1))                               then
          a1(i,6)               =  a1(i,1)
          a1(i,5)               =  a1(i,3) - x - a1(i,6)
        else
          a1(i,6)               =  a1(i,3) - x
          a1(i,5)               =  0.
        endif
        a1(i,3)                 =  x
        a1(i,5)                 =  0.5 * (a1(i,5) + y)
        a1(i,6)                 =  0.5 * (a1(i,6) + z)
   80 continue
!
!----------------------------------------------------------------------
!     determine the maximum portion in a cloud block
!     determine the maximum number for ncd and ncu, for iteration in
!     longwave
!----------------------------------------------------------------------
!
      do 105 k = 1, lay
        km1 = k - 1
        ncum(k)                 =  0
        ncdm(k)                 =  0
        do 100 i = il1, il2
          if (ncd(i,k) .gt. 1)                                      then
            cldmin(i,k)         =  min (cldmin(i,km1), cldmin(i,k))
          endif
!
          ncum(k)               =  max (ncu(i,k), ncum(k))
          ncdm(k)               =  max (ncd(i,k), ncdm(k))
  100   continue
  105 continue
!
      return
      end
