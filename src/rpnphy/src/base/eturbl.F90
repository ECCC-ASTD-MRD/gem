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

!/@*
subroutine ETURBL13(EN,ENOLD,ZN,ZD,RIF,TURBREG,RIG,SHR2,GAMA,HOL,FN, &
     GAMAL,U,V,T,TE,TVE,Q,QCE,QE,PS,S,SE, &
     TAU,KOUNT,GAMAQ,KT,Z,GZMOM,P_PROF,FRV,XH, &
     DXDY,TRNCH,N,NK,Z0)
   !#TODO: never  used: QL, H, LH, TS, CCS, IT
   use, intrinsic :: iso_fortran_env, only: INT64
   use tdpack, only: CAPPA, DELTA, KARMAN
   use series_mod, only: series_xst
   use phy_options
   use phy_status, only: phy_error_L, PHY_OK
   use mixing_length, only: ml_blend,ml_calc_blac,ml_calc_boujo
   use pbl_stabfunc, only: psf_stabfunc
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
   integer TRNCH,N,NK
   real EN(N,NK),ENOLD(N,NK),ZN(N,NK),ZD(N,NK),RIF(N,NK),TURBREG(N,NK),RIG(N,NK),SHR2(N,NK)
   real GAMA(N,NK),FN(N,NK),XH(N)
   real HOL(N),U(N,NK),V(N,NK)
   real GAMAL(N,NK),DXDY(N)
   real T(N,NK),TE(N,NK),TVE(N,NK),Q(N,NK),QCE(N,NK),QE(N,NK),PS(N)
   real S(n,NK),SE(n,NK),P_PROF(NK)
   real TAU
   integer KOUNT
   real LMN,Z0(N)
   real FRV(N)
   real KT(N,NK),GAMAQ(N,NK)
   real Z(N,NK),GZMOM(N,NK)
   real EXP_TAU

   !@Object predict EN(turbulent energy) and ZN(mixing length)

   !@Author J. Cote (RPN 1983)

   !@Revision
   ! 001      J. Cote RPN(Nov 1984)SEF version documentation
   ! 002      M. Lepine  -  RFE model code revision project (Feb 87)
   !                      -  Remove COMMON WKL2D1 and pass the
   !                         work field in parameter
   ! 003      J.Mailhot RPN(Sep 1985) Series (RIF,BILAN EN)
   ! 004      J.Mailhot RPN(Oct 1985) Scaling (ZN,ZE,C)
   ! 005      J.Mailhot RPN(Oct 1985) Add countergradient term
   !               Adaptation to revised code G.Pellerin (Oct87)
   ! 006      J. Mailhot-G.Pellerin (Nov 87) - Correction to the
   !                vertical diffusion of the turbulent energy
   ! 007      J.Mailhot-G.Pellerin (Apr88)
   !               Return to the old formula for stable LAMBDA
   !               (with relaxation term (noise))
   !               countergradient term to zero.
   ! 008      MJ L'Heureux  (Mar89) Initializations of fields
   !                                at KOUNT=0
   ! 009      R.Benoit (Mar89)   -Y. Delage (May89)
   !               Revision of vertical diffusion
   ! 010      Y. Delage (Jan90)
   !                Return DIAGSF and ETURBL coherents with FLXSRF for
   !                the calculation of unstable diffusion coefficients
   ! 011      N. Brunet  (May90)
   !                Standardization of thermodynamic functions
   ! 012      J.Mailhot RPN(Feb 1990) Shallow convection (GELEYN)
   ! 013      G.Pellerin(August90)Adaptation to thermo functions
   ! 014      Y. Delage  (Nov 1990) Options of shallow convection
   ! 015      Y. Delage  (Nov1990)
   !                  Removal OFA,EA,PRI and BETA
   !                  Replace WC and HOL by ILMO
   ! 016      N. Brunet  (May91)
   !                New version of thermodynamic functions
   !                and file of constants
   ! 017      B. Bilodeau  (July 1991)- Adaptation to UNIX
   !
   ! 018      C. Girard (Nov 1992)
   !          Modification of the shallow convection:
   !          - end of GELHU option
   !          - new significance for GELEYN option
   !          Modification to definitions:
   !          - neutral mixing length
   !          - stability functions
   !          - parameters XX=0., X1=.14
   ! 019      B. Bilodeau (May 1994) - New physics interface
   ! 020      G. Pellerin (Nov 1994) - New surface layer formulation
   ! 021      G. Pellerin (Jan 1995) - Modifier l'extraction de LE
   ! 022      G. Pellerin (Jun 1995) - Revert to original rigrad for
   !                          computation of unstable boundary layer
   ! 023      B. Bilodeau (Nov 1995) - Replace VK by KARMAN
   ! 024      B. Bilodeau and J. Mailhot (Jan 1996) -
   !          Eliminate divisions by zero
   ! 024      R. Sarrazin (Jan 1996) - Carry boundary layer pointer in KCL
   ! 025      C. Girard (Fev 1996) - Introduce different options for
   !             shalow convection - GELEYN,CONRES,SHALOW
   ! 026      A-M.Leduc (Sept 2002) - add QC in arguments and remove ISHLCVT
   !                                  eturbl4--->eturbl5.
   !                                  Add X1 calculation for call to mixlen1.
   ! 027      J. Mailhot (Mar 2003) - TKE advection.
   ! 028      A. Plante (June 2003) - IBM conversion
   !             - Replace call to CVMG* by if-else statements
   !             - call to exponen4 (to calculate power function '**')
   !             - call to vslog routine (from massvp4 library)
   !             - constants precomputations
   !             - @PROCESS STRICT compilation option added
   ! 029      S. Belair  (Mar 2003) - Add F(ZD) in arguments ...> eturbl6
   !             - Use time filter for the Bougeault-Lacarrere mixing length.
   !
   ! 030      A. Plante (July 2003) - Correct bug in IBM conversion :
   !             - Virtual temperature was incorrect for MIXLEN1. This was
   !               a problem only if longmel = BOUJO
   ! 031      B. Bilodeau (Aug 2003) - exponen4 replaced by vspow
   !                                   call to mixlen2
   ! 032      A-M. Leduc (March 2004) - add arguments S ans PS to MIXLEN2--->
   !                                    MIXLEN3
   ! 033      Y. Delage (Sept 2004) - Introduce log-linear profile mixing length for
   !                                    near-neutral cases.  Optimisation of KT.
   ! 034      L. Spacek (Dec 2007) - add "vertical staggering" option
   !                                 change the name to eturbl7
   ! 035      L. Spacek    (Sep 2011)   Eliminate obsolete options

   !@Arguments
   !          - Input/Output -
   ! EN       turbulent energy
   ! ZN       mixing length of the turbulence
   !
   !          - Output -
   ! RIF      flux Richardson number
   ! RIG      gradient Richardson number
   ! SHR2     square of wind shear
   !
   !          - Input -
   ! ENOLD    turbulent energy (at time -)
   ! ZNOLD    mixing length of the turbulence (at time -)
   ! GAMA     countergradient term in the transport coefficient of
   !          Theta and Q
   ! HOL      inverse of length of Monin-Obokhov
   ! FN       cloud fraction
   ! U        east-west component of wind
   ! V        north-south component of wind
   ! T        temperature
   ! TVE      virtual temperature on 'E' levels
   ! Q        specific humidity
   ! QC       cloud water
   ! QE       specific humidity on 'E' levels
   ! XH       convective velocity scale (w*)
   !          - Input -
   ! PS
   ! S        sigma level
   ! SE       sigma level for turbulent energy
   ! DSGDZ    sigma intervals
   ! TAU      timestep
   ! KOUNT    index of timestep
   ! KT       ratio of KT on KM (real KT calculated in DIFVRAD)
   ! Z        height of sigma level
   ! GZMOM    height of sigma momentum levels
   !          - Input -
   ! TRNCH    number of the slice
   ! N        horizontal dimension
   ! M        1st dimension of T, Q, U, V
   ! NK       vertical dimension
   ! Z0       roughness length

   !@Notes
   !          EN and ZN contain the values at time T and input U
   !          and V are wind images. C and ZE are over-written.
   !          Refer to J.Mailhot and R.Benoit JAS 39 (1982)Pg2249-2266
   !          and Master thesis of J.Mailhot.
   !*@/

#include "clefcon.cdk"
#include "surface.cdk"
#include "machcon.cdk"
   include "phyinput.inc"

   real, dimension(N) :: XB
   real, dimension(N,2) :: WK
   real, dimension(N,NK) :: WORK,ZE,C,X,DSGDZ,X1,FM,FH
   real, dimension(N,4*NK) :: B

   !***********************************************************************

   !     temporary variables used to convert a #@$%!& CVMG.. expression

   real yuk1,yuk2,dtfac, r4tmp

   real, parameter :: EPSILON_B=1.E-8,PETIT=1.E-6,LMDA=200.
   real ZNOLD(N,NK),beta_sfc(n)
   real SC,EXP_EXPLIM,TAUINV
   real, dimension(n,nk) :: zn_blac,zn_boujo,blend_hght
   real, dimension(n,nk,3) :: w_cld
   integer J,K,STAT
   integer NKE
   integer, dimension(N) :: SLK
   integer, external :: neark

   EXP_EXPLIM=exp(-EXPLIM)
   TAUINV=1.0/TAU
   w_cld = 0.

   NKE=NK

   !     SHR2  =  ( D VENT / D Z ) ** 2

   call ABSDVDZ3(SHR2,U,V,TVE,SE,DSGDZ,S,N,N,NK)

   do k=1,NKE
      do j=1,N
         SHR2(j,k) = SHR2(j,k) + PETIT
      end do
   end do

   !     RIG ( NOMBRE DE RICHARDSON GRADIENT)

   if (PBL_RIBKG) then
      call rigrad3(shr2,tve,qe,zn,ps,se,z,rig,gama,gamaq,PBL_RIBKG,n,nk)
   else
      call rigrad1b(rig, gama, gamaq, xb, shr2, t, tve, q, qe, &
           s, wk, n, n, nk)
   endif

   call series_xst(RIG, 'RI', trnch)

   !           AJOUT DE L'EFFET DE LA CONVECTION RESTREINTE

   do k=1,NKE
      do j=1,N
         FN(j,k) = 0.
         GAMAL(j,k) = 0.
         WORK(j,k) = 0.
      end do
   end do

   if( PBL_SHAL == 'CONRES' ) then
      call CONRES1(RIG,GAMA,GAMAQ,WORK,T,TVE,Q,QE,PS, &
           HOL,S,SE,SHR2,WK,N,N,NK)
   endif

   call series_xst(RIG, 'RM', trnch)


   !                               CALCUL DE LA LONGUEUR DE MELANGE


   ZNOLD(:,:) = ZN(:,:)


   !                               A) BLACKADAR (1962)

   ! Compute the PBL stability functions
   stat = psf_stabfunc(rig,z,fm,fh, &
        blend_bottom=pbl_slblend_layer(1),blend_top=pbl_slblend_layer(2))
   kt = fm/fh  !temporary storage for the inverse Prandtl number
   if (stat /= PHY_OK) then
      call physeterror('eturbl', 'error returned by PBL stability functions')
      return
   endif

   ! Compute the Blackadar (1962) mixing length,
   ! including high-resolution (<850m) reduction of asymtotic value
   do k=1,nke
      do j=1,n
         lmn = min(KARMAN*(z(j,k)+z0(j)),min(0.23*sqrt(dxdy(j)),LMDA))
         zn_blac(j,k) = lmn / fm(j,k)
      enddo
   enddo

   !     RIF ( NOMBRE DE RICHARDSON DE FLUX)

   RIF = RIG
   call RIFLUX2(RIF,KT,N,NK)

   !     APPLY HYSTERESIS BASED ON TURBULENCE REGIME
   stat = neark(se,ps,3000.,n,nk,slk) !determine "surface layer" vertical index
   if (kount == 0) then
      INIT_TURB: if (.not.any('turbreg'==phyinread_list_s(1:phyinread_n))) then
         do k=1,nk
            do j=1,n
               if (k <= slk(j)) then
                  if (RIF(j,k) > PBL_RICRIT(1)) then
                     TURBREG(j,k) = LAMINAR
                  else
                     TURBREG(j,k) = TURBULENT
                  endif
               else
                  TURBREG(j,k) = TURBULENT
               endif
            enddo
         enddo
      endif INIT_TURB
   endif
   do k=1,nk             !do not apply to the lowest levels
      if (p_prof(k) < 60000.) cycle
      do j=1,n
         ABOVE_SFCLAYER: if (k <= slk(j)) then
            if (RIF(j,k) < PBL_RICRIT(1)) then
               TURBREG(j,k) = TURBULENT
            elseif (RIF(j,k) > PBL_RICRIT(2)) then
               TURBREG(j,k) = LAMINAR
            endif
            if (RIF(j,k) > PBL_RICRIT(1) .and. nint(TURBREG(j,k)) == LAMINAR) RIF(j,k) = max(RIF(j,k),1.)
            if (RIF(j,k) < PBL_RICRIT(2) .and. nint(TURBREG(j,k)) == TURBULENT) RIF(j,k) = min(RIF(j,k),1.)
         endif ABOVE_SFCLAYER
      enddo
   enddo

   call series_xst(RIF, 'RF', trnch)

   !                                Calculate the mixing length
   !                                according to Bougeault and
   !                                Lacarrere (1989)

   IF_BOUJO: if (any(longmel == (/'BOUJO   ', 'TURBOUJO'/))) then
      do K=1,NK
         do J=1,N
            !# Virtual potential temperature (THV)
            X1(J,K)=TE(J,K)*(1.0+DELTA*QE(J,K)-QCE(J,K))* (SE(J,K)**(-CAPPA))
         end do
      end do
      if (ml_calc_boujo(zn_boujo, x1, enold, w_cld, z, se, ps) /= PHY_OK) then
         call physeterror('eturbl', 'error returned by B-L mixing length estimate')
         return
      endif
      blend_hght(:,:nk-1) = gzmom(:,2:)
      blend_hght(:,nk) = 0.
      zn_boujo = max(zn_boujo,1.)
      if (ml_blend(zn, zn_blac, zn_boujo, blend_hght, se, ps) /= PHY_OK) then
         call physeterror('eturbl', 'error returned by mixing length blending')
         return
      endif
      if (longmel == 'TURBOUJO') then
         where (nint(turbreg) == LAMINAR) zn = zn_blac ! Use local (Blackadar) mixing length for laminar flows
      endif
   else
      zn = zn_blac
   endif IF_BOUJO

   !     No time filtering at kount=0 since this step is used for initialization only
   if(KOUNT == 0)then
      if (any('zn'==phyinread_list_s(1:phyinread_n))) zn(1:n,1:nke) = znold(1:n,1:nke)
   else
      IF_ZN_TIME_RELAX: if (PBL_ZNTAU > 0.) then
         EXP_TAU = exp(-TAU/PBL_ZNTAU)
         do K=1,NKE
            do J=1,N
               ZN(J,K) = ZN(J,K) + (ZNOLD(J,K)-ZN(J,K))*EXP_TAU
            end do
         end do
      endif IF_ZN_TIME_RELAX
   endif

   if (any(longmel == (/'BOUJO   ', 'TURBOUJO'/))) then
      do K=1,NKE
         do J=1,N
            ZE(J,K) = ZN(J,K) * ( 1. - min( RIF(J,K) , 0.4) ) &
                 / ( 1. - 2.*min( RIF(J,K) , 0.4) )
            ZE(J,K) = max ( ZE(J,K) , 1.E-6 )
         end do
      end do
   else
      do K=1,NKE
         do J=1,N
            ZE(J,K)=max(ZN(J,K),1.E-6)
         end do
      end do
   end if

   if (PBL_DISS == 'LIM50') ze(:,1:nke) = min(ze(:,1:nke),50.)

   ZD(:,1:NKE) = ZE(:,1:NKE)

   call series_xst(ZN, 'L1', trnch)
   call series_xst(ZE, 'L2', trnch)

   call series_xst(ZE, 'LE', trnch)

   do K=1,NKE
      do J=1,N
         ZE(J,K)=1./ZE(J,K)
      enddo
   enddo


   !     CALCUL DES TERMES ALGEBRIQUES DE L'EQUATION
   !             DE L'ENERGIE TURBULENTE

   ! B       - PRODUCTION MECANIQUE ET THERMIQUE DE L'ENERGIE TURBULENTE
   !           PEUT ETRE NEGATIVE OU POSITIVE
   ! C       - DISSIPATION VISQUEUSE DE L'ENERGIE TURBULENTE > 0.0

   do K=1,NKE
      do J=1,N
         C(J,K)=BLCONST_CE*ZE(J,K)

         ZE(J,K)=(SHR2(J,K)-PETIT)*ZN(J,K)*BLCONST_CK/(C(J,K)+PETIT)
         X(J,K)=-SHR2(J,K)*RIF(J,K)*ZN(J,K)*BLCONST_CK/(C(J,K)+PETIT)
         B(J,K)=(C(J,K)+PETIT)*(ZE(J,K)+X(J,K))
      enddo
   enddo

   if(KOUNT.eq.0)then

      !     SOLUTION STATIONNAIRE

      !     ON INITIALISE EN
      !     STATION EST ENLEVE  (+PRECALCUL DE EN)

      do K=1,NK
         do J=1,N
            X(J,K)=0.0
         enddo
      enddo

      call series_xst(X, 'EM', trnch)
      call series_xst(X, 'EB', trnch)
      call series_xst(X, 'ED', trnch)
      call series_xst(X, 'ET', trnch)
      call series_xst(X, 'ER', trnch)

   else

      !     SOLUTION DE LA PARTIE ALGEBRIQUE DE L'EQUATION
      !               DE L'ENERGIE TURBULENTE

      DO4: do K=1,NKE
         do J=1,N

            if(abs(B(J,K)) < EPSILON_B) then
               C(J,K)= C(J,K)*TAU
               B(J,K)=0.
            else
               C(J,K)= sqrt(abs(C(J,K)/B(J,K)))
               B(J,K)=min(B(J,K)*C(J,K)*TAU,EXPLIM)
            endif

            if(B(J,K) > epsilon(B)) then
               yuk1 = -1.0+2.0/(1.0+exp(-AMIN1(abs(B(J,K)),174.))*(-1.0+2.0/ &
                    (1.0+sqrt(EN(J,K))*C(J,K))))
               B(J,K) = yuk1
            else if(B(J,K) < -epsilon(B)) then
               yuk1 = tan(min(TANLIM,max( atan( sqrt(EN(J,K))*C(J,K) ) &
                    +0.5*B(J,K) , 0.0 ) ) )
               B(J,K) = yuk1
            else
               yuk1 = sqrt(EN(J,K))*C(J,K)/(1.+0.5*sqrt(EN(J,K))*C(J,K))
               B(J,K) = yuk1
            endif

            if(abs(C(J,K)) < epsilon(C)) then
               B(J,K)=EN(J,K)
            else
               B(J,K)=(B(J,K)/C(J,K))**2
            endif
            if(B(J,K)-PETIT .lt. 0.) B(J,K)=ETRMIN
            C(J,K)=ZE(J,K)+X(J,K)

            !              TERMES DE PRODUCTION MECANIQUE ET THERMIQUE NULS SI EN=C
            if ((EN(J,K)-C(J,K)).ne.0.0) then
               yuk2=abs((B(J,K)-C(J,K)) / (EN(J,K)-C(J,K)))
               r4tmp = log(max(yuk2,exp_explim))
            else
               r4tmp = 0.  !# log(1.0)
            endif

            !     TERME DE PRODUCTION MECANIQUE

            ZE(J,K)=-ZE(J,K)*r4tmp *tauinv

            !     TERME DE PRODUCTION THERMIQUE

            X(J,K)=-X(J,K)  *r4tmp *tauinv

            !     TERME DE DISSIPATION VISQUEUSE
            C(J,K)=-X(J,K)-ZE(J,K)+(B(J,K)-EN(J,K))*TAUINV
         enddo
      enddo DO4

      do J=1,N
         ZE(J,NK)=0.0
         X (J,NK)=0.0
         C (J,NK)=0.0
         EN(J,NK)=0.0
      enddo

      call series_xst(ZE, 'EM', trnch)
      call series_xst(X,  'EB', trnch)
      call series_xst(C,  'ED', trnch)

      !     SOLUTION DE LA PARTIE DIFFUSIVE (E-F) DE L'EQUATION
      !                   DE L'ENERGIE TURBULENTE

      do K=1,NK
         do J=1,N
            SC=pbl_tkediff*CLEFAE*BLCONST_CK
            !     (E*-EN)/TAU
            X(J,K)=X(J,K)+C(J,K)+ZE(J,K)
            ZE(J,K)=B(J,K)
            C(J,K)=SC*ZN(J,K)*sqrt(ENOLD(J,K))*DSGDZ(j,k)**2
         enddo
      enddo

      !     METTRE X1 A ZERO
      do K=1,NK
         do J=1,N
            X1(J,K)=0
         enddo
      enddo
      !     C CONTIENT K(E) ET X1 CONTIENT ZERO

      if (pbl_zerobc) then
         ZE(:,NK) = 0.
         C(:,NK) = 0.
         beta_sfc = 0.
         XB = 0.
      else
         XB = BLCONST_CU*FRV**2 + BLCONST_CW*XH**2
         ZE(:,NK) = XB
         beta_sfc = XH
      endif
      dtfac = 1.
      if (pbl_tkediff2dt) dtfac = 2.

      call DIFUVDFj1 (EN,ZE,C,X1,X1,X1,XB,beta_sfc,S,SE,dtfac*TAU,4,1.,N,N,N,NKE)
      if (phy_error_L) return

      !     NOUVEAU EN
      en = max(ETRMIN,ze+tau*en)

   endif

   return
end subroutine ETURBL13
