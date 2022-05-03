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
!**S/P utco2
!
      real function utCO2 (ii,jj,ix,nx,NN,NK,SH,S,SC,DEL)
!
      implicit none
!!!#include <arch_specific.hf>
      integer ii, jj, nx, ix
      INTEGER NN,NK
      REAL DEL(NK),SC(NN),S(NN),SH(nx,NK)
!
!Author
!          R. Benoit (Aug 93) -- from L Garand's CO2INFO s/r
!
!Object
!          to calculate the quantities of CO2 and the
!          transmissivity from level ii to level jj.  Local sigma.
!          calculation done for an entire line (nx)
!
!Arguments
!
!          - Output -
! utco2    depending upon value of jj, the function returns:
!          UCO2  if jj=0
!          TCO2  if jj>0
!
! UCO2     amount of CO2 in each layer of thickness DEL (multiply by
!          PS**2 to get kg/metres squared)
! TCO2     precalculated transmissivity of CO2 from level to level
!          (EXP(-PS*TCO2) = matrix of transmission.  The upper
!          triangle of TCO2 is used for the (strong) central band.
!          The lower triangle of TCO2 is used for the average of the
!          right and left wings)
!
!          - Input -
! ii       index of first  level
! jj       index of second level
! ix       index of single sigma profile to use in sh
! nx       horizontal dimension of sh
! NN       number of levels (NK+1)
! NK       number of layers
! SH       sigma levels at the centre of layers
! S        work space (sigma levels at the borders of the layers)
! SC       work space
! DEL      work space (sigma thickness from level to level)
!
!PARAMETERS
!
      REAL A1D,A1G,A2D,A2G,AWING,ECO2,A1C,A2C,QCO2
      PARAMETER  (ECO2=1.00)
      PARAMETER (A1C=198.0)
      PARAMETER (A2C=0.9768)
      PARAMETER (QCO2=5.01E-4)
      PARAMETER (A1D=4.035)
      PARAMETER (A2D=0.8224)
      PARAMETER (A1G=5.439)
      PARAMETER (A2G=0.9239)
      parameter (AWING= (A1G*A2G + A1D*A2D)/2. )
!  A1D ET A2D SONT LES PARAMETRES DE L'AILE DROITE DU CO2
!  A1G ET A2G """"""""""""""""""""""""""""" GAUCHE """""""
!  A1C ET A2C """"""""""""""""""" DE LA BANDE CENTRALE (FORTE) DU CO2
!*
!     AWING= (A1G*A2G + A1D*A2D)/2.
!     PARAMETRE D'ABSORTION MOYEN POUR LES DEUX AILES
      integer I,JK,L,J,itrapez
      REAL ELSA,Z,ZU,TRAPEZ2
!
      SC(NN)=QCO2
      S(NN)=1.
      S(1)=2.*SH(ix,1)-((SH(ix,1)+SH(ix,2))/2.)
!  CETTE DEFINITION DU PREMIER NIVEAU DE FLUX DOIT ETRE LA MEME
!  QUE DANS LE CODE DE RADIATION
      S(1)=AMAX1(S(1),0.0003)
      S(NN)=1.
      DO 25 I=2,NK
      S(I)=(SH(ix,I)+SH(ix,I-1))/2.
      DEL(I-1)=S(I)-S(I-1)
25    CONTINUE
      DEL(NK)=1.-S(NK)
      DO 10 I=1,NK
      SC(I)=QCO2*S(I)**ECO2
  10  CONTINUE
      ELSA=1.66
      Z=1./(101325.*9.80616)
      i=ii
      j=jj
!
!     commented lines below are from s/r co2info
!
!     DO 79 I=1,NN
!     TCO2(I,I)=1.
      uTCO2=1.
      if (i.eq.j) return
      JK=I+1
      if (j.eq.0) j=jk
!     IF(I.EQ.NN) return
!     DO 79 J=JK,NN
!     L=J-I+1
      L=abs(J-I)+1
      i=min(i,j)
!     -----------------------------------------
!     REAL FUNCTION TRAPEZ2(DEL,F,N,NM)
!     TRAPEZ2=0.
!     DO 10 I=1,NM
! 10  TRAPEZ2=TRAPEZ2+(F(I)+F(I+1))/2.*DEL(I)
!     -----------------------------------------
      trapez2=0.
      do 11 itrapez=1,l-1
         trapez2=trapez2+ &
              (sc(i+itrapez-1)+sc(i+itrapez))*del(i+itrapez-1)
 11   continue
!     ZU=Z*TRAPEZ2(DEL(I),SC(I),L,L-1)/2.*ELSA
      ZU=Z*TRAPEZ2/2.*ELSA
      if (jj.eq.0) then
         utco2=zu
      elseif (jj.gt.ii) then
!     TCO2(I,J)=SQRT(A1C*A2C*ZU)
         utco2=SQRT(A1C*A2C*ZU)
      elseif (jj.lt.ii) then
!     TCO2(J,I)=SQRT(AWING*ZU)
         utco2=SQRT(AWING*ZU)
      endif
! 79  IF(J.EQ.JK)UCO2(I)=ZU
! 42  CONTINUE
      RETURN
      END
