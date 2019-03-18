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

      SUBROUTINE NUAGES2 ( CH , CM , CL , C3D , &
                           BASE, Q , T , PS, SHCL, ILMO, S, &
                           TRNCH, N, M, NK, ITASK, SATUCO, STRCLD)
      use tdpack
      use series_mod, only: series_xst
      implicit none
#include <arch_specific.hf>
      INTEGER N,M,NK,ITASK,IERROR
      REAL CH(N),CM(N),CL(N),C3D(N,NK),Q(M,NK)
      REAL T(M,NK),PS(N),SHCL(N),ILMO(N),S(n,NK)
      REAL BASE (N)
      INTEGER TRNCH
      LOGICAL SATUCO,STRCLD

!@Author R. Benoit RPN (April 1984)
!@Revision
! 001      J. Cote RPN(Nov 1984)SEF version documentation
! 002      R. Benoit RPN(April 1985) Remove clouds in unstable
!          boundary layer
! 003      R Benoit (May85) Inverse of calculations and base layer of
!          condensation
! 004      M. Lepine  -  RFE model code revision project (Feb 87)
! 005      J.Mailhot(Mar 87)Base threshold of COND= KUO
! 006      G.Pellerin(Nov 87) Adaptation to revised code
! 007      G.Pellerin(Jan 90) Adaptation to version D4P6
! 008      G.Pellerin(Mar 90) Adaptation to version D5P7
! 009      N. Brunet (May 90) Standardization of thermodynamic
!          functions
! 010      Y. Delage (Nov 1990) Replace WC by ILMO
! 011      C. Girard(Nov 90)
!               Substantial modification to the CLOUD parameter
! 012      N. Brunet  (May91)
!               New version of thermodynamic functions
!               and file of constants
! 013      B. Bilodeau  (July 1991)- Adaptation to UNIX
! 014      C. Girard (Nov 1992) New modification to the
!          definition
! 015      R. Benoit (Aug 93) Local Sigma
! 016      B. Bilodeau (May 1994) - New physics interface
! 017      R. Sarrazin (Summer 1995) - Correct bug for CL
! 018      B. Dugas (Sep 1996) - Add option to eliminate
!          stratospheric clouds
! 019      J.P. Toviessi (June 2003) - Revove CVMG functions
!
!@Object calculate simplified cloud cover
!
!@Arguments
!          - Output -
! CH       high altitude cloud fraction (0 to 1)
! CM       medium altitude cloud fraction (0 to 1)
! CL       low altitude cloud fraction (0 to 1)
! C3D      3-dimensional cloud field
! BASE     sigma base of condensation layer (+/-)
!
!          - Input -
! Q        specific humidity
! T        temperature
! PS       surface pressure
! SHCL     sigma height of the boundary layer
! ILMO     inverse of the length of Monin-Obukhov
! S        sigma levels
! TRNCH    index of the vertical plane(NI*NK) for which
!          calculations are to be done
! N        horizontal dimension
! M        1st dimension of T and Q
! NK       vertical dimension
! ITASK    task number for multi-tasking
! SATUCO   .TRUE. if water/ice phase for saturation
!          .FALSE. if water phase only for saturation
! STRCLD   .TRUE. if stratospheric clouds are acceptable
!          .FALSE. otherwise

#include "clefcon.cdk"

      REAL SST,SCL,SH,SM
      REAL F,SIG,NEBUL,U,SQRT3
      INTEGER J,K,KH,KM,KL
      INTEGER K1
      LOGICAL REGULAR
      logical flag
      REAL EPS , VAL , HM
      SAVE SST,SCL,SH,SM
!
!***********************************************************************
!     AUTOMATIC ARRAYS
!***********************************************************************
!
      REAL, dimension(NK) :: ftmp
      REAL, dimension(N,NK) :: ftmp1
!
!***********************************************************************
!
      REAL tmpNEBUL,tmp
!
      DATA SST , SCL , SH , SM / 0.225 , 0.905 , 0.395 , 0.710 /
!
!     TOPC and MINQ are the minimum values of pressure and
!     specific humidity at which the routine stops producing
!     clouds when STRCLD is set to .FALSE.
#include "nocld.cdk"
!
!     ANCIENNE FORMULATION (REVISION 10)
!     F(SIG)=MIN(.98,MAX(.8,(2.+SIG)/3.))
!     NEBUL(U,SIG)=MAX(0.0,MIN(1.0,(U-F(SIG))/(1.0-F(SIG))))**2
!
!     NOUVELLE FORMULATION (REVISION 12)
      tmpNEBUL(U,F)=MAX(0.0,MIN(1.0,(U-F)/((1.0-F))))
      NEBUL(U,F)=tmpNEBUL(U,F)*tmpNEBUL(U,F)
!
!
      SQRT3=SQRT(3.)
!
!
!  FAIRE D'ABORD U (=HUM.REL.) DANS C3D
!
      do K = 1, NK
        ftmp(K) = S(1,K)
        do J = 1,N
         ftmp1(J,K) = S(J,K)*PS(J)
        enddo
      enddo

      call mfohr4 ( C3D(1,1),Q(1,1),T(1,1),ftmp1,N,NK,N,SATUCO)
      call series_xst(c3d, 'hr', trnch)
!
!  BASE DE COUCHE DE CONDENSATION (SI EXISTE)
!  HM = SEUIL DE DEBUT DE CONDENSATION
!
      HM=0.9
!
!  TROUVER 1ER NIVEAU EN MONTANT OU U>HM . METTRE DANS CL
!
      DO 2 J = 1 ,N
2     CL (J) = 0.
!
      DO  K = NK-1 , 1 , -1
         DO  J = 1 , N
           if (C3D(J,K).GT.HM .AND. CL(J).EQ.0.) CL (J) = FLOAT(K)
         ENDDO
      ENDDO
!
      EPS = 1.E-12
      DO 4 J = 1 , N
         K1 = NINT (CL(J) )
         BASE (J) = S (j,1)

         if (K1.EQ.NK) BASE (J) = S(j,NK)
         REGULAR = K1.GE.1 .AND. K1.LT.NK
         K1 = MIN (NK-1 , MAX (K1,1) )
         VAL = S(j,K1) + (S(j,K1+1)-S(j,K1)) * (C3D(J,K1)-HM) &
                 / MAX ( EPS , C3D(J,K1)-C3D(J,K1+1))
         VAL = MIN ( VAL , S(j,K1+1) )
         if (REGULAR) BASE (J) = VAL
!
!  DONNER SIGNE A BASE . + SI W*=0 , - SI W*>0 .
!
           BASE(J) = SIGN( BASE(J), ILMO(J) )
4     CONTINUE
!
      DO 5 K=1,NK
         DO 5 J=1,N
!
!  LES NUAGES HORS DE (SST,SCL) NE SONT PAS OTES DE C3D.
!  SI STRCLD EST FAUX, IL N'Y A PAS DE NUAGES AU-DESSUS
!  DE TOPC OU BIEN SI Q EST PLUS PETIT QUE MINQ.
!                    ilmo .lt.0  instable
       if (T(J,K).LT.238.) then
           F= .95
       else
           F= .99
       endif
       if ( S(j,K).GE.SHCL(J).and. ilmo(j) .lt.0.) then
          C3D(J,K) = 0.
       else
          C3D(J,K) =  NEBUL(C3D(J,K),F)
       endif
       if (.NOT.STRCLD.AND.(S(J,K)*PS(J).LT.TOPC.OR.Q(J,K).LE.MINQ)) &
            C3D(J,K) = 0.
5     CONTINUE
!
      DO 10 J=1,N
         CH(J)=0.0
         CM(J)=0.0
   10    CL(J)=0.0
!
      do 50 k=1,nk
         do  j=1,n
            flag = (s(j,k).lt.shcl(j).or.ilmo(j).ge.0.)
            if (s(j,k).ge.sst.and.s(j,k).lt.sh) then
               if (flag) ch(j)= max(c3d(j,k) , ch(j))
            elseif (s(j,k).ge.sst.and.s(j,k).lt.sm) then
               if (flag) cm(j)= max (c3d(j,k) , cm(j))
               elseif (s(j,k).ge.sst.and.s(j,k).lt.scl) then
               if (flag) cl(j)=  max (c3d(j,k) , cl(j))
            endif
         enddo

 50   continue
!
      return
      end
