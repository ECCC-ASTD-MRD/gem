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
!**FUNCTION INDSERI - RANG DU NOM DANS LA LISTE
!
      FUNCTION INDSERI ( NOM , LISTE , N )
      implicit none
!!!#include <arch_specific.hf>
      character(len=*) NOM, LISTE(*)
      INTEGER INDSERI, N
      INTEGER I
!
!Author
!          R. Benoit
!
!Revision
! 001      V.Alex.(Feb 87)- Change in loop 10, documentation
! 002      B. Bilodeau  (July 1991)- Adaptation to UNIX
!
!Object
!          to return the index in the list of names(LISTE) pointing
!          to the name(NOM)
!
!Arguments
!
!          - Input -
! NOM      name to look for in the LISTE
! LISTE    list of names
! N        length of LISTE
!
!*
      INDSERI = 0
      DO 10 I=1,N
         IF ( NOM.EQ.LISTE(I) ) THEN
           INDSERI=I
           RETURN
         ENDIF
   10 CONTINUE
      RETURN
      END
