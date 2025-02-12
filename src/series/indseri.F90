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
