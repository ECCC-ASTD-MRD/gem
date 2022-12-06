      SUBROUTINE XIT(NAME,N)
! 
!     * OCT 01/92 - E.CHAN. (CHANGE STOP 1 TO STOP)
!     * JUN 10/91 - E.CHAN. (TRANSLATE HOLLERITH LITERALS AND 
!     *                      DIMENSION STRINGS) 
! 
!     * OCT 10/78 - J.D.HENDERSON.
!     * TERMINATES A PROGRAM BY PRINTING THE PROGRAM NAME AND 
!     * A LINE ACROSS THE PAGE FOLLOWED BY A NUMBER N. 
! 
!     * N.GE.0 IS FOR A NORMAL END. THE LINE IS DASHED. 
!     * NORMAL ENDS TERMINATE WITH   STOP.
! 
!     * N.LT.0 IS FOR AN ABNORMAL END. THE LINE IS DOTTED.
!     * IF N IS LESS THAN -100 THE PROGRAM SIMPLY TERMINATES. 
!     * OTHERWISE IF N IS LESS THAN ZERO THE PROGRAM ABORTS.
! 
      CHARACTER*(*) NAME
      CHARACTER*8   NAME8, DASH, STAR
! 
      DATA DASH /'--------'/, STAR /'********'/ 
!---------------------------------------------------------------------
! 
      NAME8 = NAME
      IF(N.GE.0) WRITE(6,6010) DASH,NAME8,(DASH,I=1,9),N 
! 
      IF(N.LT.0) WRITE(6,6010) STAR,NAME8,(STAR,I=1,9),N 
! 
      IF ( N.GE.0 .OR. N.LT.-100 ) THEN
	STOP 
      ELSE
	CALL ABORT
      ENDIF
! 
!---------------------------------------------------------------------
 6010 FORMAT('0',A8,'  END  ',A8,9A8,I8)
      END   
