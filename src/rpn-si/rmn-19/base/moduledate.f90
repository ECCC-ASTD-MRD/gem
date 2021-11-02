!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
!============================================================================
!                       THREAD SAFE ROUTINES 
!     (they call the original ones inside a OpenMP critical region)
!     the original routine names have been deliberately mangled
!============================================================================
!     M.Valin 2016/12/08  initial release of the thread safe cover routines
!
      subroutine date_thread_lock(lock)  ! .true. attempt to acquire lock, .false. release lock
      IMPLICIT NONE
      logical, intent(IN) :: lock
      integer, save :: owner_thread = 0
      call set_user_lock(owner_thread,lock)
      return
      end subroutine date_thread_lock

      subroutine INCDATi(idate1,idate2,nhours)
      IMPLICIT NONE
      integer :: idate1,idate2
      real *8 :: nhours
      call date_thread_lock(.true.)
      call IDNACTi(idate1,idate2,nhours)
      call date_thread_lock(.false.)
      end subroutine INCDATi

      subroutine INCDATr(idate1,idate2,nhours)
      IMPLICIT NONE
      integer :: idate1,idate2
      real *8 :: nhours
      call date_thread_lock(.true.)
      call IDNACTr(idate1,idate2,nhours)
      call date_thread_lock(.false.)
      end subroutine INCDATr

      subroutine DIFDATi(idate1,idate2,nhours)
      IMPLICIT NONE
      integer :: idate1,idate2
      real *8 :: nhours
      call date_thread_lock(.true.)
      call DDIAFTi(idate1,idate2,nhours)
      call date_thread_lock(.false.)
      end subroutine DIFDATi

      subroutine DIFDATr(idate1,idate2,nhours)
      IMPLICIT NONE
      integer :: idate1,idate2
      real *8 :: nhours
      call date_thread_lock(.true.)
      call DDIAFTr(idate1,idate2,nhours)
      call date_thread_lock(.false.)
      end subroutine DIFDATr

      INTEGER FUNCTION newdate(DAT1,DAT2,DAT3,MODE)
      IMPLICIT NONE
      integer :: DAT1,DAT2(*),DAT3,MODE
      integer, external :: naetwed
      call date_thread_lock(.true.)
      newdate = naetwed(DAT1,DAT2,DAT3,MODE)
      call date_thread_lock(.false.)
      end FUNCTION newdate

      INTEGER FUNCTION IDATMG2(IDATE)
      IMPLICIT NONE
      integer idate(14)
      integer, external :: itdmag2
      call date_thread_lock(.true.)
      IDATMG2 = itdmag2(IDATE)
      call date_thread_lock(.false.)
      end function IDATMG2

      subroutine DATMGP2(IDATE)
      IMPLICIT NONE
      integer idate(14)
      call date_thread_lock(.true.)
      call dmagtp2(IDATE)
      call date_thread_lock(.false.)
      end subroutine DATMGP2

      subroutine NewDate_Options( value,command )
      IMPLICIT NONE
      character*(*) value,command
      call date_thread_lock(.true.)
      call NewDate_Options_int( value,command )
      call date_thread_lock(.false.)
      end subroutine NewDate_Options

      subroutine Get_Calendar_Status( NoLeapYears,CcclxDays )
      IMPLICIT NONE
      logical ::  NoLeapYears,CcclxDays
      call date_thread_lock(.true.)
      call Get_Calendar_Status_int( NoLeapYears,CcclxDays )
      call date_thread_lock(.false.)
      end subroutine Get_Calendar_Status

      integer function Calendar_Adjust(tdate1,tdate2,true_date_mode,adding)
      IMPLICIT NONE
      integer, external :: Calendar_Adjust_int
      integer :: tdate1,tdate2
      character(len=1) true_date_mode
      logical :: adding
      call date_thread_lock(.true.)
      Calendar_Adjust = Calendar_Adjust_int(tdate1,tdate2,true_date_mode,adding)
      call date_thread_lock(.false.)
      end function Calendar_Adjust
                                       
      integer function CcclxDays_Adjust(tdate1,tdate2,true_date_mode,adding)
      IMPLICIT NONE
      integer, external :: CcclxDays_Adjust_int
      integer :: tdate1,tdate2 ! input TrueDates
      character(len=1) true_date_mode ! (B)asic or (E)xtended TrueDates
      logical :: adding ! operating mode (T=incadtr, F=difdatr)
      call date_thread_lock(.true.)
      CcclxDays_Adjust = CcclxDays_Adjust_int(tdate1,tdate2,true_date_mode,adding)
      call date_thread_lock(.false.)
      end function CcclxDays_Adjust

      integer function LeapYear_Adjust(tdate1,tdate2,true_date_mode,adding)
      IMPLICIT NONE
      integer, external :: LeapYear_Adjust_int
      logical :: adding
      character(len=1) true_date_mode ! (B)asic or (E)xtended true dates
      integer :: tdate1,tdate2
      call date_thread_lock(.true.)
      LeapYear_Adjust = LeapYear_Adjust_int(tdate1,tdate2,true_date_mode,adding)
      call date_thread_lock(.false.)
      end function LeapYear_Adjust

      subroutine Ignore_LeapYear()
      IMPLICIT NONE
      call date_thread_lock(.true.)
      call Ignore_LeapYear_int
      call date_thread_lock(.false.)
      end subroutine Ignore_LeapYear

      subroutine Accept_LeapYear()
      IMPLICIT NONE
      call date_thread_lock(.true.)
      call Accept_LeapYear_int
      call date_thread_lock(.false.)
      end subroutine Accept_LeapYear

      subroutine Get_LeapYear_Status(no_leap_year_status)
      IMPLICIT NONE
      logical :: no_leap_year_status
      call date_thread_lock(.true.)
      call Get_LeapYear_Status_int(no_leap_year_status)
      call date_thread_lock(.false.)
      end subroutine Get_LeapYear_Status
!============================================================================
!     END OF THREAD SAFE ROUTINES
!============================================================================
!   the original names of the following routines have been altered because of
!   the above mentioned thread safe routines
!   internal calls use the mangled internal names
!============================================================================
!    environment variable NEWDATE_OPTIONS usage syntax
!
!   export NEWDATE_OPTIONS="[debug][,][year=360_day|365_day|gregorian][,][debug]"
!
!   examples of usage:
!
!   export NEWDATE_OPTIONS="debug"
!   export NEWDATE_OPTIONS="year=360_day"
!   export NEWDATE_OPTIONS="debug,year=360_day"
!   export NEWDATE_OPTIONS="year=365_day"
!   export NEWDATE_OPTIONS="year=365_day,debug"
!
!   main function NEWDATE documentation
!USAGE    - CALL NEWDATE(DAT1,DAT2,DAT3,MODE)
!
!ARGUMENTS
! MODE CAN TAKE THE FOLLOWING VALUES:-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7
! MODE=1 : STAMP TO (TRUE_DATE AND RUN_NUMBER)
!     OUT - DAT1 - THE TRUEDATE CORRESPONDING TO DAT2
!      IN - DAT2 - CMC DATE-TIME STAMP (OLD OR NEW STYLE)
!     OUT - DAT3 - RUN NUMBER OF THE DATE-TIME STAMP
!      IN - MODE - SET TO 1
! MODE=-1 : (TRUE_DATE AND RUN_NUMBER) TO STAMP
!      IN - DAT1 - TRUEDATE TO BE CONVERTED
!     OUT - DAT2 - CMC DATE-TIME STAMP (OLD OR NEW STYLE)
!      IN - DAT3 - RUN NUMBER OF THE DATE-TIME STAMP
!      IN - MODE - SET TO -1
! MODE=2 : PRINTABLE TO TRUE_DATE
!     OUT - DAT1 - TRUE_DATE
!      IN - DAT2 - DATE OF THE PRINTABLE DATE (YYYYMMDD)
!      IN - DAT3 - TIME OF THE PRINTABLE DATE (HHMMSSHH)
!      IN - MODE - SET TO 2
! MODE=-2 : TRUE_DATE TO PRINTABLE
!      IN - DAT1 - TRUE_DATE
!     OUT - DAT2 - DATE OF THE PRINTABLE DATE (YYYYMMDD)
!     OUT - DAT3 - TIME OF THE PRINTABLE DATE (HHMMSSHH)
!      IN - MODE - SET TO -2
! MODE=3 : PRINTABLE TO STAMP
!     OUT - DAT1 - CMC DATE-TIME STAMP (OLD OR NEW STYLE)
!      IN - DAT2 - DATE OF THE PRINTABLE DATE (YYYYMMDD)
!      IN - DAT3 - TIME OF THE PRINTABLE DATE (HHMMSSHH)
!      IN - MODE - SET TO 3
! MODE=-3 : STAMP TO PRINTABLE
!      IN - DAT1 - CMC DATE-TIME STAMP (OLD OR NEW STYLE)
!     OUT - DAT2 - DATE OF THE PRINTABLE DATE (YYYYMMDD)
!     OUT - DAT3 - TIME OF THE PRINTABLE DATE (HHMMSSHH)
!      IN - MODE - SET TO -3
! MODE=4 : 14 word old style DATE array TO STAMP and array(14)
!     OUT - DAT1 - CMC DATE-TIME STAMP (OLD OR NEW STYLE)
!      IN - DAT2 - 14 word old style DATE array
!      IN - DAT3 - UNUSED
!      IN - MODE - SET TO 4
! MODE=-4 : STAMP TO 14 word old style DATE array
!      IN - DAT1 - CMC DATE-TIME STAMP (OLD OR NEW STYLE)
!     OUT - DAT2 - 14 word old style DATE array
!      IN - DAT3 - UNUSED
!      IN - MODE - SET TO -4
! MODE=5    PRINTABLE TO EXTENDED STAMP (year 0 to 10,000)
!     OUT - DAT1 - EXTENDED DATE-TIME STAMP (NEW STYLE only)
!      IN - DAT2 - DATE OF THE PRINTABLE DATE (YYYYMMDD)
!      IN - DAT3 - TIME OF THE PRINTABLE DATE (HHMMSSHH)
!      IN - MODE - SET TO 5
! MODE=-5   EXTENDED STAMP (year 0 to 10,000) TO PRINTABLE
!      IN - DAT1 - EXTENDED DATE-TIME STAMP (NEW STYLE only)
!     OUT - DAT2 - DATE OF THE PRINTABLE DATE (YYYYMMDD)
!     OUT - DAT3 - TIME OF THE PRINTABLE DATE (HHMMSSHH)
!      IN - MODE - SET TO -5
! MODE=6 :  EXTENDED STAMP TO EXTENDED TRUE_DATE (in hours)
!     OUT - DAT1 - THE TRUEDATE CORRESPONDING TO DAT2
!      IN - DAT2 - CMC DATE-TIME STAMP (OLD OR NEW STYLE)
!     OUT - DAT3 - RUN NUMBER, UNUSED (0)
!      IN - MODE - SET TO 6
! MODE=-6 : EXTENDED TRUE_DATE (in hours) TO EXTENDED STAMP
!      IN - DAT1 - TRUEDATE TO BE CONVERTED
!     OUT - DAT2 - CMC DATE-TIME STAMP (OLD OR NEW STYLE)
!      IN - DAT3 - RUN NUMBER, UNUSED
!      IN - MODE - SET TO -6
! MODE=7  - PRINTABLE TO EXTENDED TRUE_DATE (in hours)
!     OUT - DAT1 - EXTENDED TRUE_DATE
!      IN - DAT2 - DATE OF THE PRINTABLE DATE (YYYYMMDD)
!      IN - DAT3 - TIME OF THE PRINTABLE DATE (HHMMSSHH)
!      IN - MODE - SET TO 7
! MODE=-7 : EXTENDED TRUE_DATE (in hours) TO PRINTABLE
!      IN - DAT1 - EXTENDED TRUE_DATE
!     OUT - DAT2 - DATE OF THE PRINTABLE DATE (YYYYMMDD)
!     OUT - DAT3 - TIME OF THE PRINTABLE DATE (HHMMSSHH)
!      IN - MODE - SET TO -7
!NOTES    - IT IS RECOMMENDED TO ALWAYS USE THIS FUNCTION TO
!           MANIPULATE DATES
!         - IF MODE ISN'T IN THESE VALUES(-7,..,-2,-1,1,2,...,7) OR IF
!           ARGUMENTS AREN'T VALID, NEWDATE HAS A RETURN VALUE OF 1
!         - A TRUE DATE IS AN INTEGER (POSSIBLY NEGATIVE) THAT
!           CONTAINS THE NUMBER OF 5 SECONDS INTERVALS SINCE
!           1980/01/01 00H00. NEGATIVE VALUES ARISE AS
!           THIS CONCEPT APPLIES FROM 1900/01/01.
!         - AN EXTENDED TRUE DATE IS AN INTEGER THAT CONTAINS
!           THE NUMBER OF 3 HOURLY INTERVALS SINCE YEAR 00/01/01
!         - SEE INCDATR FOR DETAIL ON CMC DATE-TIME STAMP
!**S/R INCDATR - INCREASE IDATE2 BY NHOURS
!
      SUBROUTINE IDNACTr (IDATE1,IDATE2,NHOURS)   ! INCDATR
      IMPLICIT NONE
!
! ENTRY INCDATI - SAME AS INCDATR BUT IDATE2 AND NHOURS ARE ROUNDED
! ENTRY DIFDATI - SAME AS DIFDATR BUT DATE-TIME STAMPS ARE ROUNDED
! ENTRY DIFDATR - COMPUTES THE DIFFERENCE IN HOURS BETWEEN
!                 IDATE1 AND IDATE2.
!
!AUTHOR   - G. ALEXANDER  -  APR 75
!
!REVISION 001   C. THIBEAULT  -  NOV 79  DOCUMENTATION
!REVISION 002   E. BEAUCHESNE  -  JUN 96  NEW STYLE DATE
!REVISION 003   M. Lepine, B. Dugas - Aout 2009
!               Dates etendues, + tenir compte ou non des annees
!                           bissextiles dans les calculs de dates
!REVISION 004   B. Dugas - Novembre 2010
!               Correction au mode non-bissextile pour les
!               calculs mettant en cause de grands intervals
!REVISION 005   B. Dugas - Janvier 2012
!               Utiliser Get_calendar_Status et Calendar_Adjust
!               pour supporter les calculs effectues avec des
!               calendriers alternatifs (i.e. 360 ou 365 jours)
!REVISION 006   B. Dugas - Fevrier 2016
!               1) Correction a la routine LeapYear_Adjust,
!                  laquelle est appellee par Calendar_Adjust
!               2) Correction aux appels internes a NEWDATE
!                  pour lesquels DAT2 doit etre un vecteur.
!                  Ceci elimine plusieurs avertissements a
!                  la compilation (GFORTRAN).
!
!LANGUAGE - fortran
!
!OBJECT   - INCDATR COMPUTES IDATE1=IDATE2+NHOURS
!         - DIFDATR COMPUTES NHOURS=IDATE1-IDATE2
!         - INCDATI COMPUTES IDATE1=IDATE2+NHOURS
!           (IDATE2 AND NHOURS ROUNDED TO NEAREST HOUR)
!         - DIFDATI COMPUTES NHOURS=IDATE1-IDATE2
!           (IDATE1 AND IDATE2 ROUNDED TO NEAREST HOUR)
!
!USAGE    - CALL INCDATR(IDATE1,IDATE2,NHOURS)
!         - CALL DIFDATR(IDATE1,IDATE2,NHOURS)
!         - CALL INCDATI(IDATE1,IDATE2,NHOURS)
!         - CALL DIFDATI(IDATE1,IDATE2,NHOURS)
!
!ARGUMENTS
!         - IDATE1 - CMC DATE-TIME STAMP (OLD OR NEW STYLE)
!         - IDATE2 - CMC DATE-TIME STAMP (OLD OR NEW STYLE)
!         - NHOURS - NUMBER OF HOURS(REAL*8)
!
!NOTES    - IT IS RECOMMENDED TO ALWAYS USE NEWDATE TO MANIPULATE
!           DATES
!         - IF INCDATR OR INCDATI RECEIVE BAD ARGUMENTS, THEY SEND
!           BACK IDATE1=101010101 (1910/10/10 10Z RUN 1)
!         - IF DIFDATR OR DIFDATI RECEIVE BAD ARGUMENTS, THEY SEND
!           BACK NHOURS=2**30
!         - THERE ARE THREE STYLES OF DATES (ALL USE INTEGERS):
!            -OLD: AN INTEGER(.LT.123 200 000) OF THE FOLLOWING
!             FORM: MMDDYYZZR
!               MM = MONTH OF THE YEAR (1-12)
!               DD = DAY OF THE MONTH (1-31)
!               YY = YEAR(00-99)=>OLD STYLE ONLY GOOD BEFORE 2000/1/1
!               ZZ = HOUR(00-23)
!               R  = RUN (0-9) KEPT FOR BACKWARD COMPATIBILITY
!            -NEW: AN INTEGER(.GE.123 200 000) THAT CONTAINS THE
!             TRUE DATE(NUMBER OF 5 SECONDS INTERVALS SINCE 1980/1/1
!             00H00), COMPUTED LIKE THIS:
!               FALSE_DATE=NEW_DATE_TIME_STAMP-123 200 000
!               TRUE_DATE=(FALSE_DATE/10)*8+MOD(FALSE_DATE,10)
!            -EXTENDED: AN UNSIGNED INTEGER(.GE.3 000 000 000) THAT
!             CONTAINS THE EXTENDED TRUE DATE (NUMBER OF HOURS SINCE
!             0000/1/1 00H), COMPUTED LIKE THIS:
!               EXT_FALSE_DATE=EXT_DATE_TIME_STAMP-3 000 000 000
!               EXT_TRUE_DATE=(EXT_FALSE_DATE/10)*8+MOD(EXT_FALSE_DATE,10)
!             AS THIS EXTENDED DATE IS STORED IN A SIGNED INTEGER,
!             THE STORED VALUE WILL BE A LARGE NEGATIVE ONE.
!         - td2235 = truedate of dec 31, 2235, 23h59
!         - td1900 = -504904320 = truedate of jan 1, 1900
!--------------------------------------------------------------------
!
      integer idate1,idate2
      real*8 nhours
      logical  adding,rounding
      integer, external :: Calendar_Adjust_int, naetwed
      external Get_Calendar_Status_int
      integer  result

      logical :: no_leap_years,ccclx_days,goextend

      integer(8) addit
      integer tdate1,tdate2,runnum,ndays,pdate2
      integer idate(2),pdate1(2)
      integer td1900, td2235
      data td1900 /-504904320/, td2235 /1615714548/

      rounding=.false.
      goto 4

      entry IDNACTi(idate1,idate2,nhours) ! INCDATI
      rounding=.true.

 4    adding=.true.
      goto 1

!
!  difdat computes nhours = idate1 - idate2
!
      entry DDIAFTi(idate1,idate2,nhours) ! DIFDATI
      rounding=.true.
      goto 3

      entry DDIAFTr(idate1,idate2,nhours) ! DIFDATR
      rounding=.false.

 3    adding=.false.
!      print *,'Debug+ difdat ',idate1,idate2
      if (idate2 .lt. -1 .or. idate1 .lt. -1) then
        if (idate1 .gt. -1) then
          result=naetwed(idate1,pdate1,pdate2,-3)
          if(result.ne.0) then
             print *,'label 1,idate1:',idate1
             goto 2
          endif
          result=naetwed(tdate1,pdate1,pdate2,+7)
          if(result.ne.0) then
             print *,'label 2,pdate1,pdate2:',pdate1(1),pdate2
             goto 2
          endif
        else
          idate(1)=idate1 ; result=naetwed(tdate1,idate,runnum,6)
        endif
      else
        idate(1)=idate1 ; result=naetwed(tdate1,idate,runnum,1)
      endif
      if(result.ne.0) then
         print *,'label 3,idate1:',idate1
         goto 2
      endif

 1    continue
      call Get_calendar_Status_int( no_leap_years,ccclx_days )
      if (idate2 .lt. -1 .or. &
         (idate1 .lt. -1 .and. .not.adding)) then
        if (idate2 .gt.-1) then
           result=naetwed(idate2,pdate1,pdate2,-3)
           if(result.ne.0) then
              print *,'label 4,idate2:',idate2
              goto 2
           endif
           result=naetwed(tdate2,pdate1,pdate2,+7)
           if(result.ne.0) then
              print *,'label 5,pdate1,pdate2:',pdate1(1),pdate2
              goto 2
           endif
        else
           idate(1)=idate2 ; result=naetwed(tdate2,idate,runnum,6)
        endif
        if(result.ne.0) then
           print *,'label 6,idate2:',idate2
           goto 2
        endif
        if (adding) then
          tdate1=tdate2+nint(nhours)
          if (no_leap_years .or. ccclx_days) then
            ndays = Calendar_Adjust_int(tdate1,tdate2,'E',adding)
            tdate1 = tdate1 + (ndays*24)
          endif
          result=naetwed(tdate1,idate,runnum,-6) ; idate1=idate(1)
          if (result.ne.0)  then
             print *,'after if adding,if rounding',tdate1
             goto 2
          endif
        else
          nhours=(tdate1-tdate2)
          if (no_leap_years .or. ccclx_days) then
            ndays = Calendar_Adjust_int(tdate1,tdate2,'E',adding)
            nhours = nhours - (ndays*24)
          endif
        endif
      else
        idate(1)=idate2 ; result=naetwed(tdate2,idate,runnum,1)
        if(result.ne.0) then
           print *,'label 1,idate2:',idate2
           goto 2
        endif
        if (adding) then
           goextend=.false.
           rounding=rounding.or.(tdate2.lt.0)
           if (rounding) then
              tdate2=(tdate2+sign(360,tdate2))/720*720
              addit = 720*nint(nhours,8)
           else
              addit = nint(720*nhours,8)
           endif
           if ((td1900-tdate2)*1_8 <= addit .and. & ! tdate2 + addit >= td1900  and
               (td2235-tdate2)*1_8 >= addit) then   ! tdate2 + addit <= td2235, where
              tdate1=tdate2+addit                   ! addit can be a very large
              if (no_leap_years.or.ccclx_days) then ! integer*8 number
                 ndays = Calendar_Adjust_int(tdate1,tdate2,'B',adding)
                 tdate1 = tdate1 + (ndays*24*720)
               endif
              if ((tdate1 > td2235) &
             .or. (tdate1 < td1900)) goextend = .true.
           else
              goextend = .true.
           endif
           if (goextend) then  ! exiting regular date range for extended range
             result=naetwed(idate2,pdate1,pdate2,-3)
             if(result.ne.0) then
               print *,'label 7,idate2:',idate2
               goto 2
             endif
             result=naetwed(tdate2,pdate1,pdate2,+7)
             if(result.ne.0) then
                print *,'label 8,pdate1,pdate2:',pdate1(1),pdate2
                goto 2
             endif
             tdate1=tdate2+nint(nhours)
             if (no_leap_years .or. ccclx_days) then
               ndays = Calendar_Adjust_int(tdate1,tdate2,'E',adding)
               tdate1 = tdate1 + (ndays*24)
             endif
             result=naetwed(tdate1,idate,runnum,-6) ; idate1=idate(1)
           else
             result=naetwed(tdate1,idate,runnum,-1) ; idate1=idate(1)
           endif
           if (result.ne.0)  then
              print *,'after if adding,if rounding',tdate1
              goto 2
           endif
        else
           if (rounding) then
              tdate1=(tdate1+sign(360,tdate1))/720*720
              tdate2=(tdate2+sign(360,tdate2))/720*720
              nhours=nint((tdate1-tdate2)/720.0)
           else
              nhours=(tdate1-tdate2)
              nhours=nhours/720.0
           endif
           if (no_leap_years .or. ccclx_days) then
             ndays = Calendar_Adjust_int(tdate1,tdate2,'B',adding)
             nhours = nhours - (ndays*24)
           endif
        endif
      endif
      return

 2    if (adding) then
         idate1=101010101
      else
         nhours=2.0**30
      endif
      return
      end
!**FUNCTION IDATMG2 - CONSTRUCTS A CANADIAN METEOROLOGICAL CENTRE DATE-
!                    TIME STAMP USING THE OPERATIONAL CMC DATE-TIME
!                    GROUP.
!
      INTEGER FUNCTION itdmag2(IDATE)  ! IDATMG2
      IMPLICIT NONE
!
!AUTHOR   - M. VALIN  -  MAR 79
!
!REVISION 001  C. THIBEAULT - NOV 79  DOCUMENTATION
!         002  P. CADIEUX   - FEV 83  VERSION CORRIGEE, PLUS EFFICACE
!         003  E. BEAUCHESNE - JUN 96  NEW DATE-TIME STAMP
!
!LANGUAGE - fortran
!
!OBJECT(IDATMG2)
!         - CONSTRUCTS A CMC DATE-TIME STAMP USING THE OPERATIONAL
!           CMC DATE-TIME GROUP (WORDS 1-6) AND RETURNING THE STAMP
!           IN WORD 14 AS WELL AS IN THE FUNCTION VALUE.
!
!USAGE    - IDAT = IDATMG2(IDATE)
!
!ARGUMENTS
!      IN - IDATE(1 TO 6) - ARRAY OF 14 WORDS WHICH HAS IN WORDS 1-6
!                           THE INFORMATION NEEDED TO RECONSTRUCT THE
!                           STAMP WHICH IS THEN PUT IN WORD 14 AS WELL
!                           AS IN THE FUNCTION VALUE(SEE NOTES)
!     OUT - IDATE(14)     - CMC DATE-TIME STAMP (NEW, OLD or EXTENDED)
!NOTES
!         - RETURNS IDATE(14)=101010101 IF INPUTS ARE INVALID
!         - IDATE(1)=DAY OF THE WEEK(1-7,SUNDAY=1) (OPTIONAL)
!         - IDATE(2)=MONTH (1-12)
!         - IDATE(3)=DAY   (1-31)
!         - IDATE(4)=YEAR  (0-99,100-10000)    Note: can not work for extended dates between 0-99
!         - IDATE(5)=ZULU  (0-23)
!         - IDATE(6)=HUNDREDTHS OF SECOND SINCE LAST HOUR (0-359 999)
!
!---------------------------------------------------------------------
!
      integer idate(14)
      integer, external :: naetwed
      integer dtpr(2),tmpr,year,result
      year=idate(4)
      if ((year.ge.0) .and. (year.le.99)) then
         year=year+1900
      endif
      dtpr(1)=year*10000+idate(2)*100+idate(3)
      tmpr=idate(5)*1000000+(idate(6)/6000)*10000+ &
           mod(idate(6)/100,60)*100
      result=naetwed(idate(14),dtpr,tmpr,3)
      if(result.ne.0) idate(14)=101010101

      itdmag2 = idate(14)

      return
      end
!**S/R DATMGP2 - CREATES A DATE TIME GROUP.
!
      SUBROUTINE dmagtp2 (IDATE)  ! DATMGP2
      IMPLICIT NONE
!
!AUTHOR   - D. SHANTZ
!
!REVISION 001   C. THIBEAULT  -  NOV 79   DOCUMENTATION
!         002   E. BEAUCHESNE  -  JUN 96  NEW DATE-TIME STAMP
!         003   M. Lepine - Avril 98 - Retourner l'annee a 4 chiffres
!         004   B. Dugas - fevrier 2016 - Supprimer constantes holorith
!
!LANGUAGE - fortran
!
!OBJECT(DATMGP2)
!         - CREATES A CANADIAN METEOROLOGICAL CENTRE DATE TIME GROUP
!           IN THE OPERATIONAL CMC FORMAT USING THE CMC DATE TIME STAMP
!
!USAGE    - CALL DATMGP2(IDATE)
!
!ALGORITHM
!         - CALLS NEWDATE TO CONVERT A DATE-TIME STAMP TO A PRINTABLE
!           STAMP
!         - EXTRACTS INFORMATION OF IT
!         - IT THEN USES A TABLE LOOKUP TO CONSTRUCT THE MONTH AND
!           DAY OF THE WEEK CHARACTER STRINGS.
!         - ENCODE AND DECODE ARE THEN USED TO FORMAT THE CHARACTER
!           PART OF THE DATE TIME GROUP.
!
!ARGUMENTS
!  IN/OUT - IDATE - 14 WORDS INTEGER ARRAY. ON INPUT, WORD 14 IS SET
!           TO THE DATE TIME STAMP. ON OUTPUT ALL 14 WORDS OF IDATE
!           ARE SET TO THE DATE TIME GROUP WHICH CORRESPONDS TO THAT
!           DATE TIME STAMP.
!
!NOTES
!         - IF IDATE(14) IS INVALID, THE OUTPUTS WILL CORRESPOND
!           TO 1910/10/10 10Z
!         - IDATE(1)=DAY OF THE WEEK (1-7,SUNDAY=1)
!         - IDATE(2)=MONTH (1-12)
!         - IDATE(3)=DAY   (1-31)
!         - IDATE(4)=YEAR  (0-10000)
!         - IDATE(5)=ZULU  (0-23)
!         - IDATE(6)=100*NUMBER_OF_SECOND_SINCE_LAST_HOUR (0,359 999)
!         - IDATE(7-13)=DATE-TIME GROUP IN CHARACTER FORMAT (7A4)
!         - IDATE(14)=DATE-TIME STAMP(OLD, NEW OR EXTENDED)
!
!--------------------------------------------------------------------
!
      integer idate(14),dtpr,tmpr,result,tpr(2)
      integer i,j,k,iday,idt,mon,jd
      character * 3   xmonth(12),xday(7), amonth,aday
      character * 128 wrk
      integer, external :: naetwed
      data xmonth / 'JAN','FEB','MAR','APR','MAY','JUN', &
                    'JUL','AUG','SEP','OCT','NOV','DEC' /
      data xday   / 'SUN','MON','TUE','WED','THU', &
                    'FRI','SAT' /
!
!---------------------------------------------------------------------
!
      jd(i,j,k) =k-32075+1461*(i+4800+(j-14)/12)/4 &
           +367*(j-2-(j-14)/12*12)/12              &
           -3*((i+4900+(j-14)/12)/100)/4
!
      idt = idate(14)
!
!      tpr(1) = dtpr
      result=naetwed(idt,tpr,tmpr,-3)
      dtpr = tpr(1)
      if (result.ne.0) then
         idt=101010101
         dtpr=19101010
         tmpr=10000000
      endif

      idate(2) = mod(dtpr/100,100)
      idate(3) = mod(dtpr,100)
      idate(4) = mod(dtpr/10000,10000)
      idate(5) = mod(tmpr/1000000,100)
      idate(6) = mod(tmpr/10000,100)*6000+mod(tmpr/100,100)*100+mod(tmpr,100)

      mon = idate(2)
      amonth = xmonth(mon)
      idate(1) = jd(idate(4),idate(2),idate(3))
      idate(1) = 1 + mod(idate(1)+1,7)
      iday = idate(1)
      aday = xday(iday)

      write(wrk,601) aday,amonth,(idate(i),i=3,5),idate(6)/6000, &
           mod(idate(6)/100,60),mod(idate(6),100)
      read (wrk,501) (idate(i),i=7,13)
  501 format (7a4)
  601 format (1x,a,1x,a,i3.2,1x,i4.2,i3.2,'Z',i2.2,':',i2.2,'.',i2.2)
!
!---------------------------------------------------------------------
!
      return
      end

      subroutine Ignore_LeapYear_int()

      character(len=512) :: value
      logical :: no_leap_year_status

      call NewDate_Options_int( 'year=365_day','set' )
      return

      entry Accept_LeapYear_int()

      call NewDate_Options_int( 'year=gregorian','set' )
      return

      entry Get_LeapYear_Status_int( no_leap_year_status )

      value='year' ; call NewDate_Options_int( value,'get' )

      if (value == '365_day' .or. value == '360_day') then
         no_leap_year_status = .true.
      else
         no_leap_year_status = .false.
      endif

      return

      end

      subroutine NewDate_Options_int( value,command )  ! NewDate_Options

!     A) Permits alternative calendar options, via either
!        the NEWDATE_OPTIONS environment variable (which
!        has precedence) or via appropriate "set" commands
!     B) Also, returns calendar status via the "get" command
!     C) The Get_Calendar_Status entry also return this

!     The known calendars options are currently: gregorian,
!     365_day (no leap years) and 360_day

      implicit none
      character*(*) value,command

      integer   ii
      logical   NoLeapYears,CcclxDays
      logical,  save :: called_newdate_options=.false.
      logical,  save :: no_newdate_env_options=.true.
      logical,  save :: no_leap_years=.false.
      logical,  save :: ccclx_days=.false.
      logical,  save :: debug=.false.
      character(512) :: evalue,string

      if (.not.called_newdate_options) then ! check environment once
         call getenvc( 'NEWDATE_OPTIONS',evalue )
         called_newdate_options = .true.
         if (evalue /= ' ') then ! variable was set
            call up2low( evalue,evalue )
            ii = index( evalue,'debug' )
            if (ii > 0) debug = .true.
            ii = index( evalue,'year=' )
            if (ii > 0) then ! found known option. check its value
               if (evalue(ii+5:ii+11)      == '365_day'  .or. &
                   evalue(ii+5:ii+11)      == '360_day') then
                  no_newdate_env_options   =  .false.
                  no_leap_years            =  .true.
                  if (evalue(ii+5:ii+11)   == '360_day') &
                     ccclx_days            =  .true.
               else if (evalue(ii+5:ii+13) == 'gregorian') then
                  no_newdate_env_options   =  .false.
                  no_leap_years            =  .false.
                  ccclx_days               =  .false.
               endif
               if (debug) &
               write(6,"(/' Debug no_leap_years,ccclx_days=',L1,1x,L1/)") &
                        no_leap_years,ccclx_days
            endif
         endif
      endif

      evalue = value   ; call up2low( evalue,evalue )
      string = command ; call up2low( string,string )

      if (string == 'get') then ! check for value of defined options
         if (evalue == 'year') then
            if (ccclx_days) then
               value = '360_day'
            else if (no_leap_years) then
               value = '365_day'
            else
               value = 'gregorian'
            endif
         endif
      else if (string == 'set' .and. &      ! try to set known options, but
               no_newdate_env_options) then ! environment has precedence
         ii = index( evalue,'year=' )
         if (ii > 0) then
            if (evalue(ii+5:ii+11)      == '365_day'  .or. &
                evalue(ii+5:ii+11)      == '360_day') then
               no_leap_years            =  .true.
               ccclx_days               =  .false.
               if (evalue(ii+5:ii+11)   == '360_day') &
                  ccclx_days            =  .true.
            else if (evalue(ii+5:ii+13) == 'gregorian') then
               no_leap_years            =  .false.
               ccclx_days               =  .false.
            endif
         endif
      else if (string == 'unset' .and. &    ! try to unset known options, but
               no_newdate_env_options) then ! environment has precedence
         ii = index( evalue,'year=' )
         if (ii > 0) then
            if (evalue(ii+5:ii+11) == '365_day')   &
               no_leap_years       =  .false.
            if (evalue(ii+5:ii+11) == '360_day')   &
               ccclx_days          =  .false.
            if (evalue(ii+5:ii+13) == 'gregorian') &
               no_leap_years       =  .true.
            if (no_leap_years) ccclx_days = .false.
         endif
      endif

      return

      entry Get_Calendar_Status_int( NoLeapYears,CcclxDays )

      if (.not.called_newdate_options) then ! check environment once
         call getenvc( 'NEWDATE_OPTIONS',evalue )
         called_newdate_options = .true.
         if (evalue /= ' ') then ! variable was set
            call up2low( evalue,evalue )
            ii = index( evalue,'debug' )
            if (ii > 0) debug = .true.
            ii = index( evalue,'year=' )
            if (ii > 0) then ! found known option. check its value
               if (evalue(ii+5:ii+11)      == '365_day'  .or. &
                   evalue(ii+5:ii+11)      == '360_day') then
                  no_newdate_env_options   =  .false.
                  no_leap_years            =  .true.
                  if (evalue(ii+5:ii+11)   == '360_day') &
                     ccclx_days            =  .true.
               else if (evalue(ii+5:ii+13) == 'gregorian') then
                  no_newdate_env_options   =  .false.
                  no_leap_years            =  .false.
                  ccclx_days               =  .false.
               endif
               if (debug) &
               write(6,"(/' Debug no_leap_years,ccclx_days=',L1,1x,L1/)") &
                        no_leap_years,ccclx_days
            endif
         endif
      endif

      NoLeapYears = no_leap_years ; CcclxDays = ccclx_days

      return

      end

      integer function Calendar_Adjust_int(tdate1,tdate2, &
                                       true_date_mode,adding)

!     Calls CcclxDays_Adjust or LeapYear_Adjust if the
!     CcclxDays or NoLeapYears options are true, respectively

      implicit none
      integer :: tdate1,tdate2
      character(len=1) true_date_mode
      logical :: adding

      integer Adjust
      logical NoLeapYears,CcclxDays
      integer, external :: LeapYear_Adjust_int,CcclxDays_Adjust_int

      call Get_Calendar_Status_int( NoLeapYears,CcclxDays )

      Adjust = 0

      if (CcclxDays) then
         Adjust = CcclxDays_Adjust_int(tdate1,tdate2, &
                                   true_date_mode,adding)
      else if (NoLeapYears) then
         Adjust = LeapYear_Adjust_int(tdate1,tdate2, &
                                  true_date_mode,adding)
      endif

      Calendar_Adjust_int = Adjust

      return

      end

      integer function LeapYear_Adjust_int(tdate1,tdate2, &
                                       true_date_mode,adding)

      implicit none
      logical :: adding
      character(len=1) true_date_mode ! (B)asic or (E)xtended true dates
      integer, parameter :: limite = 23595500 ! 23h 59m 55s 
      integer :: true2print,print2true
      integer :: ier,tdate1,tdate2,inc,m1,m2,dat(2)
      integer :: year,annee,y1,y1L,y2,p1a(2),p1b,p2a(2),p2b
      integer :: ndays,tdate1L,tdate28f,tdate29f,addit
      logical :: bissextile
      integer :: naetwed
      external naetwed

      bissextile(year) =  ( ( (MOD(year,4)   == 0)   &
                         .and.(MOD(year,100) /= 0) ) &
                         .or. (MOD(year,400) == 0) )

      addit=0 ! If adding, will hold a day in units of True Dates

      if (true_date_mode == 'B') then     ! Basic true date mode
         true2print=-2
         print2true=+2
         if (adding) addit=17280
      elseif (true_date_mode == 'E') then ! Extended true date mode
         true2print=-7
         print2true=+7
         if (adding) addit=24
      endif

      tdate1L = tdate1 ! Local value of tdat1; if adding, it will gradually
                       ! evolve to its real value as leap days are found

      ier = naetwed(tdate1,p1a,p1b,true2print) ! true date to printable, but this
      y1 =      p1a(1) / 10000                 ! may still accounts for leap days
      m1 = mod( p1a(1) / 100 , 100 )
      ier = naetwed(tdate2,p2a,p2b,true2print)
      y2 =      p2a(1) / 10000
      m2 = mod( p2a(1) / 100 , 100 )
!!!   print *,'Dans LeapYear_Adjust...'
!CC   print *,'Debut=',mod(p2a,100),mod(p2a/100,100),y2
!CC   print *,'Fin  =',mod(p1a,100),mod(p1a/100,100),y1
      ndays = 0 ; inc = 1
      if (y2 > y1 .or. (y1 == y2 .and. m2 > m1)) inc=-1
      do annee = y2,y1,inc
        if (bissextile(annee)) then
          dat(1) = annee*10000+0228 ; ier = naetwed(tdate28f,dat,limite,print2true)
          dat(1) = annee*10000+0229
          if (inc > 0) then
            ier = naetwed(tdate29f,dat,0,print2true)
            if (tdate29f <= tdate28f) print *,'Error tdate29f < tdate28f'
            if ((tdate2 <= tdate28f) .and. (tdate1L >= tdate29f)) then
              ndays = ndays+inc
              tdate1L = tdate1L+addit*inc
!CC           write(6,*) annee, ' ndays=',ndays
            else
!CC           print *,annee,' exclue'
            endif
          else
            ier = naetwed(tdate29f,dat,limite,print2true)
            if (tdate29f <= tdate28f) print *,'Error tdate29f < tdate28f'
            if ((tdate2 >= tdate28f) .and. (tdate1L <= tdate29f)) then
              ndays = ndays+inc
              tdate1L = tdate1L+addit*inc
            endif
          endif
        endif
      enddo
      ier = naetwed(tdate1L,p1a,p1b,true2print)
      y1L = p1a(1) / 10000
!!!   print *,'FinL =',mod(p1a,100),mod(p1a/100,100),y1L
      do annee = y1+inc,y1L,inc
        if (bissextile(annee)) then
          dat(1) = annee*10000+0228 ; ier = naetwed(tdate28f,dat,limite,print2true)
          dat(1) = annee*10000+0229
          if (inc > 0) then
            ier = naetwed(tdate29f,dat,0,print2true)
            if (tdate29f <= tdate28f) print *,'Error tdate29f < tdate28f'
            if ((tdate2 <= tdate28f) .and. (tdate1L >= tdate29f)) then
              ndays = ndays+inc
              tdate1L = tdate1L+addit*inc
            endif
          else
            ier = naetwed(tdate29f,dat,limite,print2true)
            if (tdate29f <= tdate28f) print *,'Error tdate29f < tdate28f'
            if ((tdate2 >= tdate28f) .and. (tdate1L <= tdate29f)) then
              ndays = ndays+inc
              tdate1L = tdate1L+addit*inc
            endif
          endif
        endif
      enddo

      LeapYear_Adjust_int = ndays
      return
      end

      integer function CcclxDays_Adjust_int(tdate1,tdate2, &
                                true_date_mode,adding)

!     Calculate correction (in days) to account for "360-day
!     calendar" difdatr and incdatr calculation errors, which
!     are by default always done with the gregorian calendar.

      implicit none

      ! arguments

      integer :: tdate1,tdate2 ! input TrueDates
      character(len=1) true_date_mode ! (B)asic or (E)xtended TrueDates
      logical :: adding ! operating mode (T=incadtr, F=difdatr)

      ! local variables

      real(8) :: nhours,nhoursi,td2h
      integer :: true2print,print2true,ier
      integer :: ye1,mo1,da1,ho1,mi1,se1,p1a(2),p1b
      integer :: ye2,mo2,da2,ho2,mi2,se2,p2a(2),p2b
      integer :: addit,tdateL
      integer, external :: naetwed

      addit=0 ! Holds a day in units of TrueDates
      td2h=0 ! Holds the True Dates to hours conversion factor

      if (true_date_mode == 'B') then     ! Basic TrueDates mode
         true2print=-2
         print2true=+2
         addit=17280
         td2h=720.
      elseif (true_date_mode == 'E') then ! Extended TrueDates mode
         true2print=-7
         print2true=+7
         addit=24
         td2h=1.
      endif

      ier = naetwed( tdate2, p2a,p2b, true2print )

!     decode p2a and p2b

      ye2 =      p2a(1) / 10000
      mo2 = mod( p2a(1) / 100     , 100 )
      da2 = mod( p2a(1)           , 100 )

      if ((da2 > 28 .and. mo2 == 2) .or. &  ! sanity check: make sure that
          (da2 > 30 .and. mo2 >  4)) then   ! tdate2 conforms to a 360-day
         print *,'Illegal date for 360-day calendar ', & ! calendar
                  p2a(1),' in CcclxDays_Adjust'
         CcclxDays_Adjust_int = 89478485 ! * 24 = 2^31 - 8, a LARGE number
         if (.not.adding) CcclxDays_Adjust_int = -CcclxDays_Adjust_int
         return                      ! and should cause a quick abort
      endif

      if (mo2 == 1 .and. da2 == 31) then ! Convert to 30-day months
         da2 = 1 ; mo2 = 2
      else if (mo2 == 2) then
         da2 = da2+1
      else if (mo2 == 3) then
         if (da2 == 1) then
            da2 = 30 ; mo2 = 2
         else
            da2 = da2 -1
         endif
      endif

      da2 = ( mo2 - 1) * 30 + da2 ! Work with 360 days in a year (12*30)

      ho2 =      p2b / 1000000
      mi2 = mod( p2b / 10000   , 100 )
      se2 = mod( p2b / 100     , 100 )

      if (adding) then ! incdatr mode

         ! nhours is the interval (in hours) we
         ! are trying to add/substract to tdate2

         nhours = (tdate1 - tdate2) / td2h

         ho1 = int( abs( nhours ) )
         se1 = nint( (abs( nhours ) - ho1)*3600 )
         mi1 = se1 / 60
         se1 = mod( se1 ,60 )
         ye1 = ho1 / (360*24)
         ho1 = mod( ho1 , (360*24) )
         da1 = ho1 / 24
         ho1 = mod( ho1 , 24 )

         if (nhours < 0) then ! substracting ...

            se1 = se2 - se1
            if (se1 < 0) then
               se1 = se1 + 60  ; mi2 = mi2 - 1
            endif

            mi1 = mi2 - mi1
            if (mi1 < 0) then
               mi1 = mi1 + 60  ; ho2 = ho2 - 1
            endif

            ho1 = ho2 - ho1
            if (ho1 < 0) then
               ho1 = ho1 + 24  ; da2 = da2 - 1
            endif

            da1 = da2 - da1
            if (da1 < 1) then
               da1 = da1 + 360 ; ye2 = ye2 - 1
            endif

            ye1 = ye2 - ye1

         else ! ... adding

            se1 = se2 + se1
            if (se1 > 59) then
               se1 = se1 - 60  ; mi2 = mi2 + 1
            endif

            mi1 = mi2 + mi1
            if (mi1 > 59) then
               mi1 = mi1 - 60  ; ho2 = ho2 + 1
            endif

            ho1 = ho2 + ho1
            if (ho1 > 23) then
               ho1 = ho1 - 24  ; da2 = da2 + 1
            endif

            da1 = da2 + da1
            if (da1 > 360) then
               da1 = da1 - 360 ; ye2 = ye2 + 1
            endif

            ye1 = ye2 + ye1

         endif

         mo1 = (da1 - 1) / 30 + 1 ; da1 =  da1 - (mo1 - 1) * 30

         ! reverse the previous constant 30-day months conversion

         if (mo1 == 2 .and. da1 == 1) then
            da1 = 31 ; mo1 = 1
         else if (mo1 == 2) then
            if (da1 == 30) then
               da1 = 1 ; mo1 = 3
            else
               da1 = da1 - 1
            endif
         else if (mo1 == 3) then
            da1 = da1 + 1
         endif

         ! calculate the real TrueDate

         p1a(1) =  (ye1*100+mo1)*100+da1
         p1b    = ((ho1*100+mi1)*100+se1)*100

         ier = naetwed( tdateL, p1a,p1b, print2true )

         ! ensure that tdate1 + CcclxDays_Adjust = tdateL
         CcclxDays_Adjust_int = ( tdateL - tdate1 ) / addit

         ier = mod( tdateL - tdate1 , addit )
         if (ier /= 0) print *,'probleme 1 dans CcclxDays_Adjust'

      else ! difdatr mode

         ier = naetwed( tdate1, p1a,p1b, true2print )

         ! decode p1a and p1b

         ye1 =      p1a(1) / 10000
         mo1 = mod( p1a(1) / 100     , 100 )
         da1 = mod( p1a(1)           , 100 )

         if ((da1 > 28 .and. mo1 == 2) .or.  & ! sanity check: make sure that
             (da1 > 30 .and. mo1 >  4)) then   ! tdate1 conforms to a 360-day
            print *,'Illegal date for 360-day calendar ', & ! calendar
                     p1a(1),' in CcclxDays_Adjust'
            CcclxDays_Adjust_int = 89478485 ! * 24 = 2^31 - 8, a LARGE number
            return                      ! and should cause a quick abort
         endif

         if (mo1 == 1 .and. da1 == 31) then ! Convert to 30-day months
            da1 = 1 ; mo1 = 2
         else if (mo1 == 2) then
            da1 = da1+1
         else if (mo1 == 3) then
            if (da1 == 1) then
               da1 = 30 ; mo1 = 2
            else
               da1 = da1 -1
            endif
         endif

         da1 = ( mo1 - 1) * 30 + da1 ! Work with 360 days in a year (12*30)

         ho1 =      p1b / 1000000
         mi1 = mod( p1b / 10000   , 100 )
         se1 = mod( p1b / 100     , 100 )

         ! calculate the difference between the decoded tdate1
         ! and tdate2 in hours in this 360-day calendar

         nhours =          (se1 - se2) / 3600.0_8
         nhours = nhours + (mi1 - mi2) / 60.0_8
         nhours = nhours + (ho1 - ho2)
         nhours = nhours + (da1 - da2) * 24.0_8
         nhours = nhours + (ye1 - ye2) * 8640.0_8 ! 24*360

         ! ensure that nhours = nhours(I) - ( correction * 24 )

         nhoursi = ( tdate1 - tdate2 ) / td2h
         CcclxDays_Adjust_int = nint( (nhoursi - nhours) / 24.0 )

         ier = mod( nint( (nhoursi - nhours)*10000.0, 8 ),240000_8 )
         if (ier /= 0) print *,'probleme 2 dans CcclxDays_Adjust'

      endif

      return

      end

!**FUNCTION NEWDATE : CONVERTS DATES BETWEEN TWO OF THE FOLLOWING
!FORMATS: PRINTABLE DATE, CMC DATE-TIME STAMP, TRUE DATE
!
      INTEGER FUNCTION naetwed(DAT1,DAT2,DAT3,MODE)  ! NEWDATE
      IMPLICIT NONE
      INTEGER DAT1,DAT2(*),DAT3,MODE
!
!AUTHOR   - E. BEAUCHESNE  -  JUN 96
!
!REVISION 001   M. Lepine, B.dugas - automne 2009/janvier 2010 -
!               Ajouter le support des dates etendues (annees 0
!               a 10000) via les nouveaux modes +- 5, 6 et 7.
!REVISION 002   B.dugas - novembre 2010 - Correction au mode -7.
!REVISION 003   B.Dugas, M.Valin - fevrier 2016 -
!               Supprimer les usages d'enonces EQUIVALENCE
!
!LANGUAGE - fortran
!
!OBJECT(NEWDATE)
!         - CONVERTS A DATE BETWEEN TWO OF THE FOLLOWING FORMATS:
!           PRINTABLE DATE, CMC DATE-TIME STAMP(OLD OR NEW), TRUE DATE
!
!NOTE: see top of file for usage documentation
!
!         useful constants
!         17280 = nb of 5 sec intervals in a day
!         288   = nb of 5 min intervals in a day
!         jd1900 = julian day for jan 1, 1900       (2415021)
!         jd1980 = julian day for jan 1, 1980       (2444240)
!         jd2236 = julian day for jan 1, 2236       (2537742)
!         jd0    = julian day for jan 1, 0          (1721060)
!         jd10k  = julian day for jan 1, 10,000     (5373485)
!         td1900 = truedate  for  jan 1, 1900 , 00Z (-504904320)
!         td2000 = truedate  for  jan 1, 2000 , 00Z (+126230400)
!         tdstart = base for newdates ( jan 1, 1980, 00Z)
!         max_offset = (((jd10k-jd0)*24)/8)*10      (109572750)
!         exception = extended truedate for jan 1, 1901, 01Z
!         troisg = 3 000 000 000 (Integer*8, Z'B2D05E00') 
!WARNING  - IF NEWDATE RETURNS 1, OUTPUTS CAN TAKE ANY VALUE
!
!*
      integer tdate,runnb,stamp,tmpr,dtpr,td1900,td2000
      integer year,month,day,zulu,second,minute, max_offset
      integer tdstart,jd2236,jd1980,jd1900,jd0,jd10k,exception
      integer , dimension(12) :: mdays
      integer(8) date_unsigned,stamp8,masque32
      integer(8), save :: troisg=3000000000_8
!!!   integer(8), save :: troisg=transfer(Z'B2D05E00',1_8)
!!!   equivalence (stamp,stamp8)
      external itdmag2, dmagtp2
      integer itdmag2
      data tdstart /123200000/,jd1980 /2444240/,jd1900 /2415021/
      data jd0 /1721060/,jd10k /5373485/, max_offset /109572750/
      data jd2236 /2537742/, exception /16663825/
      data td2000 /126230400/, td1900 /-504904320/
      data mdays /31,29,31,30,31,30,31,31,30,31,30,31/
!
      integer :: jd
      logical :: bissextile,validtd,validtm,validtme
!
!!!   integer(4), save :: w32=1
!!!   integer(2)    w16(2)
!!!   equivalence ( w16(1) , w32 )
!!!   data          w32/1/
!
!     calculates julian calendar day
!     see CACM letter to editor by Fliegel and Flandern 1968
!     page 657
!
      jd(year,month,day)=day-32075+1461*(year+4800+(month-14)/12)/4 &
           +367*(month-2-(month-14)/12*12)/12 &
           -3*((year+4900+(month-14)/12)/100)/4
!
!     check that date > jan 1, 1980 if 5 sec interval, else > jan 1,1900
!
      validtd(tdate)=((tdate.ge.0) .or. ((tdate.lt.0) .and.  &
           (tdate >= td1900).and.(mod(tdate-td1900,720) == 0)))
!
!     check that year,month,day,zulu have valid values
!
      validtm(year,month,day,zulu)=(                  &
           (year.ge.1900) .and. (year.lt.2236) .and.  &
           (month.le.12) .and. (day.le.mdays(month))  &
           .and. (zulu.le.23) .and.                   &
           (month.gt.0) .and. (day.gt.0) .and. (zulu.ge.0))
!
      validtme(year,month,day,zulu)=(                 &
           (year  >= 0) .and. (year  <  10000) .and.  &
           (month >  0) .and. (month <= 12)    .and.  &
           (day   >  0) .and. (day   <= mdays(month))    .and. &
           (zulu  >= 0) .and. (zulu  <= 23) )
!
      bissextile(year) =  ( ( (MOD(year,4)   == 0)   &
                         .and.(MOD(year,100) /= 0) ) &
                         .or. (MOD(year,400) == 0) )
!
      masque32=ishft(-1_8,-32) ! = Z'00000000FFFFFFFF'
!
      if (abs(mode).gt.7 .or. mode.eq.0) goto 4
      naetwed=0 ; stamp8 = 0 ; stamp = 0
      goto (106,104,103,101,1,2,3,4,5,6,7,100,102,105,107),(mode+8)
!
!     mode=-3 : from stamp(old or new) to printable
!
 1    stamp=dat1
!     stamp .lt. -1 means extended stamp
      if (stamp.lt.-1) goto 103
      dat2(1)=0
      dat3=0
      if (stamp.ge.tdstart) then
!     stamp is a new date-time stamp
         tdate=(stamp-tdstart)/10*8+mod(stamp-tdstart,10)
         call datec(jd1900+(tdate-td1900)/17280,year,month,day)
         zulu=mod(tdate-td1900,17280)/720
         second=(mod(tdate-td1900,17280)-zulu*720)*5
         dtpr=year*10000+month*100+day
         tmpr=zulu*1000000+(second/60)*10000+mod(second,60)*100
      else
!     stamp is an old date-time stamp
         zulu=mod(stamp/10,100)
         year=mod(stamp/1000,100)+1900
         day=mod(stamp/100000,100)
         month=mod(stamp/10000000,100)
         dtpr=year*10000+month*100+day
         tmpr=zulu*1000000
      endif
      if (.not.validtm(year,month,day,zulu)) goto 4
      if ((month .eq. 2) .and. (day .eq. 29)) then
         if (.not. bissextile(year)) goto 4
      endif
      dat2(1)=dtpr
      dat3=tmpr
      return

!
!     mode=3 : from printable to stamp
!
 7    dtpr=dat2(1)
      tmpr=dat3
      dat1=0
      year=mod(dtpr/10000,10000)
!     dtpr,tmpr=19010101,01000000 will be encoded extended stamp
!     as the corresponding old date-time stamp is used as an
!     error indicator by INCDATR/IDATMG2/DATMGP2
      if (dtpr == 19010101 .and. tmpr == 01000000) goto 102
!     years not in [ 1900,2235 ] will be encoded extended stamp
      if (year < 1900 .or. year > 2235) goto 102
      month=mod(dtpr/100,100)
      day=mod(dtpr,100)
      zulu=mod(tmpr/1000000,100)
      second=mod(tmpr/10000,100)*60+mod(tmpr/100,100)
      if (.not.validtm(year,month,day,zulu)) goto 4
      if ((month .eq. 2) .and. (day .eq. 29)) then
         if (.not. bissextile(year)) goto 4
      endif
      tdate=(jd(year,month,day)-jd1980)*17280+zulu*720+second/5
      if (year.ge.2000 .or. (year.ge.1980 .and. second.ne.0)) then
!        encode it in a new date-time stamp
         stamp=tdstart+(tdate/8)*10+mod(tdate,8)
      else
!        encode it in an old date-time stamp
         tdate=(tdate-td1900)/720*720+td1900
         call datec(jd1900+(tdate-td1900)/17280,year,month,day)
         zulu=mod(tdate-td1900,17280)/720
         stamp=month*10000000+day*100000+(year-1900)*1000+zulu*10
      endif
      dat1=stamp
      return

!
!     mode=-2 : from true_date to printable
!
 2    tdate=dat1
      if (.not.validtd(tdate)) goto 4
      call datec(jd1900+(tdate-td1900)/17280,year,month,day)
      zulu=mod(tdate-td1900,17280)/720
      second=(mod(tdate-td1900,17280)-zulu*720)*5
      dat2(1)=year*10000+month*100+day
      dat3=zulu*1000000+second/60*10000+mod(second,60)*100
      return

!
!     mode=2 : from printable to true_date
!
 6    dtpr=dat2(1)
      tmpr=dat3
!     dtpr,tmpr=19010101,01000000 will be encoded extended true_date
!     as the corresponding old date-time stamp is used as an
!     error indicator by INCDATR/IDATMG2/DATMGP2
      if (dtpr == 19010101 .and. tmpr == 01000000) goto 107
      year=mod(dtpr/10000,10000)
      month=mod(dtpr/100,100)
      day=mod(dtpr,100)
      zulu=mod(tmpr/1000000,100)
      second=mod(tmpr/10000,100)*60+mod(tmpr/100,100)
      if (.not.validtm(year,month,day,zulu)) goto 4
      if ((month .eq. 2) .and. (day .eq. 29)) then
         if (.not. bissextile(year)) goto 4
      endif
      dat1=(jd(year,month,day)-jd1980)*17280+zulu*720+second/5
      return

!
!     mode=-1 : from (true_date and run_number) to stamp
!
 3    tdate=dat1
      runnb=dat3
 33   if((runnb.gt.9) .or. (.not.validtd(tdate))) goto 4
!     use new stamp if > jan 1, 2000 or fractional hour
      if (tdate.ge.td2000 .or. mod(tdate,720).ne.0) then
!     encode it in a new date-time stamp, ignore run nb
         stamp=tdstart+(tdate/8)*10+mod(tdate,8)
      else
!     encode it in an old date-time stamp
         call datec(jd1900+(tdate-td1900)/17280,year,month,day)
         tdate=(tdate-td1900)/720*720+td1900
         zulu=mod(tdate-td1900,17280)/720
         stamp=month*10000000+day*100000+(year-1900)*1000+zulu*10 &
              +runnb
      endif
      dat2(1)=stamp
      return

!
!     mode=1 : from stamp(old or new) to (true_date and run_number)
!
 5    stamp=dat2(1)
      if (stamp.ge.tdstart) then
!     stamp is a new date-time stamp
         tdate=(stamp-tdstart)/10*8+mod(stamp-tdstart,10)
         runnb=0
      else if (stamp .lt. -1) then
        print *,'newdate error: mode 1, negative stamp'
        goto 4
      else
!     stamp is an old date-time stamp
         runnb=mod(stamp,10)
         zulu=mod(stamp/10,100)
         year=mod(stamp/1000,100)+1900
         day=mod(stamp/100000,100)
         month=mod(stamp/10000000,100)
         tdate=(jd(year,month,day)-jd1980)*17280+zulu*720
      endif
      if (.not.validtd(tdate)) goto 4
      dat1=tdate
      dat3=runnb
      return

!
!     mode=4 : from 14 word old style DATE array TO STAMP and array(14)
!
100   dat1=itdmag2(dat2)
      return
!
!     mode=-4 : from STAMP TO 14 word old style DATE array
!
101   dat2(14)=dat1
      call dmagtp2(dat2)
      return
!
!     mode=5 : from printable to extended stamp
!
102   continue
      dtpr=dat2(1)
      tmpr=dat3
      dat1=0
      year=mod(dtpr/10000,10000)
      month=mod(dtpr/100,100)
      day=mod(dtpr,100)
      zulu=mod(tmpr/1000000,100)
      minute=mod(tmpr/10000,100)
      if (.not.validtme(year,month,day,zulu)) goto 4
      if ((month .eq. 2) .and. (day .eq. 29)) then
         if (.not. bissextile(year)) goto 4
      endif
      tdate=jd(year,month,day)
      if (tdate < jd0 .or. tdate >= jd10k) then
         print *,'newdate error: date outside of supported range, ', &
                 'date =',dtpr
         goto 4
      endif
      tdate=(tdate-jd0)*24+zulu+minute/60
!        encode it in a new date-time stamp
      stamp=(tdate/8)*10+mod(tdate,8)
      date_unsigned=stamp + troisg
      ! implicit INTEGER(8) to INTEGER(4) transfer
      dat1=date_unsigned
      return
!
!     mode=-5 : from extended stamp to printable
!
103   continue
      stamp8=dat1 ; stamp8=and(masque32,stamp8)
      dat2(1)= 0 ; dat3=0
      date_unsigned = stamp8
!!!   if (w16(1).eq.0) date_unsigned = ishft(stamp8,-32)
      if (date_unsigned <  troisg .or. &
          date_unsigned >= troisg + max_offset) then
        print *,'newdate error: invalid stamp for mode -5, stamp=',stamp
        goto 4
      endif
      stamp=date_unsigned - troisg
      tdate=stamp/10*8+mod(stamp,10)
      call datec(jd0+tdate/24,year,month,day)
      zulu=mod(tdate,24)
      minute=0
      if (.not.validtme(year,month,day,zulu)) goto 4
      if ((month .eq. 2) .and. (day .eq. 29)) then
         if (.not. bissextile(year)) goto 4
      endif
      dat2(1)=year*10000+month*100+day
      dat3=zulu*1000000+minute*10000
      return
!
!     mode=-6 : from extended true date to stamp
!
104   continue
      tdate=dat1
      if (tdate == exception .or.      &      ! 1901010101
         (tdate/24+jd0) <  jd1900 .or. &
         (tdate/24+jd0) >= jd2236) then ! extended stamp
         stamp=(tdate/8)*10+mod(tdate,8)
         date_unsigned=stamp + troisg
         ! implicit INTEGER(8) to INTEGER(4) transfer
         dat2(1)=date_unsigned
      else                               ! (new or old) stamp
         runnb=0
         zulu=mod(tdate,24)
         tdate=(tdate/24+jd0-jd1980)*17280+zulu*720
         goto 33
      endif
      return
!
!     mode=6 : from stamp to extended true date
!
105   continue
      stamp=dat2(1)
      if (stamp .lt. -1) then
        stamp8=stamp ; stamp8=and(masque32,stamp8)
        dat1=0
        dat3=0
        date_unsigned = stamp8
!!!     if (w16(1).eq.0) date_unsigned = ishft(stamp8,-32)
        if (date_unsigned < troisg .or. &
           date_unsigned > troisg + max_offset) then
           print *,'newdate error: invalid stamp for mode -6, ', &
                   'stamp=',stamp
           goto 4
        endif
        stamp=date_unsigned - troisg
        tdate=stamp/10*8+mod(stamp,10)
        dat1=tdate
        dat3=0
      else
        if (stamp.ge.tdstart) then
!     stamp is a new date-time stamp
          tdate=(stamp-tdstart)/10*8+mod(stamp-tdstart,10)
          call datec(jd1900+(tdate-td1900)/17280,year,month,day)
          zulu=mod(tdate-td1900,17280)/720
          runnb=0
          tdate=(jd(year,month,day)-jd0)*24+zulu
!         print *,'Debug stamp > tdstart tdate=',tdate,
!     %           ' zulu=',zulu
        else
!     stamp is an old date-time stamp
           runnb=mod(stamp,10)
           zulu=mod(stamp/10,100)
           year=mod(stamp/1000,100)+1900
           day=mod(stamp/100000,100)
           month=mod(stamp/10000000,100)
           tdate=(jd(year,month,day)-jd0)*24+zulu
!           print *,'Debug old date stamp tdate=',tdate
        endif
        if (.not.validtd(tdate)) goto 4
        dat1=tdate
        dat3=runnb
      endif
      return
!
!     mode=-7 : from extended true_date to printable
!
106   tdate=dat1
      if (.not.validtd(tdate)) goto 4
      call datec(jd0+tdate/24,year,month,day)
      zulu=mod(tdate,24)
      if (.not.validtme(year,month,day,zulu)) goto 4
      if ((month .eq. 2) .and. (day .eq. 29)) then
         if (.not. bissextile(year)) goto 4
      endif
      minute=0
      dat2(1)=year*10000+month*100+day
      dat3=zulu*1000000+minute*10000
      return
!
!     mode=7 : from printable to extended true_date
!
107   dtpr=dat2(1)
      tmpr=dat3
      year=mod(dtpr/10000,10000)
      if (year < 0 .or. year >= 10000) then
         print *,'newdate error: date outside of supported range, ', &
                 'date =',dtpr
         goto 4
      endif
      month=mod(dtpr/100,100)
      day=mod(dtpr,100)
      zulu=mod(tmpr/1000000,100)
      second=mod(tmpr/10000,100)*60+mod(tmpr/100,100)
      if (.not.validtme(year,month,day,zulu)) goto 4
      if ((month .eq. 2) .and. (day .eq. 29)) then
         if (.not. bissextile(year)) goto 4
      endif
      dat1=(jd(year,month,day)-jd0)*24+zulu
      return
!
!     error: bad mode or bad arguments
!
 4    naetwed=1
      return
      end

