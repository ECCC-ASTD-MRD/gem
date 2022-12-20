!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

!#TODO: rewrite in C with iso_c itf for F90

!/@*
module timestr_mod
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use clib_itf_mod, only: clib_toupper
   use str_mod, only: str_toreal,str_normalize
   use mu_jdate_mod, only: jdate_from_cmc, jdate_year, jdate_month
   implicit none
   private
   !@objective
   !@description
   ! Public functions
   public :: timestr_parse,timestr_check,timestr_isstep,timestr2sec,timestr2step, &
        timestr_prognum,timestr_unitfact,timestr_default_set
   ! Public constants
   integer,public :: TIMESTR_NO_MATCH = 0
   integer,public :: TIMESTR_MATCH = 1
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   interface timestr_parse
      module procedure timestr_parse_multi
      module procedure timestr_parse_1
   end interface timestr_parse

   interface timestr2sec
      module procedure timestr2sec_r4
      module procedure timestr2sec_r8
      module procedure timestr2sec_i4
      module procedure timestr2sec_i8
   end interface timestr2sec

   interface timestr2step
      module procedure timestr2step_r4
      module procedure timestr2step_r8
      module procedure timestr2step_i4
      module procedure timestr2step_i8
   end interface timestr2step

   interface timestr_prognum
      module procedure timestr_prognum_int
      module procedure timestr_prognum_int_jdate
   end interface timestr_prognum

   interface timestr_isstep
      module procedure timestr_isstep_str
      module procedure timestr_isstep_int
      module procedure timestr_isstep_str_jdate
      module procedure timestr_isstep_int_jdate
   end interface timestr_isstep

   real(REAL64),parameter :: EPSILON_8 = tiny(1.D0)
   real,parameter :: EPSILON_4 = tiny(1.)

   character(len=64) :: m_timestr_default_S = 'steps,-1'

contains

   !/@*
   subroutine timestr_default_set(F_default_S)
      !@objective Define default timestring to use if provided one is empty
      implicit none
      !@arguments
      character(len=*), intent(in)  :: F_default_S
      !*@/
      m_timestr_default_S = F_default_S
   end subroutine timestr_default_set


   !/@*
   function timestr_parse_multi(F_values,F_loop,F_hasloop_L,F_units_S,F_timestr_S,F_default_S) result(F_nvals)
      !@objective Parse time string into number and units
      implicit none
      !@arguments
      real,             intent(out) :: F_values(:)
      real,             intent(out) :: F_loop(3)
      logical,          intent(out) :: F_hasloop_L
      character(len=*), intent(out) :: F_units_S
      character(len=*), intent(in)  :: F_timestr_S
      character(len=*), intent(in), optional  :: F_default_S
      !@return
      integer :: F_nvals !# Number of values (-1 on error)
      !@description
      ! The time string takes the following formats
      !
      ! Old format:
      !    999U
      !    Where:
      !       999: time value, any fortran accepted real number
      !       U  : units of given "time value" (optional), one of:
      !            P for steps (if dt not provided, equivalent to S)
      !            S for seconds
      !            M for minutes
      !            H for hours
      !            D for days
      !       If U if not provided a default unit is used,
      !          F_default_S or 'hours' is not provided
      !
      ! New formats, accepted formats:
      !    VALUE
      !    UNITS,VALUE
      !    UNITS,VALUELIST
      !    UNITS,[VALUELIST]
      !    UNITS,<VALUE0,VALUEEND,INTERVAL>
      !    TODO: accept mix of the above
      !    Where :
      !    UNITS     : units of given time value (optional)
      !                accepted values (case insensitive):
      !                     steps, seconds, minutes, hours, days, months
      !                Only the first 3 letters are checked
      !                If not provided default units are used,
      !                     F_default_S or 'hours' is not provided
      !    VALUE     : time value, any fortran accepted real number
      !    VALUELIST : list of time values, any fortran accepted real number
      !    VALUE0,VALUEEND,INTERVAL : loop from VALUE0 to VALUEEND with strides of INTERVAL
      !
      ! Returned F_units_S is expressed as 3 lettres like in new format
      !*@/
      integer,parameter :: NMAX0 = 1024
      integer :: mylen,istat,nmax,ii
      character(len=64) :: nbr_S,units_S,default_S,timestr_S,parts_S(NMAX0)

      F_nvals  = RMN_ERR
      F_values = -1.
      F_loop(:) = -1.
      F_hasloop_L = .false.
      F_units_S = ' '
      nmax = min(NMAX0,size(F_values))

      timestr_S = F_timestr_S
      mylen = len_trim(timestr_S)
      if (mylen <= 0) then
         timestr_S = m_timestr_default_S
         mylen = len_trim(timestr_S)
         if (mylen <= 0) return
      endif

      default_S = 'HOURS'
      if (present(F_default_S)) then
         default_S = F_default_S
         istat = clib_toupper(default_S)
      endif

      istat = clib_toupper(timestr_S)
      call str_split(units_S,nbr_S,timestr_S,',')

      !# old style, no comma
      if (units_S == timestr_S) then
         units_S = timestr_S(mylen:mylen)
         select case(units_S(1:1))
         case('P') !pas
            units_S = 'STEPS'
         case('S')
            units_S = 'SECONDS'
         case('M')
            units_S = 'MINUTES'
         case('H')
            units_S = 'HOURS'
         case('D')
            units_S = 'DAYS'
         case default
            !TODO: more clever check for other cases, raise error if not known
            units_S = ' '
            mylen = mylen + 1
         end select
         if (mylen <= 1) return
         nbr_S = timestr_S(1:mylen-1)
      endif
      if (units_S == ' ') units_S = default_S

      if (.not.any(units_S(1:3) == (/'STE','SEC','MIN','HOU','DAY','MON'/))) then
         call msg(MSG_WARNING,'(timestr_parse) Invalid time string: '//trim(F_timestr_S))
         return
      endif
      F_units_S = units_S(1:3)

      call str_normalize(nbr_S)
      mylen = len_trim(nbr_S)
      if (mylen <= 0)  then
         call msg(MSG_WARNING,'(timestr_parse) Invalid time string: '//trim(F_timestr_S))
         return
      endif

      parts_S(:) = ' '
      if (nbr_S(1:1) == '<' .and. nbr_S(mylen:mylen) == '>') then
         nbr_S = nbr_S(2:mylen-1)
         call str_split2list(parts_S,nbr_S,',',3)
         istat = str_toreal(F_loop(1),parts_S(1))
         istat = min(str_toreal(F_loop(2),parts_S(2)),istat)
         F_loop(3) = 1.
         if (parts_S(3) /= ' ') &
              istat = min(str_toreal(F_loop(3),parts_S(3)),istat)
         if (.not.RMN_IS_OK(istat))  then
            call msg(MSG_WARNING,'(timestr_parse) Invalid time string: '//trim(F_timestr_S))
            return
         endif
         F_hasloop_L = .true.
         F_nvals = 0
      else
         if (nbr_S(1:1) == '[' .and. nbr_S(mylen:mylen) == ']') then
            nbr_S = nbr_S(2:mylen-1)
         endif
         call str_split2list(parts_S,nbr_S,',',min(nmax+1,NMAX0))
         istat = RMN_OK
         ii = 1
         do while (RMN_IS_OK(istat) .and. ii <= nmax .and. parts_S(ii) /= ' ')
            istat = str_toreal(F_values(ii),parts_S(ii))
            ii = ii + 1
         enddo
         if (parts_S(min(ii,NMAX0)) /= ' ') then
            call msg(MSG_WARNING,'(timestr_parse) Too many values to unpack in time string: '//trim(F_timestr_S))
            return
         endif
         if (.not.RMN_IS_OK(istat))  then
            call msg(MSG_WARNING,'(timestr_parse) Invalid time string: '//trim(F_timestr_S))
            return
         endif
         F_nvals = ii - 1
      endif

     return
  end function timestr_parse_multi


   !/@*
   function timestr_parse_1(F_nbr,F_units_S,F_timestr_S,F_default_S) result(F_status)
      !@objective Parse time string into number and units
      implicit none
      !@arguments
      real,             intent(out) :: F_nbr
      character(len=*), intent(out) :: F_units_S
      character(len=*), intent(in)  :: F_timestr_S
      character(len=*), intent(in), optional  :: F_default_S
      !@return
      integer :: F_status
      !*@/
      real :: values(1),loop(3)
      logical :: hasloop_L

      F_nbr     = -1.
      if (present(F_default_S)) then
         F_status  = timestr_parse_multi(values,loop,hasloop_L,F_units_S,F_timestr_S,F_default_S)
      else
         F_status  = timestr_parse_multi(values,loop,hasloop_L,F_units_S,F_timestr_S)
      endif
      if (F_status /= 1 .or. hasloop_L) F_status = RMN_ERR
      F_nbr = values(1)

      return
   end function timestr_parse_1


   !/@*
   function timestr_check(F_timestr_S) result(F_status)
      !@objective Check if timestr is valid
      implicit none
      !@arguments
      character(len=*), intent(in)  :: F_timestr_S
      !@return
      integer :: F_status
      !*@/
      character(len=64) :: units_S
      real :: values(128),loop(3)
      logical :: hasloop_L
      F_status  = timestr_parse_multi(values,loop,hasloop_L,units_S,F_timestr_S)
      return
   end function timestr_check


   !/@*
   function timestr2sec_r8(F_sec,F_timestr_S,F_dt,F_default_S) result(F_status)
      !@objective Parse time string and return number of seconds equivalent
      implicit none
      !@arguments
      real(REAL64),    intent(out) :: F_sec
      character(len=*), intent(in)  :: F_timestr_S
      real(REAL64),    intent(in), optional  :: F_dt !# Timestep length [sec]
      character(len=*), intent(in), optional  :: F_default_S
      !@return
      integer :: F_status
      !*@/
      real :: nunits
      real(REAL64) :: dt,fact_8
      character(len=64) :: units_S

      F_sec = 0.
      if (present(F_default_S)) then
         F_status = timestr_parse(nunits,units_S,F_timestr_S,F_default_S)
      else
         F_status = timestr_parse(nunits,units_S,F_timestr_S)
      endif
      if (.not.RMN_IS_OK(F_status)) return

      dt = 1.
      if (present(F_dt)) dt = F_dt

      fact_8 = timestr_unitfact(units_S,dt) * dt
!!$      if (F_units_S(1:3) == 'MON') !#TODO: need to know what month/year
      if (fact_8 < 0.) then
         F_status = RMN_ERR
         return
      endif

      F_sec = dble(nunits) * fact_8

      return
   end function timestr2sec_r8


   !/@*
   function timestr2sec_r4(F_sec,F_timestr_S,F_dt,F_default_S) result(F_status)
      !@objective Parse time string and return number of seconds equivalent
      implicit none
      !@arguments
      real,    intent(out) :: F_sec
      character(len=*), intent(in)  :: F_timestr_S
      real(REAL64),    intent(in), optional  :: F_dt !# Timestep length [sec]
      character(len=*), intent(in), optional  :: F_default_S
      !@return
      integer :: F_status
      !*@/
      real(REAL64) :: dt,sec_8

      F_sec = 0.
      dt = 1.
      if (present(F_dt)) dt = F_dt
      if (present(F_default_S)) then
         F_status = timestr2sec_r8(sec_8,F_timestr_S,dt,F_default_S)
      else
         F_status = timestr2sec_r8(sec_8,F_timestr_S,dt)
      endif
      if (.not.RMN_IS_OK(F_status)) return
      F_sec = real(sec_8)

      return
   end function timestr2sec_r4


   !/@*
   function timestr2sec_i4(F_sec,F_timestr_S,F_dt,F_default_S) result(F_status)
      !@objective Parse time string and return number of seconds equivalent
      implicit none
      !@arguments
      integer,    intent(out) :: F_sec
      character(len=*), intent(in)  :: F_timestr_S
      real(REAL64),    intent(in), optional  :: F_dt !# Timestep length [sec]
      character(len=*), intent(in), optional  :: F_default_S
      !@return
      integer :: F_status
      !*@/
      real(REAL64) :: dt,sec_8

      F_sec = 0.
      dt = 1.
      if (present(F_dt)) dt = F_dt
      if (present(F_default_S)) then
         F_status = timestr2sec_r8(sec_8,F_timestr_S,dt,F_default_S)
      else
         F_status = timestr2sec_r8(sec_8,F_timestr_S,dt)
      endif
      if (.not.RMN_IS_OK(F_status)) return
      F_sec = int(sec_8)

      return
   end function timestr2sec_i4


   !/@*
   function timestr2sec_i8(F_sec,F_timestr_S,F_dt,F_default_S) result(F_status)
      !@objective Parse time string and return number of seconds equivalent
      implicit none
      !@arguments
      integer(INT64), intent(out) :: F_sec
      character(len=*), intent(in)  :: F_timestr_S
      real(REAL64),    intent(in), optional  :: F_dt !# Timestep length [sec]
      character(len=*), intent(in), optional  :: F_default_S
      !@return
      integer :: F_status
      !*@/
      real(REAL64) :: dt,sec_8

      F_sec = 0.
      dt = 1.
      if (present(F_dt)) dt = F_dt
      if (present(F_default_S)) then
         F_status = timestr2sec_r8(sec_8,F_timestr_S,dt,F_default_S)
      else
         F_status = timestr2sec_r8(sec_8,F_timestr_S,dt)
      endif
      if (.not.RMN_IS_OK(F_status)) return
      F_sec = sec_8

      return
   end function timestr2sec_i8


   !/@*
   function timestr2step_r8(F_nstep,F_timestr_S,F_dt,F_default_S) result(F_status)
      !@objective Parse time string and return number of step equivalent
      implicit none
      !@arguments
      real(REAL64),    intent(out) :: F_nstep
      character(len=*), intent(in)  :: F_timestr_S
      real(REAL64),    intent(in)  :: F_dt !# Timestep length [sec]
      character(len=*), intent(in), optional  :: F_default_S
      !@return
      integer :: F_status
      !*@/
      real(REAL64) :: sec_8

      F_nstep  = -1
      if (present(F_default_S)) then
         F_status = timestr2sec_r8(sec_8,F_timestr_S,F_dt,F_default_S)
      else
         F_status = timestr2sec_r8(sec_8,F_timestr_S,F_dt)
      endif
      if (.not.RMN_IS_OK(F_status)) return

      if (abs(F_dt) < EPSILON_8) then
         call msg(MSG_ERROR,"(timestr2step) Cannot use provided timestep == 0.")
         F_status = RMN_ERR
         return
      endif
      F_nstep = sec_8 / F_dt

      return
   end function timestr2step_r8


   !/@*
   function timestr2step_r4(F_nstep,F_timestr_S,F_dt,F_default_S) result(F_status)
      !@objective Parse time string and return number of step equivalent
      implicit none
      !@arguments
      real,             intent(out) :: F_nstep
      character(len=*), intent(in)  :: F_timestr_S
      real(REAL64),    intent(in)  :: F_dt !# Timestep length [sec]
      character(len=*), intent(in), optional  :: F_default_S
      !@return
      integer :: F_status
      !*@/
      real(REAL64) :: sec_8

      F_nstep  = -1
      if (present(F_default_S)) then
         F_status = timestr2sec_r8(sec_8,F_timestr_S,F_dt,F_default_S)
      else
         F_status = timestr2sec_r8(sec_8,F_timestr_S,F_dt)
      endif
      if (.not.RMN_IS_OK(F_status)) return

      if (abs(F_dt) < EPSILON_8) then
         call msg(MSG_ERROR,"(timestr2step) Cannot use provided timestep == 0.")
         F_status = RMN_ERR
         return
      endif
      F_nstep = sec_8 / F_dt

      return
   end function timestr2step_r4


   !/@*
   function timestr2step_i4(F_nstep,F_timestr_S,F_dt,F_default_S) result(F_status)
      !@objective Parse time string and return number of step equivalent
      implicit none
      !@arguments
      integer,          intent(out) :: F_nstep
      character(len=*), intent(in)  :: F_timestr_S
      real(REAL64),    intent(in)  :: F_dt !# Timestep length [sec]
      character(len=*), intent(in), optional  :: F_default_S
      !@return
      integer :: F_status
      !*@/
      real(REAL64) :: sec_8

      F_nstep  = -1
      if (present(F_default_S)) then
         F_status = timestr2sec_r8(sec_8,F_timestr_S,F_dt,F_default_S)
      else
         F_status = timestr2sec_r8(sec_8,F_timestr_S,F_dt)
      endif
      if (.not.RMN_IS_OK(F_status)) return

      if (abs(F_dt) < EPSILON_8) then
         call msg(MSG_ERROR,"(timestr2step) Cannot use provided timestep == 0.")
         F_status = RMN_ERR
         return
      endif
      F_nstep = sec_8 / F_dt

      return
   end function timestr2step_i4


   !/@*
   function timestr2step_i8(F_nstep,F_timestr_S,F_dt,F_default_S) result(F_status)
      !@objective Parse time string and return number of step equivalent
      implicit none
      !@arguments
      integer(INT64), intent(out) :: F_nstep
      character(len=*), intent(in)  :: F_timestr_S
      real(REAL64),    intent(in)  :: F_dt !# Timestep length [sec]
      character(len=*), intent(in), optional  :: F_default_S
      !@return
      integer :: F_status
      !*@/
      real(REAL64) :: sec_8

      F_nstep  = -1
      if (present(F_default_S)) then
         F_status = timestr2sec_r8(sec_8,F_timestr_S,F_dt,F_default_S)
      else
         F_status = timestr2sec_r8(sec_8,F_timestr_S,F_dt)
      endif
      if (.not.RMN_IS_OK(F_status)) return

      if (abs(F_dt) < EPSILON_8) then
         call msg(MSG_ERROR,"(timestr2step) Cannot use provided timestep == 0.")
         F_status = RMN_ERR
         return
      endif
      F_nstep = sec_8 / F_dt

      return
   end function timestr2step_i8


   !/@*
   function timestr_prognum_int(F_prognum,F_units_S,F_interval,F_dateo,F_dt,F_step,F_maxstep) result(F_status)
      !@objective Return interger for program-number
      implicit none
      !@arguments
      integer,          intent(out) :: F_prognum
      character(len=*), intent(in)  :: F_units_S
      real,             intent(in)  :: F_interval
      integer,          intent(in)  :: F_dateo    !# Date of origin [CMC date]
      integer,          intent(in)  :: F_step
      real,             intent(in)  :: F_dt       !# Timestep length [sec]
      integer,optional, intent(in)  :: F_maxstep
      !@return
      integer :: F_status
      !*@/
      integer(INT64) :: jdateo
      jdateo = jdate_from_cmc(F_dateo)
      if (present(F_maxstep)) then
         F_status = timestr_prognum_int_jdate(F_prognum,F_units_S,F_interval,jdateo,F_dt,F_step,F_maxstep)
      else
         F_status = timestr_prognum_int_jdate(F_prognum,F_units_S,F_interval,jdateo,F_dt,F_step)
      endif
      return
   end function timestr_prognum_int


   !/@*
   function timestr_prognum_int_jdate(F_prognum,F_units_S,F_interval,F_jdateo,F_dt,F_step,F_maxstep) result(F_status)
      !@objective Return interger for program-number
      implicit none
      !@arguments
      integer,          intent(out) :: F_prognum
      character(len=*), intent(in)  :: F_units_S
      real,             intent(in)  :: F_interval
      integer(INT64), intent(in)  :: F_jdateo   !# Date of origin [Julian Sec]
      integer,          intent(in)  :: F_step
      real,             intent(in)  :: F_dt       !# Timestep length [sec]
      integer,optional, intent(in)  :: F_maxstep
      !@return
      integer :: F_status
      !*@/
      real(REAL64), parameter :: EPSILON_8b = 1.0D-12
      integer :: interval, mystep,m0,m1,y0,y1,nmonths
      integer(INT64) :: jdatev2, istep, dt
      real(REAL64) :: fact_8

      F_status  = RMN_ERR
      F_prognum = 0

      if (nint(F_interval) <= 0.) return

      fact_8 = timestr_unitfact(F_units_S,dble(F_dt))
      if (F_units_S(1:3) == 'MON') fact_8 = 1.D0
      if (fact_8 < 0.) then
         call msg(MSG_ERROR,"(timestr_prognum) Unknown units: "//trim(F_units_S))
         return
      endif

      interval = nint(F_interval*fact_8)
      if (interval == 0) then
         call msg(MSG_ERROR,"(timestr_prognum) Cannot use interval == 0")
         return
      endif

      if (F_units_S(1:3) == 'MON') then

         istep = F_step ; dt = F_dt
         jdatev2 = F_jdateo + (istep * dt)
         m0 = jdate_month(F_jdateo)
         y0 = jdate_year(F_jdateo)
         m1 = jdate_month(jdatev2)
         y1 = jdate_year(jdatev2)

         nmonths = (y1-y0)*12 + m1-m0
         F_prognum = (nmonths/interval + min(0,mod(nmonths,interval)))*interval

      else

         mystep = F_step/interval
         if (F_step>0) mystep= mystep + min(1,mod(F_step,interval))
         mystep = mystep * interval
         if (present(F_maxstep)) mystep = min(mystep,F_maxstep)
         F_prognum = ceiling(dble(mystep)/fact_8 - EPSILON_8b)

      endif

      F_status = RMN_OK
      return
   end function timestr_prognum_int_jdate


   !/@*
   function timestr_isstep_str(F_timestr_S,F_dateo,F_dt,F_step,F_maxstep,F_default_S) result(F_status)
      !@objective Return TIMESTR_MATCH if step match F_timestr_S or last step of end month
      implicit none
      !@arguments
      character(len=*), intent(in)  :: F_timestr_S
      integer,          intent(in)  :: F_dateo,F_step
      real,             intent(in)  :: F_dt
      integer,          intent(in), optional  :: F_maxstep
      character(len=*), intent(in), optional  :: F_default_S
      !@return
      integer :: F_status
      !*@/
      character(len=64) :: units_S
      real :: values(128),loop(3)
      logical :: hasloop_L

      if (present(F_default_S)) then
         F_status  = timestr_parse_multi(values,loop,hasloop_L,units_S,F_timestr_S,F_default_S)
      else
         F_status  = timestr_parse_multi(values,loop,hasloop_L,units_S,F_timestr_S)
      endif
      if (.not.RMN_IS_OK(F_status)) return

      if (hasloop_L .or. F_status /= 1) then
         !#TODO: implement for other time specifications
         F_status = RMN_ERR
         call msg(MSG_ERROR,'(timestr_isstep) Not yet implemented for loops and lists')
         return
      endif

      if (present(F_maxstep)) then
         F_status  = timestr_isstep_int(units_S,values(1),F_dateo,F_dt,F_step,F_maxstep)
      else
         F_status  = timestr_isstep_int(units_S,values(1),F_dateo,F_dt,F_step)
      endif

      return
   end function timestr_isstep_str


   !/@*
   function timestr_isstep_str_jdate(F_timestr_S,F_jdateo,F_dt,F_step,F_maxstep,F_default_S) result(F_status)
      !@objective Return TIMESTR_MATCH if step match F_timestr_S or last step of end month
      implicit none
      !@arguments
      character(len=*), intent(in)  :: F_timestr_S
      integer(INT64), intent(in)  :: F_jdateo
      integer,          intent(in)  :: F_step
      real,             intent(in)  :: F_dt
      integer,          intent(in), optional  :: F_maxstep
      character(len=*), intent(in), optional  :: F_default_S
      !@return
      integer :: F_status
      !*@/
      character(len=64) :: units_S
      real :: values(128),loop(3)
      logical :: hasloop_L

      if (present(F_default_S)) then
         F_status  = timestr_parse_multi(values,loop,hasloop_L,units_S,F_timestr_S,F_default_S)
      else
         F_status  = timestr_parse_multi(values,loop,hasloop_L,units_S,F_timestr_S)
      endif
      if (.not.RMN_IS_OK(F_status)) return

      if (hasloop_L .or. F_status /= 1) then
         !#TODO: implement for other time specifications
         F_status = RMN_ERR
         call msg(MSG_ERROR,'(timestr_isstep) Not yet implemented for loops and lists')
         return
      endif

      if (present(F_maxstep)) then
         F_status  = timestr_isstep_int_jdate(units_S,values(1),F_jdateo,F_dt,F_step,F_maxstep)
      else
         F_status  = timestr_isstep_int_jdate(units_S,values(1),F_jdateo,F_dt,F_step)
      endif

      return
   end function timestr_isstep_str_jdate


  !/@*
   function timestr_isstep_int(F_units_S,F_interval,F_dateo,F_dt,F_step,F_maxstep) result(F_status)
      !@objective Return TIMESTR_MATCH if step match interval or last step of end month
      implicit none
      !@arguments
      character(len=*), intent(in)  :: F_units_S
      real,             intent(in)  :: F_interval
      integer,          intent(in)  :: F_dateo,F_step
      real,             intent(in)  :: F_dt
      integer,optional, intent(in)  :: F_maxstep
      !@return
      integer :: F_status
      !*@/
      integer :: interval,istat,prognum,prognum1,maxstep,gap
      real(REAL64) :: fact_8, ris_8, seconds_8

      F_status  = RMN_ERR
      maxstep = F_step+9999
      if (present(F_maxstep)) then
         if (F_step >= F_maxstep) then
            F_status = TIMESTR_MATCH
            return
         endif
         maxstep = F_maxstep
      endif

      if (nint(F_interval) <= 0.) return

      IF_MON: if (F_units_S(1:3) /= 'MON') then

         fact_8 = timestr_unitfact(F_units_S,dble(F_dt))
         if (fact_8 < 0.) then
            call msg(MSG_ERROR,"(timestr_isstep) Unknown units: "//trim(F_units_S))
            return
         endif

         interval = nint(F_interval*fact_8)
         F_status = TIMESTR_NO_MATCH
         if (mod(F_step,interval) == 0) F_status = TIMESTR_MATCH

      else !# IF_MON

         if (abs(F_dt) < EPSILON_4) then
            call msg(MSG_ERROR,"(timestr_isstep) Cannot use provided timestep == 0.")
            return
         endif

         gap = 1800/F_dt
         gap = gap + 1


         istat = timestr_prognum_int(prognum,F_units_S,F_interval,F_dateo,F_dt,F_step-gap,maxstep)
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_ERROR,"(timestr_isstep) Unknown units: "//trim(F_units_S))
            return
         endif
         istat = timestr_prognum_int(prognum1,F_units_S,F_interval,F_dateo,F_dt,F_step,maxstep)
         F_status = TIMESTR_NO_MATCH

         seconds_8 = dble(F_step) * Dble(F_dt)
         ris_8 = dmod(seconds_8,3600.d0)
         if (prognum1 > prognum) then
            if(ris_8 >= 1800.d0 .or. F_step == 1) then
               F_status = TIMESTR_NO_MATCH
            else
               F_status = TIMESTR_MATCH
            endif
         endif

      endif IF_MON

      return
   end function timestr_isstep_int


   !/@*
   function timestr_isstep_int_jdate(F_units_S,F_interval,F_jdateo,F_dt,F_step,F_maxstep) result(F_status)
      !@objective Return TIMESTR_MATCH if step match interval or last step of end month
      implicit none
      !@arguments
      character(len=*), intent(in)  :: F_units_S
      real,             intent(in)  :: F_interval
      integer(INT64), intent(in)  :: F_jdateo
      integer,          intent(in)  :: F_step
      real,             intent(in)  :: F_dt
      integer,optional, intent(in)  :: F_maxstep
      !@return
      integer :: F_status
      !*@/
      integer :: interval,istat,prognum,prognum1,maxstep
      real(REAL64) :: fact_8

      F_status  = RMN_ERR
      maxstep = F_step+9999
      if (present(F_maxstep)) then
         if (F_step >= F_maxstep) then
            F_status = TIMESTR_MATCH
            return
         endif
         maxstep = F_maxstep
      endif

      if (nint(F_interval) <= 0.) return

      IF_MON: if (F_units_S(1:3) /= 'MON') then

         fact_8 = timestr_unitfact(F_units_S,dble(F_dt))
         if (fact_8 < 0.) then
            call msg(MSG_ERROR,"(timestr_isstep) Unknown units: "//trim(F_units_S))
            return
         endif

         interval = nint(F_interval*fact_8)
         F_status = TIMESTR_NO_MATCH
         if (mod(F_step,interval) == 0) F_status = TIMESTR_MATCH

      else !# IF_MON

         istat = timestr_prognum(prognum,F_units_S,F_interval,F_jdateo,F_dt,F_step,maxstep)
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_ERROR,"(timestr_isstep) Unknown units: "//trim(F_units_S))
            return
         endif
         istat = timestr_prognum(prognum1,F_units_S,F_interval,F_jdateo,F_dt,F_step+1,maxstep)
         F_status = TIMESTR_NO_MATCH
         if (prognum1 > prognum)  F_status = TIMESTR_MATCH
!!$         if (prognum1 > prognum .and. F_step/=1)  F_status = TIMESTR_MATCH

      endif IF_MON

      return
   end function timestr_isstep_int_jdate


   !/@*
   function timestr_unitfact(F_units_S,F_dt_8) result(F_fact_8)
      !@objective Return factor to convert to Steps
      implicit none
      !@arguments
      character(len=*),       intent(in)  :: F_units_S
      real(REAL64),optional, intent(in)  :: F_dt_8
      !@return
      real(REAL64) :: F_fact_8
      !*@/
      real(REAL64) :: dt_8

      F_fact_8 = -1.D0
      dt_8   = 1.D0
      if (present(F_dt_8)) dt_8 = F_dt_8
      if (abs(dt_8) < EPSILON_8) then
         call msg(MSG_ERROR,"(timestr_unitfact) Cannot use provided timestep == 0.")
         return
      endif

      select case(F_units_S(1:3))
      case('STE')
         F_fact_8 = 1.D0
      case('SEC')
         F_fact_8 = 1.D0/dt_8
      case('MIN')
         F_fact_8 = 60.D0/dt_8
      case('HOU')
         F_fact_8 = 3600.D0/dt_8
      case('DAY')
         F_fact_8 = 86400.D0/dt_8
!!$      case('MON')
!!$         F_fact_8 = 1.D0
!!$      case default
!!$         call msg(MSG_ERROR,"(timestr) Unknown units: "//trim(F_units_S))
!!$         return
      end select
      return
   end function timestr_unitfact


   !#TODO: timestr_isStep(): isstep, next-interval-len, next-interval-step/time/date?, see outcfg_time

end module timestr_mod
