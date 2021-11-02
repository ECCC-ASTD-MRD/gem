! RMNLIB - Library of useful routines for C and FORTRAN programming
! Copyright (C) 1975-2019  Division de Recherche en Prevision Numerique
!                          Environnement Canada
! 
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation,
! version 2.1 of the License.
! 
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the
! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! Boston, MA 02111-1307, USA.

#if defined(SELF_TEST)
program test
  implicit none
  character(len=64) :: timezone
  integer :: yyyy, mo, dd, hh, mi, ss, ms
  integer :: yyyymmdd, hhmmss00

  CALL GET_ENVIRONMENT_VARIABLE('TZ',timezone)
  write(*,*) 'timezone = ',trim(timezone)
  call world_date_and_time(yyyy, mo, dd, hh, mi, ss, ms, .true.)
  write(*,100) 'GMT   =',yyyy, mo, dd, hh, mi, ss, ms
  call world_date_and_time(yyyy, mo, dd, hh, mi, ss, ms, .false.)
  write(*,100) 'LOCAL =',yyyy, mo, dd, hh, mi, ss, ms
  call system_time(yyyymmdd,hhmmss00)
  write(*,*) 'system_time =',yyyymmdd,hhmmss00
100 format(A,8I6)
stop
end
#endif
subroutine system_time(yyyymmdd,hhmmss00)   ! return GMT date and time
  implicit none
  integer, intent(OUT) :: yyyymmdd, hhmmss00

  integer :: yyyy, mo, dd, hh, mi, ss, ms
  call world_date_and_time(yyyy, mo, dd, hh, mi, ss, ms, .true.)
  yyyymmdd = yyyy*10000 + mo*100 + dd
  hhmmss00 = hh*1000000 + mi*10000 + ss*100

  return
end subroutine system_time

subroutine world_date_and_time(yyyy, mo, dd, hh, mi, ss, ms, gmt)  ! return GMT or local date and time
  implicit none
  integer, intent(OUT) :: yyyy, mo, dd, hh, mi, ss, ms             ! year, month, day, hours, minutes, seconds, milliseconds
  logical, intent(IN)  :: gmt     ! if .true. get GMT time, else get local time

  integer, dimension(8) :: values
  integer :: jd
  integer :: sgn
  integer :: hrs, mins, daybump

  call date_and_time(VALUES=values)

  if(gmt) then  ! account for time offset if we want GMT
    sgn = 1
    if(values(4) < 0) sgn = -1       ! sign of time offset
    hrs = abs(values(4)) / 60        ! offset hours
    mins = mod(abs(values(4)),60)    ! offset minutes

    values(6) = values(6) - sgn * mins  ! correction for minutes
    values(5) = values(5) - sgn * hrs   ! correction for hours

    if(values(6) >= 60) then            ! minutes > 60, fix hours and minutes
      values(5) = values(5) + 1         ! + 1 hour
      values(6) = values(6) - 60        ! fix minutes
    endif
    if(values(6) < 0) then              ! minutes < 0, fix hours and minutes
      values(5) = values(5) - 1         ! - 1 hour
      values(6) = values(6) + 60        ! fix minutes
    endif

    daybump = 0
    if(values(5) >= 24) then            ! hours >= 24, fix hours, compute day fix
      daybump = 1
      values(5) = values(5) - 24
    endif
    if(values(5) < 0) then              ! hours < 0, fix hours, compute day fix
      daybump = -1
      values(5) = values(5) + 24
    endif

    if(daybump .ne. 0) then             ! must add or substract one day
      call JDATEC(JD,values(1),values(2),values(3))  ! julian day from y/m/d
      jd = jd + daybump
      call DATEC(JD,values(1),values(2),values(3))   ! y/m/d from julian day
    endif
  endif

  yyyy = values(1)    ! year
  mo   = values(2)    ! month
  dd   = values(3)    ! day
  hh   = values(5)    ! hours
  mi   = values(6)    ! minutes
  ss   = values(7)    ! seconds
  ms   = values(8)    ! milliseconds

  return
contains    ! jdatec and datec included to make this self contained

! NOTES    - ALGORITHM COPIED FROM "COMMUNICATIONS OF THE ACM" (1968), PAGE 657.
!          - IT COVERS A PERIOD OF 7980 YEARS WITH DAY 1 STARTING
!            AT YEAR=-4713, MONTH=11, DAY=25.

  SUBROUTINE JDATEC(JD,I,J,K)  ! convert 3 integers (year, month, day) into julian day
    integer, intent(IN) :: I,J,K
    integer, intent(OUT) :: jd
    JD        = K-32075+1461*(I+4800+(J-14)/12)/4   &
      &             +367*(J-2-(J-14)/12*12)/12-3   &
      &             *((I+4900+(J-14)/12)/100)/4

  return
  end
  SUBROUTINE DATEC(JD,I,J,K) ! convert julian day into 3 integers (year, month, day)
    integer, intent(OUT) :: I,J,K
    integer, intent(IN) :: jd
    integer :: L, N
    L= JD+68569
    N= 4*L/146097
    L= L-(146097*N+3)/4
    I= 4000*(L+1)/1461001
    L= L-1461*I/4+31
    J= 80*L/2447
    K= L-2447*J/80
    L= J/11
    J= J+2-12*L
    I= 100*(N-49)+I+L
  return
  end
end subroutine world_date_and_time
