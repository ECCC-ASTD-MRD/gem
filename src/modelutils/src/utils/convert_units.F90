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

!/@*
module convert_units_mod
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use clib_itf_mod, only: clib_tolower
   use tdpack_const, only: GRAV, KNAMS, TCDK
   implicit none
   !@author Stephane Chamberland, 2012-05
   !@Objective Var/Fields units converter
   private
   public :: convert_units_add,convert_units_get,convert_units
   !@Description
   !  data = postopr((preopr(data) + add0) * mul + add1)
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   interface convert_units
      module procedure convert_units_r4_2d
      module procedure convert_units_r4_3d
   end interface

   type :: convert_units_T
      character(len=32) :: in_S,out_S
      character(len=32) :: preopr_S,postopr_S
      real :: mul, add0, add1
   end type convert_units_T

   integer,parameter :: MAXUNITS = 128
   character(len=*),parameter  :: KNOWN_OPR(4) = (/ &
        'exp ','ln  ','sqrt','^2  '&
        /)
   integer,save :: m_nconv = 0
   type(convert_units_T) :: m_conv(MAXUNITS)

contains

   !/@*
   subroutine convert_units_init()
      implicit none
      !*@/
      real,parameter :: MB2PA =  100.
      real,parameter :: DAM2M =  10.
      real,parameter :: CM2M =  0.01
      logical,save :: is_init_L = .false.
      integer :: istat
      real :: deg2rad
      !---------------------------------------------------------------
      if (is_init_L) return
      is_init_L = .true.

      deg2rad = acos(-1.)/180.
      !# units_in_S, units_out_S, mul, add0, add1, preopr_S, postopr_S
      !# data =  postopr((preopr(data) + add0)*mul + add1)
      istat = convert_units_add(&
           'mb','pa',MB2PA,0.,0.,' ',' '&
           )
      istat = convert_units_add(&
           'hpa','pa',MB2PA,0.,0.,' ',' '&
           )
      istat = convert_units_add(&
           'dam','m',DAM2M,0.,0.,' ',' '&
           )
      istat = convert_units_add(&
           'dam','m2/s2',DAM2M*GRAV,0.,0.,' ',' '&
           )
      istat = convert_units_add(&
           'c','k',1.,TCDK,0.,' ',' '&
           )
      istat = convert_units_add(&
           'kt','m/s2',KNAMS,0.,0.,' ',' '&
           )
      istat = convert_units_add(&
           'deg','rad',deg2rad,0.,0.,' ',' '&
           )
      istat = convert_units_add(&
           'cm','m',CM2M,0.,0.,' ',' '&
           )
      !---------------------------------------------------------------
      return
   end subroutine convert_units_init


   !/@*
   function convert_units_add(F_units_in_S,F_units_out_S,F_mul,F_add0,F_add1,F_preopr_S,F_postopr_S) result(F_istat)
      implicit none
      !@arguments
      character(len=*),intent(in) :: F_units_in_S,F_units_out_S
      character(len=*),intent(in) :: F_preopr_S,F_postopr_S
      real,intent(in) :: F_mul,F_add0,F_add1
      !@return
      integer :: F_istat
      !@Description
      !  data = postopr((preopr(data) + add0) * mul + add1)
      !*@/
      character(len=256) :: msg_S
      character(len=32) :: units_in_S,units_out_S,preopr_S,postopr_S
      integer :: istat
      !---------------------------------------------------------------
      F_istat = RMN_ERR
      if (m_nconv >= MAXUNITS-1) then
         call msg(MSG_WARNING,'(convert_units) Cannot add new units, reached max')
         return
      endif
      if (F_units_in_S == ' ' .or. F_units_out_S == ' ') then
         call msg(MSG_WARNING,'(convert_units) Cannot add units, invalid unit names:'//trim(F_units_in_S)//' : '//trim(F_units_out_S))
         return
      endif
      if (abs(F_mul) < tiny(1.)) then
         call msg(MSG_WARNING,'(convert_units) Cannot add, mul factor must be != 0')
         return
      endif
      preopr_S = F_preopr_S
      postopr_S = F_postopr_S
      istat = clib_tolower(preopr_S)
      istat = clib_tolower(postopr_S)
      if (.not.(&
           (any(preopr_S == KNOWN_OPR)  .or. preopr_S == ' ') .and. &
           (any(postopr_S == KNOWN_OPR) .or. postopr_S == ' ') &
           )) then
         call msg(MSG_WARNING,'(convert_units) Cannot add, unknown pre/post operators:'//trim(preopr_S)//' : '//trim(postopr_S))
         return
      endif
      F_istat = RMN_OK
      units_in_S = F_units_in_S
      units_out_S = F_units_out_S
      istat = clib_tolower(units_in_S)
      istat = clib_tolower(units_out_S)
      m_nconv = m_nconv + 1
      m_conv(m_nconv)%in_S = units_in_S
      m_conv(m_nconv)%out_S = units_out_S
      m_conv(m_nconv)%preopr_S = preopr_S
      m_conv(m_nconv)%postopr_S = postopr_S
      m_conv(m_nconv)%mul = F_mul
      m_conv(m_nconv)%add0 = F_add0
      m_conv(m_nconv)%add1 = F_add1
      write(msg_S,'(a,3(1pe13.6,a))') '(convert_units) Add: ['// &
           units_in_S(1:8)//'] ==> ['//units_out_S(1:8)//'] data = '// &
           trim(postopr_S)//'('//trim(preopr_S)//'(data) + ',F_add0, &
           ') * ',F_mul,' + ',F_add1,')'
      call msg(MSG_INFOPLUS,msg_S)
      !---------------------------------------------------------------
      return
   end function convert_units_add


   !/@*
   function convert_units_get(F_units_in_S,F_units_out_S,F_mul,F_add0,F_add1,F_preopr_S,F_postopr_S) result(F_istat)
      implicit none
      !@arguments
      character(len=*),intent(in) :: F_units_in_S,F_units_out_S
      character(len=*),intent(out) :: F_preopr_S,F_postopr_S
      real,intent(out) :: F_mul,F_add0,F_add1
      !@return
      integer :: F_istat
      !@Description
      !  data = postopr((preopr(data) + add0) * mul + add1)
      !*@/
      character(len=32) :: units_in_S,units_out_S
      integer :: istat,idx,n
      logical :: reverse_L
      !---------------------------------------------------------------
      F_istat = RMN_ERR
      F_mul = 1.
      F_add0 = 0.
      F_add1 = 0.
      F_preopr_S = ' '
      F_postopr_S = ' '
      if (F_units_in_S == ' ' .or. F_units_out_S == ' ') return
      call convert_units_init()
      units_in_S = F_units_in_S
      units_out_S = F_units_out_S
      istat = clib_tolower(units_in_S)
      istat = clib_tolower(units_out_S)
      idx = -1
      reverse_L = .false.
      do n=1,m_nconv
         if (m_conv(n)%in_S == units_in_S .and. &
              m_conv(n)%out_S == units_out_S) then
            idx = n
            exit
         endif
         if (m_conv(n)%in_S == units_out_S .and. &
              m_conv(n)%out_S == units_in_S) then
            reverse_L = .true.
            idx = n
            exit
         endif
      enddo
      if (idx>0) then
         F_istat = RMN_OK
         if (reverse_L) then
            F_mul = 1./min(m_conv(n)%mul,huge(1.))
            F_add0 = -m_conv(n)%add1
            F_add1 = -m_conv(n)%add0
            select case (F_preopr_S)
            case('exp')
               F_preopr_S = 'ln'
            case('ln')
               F_preopr_S = 'exp'
            case('sqrt')
               F_preopr_S = '^2'
            case('^2')
               F_preopr_S = 'sqrt'
            end select
            select case (F_postopr_S)
            case('exp')
               F_postopr_S = 'ln'
            case('ln')
               F_postopr_S = 'exp'
            case('sqrt')
               F_postopr_S = '^2'
            case('^2')
               F_postopr_S = 'sqrt'
            end select
         else
            F_mul = m_conv(n)%mul
            F_add0 = m_conv(n)%add0
            F_add1 = m_conv(n)%add1
            F_preopr_S = m_conv(n)%preopr_S
            F_postopr_S = m_conv(n)%postopr_S
         endif
      endif
      !---------------------------------------------------------------
      return
   end function convert_units_get


   !/@*
   function convert_units_r4_2d(F_data,F_units_in_S,F_units_out_S) result(F_istat)
      implicit none
      !@arguments
      character(len=*),intent(in) :: F_units_in_S,F_units_out_S
      real,intent(inout) :: F_data(:,:)
      !@return
      integer :: F_istat
      !*@/
      real :: mul, add0, add1
      character(len=32) :: preopr_S,postopr_S
      !---------------------------------------------------------------
      F_istat = convert_units_get(F_units_in_S,F_units_out_S,mul,add0,add1,preopr_S,postopr_S)
      if (RMN_IS_OK(F_istat)) then
         call msg(MSG_INFOPLUS,'(convert_units) From: ['//trim(F_units_in_S)//'] To: ['//trim(F_units_out_S)//']')
         select case (preopr_S)
         case('exp')
            F_data = exp(F_data)
         case('ln')
            F_data = log(F_data)
         case('sqrt')
            F_data = sqrt(F_data)
         case('^2')
            F_data = F_data**2
         end select

         F_data = (F_data + add0) * mul + add1

         select case (postopr_S)
         case('exp')
            F_data = exp(F_data)
         case('ln')
            F_data = log(F_data)
         case('sqrt')
            F_data = sqrt(F_data)
         case('^2')
            F_data = F_data**2
         end select
      endif
      !---------------------------------------------------------------
      return
   end function convert_units_r4_2d


   !/@*
   function convert_units_r4_3d(F_data,F_units_in_S,F_units_out_S) result(F_istat)
      implicit none
      !@arguments
      character(len=*),intent(in) :: F_units_in_S,F_units_out_S
      real,intent(inout) :: F_data(:,:,:)
      !@return
      integer :: F_istat
      !*@/
      real :: mul, add0, add1
      character(len=32) :: preopr_S,postopr_S
      !---------------------------------------------------------------
      F_istat = convert_units_get(F_units_in_S,F_units_out_S,mul,add0,add1,preopr_S,postopr_S)
      if (RMN_IS_OK(F_istat)) then
         call msg(MSG_INFOPLUS,'(convert_units) From: ['//trim(F_units_in_S)//'] To: ['//trim(F_units_out_S)//']')
         select case (preopr_S)
         case('exp')
            F_data = exp(F_data)
         case('ln')
            F_data = log(F_data)
         case('sqrt')
            F_data = sqrt(F_data)
         case('^2')
            F_data = F_data**2
         end select

         F_data = (F_data + add0) * mul + add1

         select case (postopr_S)
         case('exp')
            F_data = exp(F_data)
         case('ln')
            F_data = log(F_data)
         case('sqrt')
            F_data = sqrt(F_data)
         case('^2')
            F_data = F_data**2
         end select
      endif
      !---------------------------------------------------------------
      return
   end function convert_units_r4_3d


end module convert_units_mod
