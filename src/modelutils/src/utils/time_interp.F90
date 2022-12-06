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
#include <rmn/msg.h>

!/@*
module time_interp_mod
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use clib_itf_mod, only: clib_tolower
   use wb_itf_mod
   use mu_jdate_mod
   use rmn_gmm
   implicit none
   private
   !@objective 
   !@author Stephane Chamberland, 2011-04
   !@revision
   !  2012-04, S.Chamberland: TIME_INTERP_NEXT
   !  2014-07, S.Chamberland: re-code
   !@description
   ! Public functions
   public :: time_interp_set, time_interp_get, time_interp_typecode, &
        time_interp_dates, time_interp_retrieve, time_interp_weight, &
        time_interp_status, time_interp_ptr
   ! Public constants
   integer,parameter,public :: TIME_INTERP_NEAREST = 1
   integer,parameter,public :: TIME_INTERP_NEAR = TIME_INTERP_NEAREST
   integer,parameter,public :: TIME_INTERP_LINEAR  = 2
   integer,parameter,public :: TIME_INTERP_LINE  = TIME_INTERP_LINEAR
   integer,parameter,public :: TIME_INTERP_STEP = 3
   integer,parameter,public :: TIME_INTERP_PREV = TIME_INTERP_STEP
   integer,parameter,public :: TIME_INTERP_PREVIOUS = TIME_INTERP_STEP
   integer,parameter,public :: TIME_INTERP_NEXT = 4
   integer,parameter,public :: TIME_INTERP_INCREMENT = 5
   integer,parameter,public :: TIME_INTERP_INCR = TIME_INTERP_INCREMENT
   integer,parameter,public :: TIME_INTERP_ANY = 6
   integer,parameter,public :: TIME_INTERP_NONE = 7
   integer,parameter,public :: TIME_INTERP_NTYPES = 7
   character(len=9),parameter,public :: &
        TIME_INTERP_TYPELIST_S(TIME_INTERP_NTYPES) = &
        (/'nearest  ', 'linear   ', 'step     ', 'next     ', &
          'increment', 'any      ', 'none     '/)

   integer,parameter,public :: TIME_INTERP_FOUND_PREV = 97
   integer,parameter,public :: TIME_INTERP_FOUND_NEXT = 98
   integer,parameter,public :: TIME_INTERP_NOT_FOUND = -99
   integer,parameter,public :: TIME_INTERP_NEED_NEXT = -98
   integer,parameter,public :: TIME_INTERP_NEED_PREV = -97

   real,   parameter,public :: TIME_INTERP_WEIGHT_FACT = 100000.

   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   interface time_interp_dates
      module procedure time_interp_dates_4
      module procedure time_interp_dates_8
   end interface time_interp_dates

   interface time_interp_status
      module procedure time_interp_status_4
      module procedure time_interp_status_8
   end interface time_interp_status

   interface time_interp_set
      module procedure time_interp_set0_4
      module procedure time_interp_set1_4
      module procedure time_interp_set0_8
      module procedure time_interp_set1_8
   end interface

   interface time_interp_get
      module procedure time_interp_get0_4
      module procedure time_interp_get0_8
      module procedure time_interp_ptr_4
      module procedure time_interp_ptr_8
      module procedure time_interp_retrieve_4
      module procedure time_interp_retrieve_8
   end interface

   interface time_interp_weight
      module procedure time_interp_weight_lin_4
      module procedure time_interp_weight_lin_8
      module procedure time_interp_weight_lin1_4
      module procedure time_interp_weight_lin1_8
      module procedure time_interp_weight_type_4
      module procedure time_interp_weight_type_8
      module procedure time_interp_weight_type1_4
      module procedure time_interp_weight_type1_8
   end interface

   interface time_interp_retrieve
      module procedure time_interp_retrieve_4
      module procedure time_interp_retrieve_8
   end interface

   interface time_interp_ptr
      module procedure time_interp_ptr_4
      module procedure time_interp_ptr_8
   end interface

   character(len=5),parameter :: TIPREFIX = '_TI_/'
   character(len=2),parameter :: TI0 = ':0'
   character(len=2),parameter :: TI1 = ':1'
   real(REAL64),parameter :: EPSILON_8 = 1.d-5
   real(REAL64),parameter :: SEC_PER_HOURS_8 = 3600.0D0
   real,parameter :: EPSILON = 1.e-5
   integer, parameter :: NVLIST = 3

contains

   !TODO-later: add interface for other array type and rank
   !TODO-later: should we prevent time_interp_set within present time range?
   !TODO-later: should we allow more than 2 time frame to be saved?
   !TODO-later: should we allow time extrapolation?
   !TODO-later: should we restrict time update to one dir at a time?
   !TODO-later: should we provide minxy,maxxy,haloxy? or meta

   !/@*
   function time_interp_typecode(F_type_S) result(F_type)
      implicit none
      character(len=*),intent(in) :: F_type_S
      integer :: F_type
      !*@/
      integer :: istat,n
      character(len=4) :: type_S
      !----------------------------------------------------------------------
      F_type = RMN_ERR
      type_S = F_type_S(1:4)
      istat = clib_tolower(type_S)
      do n=1,TIME_INTERP_NTYPES
         if (type_S(1:4) == TIME_INTERP_TYPELIST_S(n)(1:4)) then
            F_type = n
            return
         endif
      enddo
      !----------------------------------------------------------------------
      return
   end function time_interp_typecode


   !/@*
   function time_interp_dates_4(F_varname_S,F_datev1,F_datev2) result(F_istat)
      character(len=*),intent(in) :: F_varname_S
      integer,intent(out) :: F_datev1,F_datev2
      integer :: F_istat
      !*@/
      integer(INT64) :: jdatev1,jdatev2
      !----------------------------------------------------------------------
      F_istat = time_interp_dates_8(F_varname_S,jdatev1,jdatev2)
      F_datev1 = jdate_to_cmc(jdatev1)
      F_datev2 = jdate_to_cmc(jdatev2)
      !----------------------------------------------------------------------
      return
   end function time_interp_dates_4


   !/@*
   function time_interp_dates_8(F_varname_S,F_datev1,F_datev2) result(F_istat)
      character(len=*),intent(in) :: F_varname_S
      integer(INT64),intent(out) :: F_datev1,F_datev2
      integer :: F_istat
      !*@/
      character(len=GMM_MAXNAMELENGTH) :: varname_S
      integer :: istat1,istat2
      type(gmm_metadata) :: mymeta
      real,pointer :: dummy3d(:,:,:)
      !----------------------------------------------------------------------
      F_istat   = TIME_INTERP_NOT_FOUND
      F_datev1  = RMN_ERR
      varname_S = priv_varname(F_varname_S,TIME_INTERP_PREV)
      !#Note: use gmm_get, gmm_getmeta is too verbose
      nullify(dummy3d)
      istat1    = gmm_get(varname_S,dummy3d,mymeta)
      if (RMN_IS_OK(istat1) .and. associated(dummy3d)) then
         F_datev1 = priv_meta2date(mymeta)
         F_istat  = TIME_INTERP_FOUND_PREV
      endif
      F_datev2  = F_datev1
      varname_S = priv_varname(F_varname_S,TIME_INTERP_NEXT)
      !#Note: use gmm_get, gmm_getmeta is too verbose
      nullify(dummy3d)
      istat2    = gmm_get(varname_S,dummy3d,mymeta)
      if (RMN_IS_OK(istat2) .and. associated(dummy3d)) then
         F_datev2  = priv_meta2date(mymeta)
         if (F_istat == TIME_INTERP_FOUND_PREV) then
            F_istat = RMN_OK
         else
            F_istat = TIME_INTERP_FOUND_NEXT
         endif
      endif
      if (F_datev1 == RMN_ERR) F_datev1 = F_datev2
      !----------------------------------------------------------------------
      return
   end function time_interp_dates_8


   !/@*
   function time_interp_weight_lin_4(F_datev,F_datev1,F_datev2) result(F_weight)
      implicit none
      integer,intent(in) :: F_datev,F_datev1,F_datev2
      real :: F_weight
      !*@/
      integer(INT64) :: jdatev,jdatev1,jdatev2
      !----------------------------------------------------------------------
      jdatev  = jdate_from_cmc(F_datev)
      jdatev1 = jdate_from_cmc(F_datev1)
      jdatev2 = jdate_from_cmc(F_datev2)
      F_weight = time_interp_weight_lin_8(jdatev,jdatev1,jdatev2)
      !----------------------------------------------------------------------
      return
   end function time_interp_weight_lin_4


   !/@*
   function time_interp_weight_lin_8(F_datev,F_datev1,F_datev2) result(F_weight)
      implicit none
      integer(INT64),intent(in) :: F_datev,F_datev1,F_datev2
      real :: F_weight
      !*@/
      integer(INT64) :: datev1,datev2
      real(REAL64) :: nhours1_8,nhours12_8,nhours2_8
      !----------------------------------------------------------------------
      F_weight = 0.
      datev1 = F_datev1
      datev2 = F_datev2
!!$      if (datev2 < 0) datev2 = datev1
!!$      if (datev1 < 0) datev1 = datev2
!!$      if (datev2 == RMN_ERR) datev2 = datev1
!!$      if (datev1 == RMN_ERR) datev1 = datev2
      if (datev2 == RMN_ERR .or. datev2 == 0) datev2 = datev1
      if (datev1 == RMN_ERR .or. datev1 == 0) datev1 = datev2
      if (datev1 == 0. .and. datev2 == 0.) return

      !TODO: check this
      nhours1_8  = dble(F_datev-datev1)/SEC_PER_HOURS_8
      nhours2_8  = dble(F_datev-datev2)/SEC_PER_HOURS_8
      nhours12_8 = dble(datev2-datev1)/SEC_PER_HOURS_8
      if (nhours1_8+EPSILON_8 < 0.) then
         if (abs(nhours12_8) < EPSILON_8) then
            F_weight = -1.
         else
            F_weight = nhours1_8/nhours12_8
         endif
      elseif (abs(nhours1_8) < EPSILON_8) then
         F_weight = 0.
!!$         if (abs(nhours12_8) < EPSILON_8) then
!!$            F_weight = 1.
!!$         endif
      elseif (abs(nhours2_8) < EPSILON_8) then
         F_weight = 1.

      elseif ((nhours2_8+EPSILON_8) < nhours12_8) then
         F_weight = (nhours12_8+nhours2_8)/nhours12_8
      else
         if (abs(nhours12_8) < EPSILON_8) then
            F_weight = 2.
         else
            F_weight = (nhours12_8+nhours2_8)/nhours12_8
        endif
      endif
      !----------------------------------------------------------------------
      return
   end function time_interp_weight_lin_8


   !/@*
   function time_interp_weight_lin1_4(F_varname_S,F_datev,F_istat) result(F_weight)
      implicit none
      character(len=*),intent(in) :: F_varname_S
      integer,intent(in) :: F_datev
      integer,intent(out) :: F_istat
      real :: F_weight
      !*@/
      integer(INT64) :: jdatev
      !----------------------------------------------------------------------
      jdatev   = jdate_from_cmc(F_datev)
      F_weight = time_interp_weight_lin1_8(F_varname_S,jdatev,F_istat)
      !----------------------------------------------------------------------
      return
   end function time_interp_weight_lin1_4


   !/@*
   function time_interp_weight_lin1_8(F_varname_S,F_datev,F_istat) result(F_weight)
      implicit none
      character(len=*),intent(in) :: F_varname_S
      integer(INT64),intent(in) :: F_datev
      integer,intent(out) :: F_istat
      real :: F_weight
      !*@/
      integer(INT64) :: datev1,datev2
      !----------------------------------------------------------------------
      F_weight = 0.
      F_istat = time_interp_dates(F_varname_S,datev1,datev2)
      if (datev1 == RMN_ERR) return
      F_weight = time_interp_weight_lin_8(F_datev,datev1,datev2)
      !----------------------------------------------------------------------
      return
   end function time_interp_weight_lin1_8


   !/@*
   function time_interp_weight_type_4(F_type,F_datev,F_datev1,F_datev2) result(F_weight)
      implicit none
      integer,intent(in) :: F_type,F_datev,F_datev1,F_datev2
      real :: F_weight
      !*@/
      integer(INT64) :: jdatev,jdatev1,jdatev2
      !----------------------------------------------------------------------
      jdatev  = jdate_from_cmc(F_datev)
      jdatev1 = jdate_from_cmc(F_datev1)
      jdatev2 = jdate_from_cmc(F_datev2)
      F_weight = time_interp_weight_type_8(F_type,jdatev,jdatev1,jdatev2)
      !----------------------------------------------------------------------
      return
   end function time_interp_weight_type_4


   !/@*
   function time_interp_weight_type_8(F_type,F_datev,F_datev1,F_datev2) result(F_weight)
      implicit none
      integer,intent(in) :: F_type
      integer(INT64),intent(in) :: F_datev,F_datev1,F_datev2
      real :: F_weight
      !*@/
      integer :: itype
      !----------------------------------------------------------------------
      F_weight = time_interp_weight_lin_8(F_datev,F_datev1,F_datev2)
      itype = max(1,min(F_type,TIME_INTERP_NTYPES))
      if (F_type /= itype .or. F_type == TIME_INTERP_LINE) return

      select case(itype)
      case(TIME_INTERP_NEAREST)
         if (F_weight+EPSILON >= 0. .and. F_weight-EPSILON <= 1.) then
            F_weight = real(nint(F_weight))
!!$         else
!!$            if (F_weight < 0.) then
!!$               F_weight = min(real(floor(F_weight)),-1.)
!!$            else
!!$               F_weight = max(real(nint(F_weight+0.5)),2.)
!!$            endif
         endif
      case(TIME_INTERP_STEP)
         if (F_weight+EPSILON >= 0. .and. F_weight-EPSILON <= 1.) then
            F_weight = 0.
            if (F_weight >= 1.) F_weight = 1.
!!$         else
         endif
      case(TIME_INTERP_NEXT)
         if (F_weight <= 0.) then
            F_weight = 0.
         else if (F_weight+EPSILON >= 0. .and. F_weight-EPSILON <= 1.) then
            F_weight = 1.
!!$         else
         endif
      end select
      !----------------------------------------------------------------------
      return
   end function time_interp_weight_type_8


   !/@*
   function time_interp_weight_type1_4(F_type,F_varname_S,F_datev,F_istat) result(F_weight)
      implicit none
      character(len=*),intent(in) :: F_varname_S
      integer,intent(in) :: F_type,F_datev
      integer,intent(out) :: F_istat
      real :: F_weight
      !*@/
      integer(INT64) :: jdatev
      !----------------------------------------------------------------------
      jdatev   = jdate_from_cmc(F_datev)
      F_weight = time_interp_weight_type1_8(F_type,F_varname_S,jdatev,F_istat)
      !----------------------------------------------------------------------
      return
   end function time_interp_weight_type1_4


   !/@*
   function time_interp_weight_type1_8(F_type,F_varname_S,F_datev,F_istat) result(F_weight)
      implicit none
      character(len=*),intent(in) :: F_varname_S
      integer,intent(in) :: F_type
      integer(INT64),intent(in) :: F_datev
      integer,intent(out) :: F_istat
      real :: F_weight
      !*@/
      integer(INT64) :: datev1,datev2
      !----------------------------------------------------------------------
      F_weight = 0.
      F_istat = time_interp_dates(F_varname_S,datev1,datev2)
      if (datev1 == RMN_ERR) return
      F_weight = time_interp_weight_type_8(F_type,F_datev,datev1,datev2)
      !----------------------------------------------------------------------
      return
   end function time_interp_weight_type1_8


   !/@*
   function time_interp_status_4(F_varname_S,F_datev,F_type) result(F_istat)
      implicit none
      character(len=*),intent(in) :: F_varname_S
      integer,intent(in) :: F_datev
      integer,intent(in),optional :: F_type
      integer :: F_istat
      !*@/
      integer(INT64) :: jdatev
      !----------------------------------------------------------------------
      jdatev  = jdate_from_cmc(F_datev)
      if (present(F_type)) then
         F_istat = time_interp_status_8(F_varname_S,jdatev,F_type)
      else
         F_istat = time_interp_status_8(F_varname_S,jdatev)
      endif
      !----------------------------------------------------------------------
      return
   end function time_interp_status_4


   !/@*
   function time_interp_status_8(F_varname_S,F_datev,F_type) result(F_istat)
      implicit none
      character(len=*),intent(in) :: F_varname_S
      integer(INT64),intent(in) :: F_datev
      integer,intent(in),optional :: F_type
      integer :: F_istat
      !*@/
      character(len=MU_JDATE_PDF_LEN) datev_S,itype_S
      character(len=512) :: msg_S
      integer :: istat,itype
      integer(INT64) :: datev1,datev2
      real :: weight
      !----------------------------------------------------------------------
      F_istat = RMN_ERR
      itype = TIME_INTERP_LINE
      if (present(F_type)) itype = max(1,min(F_type,TIME_INTERP_NTYPES))
      itype_S = TIME_INTERP_TYPELIST_S(itype)
      datev_S = jdate_to_print(F_datev)

      istat = time_interp_dates(F_varname_S,datev1,datev2)
      weight = 0.
      if (datev1 /= RMN_ERR) &
           weight = time_interp_weight_type_8(itype,F_datev,datev1,datev2)

!!$      write(msg_S,'(a,i,f9.6,a)') ' Status0/weight0: ',istat,weight,' '// &
!!$      trim(F_varname_S)//' at '//trim(datev_S)//' [type='//trim(itype_S)//']'
      write(msg_S, '(a,1x,i0,1x,g0,1x,a)') 'Status0/weight0: ',istat,weight,' '//trim(F_varname_S)// &
           ' at '//trim(datev_S)//' [type='//trim(itype_S)//']'
      call msg(MSG_INFOPLUS,'(time_interp)'//trim(msg_S))

      select case(istat)
      case(RMN_OK)
         if (weight+EPSILON < 0.) then
            F_istat = TIME_INTERP_NEED_PREV
!!$            if (itype == TIME_INTERP_NEXT) F_istat = 0
         elseif (weight-EPSILON > 1.) then
            F_istat = TIME_INTERP_NEED_NEXT
!!$            if (itype == TIME_INTERP_STEP) F_istat = TIME_INTERP_WEIGHT_FACT
         elseif (abs(weight) < EPSILON) then
            F_istat = 0
         elseif (abs(weight-1.) < EPSILON) then
            F_istat = TIME_INTERP_WEIGHT_FACT
         else
            F_istat = nint(weight*TIME_INTERP_WEIGHT_FACT)
         endif
      case(TIME_INTERP_FOUND_PREV)
         F_istat = TIME_INTERP_NEED_NEXT
         if (abs(weight) < EPSILON .and. &
!!$              itype /= TIME_INTERP_NEXT .and. itype /= TIME_INTERP_INCR) &
              itype /= TIME_INTERP_INCR) &
              F_istat = 0
      case(TIME_INTERP_FOUND_NEXT)
         F_istat = TIME_INTERP_NEED_PREV
         if (abs(weight-1) < EPSILON .and. &
               itype /= TIME_INTERP_STEP .and. itype /= TIME_INTERP_INCR) &
             F_istat = TIME_INTERP_WEIGHT_FACT
      case default
         F_istat = TIME_INTERP_NOT_FOUND
         if (itype == TIME_INTERP_NEXT) F_istat = TIME_INTERP_NEED_NEXT
     end select
 
!!$      write(msg_S,'(a,i,f9.6,a)') ' Status/weight: ',F_istat,weight,' '// &
!!$     trim(F_varname_S)//' at '//trim(datev_S)//' [type='//trim(itype_S)//']'
      write(msg_S, '(a,1x,i0,1x,g0,1x,a)') 'Status/weight: ',F_istat,weight,' '//trim(F_varname_S)// &
           ' at '//trim(datev_S)//' [type='//trim(itype_S)//']'
      call msg(MSG_INFOPLUS,'(time_interp)'//trim(msg_S))
      !----------------------------------------------------------------------
      return
   end function time_interp_status_8


   !/@*
   function time_interp_set0_4(F_data,F_varname_S,F_datev,F_in_restart_L,F_flag) result(F_istat)
      implicit none
      real,pointer :: F_data(:,:,:)
      character(len=*),intent(in) :: F_varname_S
      integer,intent(in) :: F_datev
      logical,intent(in),optional :: F_in_restart_L
      integer,intent(in),optional :: F_flag
      integer :: F_istat
      !*@/
      integer :: myflag
      integer(INT64) :: jdatev
      !----------------------------------------------------------------------
      jdatev  = jdate_from_cmc(F_datev)
      myflag = 0
      if (present(F_flag)) myflag = F_flag
      if (present(F_in_restart_L)) then
         F_istat = time_interp_set0_8(F_data,F_varname_S,jdatev,F_in_restart_L,myflag)
      else
         F_istat = time_interp_set0_8(F_data,F_varname_S,jdatev,F_flag=myflag)
      endif
      !----------------------------------------------------------------------
      return
   end function time_interp_set0_4


   !/@*
   function time_interp_set0_8(F_data,F_varname_S,F_datev,F_in_restart_L,F_flag) result(F_istat)
      implicit none
      real,pointer :: F_data(:,:,:)
      character(len=*),intent(in) :: F_varname_S
      integer(INT64),intent(in) :: F_datev
      logical,intent(in),optional :: F_in_restart_L
      integer,intent(in),optional :: F_flag
      integer :: F_istat
      !*@/
      integer :: myflag
      !----------------------------------------------------------------------
      myflag = 0
      if (present(F_flag)) myflag = F_flag
      if (present(F_in_restart_L)) then
         F_istat = time_interp_set1_8(F_data,F_varname_S,F_datev,' ',' ',F_in_restart_L,F_flag=myflag)
      else
         F_istat = time_interp_set1_8(F_data,F_varname_S,F_datev,' ',' ',F_flag=myflag)
      endif
      !----------------------------------------------------------------------
      return
   end function time_interp_set0_8


   !/@*
   function time_interp_set1_4(F_data,F_varname_S,F_datev,F_vgrid_S,F_sfcfld_S, &
        F_in_restart_L,F_flag) result(F_istat)
      implicit none
      real,pointer :: F_data(:,:,:)
      character(len=*),intent(in) :: F_varname_S,F_vgrid_S,F_sfcfld_S
      integer,intent(in) :: F_datev
      logical,intent(in),optional :: F_in_restart_L
      integer,intent(in),optional :: F_flag
      integer :: F_istat
      !*@/
      integer :: myflag
      integer(INT64) :: jdatev
      !----------------------------------------------------------------------
      jdatev  = jdate_from_cmc(F_datev)
      myflag = 0
      if (present(F_flag)) myflag = F_flag
      if (present(F_in_restart_L)) then
         F_istat = time_interp_set1_8(F_data,F_varname_S,jdatev,F_vgrid_S,F_sfcfld_S,F_in_restart_L,myflag)
      else
         F_istat = time_interp_set1_8(F_data,F_varname_S,jdatev,F_vgrid_S,F_sfcfld_S,F_flag=myflag)
      endif
      !----------------------------------------------------------------------
      return
   end function time_interp_set1_4


   !/@*
   function time_interp_set1_8(F_data,F_varname_S,F_datev,F_vgrid_S, &
        F_sfcfld_S,F_in_restart_L,F_flag,F_sfcfld2_S) result(F_istat)
      implicit none
      real,pointer :: F_data(:,:,:)
      character(len=*),intent(in) :: F_varname_S,F_vgrid_S,F_sfcfld_S
      integer(INT64),intent(in) :: F_datev
      logical,intent(in),optional :: F_in_restart_L
      integer,intent(in),optional :: F_flag
      character(len=*),intent(in),optional :: F_sfcfld2_S
      integer :: F_istat
      !*@/
      integer :: istat,istat2,rstst_flag,myflag
      real :: weight
      character(len=MU_JDATE_PDF_LEN) datev_S
      character(len=512) :: msg_S, sfcfld2_S
      !----------------------------------------------------------------------
      F_istat = RMN_ERR
      if (trim(F_varname_S)==' ' .or. .not.associated(F_data)) return
      rstst_flag = 0
      if (present(F_in_restart_L)) then
         if (F_in_restart_L) rstst_flag = GMM_FLAG_RSTR
      endif
      myflag = 0
      if (present(F_flag)) myflag = F_flag
      sfcfld2_S = ' '
      if (present(F_sfcfld2_S)) sfcfld2_S = F_sfcfld2_S

      datev_S = jdate_to_print(F_datev)
      msg_S = ' '
!!$      write(msg_S,'(a,2e14.6,a)') trim(F_varname_S)//' valid at '//trim(datev_S)// &
!!$           [minmax=',minval(F_data),maxval(F_data),'] vgd='//trim(F_vgrid_S)// &
!!$           '; rfld='//trim(F_sfcfld_S)//' '//trim(sfcfld2_S)
      write(msg_S,'(a)') trim(F_varname_S)//' valid at '//trim(datev_S)//' vgd='// &
           trim(F_vgrid_S)//'; rfld='//trim(F_sfcfld_S)//' '//trim(sfcfld2_S)

      weight = time_interp_weight(F_varname_S,F_datev,istat)
      select case(istat)
      case(RMN_OK)
         if (weight-EPSILON > 1.) then
            istat2 = priv_shuffle(F_varname_S)
            F_istat = priv_setdata(TIME_INTERP_NEXT,F_varname_S,F_datev,F_vgrid_S, &
                 F_sfcfld_S,rstst_flag,myflag,F_data, sfcfld2_S)
         else
            if (weight < 0.) then
               call msg(MSG_INFOPLUS,'(time_interp) Backward set not yet implemented: '//trim(msg_S))
               !#TODO: allow backward set
            else
               call msg(MSG_INFOPLUS,'(time_interp) Set cannot replace a value within '// &
                    'existing time range: '//trim(msg_S))
            endif
         endif
      case(TIME_INTERP_FOUND_PREV)
         if (abs(weight) > EPSILON) then
            F_istat = priv_setdata(TIME_INTERP_NEXT,F_varname_S,F_datev,F_vgrid_S, &
                 F_sfcfld_S,rstst_flag,myflag,F_data, sfcfld2_S)
            if (weight < 0.) istat2 = priv_shuffle(F_varname_S)
         else
            call msg(MSG_INFOPLUS,'(time_interp) Set cannot replace a value within '//&
                 'existing time range: '//trim(msg_S))
         endif
      case(TIME_INTERP_FOUND_NEXT)
         continue
      case default
            F_istat = priv_setdata(TIME_INTERP_PREV,F_varname_S,F_datev,F_vgrid_S, &
                 F_sfcfld_S,rstst_flag,myflag,F_data, sfcfld2_S)
      end select

      if (RMN_IS_OK(F_istat)) then
         call msg(MSG_INFOPLUS,'(time_interp) Set: '//trim(msg_S))
      else
         call msg(MSG_INFOPLUS,'(time_interp) Set ERROR: '//trim(msg_S))
      endif
      !----------------------------------------------------------------------
      return
   end function time_interp_set1_8


   !/@*
   function time_interp_retrieve_4(F_varname_S,F_next_prev,F_datev, &
        F_meta,F_data,F_vgrid_S,F_sfcfld_S,F_flag) result(F_istat)
      implicit none
      character(len=*),intent(in) :: F_varname_S
      integer,intent(in) :: F_next_prev
      integer,intent(out) :: F_datev
      type(gmm_metadata),intent(out),optional :: F_meta
      real,pointer, optional :: F_data(:,:,:)
      character(len=*),intent(out),optional :: F_vgrid_S,F_sfcfld_S
      integer,intent(out),optional :: F_flag
      integer :: F_istat
      !*@/
      integer(INT64) :: jdatev
      type(gmm_metadata) :: meta
      real, pointer :: data(:,:,:)
      character(len=GMM_MAXNAMELENGTH) :: vgrid_S, sfcfld_S
      integer :: flag
      !----------------------------------------------------------------------
      nullify(data)
      F_istat = time_interp_retrieve(F_varname_S, F_next_prev, jdatev, &
           meta, data, vgrid_S, sfcfld_S, flag)
      F_datev = jdate_to_cmc(jdatev)
      if (present(F_meta)) F_meta = meta
      if (present(F_data)) F_data => data
      if (present(F_vgrid_S)) F_vgrid_S = vgrid_S
      if (present(F_sfcfld_S)) F_sfcfld_S = sfcfld_S
      if (present(F_flag)) F_flag = flag
      !----------------------------------------------------------------------
      return
   end function time_interp_retrieve_4


   !/@*
   function time_interp_retrieve_8(F_varname_S,F_next_prev,F_datev,F_meta,F_data, &
        F_vgrid_S,F_sfcfld_S,F_flag, F_sfcfld2_S) result(F_istat)
      implicit none
      character(len=*),intent(in) :: F_varname_S
      integer,intent(in) :: F_next_prev
      integer(INT64),intent(out) :: F_datev
      type(gmm_metadata),intent(out),optional :: F_meta
      real,pointer, optional :: F_data(:,:,:)
      character(len=*),intent(out),optional :: F_vgrid_S, F_sfcfld_S, F_sfcfld2_S
      integer,intent(out),optional :: F_flag
      integer :: F_istat
      !*@/
      character(len=512) :: msg_S,dummy_S
      character(len=GMM_MAXNAMELENGTH) :: varname_S,vlist0_S(NVLIST),datev_S
      integer :: istat,nitems
      type(gmm_metadata) :: mymeta
      real,pointer :: dummy3d(:,:,:)
      !----------------------------------------------------------------------
      F_istat = TIME_INTERP_NOT_FOUND
      F_datev = 0
      if (present(F_meta)) F_meta = GMM_NULL_METADATA
      if (present(F_vgrid_S)) F_vgrid_S = ' '
      if (present(F_sfcfld_S)) F_sfcfld_S = ' '
      if (present(F_sfcfld2_S)) F_sfcfld2_S = ' '
      if (present(F_flag)) F_flag = 0

      varname_S = priv_varname(F_varname_S,F_next_prev)
      nullify(dummy3d)
      if (present(F_data)) then
         nullify(F_data)
         istat = gmm_get(varname_S,F_data,mymeta)
         if (associated(F_data)) dummy3d => F_data
      else
         !#use gmm_get, gmm_getmeta is too verbose
         istat = gmm_get(varname_S,dummy3d,mymeta)
      endif
      if (associated(dummy3d)) F_istat = RMN_OK

      dummy_S = 'PREV'//' '//trim(F_varname_S)
      if (F_next_prev==TIME_INTERP_NEXT) &
           dummy_S = 'NEXT'//' '//trim(F_varname_S)
      if (.not.RMN_IS_OK(F_istat)) then
         msg_S  = trim(dummy_S)//' [NOT FOUND]'
      else
         F_datev = priv_meta2date(mymeta)
         datev_S = ' ' ; msg_S = ' '
         datev_S = jdate_to_print(F_datev)
!!$         write(msg_S,'(a,2e14.6,a)') trim(dummy_S)//' valid at '//trim(datev_S)// &
!!$             ' [minmax=',minval(dummy3d),maxval(dummy3d),']'
         write(msg_S, '(a)') trim(dummy_S)//' valid at '//trim(datev_S)
      endif
      call msg(MSG_INFOPLUS,'(time_interp) Retrieve: '//trim(msg_S))
      if (.not.RMN_IS_OK(F_istat)) return

      if (present(F_meta)) F_meta = mymeta
      if (present(F_vgrid_S) .or. present(F_sfcfld_S) .or. present(F_sfcfld2_S)) then
         istat = wb_get(varname_S,vlist0_S,nitems)
         if (RMN_IS_OK(istat)) then
            if (present(F_vgrid_S)) F_vgrid_S = vlist0_S(1)
            if (present(F_sfcfld_S)) F_sfcfld_S = vlist0_S(2)
            if (present(F_sfcfld2_S)) F_sfcfld2_S = vlist0_S(3)
!!$         else
!!$            F_istat = TIME_INTERP_NOT_FOUND
         endif
      endif
      if (present(F_flag)) F_flag = priv_meta2flag(mymeta)
      !----------------------------------------------------------------------
      return
   end function time_interp_retrieve_8


   !/@*
   function time_interp_ptr_4(F_data,F_data0,F_data1,F_datev,F_datev0,F_datev1, &
        F_type,F_incr_len,F_varname_S) result(F_istat)
      implicit none
      real,pointer :: F_data(:,:,:),F_data0(:,:,:),F_data1(:,:,:)
      integer,intent(in) :: F_datev,F_datev0,F_datev1
      integer,intent(in),optional :: F_type
      integer,intent(in),optional :: F_incr_len !# Lenght of increment interval [sec]
      character(len=*),intent(in),optional :: F_varname_S
      integer :: F_istat
      !*@/
      integer(INT64) :: jdatev,jdatev0,jdatev1
      integer :: itype, incr_len
      character(len=GMM_MAXNAMELENGTH) :: varname_S
      !----------------------------------------------------------------------
      itype = TIME_INTERP_LINE
      incr_len = 3600
      varname_S = ' '
      if (present(F_type)) itype = max(1,min(F_type,TIME_INTERP_NTYPES))
      if (present(F_incr_len)) incr_len = F_incr_len
      if (present(F_varname_S)) varname_S = F_varname_S
      jdatev  = jdate_from_cmc(F_datev) !#TODO: if /= RMN_ERR?
      jdatev0 = jdate_from_cmc(F_datev0)
      jdatev1 = jdate_from_cmc(F_datev1)
      F_istat = time_interp_ptr_8(F_data,F_data0,F_data1,jdatev,jdatev0,jdatev1,itype,incr_len,varname_S)
       !----------------------------------------------------------------------
      return
   end function time_interp_ptr_4


   !/@*
   function time_interp_ptr_8(F_data,F_data0,F_data1,F_datev,F_datev0,F_datev1, &
        F_type,F_incr_len,F_varname_S) result(F_istat)
      implicit none
      real,pointer :: F_data(:,:,:),F_data0(:,:,:),F_data1(:,:,:)
      integer(INT64),intent(in) :: F_datev,F_datev0,F_datev1
      integer,intent(in),optional :: F_type
      integer,intent(in),optional :: F_incr_len !# Lenght of increment interval [sec]
      character(len=*),intent(in),optional :: F_varname_S
      integer :: F_istat
      !*@/
      logical :: ok_L
      integer :: istat,itype,lijk(3),uijk(3)
      real :: weight,incr_len
      real(REAL64) :: nhours1_8
      character(len=16) datev_S,varname_S
      character(len=512) :: itype_S,msg_S
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG, '(time_interp) get ptr [BEGIN]')
      F_istat = RMN_ERR
      itype = TIME_INTERP_LINE
      incr_len = 1.
      varname_S = ' '
      if (present(F_type)) itype = max(1,min(F_type,TIME_INTERP_NTYPES))
      if (present(F_incr_len)) incr_len = real(max(1,F_incr_len))/3600.
      if (present(F_varname_S)) varname_S = F_varname_S
      itype_S = TIME_INTERP_TYPELIST_S(itype)

      weight = time_interp_weight(itype,F_datev,F_datev0,F_datev1)
      if (weight < -EPSILON .or. weight > 1.+EPSILON) then
         call msg(MSG_DEBUG,'(time_interp) ERROR - Get ptr ('//trim(itype_S)//'): Does not extrapolate')
         F_istat = TIME_INTERP_NEED_NEXT
         if (weight < -EPSILON) F_istat = TIME_INTERP_NEED_PREV
         return
      endif

      if (.not.associated(F_data)) then
         if (associated(F_data0)) then
            lijk = lbound(F_data0) ; uijk = ubound(F_data0)
            allocate(F_data(lijk(1):uijk(1),lijk(2):uijk(2),lijk(3):uijk(3)),stat=istat)
         elseif (associated(F_data1)) then
            lijk = lbound(F_data1) ; uijk = ubound(F_data1)
            allocate(F_data(lijk(1):uijk(1),lijk(2):uijk(2),lijk(3):uijk(3)),stat=istat)
         else
            call msg(MSG_DEBUG,'(time_interp) ERROR Get ptr - Only got null pointer as input')
            return
         endif
      else
         ok_L = .true.
         if (associated(F_data0)) then
            if (.not.all(shape(F_data)==shape(F_data0))) ok_L = .false.
         endif
         if (associated(F_data1)) then
            if (.not.all(shape(F_data)==shape(F_data1))) ok_L = .false.
         endif
         if (.not.ok_L) then
            call msg(MSG_DEBUG,'(time_interp) ERROR - Get ptr ('//trim(itype_S)// &
                 '):  Cannot get values, array bounds not compatible')
            return
         endif
      endif

      if (.not.associated(F_data)) then
         call msg(MSG_DEBUG,'(time_interp) ERROR Get ptr - Probleme allocating mem')
         return
      endif

      if (abs(weight) < EPSILON .and. itype /= TIME_INTERP_INCR) then
         if (.not.associated(F_data0)) then
            F_istat = TIME_INTERP_NEED_PREV
            return
         endif
         F_data = F_data0
      else if (abs(weight-1.) < EPSILON .and. itype /= TIME_INTERP_INCR) then
         if (.not.associated(F_data1)) then
            F_istat = TIME_INTERP_NEED_NEXT
            call msg(MSG_DEBUG,'(time_interp) get_ptr - NEED_NEXT')
            return
         endif
         F_data = F_data1
      else
         if (.not.(associated(F_data0).and.associated(F_data1))) then
            call msg(MSG_DEBUG,'(time_interp) ERROR Get ptr - Missing input data (null ptr)')
            return
         endif
         if (.not.all(shape(F_data0)==shape(F_data1))) then
            call msg(MSG_DEBUG,'(time_interp) ERROR - Get ('//trim(itype_S)// &
                 '):  Cannot get values, array bounds not compatible')
            return
         endif
         if (itype == TIME_INTERP_INCR) then
            !TODO: for intervals going outside datev0:datev1 need a 3rd date
            nhours1_8 = abs(dble(F_datev0 - F_datev1)/SEC_PER_HOURS_8)
            F_data = (F_data1 - F_data0) / abs(real(nhours1_8)/incr_len)
         else
            F_data = (1.-weight)*F_data0 + weight*F_data1
         endif
      endif

      if (varname_S /= ' ') then
         datev_S = jdate_to_print(F_datev)
         msg_S = ' '
         if (itype /= TIME_INTERP_INCR) then
            write(msg_S,'(a,a,2e14.6,a,f9.6,a,f9.6,a)') trim(varname_S),' valid at '// &
                 trim(datev_S)//' [minmax=',minval(F_data),maxval(F_data), &
                 '] (data=',1.-weight,'*data0 + ',weight,'*data1)'
         else
            write(msg_S, '(a,a,2e14.6,a)') trim(varname_S),' valid at '//trim(datev_S)// &
                 ' [minmax=',minval(F_data),maxval(F_data),']'
         endif
         call msg(MSG_INFOPLUS,'(time_interp) Get ('//trim(itype_S)//'): '//trim(msg_S))
      endif
      F_istat = max(RMN_OK, nint(abs(weight)*TIME_INTERP_WEIGHT_FACT))
      write(msg_S, '(i0)') F_istat
      call msg(MSG_DEBUG, '(time_interp) get ptr [END] '//trim(msg_S))
      !----------------------------------------------------------------------
      return
   end function time_interp_ptr_8


   !#TODO: below


   !/@*
   function time_interp_get0_4(F_data,F_varname_S,F_datev,F_type,F_incr_len) result(F_istat)
      implicit none
      real,pointer :: F_data(:,:,:)
      character(len=*),intent(in) :: F_varname_S
      integer,intent(in) :: F_datev
      integer,intent(in),optional :: F_type
      integer,intent(in),optional :: F_incr_len !# Lenght of increment interval [sec]
      integer :: F_istat
      !*@/
      integer(INT64) :: jdatev
      integer :: itype, incr_len
      !----------------------------------------------------------------------
      itype = TIME_INTERP_LINE
      incr_len = 3600
      if (present(F_type)) itype = max(1,min(F_type,TIME_INTERP_NTYPES))
      if (present(F_incr_len)) incr_len = F_incr_len
      jdatev  = jdate_from_cmc(F_datev) !#TODO: if /= RMN_ERR?
      F_istat = time_interp_get0_8(F_data,F_varname_S,jdatev,itype,incr_len)
      !----------------------------------------------------------------------
      return
   end function time_interp_get0_4


   !/@*
   function time_interp_get0_8(F_data,F_varname_S,F_datev,F_type,F_incr_len) result(F_istat)
      implicit none
      real,pointer :: F_data(:,:,:)
      character(len=*),intent(in) :: F_varname_S
      integer(INT64),intent(in) :: F_datev
      integer,intent(in),optional :: F_type
      integer,intent(in),optional :: F_incr_len !# Lenght of increment interval [sec]
      integer :: F_istat
      !*@/
      real,pointer :: data0(:,:,:),data1(:,:,:)
      type(gmm_metadata) :: meta
      integer(INT64) :: datev0,datev1
      integer :: istat,itype,incr_len
      real :: weight,incr_len_r
      !----------------------------------------------------------------------
      itype = TIME_INTERP_LINE
      if (present(F_type)) itype = max(1,min(F_type,TIME_INTERP_NTYPES))
      incr_len_r = 0.
      if (present(F_incr_len)) incr_len_r = real(F_incr_len)/3600.
      incr_len = 0
      if (present(F_incr_len)) incr_len = F_incr_len

      F_istat = time_interp_status(F_varname_S,F_datev,itype)
      if (.not.RMN_IS_OK(F_istat)) return
      weight = real(F_istat)/TIME_INTERP_WEIGHT_FACT

      nullify(data0,data1) ; datev0 = 0 ; datev1 = 0
      istat = time_interp_retrieve(F_varname_S,TIME_INTERP_PREV,datev0,meta,data0)
      istat = time_interp_retrieve(F_varname_S,TIME_INTERP_NEXT,datev1,meta,data1)
      !TODO: may want to swap data1,data2 depending on itype.and.weight
      if (F_istat /= TIME_INTERP_WEIGHT_FACT .and. .not.associated(data0)) F_istat = RMN_ERR
      if (F_istat /= 0 .and. .not.associated(data1)) F_istat = RMN_ERR
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_ERROR,'(time_interp) Problem retrieving saved data')
         return
      endif

      F_istat = time_interp_ptr(F_data,data0,data1,F_datev,datev0,datev1,itype,incr_len,F_varname_S)
      !----------------------------------------------------------------------
      return
   end function time_interp_get0_8


   !==== Private Functions =================================================


   function priv_varname(F_inname_S,F_next_prev) result(F_outname_S)
      implicit none
      character(len=*),intent(in) :: F_inname_S
      integer,intent(in) :: F_next_prev
      character(len=GMM_MAXNAMELENGTH) :: F_outname_S
      integer :: istat
      F_outname_S = F_inname_S
      istat = clib_tolower(F_outname_S) !# Note: GMM is NOT case independent
      if (F_next_prev == TIME_INTERP_NEXT) then
         F_outname_S = TIPREFIX//trim(F_outname_S)//TI1
      else
         F_outname_S = TIPREFIX//trim(F_outname_S)//TI0
      endif
   end function priv_varname


   function priv_meta2date(F_meta) result(F_date)
      implicit none
      type(gmm_metadata),intent(in) :: F_meta
      integer(INT64) :: F_date
      F_date = F_meta%a%uuid2
   end function priv_meta2date

   function priv_meta2flag(F_meta) result(F_flag)
      implicit none
      type(gmm_metadata),intent(in) :: F_meta
      integer :: F_flag
      F_flag = F_meta%a%uuid1
   end function priv_meta2flag


   function priv_date2meta(F_meta,F_date,F_flag) result(F_meta2)
      implicit none
      type(gmm_metadata),intent(in) :: F_meta
      integer(INT64),intent(in) :: F_date
      integer,intent(in) :: F_flag
      type(gmm_metadata) :: F_meta2
      F_meta2 = F_meta
      F_meta2%a%uuid1 = F_flag
      F_meta2%a%uuid2 = F_date
   end function priv_date2meta


   function priv_setmeta(F_data,F_date,F_flag) result(F_meta)
      implicit none
      real,pointer :: F_data(:,:,:)
      integer(INT64),intent(in) :: F_date
      integer,intent(in) :: F_flag
      type(gmm_metadata) :: F_meta
      integer :: ldims(3),udims(3),ii
      F_meta = GMM_NULL_METADATA
      F_meta%a%uuid1 = F_flag
      F_meta%a%uuid2 = F_date
      ldims = lbound(F_data)
      udims = ubound(F_data)
      do ii = 1,3
         F_meta%l(ii) = gmm_layout(ldims(ii),udims(ii),0,0,udims(ii))
      enddo
   end function priv_setmeta


   function priv_shuffle(F_varname_S) result(F_istat)
      implicit none
      character(len=*),intent(in) :: F_varname_S
      character(len=GMM_MAXNAMELENGTH) :: list_S(2),vlist1_S(NVLIST),vlist2_S(NVLIST)
      integer :: F_istat,istat,nitems
      list_S(1) = priv_varname(F_varname_S,TIME_INTERP_PREV)
      list_S(2) = priv_varname(F_varname_S,TIME_INTERP_NEXT)
      istat = wb_get(list_S(1),vlist1_S,nitems)
      istat = wb_get(list_S(2),vlist2_S,nitems)
      F_istat = wb_put(list_S(1),vlist2_S,WB_REWRITE_MANY)
      F_istat = min(wb_put(list_S(2),vlist1_S,WB_REWRITE_MANY),F_istat)
      istat = gmm_shuffle(list_S)
   end function priv_shuffle


   !/@*
   function priv_setdata(F_next_prev,F_varname_S,F_datev,F_vgrid_S,F_sfcfld_S, &
        F_rstst_flag,F_myflag,F_data, F_sfcfld2_S) result(F_istat)
      implicit none
      real,pointer :: F_data(:,:,:)
      character(len=*),intent(in) :: F_varname_S,F_vgrid_S,F_sfcfld_S, F_sfcfld2_S
      integer,intent(in) :: F_next_prev
      integer(INT64),intent(in) :: F_datev
      integer,intent(in) :: F_rstst_flag,F_myflag
      integer :: F_istat
      !*@/
      type(gmm_metadata) :: meta0
      character(len=GMM_MAXNAMELENGTH) :: varname_S,varname2_S,vlist_S(NVLIST)
      logical :: ok_L
      integer :: istat,ldims(3),udims(3),ii
      real,pointer :: data(:,:,:)
      !----------------------------------------------------------------------
      F_istat = RMN_ERR
      if (.not.associated(F_data)) return
      varname_S = priv_varname(F_varname_S,F_next_prev)

      nullify(data)
      istat = gmm_get(varname_S,data,meta0)
      if (.not.(RMN_IS_OK(istat) .and. associated(data))) then
         if (F_next_prev == TIME_INTERP_NEXT) then
            varname2_S = priv_varname(F_varname_S,TIME_INTERP_PREV)
            istat = gmm_get(varname2_S,data,meta0)
            nullify(data)
         endif
         if (.not.RMN_IS_OK(istat)) &
              meta0 = priv_setmeta(F_data,F_datev,F_myflag)
      endif

      ldims = lbound(F_data)
      udims = ubound(F_data)
      ok_L = .true.
      do ii = 1,3
         if  (ldims(ii) /= meta0%l(ii)%low .or. &
              udims(ii) /= meta0%l(ii)%high) ok_L = .false.
      end do
      if (.not.ok_L) then
         call msg(MSG_ERROR,'(time_interp) Set: Cannot set values, array bounds '// &
              'not compatible wiht previous set for: '//trim(F_varname_S))
         return
      endif
      if (.not.associated(data)) &
           istat = gmm_create(varname_S,data,meta0,F_rstst_flag)
      if (.not.associated(data)) then
         call msg(MSG_ERROR,'(time_interp) Set: Problem getting pointer for: '//trim(F_varname_S))
         return
      endif

      meta0 = priv_date2meta(meta0,F_datev,F_myflag)
      istat = gmm_updatemeta(varname_S,meta0)
      data = F_data
      vlist_S(1) = F_vgrid_S
      vlist_S(2) = F_sfcfld_S
      vlist_S(3) = F_sfcfld2_S
      istat = wb_put(varname_S,vlist_S,WB_REWRITE_MANY)
      F_istat = RMN_OK
      !----------------------------------------------------------------------
      return
   end function priv_setdata

end module time_interp_mod
