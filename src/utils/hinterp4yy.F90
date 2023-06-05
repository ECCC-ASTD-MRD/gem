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
module hinterp4yy_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_tolower
   use ezgrid_mod
   implicit none
   private
   !@objective 
   !@author Stephane Chamberland,2011-04
   !@description
   ! Public functions
   public :: hinterp4yy2d
   ! Public constants
   character(len=8),parameter,public :: HINTERP4YY_NEAREST = 'nearest'
   character(len=8),parameter,public :: HINTERP4YY_LINEAR = 'linear'
   character(len=8),parameter,public :: HINTERP4YY_CUBIC = 'cubic'
   integer,parameter,public :: HINTERP4YY_OK = 0
   integer,parameter,public :: HINTERP4YY_NONE = 1
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   interface hinterp4yy2d
      module procedure hinterp4yy2d_scalar
      module procedure hinterp4yy2d_vect
   end interface


contains

   !/@*
   function hinterp4yy2d_scalar(F_outdata,F_indata,F_k,F_ingridid,F_outgridid,F_coregridid,F_h_int_S,F_varname_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      real,dimension(:,:,:),pointer :: F_outdata,F_indata
      integer,intent(in) :: F_k,F_ingridid,F_outgridid,F_coregridid
      character(len=*),intent(in) :: F_h_int_S,F_varname_S
      !@return
      integer :: F_istat
      !*@/
      real,dimension(:,:,:),pointer :: outdata2,indata2
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(hinterp4yy2d_scalar) [BEGIN]')
      nullify(outdata2,indata2)
      F_istat = hinterp4yy2d_vect(F_outdata,outdata2,F_indata,indata2,F_k,F_ingridid, &
           F_outgridid,F_coregridid,F_h_int_S,F_varname_S,' ')
      call msg(MSG_DEBUG,'(hinterp4yy2d_scalar) [END]')
      !------------------------------------------------------------------
      return
   end function hinterp4yy2d_scalar


   !/@*
   function hinterp4yy2d_vect(F_outdata,F_outdata2,F_indata,F_indata2,F_k,F_ingridid, &
        F_outgridid,F_coregridid,F_h_int_S,F_varname_S,F_varname2_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      real,dimension(:,:,:),pointer :: F_outdata,F_outdata2,F_indata,F_indata2
      integer,intent(in) :: F_k,F_ingridid,F_outgridid,F_coregridid
      character(len=*),intent(in) :: F_h_int_S,F_varname_S,F_varname2_S
      !@return
      integer :: F_istat
      !*@/
      integer :: gridsetid,subgridid
      character(len=16) :: h_int_S, onesubgrid_S
      integer :: istat
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(hinterp4yy2d_vect) [BEGIN]')
      F_istat = RMN_ERR

      if (.not.associated(F_indata) .or. .not.associated(F_outdata)) then
         call msg(MSG_WARNING,'(hinterp4yy2d) Cannot hinterpolate, Pointer not associated for: '//trim(F_varname_S)//' '//trim(F_varname2_S))
         return
      endif

      h_int_S = F_h_int_S
      call priv_ezparams(gridsetid,subgridid,onesubgrid_S,h_int_S, &
           F_ingridid,F_outgridid,F_coregridid)

      !#TODO: rm temporary patch after samegrid fix
!!$      if (h_int_S(1:4) == 'near' .and. &
!!$           .not.any(F_h_int_S(1:1) == (/'n','N'/))) then
!!$         h_int_S = F_h_int_S
!!$!#         onesubgrid_S = 'NO'
!!$      endif

!!$      if (h_int_S(1:4) == 'near' .and. &
!!$           .not.any(F_h_int_S(1:1) == (/'n','N'/))) then
!!$         print *,'(hinterp4yy) Changed from '//trim(F_h_int_S)// &
!!$              ' to Nearest horizontal interpolation for: '//trim(F_varname_S )//' '// &
!!$              trim(F_varname2_S )
!!$      else
!!$         print *,'(hinterp4yy) Not Changed '//trim(F_h_int_S)// &
!!$              ' horizontal interpolation for: '//trim(F_varname_S )//' '//trim(F_varname2_S )
!!$      endif

!!$      print *,'hinterp4yy2d',subgridid,';',trim(onesubgrid_S),';',trim(h_int_S)

      if (subgridid >= 0) &
           istat = ezsetival('SUBGRIDID',subgridid)
!!$         !istat = ezsetopt('EXTRAP_DEGREE','ABORT')
!!$         istat = ezsetopt('EXTRAP_DEGREE','VALUE')
!!$         istat = ezsetval('EXTRAP_VALUE',huge(1.))
      istat = ezsetopt('INTERP_DEGREE',trim(h_int_S))
      istat = min(ezsetopt('USE_1SUBGRID',trim(onesubgrid_S)),istat)

      if (h_int_S(1:4) == 'near') &
           call msg(MSG_INFO,'(hinterp4yy) Nearest horizontal interpolation for: '//trim(F_varname_S )//' '//trim(F_varname2_S ))
      if (h_int_S(1:4) /= 'near') &
           call msg(MSG_INFOPLUS,'(hinterp4yy) Horizontal Interpolation for: '//trim(F_varname_S )//' '//trim(F_varname2_S ))

      if (F_varname2_S == ' ' .or. .not.associated(F_indata2) .or. .not.associated(F_outdata2)) then
         F_istat = ezsint(F_outdata(:,:,F_k),F_indata)
      else
         if (subgridid < 0) then
!!$            print *,'[input] hInterp vect fields: '//trim(F_varname_S )//' '//trim(F_varname2_S ) ; call flush(6)
            F_istat = ezuvint(F_outdata(:,:,F_k),F_outdata2(:,:,F_k),F_indata,F_indata2)
         else
            call msg(MSG_INFOPLUS, &
                 '(hinterp4yy2d_vect) Scalar Hinterp done for vect fields on samegrid for: '// &
                 trim(F_varname_S )//' '//trim(F_varname2_S ))
            F_istat = ezsint(F_outdata(:,:,F_k),F_indata)
            F_istat = min(ezsint(F_outdata2(:,:,F_k),F_indata2),F_istat)
        endif
      endif

      if (RMN_IS_OK(F_istat)) then
!!$         if (any(F_outdata(:,:,F_k) == huge(1.))) then
!!$            call msg(MSG_WARNING,'(hinterp4yy) Extrapolation occured for: '//trim(F_varname_S))
!!$            F_istat = RMN_ERR
!!$         else
            F_istat = HINTERP4YY_OK
            if (subgridid >= 0) F_istat = HINTERP4YY_NONE
!!$         endif
      endif
      call msg(MSG_DEBUG,'(hinterp4yy2d_vect) [END]')
      !------------------------------------------------------------------
      return
   end function hinterp4yy2d_vect


   !==== Private Functions =================================================


   !/@*
   subroutine priv_ezparams(F_gridsetid,F_subgridid,F_onesubgrid_S,F_h_int_S, &
        F_ingridid,F_outgridid,F_coregridid)
      integer,intent(out) :: F_gridsetid,F_subgridid
      character(len=*),intent(out) :: F_onesubgrid_S
      character(len=*),intent(inout) :: F_h_int_S
      integer,intent(in) :: F_ingridid,F_outgridid,F_coregridid
      !*@/
      integer,save :: gridsetid=-1,outgridid=-1,ingridid=-1,subgridid=-1
      character(len=16),save :: h_int0_S = ' '
      character(len=16),save :: onesubgrid_S = 'NO'
      integer :: istat
      !------------------------------------------------------------------
      F_gridsetid = gridsetid
      F_subgridid = subgridid
      F_onesubgrid_S = onesubgrid_S
      istat = clib_tolower(F_h_int_S)
!!$      if (F_outgridid == outgridid .and. F_ingridid == ingridid &
!!$           .and. gridsetid >= 0 .and. F_h_int_S == h_int0_S) then
!!$         if (subgridid >= 0) F_h_int_S = 'nearest'
!!$         return
!!$      endif
      h_int0_S = F_h_int_S
      outgridid = F_outgridid
      ingridid = F_ingridid
      gridsetid = ezdefset(outgridid, ingridid)
      if (gridsetid < 0) then
         call msg(MSG_WARNING,'(hinterp4yy2d) Unable to set an interpolation gridset')
         return
      endif

      onesubgrid_S='NO'
      subgridid=-1
      subgridid = ezgrid_subcolocated(ingridid, F_coregridid)
      if (subgridid >= 0)  then
         F_h_int_S = 'nearest'
         onesubgrid_S = 'YES'
      endif
      F_gridsetid = gridsetid
      F_subgridid = subgridid
      F_onesubgrid_S = onesubgrid_S
      !------------------------------------------------------------------
      return
   end subroutine priv_ezparams

end module hinterp4yy_mod
