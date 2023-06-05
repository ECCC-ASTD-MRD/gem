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

!/@
subroutine test_time_interp2()
use iso_c_binding
   use testutils
   use time_interp_mod
   use rmn_gmm
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04
!@/
#include <rmnlib_basics.hf>
   include "rpn_comm.inc"
   real,parameter :: EPSILON = 1.e-5
   character(len=512) :: dateo_S,varname_S,dummy_S,vgrid_S,sfcfld_S
   integer :: istat,datev,datev1,dateo,datev2,datevm1,dt,istep,nsteps,itype,datevtmp,myproc
   real :: weight
   real(REAL64) :: nhours_8,nhours2_8
   real,dimension(4,3,2),target :: data0v,data1v
   real,pointer,dimension(:,:,:):: data,data0,data1
   logical :: ok_L
   ! ---------------------------------------------------------------------
!!$   myproc = testutils_initmpi()
   call testutils_verbosity()
   call testutils_set_name('test_time_interp')
   data0 => data0v
   data1 => data1v

   data0 = 2.
   data1 = 4.
   varname_S = 'TT'
!!$   print *,'2 data1:',minval(data1),maxval(data1)

   do itype = 1,TIME_INTERP_NTYPES
      istat = time_interp_typecode(TIME_INTERP_TYPELIST_S(itype))
      call testutils_assert_eq(istat,itype,'typecode '//trim(TIME_INTERP_TYPELIST_S(itype)))
   enddo

   dateo=0
   dateo_S = '20090427.000000'
   call datp2f(dateo,dateo_S)
   nhours_8 = 2.D0
   datev1=dateo ;datev=dateo ; datev2=dateo ; datevm1=dateo
   call incdatr(datev1,dateo,nhours_8)
   call incdatr(datev,dateo,nhours_8/2.D0)
   call incdatr(datev2,dateo,nhours_8*2.D0)
   call incdatr(datevm1,dateo,-nhours_8)

   nullify(data)
   istat = time_interp_get(data,varname_S,dateo)
   print *,'dateo',istat,associated(data) ; call flush(6)
   nullify(data)
   istat = time_interp_get(data,varname_S,datev1)
   print *,'datev1',istat,associated(data) ; call flush(6)
   nullify(data)
   istat = time_interp_get(data,varname_S,datev)
   print *,'datev',istat,associated(data) ; call flush(6)
   if (istat == TIME_INTERP_NOT_FOUND) then
      istat = time_interp_set(data0,varname_S,dateo,'vgrid0','sfcfld0')
   elseif (istat == TIME_INTERP_NEED_NEXT) then
      istat = time_interp_set(data1,varname_S,datev1)
   endif

   nullify(data)
   istat = time_interp_get(data,varname_S,dateo)
   print *,'dateo',istat,associated(data) ; call flush(6)
   nullify(data)
   istat = time_interp_get(data,varname_S,datev1)
   print *,'datev1',istat,associated(data) ; call flush(6)
   nullify(data)
   istat = time_interp_get(data,varname_S,datev)
   print *,'datev',istat,associated(data) ; call flush(6)
   if (istat == TIME_INTERP_NOT_FOUND) then
      istat = time_interp_set(data0,varname_S,dateo,'vgrid0','sfcfld0')
   elseif (istat == TIME_INTERP_NEED_NEXT) then
      istat = time_interp_set(data1,varname_S,datev1)
   endif

   nullify(data)
   istat = time_interp_get(data,varname_S,dateo)
   print *,'dateo',istat,associated(data) ; call flush(6)
   nullify(data)
   istat = time_interp_get(data,varname_S,datev1)
   print *,'datev1',istat,associated(data) ; call flush(6)
   nullify(data)
   istat = time_interp_get(data,varname_S,datev)
   print *,'datev',istat,associated(data) ; call flush(6)
   if (istat == TIME_INTERP_NOT_FOUND) then
      istat = time_interp_set(data0,varname_S,dateo,'vgrid0','sfcfld0')
   elseif (istat == TIME_INTERP_NEED_NEXT) then
      istat = time_interp_set(data1,varname_S,datev1)
   endif

!!$   call rpn_comm_barrier(RPN_COMM_WORLD, istat)
!!$   call rpn_comm_finalize(istat)
   ! ---------------------------------------------------------------------
   return
end subroutine test_time_interp2
