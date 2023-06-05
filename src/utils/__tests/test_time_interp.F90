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
subroutine test_time_interp()
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
   integer :: istat,datev,datev1,dateo,datev2,datevm1,dt,istep,nsteps,itype,datevtmp,myproc,tinterp_type,trials,istep_prev,istep_next
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
   call testutils_assert_ok(istat == TIME_INTERP_NOT_FOUND,'dateo TIME_INTERP_NOT_FOUND')

   istat = time_interp_set(data0,varname_S,dateo,'vgrid0','sfcfld0')
   call testutils_assert_ok(RMN_IS_OK(istat),'set dateo')

   nullify(data)
   istat = time_interp_get(data,trim(varname_S)//'_none_',dateo)
   call testutils_assert_ok(istat == TIME_INTERP_NOT_FOUND,'var TIME_INTERP_NOT_FOUND')

   istat = time_interp_retrieve(varname_S,TIME_INTERP_PREV,datevtmp,F_data=data,F_vgrid_S=vgrid_S,F_sfcfld_S=sfcfld_S)
   call testutils_assert_ok(RMN_IS_OK(istat),'retrieve prev')
   call testutils_assert_ok(datevtmp==dateo,'retrieve prev values date')
   call testutils_assert_ok((vgrid_S=='vgrid0'.and.sfcfld_S=='sfcfld0'),'retrieve prev values str:'//trim(vgrid_S)//':'//trim(sfcfld_S)//':')
   call testutils_assert_ok(all(data==data0),'retrieve prev values data')

   istat = time_interp_retrieve(varname_S,TIME_INTERP_NEXT,datevtmp,F_data=data,F_vgrid_S=vgrid_S,F_sfcfld_S=sfcfld_S)
   call testutils_assert_ok(.not.RMN_IS_OK(istat),'retrieve next')

   weight = time_interp_weight(varname_S,dateo,istat)
   call testutils_assert_ok(RMN_IS_OK(istat),'weight dateo 0')
   call testutils_assert_ok(weight==0.,'weight dateo 0 value')

   istat = time_interp_status(varname_S,dateo)
   call testutils_assert_ok(RMN_IS_OK(istat),'status dateo 0')

   !#------------------------------

   allocate(data(24,1,1),stat=istat)
   data = 0.
   istat = time_interp_set(data,varname_S,datev1)   
   call testutils_assert_ok(.not.RMN_IS_OK(istat),'set wrong dims')   

   istat = time_interp_retrieve(varname_S,TIME_INTERP_PREV,datevtmp,F_data=data,F_vgrid_S=vgrid_S,F_sfcfld_S=sfcfld_S)
   call testutils_assert_ok(RMN_IS_OK(istat),'retrieve prev post')
   call testutils_assert_ok(datevtmp==dateo,'retrieve prev post values date')
   call testutils_assert_ok((vgrid_S=='vgrid0'.and.sfcfld_S=='sfcfld0'),'retrieve prev post values str:'//trim(vgrid_S)//':'//trim(sfcfld_S)//':')
   call testutils_assert_ok(all(data==data0),'retrieve prev post values data')

   istat = time_interp_retrieve(varname_S,TIME_INTERP_NEXT,datevtmp,F_data=data,F_vgrid_S=vgrid_S,F_sfcfld_S=sfcfld_S)
   call testutils_assert_ok(.not.RMN_IS_OK(istat),'retrieve next post')

   weight = time_interp_weight(varname_S,dateo,istat)
   call testutils_assert_ok(RMN_IS_OK(istat),'weight dateo 0 post')
   call testutils_assert_ok(weight==0.,'weight dateo 0 value post')

   istat = time_interp_status(varname_S,dateo)
   call testutils_assert_ok(RMN_IS_OK(istat),'status dateo 0 post')

   nullify(data)
   istat = time_interp_get(data,varname_S,dateo)
   call testutils_assert_ok(RMN_IS_OK(istat),'dateo 0')
   ok_L = .false.
   if (associated(data)) ok_L = any(abs(data-data0) < EPSILON)
   call testutils_assert_ok(ok_L,'dateo 0 values')

   nullify(data)
   istat = time_interp_get(data,varname_S,datev)
   call testutils_assert_ok(istat == TIME_INTERP_NEED_NEXT,'datev TIME_INTERP_NEED_NEXT')

   !#------------------------------

   istat = time_interp_set(data1,varname_S,datev1,'vgrid1','sfcfld1')
   call testutils_assert_ok(RMN_IS_OK(istat),'set datev1')

   istat = time_interp_retrieve(varname_S,TIME_INTERP_PREV,datevtmp,F_data=data,F_vgrid_S=vgrid_S,F_sfcfld_S=sfcfld_S)
   call testutils_assert_ok(RMN_IS_OK(istat),'retrieve prev 2')
   call testutils_assert_ok(datevtmp==dateo,'retrieve prev 2 values date')
   call testutils_assert_ok((vgrid_S=='vgrid0'.and.sfcfld_S=='sfcfld0'),'retrieve prev 2 values str:'//trim(vgrid_S)//':'//trim(sfcfld_S)//':')
   call testutils_assert_ok(all(data==data0),'retrieve prev 2 values data')

   istat = time_interp_retrieve(varname_S,TIME_INTERP_NEXT,datevtmp,F_data=data,F_vgrid_S=vgrid_S,F_sfcfld_S=sfcfld_S)
   call testutils_assert_ok(RMN_IS_OK(istat),'retrieve next 2')
   call testutils_assert_ok(datevtmp==datev1,'retrieve next 2 values date')
   call testutils_assert_ok((vgrid_S=='vgrid1'.and.sfcfld_S=='sfcfld1'),'retrieve next 2 values str:'//trim(vgrid_S)//':'//trim(sfcfld_S)//':')
   call testutils_assert_ok(all(data==data1),'retrieve next 2 values data')

   weight = time_interp_weight(varname_S,dateo,istat)
   call testutils_assert_ok(RMN_IS_OK(istat),'weight dateo 0 2')
   call testutils_assert_ok(weight==0.,'weight dateo 0 value 2')

   weight = time_interp_weight(varname_S,datev,istat)
   call testutils_assert_ok(RMN_IS_OK(istat),'weight datev 0 2')
   call testutils_assert_ok(weight==0.5,'weight datev 0 value 2')

   weight = time_interp_weight(varname_S,datev1,istat)
   call testutils_assert_ok(RMN_IS_OK(istat),'weight datev1 0 2')
   call testutils_assert_ok(weight==1.,'weight datev1 0 value 2')

   istat = time_interp_status(varname_S,dateo)
   call testutils_assert_ok(RMN_IS_OK(istat),'status dateo 0 2')


   nullify(data)
   istat = time_interp_get(data,varname_S,dateo)
   call testutils_assert_ok(RMN_IS_OK(istat),'dateo 1')   
   ok_L = .false.
   if (associated(data)) ok_L = any(abs(data-data0) < EPSILON)
   call testutils_assert_ok(ok_L,'dateo 1 values')

   nullify(data)
   istat = time_interp_get(data,varname_S,datev)   
   call testutils_assert_ok(RMN_IS_OK(istat),'datev')   
   ok_L = .false.
   if (associated(data)) ok_L = any(abs(data-(0.5*(data0+data1))) < EPSILON)
   call testutils_assert_ok(ok_L,'datev values')


   nullify(data)
   istat = time_interp_get(data,varname_S,datev1)   
   call testutils_assert_ok(RMN_IS_OK(istat),'datev1')   
   ok_L = .false.
   if (associated(data)) ok_L = any(abs(data-data1) < EPSILON)
   call testutils_assert_ok(ok_L,'datev1 values')

   istat = time_interp_get(data,varname_S,datev1)   
   call testutils_assert_ok(RMN_IS_OK(istat),'datev1 allocated')   
   ok_L = .false.
   if (associated(data)) ok_L = any(abs(data-data1) < EPSILON)
   call testutils_assert_ok(ok_L,'datev1 values allocated')

   !#------------------------------

   allocate(data(24,1,1),stat=istat)
   istat = time_interp_get(data,varname_S,datev1)   
   call testutils_assert_ok(.not.RMN_IS_OK(istat),'datev1 wrong dims')   

   nullify(data)
   istat = time_interp_get(data,varname_S,datev2)
   call testutils_assert_ok(istat == TIME_INTERP_NEED_NEXT,'datev2 TIME_INTERP_NEED_NEXT')

   nullify(data)
   istat = time_interp_get(data,varname_S,datevm1)   
   call testutils_assert_ok(istat == TIME_INTERP_NEED_PREV,'datevm1 TIME_INTERP_NEED_PREV')


   !#------------------------------

   call incdatr(datev,dateo,nhours_8/3.D0)


   nullify(data)
   istat = time_interp_get(data,varname_S,datev,TIME_INTERP_NEAR)   
   call testutils_assert_ok(RMN_IS_OK(istat),'datev NEAR 0')
   ok_L = .false.
   if (associated(data)) ok_L = any(abs(data-data0) < EPSILON)
   call testutils_assert_ok(ok_L,'datev NEAR 0 values')


   nullify(data)
   istat = time_interp_get(data,varname_S,datev,TIME_INTERP_STEP)   
   call testutils_assert_ok(RMN_IS_OK(istat),'datev STEP 0')
   ok_L = .false.
   if (associated(data)) ok_L = any(abs(data-data0) < EPSILON)
   call testutils_assert_ok(ok_L,'datev STEP 0 values')


   nullify(data)
   istat = time_interp_get(data,varname_S,datev,TIME_INTERP_NEXT)   
   call testutils_assert_ok(RMN_IS_OK(istat),'datev NEXT 0')
   ok_L = .false.
   if (associated(data)) ok_L = any(abs(data-data1) < EPSILON)
   call testutils_assert_ok(ok_L,'datev NEXT 0 values')


   nullify(data)
   istat = time_interp_get(data,varname_S,datev,TIME_INTERP_INCR,3600)   
   call testutils_assert_ok(RMN_IS_OK(istat),'datev INCR 0')
   ok_L = .false.
   if (associated(data)) ok_L = any(abs(data-abs(data1-data0)/2.) < EPSILON)
   call testutils_assert_ok(ok_L,'datev INCR 0 values')



   call incdatr(datev,dateo,2.d0*nhours_8/3.D0)
   nullify(data)
   istat = time_interp_get(data,varname_S,datev,TIME_INTERP_NEAR)   
   call testutils_assert_ok(RMN_IS_OK(istat),'datev NEAR 1')
   ok_L = .false.
   if (associated(data)) ok_L = any(abs(data-data1) < EPSILON)
   call testutils_assert_ok(ok_L,'datev NEAR 1 values')

   nullify(data)
   istat = time_interp_get(data,varname_S,datev,TIME_INTERP_STEP)   
   call testutils_assert_ok(RMN_IS_OK(istat),'datev STEP 1')
   ok_L = .false.
   if (associated(data)) ok_L = any(abs(data-data0) < EPSILON)
   call testutils_assert_ok(ok_L,'datev STEP 1 values')


   nullify(data)
   istat = time_interp_get(data,varname_S,datev,TIME_INTERP_NEXT)   
   call testutils_assert_ok(RMN_IS_OK(istat),'datev NEXT 1')
   ok_L = .false.
   if (associated(data)) ok_L = any(abs(data-data1) < EPSILON)
   call testutils_assert_ok(ok_L,'datev NEXT 1 values')


   nullify(data)
   istat = time_interp_get(data,varname_S,datev,TIME_INTERP_INCR,3600)   
   call testutils_assert_ok(RMN_IS_OK(istat),'datev INCR 1')
   ok_L = .false.
   if (associated(data)) ok_L = any(abs(data-abs(data1-data0)/2.) < EPSILON)
   call testutils_assert_ok(ok_L,'datev INCR 1 values')


   istat = time_interp_retrieve(varname_S,TIME_INTERP_PREV,datevtmp,F_vgrid_S=vgrid_S,F_sfcfld_S=sfcfld_S)
   call testutils_assert_ok(RMN_IS_OK(istat),'retrieve prev')
   call testutils_assert_ok((vgrid_S=='vgrid0'.and.sfcfld_S=='sfcfld0'),'retrieve prev values')

   istat = time_interp_retrieve(varname_S,TIME_INTERP_NEXT,datevtmp,F_vgrid_S=vgrid_S,F_sfcfld_S=sfcfld_S)
   call testutils_assert_ok(RMN_IS_OK(istat),'retrieve next')
   call testutils_assert_ok((vgrid_S=='vgrid1'.and.sfcfld_S=='sfcfld1'),'retrieve next values')


   tinterp_type = TIME_INTERP_LINE
   varname_S = 'P0'
   dt = 3600
   nsteps = 4
   istep_prev = -99
   istep_next = -99
   do istep=0,nsteps
      write(dummy_S,'(a,i2)') ' istep=',istep
      
      nhours_8 = dble(dt*istep)/3600.D0
      call incdatr(datev,dateo,nhours_8)

      trials = 0
      DO_TRIALS: do
         istat = time_interp_status(varname_S,datev,tinterp_type)
         if (istat==TIME_INTERP_NOT_FOUND .or. &
              istat==TIME_INTERP_NEED_PREV) then
            istep_prev = (istep-1)*2
            nhours2_8 = dble(dt*istep_prev)/3600.D0
            call incdatr(datev2,dateo,nhours2_8)
            allocate(data(4,3,2),stat=istat)
            data = real(istep_prev)
            istat = time_interp_set(data,varname_S,datev2)
            call testutils_assert_ok(RMN_IS_OK(istat),'set'//trim(dummy_S))
!!$            print *,'test_time_interp se-',istep,istep_prev,maxval(data)
            deallocate(data,stat=istat)         
         else if (istat==TIME_INTERP_NEED_NEXT) then
            if (istep_next /= -99) istep_prev = istep_next
            istep_next = istep*2
            nhours2_8 = dble(dt*istep_next)/3600.D0
            call incdatr(datev2,dateo,nhours2_8)
            allocate(data(4,3,2),stat=istat)
            data = real(istep_next)
            istat = time_interp_set(data,varname_S,datev2)
            call testutils_assert_ok(RMN_IS_OK(istat),'set'//trim(dummy_S))
!!$            print *,'test_time_interp se+',istep,istep_next,maxval(data)
            deallocate(data,stat=istat)
         else
            exit
         endif
         trials = trials + 1
         if (trials > 2) exit
      enddo DO_TRIALS

      istat = time_interp_retrieve(varname_S,TIME_INTERP_PREV,datevtmp,F_data=data)
      call testutils_assert_ok(RMN_IS_OK(istat),'retrieve prev'//trim(dummy_S))
      ok_L = .false.
      data1 = real(istep_prev)
      if (associated(data)) ok_L = any(abs(data-data1) < EPSILON)
      call testutils_assert_ok(ok_L,'retrieve prev values'//trim(dummy_S))
!!$      print *,'test_time_interp re-',istep,istep_next,maxval(data)
            
      istat = time_interp_retrieve(varname_S,TIME_INTERP_NEXT,datevtmp,F_data=data)
      call testutils_assert_ok(RMN_IS_OK(istat),'retrieve next'//trim(dummy_S))
      ok_L = .false.
      data1 = real(istep_next)
      if (associated(data)) ok_L = any(abs(data-data1) < EPSILON)
      call testutils_assert_ok(ok_L,'retrieve next values'//trim(dummy_S))
!!$      print *,'test_time_interp re+',istep,istep_next,maxval(data)

      nullify(data)
      istat = time_interp_get(data,varname_S,datev)
      call testutils_assert_ok(RMN_IS_OK(istat),'get'//trim(dummy_S))
      ok_L = .false.
      data1 = real(istep)
      if (associated(data)) ok_L = any(abs(data-data1) < EPSILON)
      call testutils_assert_ok(ok_L,'datev values'//trim(dummy_S))
   enddo




   varname_S = 'HU'
   do istep=nsteps,0,-1
      write(dummy_S,'(a,i2)') ' backward istep=',istep
      
      nhours_8 = dble(dt*istep)/3600.D0
      call incdatr(datev,dateo,nhours_8)

      nullify(data)
      istat = time_interp_get(data,varname_S,datev)

      if (istat==TIME_INTERP_NOT_FOUND .or. &
           istat==TIME_INTERP_NEED_NEXT .or. &
           istat==TIME_INTERP_NEED_PREV) then
         nhours2_8 = dble(dt*istep*2)/3600.D0
         call incdatr(datev2,dateo,nhours2_8)
         allocate(data(4,3,2),stat=istat)
         data = real(istep*2)
         istat = time_interp_set(data,varname_S,datev2)
         call testutils_assert_ok(RMN_IS_OK(istat),'set'//trim(dummy_S))

         nullify(data)
         istat = time_interp_get(data,varname_S,datev)
         call testutils_assert_ok(RMN_IS_OK(istat),'get'//trim(dummy_S))
      else
         call testutils_assert_ok(RMN_IS_OK(istat),'get'//trim(dummy_S))
      endif
      ok_L = .false.
      data1 = real(istep)
      if (associated(data)) ok_L = any(abs(data-data1) < EPSILON)

      call testutils_assert_ok(ok_L,'datev values'//trim(dummy_S))
   enddo
   ! ---------------------------------------------------------------------
   return
end subroutine test_time_interp
