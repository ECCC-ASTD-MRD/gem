!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer, 
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms 
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer 
!version 3 or (at your option) any later version that should be found at: 
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html 
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software; 
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec), 
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------
#include <rmn/msg.h>

!/@*
subroutine test_timestr()
   use testutils
   use timestr_mod
   use cmcdate_mod
   implicit none
   !@author S.Chamberland, 2011-09
   !*@/
#include <rmnlib_basics.hf>
   logical :: ok_L,hasloop_L
   integer :: istat,istat0,nstep,nbr_i,dateo,step,step2,maxstep,prognum
   integer(INT64) :: nbr_i8
   real :: nbr,nbr_r4,values(8),loop(3),interval,interval2,dt
   real(REAL64) :: nbr_8,dt_8,nbr_r8
   real(REAL64) :: dt2_8
   character(len=64) :: tmp_S,units_S
   !---------------------------------------------------------------------
   call testutils_verbosity()
   call testutils_set_name('timestr_mod')

   !# valid old
   istat = timestr_parse(nbr,units_S,' ')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'0')
   call testutils_assert_eq(nbr,    -1.,'0b')
   call testutils_assert_eq(units_S,'STE','0c')

   istat = timestr_parse(nbr,units_S,'132P')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'1')
   call testutils_assert_eq(nbr,    132.,'1b')
   call testutils_assert_eq(units_S,'STE','1c')

   istat = timestr_check('132P')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'1d')

   istat = timestr_parse(nbr,units_S,'132.P')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'2')
   call testutils_assert_eq(nbr,    132.,'2b')
   call testutils_assert_eq(units_S,'STE','2c')

   istat = timestr_parse(nbr,units_S,'132.S')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'3')
   call testutils_assert_eq(nbr,    132.,'')
   call testutils_assert_eq(units_S,'SEC','')

   istat = timestr_parse(nbr,units_S,'132.M')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'4')
   call testutils_assert_eq(nbr,    132.,'')
   call testutils_assert_eq(units_S,'MIN','')

   istat = timestr_parse(nbr,units_S,'132.H')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'5')
   call testutils_assert_eq(nbr,    132.,'')
   call testutils_assert_eq(units_S,'HOU','')

   istat = timestr_parse(nbr,units_S,'132.D')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'6')
   call testutils_assert_eq(nbr,    132.,'')
   call testutils_assert_eq(units_S,'DAY','')

   istat = timestr_parse(nbr,units_S,'132.')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'7')
   call testutils_assert_eq(nbr,    132.,'')
   call testutils_assert_eq(units_S,'HOU','')

   istat = timestr_parse(nbr,units_S,'132.','sec')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'8')
   call testutils_assert_eq(nbr,    132.,'')
   call testutils_assert_eq(units_S,'SEC','')

   !# valid new
   istat = timestr_parse(nbr,units_S,'step,132')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'9')
   call testutils_assert_eq(nbr,    132.,'')
   call testutils_assert_eq(units_S,'STE','')

   istat = timestr_parse(nbr,units_S,'sec, 132.')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'10')
   call testutils_assert_eq(nbr,    132.,'')
   call testutils_assert_eq(units_S,'SEC','')

   istat = timestr_parse(nbr,units_S,'min ,132.')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'11')
   call testutils_assert_eq(nbr,    132.,'')
   call testutils_assert_eq(units_S,'MIN','')

   istat = timestr_parse(nbr,units_S,'hour,132.')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'12')
   call testutils_assert_eq(nbr,    132.,'')
   call testutils_assert_eq(units_S,'HOU','')

   istat = timestr_check('hour,132.')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'12d')

   istat = timestr_parse(nbr,units_S,'day, 132.')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'13')
   call testutils_assert_eq(nbr,    132.,'')
   call testutils_assert_eq(units_S,'DAY','')

   istat = timestr_parse(nbr,units_S,'month,132.')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'14')
   call testutils_assert_eq(nbr,    132.,'')
   call testutils_assert_eq(units_S,'MON','')

   istat = timestr_parse(nbr,units_S,'132.')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'15')
   call testutils_assert_eq(nbr,    132.,'')
   call testutils_assert_eq(units_S,'HOU','')

   istat = timestr_parse(nbr,units_S,',132.','day')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'16')
   call testutils_assert_eq(nbr,    132.,'')
   call testutils_assert_eq(units_S,'DAY','')

   istat = timestr_check('132.')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'16d')

   istat = timestr_parse(nbr,units_S,'132.','month')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'17')
   call testutils_assert_eq(nbr,    132.,'')
   call testutils_assert_eq(units_S,'MON','')

   !# full
   istat = timestr_parse(values,loop,hasloop_L,units_S,'min ,132.,133.')
   call testutils_assert_eq(istat,2,'17.1')
   call testutils_assert_eq(values(1),132.,'')
   call testutils_assert_eq(values(2),133.,'')
   call testutils_assert_eq(hasloop_L,.false.,'')
   istat = timestr_parse(values,loop,hasloop_L,units_S,'min ,[132.,133.]')
   call testutils_assert_eq(istat,2,'17.2')
   call testutils_assert_eq(values(1),132.,'')
   call testutils_assert_eq(values(2),133.,'')
   call testutils_assert_eq(hasloop_L,.false.,'')
   istat = timestr_parse(values,loop,hasloop_L,units_S,'min ,<0,133.,1>')
   call testutils_assert_eq(istat,0,'17.3')
   call testutils_assert_eq(loop(1),0.,'')
   call testutils_assert_eq(loop(2),133.,'')
   call testutils_assert_eq(loop(3),1.,'')
   call testutils_assert_eq(hasloop_L,.true.,'')

   !# invalid
   istat = timestr_check('132.N')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'18a')

   istat = timestr_parse(nbr,units_S,'132.N')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'18')
   istat = timestr_parse(nbr,units_S,'132.Month')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'19')
   istat = timestr_parse(nbr,units_S,'13.2.M')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'20')
   istat = timestr_parse(nbr,units_S,'m,132.')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'21')
   istat = timestr_parse(nbr,units_S,'sec132.')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'22')
   istat = timestr_parse(nbr,units_S,'sec,13.2.')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'23')
   istat = timestr_parse(nbr,units_S,'Month')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'24')
   istat = timestr_parse(nbr,units_S,'Month,')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'25')

   !# full invalid
   istat = timestr_parse(nbr,units_S,'min ,132.,133.')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'25.1')
   istat = timestr_parse(nbr,units_S,'min ,[132.,133.]')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'25.2')
   istat = timestr_parse(nbr,units_S,'min ,<0,133.,1>')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'25.3')


   !# timestr2sec
   istat = timestr2sec(nbr,'2')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'26')
   call testutils_assert_eq(nbr,    7200.,'')

   istat = timestr2sec(nbr_8,'2')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'26_r8')
   call testutils_assert_eq(nbr_8==7200.,.true.,'')

   istat = timestr2sec(nbr_i,'2')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'26_i4')
   call testutils_assert_eq(nbr_i==7200,.true.,'')

   istat = timestr2sec(nbr_i8,'2')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'26_i8')
   call testutils_assert_eq(nbr_i8==7200,.true.,'')


   istat = timestr2sec(nbr,'2.M')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'27')
   call testutils_assert_eq(nbr,    120.,'')

   istat = timestr2sec(nbr,'2P')
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'28')
   call testutils_assert_eq(nbr,    2.,'')
   
   dt_8 = 10.
   dt2_8 = 10.
   istat = timestr2sec(nbr,'2P',F_dt=dt_8)
   istat = timestr2sec(nbr,'2P',F_dt=dt2_8)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'29')
   call testutils_assert_eq(nbr,    20.,'')

   istat = timestr2sec(nbr,'2N')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'30')


   !# timestr2step
   dt_8 = 30.
   dt2_8 = 30.
   istat = timestr2step(nstep,'2',dt2_8)
   istat = timestr2step(nstep,'2',dt_8)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'31')
   call testutils_assert_eq(nstep,    240,'')

   istat = timestr2step(nbr_i8,'2',dt_8)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'31_i8')
   call testutils_assert_eq(nbr_i8==240,.true.,'')

   istat = timestr2step(nbr_r4,'2',dt_8)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'31_r4')
   call testutils_assert_eq(nbr_r4==240.,.true.,'')

   istat = timestr2step(nbr_r8,'2',dt_8)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'31_r8')
   call testutils_assert_eq(nbr_r8==240.,.true.,'')

   istat = timestr2step(nstep,'2.M',dt_8)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'32')
   call testutils_assert_eq(nstep,    4,'')

   istat = timestr2step(nstep,'2P',dt_8)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'33')
   call testutils_assert_eq(nstep,    2,'')

   istat = timestr2step(nstep,'2N',dt_8)
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'34')


   !# timestr_prognum_int
   dateo   = cmcdate_fromprint('20150430.000000')
   dt      = 1800.
   maxstep = 12

   units_S  = 'STE'
   interval = 1.
   interval2 = interval
   do step=0,15
      step2 = ceiling(min(step,maxstep)/interval2)*interval
      istat = timestr_prognum(prognum,units_S,interval,dateo,dt,step,maxstep)
      call testutils_assert_eq(RMN_IS_OK(istat),.true.,'40')
      call testutils_assert_eq(prognum,step2,'')
   enddo

   units_S  = 'STE'
   interval = 2.
   interval2 = interval
   do step=0,14
      step2 = ceiling(min(step,maxstep)/interval2)*interval
      istat = timestr_prognum(prognum,units_S,interval,dateo,dt,step,maxstep)
      call testutils_assert_eq(RMN_IS_OK(istat),.true.,'41')
      call testutils_assert_eq(prognum,step2,'')
   enddo

   units_S  = 'HOU'
   interval = 1.
   interval2 = interval*3600.
   do step=0,13
      step2 = ceiling(min(step,maxstep)*dt/interval2)*interval
      istat = timestr_prognum(prognum,units_S,interval,dateo,dt,step,maxstep)
      call testutils_assert_eq(RMN_IS_OK(istat),.true.,'42')
      call testutils_assert_eq(prognum,step2,'')
   enddo

   units_S  = 'HOU'
   interval = 2.
   interval2 = interval*3600.
   do step=0,14
      step2 = ceiling(min(step,maxstep)*dt/interval2)*interval
      istat = timestr_prognum(prognum,units_S,interval,dateo,dt,step,maxstep)
      call testutils_assert_eq(RMN_IS_OK(istat),.true.,'43')
      call testutils_assert_eq(prognum,step2,'')
   enddo

   units_S  = 'MON'
   interval = 1.
   do step=0,2
      step2 = 0
      istat = timestr_prognum(prognum,units_S,interval,dateo,dt,step,maxstep)
      call testutils_assert_eq(RMN_IS_OK(istat),.true.,'44')
      call testutils_assert_eq(prognum,step2,'')
   enddo
   do step=47,49
      step2 = 0
      if (step >= 48) step2 = 1
      istat = timestr_prognum(prognum,units_S,interval,dateo,dt,step,maxstep)
      call testutils_assert_eq(RMN_IS_OK(istat),.true.,'45')
      call testutils_assert_eq(prognum,step2,'')
   enddo


   !# timestr_isstep
   dateo   = cmcdate_fromprint('20150430.000000')
   dt      = 1800.
   maxstep = 12
   tmp_S  = 'STEP,2'
   do step=0,15
      istat0 = TIMESTR_NO_MATCH
      if (step >= maxstep .or. mod(step,2) == 0) istat0 = TIMESTR_MATCH
      istat = timestr_isstep(tmp_S,dateo,dt,step,maxstep)
      call testutils_assert_eq(istat,istat0,'50')
   enddo

   maxstep = 51
   tmp_S  = 'MONTH,1'
   do step=0,2
      istat = timestr_isstep(tmp_S,dateo,dt,step,maxstep)
      call testutils_assert_eq(istat,TIMESTR_NO_MATCH,'51')
   enddo
   do step=46,51
      istat0 = TIMESTR_NO_MATCH
      if (step == 47 .or. step >=maxstep) istat0 = TIMESTR_MATCH
      istat = timestr_isstep(tmp_S,dateo,dt,step,maxstep)
      call testutils_assert_eq(istat,istat0,'51b')
   enddo
   do step=46,51
      istat0 = TIMESTR_NO_MATCH
      if (step == 47) istat0 = TIMESTR_MATCH
      istat = timestr_isstep(tmp_S,dateo,dt,step)
      call testutils_assert_eq(istat,istat0,'51c')
   enddo

   !---------------------------------------------------------------------
   return
end subroutine test_timestr
