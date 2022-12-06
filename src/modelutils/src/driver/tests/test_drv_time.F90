!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!-------------------------------------------------------------------------- 
!/@*
subroutine test_drv_time()
   use wb_itf_mod
   use testutils
   use cmcdate_mod
   use drv_time_mod, only: drv_time_config, drv_time_init, drv_time_increment
   implicit none
   !@objective Test drv_time_mod module functions
   !@author Stephane Chamberland, Feb 2012
!*@/
#include <rmnlib_basics.hf>

   integer, parameter :: NB_TIME_LEVELS = 2

   integer :: istat, step,step1,step2,dateo,datev,datev2
   logical :: is_last, is_chkpt,is_stat
   real(8) :: dt_8,nhours_8
   character(len=RMN_PATH_LEN) :: basename_S,msg_S,dateo_S
   !---------------------------------------------------------------------
   call testutils_verbosity()
   basename_S = 'test_drv_time'
   call testutils_set_name(basename_S)

   istat = drv_time_config(basename_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'drv_time_config status')
   istat = drv_time_init(dateo_S,dt_8,step,is_chkpt,is_last,NB_TIME_LEVELS)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'drv_time_init status')
   call testutils_assert_eq(dateo_S,'20080101.000000','drv_time_init dateo')
   call testutils_assert_eq(dt_8,75.D0,'drv_time_init dt_8')
   call testutils_assert_eq(step,0,'drv_time_init step')
   call testutils_assert_eq(is_chkpt,.false.,'drv_time_init is_chkpt')
   call testutils_assert_eq(is_last,.false.,'drv_time_init is_last')

   dateo = cmcdate_fromprint(dateo_S)

   step1 = 0
   istat = RMN_OK
   DO_STEP: do while (.not.is_last .and. istat == RMN_OK)
      step1 = step1 + 1
      write(msg_S,'(i2.2)') step1
      call drv_time_increment(step,is_chkpt,is_last,is_stat,datev)
      call testutils_assert_eq(step,step1,'increment step:'//trim(msg_S))
      call testutils_assert_eq(is_chkpt,mod(step1,3)==0,'increment is_chkpt:'//trim(msg_S))
      call testutils_assert_eq(is_stat,mod(step1,2)==0,'increment is_stat:'//trim(msg_S))
      call testutils_assert_eq(is_last,step1==6,'increment is_last:'//trim(msg_S))
      nhours_8 = dble(step1) * dt_8
      call incdatr(datev2,dateo,nhours_8)
      call testutils_assert_eq(datev2,datev,'increment datev:'//trim(msg_S))

      !TODO: test drv_time_shuffle
      !TODO: test drv_time_shuffle w/ copy
   end do DO_STEP
   call testutils_assert_eq(step,6,'increment step')

   call stop_mpi(RMN_OK,'test_drv_time',' ')
   !---------------------------------------------------------------------
   return
end subroutine test_drv_time
