!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!--------------------------------------------------------------------------

!/@*
module drv_time_mod
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use clib_itf_mod, only: clib_tolower
   use wb_itf_mod
   use config_mod
   use cmcdate_mod
   use sort_mod
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <mu_gmm.hf>
#include <msg.h>
   private
   !@objective Manage time steping
   !@author
   !  Michel Desgagne, Feb 2008
   !  Ron McTaggart-Cowan, Feb 2008
   !  Stephane Chamberland, Feb 2008
   !@revisions
   !  2012-02, Stephane Chamberland: RPNPhy offline
   !  2012-03, Stephane Chamberland: Generalize time levels flags
   !@public_functions
   public :: drv_time_config, drv_time_init, drv_time_info, drv_time_increment, &
        drv_time_shuffle_name,drv_time_shuffle,drv_time_shuffle_list
   !@public_vars
   character(len=15), public :: time_run_start !- YYYYMMDD.hhmmss
   real(REAL64), public,save :: time_dt       !- time step [sec]
   integer, public,save :: time_stepno_stat    !- echo stats every time_stepno_stat steps
   !@public_params
   integer,parameter,public :: DRV_TIME_MODE_SHUFFLE = 0
   integer,parameter,public :: DRV_TIME_MODE_COPY = 1
   character(len=1),parameter,public :: DRV_TIME_FLAG_2LVL(2) = (/'+','-'/)
   character(len=1),parameter,public :: DRV_TIME_FLAG_3LVL(3) = (/'+','0','-'/)
   !@description
   !  This module read time settings and inititalize time values
   !  It also offer function to increment/decrement the time step number
   !  while keeping GMM labels in sync [gmm_shuffle]
   !  Some useful time values are put in the global WB for other component 
   !  to access [read only]
   !  time_run_start : Initial start date [r8], YYYYMMDD.hhmmss
   !  time_dt        : Time step [r8], seconds
   !  time_stepno    : Present Step number [I]
   !  time_datetime  : Actual Step date [r8], YYYYMMDD.hhmmss
!*@/

   !- Private module vars
   integer,parameter :: MAX_FIELDS  = 2048
   integer,parameter :: TIME_FLAGS_MAX = 8

   integer,save :: time_stepno       !- current step number
   integer,save :: time_nb_levels    !- number of time levels [dyn dependent]
   integer,save :: time_stepno_tot   !- tot number of steps
   integer,save :: time_stepno_chkpt !- wrtie a rstrt file every time_stepno_chkpt steps
   integer,save :: time_stepno_rstrt !- restart every time_stepno_rstrt steps
   integer,save :: time_dateo        !- Date of run start [CMC datetime stamp]
   integer,save :: stepno_last       !- 
   integer,save :: m_swap_mode       !- Copy or shuffle tie swap mode
   logical,save :: m_init_L = .false.

contains

   !/@*
   function drv_time_config(F_cfg_basename_S) result(F_istat)
      implicit none
      !@objective Read time config from file to WB
      !@arguments
      character(len=*),intent(in) :: F_cfg_basename_S  !- Name of the config file [w/o extension]
      !@return
      integer :: F_istat
      !@author
      !  Stephane Chamberland, Feb 2008
      !*@/
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'[BEGIN] drv_time_config')
      F_istat = RMN_OK
      if (m_init_L) return
      F_istat = config_read(F_cfg_basename_S,'time_cfgs')
      call msg(MSG_DEBUG,'[END] drv_time_config')
      !---------------------------------------------------------------------
      return
   end function drv_time_config


   !/@*
   function drv_time_init(F_runstrt_S,F_dt_8,F_stepno,F_is_chkpt_L,F_is_last_L,F_ntime_levels,F_datev,F_swap_mode) result(F_istat)
      implicit none
      !@objective Initialize time config
      !@arguments
      character(len=*),intent(out) :: F_runstrt_S
      real(REAL64),intent(out) :: F_dt_8
      integer,intent(out) :: F_stepno        !- current step
      logical,intent(out) :: F_is_chkpt_L  !- .T. if it's a checkpoint step
      logical,intent(out) :: F_is_last_L   !- .T. if it's the last step [restart of run end]
      integer,intent(in),optional :: F_ntime_levels !may want to provide a callback fn instead
      integer,intent(out),optional :: F_datev !- Actual valid date/time at F_stepno [CMC datetime stamp]
      integer,intent(in),optional :: F_swap_mode
      !TODO-FR: add time_flags optional arg
      !@return
      integer :: F_istat
      !@author
      !  Michel Desgagne, Feb 2008
      !  Ron McTaggart-Cowan, Feb 2008
      !  Stephane Chamberland, Feb 2008
      !*@/
      integer :: istat,datev
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'[BEGIN] drv_time_init')
      istat = RMN_OK
      if (m_init_L) return
      m_init_L = .true.

      !- Handle 1st time slice [no Restart mode]
      istat = wb_get('time_stepno',time_stepno)
      if (.not.RMN_IS_OK(istat)) then
         time_stepno = 0
         time_nb_levels = 2
         if (present(F_ntime_levels)) time_nb_levels = F_ntime_levels
         istat = wb_put('time_cfgs/time_nb_levels',time_nb_levels)
      else
         istat = wb_get('time_cfgs/time_nb_levels',time_nb_levels)
      endif

      !- Init module vars
      istat = min(wb_get('time_cfgs/step_runstrt_S',time_run_start),istat)
      istat = min(wb_get('time_cfgs/step_dt',time_dt),istat)
      istat = min(wb_get('time_cfgs/step_total',time_stepno_tot),istat)
      istat = min(wb_get('time_cfgs/step_rsti',time_stepno_rstrt),istat)
      istat = min(wb_get('time_cfgs/step_bkup',time_stepno_chkpt),istat)
      istat = min(wb_get('time_cfgs/step_gstat',time_stepno_stat),istat)

      time_dateo = cmcdate_fromprint(time_run_start)

      F_runstrt_S = time_run_start
      F_dt_8 = time_dt
      F_stepno = time_stepno
      F_is_chkpt_L = .false.
      if (time_stepno_rstrt <= 0) time_stepno_rstrt = time_stepno_tot
      stepno_last = min(time_stepno + time_stepno_rstrt, time_stepno_tot)
      F_is_last_L = (time_stepno >= stepno_last)

      call drv_time_info(F_datev=datev)
      if (present(F_datev)) F_datev = datev
 
      m_swap_mode = DRV_TIME_MODE_SHUFFLE
      if (present(F_swap_mode)) m_swap_mode = F_swap_mode

      !- put usefull time values in global WB for other components to access
      istat = min(wb_put('time_run_start',time_run_start),istat)
      istat = min(wb_put('time_dt',time_dt),istat)
      istat = min(wb_put('time_stepno',time_stepno,WB_REWRITE_MANY),istat)
      istat = min(wb_put('time_datev',datev,WB_REWRITE_MANY),istat)
      F_istat = istat
      if (RMN_IS_OK(F_istat)) then
         call msg(MSG_INFO,'(drv_time) Initialisation OK')
      else
         m_init_L = .false.
         call msg(MSG_ERROR,'(drv_time) Problem in Initialisation')
      endif
      call msg(MSG_DEBUG,'[END] drv_time_init')
      !---------------------------------------------------------------------
      return
   end function drv_time_init


   !/@*
   subroutine drv_time_increment(F_stepno,F_is_chkpt_L,F_is_last_L,F_is_stat_L,F_datev,F_timeFlags_S,F_delta_step,F_swap_mode)
      implicit none
      !@objective Increase time step and update GMM fields names
      !@arguments
      integer,intent(out) :: F_stepno       !- current step
      logical,intent(out) :: F_is_chkpt_L  !- .T. if it's a checkpoint step
      logical,intent(out) :: F_is_last_L   !- .T. if it's the last step [restart of run end]
      logical,intent(out),optional :: F_is_stat_L   !- .T. if it's time to do print stats
      integer,intent(out),optional :: F_datev !- Actual valid date/time at F_stepno [CMC datetime stamp]
      character(len=*),intent(in),optional :: F_timeFlags_S(:)
      integer,intent(in),optional :: F_delta_step
      integer,intent(in),optional :: F_swap_mode
     !@author
      !  Michel Desgagne, Feb 2008
      !  Ron McTaggart-Cowan, Feb 2008
      !  Stephane Chamberland, Feb 2008
   !*@/
      logical :: is_stat_L
      integer :: istat,nlvl,datev,delta_step,nn,swap_mode
      character(len=2) :: timeFlags_S(TIME_FLAGS_MAX)
      character(len=16) :: datev_S
      character(len=256) :: msg_S
      !---------------------------------------------------------------------
      if (.not.m_init_L) then
         call msg(MSG_ERROR,'(drv_time) Not Initialised')
         return
      endif

      delta_step = 1
      if (present(F_delta_step)) delta_step = F_delta_step
      time_stepno = time_stepno + delta_step

      if (present(F_timeFlags_S)) then
         nlvl = min(size(F_timeFlags_S),TIME_FLAGS_MAX)
         timeFlags_S(1:nlvl) = F_timeFlags_S(1:nlvl)
      else
         if (time_nb_levels == 3) then
            nlvl = min(size(DRV_TIME_FLAG_3lvl),TIME_FLAGS_MAX)
            timeFlags_S(1:nlvl) = DRV_TIME_FLAG_3lvl(1:nlvl)
         else !if (time_nb_levels == 2) then
            nlvl = min(size(DRV_TIME_FLAG_2lvl),TIME_FLAGS_MAX)
            timeFlags_S(1:nlvl) = DRV_TIME_FLAG_2lvl(1:nlvl)
         endif
      endif

      swap_mode = m_swap_mode
      if (present(F_swap_mode)) swap_mode = F_swap_mode

      if (delta_step /= abs(delta_step)) call reverse(timeFlags_S(1:nlvl))
      do nn=1,abs(delta_step)
         call drv_time_shuffle(timeFlags_S(1:nlvl),swap_mode)
      enddo

      call drv_time_info(F_stepno,F_is_chkpt_L,F_is_last_L,is_stat_L,datev)
      if (present(F_is_stat_L)) F_is_stat_L = is_stat_L
      if (present(F_datev)) F_datev = datev

      datev_S = cmcdate_toprint(datev)
      write(msg_S,'("Step=",i4.4," (datev=",a,") (ckpt,last,stat)=(",l2,l2,l2,")")') &
           F_stepno,datev_S(1:15),F_is_chkpt_L,F_is_last_L,is_stat_L
      call msg(MSG_INFO,'(drv_time) Increment: '//trim(msg_S))

      !- update useful time values in global WB for other components to access
      istat = wb_put('time_stepno',time_stepno,WB_REWRITE_MANY)
      istat = wb_put('time_datev',datev,WB_REWRITE_MANY)
      !---------------------------------------------------------------------
      return
   end subroutine drv_time_increment


   !/@*
   subroutine drv_time_info(F_stepno,F_is_chkpt_L,F_is_last_L,F_is_stat_L,F_datev,F_dateo,F_dt, &
        F_nb_time_levels,F_swap_mode)
      implicit none
      !@objective Increase time step and update GMM fields names
      !@arguments
      integer,intent(out),optional :: F_stepno       !- current step
      logical,intent(out),optional :: F_is_chkpt_L  !- .T. if it's a checkpoint step
      logical,intent(out),optional :: F_is_last_L   !- .T. if it's the last step [restart of run end]
      logical,intent(out),optional :: F_is_stat_L   !- .T. if it's time to do print stats
      integer,intent(out),optional :: F_datev !- Actual valid date, CMC datetime stamp
      integer,intent(out),optional :: F_dateo !- Date of origine, CMC datetime stamp
      real,intent(out),optional :: F_dt !- Time step length [sec]
      integer,intent(out),optional :: F_nb_time_levels
      integer,intent(out),optional :: F_swap_mode
      !@author
   !*@/
      real(REAL64) :: nhours_8
      !---------------------------------------------------------------------
      if (.not.m_init_L) then
         call msg(MSG_ERROR,'(drv_time) Not Initialised')
         return
      endif

      if (present(F_stepno)) F_stepno = time_stepno
      if (present(F_is_chkpt_L)) F_is_chkpt_L = &
           (time_stepno_chkpt >= 0 .and. &
           (mod(time_stepno,time_stepno_chkpt) == 0 .or. &
           mod(time_stepno,time_stepno_rstrt) == 0))
      if (present(F_is_last_L)) F_is_last_L = (time_stepno >= stepno_last)
      if (present(F_is_stat_L)) F_is_stat_L = &
           (time_stepno_stat >= 0 .and. &
           (mod(time_stepno,time_stepno_stat) == 0))
      if (present(F_dateo)) F_dateo = time_dateo
      if (present(F_datev)) then
         nhours_8 = (dble(time_stepno) * dble(time_dt))/3600.d0
         call incdatr(F_datev,time_dateo,nhours_8)
      endif
      if (present(F_dt)) F_dt = time_dt
      if (present(F_nb_time_levels)) F_nb_time_levels = time_nb_levels
      if (present(F_swap_mode)) F_swap_mode = m_swap_mode
      !---------------------------------------------------------------------
      return
   end subroutine drv_time_info


   !/@*
   subroutine drv_time_shuffle_name(F_label_S,F_timeFlags_S)
      implicit none
      !@objective Time Shuffle a GMM label time flag
      !@arguments
      character(len=*),intent(inout) :: F_label_S
      character(len=*),intent(in) :: F_timeFlags_S(:)
      !@author
      !  Michel Desgagne, Feb 2008
      !  Ron McTaggart-Cowan, Feb 2008
      !  Stephane Chamberland, Feb 2008
   !*@/
      integer :: i,nlvl,i2,i3
      character(len=GMM_MAXNAMELENGTH):: timeflag
      !---------------------------------------------------------------------
      i = index(F_label_S,':')
      if (i == 0) return !- pattern not found, no timeflag, no possible shuffle

      timeflag = F_label_S(i+1:)
 
      nlvl = size(F_timeFlags_S)
      i3 = -1
      do i2 = 1,nlvl
         if (timeflag == F_timeFlags_S(i2)) then
            i3 = i2+1
            if (i3 > nlvl) i3 = 1
            F_label_S = F_label_S(1:i)//F_timeFlags_S(i3)
            exit
         endif
      enddo
!!$      if (i3 < 0) error
      !---------------------------------------------------------------------
      return
   end subroutine drv_time_shuffle_name


   !/@*
   function drv_time_shuffle_list(F_namelist_S,F_timeFlags_S) result(F_nflds)
      implicit none
      !@objective Return list of base var name to shuffle (w/o timeflag)
      !@arguments
      character(len=*),intent(out) :: F_namelist_S(:)
      character(len=*),intent(in) :: F_timeFlags_S(:)
      !@return
      integer :: F_nflds
      !@author
      !  Michel Desgagne, Feb 2008
      !  Ron McTaggart-Cowan, Feb 2008
      !  Stephane Chamberland, Feb 2008
   !*@/
      character(len=GMM_MAXNAMELENGTH) :: name1,name2,timeflag
      character(len=GMM_MAXNAMELENGTH) :: keys(MAX_FIELDS),timeflags0(8)
      integer :: istat,nkeys,ii,jj,nlvl,nflds
      !---------------------------------------------------------------------
      F_nflds = 0

      !# Note: GMM is case sensitive, use all lowercase
      nlvl = size(F_timeFlags_S)
      timeflags0(1:nlvl) = F_timeFlags_S(1:nlvl)
      do ii=1,nlvl
         istat = clib_tolower(timeflags0(ii))
      enddo

      !# Get list of fields matching timeflags
      nkeys = gmm_keys(keys)
      do ii = 1,nkeys
         istat = clib_tolower(keys(ii))
         jj = index(keys(ii),':')
         if (jj == 0) continue !- no timeflag, ignore
         name1 = keys(ii)
         timeflag = name1(jj+1:)
         if (any(timeflags0(1:nlvl) == timeflag)) then
            name2 = name1(1:jj)
            if (.not.any(F_namelist_S(1:F_nflds) == name2)) then
               F_nflds = min(F_nflds + 1,size(F_namelist_S))
               F_namelist_S(F_nflds) = name2
            endif
         endif
      end do

      !# Keep only fields matching all time flags
      nflds = 0
      do ii = 1,F_nflds
         istat = 0
         do jj=1,nlvl
            if (any(keys(:) == trim(F_namelist_S(ii))//trim(timeflags0(jj)))) &
                 istat=istat+1
         enddo
         if (istat == nlvl) then
            nflds = nflds + 1
            F_namelist_S(nflds) = F_namelist_S(ii)
         endif
      enddo
      F_nflds = nflds
      !---------------------------------------------------------------------
      return
   end function drv_time_shuffle_list


   !/@*
   subroutine drv_time_shuffle(F_timeFlags_S,F_swap_mode)
      implicit none
      !@objective Time Shuffle lables in GMM
      !@arguments
      character(len=*),intent(in) :: F_timeFlags_S(:)
      integer,intent(in),optional :: F_swap_mode
      !@author
      !  Michel Desgagne, Feb 2008
      !  Ron McTaggart-Cowan, Feb 2008
      !  Stephane Chamberland, Feb 2008
   !*@/
      character(len=GMM_MAXNAMELENGTH) :: namelist(TIME_FLAGS_MAX)
      character(len=GMM_MAXNAMELENGTH) :: flds(MAX_FIELDS),timeflags0(8)
      integer :: ii,jj,istat,nflds,nlvl,swap_mode
      !---------------------------------------------------------------------
      if (.not.m_init_L) then
         call msg(MSG_ERROR,'(drv_time) Not Initialised')
         return
      endif

      swap_mode = m_swap_mode
      if (present(F_swap_mode)) swap_mode = F_swap_mode

      nlvl = size(F_timeFlags_S)
      timeflags0(1:nlvl) = F_timeFlags_S(1:nlvl)
      do ii=1,nlvl
         istat = clib_tolower(timeflags0(ii))
      enddo

      nflds = drv_time_shuffle_list(flds,timeflags0(1:nlvl))

      do ii = 1,nflds
         namelist(1) = trim(flds(ii))//trim(timeflags0(1))
         do jj = 2,nlvl
            namelist(jj) = namelist(jj-1)
            call drv_time_shuffle_name(namelist(jj),timeflags0(1:nlvl))
         enddo
         if (nlvl == 2) then
            call msg(MSG_INFOPLUS,'(drv_time) shuffle: '//trim(namelist(1))//' <-> '//trim(namelist(2)))
         else if (nlvl == 3) then
            call msg(MSG_INFOPLUS,'(drv_time) shuffle: '//trim(namelist(1))//' <- '//trim(namelist(2))//' <- '//trim(namelist(3)))
         else
            call msg(MSG_INFOPLUS,'(drv_time) shuffle: '//trim(namelist(1))//' <-> ... <->'//trim(namelist(nlvl)))
         endif
         if (swap_mode == DRV_TIME_MODE_SHUFFLE) then
            istat = gmm_shuffle(namelist(1:nlvl))
         else
            istat = priv_shuffle_copy(namelist(1:nlvl))
         endif
      end do
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_ERROR,'(drv_time) Problem doing shuffle')
      endif
      !---------------------------------------------------------------------
      return
   end subroutine drv_time_shuffle


   !/@*
   function priv_shuffle_copy(F_varlist_S) result(F_istat)
      implicit none
      !@objective Initialize time config
      !@arguments
      character(len=*),intent(in) :: F_varlist_S(:)
      !@return
      integer :: F_istat
      !*@/
      integer :: ivar
      real,pointer :: ptr2d1(:,:),ptr3d1(:,:,:),ptr2d2(:,:),ptr3d2(:,:,:)
      type(gmm_metadata) :: mymeta,mymeta2
      !---------------------------------------------------------------------
      F_istat = RMN_OK
      do ivar = 1,size(F_varlist_S)-1
         F_istat = gmm_getmeta(F_varlist_S(ivar),mymeta)
         F_istat = min(gmm_getmeta(F_varlist_S(ivar+1),mymeta2),F_istat)
         if (F_istat /= 0) then
            call msg(MSG_ERROR,'(drv_time) Cannot copy fields, not same size/rank')
            return
         endif
         if (mymeta%l(3)%n /= mymeta2%l(3)%n) then
            call msg(MSG_ERROR,'(drv_time) Cannot copy fields, not same size/rank')
            F_istat = RMN_ERR
            return
         endif
         if (mymeta%l(3)%n == 0) then
            nullify(ptr2d1,ptr2d2)
            F_istat = gmm_get(F_varlist_S(ivar),ptr2d1)
            F_istat = min(gmm_get(F_varlist_S(ivar+1),ptr2d2),F_istat)
            if (.not.(RMN_IS_OK(F_istat).and.associated(ptr2d1).and.associated(ptr2d2))) then
               call msg(MSG_ERROR,'(drv_time) Cannot copy fields, Problem getting data pointers')
               F_istat = RMN_ERR
               return
            endif
            if (any(lbound(ptr2d1) /= lbound(ptr2d2)) .or. &
                 any(ubound(ptr2d1) /= ubound(ptr2d2))) then
               call msg(MSG_ERROR,'(drv_time) Cannot copy fields, size mismatch')
               F_istat = RMN_ERR
               return
            endif
            ptr2d1(:,:) = ptr2d2(:,:)
         else
            nullify(ptr3d1,ptr3d2)
            F_istat = gmm_get(F_varlist_S(ivar),ptr3d1)
            F_istat = min(gmm_get(F_varlist_S(ivar+1),ptr3d2),F_istat)
            if (.not.(RMN_IS_OK(F_istat).and.associated(ptr3d1).and.associated(ptr3d2))) then
               call msg(MSG_ERROR,'(drv_time) Cannot copy fields, Problem getting data pointers')
               F_istat = RMN_ERR
               return
            endif
            if (any(lbound(ptr3d1) /= lbound(ptr3d2)) .or. &
                 any(ubound(ptr3d1) /= ubound(ptr3d2))) then
               call msg(MSG_ERROR,'(drv_time) Cannot copy fields, size mismatch')
               F_istat = RMN_ERR
               return
            endif
            ptr3d1(:,:,:) = ptr3d2(:,:,:)
         endif
      enddo
      !---------------------------------------------------------------------
      return
   end function priv_shuffle_copy

end module drv_time_mod
