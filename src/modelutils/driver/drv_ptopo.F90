!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!-------------------------------------------------------------------------- 

!/@*
module drv_ptopo_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use iso_c_binding
   use rpn_comm_itf_mod
   use clib_itf_mod, only: clib_chdir, clib_getenv
   use wb_itf_mod
   use ptopo_utils
   use drv_path_mod
   implicit none
   private
   !@objective Initialization of execution grid and processor topology
   !@author
   !  Michel Desgagne, Feb 2008
   !  Ron McTaggart-Cowan, Feb 2008
   !  Stephane Chamberland, Feb 2008
   !@revisions
   !  2008-09, Stephane Chamberland: allow multi grid topology
   !  2010-06, Stephane Chamberland: split in smaller pieces
   !  2012-02, Stephane Chamberland: update for maestro/taks_setup,multi domains
   !@public_functions
   public :: drv_ptopo_init,drv_ptopo_terminate
   !@description
   !  Initialization of execution grid and processor topology
   !
   !  VARIABLES ASSOCIATED WITH LOGICAL PARALLEL PROCESSOR TOPOLOGY       |
   !                                                                      |
   !     along Y                                                          |
   !        .           .                    .                            |
   !        .           .                    .                            |
   !   +-----------+-----------+     +---------------+                    |
   !   | (0,myrow) | (1,myrow) |.....| (mycol,myrow) |.....               |
   !   +-----------+-----------+     +---------------+                    |
   !        .           .                    .                            |
   !        .           .                    .                            |
   !   +-----------+-----------+     +---------------+                    |
   !   |   (0,2)   |   (1,2)   |.....|   (mycol,2)   |.....               |
   !   +-----------+-----------+     +---------------+                    |
   !   |   (0,1)   |   (1,1)   |.....|   (mycol,1)   |.....               |
   !   +-----------+-----------+     +---------------+                    |
   !   |   (0,0)   |   (1,0)   |.....|   (mycol,0)   |..... along X       |
   !   +-----------+-----------+     +---------------+                    |
!*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>

   !- privates vars and parameters
   character(len=*),parameter :: status_file_name_S = 'status_drv.dot'
   character(len=*),parameter :: PROCESS_S = 'process'


contains


   !/@*
    function drv_ptopo_init(F_ngrids) result(F_istat)
      implicit none
      !@objective Initialize processor topology
      !@argument
      integer,intent(in),optional :: F_ngrids
      !@returns
      integer :: F_istat !- Error Status
      !@author
      !  Michel Desgagne, Feb 2008
      !  Ron McTaggart-Cowan, Feb 2008
      !  Stephane Chamberland, Feb 2008
      !@revisions
      !  2012-02, Stephane Chamberland: update for maestro/taks_setup,multi domains
   !*@/
      logical,parameter :: DO_BCAST_L = .true.

      logical,save :: ptopo_is_init_L = .false.

      integer :: istat,ngrids,npex,npey,igrid,myproc,numproc,nblocx,nblocy,mydomain,ndomains
      character(len=RMN_PATH_LEN) :: tmp_S

      integer,external :: rpn_comm_init_multi_level
      external :: drv_ptopo_ndoms,drv_ptopo_p0 !- call back s/r to define: ndomains, mydomain, npex,npey
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'[BEGIN] drv_ptopo_init')
      F_istat = RMN_ERR
      if (ptopo_is_init_L) return
      F_istat = RMN_OK

      !- RPN_COMM [MPI] and OpenMP Initialization
      call rpn_comm_mydomain(drv_ptopo_ndoms, mydomain) !- UM_EXEC_NDOMAINS
      ndomains = 1
      istat = wb_get('ptopo_cfgs/ndomains',ndomains)
!!$      call rpn_comm_bcast(ndomains,1,RPN_COMM_INTEGER,RPN_COMM_MASTER, &
!!$           RPN_COMM_GRID,istat) !TODO-later: is this needed? is drv_ptopo_ndoms called on all pe?
      ngrids = 1
      if (present(F_ngrids)) ngrids = F_ngrids
      istat = wb_put('ptopo_cfgs/ngrids',ngrids)

      F_istat = drv_path_set(mydomain)
      if (.not.RMN_IS_OK(F_istat)) return

      npex = 0
      npey = 0
      igrid = rpn_comm_init_multi_level(drv_ptopo_p0,myproc,numproc,npex,npey,ndomains,ngrids)
!!$      igrid = rpn_comm_init_multigrid(drv_ptopo_p0,myproc,numproc,npex,npey,ngrids)
      F_istat = min(F_istat,igrid)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_ERROR,'(drv_ptopo_init) Problem in rpn_comm_init_multi_level')
         return
      endif
 
      call ptopo_init_var(ndomains,mydomain,ngrids,igrid)
 
      if (ptopo_grid_ipe == RPN_COMM_MASTER) then
         !NOTE: wb_broadcast broken, need to get from pe0 then bcast
         istat = wb_get('ptopo_cfgs/nblocx',nblocx)
         istat = wb_get('ptopo_cfgs/nblocy',nblocy)
      endif
      istat = ptopo_bloc_set(nblocx,nblocy,DO_BCAST_L)
      istat = ptopo_openMP_set()

      istat = drv_path_set_basedir()
      istat = drv_path_set(mydomain,ngrids,igrid,ptopo_grid_ipex,ptopo_grid_ipey)
      istat = clib_chdir(trim(drv_path_work_S))

      write(tmp_S,'(a,4i6)') '(drv_ptopo_init) idomain,igrid,ipex,ipey:',&
           ptopo_world_idom,ptopo_dom_igrid,ptopo_grid_ipex,ptopo_grid_ipey
      call msg(MSG_INFO,tmp_S)

      ptopo_is_init_L = (RMN_IS_OK(F_istat))
      call msg(MSG_DEBUG,'[END] drv_ptopo_init')
      !---------------------------------------------------------------------
      return
   end function drv_ptopo_init


!!$   !/@*
!!$   function ptopo_setdir() result(F_istat)
!!$      implicit none
!!$      !@objective Create process dirs
!!$      !@returns
!!$      integer :: F_istat !- Error Status
!!$      !@author
!!$      !  Michel Desgagne, Feb 2008
!!$      !  Ron McTaggart-Cowan, Feb 2008
!!$      !  Stephane Chamberland, Feb 2008
!!$   !*@/
!!$      integer :: istat,iun
!!$      character(len=RMN_PATH_LEN) :: localdir
!!$      !---------------------------------------------------------------------
!!$      F_istat = clib_getcwd(ptopo_basedir)
!!$      write(localdir,'(A10,I2.2,A1,I2.2,A1,I2.2)') './'//trim(PROCESS_S)//'/',ptopo_mycol,'-',ptopo_myrow,'-',ptopo_igrid
!!$      if (.not.RMN_IS_OK(clib_isdir(trim(localdir)))) then
!!$          istat = clib_mkdir(trim(PROCESS_S))
!!$          istat = clib_mkdir(trim(localdir))
!!$      endif
!!$
!!$      !- Create a status file
!!$      if (ptopo_myproc == RPN_COMM_MASTER) then
!!$         iun = 0
!!$         istat = fnom(iun, status_file_name_S, 'FMT',0)
!!$         write(iun,'(a)') '_status_drv=ABORT;'
!!$         call flush(iun)
!!$         istat = fclos(iun)
!!$      endif
!!$
!!$      F_istat = min(clib_chdir(trim(localdir)),F_istat)
!!$      !---------------------------------------------------------------------
!!$      return
!!$   end function ptopo_setdir


   !/@*
   subroutine drv_ptopo_terminate(F_status_S)
      implicit none
      !@objective Finalize... the MPI way
      !@arguments
      character(len=*), intent(in), optional :: F_status_S
      !@author
      !  Michel Desgagne, Feb 2008
      !  Ron McTaggart-Cowan, Feb 2008
      !  Stephane Chamberland, Feb 2008
   !*@/
      logical,parameter :: IS_BEGIN_L = .true.
      character(len=64) :: status_S,model_name_S
      integer :: istat
      !---------------------------------------------------------------------
      status_S = 'ED'
      if (present(F_status_S)) status_S = F_status_S
      call write_status_file3('_status='//trim(status_S))
      call close_status_file3()
  
      call drv_print_banner(.not.IS_BEGIN_L)

      istat = wb_get('ATM_MODEL_NAME',model_name_S)
      if (.not.RMN_IS_OK(istat)) &
           istat = clib_getenv('ATM_MODEL_NAME',model_name_S)
      if (.not.RMN_IS_OK(istat)) model_name_S = '(unknown name)'
      call stop_mpi(RMN_OK,trim(model_name_S),'Normal Ending')
!!$      call rpn_comm_barrier(RPN_COMM_WORLD, istat)
!!$      call rpn_comm_finalize(istat)
      !---------------------------------------------------------------------
      return
   end subroutine drv_ptopo_terminate


end module drv_ptopo_mod
