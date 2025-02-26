
module series_mod
   use series_options, only: series_paused_L, series_fileid, PKGNAME_S
   use series_nml_mod, only: series_nml, SER_NML_ERR, SER_NML_NONE, SER_NML_OK
   use series_init_mod, only: series_init
   use series_step_mod, only: series_stepinit, series_stepend
   use series_xst_mod, only: series_xst, series_isstep, series_isvar
   implicit none
   private

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   public :: series_nml, series_init, series_terminate, series_stepinit, &
        series_stepend, series_pause, series_resume, series_isstep, &
        series_isvar, series_xst
   public :: SER_NML_ERR, SER_NML_NONE, SER_NML_OK

contains

   !/@*
   function series_pause() result(F_istat)
      implicit none
      !@returns
      integer :: F_istat
      !*@/
      !---------------------------------------------------------------
      F_istat = RMN_OK
      series_paused_L = .true.

      !#TODO: Append doesn't seem to work, only close at the end
!!$      if (series_fileid /= -1) then
!!$         call msg(MSG_INFO, PKGNAME_S//'Close time series file')
!!$         istat = fclos(series_fileid)
!!$         series_fileid = -1
!!$      endif
      !---------------------------------------------------------------
      return
   end function series_pause


   !/@*
   function series_resume() result(F_istat)
      implicit none
      !@returns
      integer :: F_istat
      !*@/
      !---------------------------------------------------------------
      F_istat = RMN_OK
      series_paused_L = .false.
      !---------------------------------------------------------------
      return
   end function series_resume
 

   !/@*
   function series_terminate() result(F_istat)
      implicit none
      !@returns
      integer :: F_istat
      !*@/
      logical, parameter :: FORCE_L = .true.
      integer :: istat
      !---------------------------------------------------------------
      !#TODO: protect from multi thread
      !#TODO: adjust for not writing at every extraction
!!$      if (mod(series_kount, series_out_nsteps) /= 0) &
!!$           call series_write(FORCE_L)

      if (series_fileid /= -1) then
         call msg(MSG_INFO, PKGNAME_S//'Close time series file')
         istat = fclos(series_fileid)
         series_fileid = -1
      endif

      F_istat = RMN_OK
      !---------------------------------------------------------------
      return
   end function series_terminate

end module series_mod
