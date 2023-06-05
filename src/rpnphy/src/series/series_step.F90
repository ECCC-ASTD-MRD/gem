!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

module series_step_mod
   use phygridmap, only: phydim_ni, phydim_nj
   use series_options
   use series_geop_mod, only: series_geop
   use series_write_mod, only: series_write_geo, series_write
   implicit none
   private

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   public :: series_stepinit, series_stepend

contains

   !/@*
   function series_stepinit(F_stepno) result(F_istat)
      implicit none
      !@object 
      !@arguments
      integer, intent(in) :: F_stepno
      !@returns
      integer :: F_istat
      !@author Andre Methot - aug 94 - v0_14
      !@revision
      ! v2_00 - Desgagne M.       - initial MPI version
      ! v3_30 - Winger K.         - correct time series handling in climate mode
      ! 2017-01, S.Chamberland    - Major revision
      !*@/
      ! ---------------------------------------------------------------
      !#TODO: protect from multi thread
      F_istat = RMN_OK
      if (.not.series_on_L) return

      F_istat = RMN_ERR

      series_end_L = (F_stepno <= P_serg_serstp)
      series_kount = F_stepno
      series_heure_8 = dble(F_stepno) * dble(series_delt)/3600.d0

      F_istat = RMN_OK
      ! ---------------------------------------------------------------
      return
   end function series_stepinit


   !/@*
   function series_stepend() result(F_istat)
      implicit none
      !@Object Write time series buffers
      !@arguments
      !@returns
      integer :: F_istat
      !*@/
      logical, parameter :: FORCE_L = .true.
      ! ---------------------------------------------------------------
      F_istat = RMN_OK
      if (.not.series_on_L) return

      !#TODO: protect from multi thread
      F_istat = RMN_ERR

      IF_KOUNT0: if (series_kount == 0) then
         call series_geop()
         call series_write_geo()
      else
         call series_write(.not.FORCE_L)
      endif IF_KOUNT0

      !#TODO: Append doesn't seem to work, only close at the end
!!$      if (series_fileid /= -1) then
!!$         call msg(MSG_INFO, PKGNAME_S//'Close time series file')
!!$         istat = fclos(series_fileid)
!!$         series_fileid = -1
!!$      endif

      F_istat = RMN_OK
      ! ---------------------------------------------------------------
      return
   end function series_stepend


   !#---- Private functions ------------------------------------------

end module series_step_mod
