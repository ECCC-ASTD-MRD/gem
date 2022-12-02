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

!/@
module fstmpio_mod
   use iso_c_binding
   use rpn_comm_itf_mod
   use fst_mod
   use fstmpio_rdhint_mod
!!$   use fstmpio_write_mod
   use ptopo_utils
   implicit none
   private
   !@objective
   !@author  Stephane Chamberland, 2017-04
   !@description
   ! Public functions
   public :: fstmpio_set_iotype, fstmpio_open, fstmpio_close, fstmpio_find, &
        fstmpio_getmeta, fstmpio_get_vgrid, fstmpio_rdhint
   public :: fstmpio_find_0, fstmpio_find_vect, fstmpio_find_3d_0, fstmpio_find_3d_vect
   public :: fstmpio_rdhint_3d_r4, fstmpio_rdhint_3d_r4_vect, fstmpio_get_hgridid

   !#TODO: , fstmpio_write

   ! Public constants
   public :: FST_READONLY, FST_FIND_LT, FST_FIND_LE, FST_FIND_NEAR, FST_FIND_GE, &
        FST_FIND_GT, FST_NPAK_DEFAULT, FST_NPAK_FULL32, FST_FIND_DIAG_T, FST_FIND_DIAG_M
!@/

#include <rmn/msg.h>
#include <rmnlib_basics.hf>
!!!#include <arch_specific.hf>

contains

   !/@*
   subroutine fstmpio_set_iotype(F_iotype)
      implicit none
      !@objective Set the type of rpn_comm communicator used in ptopo for io
      !           (PTOPO_BLOC or PTOPO_IODIST)
      !@arguments
      integer, intent(in) :: F_iotype
      !*@/
      !----------------------------------------------------------------------
      ptopo_iotype = F_iotype
      !----------------------------------------------------------------------
      return
   end subroutine fstmpio_set_iotype


   !/@
   function fstmpio_open(F_filename_S, F_readonly_L, F_dir_ok_L) result(F_fileid)
      implicit none
      !@objective Open rpn std file
      !@arguments
      character(len=*), intent(in) :: F_filename_S
      logical, intent(in), optional :: F_readonly_L
      logical, intent(in), optional :: F_dir_ok_L
      !@author
      !@return
      integer :: F_fileid
      !@/
      logical :: readonly_L, dir_ok_L, isiomaster_L, isiope_L
      integer :: istat, mysize, mydata(1), comm_ipe_io_master
      character(len=64) :: communicator_S
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG, '(fstmpio) open [BEGIN]')

      istat = ptopo_get_io_params(isiomaster_L, isiope_L, comm_ipe_io_master, communicator_S)

      F_fileid = RMN_ERR
      if (isiope_L) then
         readonly_L = .not.FST_READONLY
         dir_ok_L = .false.
         if (present(F_readonly_L)) readonly_L = F_readonly_L
         if (present(F_dir_ok_L)) dir_ok_L = F_dir_ok_L
         F_fileid = fst_open(F_filename_S, readonly_L, dir_ok_L)
      endif
      mydata(1) = F_fileid
      mysize = size(mydata)
      !#TODO: reduce min on io PEs to check if all io PEs open are ok
      call rpn_comm_bcast(mydata, mysize, RPN_COMM_INTEGER, comm_ipe_io_master, &
           trim(communicator_S), istat)
      if (.not.isiope_L) F_fileid = mydata(1)
      call msg(MSG_DEBUG, '(fstmpio) open [END]')
      ! ---------------------------------------------------------------------
      return
   end function fstmpio_open


   !/@
   function fstmpio_close(F_fileid) result(F_istat)
      implicit none
      !@objective Close rpn std file
      !@arguments
      integer, intent(in) :: F_fileid
      !@author
      !@return
      integer :: F_istat
      !@/
      logical :: isiomaster_L, isiope_L
      integer :: istat, comm_ipe_io_master
      character(len=64) :: communicator_S
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG, '(fstmpio) close [BEGIN]')

      istat = ptopo_get_io_params(isiomaster_L, isiope_L, comm_ipe_io_master, communicator_S)

      F_istat = RMN_OK
      if (isiope_L) then
         F_istat = fst_close(F_fileid)
      endif
      call msg(MSG_DEBUG, '(fstmpio) close [END]')
      ! ---------------------------------------------------------------------
      return
   end function fstmpio_close


end module fstmpio_mod
