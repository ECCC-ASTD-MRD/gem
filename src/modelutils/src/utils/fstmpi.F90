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
module fstmpi_mod
   use iso_c_binding
   use rpn_comm_itf_mod
   use fst_mod
   use fstmpi_read_mod
   use fstmpi_write_mod
   use ptopo_utils
   implicit none
   private
   !@objective 
   !@author  Stephane Chamberland, 2011-06
   !@description
   ! Public functions
   public :: fstmpi_open, fstmpi_close, fstmpi_find, fstmpi_read, fstmpi_rdhint, &
        fstmpi_rdhint_3d_r4, fstmpi_rdhint_3d_r4_vect, fstmpi_write, &
        fstmpi_getmeta, fstmpi_get_gridid, fstmpi_get_hgridid, fstmpi_get_vgrid

   ! Public constants
   public :: FST_READONLY,FST_FIND_LT,FST_FIND_LE,FST_FIND_NEAR,FST_FIND_GE,FST_FIND_GT,FST_NPAK_DEFAULT,FST_NPAK_FULL32
!@/

#include <rmn/msg.h>
#include <rmnlib_basics.hf>
!!!#include <arch_specific.hf>

contains

   !/@
   function fstmpi_open(F_filename_S,F_readonly_L,F_dir_ok_L) result(F_fileid)
      implicit none
      !@objective Open rpn std file
      !@arguments
      character(len=*), intent(in) :: F_filename_S
      logical, intent(in),optional :: F_readonly_L
      logical, intent(in),optional :: F_dir_ok_L
      !@author
      !@return
      integer :: F_fileid
      !@/
      logical :: readonly_L,dir_ok_L
      integer :: istat,mysize,mydata(1)
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(fstmpi) open [BEGIN]')
      call ptopo_init_var()
      F_fileid = RMN_ERR
      if (ptopo_isblocmaster_L) then
         readonly_L = .not.FST_READONLY
         dir_ok_L = .false.
         if (present(F_readonly_L)) readonly_L = F_readonly_L
         if (present(F_dir_ok_L)) dir_ok_L = F_dir_ok_L
         F_fileid = fst_open(F_filename_S,readonly_L,dir_ok_L)
      endif
      mydata(1) = F_fileid
      mysize = size(mydata)
      call rpn_comm_bcast(mydata, mysize, RPN_COMM_INTEGER, RPN_COMM_MASTER, RPN_COMM_BLOC_COMM,istat)
      F_fileid = mydata(1)
      call msg(MSG_DEBUG,'(fstmpi) open [END]')
      ! ---------------------------------------------------------------------
      return
  end function fstmpi_open


   !/@
   function fstmpi_close(F_fileid) result(F_istat)
      implicit none
      !@objective Close rpn std file
      !@arguments
      integer, intent(in) :: F_fileid
      !@author
      !@return
      integer :: F_istat
      !@/
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(fstmpi) close [BEGIN]')
      call ptopo_init_var()
      F_istat = RMN_OK
      if (ptopo_isblocmaster_L) then
         F_istat = fst_close(F_fileid)
      endif
      call msg(MSG_DEBUG,'(fstmpi) close [END]')
      ! ---------------------------------------------------------------------
      return
   end function fstmpi_close


end module fstmpi_mod
