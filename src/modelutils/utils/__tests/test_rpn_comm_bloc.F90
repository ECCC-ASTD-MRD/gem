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

subroutine test_rpn_comm_bloc()
   use, intrinsic :: iso_fortran_env, only: INT64
   use iso_c_binding
   implicit none
#include <clib_interface_mu.hf>
   include "rpn_comm.inc"
   external :: test_rpn_comm_bloc_idoms,test_rpn_comm_bloc_p0
   integer :: mydomain,myproc,numproc,npex,npey,ngrids,igrid,nblocx,nblocy
   integer :: istat,bloc_npe,bloc_ipe,bloc_ipe2,bloc_ipex,bloc_ipey
   character(len=80) :: tmp_S
   ! ---------------------------------------------------------------------
   call rpn_comm_mydomain(test_rpn_comm_bloc_idoms,mydomain)

   ngrids = 1
   npex = 0
   npey = 0
   igrid = rpn_comm_init_multigrid(test_rpn_comm_bloc_p0,myproc,numproc,npex,npey,ngrids)

   nblocx = npex
   nblocy = npey
   istat = clib_getenv(trim('MPI_NBLOCX'),tmp_S)
   if (istat >= 0) then
      read(tmp_S,fmt=*,iostat=istat) nblocx
      if (istat /= 0) nblocx=npex
   endif
   istat = clib_getenv(trim('MPI_NBLOCY'),tmp_S)
   if (istat >= 0) then
      read(tmp_S,fmt=*,iostat=istat) nblocy
      if (istat /= 0) nblocy=npey
   endif
   
   print *,'rpn_comm_bloc(',nblocx,nblocy,')'
   istat = rpn_comm_bloc(nblocx,nblocy)

   call rpn_comm_size(RPN_COMM_BLOC_COMM,bloc_npe,istat)
   call rpn_comm_rank(RPN_COMM_BLOC_COMM,bloc_ipe,istat)
   call rpn_comm_bloctopo(bloc_ipe2,bloc_ipex,bloc_ipey,nblocx,nblocy)
   print *,bloc_npe,':',nblocx,nblocy,':',bloc_ipe,'=',bloc_ipe2,':',bloc_ipex,bloc_ipey
!!$P=2x1 B=2x1
!!$0-1: 1 : 2 1 : 0 = 0 : 0 0
!!$P=2x1 B=1x1
!!$0: 2 : 1 1 : 0 = 0 : 0 0
!!$1: 2 : 1 1 : 1 = 1 : 1 0
!!$P=1x2 B=1x2
!!$0-1: 1 : 1 2 : 0 = 0 : 0 0
!!$P=1x2 B=1x1
!!$0: 2 : 1 1 : 0 = 0 : 0 0
!!$1: 2 : 1 1 : 1 = 1 : 0 1

   call rpn_comm_finalize(istat)
   ! ---------------------------------------------------------------------
   return
end subroutine test_rpn_comm_bloc


subroutine test_rpn_comm_bloc_p0(F_npex,F_npey)
   implicit none
   integer,intent(out) :: F_npex,F_npey
#include <clib_interface_mu.hf>
   integer :: istat
   character(len=80) :: tmp_S
   ! ---------------------------------------------------------------------
   F_npex = 1
   F_npey = 1
   istat = clib_getenv(trim('MPI_NPEX'),tmp_S)
   if (istat >= 0) then
      read(tmp_S,fmt=*,iostat=istat) F_npex
      if (istat /= 0) F_npex = 1
   endif
   istat = clib_getenv(trim('MPI_NPEY'),tmp_S)
   if (istat >= 0) then
      read(tmp_S,fmt=*,iostat=istat) F_npey
      if (istat /= 0) F_npey = 1
   endif
   ! ---------------------------------------------------------------------
   return
end subroutine test_rpn_comm_bloc_p0


subroutine test_rpn_comm_bloc_idoms(F_ndomains,F_offset,F_istat)
   implicit none
   integer, intent(out) :: F_ndomains,F_offset,F_istat
   ! ---------------------------------------------------------------------
   F_istat = 0
   F_ndomains = 1
   F_offset = 0
   ! ---------------------------------------------------------------------
   return
end subroutine test_rpn_comm_bloc_idoms
