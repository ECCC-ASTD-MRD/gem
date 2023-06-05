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
subroutine test_inbloc_mpi()
   use testutils
   implicit none
   !@objective Test that rpn_comm blocking can be changed back and forth
   !@author Stephane Chamberland, 2011-08
   !@/
#include <rmnlib_basics.hf>
#include <rpn_comm.hf>
   integer,parameter :: NSTEPS = 1024
   integer :: npex,npey,myproc,medomm,mex,mey,nblocx,nblocy,ismaster,mymaster,   &
        mybloc,myblocx,myblocy,blocme,istat,blocrank0,blocrank,nblocx2,nblocy2,ninblocx,ninblocy,istep
   character(len=12) :: domname_S
   character(len=256) :: msg_S
   ! ---------------------------------------------------------------------
   myproc = testutils_initmpi()

   !- get original bloc topology
   call rpn_comm_carac(npex,npey,myproc, &
        medomm,mex,mey,nblocx,nblocy,ismaster, mymaster,   &
        mybloc, myblocx,myblocy,blocme,domname_S)
   call rpn_comm_rank(RPN_COMM_BLOC_COMM,blocrank0,istat)
!!$   write(msg_S,*) myproc,blocrank0,':',npex,npey,':',nblocx,nblocy
!!$   call msg(MSG_INFO,msg_S)

   do istep=1,NSTEPS
      print *,'step=',istep ; call flush(6)

   !- reset to inbloc topology
   ninblocx = npex
   ninblocy = npey
   istat = rpn_comm_bloc(ninblocx,ninblocy)

   !- check
   call rpn_comm_carac(npex,npey,myproc, &
        medomm,mex,mey,nblocx2,nblocy2,ismaster, mymaster,   &
        mybloc, myblocx,myblocy,blocme,domname_S)
   call testutils_assert_ok(nblocx2==ninblocx.and.nblocy2==ninblocy,'test_inbloc_mpi','new nblocxy')

   call rpn_comm_rank(RPN_COMM_BLOC_COMM,blocrank,istat)
   call testutils_assert_ok(blocrank==0,'test_inbloc_mpi','new nblocxy rank')

   !- return to original bloc topo
   istat = rpn_comm_bloc(nblocx,nblocy)

   !- check
   call rpn_comm_carac(npex,npey,myproc, &
        medomm,mex,mey,nblocx2,nblocy2,ismaster, mymaster,   &
        mybloc, myblocx,myblocy,blocme,domname_S)
   call testutils_assert_ok(nblocx2==nblocx.and.nblocy2==nblocy,'test_inbloc_mpi','reset nblocxy')

   call rpn_comm_rank(RPN_COMM_BLOC_COMM,blocrank,istat)
   call testutils_assert_ok(blocrank==blocrank0,'test_inbloc_mpi','reset nblocxy rank')

enddo

   call rpn_comm_finalize(istat)
   ! ---------------------------------------------------------------------
   return
end subroutine test_inbloc_mpi
