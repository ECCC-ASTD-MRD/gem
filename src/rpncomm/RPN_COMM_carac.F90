!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! !                          Environnement Canada
! !
! ! This library is free software; you can redistribute it and/or
! ! modify it under the terms of the GNU Lesser General Public
! ! License as published by the Free Software Foundation,
! ! version 2.1 of the License.
! !
! ! This library is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! ! Lesser General Public License for more details.
! !
! ! You should have received a copy of the GNU Lesser General Public
! ! License along with this library; if not, write to the
! ! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! ! Boston, MA 02111-1307, USA.
! !/
!InTf!
!!subroutine RPN_COMM_carac(npex,npey,me,medomm,mex,mey,sizex,sizey,ismaster, mymaster, mybloc, myblocx,myblocy,blocme,domname) !InTf!
      subroutine RPN_COMM_carac(npex,npey,me,medomm,mex,mey,sizex,sizey,&
          ismaster, mymaster, mybloc, myblocx,myblocy,blocme,domname)
! npex,npey: Number of PE along x and y
! me       : Rank in the context of communicator "ALL"
! medomm   : Rank in the context of communicator "DOMM"
! mex,mey  : x and y coordinates in domain
! sizex,sizey: Size of blocks along x and y axis
! ismaster : Equals 1 if PE is included in communicator "BLOCMASTER", 0 else
! mymaster : Rank (relative to the domain) of the master of the PE's block
! mybloc   : Rank of PE's block
! myblocx,myblocy: Coordinate of PE's block along x and y
! blocme   : Rank of PE relative to its block
! domname  : Domain name of the current PE
      use rpn_comm
      implicit none                                                             !InTf!
      integer, intent(out) :: npex,npey,me,mex,mey,sizex,sizey,ismaster         !InTf!
      integer, intent(out) :: mymaster, mybloc, myblocx,myblocy,blocme          !InTf!
      integer, intent(out) :: medomm                                            !InTf!
      character(len=*), intent(out) :: domname                                  !InTf!
!arguments
!     I nblocx, nblocy: number of blocks on the subgrid in x-y direction
!     O RPN_COMM_bloc : error status (-1 if error, else 0)
!      include 'rpn_comm.h'
!      include 'mpif.h'
!
      npex     = pe_nx
      npey     = pe_ny
      me       = pe_me
      medomm   = pe_medomm
      mex      = pe_mex
      mey      = pe_mey
      sizex    = BLOC_sizex
      sizey    = BLOC_sizey
      ismaster = BLOC_master
      mymaster = BLOC_corner
      mybloc   = BLOC_mybloc
      myblocx  = BLOC_myblocx
      myblocy  = BLOC_myblocy
      blocme   = BLOC_me
      domname  = 'DOM1'
      return
      end subroutine RPN_COMM_carac                                              !InTf!
