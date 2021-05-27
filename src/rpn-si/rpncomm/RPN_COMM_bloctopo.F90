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
      subroutine RPN_COMM_bloctopo(blocme,blocmex,blocmey,blocsizex,blocsizey)   !InTf!
	use rpn_comm
      implicit none                                                              !InTf!
      integer, intent(out) :: blocme,blocmex,blocmey                             !InTf!
      integer, intent(out) :: blocsizex, blocsizey                               !InTf!
!      include 'mpif.h'
!

      blocsizex    = BLOC_sizex
      blocsizey    = BLOC_sizey
      blocme   = BLOC_me
      blocmex = pe_mex-BLOC_myblocx*(pe_nx/BLOC_sizex)
      blocmey = pe_mey-BLOC_myblocy*(pe_ny/BLOC_sizey)
      
      return
      end subroutine RPN_COMM_bloctopo                                            !InTf!
