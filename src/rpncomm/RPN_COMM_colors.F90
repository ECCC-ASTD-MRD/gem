!/* RPN_COMM - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
!InTf!
      integer function RPN_COMM_colors(comm)         !InTf!
      use rpn_comm
      implicit none                                  !InTf!
      character(len=*) :: comm                       !InTf!
      RPN_COMM_colors=-1
      if (comm(1:9).eq.'ALLGRIDS') RPN_COMM_colors=my_colors(1)
      if (comm(1:5).eq.'WORLD')    RPN_COMM_colors=my_colors(1)
      if (comm(1:9).eq.'MULTIGRID')RPN_COMM_colors=my_colors(2)
      if (comm(1:4).eq.'GRID')     RPN_COMM_colors=my_colors(3)
      return
      end function RPN_COMM_colors                    !InTf!
