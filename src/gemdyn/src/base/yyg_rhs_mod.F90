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

module yyg_rhs
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save
!______________________________________________________________________
!                                                                      |
!  VARIABLES ASSOCIATED WITH Ru,Rv Solver for Yin-Yang communication   |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! Rhsx_recv_len(M)   | Number of values to receive from each PE for    |
!                    | its West Pilot area  (M=Ptopo_numproc)          |
! Rhsx_recv_i(*,M)   | local gridpoint I to receive value from PE(*)   |
! Rhsx_recv_j(*,M)   | local gridpoint J to receive value from PE(*)   |
! Rhsx_send_len(*)   | Number of values to send to  PE (*) for West    |
! Rhsx_send_imx(*,M) | closest I gridpoint on the other panel to find  |
!                    | the value for Rhsx_sendw_xxr,Rhsx_sendw_yyr     |
! Rhsx_send_imy(*,M) | closest J gridpoint on the other panel to find  |
!                    | the value for  Rhsx_sendw_xxr,Rhsx_sendw_yyr    |
! Rhsx_send_xxr(*,M) | longitude in the other panel to find the value  |
!                    | for receiving panel                             |
! Rhsx_send_yyr(*,M) | latitude in the other panel to find the value   |
! Rhsx_send_sten(*,M)|position J or I to find the stencil             |
!______________________________________________________________________|
!
!Declarations for Ru,Rv variables (on Phi grid)

   integer  :: Rhsx_send_all,Rhsx_recv_all,Rhsx_sendmaxproc,Rhsx_recvmaxproc
   integer, dimension (:), allocatable :: &
               Rhsx_sendproc,  Rhsx_recvproc, &
               Rhsx_recv_len, Rhsx_send_len, &
               Rhsx_recv_adr, Rhsx_send_adr, &
               Rhsx_recv_i   ,Rhsx_recv_j   ,Rhsx_send_imx,Rhsx_send_imy
   real(kind=REAL64),  dimension (:), allocatable :: &
               Rhsx_send_sten, Rhsx_send_xxr, Rhsx_send_yyr

end module yyg_rhs
