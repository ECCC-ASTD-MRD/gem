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
      integer function RPN_COMM_topo_xy(nig,njg,mini,maxi,minj,maxj,nil,njl,haloi,haloj,peri,perj) !InTf!
      use rpn_comm
      implicit none                                               !InTf!
      integer, intent(IN) :: nig,njg,haloi,haloj                  !InTf!
      integer, intent(OUT) :: nil,njl,mini,maxi,minj,maxj         !InTf!
      logical, intent(IN) ::  peri,perj                           !InTf!
!
!	include 'rpn_comm.h'
!	include 'mpif.h'
!
      mini = 1-haloi
      minj = 1-haloj
      nil = (nig + pe_nx - 1)/pe_nx
      njl = (njg + pe_ny - 1)/pe_ny
!      if (pe_me .eq.0) then
!	  print *,' RPN_COMM_topo: nig,njg=',nig,njg
!	  print *,' RPN_COMM_topo: nil,njl=',nil,njl
!	  print *,' RPN_COMM_topo: pe_nx,pe_ny=',pe_nx,pe_ny
!      endif
      maxi = nil + haloi + 1 - mod(nil,2)
      maxj = njl + haloj + 1 - mod(njl,2)
      nil = min(nil,nig-nil*(pe_mex))
      njl = min(njl,njg-njl*(pe_mey))
      RPN_COMM_topo_xy = MPI_SUCCESS
      if(nil.le.0 .or. njl.le.0) RPN_COMM_topo_xy = MPI_ERROR
      return
      end function RPN_COMM_topo_xy !InTf!
