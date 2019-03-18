!/* RPN_COMM - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2012  Division de Recherche en Prevision Numerique
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
!**function RPN_COMM_topo_2 compute the local bounds for a distributed array
!InTf!
      integer function RPN_COMM_topo_2(nxg,minx,maxx,nxl,nxlmax,halo,nx0,alongx,fill,relax,abort) !InTf!
      use rpn_comm
      implicit none                                                         !InTf!
      integer, intent(in) :: nxg,halo,relax                                 !InTf!
      logical, intent(in) :: alongx,fill,abort                              !InTf!
      integer, intent(out):: minx,maxx,nxl,nxlmax,nx0                       !InTf!
!
!arguments
!  I	nxg	global dimension of data (along selected X or Y axis)
!  O	minx,maxx
!		dimensionsing bounds that should be used for local array (minx:maxx)
!               "tile" own data will be (1:nxl-nx0+1)
!  O    nxl     end of the current "tile" in global space
!  O    nxlmax  dimension of the largest "tile"
!  I    halo    size of the halo area (along  selected X or Y axis)
!  O    nx0     start of the current "tile" in global space
!  I    alongx .true. if distributing along X axis, .false. if distributing along Y axis
!  I    fill   .true. if dimension (usually along X) should be optimized for current machine
!  I    relax  currently 0(strict) 1, 2, or 3. this value is passed to function RPN_COMM_limit_2 that
!              contains the actual splitting algorithm. see that function documentation for the
!              actual meaning of this flag.
!  I    abort  .true.  : abort run if split is unsuccessful
!              .false. : return error code if unsuccessful
!notes
!       this function will return 0 in case of success and will either terminate execution in case of failure.
!       or return the error code from RPN_COMM_limit_2 depending on the value of the abort argument
!       the number of PEs along the selected axis as well as the ordinal along that axis of this PE are
!       obtained from module rpn_comm
!*
!
        integer gmin, gmax, ierr
        integer count(pe_nx + pe_ny)
        integer depl(pe_nx + pe_ny)
      external RPN_COMM_limit_2
      integer RPN_COMM_limit_2
!
        if (alongx) then    ! distribute along X axis
          ierr = RPN_COMM_limit_2(pe_mex, pe_nx, 1, nxg, gmin, gmax, count, depl,relax)
        else                ! distribute along Y axis
          ierr = RPN_COMM_limit_2(pe_mey, pe_ny, 1, nxg, gmin, gmax, count, depl,relax)
        end if
!
      if(ierr.ne.0) then
         write(rpn_u,*) 'RPN_COMM_topo: invalid distribution, ABORT'
           if(abort) then  ! this may have brutal consequences
            call RPN_COMM_finalize(ierr)
            stop
           else
              RPN_COMM_topo_2 = ierr    ! error code from RPN_COMM_limit_2
              return
           endif
      endif

        nx0 = gmin              ! "tile" covers from nx0 to nxl in global space
        minx = 1 - halo         ! "tile" dimension along this axis should be minx:maxx
        nxl = gmax - gmin + 1
        nxlmax = count(1)       ! first "tile" has the largest dimension (but may not be the only one)
        maxx = nxlmax + halo
        if(fill) maxx = maxx + 1 - mod(nxlmax,2)   ! this needs to be revised for cache based machines
!
        RPN_COMM_topo_2 = 0     ! SUCCESS
        return
        end  function RPN_COMM_topo_2                     !InTf!
!    this is the old function, it calls the newer RPN_COMM_topo_2 forcing 
!    the strict distribution mode used previously and abort in case of error
!     kept for compatibility with older versions of this library
!InTf!
!!integer function RPN_COMM_topo(nxg,minx,maxx,nxl,nxlmax,halo,nx0,alongx,fill) !InTf!
      integer function RPN_COMM_topo(nxg,minx,maxx,nxl,nxlmax,halo,nx0,alongx,fill)
! Generate needed information about local tile along a specified axis.
! The input is the total number of point to divide and the size of the halo.
! The function will split the domain depending of the topology
! induced by the PEs.
      implicit none                                                                !InTf!
      integer, intent(in) :: nxg,halo                                              !InTf!
      logical, intent(in) :: alongx,fill                                           !InTf!
      integer, intent(out):: minx,maxx,nxl,nxlmax,nx0                              !InTf!
        external RPN_COMM_topo_2
        integer RPN_COMM_topo_2
        RPN_COMM_topo=RPN_COMM_topo_2(nxg,minx,maxx,nxl,nxlmax,halo,nx0,alongx,fill,0,.true.)
        return
        end function RPN_COMM_topo                                                 !InTf!
