!/* RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! !                          Environnement Canada
! *
! ! This library is free software; you can redistribute it and/or
! ! modify it under the terms of the GNU Lesser General Public
! ! License as published by the Free Software Foundation,
! ! version 2.1 of the License.
! *
! ! This library is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! ! Lesser General Public License for more details.
! *
! ! You should have received a copy of the GNU Lesser General Public
! ! License along with this library; if not, write to the
! ! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! ! Boston, MA 02111-1307, USA.
! */
!
!InTf!
      subroutine RPN_COMM_set_petopo(sizx,sizy)           !InTf!
!     set PE block topology
!     assemble PEs in blocks of sizx by sizy tiles
!
!     default behavior is 1 by 1 
!     (along X first then along Y, see horizontal fill example)
!
!     PEs within a block : along X first, then along Y
!     blocks on nodes : X direction first, then Y direction
!     (a block is sizx PEs by sizy PEs)
!
! examples: 
!          8 by 8 model grid, sizx=4, sizy=4,  (block fill)
!          16 processes per node, 64 processes (00-63)
!
!  pe_mey        (node 2)               (node 3)
!         +----+----+----+----+  +----+----+----+----+
!    7    | 44 | 45 | 46 | 47 |  | 60 | 61 | 62 | 63 |
!         +----+----+----+----+  +----+----+----+----+
!    6    | 40 | 41 | 42 | 43 |  | 56 | 57 | 58 | 59 |
!         +----+----+----+----+  +----+----+----+----+
!    5    | 36 | 37 | 38 | 39 |  | 52 | 53 | 54 | 55 |
!         +----+----+----+----+  +----+----+----+----+
!    4    | 32 | 33 | 34 | 35 |  | 48 | 49 | 50 | 51 |
!         +----+----+----+----+  +----+----+----+----+
!
!                 (node 0)              (node 1)
!         +----+----+----+----+  +----+----+----+----+
!    3    | 12 | 13 | 14 | 15 |  | 28 | 29 | 30 | 31 |
!         +----+----+----+----+  +----+----+----+----+
!    2    | 08 | 09 | 10 | 11 |  | 24 | 25 | 26 | 27 |
!         +----+----+----+----+  +----+----+----+----+
!    1    | 04 | 05 | 06 | 07 |  | 20 | 21 | 22 | 23 |
!         +----+----+----+----+  +----+----+----+----+
!    0    | 00 | 01 | 02 | 03 |  | 16 | 17 | 18 | 19 |
!         +----+----+----+----+  +----+----+----+----+
!            0    1    2    3    4    5    6    7    pe_mex
!   
!          8 by 8 model grid, sizx=2, sizy=4,  (block fill)
!          16 processes per node, 64 processes (00-63)
!
!  pe_mey        (node 2)               (node 3)
!         +----+----+----+----+  +----+----+----+----+
!    7    | 38 | 39 | 46 | 47 |  | 54 | 55 | 62 | 63 |
!         +----+----+----+----+  +----+----+----+----+
!    6    | 36 | 37 | 44 | 45 |  | 52 | 53 | 60 | 61 |
!         +----+----+----+----+  +----+----+----+----+
!    5    | 34 | 35 | 42 | 43 |  | 50 | 51 | 58 | 59 |
!         +----+----+----+----+  +----+----+----+----+
!    4    | 32 | 33 | 40 | 41 |  | 48 | 49 | 56 | 57 |
!         +----+----+----+----+  +----+----+----+----+
!
!                 (node 0)              (node 1)
!         +----+----+----+----+  +----+----+----+----+
!    3    | 06 | 07 | 14 | 15 |  | 22 | 23 | 30 | 31 |
!         +----+----+----+----+  +----+----+----+----+
!    2    | 04 | 05 | 12 | 13 |  | 20 | 21 | 28 | 29 |
!         +----+----+----+----+  +----+----+----+----+
!    1    | 02 | 03 | 10 | 11 |  | 18 | 19 | 26 | 27 |
!         +----+----+----+----+  +----+----+----+----+
!    0    | 00 | 01 | 08 | 09 |  | 16 | 17 | 24 | 25 |
!         +----+----+----+----+  +----+----+----+----+
!            0    1    2    3    4    5    6    7    pe_mex
!   
!          8 by 8 model grid, sizx=2, sizy=2,  (block fill)
!          16 processes per node, 64 processes (00-63)
!
!  pe_mey
!         +----+----+----+----+----+----+----+----+
!    7    | 50 | 51 | 54 | 55 | 58 | 59 | 62 | 63 |
!         +----+----+----+----+----+----+----+----+ (node 3)
!    6    | 48 | 49 | 52 | 53 | 56 | 57 | 60 | 61 |
!         +----+----+----+----+----+----+----+----+
!
!         +----+----+----+----+----+----+----+----+
!    5    | 34 | 35 | 38 | 39 | 41 | 43 | 46 | 47 |
!         +----+----+----+----+----+----+----+----+ (node 2)
!    4    | 32 | 33 | 36 | 37 | 40 | 42 | 44 | 45 |
!         +----+----+----+----+----+----+----+----+
!
!         +----+----+----+----+----+----+----+----+
!    3    | 18 | 19 | 22 | 23 | 26 | 27 | 30 | 31 |
!         +----+----+----+----+----+----+----+----+ (node 1)
!    2    | 16 | 17 | 20 | 21 | 24 | 25 | 28 | 29 |
!         +----+----+----+----+----+----+----+----+
!
!         +----+----+----+----+----+----+----+----+
!    1    | 02 | 03 | 06 | 07 | 10 | 11 | 14 | 15 |
!         +----+----+----+----+----+----+----+----+ (node 0)
!    0    | 00 | 01 | 04 | 05 | 08 | 09 | 12 | 13 |
!         +----+----+----+----++----+----+----+----+
!            0    1    2    3    4    5    6    7    pe_mex
!   
!          4 by 4 model grid, sizx=4, sizy=1, (horizontal fill)
!          4 processes per node, 16 processes (00-15)
!
!  pe_mey 
!         +----+----+----+----+
!    3    | 12 | 13 | 14 | 15 | (node 3)
!         +----+----+----+----+
!         +----+----+----+----+
!    2    | 08 | 09 | 10 | 11 | (node 2)
!         +----+----+----+----+
!         +----+----+----+----+
!    1    | 04 | 05 | 06 | 07 | (node 1)
!         +----+----+----+----+
!         +----+----+----+----+
!    0    | 00 | 01 | 02 | 03 | (node 0)
!         +----+----+----+----+ 
!            0    1    2    3     pe_mex
!   
!          4 by 4 model grid, sizx=1, sizy=4, (vertical fill)
!          4 processes per node, 16 processes (00-15)
!
!           (0)     (1)     (2)     (3)   (node) 
!         +----+  +----+  +----+  +----+
!    3    | 03 |  | 07 |  | 11 |  | 15 |
!         +----+  +----+  +----+  +----+
!    2    | 02 |  | 06 |  | 10 |  | 14 |
!         +----+  +----+  +----+  +----+
!    1    | 01 |  | 05 |  | 09 |  | 13 |
!         +----+  +----+  +----+  +----+
!    0    | 00 |  | 04 |  | 08 |  | 12 |
!  pe_mey +----+  +----+  +----+  +----+ 
!            0       1       2       3    pe_mex

      use rpn_comm
      implicit none                      !InTf!
!      include 'mpif.h'
      integer, intent(IN) :: sizx,sizy   !InTf!

      deltai = sizx
      deltaj = sizy

      return
      end subroutine RPN_COMM_set_petopo                  !InTf!
!InTf!
      integer function RPN_COMM_petopo(pex,pey)          !InTf!
      use rpn_comm
      implicit none                                      !InTf!
!      include 'mpif.h'
      integer, intent(IN) :: pex,pey                     !InTf!

      integer count, ierr,i,j,i0,j0

      pe_nx=pex
      pe_ny=pey
      deltai = min(pe_nx,max(1,deltai))
      deltaj = min(pe_ny,max(1,deltaj))

      RPN_COMM_petopo = -1

      if(pex*pey .ne. pe_dommtot) then
        write(rpn_u,*) pe_me
        write(rpn_u,*) &
     &           'RPN_COMM_petopo: invalid distribution, pe_nx*pe_ny'
        write(rpn_u,*) '                 is not equal to pe_dommtot ',&
     &            pex*pey,pe_dommtot
        return
      endif

      count = pe_pe0
      
      if(allocated(pe_id))      deallocate(pe_id)
      if(allocated(pe_xtab))    deallocate(pe_xtab)
      if(allocated(pe_ytab))    deallocate(pe_ytab)
      if(allocated(ord_tab))    deallocate(ord_tab)

      allocate(pe_id(-1:pe_nx,-1:pe_ny))
      allocate(pe_xtab(0:pe_tot),pe_ytab(0:pe_tot))
      allocate(ord_tab(0:pe_tot))

      do i=0,pe_tot
         pe_xtab(i)=-1
         pe_ytab(i)=-1
         ord_tab(i)=-1
      enddo
      ord_max=0

!     deltai = deltaj = 1          : as before, X increasing, then Y increasing (horizontal stripes, bottom to top)
!     deltai = 1, deltaj = pe_ny   : Y increasing, then X increasing (vertical stripes, left ro right)
!     otherwise                    : blocks of deltai by deltaj, X increasing, then Y increasing 
!                                    blocks going left to right, then bottom to top
!     in the latter case, deltai * deltaj = number of PEs per 
      do j0 = 0 , pe_ny-1 , deltaj   ! distribute by blocks of size deltai by deltaj
      do i0 = 0 , pe_nx-1 , deltai
         do j=j0,min(pe_ny-1 , j0+deltaj-1)
         do i=i0,min(pe_nx-1 , i0+deltai-1)
            if(count == pe_me) then  ! this is me, get pe_mex and pe_mey
               pe_mex = i
               pe_mey = j
            endif
            pe_id(i,j)=count-pe_pe0
            pe_xtab(count)=i
            pe_ytab(count)=j
            ord_tab(ord_max)=pe_pe0+ord_max
            ord_max=ord_max+1
            count=count+1
         enddo
         enddo
      enddo
      enddo
      ord_max = ord_max -1
!     
!     fill processor topology matrix, including a border used
!     to find "neighbors" in case of periodic domain
!     
      do j=-1,pe_ny
         pe_id(-1,j) = pe_id(pe_nx-1,j)
         pe_id(pe_nx,j) = pe_id(0,j)
      enddo
      do i=-1,pe_nx
         pe_id(i,-1) = pe_id(i,pe_ny-1)
         pe_id(i,pe_ny) =  pe_id(i,0)
      enddo
      if(pe_me.eq.pe_pe0 .and. diag_mode .ge.2)then
         print *,'PE MATRIX :'
         do j=pe_ny,-1,-1
	    print 101,(pe_id(i,j),i=-1,pe_nx)
 101        format(30I4)
         enddo
         print *,'PE_xtab :'
         print 101,(pe_xtab(i),i=0,pe_tot-1)
         print *,'PE_ytab :'
         print 101,(pe_ytab(i),i=0,pe_tot-1)
         print *,'ordinals table'
         print 101,(ord_tab(i),i=0,ord_max)
      endif
!     
!     compute position of PE in the computational grid
!     a value of -1 for X or Y means "OUT OF GRID"
!     

!      pe_mex=mod(pe_me-pe_pe0,pe_nx)  ! this is no longer true
!      pe_mey=(pe_me-pe_pe0)/pe_nx     ! this is no longer true
      pe_extra=0
!     
!     are we on a boundary ?
!     
      bnd_south=pe_mey.eq.0
      bnd_north=pe_mey.eq.pe_ny-1
      bnd_west=pe_mex.eq.0
      bnd_east=pe_mex.eq.pe_nx-1    
!     
!     split communication domain into rows
!     
      pe_myrow=MPI_COMM_NULL
!
      call MPI_COMM_SPLIT(pe_indomm,pe_mey+1,pe_mex+1,pe_myrow,ierr)	     
!
      call MPI_COMM_rank(pe_myrow,i,ierr)
      pe_mycol=MPI_COMM_NULL
      call MPI_COMM_SPLIT(pe_indomm,pe_mex+1,pe_mey+1,pe_mycol,ierr)
      RPN_COMM_petopo = 0
      return
      end function RPN_COMM_petopo                              !InTf!
!
! ========================================================================================
!
!InTf!
      integer function RPN_COMM_get_pe(x,y,grd,sgrd,communicator) !InTf!
!
!     get PE ordinal in grid/supergrid/domain
!     x,y  : coordinates of PE in its own grid ( 0 -> pe_nx-1 , 0 -> pe_ny-1)
!     grd  : grid number the PE belongs to (-1 = this PE's grid)
!     sgrd : supergrid number the PE belongs to (-1 = this PE's supergrid)
!     communicator : may be 'GRID', 'SUPERGRID', or 'ALLGRIDS'
!          communicator                value returned by function
!          'GRID'      : coordinates of PE x,y in current GRID (grd, sgrd must be -1)
!          'SUPERGRID' : coordinates of PE x,y in current SUPERGRID (sgrd must be -1)
!          'ALLGRIDS'  : coordinates of PE x,y in current domain (set of SUPERGRIDs)
!
      use rpn_comm
      implicit none                                                !InTf!
      integer, intent(IN) :: x,y,grd,sgrd                          !InTf!
      character (len=*), intent(IN) :: communicator                   !InTf!

      integer ordinal

      RPN_COMM_get_pe = -1
      if(x >= pe_nx .or. y >= pe_ny .or. x < 0 .or. y < 0)  return    ! bad PE in grid coordinates
      if(grd < -1 .or. grd > pe_tot_multi_grid/pe_tot_grid-1) return  ! bad grid number
      if(sgrd < -1 .or. sgrd > pe_tot_a_domain/pe_tot_multi_grid) &
     &                                                      return    ! bad multigrid number

      if(trim(communicator) == 'GRID') then        ! PE ordinal in own grid
         ordinal = pe_id(x,y)
      endif
      if(trim(communicator) == 'MULTIGRID') then   ! PE ordinal in own supergrid
         if(grd == -1) then
            ordinal = pe_id(x,y) + (my_colors(3))*pe_tot_grid   ! add own grid offset if grd = -1
         else
            ordinal = pe_id(x,y) + (grd-1)*pe_tot_grid          ! add grid offset
         endif
      endif
      if(trim(communicator) == 'ALLGRIDS') then     ! PE ordinal in own domain
         ordinal = pe_id(x,y)                               ! PE ordinal in own grid
         if(grd == -1) then
            ordinal = ordinal + my_colors(3)*pe_tot_grid    ! add own grid offset if grd = -1
         else
            ordinal = ordinal + (grd-1)*pe_tot_grid         ! add grid offset
         endif
         if(sgrd == -1) then
            ordinal = ordinal + my_colors(2)*pe_tot_grid    ! add own supergrid offset if grd = -1
         else
            ordinal = ordinal + (sgrd-1)*pe_tot_multi_grid  ! add supergrid offset
         endif
      endif

      RPN_COMM_get_pe = ordinal
      return
      end function RPN_COMM_get_pe                              !InTf!
