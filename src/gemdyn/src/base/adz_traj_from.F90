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

      subroutine adz_traj_from (F_reqs,F_nr,F_export,F_nptr,F_dimr)
      use ISO_C_BINDING
      use, intrinsic :: iso_fortran_env
      use adz_mem
      use ptopo
      implicit none

      integer, intent(OUT) :: F_nr
      integer, intent(IN ) :: F_nptr, F_dimr
      real, dimension(F_dimr) , intent(OUT) :: F_reqs
      type(ADZ_SLOD), intent(IN) :: F_export

      integer :: i,cnt,err,nreq,err_get
      integer, dimension(4,0:Ptopo_world_numproc) :: lis_req
!
!     ---------------------------------------------------------------
!
      ORIGIN_DATATYPE = REAL_DATATYPE
      TARGET_DATATYPE = REAL_DATATYPE

      call MPI_Win_fence(0, F_export%winreqs, err)
      call rpn_comm_barrier ("MULTIGRID", err)
      call MPI_Win_fence(0, F_export%wintraj, err)

! get requested positions from others
      cnt=0 ; F_nr=0 ; nreq= 0 ; lis_req(4,0)= 0; err_get= 0
      do i=1,2*Ptopo_world_numproc,2
         if (F_export%from(i) > 0) then
            cnt=cnt+1
            ORIGIN_COUNT = F_export%from(i) * 3
            TARGET_COUNT = ORIGIN_COUNT
            TARGET_RANK = i/2
            TARGET_DISP = F_export%from(i+1)
            nreq= nreq+1
            lis_req(1,nreq) = TARGET_RANK
            lis_req(2,nreq) = F_export%from(i)
            lis_req(3,nreq) = lis_req(4,nreq-1)
            lis_req(4,nreq) = lis_req(4,nreq-1) + &
                              F_nptr*F_export%from(i)
            if (cnt+TARGET_COUNT-1>F_dimr) then
               err_get=-1
            else
            call MPI_Get(F_reqs(cnt), ORIGIN_COUNT, ORIGIN_DATATYPE,&
                              TARGET_RANK,TARGET_DISP, TARGET_COUNT,&
                              TARGET_DATATYPE, F_export%wintraj, err)
            F_nr = F_nr + F_export%from(i)
            endif
            cnt= cnt + TARGET_COUNT - 1
         endif
      end do
      call MPI_Win_fence(0, F_export%wintraj, err)

! clip the requests
      do i=1, F_nr*3, 3
         F_reqs(i  )= min(max(dble(F_reqs(i  )),&
                      Adz_iminposx),Adz_imaxposx)
         F_reqs(i+1)= min(max(dble(F_reqs(i+1)),&
                      Adz_iminposy),Adz_imaxposy)
      end do
                   
! send adresses of interpolated requests to tell neighbors
! where to fetch results within Adz_cor
                              
      Adz_offs= -999
      ORIGIN_COUNT = 2
      TARGET_COUNT = ORIGIN_COUNT
      ORIGIN_DATATYPE = INTEGER_DATATYPE
      TARGET_DATATYPE = INTEGER_DATATYPE
      
      call MPI_Win_fence(0, Adz_offs_win, err)
      call gem_error(err_get,'Adz_traj_from',&
           'NOT ENNOUGH MEMORY - Adz_MAX_MPI_OS_SIZE for F_reqs')
      
      do i=1, nreq
         TARGET_RANK = lis_req(1,i)
         TARGET_DISP = (Ptopo_world_myproc*2)
         call MPI_Put( lis_req(2,i), ORIGIN_COUNT, ORIGIN_DATATYPE,&
                               TARGET_RANK,TARGET_DISP, TARGET_COUNT,&
                                    TARGET_DATATYPE,Adz_offs_win, err)
      end do
!     
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_traj_from
