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
      subroutine INs_wait (F_dateV_S,F_client)
      use, intrinsic :: iso_fortran_env
      use iso_c_binding
      use omp_timing
      use IOs
      use INs
      implicit none
      
      character(len=16), intent(OUT) :: F_dateV_S
      integer, intent(OUT) :: F_client

      character(len=16) :: msg
      integer :: ierr,tag,i,indx1,indx2,cnt
!     
!--------------------------------------------------------------------
!
      if (allocated(INs_list_S)) deallocate (INs_list_S)
      INs_nreq= -1
      
    !  call clock ( Lun_out, 'Blocking rcv1', .false. )
      if (INs_1o1_L) then
         tag=801
         call MPI_recv ( INs_nreq, 1, MPI_INTEGER,INs_gem1o1,&
                         tag,INs_GEM_COMM,MPI_STATUSES_IGNORE,ierr)
         if (INs_nreq>0) then
            allocate (INs_list_S(INs_nreq))
            call MPI_recv (INs_list_S, INs_nreq*len(INs_list_S),&
                           MPI_CHARACTER,INs_gem1o1,tag+1      ,&
                            INs_GEM_COMM,MPI_STATUSES_IGNORE,ierr)
         endif
      endif
     ! call clock ( Lun_out, 'Blocking rcv1 ... DONE', .false. )
      call MPI_bcast (INs_nreq,1,MPI_INTEGER,0,&
                      MY_WORLD_COMM,ierr)
      if (INs_nreq<1) return
                      
      if (.not.INs_1o1_L) allocate (INs_list_S(INs_nreq))
      call MPI_bcast (INs_list_S,INs_nreq*len(INs_list_S),&
                      MPI_CHARACTER,0,MY_WORLD_COMM,ierr)

! Cracking the request list
      
      indx1 = index(INs_list_S(1),":")
      indx2 = index(INs_list_S(1)(indx1+1:),":")+indx1
      if ((indx1==0).or.(indx2==0)) then
         print*, 'Invalid FORMAT in request list item 1',&
                 '- will do nothing'
         INs_nreq= 0
         return
      endif
      F_dateV_S= INs_list_S(1)(indx1+1:indx2-1)
      read(INs_list_S(1)(indx2+1:),'(i4)') Inp_rtag
      F_client= Inp_rtag
      if (F_client==1001) msg='NESTING'
      if (F_client==2001) msg='IAU'
      if (F_client==3001) msg='SPN'

      if (Lun_out > 0) write(lun_out,9000) trim(msg), trim(F_dateV_S),&
                                           Inp_rtag

! Client list of requests in INs_list_S is meant to be very flexible
! while the ubound(SRL) is set to INs_maxreqs in INs_init.
! Client must provide the proper maximum at init stage.
      cnt=0
      do i=2,INs_nreq
         indx1 = index(INs_list_S(i),":")
         indx2 = index(INs_list_S(i)(indx1+1:),":")+indx1  
         if ((indx1==0).or.(indx2==0)) then
            print*, 'Invalid FORMAT in request list item: ',&
                     i, ' Ignoring item'
            cycle
         endif
         cnt=cnt+1
         if (cnt>INs_maxreqs) then
            print*, 'INs_wait: NOT enough storage in SRL - ABORT'
            stop
         endif
         SRL(cnt)%vname(1) =INs_list_S(i)(1:indx1-1)
         SRL(cnt)%stag     =INs_list_S(i)(indx1+1:indx2-1)
         SRL(cnt)%src      =INs_list_S(i)(indx2+1:)
      end do
      INs_nreq= cnt
      deallocate (INs_list_S)

 9000 format(/,' TREATING ',a,' INPUT DATA VALID AT: ',a,' on tag: ',i4,&
      /,' ===============================================')
!     
!--------------------------------------------------------------------
!
      return
      end subroutine INs_wait  
