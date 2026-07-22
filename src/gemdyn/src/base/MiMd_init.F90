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

      subroutine MiMd_init (F_component_S, F_color, F_nc, F_COMM,&
                            F_cnum, F_cme)
      use ISO_C_BINDING
      use MiMd
      implicit none
      
      character(len=*), intent(IN) :: F_component_S
      integer, intent(IN ) :: F_nc
      integer, intent(IN ) :: F_color(F_nc)
      integer, intent(OUT) :: F_COMM(F_nc), F_cnum, F_cme
      
      include 'mpif.h'

      integer :: i,j,k,ierr,color,numproc,myproc,required,provided
      integer :: lcl_comm,nb_color,me_color,me(4)
      integer :: indx_color,ordinal_in_set,comm
      integer, dimension(:  ), allocatable :: inter_connect
      integer, dimension(:,:), allocatable :: inter_color_comm
!     
!--------------------------------------------------------------------
!
!MPI_THREAD_SINGLE: Only one thread will execute. 
!MPI_THREAD_FUNNELED: The process may be multi-threaded, but only the main thread will make MPI calls (all MPI calls are funneled to the main thread). 
!MPI_THREAD_SERIALIZED: The process may be multi-threaded, and multiple threads may make MPI calls, but only one at a time: MPI calls are not made concurrently from two distinct threads (all MPI calls are serialized). 
!MPI_THREAD_MULTIPLE: Multiple threads may call MPI, with no restrictions.
!      required = MPI_THREAD_MULTIPLE
      required = MPI_THREAD_SINGLE
      call MPI_Init_thread (required, provided, ierr)
      if (provided /= required ) then
         if (MiMd_Wmyproc==0) write (6,'(/3x,a/)')&
         'FAILED in MPI_Init_thread: your system does NOT'//&
         'support MPI_THREAD_MULTIPLE -ABORT-'
         call MPI_finalize (ierr)
         stop
      endif

      call gem_init_appmpi (F_COMM)

      MiMd_gemworld= F_COMM(1)

      call MPI_COMM_size (MiMd_gemworld,MiMd_Wnumproc,ierr)
      call MPI_COMM_rank (MiMd_gemworld,MiMd_Wmyproc ,ierr)
      
      call MPI_Comm_split (MiMd_gemworld, F_color(1), MiMd_Wmyproc, F_COMM(1), ierr)
      call MPI_barrier(MiMd_gemworld,ierr)

      allocate (partners(4,MiMd_Wnumproc),names_S(MiMd_Wnumproc))
      
      call MPI_COMM_size (F_COMM(1),numproc,ierr)
      call MPI_COMM_rank (F_COMM(1),myproc ,ierr)

      F_cnum= numproc ; F_cme= myproc
      color= MPI_UNDEFINED
      if (myproc==0) then
         color=1
      endif
      call MPI_Comm_split (MiMd_gemworld, color, MiMd_Wmyproc, lcl_comm, ierr)
      me_color = -99
      if(lcl_comm .ne. MPI_COMM_NULL) then 
         call MPI_COMM_size (lcl_comm,nb_color,ierr)
         call MPI_COMM_rank (lcl_comm,me_color,ierr)
      endif
      
      me(1) = F_color(1)
      me(2) = MiMd_Wmyproc
      me(3) = myproc
      me(4) = numproc
      call MPI_barrier(MiMd_gemworld,ierr)
      call MPI_Allgather(F_component_S,256,MPI_CHARACTER,names_S,256,MPI_CHARACTER,MiMd_gemworld,ierr)
      call MPI_Allgather(me,4,MPI_INTEGER,partners,4,MPI_INTEGER,MiMd_gemworld,ierr)
      color=-1 ;  MiMd_ncolors=0
      do i=1, MiMd_Wnumproc
         if ( partners(1,i)/=color) then
            MiMd_ncolors= MiMd_ncolors+1
            MiMd_world(MiMd_ncolors)%color  = partners(1,i)
            MiMd_world(MiMd_ncolors)%rank   = partners(4,i)
            MiMd_world(MiMd_ncolors)%name_S = names_S (  i)
            if (MiMd_world(MiMd_ncolors)%color==F_color(1)) indx_color=MiMd_ncolors
            color= partners(1,i)
         end if
      end do
      do j=1, MiMd_ncolors
         do i=1, MiMd_Wnumproc
            if ( partners(1,i)==MiMd_world(j)%color) then
               if ( partners(3,i)==0) MiMd_world(j)%wpe0=partners(2,i)
            endif
         enddo
      enddo
      allocate (inter_color_comm(MiMd_ncolors,MiMd_ncolors));inter_color_comm=-99
      if (myproc==0) then
         allocate (inter_connect(MiMd_ncolors)) ; inter_connect=-1
         do i=2,F_nc
            do j=1, MiMd_ncolors
               if (MiMd_world(j)%color == F_color(i)) inter_connect(j)=F_color(i)
            end do
         end do
         call MPI_Allgather(inter_connect,MiMd_ncolors,MPI_INTEGER,inter_color_comm,MiMd_ncolors,MPI_INTEGER,lcl_comm,ierr)
      endif
      call MPI_bcast (inter_color_comm, size(inter_color_comm), MPI_INTEGER,0,F_COMM(1),ierr)

      F_COMM(2:)= MPI_COMM_NULL
      do j=1, MiMd_ncolors-1
         do i=1+j, MiMd_ncolors
            if (inter_color_comm(i,j)>=0) then
               color= MPI_UNDEFINED ; ordinal_in_set= -1
               if ((j==indx_color).or.(i==indx_color)) then
                  color=1
                  ordinal_in_set=MiMd_Wmyproc
               endif
               call MPI_Comm_split (MiMd_gemworld, color, ordinal_in_set, comm, ierr)
               call MPI_barrier(MiMd_gemworld,ierr)
               if (comm/=MPI_COMM_NULL) then
                  color= MPI_UNDEFINED ; ordinal_in_set= -1
                  if (myproc==0) then
                     color=1
                     ordinal_in_set=MiMd_Wmyproc
                  endif
                  call MPI_Comm_split (comm, color, ordinal_in_set, comm, ierr)
                  do k=2,F_nc
                     if ((F_color(k)==MiMd_world(i)%color).or.(F_color(k)==MiMd_world(j)%color)) F_COMM(k)= comm
                  enddo
               endif
            endif
         end do
      end do

      call check_host (F_component_S,F_color(1)) 
!     
!--------------------------------------------------------------------
!
      return

contains

      subroutine check_host (F_component_S,F_color)
      use iso_c_binding
      implicit none

      character(len=*), intent(IN) :: F_component_S
      integer, intent(IN) :: F_color

      interface
         integer(C_INT) function f_gethostid()BIND(C,name='gethostid')
         import :: C_INT
         end function f_gethostid
      end interface

      logical :: not_same_host=.false. ,shared_host=.false.
      integer i,myhost,err
      integer, dimension(2,0:MiMd_Wnumproc-1) :: host_list, glist
!     
!--------------------------------------------------------------------
!
      host_list= 0
      host_list(1,MiMd_Wmyproc) = F_color
      host_list(2,MiMd_Wmyproc) = abs(f_gethostid())
      call MPI_ALLREDUCE ( host_list, glist, 2*MiMd_Wnumproc,&
                     MPI_INTEGER,MPI_SUM,MiMd_gemworld,err)

      do i=0,MiMd_Wnumproc-1
         if (glist(1,i)==F_color) then
            myhost=glist(2,i)
            exit
         endif
      end do
      if (F_color .ne. 1) then
         print*, 'Checking hosts assigned to server ',trim(F_component_S)
         do i=0,MiMd_Wnumproc-1
            if ((glist(1,i)==F_color).and.(glist(2,i)/=myhost)) not_same_host=.true.
            if ((glist(1,i)==1).and.(glist(2,i)==myhost)) shared_host=.true.
         end do
         if (not_same_host) then
            write (6,987) ' CRITICAL ERROR: server ',trim(F_component_S), ' is assigned to multiple host '
         endif
         if (shared_host) then
            write (6,987) ' WARNING ERROR: server ',trim(F_component_S), ' is sharing hosts with GEM '
         endif
      endif
 987  format (//20("#"),3a,20("#")//)
!     
!--------------------------------------------------------------------
!
      return
      end subroutine check_host

      end subroutine MiMd_init
