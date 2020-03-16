#!/bin/bash
set -ex
# for the following mpirun
ulimit -s 512000
cat <<EOT >make_mpif_include.f90
program make_mpif_includes
include 'mpif.h'
EOT
for i in MPI_2DOUBLE_PRECISION MPI_2INTEGER MPI_2REAL MPI_DOUBLE_PRECISION MPI_INTEGER \
         MPI_INTEGER4 MPI_INTEGER8 MPI_LOGICAL MPI_REAL MPI_REAL4 MPI_REAL8 MPI_SUCCESS
do
  echo "print 100,'#define $i ',$i" >>make_mpif_include.f90
done
cat <<EOT >>make_mpif_include.f90
100 format(A,I10)
end program
EOT
#s.f90 -mpi -o make_mpif_include make_mpif_include.f90
MPIRUN=mpirun
Fortran="s.f90 -mpi"
Ccompiler="s.cc -mpi"
if which aprun 1>/dev/null 2>/dev/null ; then
  MPIRUN=aprun
  Fortran=ftn
  Ccompiler=cc
fi
${Fortran} -o make_mpif_include make_mpif_include.f90
rm -f make_mpif_include.f90
if which poe 1>/dev/null 2>/dev/null ; then
  export MP_PROCS=1
  export MP_HOSTFILE=$TMPDIR/MP_HOSTFILE_$$
  echo localhost>$MP_HOSTFILE
  unset MP_LABELIO
  ./make_mpif_include | grep '#define'  >mpi_stub.h
  rm -f $MP_HOSTFILE
else
#  ${MPIRUN} -n 1 ./make_mpif_include | grep '#define' >mpi_stub.h
  ./make_mpif_include | grep '#define' >mpi_stub.h
fi
rm -f make_mpif_include.f90 make_mpif_include
if [[ "$1" == fortran || "$1" == all ]] ; then
  echo "CREATING: rpn_comm_fortran_stubs.F90"
  cat <<EOT >rpn_comm_fortran_stubs.F90
#include "mpi_stub.h"
      subroutine MPI_abort
      write(6,*) 'MPI_abort is called, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_allgather (a, cnt, b, c, cnt2, d, e, ierr )
        implicit none
        !include 'mpi_stub.h'
        integer MPI_STUBS_length
        integer a(*),cnt,b,c(*),cnt2,d,e,ierr,i
        do i=1,min(cnt,cnt2)*MPI_STUBS_length(b)
           c(i)=a(i)
        enddo
        ierr=MPI_SUCCESS
        return
      end
!
      subroutine MPI_allgatherv
      write(6,*) 'MPI_allgatherv not authorized in non-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_allreduce(send,recv,ni,l,m,n,ierr)
      implicit none
!
      !include 'mpi_stub.h'
!
      integer MPI_STUBS_length
      integer i,l,m,n,ierr
      integer send(*),recv(*),ni
      do i=1,ni*MPI_STUBS_length(l)
        recv(i)=send(i)
      enddo
      ierr=MPI_SUCCESS
      return
      end
!
      subroutine MPI_alltoall
      write(6,*)'MPI_alltoall not authorized in non-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_alltoallv
      write(6,*)'MPI_alltoallv not authorized in non-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_barrier(i,ierr)
      implicit none
      !include 'mpi_stub.h'
      integer i,ierr
      ierr=MPI_SUCCESS
      return
      end
!
      subroutine MPI_bcast(i,j,k,l,m,ierr)
      implicit none
!
      !include 'mpi_stub.h'
!
      integer i,j,k,l,m,ierr
      ierr=MPI_SUCCESS
      return
      end
!
      subroutine mpi_comm_create(pe_indomm,pe_gr_blocmaster, pe_blocmaster, ierr)
      implicit none
      !include 'mpi_stub.h'
      integer pe_indomm,pe_gr_blocmaster,  pe_blocmaster, ierr
      pe_blocmaster = 0
      ierr = MPI_SUCCESS
      return
      end
!
      subroutine MPI_comm_group(i,group,ierr)
      implicit none
      integer i,group,ierr
!
      !include 'mpi_stub.h'
!
      group=-1
      ierr=MPI_SUCCESS
      return
      end
!
      subroutine MPI_COMM_RANK(comm,pe_me,ierr)
      implicit none
!
      !include 'mpi_stub.h'
!
      integer comm,pe_me,ierr
      ierr=MPI_SUCCESS
      pe_me=0
      return
      end
!
      subroutine MPI_COMM_SIZE(comm,pe_tot,ierr)
      implicit none
!
      !include 'mpi_stub.h'
!
      integer comm,pe_tot,ierr
      ierr=MPI_SUCCESS
      pe_tot=1
      return
      end
!
      subroutine MPI_comm_split(i,j,k,newcomm,ierr)
      implicit none
      integer i,j,k,newcomm,ierr
!
      !include 'mpi_stub.h'
!
      newcomm=-1
      ierr=MPI_SUCCESS
      return
      end
!
      subroutine MPI_finalize
      return
      end
!
      subroutine MPI_gather(a, cnt, b, c, cnt2, d, e, f,ierr )
      implicit none
      !include 'mpi_stub.h'
      integer MPI_STUBS_length
      integer a(*),cnt,b,c(*),cnt2,d,e,f,ierr,i
      do i=1,min(cnt,cnt2)*MPI_STUBS_length(b)
         c(i)=a(i)
      enddo
      ierr=MPI_SUCCESS
      return
      end
!
      subroutine MPI_gatherv(a, cnt, b, c, cnt2s, displ, d, e, f,ierr )
      implicit none
      !include 'mpi_stub.h'
      integer MPI_STUBS_length
      integer a(*),cnt,b,c(*),cnt2s(*),cnt2,d,e,f,ierr,i,displ(*)
        cnt2=cnt2s(1)
      do i=1,min(cnt,cnt2)*MPI_STUBS_length(b)
         c(i+displ(1))=a(i+displ(1))
      enddo
      ierr=MPI_SUCCESS
      return
      end
!
      subroutine MPI_get_count
      write(6,*) 'MPI_get_count not authorized in non-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_GROUP_incl(pe_gr_wcomm, pe_dommtot,proc_indomm,pe_gr_indomm,ierr)
      implicit none
      !include 'mpi_stub.h'
      integer pe_gr_wcomm, pe_dommtot,proc_indomm,pe_gr_indomm,ierr
      pe_gr_indomm = 0
      ierr = MPI_SUCCESS
      return
      end
!
      subroutine MPI_GROUP_RANK(comm,pe_me,ierr)
      implicit none
!
      !include 'mpi_stub.h'
!
      integer comm,pe_me,ierr
      ierr=MPI_SUCCESS
      pe_me=0
      return
      end
!
      subroutine MPI_GROUP_SIZE(comm,pe_tot,ierr)
      implicit none
!
      !include 'mpi_stub.h'
!
      integer comm,pe_tot,ierr
      ierr=MPI_SUCCESS
      pe_tot=1
      return
      end
!
      subroutine MPI_init(ierr)
      implicit none
      integer pe_tot
      save pe_tot
!
      !include 'mpi_stub.h'
!
      integer ierr
      data pe_tot / -1/
      ierr=0
!     if(pe_tot .ne. -1) call ABORT
      pe_tot=0
      return
      end
!
      subroutine MPI_INITIALIZED(mpi_started,ierr)
      implicit none
!
      !include 'mpi_stub.h'
!
      logical mpi_started
      integer ierr
!
      ierr=MPI_SUCCESS
      mpi_started=.false.
      return
      end
!
      subroutine MPI_irecv
      write(6,*) 'MPI_irecv not authorized in non-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_isend
      write(6,*) 'MPI_isend not authorized in non-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_recv
      write(6,*) 'MPI_recv not authorized in non-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_reduce(send,recv,ni,l,m,n,o,ierr)
      implicit none
!
      !include 'mpi_stub.h'
!
      integer MPI_STUBS_length
      integer i,l,m,n,o,ierr
      integer send(*),recv(*),ni
      do i=1,ni*MPI_STUBS_length(l)
        recv(i)=send(i)
      enddo
      ierr=MPI_SUCCESS
      return
      end
!
      subroutine MPI_scatterv
      write(6,*) 'MPI_scatterv not authorized in non-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_send
      write(6,*) 'MPI_send not authorized in non-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_sendrecv
      write(6,*)  'MPI_sendrecv not authorized in non-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_ssend
      write(6,*) 'MPI_ssend not authorized in non-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_comm_free
      write(6,*) 'MPI_comm_free not authorized in non-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_type_get_extent(dtyp_m,extent,ierr)  ! replaces MPI_type_extent
      integer,intent(IN) :: dtyp_m
      !integer, intent(out) :: ierr
      integer :: ierr
      integer :: extent
      write(6,*) 'MPI_type_get_extent not authorized in non-mpi mode, ABORT'
      call ABORT
      ierr = -1
      return
      end
      subroutine MPI_type_extent(dtyp_m,extent,ierr)  ! old, deprecated entry
      integer,intent(IN) :: dtyp_m
      !integer, intent(out) :: ierr
      integer :: ierr
      integer :: extent
      write(6,*) 'MPI_type_extent not authorized in non-mpi mode, ABORT'
      call ABORT
      ierr = -1
      return
      end
!
      integer function MPI_STUBS_length(itype)
      implicit none
      !include 'mpi_stub.h'
      integer itype

      MPI_STUBS_length=0

      if(itype .eq. MPI_DOUBLE_PRECISION) MPI_STUBS_length=2
      if(itype .eq. MPI_2DOUBLE_PRECISION) MPI_STUBS_length=2
      if(itype .eq. MPI_REAL8) MPI_STUBS_length=2
      if(itype .eq. MPI_INTEGER8) MPI_STUBS_length=2

      if(itype .eq. MPI_LOGICAL) MPI_STUBS_length=1

      if(itype .eq. MPI_INTEGER) MPI_STUBS_length=1
      if(itype .eq. MPI_INTEGER4) MPI_STUBS_length=1
      if(itype .eq. MPI_2INTEGER) MPI_STUBS_length=1
      if(itype .eq. MPI_REAL) MPI_STUBS_length=1
      if(itype .eq. MPI_REAL4) MPI_STUBS_length=1
      if(itype .eq. MPI_2REAL) MPI_STUBS_length=1

        if(MPI_STUBS_length.eq.0)then
        write(6,*)'MPI_STUBS_length ERROR: unrecognized type'
        call abort
      endif
      return
      end
!
      subroutine MPI_wait
      return
      end
!
      subroutine MPI_waitall
      return
      end
!
      REAL*8 function  MPI_wtick()
      MPI_wtick = 1.0E-9
      return
      end
      REAL*8 function  PMPI_wtick()
      PMPI_wtick = 1.0E-9
      return
      end
!
      REAL*8 function  MPI_wtime()
      real *8, save :: dummy_time=1.0E-9
      MPI_wtime=dummy_time
      dummy_time=dummy_time+1.0E-9
      return
      end
      REAL*8 function  PMPI_wtime()
      real *8, save :: dummy_time=1.0E-9
      PMPI_wtime=dummy_time
      dummy_time=dummy_time+1.0E-9
      return
      end
!
EOT
  ${Fortran} -c rpn_comm_fortran_stubs.F90
fi
if [[ "$1" == c || "$1" == all ]] ; then
  echo "CREATING: rpn_comm_c_stubs.c"
  cat <<EOT >rpn_comm_c_stubs.c
#include <stdio.h>
#include <stdlib.h>
#include "mpi_stub.h"
int MPI_Comm_rank(int comm, int *rank){
   *rank = 0;
   return(MPI_SUCCESS);
}
int MPI_Comm_size(int comm, int *size){
   *size = 1;
   return(MPI_SUCCESS);
}
int MPI_Allgather(void *outx, int nout, int outtype, void *inx, int nin, int intype, int comm){
   int *out = outx;
   int *in = inx;
/*   if( nin != 1 || nout != 1 || outtype != MPI_INTEGER || intype != MPI_INTEGER ) {  */
   if( nin != 1 || nout != 1 ) {
/*      fprintf(stderr,"ERROR: MPI_Allgather, number of elements is not one or type is not MPI_INTEGER\n");  */
      fprintf(stderr,"ERROR: MPI_Allgather, number of elements is not one \n");
      exit(1);
   *out = *in;
   }
   return(MPI_SUCCESS);
}
EOT
  ${Ccompiler} -c rpn_comm_c_stubs.c
fi
if [[ "$1" == "none" ]] ; then
echo "CREATING: rpn_comm_fortran_stubs.F90 and rpn_comm_c_stubs.c with NO stubs"
cat <<EOT >rpn_comm_fortran_stubs.f90
      subroutine no_mpi_ftn_stubs
      return
      end
EOT
  ${Fortran} -c rpn_comm_fortran_stubs.F90
cat <<EOT >rpn_comm_c_stubs.c
int no_mpi_c_stubs(){
   return(0);
}
EOT
  ${Ccompiler} -c rpn_comm_c_stubs.c
fi

