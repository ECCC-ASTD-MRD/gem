#include <stdio.h>
#include <stdlib.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <mpi.h>
#include <unistd.h>


/* this function needs an explicit fortran interface using ISO_C_BINDING */
/* because it returns a C pointer */

#ifdef MUST_NEVER_BE_TRUE
!InTf!
        function rpn_comm_shmget(comm,size) result(where)     !InTf!
        use ISO_C_BINDING                                     !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN) :: comm                  !InTf!   ! RPN_COMM communicator
          integer(C_INT), intent(IN) :: size                  !InTf!   ! size in bytes of shared memory area
          type(C_PTR) :: where                                !InTf!   ! pointer to shared memory area
        end function rpn_comm_shmget                          !InTf!
#endif

void *f_rpn_comm_shmget(int comm, unsigned int shm_size)
{
  size_t size=shm_size;
  int id;
  struct shmid_ds shm_buf;
  void *ptr;
  int ierr, myrank, myhost, myhost2;
  MPI_Fint f_comm=comm;
  MPI_Comm c_comm;
  

  c_comm = MPI_Comm_f2c(f_comm);  /* translate Fortran communicator into C communicator */

  myhost=gethostid();
  ierr=MPI_Allreduce(&myhost,&myhost2,1,MPI_INTEGER,MPI_BOR,c_comm); /* boolean OR of hostid from all members of this comunicator */
  if(myhost != myhost2){
    fprintf(stderr,"ERROR: (rpn_comm_shmget) processes are not all on same node \n");
    return NULL;
    }
  ierr=MPI_Comm_rank(c_comm,&myrank);
  if(myrank == 0) id=shmget(IPC_PRIVATE,size,IPC_CREAT|S_IRUSR|S_IWUSR);  /* rank 0 allocates shared memory */
  ierr=MPI_Bcast(&id,1,MPI_INTEGER,0,c_comm);                             /* all processes get id */
  printf("id = %d %x \n",id,id);
  ptr=shmat(id,NULL,0);                                                   /* all processes attach memory segment */
  if(ptr == NULL) printf("ERROR: (rpn_comm_shmget) got a null pointer from shmat\n");
  ierr=MPI_Barrier(c_comm);                                               /* all processes hgave attached the segment */
//  printf("addr = %16p\n",ptr);
  if(myrank == 0) shmctl(id,IPC_RMID,&shm_buf);      /* mark segment for deletion to make sure it is released when all processes terminate */
  return ptr;                                        /* return pointer tio shared memory area */
}

