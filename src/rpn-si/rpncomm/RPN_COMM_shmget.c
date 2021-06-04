#include <stdio.h>
#include <stdint.h>
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
        function rpn_comm_shmget(comm,size) result(where) bind(C,name='RPNCOMM_Shmget')     !InTf!
        import :: C_INT, C_PTR                                !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: comm           !InTf!   ! RPN_COMM Fortran communicator
          integer(C_INT), intent(IN), value :: size           !InTf!   ! size in bytes of shared memory area
          type(C_PTR) :: where                                !InTf!   ! pointer to shared memory area
        end function rpn_comm_shmget                          !InTf!
!InTf!
        function rpn_comm_numasplit(comm, mode, numacomm, numarank) result(size) bind(C,name='RPNCOMM_Numa_split')     !InTf!
        import :: C_INT                                       !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: comm           !InTf!   ! RPN_COMM Fortran communicator to split
          integer(C_INT), intent(IN), value :: mode           !InTf!   ! 0 : split by node, 1 : split by socket
          integer(C_INT), intent(OUT) :: numacomm             !InTf!   ! new communicator
          integer(C_INT), intent(OUT) :: numarank             !InTf!   ! rank in new communicator
          integer(C_INT) :: size                              !InTf!   ! size of new communicator
        end function rpn_comm_numasplit                       !InTf!
#endif

static int socket = -1;
static int core   = -1;

static void get_core_and_socket(){  // get socket/core info if available
  if(socket >= 0) return ;          // already initialized
#if defined(__x86_64__) &&  defined(__linux__)
   uint32_t lo, hi, c;
   // rdtscp instruction
   // EDX:EAX contain TimeStampCounter
   // ECX contains IA32_TSC_AUX[31:0] (MSR_TSC_AUX value set by OS, lower 32 bits contain socket+core)
   __asm__ volatile("rdtscp" : "=a" (lo), "=d" (hi), "=c" (c));
    socket = (c>>12) & 0xFFF;
    core   =  c      & 0xFFF;
#else
  socket = 0xFFFF;                     // no socket/core info available
  core   = 0xFFFF;
#endif
}

void *RPNCOMM_Shmget(int comm, unsigned int shm_size)
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

typedef struct{
  MPI_Fint fcom;
  MPI_Fint ncom;
  int mode;
  int rank;
  int size;
}cacheentry;
static int nused = 0;
#define MAX_CACHE 16
static cacheentry cache[MAX_CACHE];

int RPNCOMM_Numa_split(MPI_Fint origcomm, int mode, MPI_Fint *numacomm, int *numarank){
  MPI_Comm comm;
  MPI_Comm Tempcomm, Numacomm, Nodecomm;
  int i, myhost, myhost0, myhost1, err, myrank, newrank, newsize;

  for(i=0 ; i<nused ; i++){
    if(cache[i].fcom == origcomm && cache[i].mode == mode){
      *numacomm = cache[i].ncom;
      *numarank = cache[i].rank;
      return cache[i].size ;
    }
  }
  // not found in cache, split original Fortran communicator
  comm = MPI_Comm_f2c(origcomm);
  myhost  = gethostid();
  myhost0 = myhost & 0x7FFFFFFF;
  myhost1 = (myhost >> 31) & 0x1;
  err = MPI_Comm_rank(comm, &myrank);
  err = MPI_Comm_split(comm,     myhost0, myrank, &Tempcomm);  // split using lower 31 bits
  err = MPI_Comm_split(Tempcomm, myhost1, myrank, &Nodecomm);  // resplit using upper bit
  err = MPI_Comm_free(&Tempcomm);                              // free unused communicator

  if(mode == 0) { // node split
    err = MPI_Comm_rank(Nodecomm, &newrank);     // rank on node
    err = MPI_Comm_size(Nodecomm, &newsize);     // size of new communicator
    *numacomm = MPI_Comm_c2f(Nodecomm);
  }else{          // NUMA(socket) split
    get_core_and_socket();
    err = MPI_Comm_split(Nodecomm, socket, myrank, &Numacomm);  // resplit node into NUMA (socket)
    err = MPI_Comm_free(&Nodecomm);                             // free unused communicator
    err = MPI_Comm_rank(Numacomm, &newrank);     // rank on NUMA/socket
    err = MPI_Comm_size(Numacomm, &newsize);     // size of new communicator
    *numacomm = MPI_Comm_c2f(Numacomm);
  }
  *numarank = newrank;

  if(nused < MAX_CACHE){   // if cache table is not full
    cache[nused].fcom = origcomm;
    cache[nused].ncom = *numacomm;
    cache[nused].mode = mode;
    cache[nused].rank = newrank;
    cache[nused].size = newsize;
    nused++;
  }
  return newsize;
}
#if defined(SELF_TEST)
int main(int argc, char **argv){
  MPI_Fint origcomm, numacomm;
  int numarank, npes, myrank;
  int *addr;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  npes = RPNCOMM_Numa_split(MPI_Comm_c2f(MPI_COMM_WORLD), 1, &numacomm, &numarank);  // socket split
  printf("size = %d, rank = %d, core = %d, socket = %d\n",npes,numarank, core, socket);
  addr = (int *) RPNCOMM_Shmget(numacomm, 32768);
  printf("address = %p\n",addr);
  MPI_Barrier(MPI_COMM_WORLD);
  if(numarank == 0) {
    addr[0] = myrank;
    addr[1] = 2345;
    addr[2] = 3456;
    addr[3] = 4567;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  printf("addr[0:3] = %d %d %d %d\n",addr[0], addr[1], addr[2], addr[3]);
  MPI_Finalize();
}
#endif
