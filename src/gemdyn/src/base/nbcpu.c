#include <stdio.h>
#include <errno.h>
#include <mpi.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <string.h>
#include <limits.h>

int cpu_per_numa() {
   /* csubich -- 27 July 2021 */

   /* To implement the single-transpose solver with the minimum amount of communication,
      we want to divide the grid up into numa-local groupings, each of which will have
      a shared-memory region.

      There is no standard-compliant way to accomplish this with just the MPI library.
      OpenMPI has a few different options for this, but they are nonstandard extensions.
      Instead, we will interrogate the system more directly.  This is Linux-specific, but
      for now it is not a limiting factor.

      One approach to this is to use the 'hwloc' library, which describes the hardware
      topolgoy directly.  However, operational requirements mean that everything needs
      to be statically linked, and there is currently no static version of this library
      on the back-end machines.

      Instead, we will use the getcpu system call to interrogate the current processor/
      NUMA allocation.  This system call takes the form:
      
      getcpu(&cpu,&numa,&unused)

      and it returns the current CPU being used by the active process (&cpu) along
      with its current NUMA node (&numa).

      By using this system call on each processor and then aggregating results, we
      can find out the effective NUMA region size.  Moreover, we can also immediately
      tell if processors have not been distributed evenly among NUMA nodes, allowing
      us to return an error flag (-1).

      Unfortunately, this system call is not conveniently wrapped for us by the glibc
      library, so we must use the direct syscall interface. */
   
   int err = 0;
   unsigned cpu, node;
   err = syscall(SYS_getcpu, &cpu, &node, 0); // make the system call

   // This system call shouldn't result in an error on a Linux system, but it may
   // if we're trying to compile on another operating system.  Note the error but
   // return somewhat gracefully, since execution can continue without the single-
   // transpose solver.
   if (err != 0) {
      fprintf(stderr,"Error %d (%s) in getcpu syscall!\n",errno,strerror(errno));
   }

   // Check whether any process has had an error, and if so return an error flag
   int glb_err = 0;
   MPI_Allreduce(&err,&glb_err,1,MPI_INT,MPI_LOR,MPI_COMM_WORLD);
   if (glb_err != 0) {
      return -1;
   }
   
   // Proceed, since the system call was successful.

   /* cpu contains the current CPU, node contains the current NUMA node.  First,
      find the maximum numa node over all processes */

   int max_node = 0;
   MPI_Allreduce(&node,&max_node,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

   /* Now, we want to count the number of processes using each NUMA node 
      *per compute node*.  To do this, first create arrays of the potential
      NUMA node membership. */
   int numa_nodes_rank[max_node+1]; //For this individual rank only
   int numa_nodes_host[max_node+1]; //For all ranks on this host
   for (int i = 0; i <= max_node; i++) {
      numa_nodes_rank[i] = 0; // initialize to 0
      numa_nodes_host[i] = 0;
   }
   numa_nodes_rank[node] = 1; // Count this rank's single NUMA node


   // Split the comm_world MPI communicator into one per host
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   MPI_Comm host_comm;
   MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,myrank,MPI_INFO_NULL,&host_comm);

   int hostrank;
   MPI_Comm_rank(host_comm,&hostrank);

   // Accumulate the NUMA membership on this host
   MPI_Allreduce(numa_nodes_rank,numa_nodes_host,max_node+1,MPI_INT,MPI_SUM,host_comm);

   int max_cpu_per_numa_host = 0;
   int min_cpu_per_numa_host = INT_MAX;
   /* Loop over each potential NUMA node to find the maximum and minium occupancy on
      this host.  The minimum is found for error checking, since we need NUMA nodes
      to fill evenly */
   for (int i = 0; i <= max_node; i++) {
      if (numa_nodes_host[i] > max_cpu_per_numa_host) {
         max_cpu_per_numa_host = numa_nodes_host[i];
      }
      if (numa_nodes_host[i] > 0 && numa_nodes_host[i] < min_cpu_per_numa_host) {
         min_cpu_per_numa_host = numa_nodes_host[i];
      }
   }

   int max_cpu_per_numa_global = 0;
   int min_cpu_per_numa_global = 0;

   /* Finally, reduce over all ranks (globally) to find the global maximum and minimum
      NUMA occupancies */

   MPI_Allreduce(&max_cpu_per_numa_host,&max_cpu_per_numa_global,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
   MPI_Allreduce(&min_cpu_per_numa_host,&min_cpu_per_numa_global,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);

   // Deallocate the host communicator to avoid leaking the handle
   MPI_Comm_free(&host_comm);

   // If the maximum and minimum differ, return -1 as an error condition
   if (max_cpu_per_numa_global != min_cpu_per_numa_global) {
      return -1;
   }

   // Otherwise, return the maximum occupancy per numa node, which we've verified is the
   // same over all hosts
   return max_cpu_per_numa_global;
   
}
