/* RPN_COMM - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2012  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/*
   unset any processor affinity for the current process
   if the FULL_UNBIND environment variable is defined
   ( FULL_UNBIND would normally be set by r.mpirun / ord_soumet )
   in the case of an MPI launch, only process 0 will print a message
   this code is only applicable to Linux for the time being
   and has been written to counter mpi implementation behaviour in some cases.
   (written specifically for use in the RPN_COMM library)

   this code is Linux only

   Michel Valin , 2011 / 11 / 02 UQAM (openmpi)
                  2011 / 11 / 20 UQAM (mpich / slurm)
*/

/* use pragma weak to create alternate FORTRAN entry points */
#pragma weak rpn_comm_unbind_process_=rpn_comm_unbind_process
#pragma weak rpn_comm_unbind_process__=rpn_comm_unbind_process


/* WARNING: some old versions of gcc do not generate the weak entry points correctly */

void rpn_comm_unbind_process__(void);
void rpn_comm_unbind_process_(void);
void rpn_comm_unbind_process(void)
{
#if defined(USE_OLD_CODE)
#if defined(linux)

#define _GNU_SOURCE
#include <sched.h>
#include <stdio.h>
#include <unistd.h>

int atoi(const char *nptr);
char *getenv(const char *name);

cpu_set_t set;
int i;
int will_print=1;
char *ompi=getenv("OMPI_COMM_WORLD_RANK");  /* openmpi */
int ncores=sysconf(_SC_NPROCESSORS_CONF);
int nthreads=1;
int nbound=0;
char *omp=getenv("OMP_NUM_THREADS");

if(omp != NULL) nthreads=atoi(omp);  /* OMP_NUM_THREADS value */

if(ompi == NULL) ompi=getenv("PMI_RANK");   /* mpich family */

if(ompi != NULL) if(0 != atoi(ompi)) will_print=0;  /* not MPI process 0, no message */

CPU_ZERO(&set);
sched_getaffinity(0,sizeof(set),&set);
i=ncores;
while(--i >=0) { if (CPU_ISSET(i,&set)){ nbound++; } }  /* how many cores are we allowed to run on ? */
  
if(getenv("FULL_UNBIND") != NULL) nbound = 0;  /* FULL_UNBIND variable defined, unbind no matter what */
  
if(nthreads > nbound) {  /* need more threads than cores we can run on , unbind everything */
  if(will_print) printf("FULL unbinding will be done, cores=%d, threads needed=%d, usable cores=%d\n",ncores,nthreads,nbound);
  CPU_ZERO(&set);
  i=ncores;
  while(--i >=0) { CPU_SET(i,&set);}
  sched_setaffinity(0,sizeof(set),&set);  /* set affinity to all cores */
}else{
  if(will_print) printf("NO unbinding will be done\n");  /* enough resources available and no forced unbind */
}

#endif
#endif

return;
}
