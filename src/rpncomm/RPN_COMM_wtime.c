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
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>

/* dummy wall time clock, pseudo clock increases at one nanosecond per call */
static double dummy_time=0.0;
static double dummy_wtime()
{
   dummy_time += 1.0E-09;
   return dummy_time;
}

static double (*fn)() = dummy_wtime ;
static double time0 = 0.0;

/* returns wall time as a 64 bit real number */
/* example FORTRAN CALL:                     */
/* real *8 time                              */
/* time = RPN_COMM_Wtime()                   */
#if defined(is_never_defined)
#pragma weak rpn_comm_wtime__=rpn_comm_wtime
#pragma weak rpn_comm_wtime_=rpn_comm_wtime
double rpn_comm_wtime__();
double rpn_comm_wtime_();

double rpn_comm_wtime()
#endif
double Rpn_comm_wtime()
{
  return (*fn)() - time0;
}

/* returns time stamp counter as a 64 bit real number */
#if defined(linux)
typedef unsigned long long ticks;
static int t0_=0;
static ticks t0;
static double cpu_mult=1.0E-9;  /* 1GHz */
static ticks getticks(void)
{
#if defined(__x86_64__)
     unsigned a, d;
     asm("cpuid");
     asm volatile("rdtsc" : "=a" (a), "=d" (d));

     return (((ticks)a) | (((ticks)d) << 32));
#else
     double MPI_Wtime();
     ticks tmp = MPI_Wtime() * 1000000.0; // microseconds
     return tmp;
#endif
}

#endif
#if defined(AIX)
static int t0_=0;
static timebasestruct_t t0;
#endif

#pragma weak rpn_comm_tsc__=rpn_comm_tsc
#pragma weak rpn_comm_tsc_=rpn_comm_tsc
double rpn_comm_tsc__();
double rpn_comm_tsc_();

double rpn_comm_tsc()
{
  double temp;
#if defined(linux)
  char *junk;
  ticks value;
  FILE *cpuInfo;
  char buffer[1024];
  char *pbuf=buffer;
  float freq=1000.0;
  if(t0_ == 0) {
    if( (cpuInfo = fopen("/proc/cpuinfo", "r")) != NULL  ) {
      while(1) {
	junk=fgets(buffer,sizeof(buffer),cpuInfo);
	if(buffer[4]=='M' && buffer[5]=='H' && buffer[6]=='z') {
	  while(*pbuf != ':') pbuf++; pbuf++;
	  sscanf(pbuf,"%f",&freq);
	  freq *= 1000000.0;
	  cpu_mult = 1.0 / freq;
	  printf("freq = %f Hz, cpu_mult=%G\n",freq,cpu_mult);
	  break;
	}
      }
      fclose(cpuInfo);
    }
    t0=getticks();
    t0_=1;
  }
  value = getticks();
  temp = value - t0;
  temp *= cpu_mult;   /* translate ticks to seconds */
#endif
#if defined(AIX)
  timebasestruct_t t;
  if(t0_ == 0) {
    read_real_time(&t0,TIMEBASE_SZ);
    time_base_to_time(&t0,TIMEBASE_SZ);
    t0_=1;
  }
  read_real_time(&t,TIMEBASE_SZ);
  time_base_to_time(&t,TIMEBASE_SZ);
  t.tb_low = t.tb_low - t0.tb_low;
  t.tb_high = t.tb_high - t0.tb_high;
  if(t.tb_low < 0){
    t.tb_high--;
    t.tb_low += 1000000000;
  }
  temp = t.tb_low;
  temp *= 1.0E-9;
  temp += t.tb_high;

#endif
#if ! defined(linux) && ! defined(AIX)
  temp = dummy_time;
  dummy_time += 1.0E-09;
#endif
return temp;
}

/* returns time of day as a 64 bit real number */

#pragma weak rpn_comm_timeofday__=rpn_comm_timeofday
#pragma weak rpn_comm_timeofday_=rpn_comm_timeofday
double rpn_comm_timeofday__();
double rpn_comm_timeofday_();
double rpn_comm_timeofday()
{
  double temp;
  struct timeval tv;
  gettimeofday(&tv,NULL);
  temp = tv.tv_usec;
  temp *= .000001;
  temp += tv.tv_sec;
  return temp;
}

/* set wall clock function used by Rpn_comm_wtime to arbitrary clock function */
/* that function MUST return a 64 bit real number (e.g MPI clock MPI_Wtime) */
/* example FORTRAN call : call RPN_COM_Wtime_set(MPI_Wtime) */
/* example FORTRAN call : call RPN_COM_Wtime_set(RPN_COMM_Timeofday)   */
/* special case :                                                      */
/*   call RPN_COM_Wtime_set(RPN_COM_Wtime_set)                         */
/*   resets the clock to the very fast dummy internal clock            */


#pragma weak rpn_comm_wtime_set__=rpn_comm_wtime_set
#pragma weak rpn_comm_wtime_set_=rpn_comm_wtime_set
void rpn_comm_wtime_set__(double (*)());
void rpn_comm_wtime_set_(double (*)());
void rpn_comm_wtime_set(double (*function)())
{
  fn = function;
  if(fn == (void *)rpn_comm_wtime_set) fn = (void *)dummy_wtime;
  time0 = (*fn)();
}

#ifdef TEST
#include <stdio.h>
#include <mpi.h>
int main(int argc,char **argv)
{
  int i;
  int ierr;
  double dummy, dummy2, dummy3, dummy4;
  ticks t1, t2, t3, t4, t5;
  
  ierr = MPI_Init(&argc, &argv);
  
  dummy = rpn_comm_tsc();
  
  t1 = getticks();
  t2 = getticks();
  t3 = getticks();
  t4 = getticks();
  t5 = getticks();
  fprintf(stdout,"%u %u %u %u \n",t2-t1,t3-t2,t4-t3,t5-t4);
  
  fprintf(stdout,"Phase 1, dummy timing function\n");
  for (i=0 ; i<5 ; i++){
    double x = Rpn_comm_wtime();
    double x1 = Rpn_comm_wtime();
    int loop = 1;
    while(x == x1) { x1 = Rpn_comm_wtime(); loop++; }
    fprintf(stdout,"TIME1= %G %d\n",x1-x,loop);
  }
  fprintf(stdout,"Phase 2, using MPI function MPI_Wtime\n");
  rpn_comm_wtime_set(MPI_Wtime);
  for (i=0 ; i<5 ; i++){
    double x = Rpn_comm_wtime();
    double x1 = Rpn_comm_wtime();
    int loop = 1;
    while(x == x1) { x1 = Rpn_comm_wtime(); loop++; }
    fprintf(stdout,"TIME2= %G %d\n",x1-x,loop);
  }
  fprintf(stdout,"Phase 3, using get_time_of_day\n");
  rpn_comm_wtime_set(rpn_comm_timeofday);
  for (i=0 ; i<5 ; i++){
    double x = Rpn_comm_wtime();
    double x1 = Rpn_comm_wtime();
    int loop = 1;
    while(x == x1) { x1 = Rpn_comm_wtime(); loop++; }
    fprintf(stdout,"TIME3= %G %d\n",x1-x,loop);
  }
  fprintf(stdout,"Phase 4, dummy timing function\n");
  rpn_comm_wtime_set((void *)rpn_comm_wtime_set);
  for (i=0 ; i<5 ; i++){
    double x = Rpn_comm_wtime();
    double x1 = Rpn_comm_wtime();
    int loop = 1;
    while(x == x1) { x1 = Rpn_comm_wtime(); loop++; }
    fprintf(stdout,"TIME4= %G %d\n",x1-x,loop);
  }
  fprintf(stdout,"Phase 5, TSC timing function\n");
  rpn_comm_wtime_set((void *)rpn_comm_tsc);
  for (i=0 ; i<5 ; i++){
    double x = Rpn_comm_wtime();
    double x1 = Rpn_comm_wtime();
    int loop = 1;
    while(x == x1) { x1 = Rpn_comm_wtime(); loop++; }
    fprintf(stdout,"TIME5= %G %d\n",x1-x,loop);
  }
  MPI_Finalize();
  return 0;
}
#endif
