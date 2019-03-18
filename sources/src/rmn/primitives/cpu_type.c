/* RMNLIB - Library of hopefully useful routines for C and FORTRAN
 * Copyright (C) 2016  Recherche en Prevision Numerique
 *                     Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
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
/****p* X86/ utilities
 * FUNCTION
 *  set of X86 utilities for timing and cpu type and properties determination
 * 
 * NOTES
 *    Fortran example:
 * 
 *      program test_feature
 *        use ISO_C_BINDING
 *        implicit none
 *      #define IN_FORTRAN_CODE
 *      #include <cpu_type.h>
 *        real(C_DOUBLE) :: s1, s2
 *        integer(C_INT64_T) :: t1, t2
 *        if(cpu_has_feature(FLAG_AVX2) .ne. 0) print *,'AVX2'
 *        if(cpu_has_feature(FLAG_AVX) .ne. 0) print *,'AVX'
 *        if(cpu_has_feature(FLAG_FMA) .ne. 0) print *,'FMA'
 *        if(cpu_has_feature(FLAG_BMI) .ne. 0) print *,'BMI'
 *        if(cpu_has_feature(FLAG_SSE4) .ne. 0) print *,'SSE4'
 *        if(cpu_has_feature(FLAG_SSE3) .ne. 0) print *,'SSE3'
 *        if(cpu_has_feature(FLAG_SSE2) .ne. 0) print *,'SSE2'
 *        status = get_fp_status_ctl()
 *        print 100,'FP status = ',status
 *        if(iand(FP_STATUS_PE,status) .ne. 0) print *,"Precision ON"
 *        if(iand(FP_STATUS_UE,status) .ne. 0) print *,"Underflow ON"
 *        if(iand(FP_STATUS_OE,status) .ne. 0) print *,"Overflow ON"
 *        if(iand(FP_STATUS_ZE,status) .ne. 0) print *,"Zero divide ON"
 *        print *,'CPU apicid =',get_cpu_id()
 *        s1 = rdtsc_seconds()
 *        s2 = rdtsc_seconds()
 *        print *,'rdtsc_seconds overhead=',s2
 *        s1 = rdtscp_seconds()
 *        s2 = rdtscp_seconds()
 *        print *,'rdtscp_seconds overhead=',s2
 *        t1 = rdtsc()
 *        t2 = rdtsc()
 *        print *,'rdtsc overhead =',t2-t1,wall_clock_seconds(t2-t1)
 *        t1 = rdtscp()
 *        t2 = rdtscp()
 *        print *,'rdtscp overhead =',t2-t1,wall_clock_seconds(t2-t1)
 *        stop
 *      end
 * 
 *    C example:
 * 
 *      #include <stdint.h>
 *      #include <stdio.h>
 *      #include <cpu_type.h>
 *      int main(int argc, char** argv)
 *      {
 *       uint64_t t1, t2;
 *       double s1, s2;
 *       printf("CPU speed: %lu Hz\n",get_cpu_freq());
 *       printf("FLAGS: ");
 *       if(cpu_has_feature(FLAG_SSE))  printf(" SSE");
 *       if(cpu_has_feature(FLAG_SSE2)) printf(" SSE2");
 *       if(cpu_has_feature(FLAG_SSE3)) printf(" SSE3");
 *       if(cpu_has_feature(FLAG_SSE4)) printf(" SSE4");
 *       if(cpu_has_feature(FLAG_AVX))  printf(" AVX");
 *       if(cpu_has_feature(FLAG_FMA))  printf(" FMA");
 *       if(cpu_has_feature(FLAG_AVX2)) printf(" AVX2");
 *       printf("\n");
 *       status = get_fp_status_ctl();
 *       printf("FPU status = %8.8x\n",status);
 *       printf("Precision status %s\n",status & FP_STATUS_PE ? "ON" : "OFF");
 *       printf("Underflow status %s\n",status & FP_STATUS_UE ? "ON" : "OFF");
 *       printf("Overflow status %s\n",status & FP_STATUS_OE ? "ON" : "OFF");
 *       printf("Zero divide status %s\n",status & FP_STATUS_ZE ? "ON" : "OFF");
 *       t1 = rdtsc();
 *       t2 = rdtsc();
 *       printf("rdtsc overhead = %lu %E\n",t2-t1,wall_clock_seconds(t2-t1));
 *       t1 = rdtscp();
 *       t2 = rdtscp();
 *       printf("rdtscp overhead = %lu %E\n",t2-t1,wall_clock_seconds(t2-t1));
 *       s1 = rdtsc_seconds();
 *       s2 = rdtsc_seconds();
 *       printf("rdtsc_seconds overhead = %E\n",s2);
 *       s1 = rdtscp_seconds();
 *       s2 = rdtscp_seconds();
 *       printf("rdtscp_seconds overhead = %E\n",s2);
 *       return(0);
 *      }
 *
 *    result of C program execution:
 * 
 *     CPU speed: 3700000000 Hz
 *     FLAGS:  SSE SSE2 SSE3 SSE4 AVX FMA AVX2
 *     rdtsc overhead = 36 9.729730E-09
 *     rdtscp overhead = 124 3.351351E-08
 *     rdtsc_seconds overhead = -4.419891E+05
 *     rdtscp_seconds overhead = -2.312973E-06
*****
*/
#include <stdint.h>
// cc -I. -DTEST_CPUID -DTEST_RDTSC -DDEBUG cpu_type.c
#include <stdio.h>

#include <cpu_type.h>

#define HAS_FLAG(flag,flagword) (flag & flagword)
#define X86_FLAG(flag) (flag & ProcessorCapabilities)

static int ProcessorCapabilities = 0 ;  /* by default, no capabilities are indicated as available */
static uint32_t vstring[12];            /* version string */
static unsigned char *cstring = (unsigned char *)vstring;  /* pointer to version string */
static uint64_t hz=0;
static double cycle=0.0;
static uint64_t last_time=0;
static int ncores=0;
static int threadpercore=1;

static void X86_cpuid(uint32_t eax, uint32_t ecx, uint32_t* regs)  /* interface to x86 cpuid instruction */
{
#if defined(__x86_64__) || defined( __i386__ )
    uint32_t ebx, edx;
#if defined(DEBUG)
    uint32_t ieax=eax, iecx=ecx;
#endif
# if defined( __i386__ ) && defined ( __PIC__ )
     /* PIC under 32-bit EBX must not be clobbered */
    __asm__ ( "movl %%ebx, %%edi \n\t cpuid \n\t xchgl %%ebx, %%edi" : "=D" (ebx),
# else
   ebx = 0;
    __asm__ ( "cpuid" : "+b" (ebx),
# endif
              "+a" (eax), "+c" (ecx), "=d" (edx) );
    regs[0] = eax; regs[1] = ebx; regs[2] = ecx; regs[3] = edx;
#if defined(DEBUG)
    printf("BEFORE: eax = %8.8x,                 ecx = %8.8x\n",ieax,iecx);
    printf("AFTER : eax = %8.8x, ebx = %8.8x, ecx = %8.8x, edx = %8.8x\n",regs[0],regs[1],regs[2],regs[3]);
#endif
#else
    regs[0] = 0; regs[1] = 0; regs[2] = 0; regs[3] = 0;
#endif
}     

static void get_cpu_capabilities()
{
#if defined(__x86_64__)
  uint32_t regs[4], eax1;;
  int j;
  float freq;

  X86_cpuid( 1, 0, regs );  /* get CPU capabilities EAX=1, ECX=0 */
  eax1 = regs[0];
  if((1 <<  0) & regs[2]) ProcessorCapabilities |= FLAG_SSE3 ;  /* SSE3   ECX bit  0 */
  if((1 << 12) & regs[2]) ProcessorCapabilities |= FLAG_FMA ;   /* FMA    ECX bit 12 */
  if((1 << 20) & regs[2]) ProcessorCapabilities |= FLAG_SSE4 ;  /* SSE4.2 ECX bit 20 */
  if((1 << 28) & regs[2]) ProcessorCapabilities |= FLAG_AVX ;   /* AVX    ECX bit 28 */
  if((1 << 25) & regs[3]) ProcessorCapabilities |= FLAG_SSE ;   /* SSE    EDX bit 25 */
  if((1 << 26) & regs[3]) ProcessorCapabilities |= FLAG_SSE2 ;  /* SSE2   EDX bit 26 */

  X86_cpuid( 0x80000002, 0, regs );  /* version string (3 calls) */
  for(j=0 ; j<4 ; j++){
    vstring[j] = regs[j];
  }
  X86_cpuid( 0x80000003, 0, regs );
  for(j=0 ; j<4 ; j++){
    vstring[j+4] = regs[j];
  }
  X86_cpuid( 0x80000004, 0, regs );
  for(j=0 ; j<4 ; j++){
    vstring[j+8] = regs[j];
  }
  j = 0 ;
  while(cstring[j++]) ;
  j = j - 3;
  while(cstring[j] != ' ') j--;
  j++;
  sscanf((const char *__restrict__)(cstring+j),"%f",&freq);
  hz = freq*1000.0 + .5;  // MHz
  hz = hz * 1000000;      // Hz
  cycle = hz;
  cycle = 1.0 / cycle;    // seconds

  X86_cpuid( 0x0B, 0, regs ); threadpercore = regs[1] & 0xFFFF;
  if(threadpercore < 1) threadpercore = 1;
//   printf("0x0B, 0,EBX=%8.8x, ECX=%8.8x, eDX=%8.8x \n",regs[1],regs[2],regs[3]);
  X86_cpuid( 0x0B, 1, regs ); ncores = regs[1] & 0xFFFF; ncores /= threadpercore;
  if(ncores < 1) ncores = 1;
//   printf("0x0B, 1,EBX=%8.8x, ECX=%8.8x, eDX=%8.8x \n",regs[1],regs[2],regs[3]);
  X86_cpuid( 0x0B, 2, regs );
//   printf("0x0B, 2,EBX=%8.8x, ECX=%8.8x, eDX=%8.8x \n",regs[1],regs[2],regs[3]);

  if((ProcessorCapabilities & FLAG_FMA) == 0) return ; /* if FMA flag not present, AVX2 will not be */

  X86_cpuid( 7, 0, regs );    /* get more CPU capabilities EAX=7, ECX=0 */
  if((1 << 5) & regs[1]) ProcessorCapabilities |= FLAG_AVX2 ;   /* AVX2  EBX bit 5 */
  if((1 << 8) & regs[1]) ProcessorCapabilities |= FLAG_BMI  ;   /* BMI2  EBX bit 8 needed to set our BMI flag */
#else
  cstring = "Unknown processor";
#endif
}

/****f* X86/get_cpu_freq
 * FUNCTION
 *  get the CPU nominal frequency
 * SYNOPSIS
    C:
     uint64_t get_cpu_freq();

    Fortran:
     interface
       function get_cpu_freq() result(freq)
         import C_INT64_T
         integer(C_INT64_T) :: freq
       end function get_cpu_freq
     end interface
 * AUTHOR
 *  M.Valin Recherche en Prevision Numerique 2016
 * ARGUMENTS
    none
 * RESULT
 *  the CPU nominal frequency in Hertz (64 bit unsigned integer) 
 *  (nonzero on X86 family cpus only)
 * 
*****
*/
#pragma weak get_cpu_freq__=get_cpu_freq
#pragma weak get_cpu_freq_=get_cpu_freq
uint64_t get_cpu_freq__(void);
uint64_t get_cpu_freq_(void);
uint64_t get_cpu_freq(void)
{
  if(ProcessorCapabilities == 0) get_cpu_capabilities();
  return( hz );
}

/****f* X86/cpu_has_feature
 * FUNCTION
 *   determine whether certain capabilities are present on the current CPU
 *   e.g. AVX/AVX2/SSE2/FMA ...
 * SYNOPSIS
    C:
     int cpu_has_feature(int feature);

    Fortran:
     interface
       function cpu_has_feature(feature) result(status)
         import C_INT
         integer(C_INT), intent(IN), value :: feature
         integer(C_INT) :: status
       end function cpu_has_feature
     end interface
 * AUTHOR
 *  M.Valin Recherche en Prevision Numerique 2016
 * ARGUMENTS
    feature : feature symbol from cpu_type (#include <cpu_type.h>)
 * RESULT
 *  1 is feature is supported by CPU, 0 otherwise
 * NOTES
 * *  #include <cpu_type.h> 
 * *  is usable by C and Fortran programs alike
 * *  Fortran programs must use
 * *  #define IN_FORTRAN_CODE
 * *  before including cpu_type.h
*****
*/
#pragma weak cpu_has_feature__=cpu_has_feature
#pragma weak cpu_has_feature_=cpu_has_feature
int cpu_has_feature__(int flag);
int cpu_has_feature_(int flag);
int cpu_has_feature(int flag)
{
  if(ProcessorCapabilities == 0) get_cpu_capabilities();
  return( X86_FLAG(flag) );
}


/****f* X86/rdtscp
 * FUNCTION
 *   read the time stamp counter WITH serializing
 * SYNOPSIS
    C:
     uint64_t rdtscp();

    Fortran:
     interface
       function rdtscp() result(tsc)
         import C_INT64_T
         integer(C_INT64_T) :: tsc
       end function rdtscp
     end interface
 * AUTHOR
 *  M.Valin Recherche en Prevision Numerique 2016
 * ARGUMENTS
    none
 * RESULT
 *  the 64 bit unsigned integer value of the time stamp counter
 * NOTES
 *  slower than rdtsc but serializing is used to minimize impact of out of order execution
 *  (always 0 on non X86 family cpus)
 *
*****
*/
#pragma weak rdtscp__=rdtscp
#pragma weak rdtscp_=rdtscp
uint64_t rdtscp__(void);
uint64_t rdtscp_(void);
uint64_t rdtscp(void) {   // version "in order" avec "serialization"
  uint32_t lo, hi;
#if defined(__x86_64__)
  __asm__ volatile ("rdtscp"
      : /* outputs   */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers  */ "%rcx");
  __asm__ volatile ("mfence");
#else
  lo = 0;
  hi = 0;
#endif
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}

/****f* X86/rdtsc
 * FUNCTION
 *   read the time stamp counter WITHOUT serializing
 *  SYNOPSIS
   C:
     uint64_t rdtsc();

   Fortran:
     interface
       function rdtsc() result(tsc)
         import C_INT64_T
         integer(C_INT64_T) :: tsc
       end function rdtsc
     end interface
 * AUTHOR
 *  M.Valin Recherche en Prevision Numerique 2016
 * ARGUMENTS
    none
 * RESULT
 *  the 64 bit unsigned integer value of the time stamp counter
 * NOTES
 *  faster than rdtscp but out of order execution might affect result as no serializing is done
 *  (always 0 on non X86 family cpus)
 *
*****
*/
#pragma weak rdtsc__=rdtsc
#pragma weak rdtsc_=rdtsc
uint64_t rdtsc__(void);
uint64_t rdtsc_(void);
uint64_t rdtsc(void) {   // version rapide "out of order"
  uint32_t lo, hi;
#if defined(__x86_64__)
  __asm__ volatile ("rdtsc"
      : /* outputs   */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers  */ "%rcx");
#else
  lo = 0;
  hi = 0;
#endif
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}

/****f* X86/wall_clock_seconds
 * FUNCTION
 *   convert CPU clock ticks into seconds
 *  SYNOPSIS
   C:
     double wall_clock_seconds(uint64_t ticks);

   Fortran:
     interface
       function wall_clock_seconds(ticks) result(seconds)
         import C_INT64_T, C_DOUBLE
         integer(C_INT64_T), intent(IN), value :: ticks
         real(C_DOUBLE) :: seconds
       end function wall_clock_seconds
     end interface
 * AUTHOR
 *  M.Valin Recherche en Prevision Numerique 2016
 * ARGUMENTS
    ticks  : cpu ticks (from rdtsc/rdtscp)
 * RESULT
 *  ticks converted into seconds
 *
*****
*/
#pragma weak wall_clock_seconds__=wall_clock_seconds
#pragma weak wall_clock_seconds_=wall_clock_seconds
double wall_clock_seconds__(uint64_t ticks) ;
double wall_clock_seconds_(uint64_t ticks) ;
double wall_clock_seconds(uint64_t ticks) {
  double t=ticks;
  if(ProcessorCapabilities == 0) get_cpu_capabilities();
  return(t*cycle);
}

/****f* X86/rdtscp_seconds
 * FUNCTION
 *   get elapsed time since last call to rdtscp_seconds or rdtsc_seconds
 *   WITH serializing
 *  SYNOPSIS
   C:
     double rdtscp_seconds();

   Fortran:
     interface
       function rdtscp_seconds() result(seconds)
         import C_DOUBLE
         real(C_DOUBLE) :: seconds
       end function rdtscp_seconds
     end interface
 * AUTHOR
 *  M.Valin Recherche en Prevision Numerique 2016
 * ARGUMENTS
    none
 * RESULT
 *  the elapsed time in seconds since the last call to rdtsc_seconds or rdtscp_seconds
 * NOTES
 *  slower than rdtsc_seconds but serializing is used to minimize impact of out of order execution
 *
*****
*/
#pragma weak rdtscp_seconds__=rdtscp_seconds
#pragma weak rdtscp_seconds_=rdtscp_seconds
double rdtscp_seconds__(void) ;
double rdtscp_seconds_(void) ;
double rdtscp_seconds(void) {
  double t;
  uint64_t now;
  now = rdtscp();
  t = now -last_time;
  last_time = now;
  if(ProcessorCapabilities == 0) get_cpu_capabilities();
  return(t*cycle);
}

/****f* X86/rdtsc_seconds
 * FUNCTION
 *   get elapsed time since last call to rdtscp_seconds or rdtsc_seconds
 *   WITHOUT serializing
 *  SYNOPSIS
   C:
     double rdtsc_seconds();

   Fortran:
     interface
       function rdtsc_seconds() result(seconds)
         import C_DOUBLE
         real(C_DOUBLE) :: seconds
       end function rdtsc_seconds
     end interface
 * AUTHOR
 *  M.Valin Recherche en Prevision Numerique 2016
 * ARGUMENTS
    none
 * RESULT
 *  the elapsed time in seconds since the last call to rdtsc_seconds or rdtscp_seconds
 * NOTES
 *  faster than rdtscp_seconds but out of order execution might affect result as no serializing is done
 *
*****
*/
#pragma weak rdtsc_seconds__=rdtsc_seconds
#pragma weak rdtsc_seconds_=rdtsc_seconds
double rdtsc_seconds__(void) ;
double rdtsc_seconds_(void) ;
double rdtsc_seconds(void) {
  double t;
  uint64_t now;
  now = rdtsc();
  t = now -last_time;
  last_time = now;
  if(ProcessorCapabilities == 0) get_cpu_capabilities();
  return(t*cycle);
}

/****f* X86/get_fp_status_ctl
 * FUNCTION
 *   get floating point status word
 *  SYNOPSIS
   C:
     int get_fp_status_ctl();

   Fortran:
     interface
       function get_fp_status_ctl() result(id)
         import C_INT
         integer(C_INT) :: id
       end function get_fp_status_ctl
     end interface
 * AUTHOR
 *  M.Valin Recherche en Prevision Numerique 2016
 * ARGUMENTS
    none
 * RESULT
 *  the floating point status word (nozero on X86 family cpus only)
*****
*/
#pragma weak get_fp_status_ctl__= get_fp_status_ctl
#pragma weak get_fp_status_ctl_= get_fp_status_ctl
int get_fp_status_ctl__(void);
int get_fp_status_ctl_(void);
int get_fp_status_ctl(void) {
  int fpstat = 0;
#if defined(__i386__) || defined(__x86_64__)
  __asm__ volatile ("stmxcsr %0" : "=m" (fpstat));
#endif
  return (fpstat);
}

/****f* X86/set_fp_status_ctl
 * FUNCTION
 *   set the floating point status and control word
 *  SYNOPSIS
   C:
     void set_fp_status_ctl(int id);

   Fortran:
     interface
       subroutine set_fp_status_ctl(id)
         import C_INT
         integer(C_INT), intent(IN), value :: id
       end subroutine set_fp_status_ctl
     end interface
 * AUTHOR
 *  M.Valin Recherche en Prevision Numerique 2016
 * ARGUMENTS
    id : integer value to store into the floating point status and control word
 * NOTES
 *  (on X86 family cpus only)
*****
*/
#pragma weak set_fp_status_ctl__=set_fp_status_ctl
#pragma weak set_fp_status_ctl_=set_fp_status_ctl
void set_fp_status_ctl__(int fpstat_in);
void set_fp_status_ctl_(int fpstat_in);
void set_fp_status_ctl(int fpstat_in) {
  int fpstat = fpstat_in;
#if defined(__i386__) || defined(__x86_64__)
  __asm__ volatile ("ldmxcsr %0" : "=m" (fpstat));
#endif
  return ;
}

#if defined(TEST_CPUID)
int get_cpu_core_thread()  /* Intel CPUs only and even in this case not always reliable */
{
  uint32_t regs[4];
  if(ProcessorCapabilities == 0) get_cpu_capabilities();
  X86_cpuid( 0x0B, 0, regs );
  return( ((regs[3]>>regs[0]) << 8) | (regs[3]-((regs[3]>>regs[0])<<regs[0])) );  /* core << 8 + thread  */
}
#endif

/****f* X86/get_cpu_hyperthreads
 * FUNCTION
 *   get the number of hyperthreads of this CPU
 *  SYNOPSIS
   C:
     int get_cpu_hyperthreads();

   Fortran:
     interface
       function get_cpu_hyperthreads() result(nhyperthreads)
         import C_INT
         integer(C_INT) :: nhyperthreads
       end function get_cpu_hyperthreads
     end interface
 * AUTHOR
 *  M.Valin Recherche en Prevision Numerique 2016
 * ARGUMENTS
    none
 * RESULT
 *  the number of hyperthreads of this CPU (not 1 on X86 family cpus only)
*****
*/
#pragma weak get_cpu_hyperthreads__=get_cpu_hyperthreads
#pragma weak get_cpu_hyperthreads_=get_cpu_hyperthreads
int get_cpu_hyperthreads__();
int get_cpu_hyperthreads_();
int get_cpu_hyperthreads()  /* Intel CPUs only */
{
  if(ProcessorCapabilities == 0) get_cpu_capabilities();
  return( threadpercore );
}

/****f* X86/get_cpu_cores
 * FUNCTION
 *   get the number of cores of this CPU
 *  SYNOPSIS
   C:
     int get_cpu_cores();

   Fortran:
     interface
       function get_cpu_cores() result(ncores)
         import C_INT
         integer(C_INT) :: ncores
       end function get_cpu_cores
     end interface
 * AUTHOR
 *  M.Valin Recherche en Prevision Numerique 2016
 * ARGUMENTS
    none
 * RESULT
 *  the number of cores of this CPU (not one on X86 family cpus only)
*****
*/
#pragma weak get_cpu_cores__=get_cpu_cores
#pragma weak get_cpu_cores_=get_cpu_cores
int get_cpu_cores__();
int get_cpu_cores_();
int get_cpu_cores()  /* Intel CPUs only */
{
  if(ProcessorCapabilities == 0) get_cpu_capabilities();
  return( ncores );
}

/****f* X86/get_cpu_id
 * FUNCTION
 *   get the unique ID of this CPU in the node
 *  SYNOPSIS
   C:
     int get_cpu_id();

   Fortran:
     interface
       function get_cpu_id() result(id)
         import C_INT
         integer(C_INT) :: id
       end function get_cpu_id
     end interface
 * AUTHOR
 *  M.Valin Recherche en Prevision Numerique 2016
 * ARGUMENTS
    none
 * RESULT
 *  the APIC id of the current cpu (nozero on X86 family cpus only)
 * NOTES
 *  to be correlated with /proc/cpuinfo to find actual CPU number
*****
*/
#pragma weak get_cpu_id__=get_cpu_id
#pragma weak get_cpu_id_=get_cpu_id
int get_cpu_id__();
int get_cpu_id_();
int get_cpu_id()  /* Intel CPUs only */
{
  uint32_t regs[4];
  X86_cpuid( 0x0B, 0, regs );
  return( regs[3] );  /* x2APIC id from EDX  */
}

#if defined(TEST_CPUID)
int main_cpuid(int argc, char** argv)
{
  int core_and_thread;
  int status;
  float r;
  uint32_t regs[4];
#if defined(DEBUG)
#endif
  printf("CPU with %d cores and %d threads/core \n",get_cpu_cores(),get_cpu_hyperthreads());
  printf("FPU status = %8.8x\n",get_fp_status_ctl());
  core_and_thread = get_cpu_core_thread();
  printf("core = %d, thread = %d\n",core_and_thread>>8,core_and_thread&0xFF);
#if defined(DEBUG)
#endif
  printf("CPU speed: %lu Hz\n",get_cpu_freq());
  printf("CPU cycle = %g sec\n",cycle);
  printf("FLAGS: ");
  if(cpu_has_feature(FLAG_SSE))  printf(" SSE");
  if(cpu_has_feature(FLAG_SSE2)) printf(" SSE2");
  if(cpu_has_feature(FLAG_SSE3)) printf(" SSE3");
  if(cpu_has_feature(FLAG_SSE4)) printf(" SSE4");
  if(cpu_has_feature(FLAG_AVX))  printf(" AVX");
  if(cpu_has_feature(FLAG_FMA))  printf(" FMA");
  if(cpu_has_feature(FLAG_AVX2)) printf(" AVX2");
  printf("\n");
  printf("CPU Version string: '%s'\n",cstring);
  status = get_fp_status_ctl();
  printf("FPU status = %8.8x\n",status);
  printf("Precision status %s\n",status & FP_STATUS_PE ? "ON" : "OFF");
  printf("Underflow status %s\n",status & FP_STATUS_UE ? "ON" : "OFF");
  printf("Overflow status %s\n",status & FP_STATUS_OE ? "ON" : "OFF");
  printf("Zero divide status %s\n",status & FP_STATUS_ZE ? "ON" : "OFF");
  printf("forcing underflow, zero divide, overflow\n");
  r = 1.0E-30; r = r*r; /* underflow */
  r = 1.0/r;            /* zero divide */
  r = 1.0E30; r = r*r;  /* overflow */
  status = get_fp_status_ctl();
  printf("FPU status = %8.8x\n",status);
  printf("Precision status %s\n",status & FP_STATUS_PE ? "ON" : "OFF");
  printf("Underflow status %s\n",status & FP_STATUS_UE ? "ON" : "OFF");
  printf("Overflow status %s\n",status & FP_STATUS_OE ? "ON" : "OFF");
  printf("Zero divide status %s\n",status & FP_STATUS_ZE ? "ON" : "OFF");
  return (0);
}
#endif

#if defined(TEST_RDTSC)
int main_rdtsc(int argc, char** argv){   // X86_64 timing demo:  cc -DSELF_TEST rdtscp.c ; ./a.out
  uint64_t t1;
  uint64_t t2;
  int i;
  int ta, tb, tc;
  float xc;
  double s1, s2, s3;
  t1=rdtscp();
  sleep(1);
  t2=rdtscp();
  t2=(t2-t1)/1000000;  // combien de megaticks dans une seconde
  tc=t2;
  xc=tc*.001;
  for(i=0;i<5;i++){
    t1=rdtscp();
    t2=rdtscp();
    ta=t2-t1;
    if(xc != 0) ta/=xc;   // conversion en nanosecondes
    t1=rdtscp();
    sleep(1);
    t2=rdtscp();
    t2=(t2-t1)/1000000;
    tc=t2;
    t1=rdtsc();
    t2=rdtsc();
    tb=t2-t1;
    if(xc != 0) tb/=xc;   // conversion en nanosecondes
    printf("%d ns, %d ns, %d ticks/us, cycles = %lu, sec = %E\n",ta,tb,tc,t2-t1,wall_clock_seconds(t2-t1));
  }
  for(i=0;i<5;i++){
    s1 = rdtsc_seconds();
    s2 = rdtsc_seconds();
    s1 = rdtscp_seconds();
    s3 = rdtscp_seconds();
    printf("rdtsc(p)_seconds overhead = %E sec %E sec\n",s2,s3);
  }
  return(0);
}
#endif
#if defined(TEST_CPUID) || defined(TEST_RDTSC)
int main(int argc, char** argv){
  int status = 0;
#if defined(TEST_RDTSC)
  status = main_rdtsc(argc,argv);
  if(status) return(status);
#endif
#if defined(TEST_CPUID)
  status = main_cpuid(argc,argv);
#endif
  return(status);
}
#endif
