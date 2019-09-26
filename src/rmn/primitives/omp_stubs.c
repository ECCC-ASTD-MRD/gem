//! Functions stubs for when OpenMP is not enabled to allow the code to compile and work
#ifndef _OPENMP
   // Define OpenMP function stubs so that GEM can build and
   // work, with reduced non-essential functionnality.

   // Create function aliases to be callable from Fortran

   #include <time.h>
   #include <stdio.h>
   #include <stdlib.h>

   //! Provide the elapsed real time as a double
   double omp_get_wtime() {
      struct timespec ts;

      if ( clock_gettime(CLOCK_REALTIME, &ts) != 0 ) {
         fprintf(stderr, "clock_gettime(CLOCK_REALTIME, &ts) system call failed!\n");
         exit(1);
      }
      return ts.tv_sec + ts.tv_nsec * 1.0e-9;
   }
   void omp_get_wtime_() __attribute__((alias ("omp_get_wtime")));


   //! Get the number of OpenMP threads
   int omp_get_num_threads() {
      return 1;
   }
   void omp_get_num_threads_() __attribute__((alias ("omp_get_num_threads")));


   //! Get the maximum number of OpenMP threads
   int omp_get_max_threads() {
      return 1;
   }
   void omp_get_max_threads_() __attribute__((alias ("omp_get_max_threads")));


   //! Get the number of current OpenMP thread
   int omp_get_thread_num() {
      return 1;
   }
   void omp_get_thread_num_() __attribute__((alias ("omp_get_thread_num")));


   //! Pretend setting the number of OpenMP threads
   void omp_set_num_threads(int dummy) {
      return;
   }
   void omp_set_num_threads_() __attribute__((alias ("omp_set_num_threads")));
#endif
