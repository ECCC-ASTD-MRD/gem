/* capability flags used by cpu_type.c */
#define FLAG_SSE    1
#define FLAG_SSE2   2
#define FLAG_SSE3   4
#define FLAG_SSE4   8
#define FLAG_AVX   16
#define FLAG_AVX2  32
#define FLAG_FMA   64
#define FLAG_BMI  128

#define FP_STATUS_RSHFT 13  
#define FP_STATUS_RNR    0
#define FP_STATUS_RDN    1
#define FP_STATUS_RUP    2
#define FP_STATUS_RZR    3

#define FP_STATUS_PM  4096
#define FP_STATUS_UM  2048
#define FP_STATUS_OM  1024
#define FP_STATUS_ZM   512
#define FP_STATUS_DM   256
#define FP_STATUS_IM   128
#define FP_STATUS_PE    32
#define FP_STATUS_UE    16
#define FP_STATUS_OE     8
#define FP_STATUS_ZE     4
#define FP_STATUS_DE     2
#define FP_STATUS_IE     1

#if defined(IN_FORTRAN_CODE)
  interface
    function cpu_has_feature(feature) result(status)
      import :: C_INT
      integer(C_INT), intent(IN), value :: feature
      integer(C_INT) :: status
    end function cpu_has_feature

    function get_cpu_id() result(id)
      import :: C_INT
      integer(C_INT) :: id
    end function get_cpu_id

    function get_cpu_freq() result(freq)
      import C_INT64_T
      integer(C_INT64_T) :: freq
    end function get_cpu_freq

    function get_cpu_cores() result(ncores)
      import C_INT
      integer(C_INT) :: ncores
    end function get_cpu_cores

    function get_cpu_hyperthreads() result(nhyperthreads)
      import C_INT
      integer(C_INT) :: nhyperthreads
    end function get_cpu_hyperthreads

    function rdtsc() result(tsc)
      import C_INT64_T
      integer(C_INT64_T) :: tsc
    end function rdtsc

    function rdtscp() result(tsc)
      import C_INT64_T
      integer(C_INT64_T) :: tsc
    end function rdtscp

    function wall_clock_seconds(ticks) result(seconds)
      import C_INT64_T, C_DOUBLE
      integer(C_INT64_T), intent(IN), value :: ticks
      real(C_DOUBLE) :: seconds
    end function wall_clock_seconds

    function rdtscp_seconds() result(seconds)
      import C_DOUBLE
      real(C_DOUBLE) :: seconds
    end function rdtscp_seconds

    function rdtsc_seconds() result(seconds)
      import C_DOUBLE
      real(C_DOUBLE) :: seconds
    end function rdtsc_seconds

    function get_fp_status_ctl() result(id)
      import C_INT
      integer(C_INT) :: id
    end function get_fp_status_ctl

    subroutine set_fp_status_ctl(id)
      import C_INT
      integer(C_INT), intent(IN), value :: id
    end subroutine set_fp_status_ctl
  end interface
#else
  int cpu_has_feature(int feature);
  int get_cpu_cores();
  int get_cpu_hyperthreads();
  int get_cpu_id(void);
  uint64_t get_cpu_freq(void);
  uint64_t rdtsc(void);
  uint64_t rdtscp(void);
  double wall_clock_seconds(uint64_t ticks);
  double rdtsc_seconds(void);
  double rdtscp_seconds(void);
  void set_fp_status_ctl(int fpstat);
  int get_fp_status_ctl(void);
#endif
