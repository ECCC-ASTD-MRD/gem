      subroutine clock ( unf, from, last )
      use, intrinsic :: iso_fortran_env
      implicit none

      ! Unit file for output
      integer, intent(in) :: unf
      ! Does something special if the value is 'CURRENT TIMESTEP'
      character(len=*), intent(in) :: from
      ! If true, print the accumulated time
      logical, intent(in) :: last

      integer, external :: get_max_rss

      character(5)  :: zone
      character date_S*8, time_S*10, jour_S*11
      logical, save :: timini_L = .false.
      real          users, systs
      real, save :: user0 = 0.0, syst0 = 0.0
      real(kind=REAL64), save ::  START = -1.d0, END, avgtime(10) = 0.d0, &
                                 ACCUM_w = 0.d0, ACCUM_u = 0.d0, ACCUM_s = 0.d0
      real(kind=REAL64) omp_get_wtime

!----------------------------------------------------------------

      if (unf <= 0) return

      if (START < 0.d0) START = omp_get_wtime()

      call date_and_time( date_S,time_S )
      jour_S  = date_S(1:4)//'/'//date_S(5:6)//'/'//date_S(7:8)//' '

      call cpu_time( users )

      call cpu_time( systs )

      END = omp_get_wtime()

      if (timini_L) then
         if (trim(from) =='CURRENT TIMESTEP') then
            avgtime(1:9) = avgtime(2:10)
            avgtime(10)  = END-START
         endif

         write(unf,1000) 'TIME: '//jour_S//time_S, END-START, users-user0, &
                          systs-syst0, get_max_rss(), from
         ACCUM_w = ACCUM_w + END - START
         ACCUM_u = ACCUM_u + users - user0
         ACCUM_s = ACCUM_s + systs - syst0
         if (last) write(unf,1001) ACCUM_w,ACCUM_u,ACCUM_s
      endif

      user0 = users
      syst0 = systs
      START = END

      timini_L = .true.

!----------------------------------------------------------------

 1000    format(/A,' W: ',1pe13.6, &
                   ' U: ',1pe13.6, &
                   ' S: ',1pe13.6, &
                  ', Mem: ',i8,' (Kbytes/PE) ',a)
 1001    format(/'ACCUMULATED TIME: W: ',1pe13.6, &
                   ' U: ',1pe13.6, &
                   ' S: ',1pe13.6)

      return
      end subroutine clock
