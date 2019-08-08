program SpeedTest
  implicit none 

  integer, parameter :: NNI=500         ! scales the arrays involved
  real,    dimension(NNI,50) :: vLevelSource, vLevelDestn, &
                                stateIn, derivIn, &
                                stateOut, derivOut
  integer, dimension(NNI,50) :: posnDestInSrc

  real :: extrapGuideDown, extrapGuideUp

  integer :: h                          ! horizontal loop index
  integer :: v                          ! vertical loop index

  integer :: clock_count_start, clock_count_end, count_rate
  real :: time_sec                      ! time, in seconds


  ! Fill in some source and destination levels, and source states and derivatives
  do v=1,50
    do h=1,NNI
      vLevelSource(h,v)=v
      vLevelDestn (h,v)=v + 0.2
      stateIn(h,v)=sin(v * 0.01)
      derivIn(h,v)=cos(v * 0.01)
    end do ! h
  end do ! v

  ! Find the destination positions in the source levels
  call Interp1D_FindPos(NNI, 50, 50, &
                        NNI, NNI, &
                        vLevelSource, &
                        posnDestInSrc, &
                        vLevelDestn &
                       )

  ! Process 200,000 horizontal points
  
  extrapGuideDown=0.              ! extrapolation not used yet
  extrapGuideUp=0.

  call system_clock(clock_count_start, count_rate)

  do h=1,200000/NNI
    call Extrap1D_LapseRate(NNI, 50, 50, &
                            NNI, NNI, &

                            vLevelSource, &
                            stateIn, &
                            derivIn, &

                            posnDestInSrc, &
                            vLevelDestn, &
                            stateOut, &
                            derivOut, &

                            extrapGuideDown, extrapGuideUp, &
                            .true., .true. &
                           )
  end do

  call system_clock(clock_count_end)
  time_sec=(clock_count_end - clock_count_start)
  time_sec=time_sec / count_rate
  write(6,*)'count_rate=', count_rate
  write(6,*)'elapsed time (in s) = ', time_sec

end program SpeedTest
