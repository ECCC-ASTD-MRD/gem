! local times array post processing code for RPN_COMM package
! will insert T1, T2, T3 into times array
      if(ntimes+2 <= ntimes_sz) then  ! timing array is not full
        times(ntimes)   = T1
        times(ntimes+1) = T2
        times(ntimes+2) = T3
        ntimes = ntimes + 3
        times(0) = ntimes
      endif


