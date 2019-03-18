        external RPN_COMM_Wtime
        real *8 RPN_COMM_Wtime
        integer, save :: ntimes=0
        integer, save :: ntimes_sz=0
        real times(0:*)
        pointer(times_,times)
        save times_
        real *8 T0        ! T0 time of entry, reference base time for T1, T2,T3
        real T1,T2,T3     ! T1-T0 after 1st memory move, T2-T0 after MPI, T3-T0 after 2nd memory move
