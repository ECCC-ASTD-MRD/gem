        subroutine rpn_comm_test_003
        implicit none
        external :: f_omp_get_max_threads_
        integer ::  f_omp_get_max_threads_
        integer junk

        call f_omp_set_num_threads_(4)
        junk=f_omp_get_max_threads_()
        print *,junk
        stop
        end

