!
!	C callable interface to FORTRAN library openMP routines/functions
!	needed within the RPN_COMM library
!
!	the purpose of this is to have predictable names to be called from C
!	and make it FORTRAN callable no matter what the binding is
!
!	no underscore version  (C style)
!
	function f_omp_get_max_threads() result(r)&
     &                 bind(C,name='f_omp_get_max_threads')
	use iso_c_binding
	implicit none
        integer(C_INT) :: r
	external :: omp_get_max_threads
	integer ::  omp_get_max_threads

	r=omp_get_max_threads()
	end function f_omp_get_max_threads
!
	subroutine f_omp_set_num_threads(nthreads)&
     &                 bind(C,name='f_omp_set_num_threads')
	use iso_c_binding
	implicit none
	integer(C_INT), intent(IN) :: nthreads
	call omp_set_num_threads(nthreads)
	return
	end subroutine f_omp_set_num_threads
!
!	one underscore version
!
	function f_omp_get_max_threads_() result(r)&
     &                 bind(C,name='f_omp_get_max_threads_')
	use iso_c_binding
	implicit none
        integer(C_INT) :: r
	external :: omp_get_max_threads
	integer ::  omp_get_max_threads

	r=omp_get_max_threads()
	end function f_omp_get_max_threads_
!
	subroutine f_omp_set_num_threads_(nthreads)&
     &                 bind(C,name='f_omp_set_num_threads_')
	use iso_c_binding
	implicit none
	integer(C_INT), intent(IN) :: nthreads
	call omp_set_num_threads(nthreads)
	return
	end subroutine f_omp_set_num_threads_
!
!	two underscores version
!
	function f_omp_get_max_threads__() result(r)&
     &                 bind(C,name='f_omp_get_max_threads__')
	use iso_c_binding
	implicit none
        integer(C_INT) :: r
	external :: omp_get_max_threads
	integer ::  omp_get_max_threads

	r=omp_get_max_threads()
	end function f_omp_get_max_threads__
!
	subroutine f_omp_set_num_threads__(nthreads)&
     &                 bind(C,name='f_omp_set_num_threads__')
	use iso_c_binding
	implicit none
	integer(C_INT), intent(IN) :: nthreads
	call omp_set_num_threads(nthreads)
	return
	end subroutine f_omp_set_num_threads__
