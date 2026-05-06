subroutine cpl_declare_exchanged_fields()
    use iris_mod
    use iso_c_binding
    use iso_fortran_env
    use rmn_gmm ! Because cplocn.cdk uses a bunch of values from gmm_mod
    use cplocn_mod
    implicit none

    integer :: ivar
    integer :: ierr
    character(len=5) :: nomvar
    real(C_FLOAT) :: level
    integer :: my_pe, my_pe_x, my_pe_y

    call rpn_comm_mype(my_pe, my_pe_x, my_pe_y)
    if(my_pe /= 0) then
        return
    endif

    level = 0.0
    write(error_unit,*) "NEMO calling Iris%Model_Provide for variable ", cplocn_iris_provides
    ierr = iris%model_provide(iris_grid, cplocn_iris_provides, level)
    if(ierr /= 0) then
        write(error_unit,*) "NEMO: Iris%Model_Provide level=", level
        return
    endif

    write(error_unit,*) "NEMO: Calling Iris%Model_Consume for variable ", cplocn_iris_consumes
    ierr = iris%model_consume(iris_grid, cplocn_iris_consumes, level)
    if(ierr /= 0) then
        write(error_unit,*) "NEMO: Iris%Model_Provide level=", level
        return
    endif
end subroutine
