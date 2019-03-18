interface
  subroutine tricublin_zyxf_beta(d, f, px, py, pz, ni, ninj,zcubic_L) bind(C,name='Tricublin_zyxf_beta')
    import :: C_DOUBLE, C_FLOAT, C_INT, C_BOOL
    logical(C_BOOL) , value ::  zcubic_L
    real(C_FLOAT), intent(OUT) :: d
    real(C_FLOAT), dimension(*), intent(IN) :: f
    real(C_DOUBLE), dimension(4), intent(IN) :: px, py, pz
    integer(C_INT), intent(IN), value :: ni, ninj
  end subroutine tricublin_zyxf_beta
  subroutine tricubic_coeffs_d(px, py, pz, x, y, z) bind(C,name='Tricubic_coeffs_d')
    import C_DOUBLE
    real(C_DOUBLE), dimension(4), intent(OUT) :: px, py, pz
    real(C_DOUBLE), intent(IN), value :: x, y, z
  end subroutine tricubic_coeffs_d
  subroutine tricublin_zyx3f_beta(d, f, px, py, pz, ni, ninj, zcubic_L) bind(C,name='Tricublin_zyx3f_beta')
    import :: C_DOUBLE, C_FLOAT, C_INT, C_BOOL
    logical(C_BOOL), value ::  zcubic_L
    real(C_FLOAT), dimension(3), intent(OUT) :: d
    real(C_FLOAT), dimension(*), intent(IN) :: f
    real(C_DOUBLE), dimension(4), intent(IN) :: px, py, pz
    integer(C_INT), intent(IN), value :: ni, ninj
  end subroutine tricublin_zyx3f_beta
end interface
