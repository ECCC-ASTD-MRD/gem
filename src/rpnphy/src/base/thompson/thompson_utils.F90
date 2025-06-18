MODULE thompson_utils

  use machine , only : kind_phys

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: print_state

CONTAINS

  !------------------------------------------------------------------
  ! print_state
  !
  ! Prints statistics for the kernel state variables
  !------------------------------------------------------------------
  SUBROUTINE print_state(msg,   &
       ext_ndiag3d,      &
       spechum,          &
       qc,               &
       qr,               &
       qi,               &
       qs,               &
       qg,               &
       ni,               &
       nr,               &
       nc,               &
       nwfa,             &
       nifa,             &
       nwfa2d,           &
       nifa2d,           &
       tgrs,             &
       prsl,             &
       phil,             &
       area,             &
       diag3d,           &
       omega,            &
       prcp,             &
       rain,             &
       graupel,          &
       ice,              &
       snow,             &
       sr,               &
       refl_10cm,        &
       phii,             &
       spp_wts_mp,       &
       spp_prt_list,     &
       spp_stddev_cutoff &
       )

    INTEGER :: i
    CHARACTER(LEN=9) :: varn

    CHARACTER(LEN=*) :: msg
    INTEGER ext_ndiag3d
    REAL(4) spechum(:,:)
    REAL(4) qc(:,:)
    REAL(4) qr(:,:)
    REAL(4) qi(:,:)
    REAL(4) qs(:,:)
    REAL(4) qg(:,:)
    REAL(4) ni(:,:)
    REAL(4) nr(:,:)
    REAL(4) nc(:,:)
    REAL(4) nwfa(:,:)
    REAL(4) nifa(:,:)
    REAL(4) nwfa2d(:)
    REAL(4) nifa2d(:)
    REAL(4) tgrs(:,:)
    REAL(4) prsl(:,:)
    REAL(4) phil(:,:)
    REAL(4) area(:)
    REAL(4) diag3d(:,:,:)
    REAL(4) omega(:,:)
    REAL(4) prcp(:)
    REAL(4) rain(:)
    REAL(4) graupel(:)
    REAL(4) ice(:)
    REAL(4) snow(:)
    REAL(4) sr(:)
    REAL(4) refl_10cm(:,:)
    REAL(4) phii(:,:)
    REAL(4) spp_wts_mp(:,:)
    REAL(4) spp_prt_list(:)
    REAL(4) spp_stddev_cutoff(:)

    WRITE(*,'(A4)') "TEST"
    WRITE(*,'(A5,A137)') "TEST ", REPEAT("=",137)
    WRITE(*,'(A5,A32)') "TEST ", msg
    WRITE(*,'(A5,A137)') "TEST ", REPEAT("=",137)
    WRITE(*,'(A5,A17,6A20)') "TEST ", "Variable", "Min", "Max", "Avg", "First", "Last", "RMS"
    WRITE(*,'(A5,A137)') "TEST ", REPEAT("-",137)

    CALL print_2d_variable("spechum", spechum)
    CALL print_2d_variable("qc", qc)
    CALL print_2d_variable("qr", qr)
    CALL print_2d_variable("qi", qi)
    CALL print_2d_variable("qs", qs)
    CALL print_2d_variable("qg", qg)
    CALL print_2d_variable("ni", ni)
    CALL print_2d_variable("nr", nr)
    CALL print_2d_variable("nc", nc)
    CALL print_2d_variable("nwfa", nwfa)
    CALL print_2d_variable("nifa", nifa)
    CALL print_1d_variable("nwfa2d", nwfa2d)
    CALL print_1d_variable("nifa2d", nifa2d)
    CALL print_2d_variable("tgrs", tgrs)
    CALL print_2d_variable("prsl", prsl)
    CALL print_2d_variable("phil", phil)
    CALL print_1d_variable("area", area)
    DO i = 1, ext_ndiag3d
        WRITE(varn, '(A,I0.2)') "diag3d_", i
        CALL print_2d_variable(varn, diag3d(:,:,i))
    ENDDO
    CALL print_2d_variable("omega", omega)
    CALL print_1d_variable("prcp", prcp)
    CALL print_1d_variable("rain", rain)
    CALL print_1d_variable("graupel", graupel)
    CALL print_1d_variable("ice", ice)
    CALL print_1d_variable("snow", snow)
    CALL print_1d_variable("sr", sr)
    CALL print_2d_variable("refl_10cm", refl_10cm)
    CALL print_2d_variable("phii", phii)
    CALL print_2d_variable("spp_wts_mp", spp_wts_mp)
    CALL print_1d_variable("spp_prt_list", spp_prt_list)
    CALL print_1d_variable("spp_stddev_cutoff", spp_stddev_cutoff)

    WRITE(*,'(A5,A137)') "TEST ", REPEAT("-",137)
    WRITE(*,'(A4)') "TEST"

  END SUBROUTINE print_state

  !------------------------------------------------------------------
  ! print_1d_variable
  !
  ! Prints statistics for a 1d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_1d_variable(name, data)

    CHARACTER(LEN=*) :: name
    REAL(4)         :: data(:), avg

    ! Note: Assumed shape array sections always start with index=1 for all
    ! dimensions
    !       So we don't have to know start/end indices here
    avg = SUM(data) / SIZE(data)
    WRITE(*,'(A5, A17,6ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), avg, data(1), &
                            data(SIZE(data,1)),            &
                            SQRT(SUM(data**2 - avg**2) / SIZE(data))

  END SUBROUTINE print_1d_variable

  !------------------------------------------------------------------
  ! print_2d_variable
  !
  ! Prints statistics for a 2d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_2d_variable(name, data)

    CHARACTER(LEN=*) :: name
    REAL(4)         :: data(:,:), avg

    ! Note: Assumed shape array sections always start with index=1 for all
    ! dimensions
    !       So we don't have to know start/end indices here
    avg = SUM(data) / SIZE(data)
    WRITE(*,'(A5, A17, 6ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), avg, data(1,1), &
                            data(SIZE(data,1), SIZE(data,2)),            &
                            SQRT(SUM(data**2 - avg**2) / SIZE(data))

  END SUBROUTINE print_2d_variable

  !------------------------------------------------------------------
  ! print_3d_variable
  !
  ! Prints statistics for a 3d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_3d_variable(name, data)

    CHARACTER(LEN=*) :: name
    REAL(4)         :: data(:,:,:)

    ! Note: Assumed shape array sections always start with index=1 for all dimensions
    !       So we do not have to know start/end indices here
    WRITE(*,'(A5,A17,5ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), data(1,1,1),  &
                            data(SIZE(data,1), SIZE(data,2), SIZE(data,3)), &
                            SQRT(SUM(data**2) / SIZE(data))

  END SUBROUTINE print_3d_variable

  !------------------------------------------------------------------
  ! print_4d_variable
  !
  ! Prints statistics for a 4d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_4d_variable(name, data)

    CHARACTER(LEN=*) :: name
    REAL(4)         :: data(:,:,:,:)

    ! Note: Assumed shape array sections always start with index=1 for all dimensions
    !       So we do not have to know start/end indices here
    WRITE(*,'(A5,A17,5ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), data(1,1,1,1),  &
                            data(SIZE(data,1), SIZE(data,2), SIZE(data,3), SIZE(data,4)), &
                            SQRT(SUM(data**2) / SIZE(data))

  END SUBROUTINE print_4d_variable


  !------------------------------------------------------------------
  ! print_5d_variable
  !
  ! Prints statistics for a 5d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_5d_variable(name, data)

    CHARACTER(LEN=*) :: name
    REAL(4)         :: data(:,:,:,:,:)

    ! Note: Assumed shape array sections always start with index=1 for all dimensions
    !       So we do not have to know start/end indices here
    WRITE(*,'(A5,A17,5ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), data(1,1,1,1,1),  &
                            data(SIZE(data,1), SIZE(data,2), SIZE(data,3), SIZE(data,4), SIZE(data,5)), &
                            SQRT(SUM(data**2) / SIZE(data))

  END SUBROUTINE print_5d_variable

  !------------------------------------------------------------------
  ! print_1d_variable
  !
  ! Prints statistics for a 1d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_1d_variable_int(name, data)

    CHARACTER(LEN=*) :: name
    INTEGER         :: data(:)

    ! Note: Assumed shape array sections always start with index=1 for all
    ! dimensions
    !       So we don't have to know start/end indices here
    WRITE(*,'(A5, A17,4I20,ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), data(1), &
                            data(SIZE(data,1)),            &
                            SQRT(REAL(SUM(data**2) / SIZE(data)))

  END SUBROUTINE print_1d_variable_int

  !------------------------------------------------------------------
  ! print_2d_variable
  !
  ! Prints statistics for a 2d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_2d_variable_int(name, data)

    CHARACTER(LEN=*) :: name
    INTEGER         :: data(:,:)

    ! Note: Assumed shape array sections always start with index=1 for all
    ! dimensions
    !       So we don't have to know start/end indices here
    WRITE(*,'(A5, A17,4I20,ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), data(1,1), &
                            data(SIZE(data,1), SIZE(data,2)),            &
                            SQRT(REAL(SUM(data**2) / SIZE(data)))

  END SUBROUTINE print_2d_variable_int

  !------------------------------------------------------------------
  ! print_3d_variable
  !
  ! Prints statistics for a 3d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_3d_variable_int(name, data)

    CHARACTER(LEN=*) :: name
    REAL(4)         :: data(:,:,:)

    ! Note: Assumed shape array sections always start with index=1 for all dimensions
    !       So we do not have to know start/end indices here
    WRITE(*,'(A5,A17,4I20,ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), data(1,1,1),  &
                            data(SIZE(data,1), SIZE(data,2), SIZE(data,3)), &
                            SQRT(REAL(SUM(data**2) / SIZE(data)))

  END SUBROUTINE print_3d_variable_int

  !------------------------------------------------------------------
  ! print_4d_variable
  !
  ! Prints statistics for a 4d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_4d_variable_int(name, data)

    CHARACTER(LEN=*) :: name
    REAL(4)         :: data(:,:,:,:)

    ! Note: Assumed shape array sections always start with index=1 for all dimensions
    !       So we do not have to know start/end indices here
    WRITE(*,'(A5,A17,4I20,ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), data(1,1,1,1),  &
                            data(SIZE(data,1), SIZE(data,2), SIZE(data,3), SIZE(data,4)), &
                            SQRT(REAL(SUM(data**2) / SIZE(data)))

  END SUBROUTINE print_4d_variable_int


  !------------------------------------------------------------------
  ! print_5d_variable
  !
  ! Prints statistics for a 5d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_5d_variable_int(name, data)

    CHARACTER(LEN=*) :: name
    REAL(4)         :: data(:,:,:,:,:)

    ! Note: Assumed shape array sections always start with index=1 for all dimensions
    !       So we do not have to know start/end indices here
    WRITE(*,'(A5,A17,4I20,ES20.10)') "TEST ", name, MINVAL(data), MAXVAL(data), data(1,1,1,1,1),  &
                            data(SIZE(data,1), SIZE(data,2), SIZE(data,3), SIZE(data,4), SIZE(data,5)), &
                            SQRT(REAL(SUM(data**2) / SIZE(data)))

  END SUBROUTINE print_5d_variable_int


END MODULE thompson_utils
