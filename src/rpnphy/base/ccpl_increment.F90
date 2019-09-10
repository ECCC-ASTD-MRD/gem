subroutine ccpl_increment(F_int,F_intprev,F_ext,F_extprev,F_intlev,F_extlev,F_n,F_nkint,F_nkext,F_name,F_interp)

  implicit none
!!!#include <arch_specific.hf>

  ! Update a column-coupled internal field based on external forcings using the form:
  !    (internal tendency) = (external tendency)
  !    F_int - F_intprev = F_ext - F_extprev
  ! expressed as,
  !    F_int = F_intprev + (F_ext - F_extprev)

  ! Argument declarations
  integer :: F_n                                !Horizontal dimension
  integer :: F_nkint                            !Number of internal component levels
  integer :: F_nkext                            !Number of external model levels
  real, dimension(F_n,F_nkint) :: F_int         !Internal state of coupled component to be updated (internal component levels)
  real, dimension(F_n,F_nkint) :: F_intprev     !State of coupled component to update from (internal component levels)
  real, dimension(F_n,F_nkext) :: F_ext         !Current state of external driver (external model levels)
  real, dimension(F_n,F_nkext) :: F_extprev     !State of external driver at last update (external model levels)
  real, dimension(F_n,F_nkint) :: F_intlev      !Pressure of internal component levels for this field (Pa)
  real, dimension(F_n,F_nkext) :: F_extlev      !Pressure of external model levels for this field (Pa)
  character(len=*) :: F_name                       !Name of field (used for extrapolation if necessary)
  character(len=*) :: F_interp                     !Order of interpolation ('linear' or 'cubic')

  ! Internal variables
  real, dimension(F_n,F_nkext) :: delta_ext
  real, dimension(F_n,F_nkint) :: delta_int

  ! Obtain driving model differences since last call to component
  delta_ext = F_ext - F_extprev

  ! Vertical interpolation to component levels of forcing fields
  call vte_intvertx3(delta_int,delta_ext,F_extlev,F_intlev,F_n,F_nkext,F_nkint,F_name,F_interp)

  ! Apply interpolated tendency to component state
  F_int = F_intprev + delta_int

  ! End of subprogram
  return
end subroutine ccpl_increment
