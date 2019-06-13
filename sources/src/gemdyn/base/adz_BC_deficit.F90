module adz_BC_deficit

      use, intrinsic :: iso_fortran_env
      implicit none
      public
      save

      integer, parameter :: MAXTR3D_ = 250

      real(kind=REAL64) :: KEEP_mass_deficit_8(MAXTR3D_)
      character(len=4) :: tracer_name(MAXTR3D_)

end module adz_BC_deficit
