module adz_BC_deficit 

      implicit none
      public
      save

      integer, parameter :: MAXTR3D_ = 250

      real*8 :: KEEP_mass_deficit_8(MAXTR3D_)
      character(len=4) :: tracer_name(MAXTR3D_)

end module adz_BC_deficit 
