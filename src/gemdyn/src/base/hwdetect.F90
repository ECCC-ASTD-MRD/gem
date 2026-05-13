module hwdetect
  use, intrinsic :: iso_c_binding
  interface
    
   integer(c_int) function cpu_per_numa(f_world_communicator) bind(C, name='cpu_per_numa')
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: f_world_communicator
   end function cpu_per_numa

  end interface
end module hwdetect
