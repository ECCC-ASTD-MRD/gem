module hwdetect
  use, intrinsic :: iso_c_binding
  interface
    
   integer(c_int) function cpu_per_numa() bind(C, name='cpu_per_numa')
      use, intrinsic :: iso_c_binding
   end function cpu_per_numa

  end interface
end module hwdetect
