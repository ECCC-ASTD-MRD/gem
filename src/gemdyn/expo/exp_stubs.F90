! Delete this file to run EXPO

subroutine exp_dynstep()
return
end

integer function exp_nml(F_namelistf_S)
character(len=*) F_namelistf_S
exp_nml = 0
return
end

subroutine exp_set_vt()
stop 'exp_set_vt : not yet implemented (stub)'
return
end

subroutine exp_hybrid ( F_hybuser, Nk )
integer Nk
real, dimension(Nk) :: F_hybuser        !user-specified hybrid coordinate values
stop 'exp_hybrid : not yet implemented (stub)'
return
end

module exp_geom
   public :: exp_geometry
   contains
   subroutine exp_geometry()
   return
   end subroutine exp_geometry
end module exp_geom
