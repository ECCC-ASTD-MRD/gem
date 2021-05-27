integer function f_iargc()
f_iargc=command_argument_count()
return
end
subroutine f_getarg(pos,value)
integer pos,len,stat
character (len=*) :: value
call get_command_argument(pos,value,len,stat)
return
end

