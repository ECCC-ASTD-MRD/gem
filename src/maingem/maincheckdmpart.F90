program maincheckdmpart
   call checkdmpart
end program maincheckdmpart

subroutine atm_model_getversion2(name_S,version_S,date_S,arch_S,compil_S,user_S,is_official_L)
   implicit none
   character(len=*),intent(out) :: name_S,version_S,date_S,arch_S,user_S,compil_S
   logical,intent(out) :: is_official_L
   name_S = "checkdmpart"
   version_S = "5.2.0-a14"
   date_S = "2022-02-16 16:41 GMT"
   arch_S = ""
   user_S = ""
   compil_S = ""
   is_official_L = .false.
   return
end subroutine atm_model_getversion2
