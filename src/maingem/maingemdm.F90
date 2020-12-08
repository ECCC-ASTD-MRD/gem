program gem
   call gemdm
end program gem

       subroutine atm_model_getversion2(name_S,version_S,date_S,arch_S,compil_S,user_S,is_official_L)
          implicit none
          character(len=*),intent(out) :: name_S,version_S,date_S,arch_S,user_S,compil_S
          logical,intent(out) :: is_official_L
          name_S = "GEMDM"
          version_S = "5.1.0-rc3"
          date_S = "2020-12-08 14:13 GMT"
          arch_S = ""
          user_S = ""
          compil_S = ""
          is_official_L = .false.
          return
       end subroutine atm_model_getversion2
