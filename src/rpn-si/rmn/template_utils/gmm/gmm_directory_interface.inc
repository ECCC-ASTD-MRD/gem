       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
