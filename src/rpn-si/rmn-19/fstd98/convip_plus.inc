      interface

         !> Translate ip kind into a 2 character string code
         FUNCTION kind_to_string(code) RESULT(string)
            integer, intent(IN) :: code
            character(len=2) :: string
         end FUNCTION kind_to_string

         integer function value_to_string(val,string,maxlen)  
            integer, intent(IN) :: maxlen
            character (len=*), intent(OUT) :: string
            real *4, intent(IN) :: val
         end function value_to_string

      end interface
