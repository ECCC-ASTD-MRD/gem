#if ! defined(RMNLIB_PLUGINS)

#define RMNLIB_PLUGINS

#if defined(ThisIsNeverDefined)
  Fortran anc C callable package to create/use runtime plugins using dynamic libraries

  step 0: 
          Fortran :
            USE ISO_C_BINDING
            #define __FORTRAN_CODE__
          Fortran and C:
            #include <plugins.h>

  step 1: Load a plugin library
          Fortran :
            type(C_PTR) :: handle
            handle = load_plugin("my_library_name.so")
          C :
            void *handle;
            handle = load_plugin("my_library_name.so");

  step 2: Get number of entry points in plugin
          Fortran:
            integer(C_INT) nsym
            nsym = plugin_n_functions(handle)
          C:
            int nsym = plugin_n_functions(handle);

  step 3: Get the name of entry point n
          Fortran:
            type(C_PTR) :: string
            integer(C_INT) :: n
            string = plugin_function_name(handle,n)
          C:
            int n;
            char *string = plugin_function_name(handle,n);

  step 4: Get address of entry point by name and call it
          Fortran:
            type(C_FUNPTR) :: faddress
            character(C_CHAR), dimension(*) :: name
            procedure(xxx), pointer :: fptr
            faddress = plugin_function(handle,name)
            call c_f_procpointer(faddress,fptr)
            .. = fptr(...arguments...)
          C:
            void *faddress;
            char *name;
            faddress = plugin_function(handle,name);
            .. = (*faddress)(...arguments...);

  step n: unload plugin
          Fortran:
            integer(C_INT) :: status
            status = unload_plugin(handle)
          C:
            int status = unload_plugin(handle);

  other : set verbosity
          Fortran:
            integer(C_INT) :: verbose
            call set_plugin_diag(verbose)
          C:
            int verbose;
            set_plugin_diag(verbose);

  NOTES:
    the presence of entry point "entry_list" is mandatory, it must 
    a NULL pointer terminated list of pointers to NULL terminated strings
    providing the names of the advertised entry points in the plugin

    function "get_symbol_number" is optional and may be used to set the
    proper values in "entry_list" (see Fortran)

----------------------- Example of C plugin -----------------------
cc -shared -fpic -o libxxx.so xxx.c

#include <stdio.h>
#include <string.h>

char *entry_list[4] = { "name1","name2","name3",NULL} ;

int name1(int arg){
printf("name1: %d\n",arg);
return(arg);
}

int name2(int arg){
printf("name2: %d\n",arg);
return(arg);
}

int name3(int arg){
printf("name3: %d\n",arg);
return(arg);
}

int get_symbol_number(){  // like fortran, function to get number of symbols, optional
  return(3);
}

----------------------- Example of Fortran plugin -----------------------
                ( needs a little more helper code than C )
fcompiler -shared -fpic -o libxxx.so xxx.F90

! beginning of routines in plugin
integer function  fn1(arg) BIND(C,name='name1f')
integer, intent(IN) :: arg
print *,'fortran name1 =',arg
fn1 = arg
return
end

integer function  fn2(arg) BIND(C,name='name2f')
integer, intent(IN) :: arg
print *,'fortran name2 =',arg
fn2 = arg
return
end
! end of routines in plugin
!
! what follows is boiler plate code
! to be adjusted by user : MAX_NAMES, MAX_NAME_LENGTH, calls to insert_in_name_table in subroutine symbols
module interop
  use ISO_C_BINDING
  implicit none
! start of user adjusted code
#define MAX_NAMES 2
#define MAX_NAME_LENGTH 8
! end of user adjusted code
  type(C_PTR), dimension(MAX_NAMES+1), save, target, BIND(C,name='entry_list') :: name_table
  character(len=1), dimension(MAX_NAME_LENGTH+1,MAX_NAMES), save, target :: names
  integer, save :: nargs
  contains
  subroutine insert_in_name_table(name)  ! add name to name table and neme pointers
    use ISO_C_BINDING
    implicit none
    character(len=*) :: name
    integer :: l
    l = len(trim(name)) + 1
    nargs = nargs + 1
    names(1:l,nargs) = transfer(trim(name)//achar(0) , names, l)
    name_table(nargs) = C_LOC(names(1,nargs))
    return
  end subroutine insert_in_name_table
  function symbols() bind(C,name='get_symbol_number') result(number)
    use ISO_C_BINDING
    implicit none
    integer(C_INT) :: number
    nargs = 0
! start of user adjusted code
    call insert_in_name_table('name1f')
    call insert_in_name_table('name2f')
! end of user adjusted code
    number = nargs   ! return number of arguments
    return
  end function symbols
end module interop

#endif

#if defined(IN_FORTRAN_CODE) || defined(__FORTRAN_CODE__)

interface

  function load_plugin(plugin_name) result(handle) BIND(C,name='load_plugin')
    import :: C_CHAR, C_INT, C_PTR
    character(C_CHAR), dimension(*), intent(IN) :: plugin_name
    type(C_PTR) :: handle
  end function load_plugin

  function unload_plugin(handle) result(status) BIND(C,name='unload_plugin')
    import :: C_INT, C_PTR
    integer(C_INT) :: status
    type(C_PTR), value :: handle
  end function unload_plugin

  function plugin_function(handle,fname) result(faddress) BIND(C,name='plugin_function')
    import :: C_PTR, C_FUNPTR, C_CHAR
    type(C_PTR), value :: handle
    character(C_CHAR), dimension(*), intent(IN) :: fname
    type(C_FUNPTR) :: faddress
  end function plugin_function

  function plugin_function_name(handle,ordinal) result(string) BIND(C,name='plugin_function_name')
    import :: C_PTR, C_INT
    type(C_PTR), value :: handle
    integer(C_INT), value :: ordinal
    type(C_PTR) :: string
  end function plugin_function_name

  function plugin_n_functions(handle) result(number) BIND(C,name='plugin_n_functions')
    import :: C_PTR, C_INT
    type(C_PTR), value :: handle
    integer(C_INT) :: number
  end function plugin_n_functions

  subroutine set_plugin_diag(level) BIND(C,name='set_plugin_diag')
    import C_INT
    integer(C_INT), value :: level
  end subroutine set_plugin_diag

  function c_strlen(string) result(length) BIND(C,name='strlen')
    import C_SIZE_T, C_PTR
    type(C_PTR), value :: string
    integer(C_SIZE_T) :: length
  end function c_strlen

end interface

#else

typedef const char * charptr;       // pointer to char string
typedef int (*fnptr)();             // pointer to function
typedef void * (*fnpptr)();             // pointer to function

#if ! defined(NO_PROTOTYPES)
void *load_plugin(const char *lib);
int unload_plugin(const void *p);
const void *plugin_function(const void *p, const char *name);
charptr *plugin_function_names(const void *p);
charptr plugin_function_name(const void *p, int ordinal);
int plugin_n_functions(const void *p);
void set_plugin_diag(int diag);
#endif

#endif

#endif