#ifndef DlInterface_DEFINED
#define DlInterface_DEFINED
#include <stdlib.h>
#include <dlfcn.h>
/*
  rmnlib interface to the dynamic loading functions
  if this module is compiled without -DLIVE
  it provides stubs that return a failure code when called
  the purpose is to have a default version in the library that does not
  necessitate -ldl at link time for applications
  if this module is compiled with -DLIVE
  it becomes a direct interface do dlopen/dlsym/dlerror/dlclose
*/
#define ERR_NOT_ACTIVE "ERROR: this is the dummy dynamic loader\n"

void *DlOpen(const char *filename, int flag);
void *DlSym(void *handle, const char *symbol);
char *DlError(void);
int DlClose(void *handle);
void DlRegister(void *open, void *sym, void *error, void *close);
#endif
