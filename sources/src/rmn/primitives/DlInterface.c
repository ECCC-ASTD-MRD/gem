#include <DlInterface.h>
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

static void *_DlOpen_(const char *filename, int flag)
{
  return(NULL);
}
static void *(*P_DlOpen)(const char *,int) = _DlOpen_;
void *DlOpen(const char *filename, int flag)
{
  return((*P_DlOpen)(filename,flag));
}

void *_DlSym_(void *handle, const char *symbol)
{
  return(NULL);
}
static void *(*P_DlSym)(void *,const char *) = _DlSym_;
void *DlSym(void *handle, const char *symbol)
{
  return ((*P_DlSym)(handle,symbol));
}

char *_DlError_(void)
{
  return(ERR_NOT_ACTIVE);
}
static char *(*P_DlError)(void) = _DlError_;
char *DlError(void)
{
  return ((*P_DlError)());
}

int _DlClose_(void *handle)
{
  return(-1);
}
static int (*P_DlClose)(void *) = NULL;
int DlClose(void *handle)
{
  return ((*P_DlClose)(handle));
}

void DlRegister(void *open, void *sym, void *error, void *close)
{
  P_DlOpen = (void *(*)(const char *,int)) open;
  P_DlSym = (void *(*)(void *,const char *)) sym;
  P_DlError = (char *(*)(void)) error;
  P_DlClose = (int (*)(void *)) close;
}
