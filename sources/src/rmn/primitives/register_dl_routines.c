#include <dlfcn.h>

void DlRegister(void *open, void *sym, void *error, void *close);

#pragma weak register_dl_routines_ = register_dl_routines
#pragma weak register_dl_routines__ = register_dl_routines
void register_dl_routines_();
void register_dl_routines__();
void register_dl_routines()
{
 DlRegister(dlopen,dlsym,dlerror,dlclose);
}
