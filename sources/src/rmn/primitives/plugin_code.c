#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dlfcn.h>

#define NO_PROTOTYPES
#include <plugins.h>

typedef struct {
  void *handle;       // blind pointer to plugin handle (from dlopen)
  char *name;         // pointer to plugin name
  charptr *symbol;    // pointer to list symbol table names (from plugin)
  fnptr *addr;        // list of functions addresses (function names in symbol) (from dlsym)
  int nentries;       // number of functions advertised in plugin (size of symbol table/address list)
  int ordinal;        // index in plugin table
} plugin;

#define MAX_PLUGINS 256
static plugin plugin_table[MAX_PLUGINS];
static int last_plugin=0;
static int verbose = 0;

int unload_plugin(plugin *p){
  if(p == NULL) return(-1) ; // bad handle
  if( (p - plugin_table) >= last_plugin) return(-1);  // out of table

  if(p->addr) free(p->addr);           // free address table if allocated
  p->addr = NULL;
  p->symbol = NULL;
  dlclose(p->handle);
  p->handle = NULL;
  p->ordinal = -1;
  p->nentries = 0;
  if(verbose) printf("INFO: plugin %s has been closed (slot %ld)\n",p->name,p-plugin_table);
  if(p->name) free(p->name);
  p->name = NULL;
  return (0);   // entry has been freed
}

// set diagnostic level
// 0 = silent
// 1 = verbose
void set_plugin_diag(int diag){
  verbose = diag;
}
//----------------------------------------------------------------------------------------
//
//  load a plugin (shared object)
//
//  lib     : character string, name of shared object
//  verbose : int, if non zero some diagnostics are printed
//
//  function return :  handle (pointer to a plugin structure)
//
//  the rules for finding "lib" are the rules followed by dlopen (see: man dlopen)
//
plugin *load_plugin(const char *lib, int diag){
  const char **temp;
  int nsym;
  int i;
  plugin *p;
  fnptr func_nb;
  int slot ;

  verbose = diag;
  for (slot=0 ; slot < MAX_PLUGINS ; slot++){
    if(plugin_table[slot].name == NULL) break;  // free slot
  }
  if(slot >= MAX_PLUGINS) {
    fprintf(stderr,"ERROR: plugin table is full\n");
    return(NULL);   // table is full
  }

  p = &plugin_table[slot];
  p->ordinal = slot;

  p->name = (char *)malloc(strlen(lib)+1);
  strncpy(p->name, lib, strlen(lib)+1);
  if(verbose) fprintf(stderr,"\nINFO: attempting to load plugin %s (slot %d)\n",p->name,slot);
  p->handle = dlopen(p->name,RTLD_LAZY);
  if( p->handle == NULL ) {
    fprintf(stderr,"ERROR: load plugin failed for %s\n",p->name);
    free(p->name);
    p->name = NULL;
    return(NULL);
  }

  func_nb = dlsym(p->handle,"get_symbol_number");         // provide get_symbol_number (may initialize entry_list) (optional)
  temp = (const char **)dlsym(p->handle,"entry_list");    // provide **entry_list  (mandatory)
  nsym = 0;
  if(func_nb) nsym = (*func_nb)();                        // optional function (mainly for Fortran usage)

  if( (temp == NULL) && (func_nb == NULL)) {
    if(verbose) fprintf(stderr,"WARNING: no function table found in plugin %s\n",p->name);
    dlclose(p->handle);     // nothing useful, close
    free(p->name);
    p->name = NULL;
    return(NULL);
  }

  if(nsym==0){
    nsym = 0 ; while(temp[nsym]) nsym++;
  }

  if (nsym == 0) {
    if(verbose) fprintf(stderr,"WARNING: no functions advertised in plugin %s\n",p->name);
    dlclose(p->handle);     // nothing useful, close
    free(p->name);
    p->name = NULL;
    return(NULL);     // no useful entry points found
  }

  p->symbol = temp;                                       // symbol table is at address of symbol entry_list
  p->nentries = nsym;
  if(verbose) fprintf(stderr,"INFO: %d functions advertised in plugin %s\n",p->nentries,p->name);
  p->addr = (void *)malloc(nsym*sizeof(void *));   // allocate address table
  if(slot >= last_plugin) last_plugin = slot + 1;

  for(i=0 ; i<nsym ; i++){
    p->addr[i] =  dlsym(p->handle,p->symbol[i]);   // fill address table
    if(verbose) fprintf(stderr,"INFO:   %p %s\n",p->addr[i],p->symbol[i]);
  }
  return(p);
}
//----------------------------------------------------------------------------------------
//
//  find number of functions in a given plugin
//
//  p  : handle obtained from load_plugin
//
//  function return :  number of symbols advertised in the shared object
//
int plugin_n_functions(const plugin *p){     // how many functions are defined in this plugin
  if(p == NULL) return(0) ; // bad handle
  if( (p - plugin_table) >= last_plugin) return(0);  // out of table

  return(p->nentries);           // address of list of names
}
//----------------------------------------------------------------------------------------
//
// get name of entry number ordinal from plugin entry name table
//
//  p       : handle obtained from load_plugin
//  ordinal : int, ordinal in entry list of desired name (first entry has ordinal 1)
//
//  function return : pointer to the name of requested entry (char*)
//
charptr plugin_function_name(const plugin *p, int ordinal){
  if(p == NULL) return(NULL) ; // bad handle
  if( (p - plugin_table) >= last_plugin) return(NULL);  // out of plugin table

  if(ordinal<1 || ordinal>p->nentries) return(NULL);  // out of range for this plugin

  return((const charptr)p->symbol[ordinal-1]);           // address of list of names
}
//----------------------------------------------------------------------------------------
//
//  get list of advertised entry names in a given plugin
//
//  p : handle obtained from load_plugin
//
//  function return : pointer to the list of names (char**)
//
charptr *plugin_function_names(const plugin *p){
  if(p == NULL) return(NULL) ; // bad handle
  if( (p - plugin_table) >= last_plugin) return(NULL);  // out of table

  return((charptr *)p->symbol);           // address of list of names
}
//----------------------------------------------------------------------------------------
//
//  get address of plugin entry name
//  if pointer to plugin is NULL, scan all known plugins
//
//  p    : handle obtained from load_plugin
//  name : char *, null terminated string, name of entry 
//
//  function return : address of requested entry
//
void *plugin_function(const plugin *p, const char *name){
  int i, j;
  void *faddr;

  if(p == NULL) {           // scan all plugins for name
    for(j=0 ; j<last_plugin ; j++){
      p = &plugin_table[j];
      for(i=0 ; i<p->nentries ; i++){
        if(strcmp(name, p->symbol[i]) == 0) return(p->addr[i]) ;
      }
    }
    return(NULL);   // nothing found
  }

  if( (p - plugin_table) >= last_plugin) return(NULL);  // out of table

  for(i=0 ; i<p->nentries ; i++){
    if(strcmp(name, p->symbol[i]) == 0) return(p->addr[i]) ;
  }
  if( (faddr=dlsym(p->handle,name)) ) return(faddr); // try for unadvertised name in plugin
  return(NULL);                         // name not found
}
