#include <stdio.h>
#include <stdlib.h>
#include <rmnlib.h>

#define MAX_ENTRIES 1024

typedef struct {
    char * name;
    char * value;
    int ncn;         /* number of caracter in name */
    int ncv;         /* number of caracter in value */
} appl_var_entry;

static appl_var_entry appl_var_table[MAX_ENTRIES];
static int In_Used;

void init_appl_var_table()
{
  int i;
  
  for (i=0; i < MAX_ENTRIES; i++) {
    appl_var_table[i].name = (char *) NULL;
    appl_var_table[i].value = (char *) NULL;
    appl_var_table[i].ncn = 0;
    appl_var_table[i].ncv = 0;
    }
  In_Used = 0;
}

void free_appl_var_table()
{
  int i;
  
  for (i=0; i < In_Used; i++) {
    if (appl_var_table[i].name) free(appl_var_table[i].name);
    if (appl_var_table[i].value) free(appl_var_table[i].value);
    appl_var_table[i].ncn = 0;
    appl_var_table[i].ncv = 0;
    }
  In_Used = 0;
}

void set_appl_var(char* name, char* value, int ln, int lv)
{
  int i, ind=-1;
  while ((ln > 0) && (name[ln-1] == ' ')) ln--;
  while ((lv > 0) && (value[lv-1] == ' ')) lv--;  
  if (ln != 0) {
    for (i = 0; i < In_Used; i++)
      if (strncmp(name,appl_var_table[i].name,ln) == 0) {
        ind = i;
        break;
        }
    }
  else {  
    for (i = 0; i < In_Used; i++)
      if (strcmp(name,appl_var_table[i].name) == 0) {
        ind = i;
        break;
       }
     }
      
  if (ind != -1)
    i = ind;
  else
    i = In_Used++;
    
  if (appl_var_table[i].name) free(appl_var_table[i].name);
  if (appl_var_table[i].value) free(appl_var_table[i].value);
  appl_var_table[i].name = malloc(ln+1);
  appl_var_table[i].value = malloc(lv+1);
  appl_var_table[i].ncn = ln;
  appl_var_table[i].ncv = lv;
  strncpy(appl_var_table[i].name,name,ln);
  strncpy(appl_var_table[i].value,value,lv);
}

int get_appl_var(char* varname,char *value, int ln, int lng)
{
  int i, ind=-1;

  while ((ln > 0) && (varname[ln-1] == ' ')) ln--;
  
  if (ln != 0) {
    for (i = 0; i < In_Used;  i++) {
/*      printf("i=%d appl_var_table[i].name=%s varname=%s\n",i,appl_var_table[i].name,varname); */
      if (strncasecmp(varname,appl_var_table[i].name,ln) == 0) {
        ind = i;
        break;
        }
      }
    }
  else {
    for (i = 0; i < In_Used; i++)
      if (strcasecmp(varname,appl_var_table[i].name) == 0) {
        ind = i;
        break;
        }
    } 
  if (ind == -1) return(0);
  strncpy(value,appl_var_table[ind].value,lng);
  return((lng >= appl_var_table[ind].ncv) ? appl_var_table[ind].ncv : -(appl_var_table[ind].ncv));
}

ftnword f77name(c_get_appl_var)(char* name, char* value, F2Cl lln, F2Cl llv)
{

  int i, lng, ln=lln, lv=llv;
  
  lng = get_appl_var(name,value,ln,lv);
  i = lng;
  while (i <= lv)
    value[i++] = ' ';          /* remove null, blank pad */
  return((ftnword) lng);
}

void f77name(c_init_appl_var_table)()
{
  init_appl_var_table();
}

void f77name(c_free_appl_var_table)()
{
  free_appl_var_table();
}

void f77name(c_set_appl_var)(char* name, char* value, F2Cl lln, F2Cl llv)
{
  int ln=lln, lv=llv;
  set_appl_var(name,value,ln,lv);
}  

