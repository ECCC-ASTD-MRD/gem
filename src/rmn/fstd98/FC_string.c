#include <rpnmacros.h>
#include <stdio.h>
#include <malloc.h>

char **allocate_string_array(int ns)
{
  int i;
  char **string_array;

  string_array = malloc((ns+1) * sizeof(char *));
  for (i=0; i <= ns; i++)
    string_array[i] = (char *) NULL;

  return(string_array);
}

void free_string_array(char **string_array)
{
  int i=0;

  while (string_array[i]) {
/*    printf("Debug free_string_array i=%d\n",i);*/
    free(string_array[i]);
    i++;
  }
  free(string_array);
}

void free_cstring(char *cstring)
{
  free(cstring);
}

void cstring_to_fstring(char *cstring, char *fstring, int nc)
{
  int i,j;

  i = 0;
  while ((*cstring) && (i < nc)) {
    *fstring = *cstring;
    fstring++;
    cstring++;
    i++;
  }

  for (j=i; j < nc; j++) {
    *fstring = ' ';
    fstring++;
  }
}

char *fstring_to_cstring(char *fstring, int nc, int rmblanks)
{
  int i;
  char *cstring, *ctemp;

  cstring = malloc(nc+1);
  ctemp = cstring;
  for (i=0; i < nc; i++) {
    *ctemp = *fstring;
    ctemp++;
    fstring++;
  }

  *ctemp = '\0';
  ctemp--;
  if (rmblanks) {
    i = nc;
    while ((i > 0) && (*ctemp == ' ')) {
      *ctemp = '\0';
      ctemp--;
    }
  }

  return(cstring);
}

char **fill_string_array(char **string_array, char *farray, int nc, int ns, int rmblanks)
{
  int i;

  for (i=0; i<ns; i++) {
    string_array[i] = fstring_to_cstring(farray,nc,rmblanks);
    farray += nc;
  }

  return(string_array);
}

void f77name(fs_to_cs)(char *fstring, int *rmblanks, int *ns, F2Cl fnc)
{
  char *cstring, *cmpstring;
  char **string_array;
  int i, nc=fnc;

  cmpstring = malloc(13);

  if (*ns == 1) {
    cstring = fstring_to_cstring(fstring,nc,*rmblanks);
    printf("Debug fs_to_cs cstring-->%s<--\n",cstring);
    strcpy(cmpstring,"Label01");
    printf("Debug fs_to_cs cmpstring-->%s<--\n",cmpstring);
    printf("Debug fs_to_cs strncmp sans blancs=%d\n",strncmp(cstring,cmpstring,13));
    strcpy(cmpstring,"Label01     ");
    printf("Debug fs_to_cs cmpstring-->%s<--\n",cmpstring);
    printf("Debug fs_to_cs strncmp avec blancs=%d\n",strncmp(cstring,cmpstring,13));
  }
  else {
    string_array = fill_string_array(allocate_string_array(*ns),fstring,nc,*ns,*rmblanks);
/*    for (i=0; i < *ns; i++)
      printf("Debug string_array[%d]-->%s<--\n",i,string_array[i]);*/
  }
}

void f77name(cs_to_fs)(char *fstring, int *nc)
{
  char *cstring;

  cstring = malloc(12);
  strcpy(cstring,"Test001");

  cstring_to_fstring(cstring,fstring,*nc);

}
