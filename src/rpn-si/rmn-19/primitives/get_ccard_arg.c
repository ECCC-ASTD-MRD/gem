#include <stdio.h>
#include <stdlib.h>
char* get_ccard_arg_m(char* nom, int ord)
{
  char keyname[80], value;
  char *pt = NULL;
  int lng, lng2;
  
  sprintf(keyname,"%s%s%04d%s","%%",nom,ord,"%%");
  /* printf("Debug keyname --->>%s<<---\n",keyname); */
  /* printf("Debug strlen(keyname)=%d\n",strlen(keyname)); */
  lng = get_appl_var(keyname,value,strlen(keyname),0);
  if (lng != 0) {
    pt = malloc(-lng);
    lng2 = get_appl_var(keyname,pt,strlen(keyname),-lng);
    pt[lng2] = '\0';
    }
  return(pt);
}
