
#if !defined(__FTN2C_HELPER__)
#define __FTN2C_HELPER__

#include <rpnmacros.h>

#define FTN2C_FSTR2CSTR(FSTR,CSTR,LFSTR,LCSTR) ftn2c_string_copy((unsigned char *)FSTR,(unsigned char *)CSTR,LFSTR,LCSTR,0)
#define FTN2C_CSTR2FSTR(CSTR,FSTR,LCSTR,LFSTR) ftn2c_string_copy((unsigned char *)CSTR,(unsigned char *)FSTR,LCSTR,LFSTR,' ')
#define FTN2C_FSTR2CSTR_A(FSTR,CSTR,LFSTR,LCSTR,NITEMS) ftn2c_fstra_cstra((unsigned char *)FSTR,(unsigned char **)CSTR,LFSTR,LCSTR,NITEMS,'\0')
#define FTN2C_CSTR2FSTR_A(CSTR,FSTR,LCSTR,LFSTR,NITEMS) ftn2c_cstra_fstra((unsigned char **)CSTR,(unsigned char *)FSTR,LCSTR,LFSTR,NITEMS,' ')

int ftn2c_string_copy(unsigned char *src, unsigned char *dest, int lsrc, int ldest, unsigned char pad);

int ftn2c_cstra_fstra(unsigned char **src, unsigned char *dest, int lsrc, int ldest, int nitems, unsigned char pad);

int ftn2c_fstra_cstra(unsigned char *src, unsigned char **dest, int lsrc, int ldest, int nitems, unsigned char pad);

#endif
