
#if !defined(__ARMNLIB__)
#define __ARMNLIB__

#include <rpnmacros.h>

int c_fnom(int *iun,char *filename,char *options,int i);
int c_fstouv(int iun,char *filename,char *options);
void c_fstvoi(int iun,char *options);
int c_fstrwd(int iun);
int c_fstinf(int iun,int *ni, int *nj,int *nk,int datev,char *etiket,
    int ip1,int ip2,int ip3,char *typvar,char *nomvar);
int c_fstsui(int iun,int *ni, int *nj,int *nk);
int c_fstinfx(int inhandle,int iun,int *ni, int *nj,int *nk,int datev,char *etiket,
    int ip1,int ip2,int ip3,char *typvar,char *nomvar);
int c_fstprm(int handle,int *dateo,int *deet,int *npas,int *ni,int *nj,int *nk,
    int *nbits,int *datyp,int *ip1,int *ip2,int *ip3,char *TYPVAR,char *NOMVAR,
    char *ETIKET,char *GRTYP,int *ig1,int *ig2,int *ig3,int *ig4,
    int *swa,int *lng,int *dltf,int *ubc,int *extra1,int *extra2,int *extra3);
int c_fstinl(int iun,int *ni, int *nj,int *nk,int datev,char *etiket,
    int ip1,int ip2,int ip3,char *typvar,char *nomvar,int *liste,int *nliste,int nmax);
int c_fstluk(void *data,int handle,int *ni, int *nj,int *nk);
int c_fst_edit_dir(int handle,int date,int deet,int npas,int ni,int nj,int nk,
    int ip1,int ip2,int ip3,char *typvar,char *nomvar,char *etiket,char *grtyp,
    int ig1,int ig2,int ig3,int ig4,int datyp);
int c_fstecr(void *data,void *work,int nbits,int iun,int dateo,int deet,int npas,
    int ni,int nj,int nk,int ip1,int ip2,int ip3,char *typvar,char *nomvar,
    char *etiket,char *grtyp,int ig1,int ig2,int ig3,int ig4,int datypdtl,int rewrit);
int c_fsteff(int handle);
int c_fstfrm(int iun);
int c_fclos(int iun);
void f77name(cxgaig)(char *grtyp,
    F77_INTEGER *ig1,F77_INTEGER *ig2,F77_INTEGER *ig3,F77_INTEGER *ig4,
    F77_REAL *xg1,F77_REAL *xg2,F77_REAL *xg3,F77_REAL *xg4);
void f77name(cigaxg)(char *grtyp,
    F77_REAL *xg1,F77_REAL *xg2,F77_REAL *xg3,F77_REAL *xg4,
    F77_INTEGER *ig1,F77_INTEGER *ig2,F77_INTEGER *ig3,F77_INTEGER *ig4);
void f77name(convip_plus)(F77_INTEGER *ipnew,F77_REAL *level,F77_INTEGER *fkind,
    F77_INTEGER *fmode,char *strg,F77_INTEGER *flag,F77_INTEGER strglen);
int f77name(newdate)(F77_INTEGER *fdat1,F77_INTEGER *fdat2,
    F77_INTEGER *fdat3,F77_INTEGER *fmode);
int f77name(difdatr)(F77_INTEGER *fdat1,F77_INTEGER *fdat2,F77_REAL8 *fnhours);
int f77name(incdatr)(F77_INTEGER *fdat1,F77_INTEGER *fdat2,F77_REAL8 *fnhours);
int c_gdll(int gdid, float *lat, float *lon);
int c_gdxyfll(int gdid,float *x,float *y,float *lat,float *lon,int n);
int c_gdllfxy(int gdid,float *lat,float *lon,float *x,float *y,int n);
int c_ezdefset(int gdid_dst,int gdid_src);
int c_ezgdef_fmem(int ni,int nj,char *grtypZ,char *grref,
    int ig1,int ig2,int ig3,int ig4,float *xs,float *ys);
int c_ezqkdef(int ni,int nj,char *grtyp,int ig1,int ig2,int ig3,int ig4,int x);
int c_ezuvint(void *newarray,void *newarray2,void *arrayin,void *arrayin2);
int c_ezsint(void *newarray,void *arrayin);
// Add get/setopt and val and ival
int c_ezgetopt(char * option, char * value);
int c_ezsetopt(char * option, char * value);
int c_ezgetval(char * option, float * value);
int c_ezsetval(char * option, float * value);
int c_ezgetival(char * option, int * value);
int c_ezsetival(char * option, int value); // This appears to pass a literal integer
int c_fstopi(char *optname,int lvl,int setget);
// Add grid release function
int c_gdrls(int gdid);

// Add wind conversion routines
int c_gdwdfuv(int gdid, float * spdllout, float * dirllout, float * uugdin,
               float * vvgdin, float * lat, float * lon, int npts);
int c_gduvfwd(int gdid, float * uugdout, float * vvgdout, float * spdllin,
               float * dirllin, float * lat, float * lon, int npts);

// Scattered point interpolation

// (lat,lon) scalar interpolation
int c_gdllsval(int gdid, float * zvals, float * zin, float * lat, float * lon, int n);
// (x,y) scalar interpolation
int c_gdxysval(int gdid, float * zvals, float * zin, float * x, float * y, int n);
// (lat,lon) vector interpolation
int c_gdllvval(int gdid, float * uuvals, float * vvvals, float * uuin, float * vvin,
               float * lat, float * lon, int n);
// (x,y) vector interpolation
int c_gdxyvval(int gdid, float * uuvals, float * vvvals, float * uuin, float * vvin,
               float * x, float * y, int n);

#endif
