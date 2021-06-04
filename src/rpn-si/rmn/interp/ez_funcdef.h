#include "ezscint.h"

#ifndef _ezfuncdef
#include "gd_key2rowcol.hc"

static _Grille **Grille  = NULL;
/*static wordint  **gr_list = NULL;*/
static _Grille **gr_list = NULL;
static wordint nGrilles = 0;
static wordint nGrillesMax = CHUNK*CHUNK;
static wordint cur_log_chunk = 7;
/*static wordint cur_log_chunk = 3;*/

static __thread wordint nsets    = 0;
static __thread wordint iset     = -1;
static __thread wordint iset_gdin = -1;
static __thread wordint iset_gdout = -1;
static __thread _gridset *gridset = NULL;
static  __thread _groptions groptions = { OUI, CUBIQUE,  MAXIMUM, NON, -1, SYM, SCALAIRE, NON, NON, OUI, 16, 0, DISTANCE, NEAREST, 0.5, 3.0, 0.0  };

static wordint log_chunks[]= {0, 1,  2, 3,    4,    5,    6,      7,     8,      9,      10,     11,        12};
static wordint primes[]    = {0, 0,  3, 7,   13,   31,   61,    127,   251,    509,    1021,   2039,      4093};
static wordint chunks[]    = {0, 0,  4, 8,   16,   32,   64,    128,   256,    512,    1024,   2048,      4096};
static wordint primes_sq[] = {0, 0,  3, 61, 251, 1021, 4093,  16381, 65521, 262139, 1048573, 4194301, 16777213};
static wordint chunks_sq[] = {0, 0, 16, 64, 256, 1024, 4096,  16384, 65536, 262144, 1048576, 4194304, 16777216};


/*****************************************************************************/
void EliminerGrille(wordint gridid);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
void f77name(ez_avg)(float *zout, float *x, float *y, int *ni_src, int *nj_src,
            float *zin, int *ni_dst, int *nj_dst, int *extension);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint LireEnrPositionnels(_Grille *gr, wordint iunit, wordint ip1, wordint ip2, wordint ip3, wordint ip4, wordint read);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
void c_llfgr(ftnfloat *lat, ftnfloat *lon, ftnfloat *x, ftnfloat *y, wordint npts,
        ftnfloat latOrigine, ftnfloat lonOrigine, ftnfloat deltaLat, ftnfloat deltaLon);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
unsigned int ez_calc_crc(int *p, int *flen,  float *ax, float *ay, int ni, int nj);
wordint ez_calclatlon(wordint gdid);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
void ez_calcntncof(wordint gdid);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint ez_calcxpncof(wordint gdid);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint ez_calcxy(wordint gdin, wordint gdout);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint c_ez_check_xpndable(wordint *extension, wordint ni, wordint nj, char grtyp, wordint ig1, wordint ig2, wordint ig3, wordint ig4);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint ez_corrval(ftnfloat *zout, ftnfloat *zin,  wordint gdin, wordint gdout);
wordint ez_corrvec(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, wordint gdin, wordint gdout);
wordint ez_corrval_ausud(ftnfloat *zout, ftnfloat *zin,  wordint gdin, wordint gdout);
wordint ez_corrval_aunord(ftnfloat *zout, ftnfloat *zin,  wordint gdin, wordint gdout);
wordint ez_corrvec_aunord(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin,  wordint gdin, wordint gdout);
wordint ez_corrvec_ausud(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin,  wordint gdin, wordint gdout);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint ez_defzones(wordint gdin, wordint gdout);/*, wordint gdin, ftnfloat *x, ftnfloat *y, wordint npts, _zone *zones);*/
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint ez_defzone_dehors(wordint gdin, ftnfloat *x, ftnfloat *y, wordint npts, _zone *zone);
wordint ez_defzone_polenord(wordint gdin, ftnfloat *x, ftnfloat *y, wordint npts, _zone *zone);
wordint ez_defzone_polesud(wordint gdin, ftnfloat *x, ftnfloat *y, wordint npts, _zone *zone);
wordint ez_defzone_nord(wordint gdin, ftnfloat *x, ftnfloat *y, wordint npts, _zone *zone);
wordint ez_defzone_sud(wordint gdin, ftnfloat *x, ftnfloat *y, wordint npts, _zone *zone);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint ez_interp(ftnfloat *zout, ftnfloat *zin, wordint gdin, wordint gdout);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint ez_poleovrw(ftnfloat *zout, wordint gdid);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
void ez_xpncof(wordint *i1,wordint *i2,wordint *j1,wordint *j2,wordint *couverture,
            wordint ni,wordint nj,char grtyp, char grref,
            wordint ig1,wordint ig2,wordint ig3,wordint ig4,wordint sym,
            ftnfloat *ax, ftnfloat *ay);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
void ez_xpnsrcgd(wordint gdid, ftnfloat *zout, ftnfloat *zin);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint c_ezfreegridset(wordint gdid, wordint index);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezdefset)(wordint *gdout, wordint *gdin);
wordint c_ezdefset(wordint gdout, wordint gdin);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezgdef)(wordint *ni, wordint *nj, char *grtyp, char *grref,
                    wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4,
                    ftnfloat *ax, ftnfloat *ay, F2Cl lengrtyp, F2Cl lengrref);
wordint c_ezgdef(wordint ni, wordint nj, char *grtyp, char *grref,
             wordint ig1, wordint ig2, wordint ig3, wordint ig4, ftnfloat *ax, ftnfloat *ay);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezgdef_ffile)(wordint *ni, wordint *nj, char *grtyp,
            wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4,
            wordint *iunit, F2Cl lengrtyp);
wordint c_ezgdef_ffile(wordint ni, wordint nj, char *grtyp,
           wordint ig1, wordint ig2, wordint ig3, wordint ig4, wordint iunit);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezgdef_fll)(wordint *ni, wordint *nj, ftnfloat *lat, ftnfloat *lon);
wordint c_ezgdef_fll(wordint ni, wordint nj,ftnfloat *lat, ftnfloat *lon);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezgdef_fmem)(wordint *ni, wordint *nj, char *grtyp, char *grref,
                    wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4,
                    ftnfloat *ax, ftnfloat *ay, F2Cl lengrtyp, F2Cl lengrref);
wordint c_ezgdef_fmem(wordint ni, wordint nj, char *grtyp, char *grref,
             wordint ig1, wordint ig2, wordint ig3, wordint ig4, ftnfloat *ax, ftnfloat *ay);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

wordint f77name(ezgdef_supergrid)(wordint *ni, wordint *nj, char *grtyp, char *grref, wordint *vercode, wordint *nsubgrids, wordint *subgrid, F2Cl lengrtyp, F2Cl lengrref);

wordint c_ezgdef_supergrid(wordint ni, wordint nj, char *grtyp, char *grref, wordint vercode, wordint nsubgrids, wordint *subgrid);
wordint c_ezgdef_yymask(_Grille *gr);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezgenpole)(ftnfloat *vpolnor, ftnfloat *vpolsud, ftnfloat *fld,
                           wordint *ni, wordint *nj, wordint *vecteur,
                           char *grtyp, wordint *hem, F2Cl lengrtyp);
wordint c_ezgenpole(ftnfloat *vpolnor, ftnfloat *vpolsud, ftnfloat *fld,
                           wordint ni, wordint nj, wordint vecteur,
                           char *grtyp, wordint hem);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezgetopt)(char *option, char *value, F2Cl lenoption, F2Cl lenvalue);
wordint c_ezgetopt(char *option, char *value);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezgetival)(char *option, ftnword *value, F2Cl lenoption);
wordint c_ezgetival(char *option, ftnword *value);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezgetval)(char *option, ftnfloat *value, F2Cl lenoption);
wordint c_ezgetval(char *option, ftnfloat *value);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezgfstp)(wordint *gdid,
         char *nomvarx, char *typvarx, char *etiketx,
         char *nomvary, char *typvary, char *etikety,
         wordint *ip1, wordint *ip2, wordint *ip3, wordint *dateo,
                     wordint *deet, wordint *npas, wordint *nbits,
         F2Cl lennomvarx, F2Cl lentypvarx, F2Cl lenetiketx,
         F2Cl lennomvary, F2Cl lentypvary, F2Cl lenetikety);
wordint c_ezgfstp(wordint gdid, char *nomvarx, char *typvarx, char *etiketx,
              char *nomvary, char *typvary, char *etikety,
              wordint *ip1, wordint *ip2, wordint *ip3, wordint *dateo, wordint *deet, wordint *npas, wordint *nbits);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezgprm)(wordint *gdid, char *grtyp, wordint *ni, wordint *nj,
             wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4, F2Cl lengrtyp);
wordint   c_ezgprm(wordint gdid, char *grtyp, wordint *ni, wordint *nj, wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezgxprm)(wordint *gdid, wordint *ni, wordint *nj, char *grtyp,
                     wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4,
                     char *grref, wordint *ig1ref, wordint *ig2ref,
                     wordint *ig3ref, wordint *ig4ref,
                     F2Cl lengrtyp, F2Cl lengrref);
wordint c_ezgxprm(wordint gdid, wordint *ni, wordint *nj,
              char *grtyp, wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4,
              char *grref, wordint *ig1ref, wordint *ig2ref, wordint *ig3ref, wordint *ig4ref);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdll)(wordint *gdid, ftnfloat *lat, ftnfloat *lon);
wordint c_gdll(wordint gdid, ftnfloat *lat, ftnfloat *lon);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezqkdef)(wordint *ni, wordint *nj, char *grtyp,
                    wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4, wordint *iunit, F2Cl lengrtyp);
wordint c_ezqkdef(wordint ni, wordint nj, char *grtyp,
             wordint ig1, wordint ig2, wordint ig3, wordint ig4, wordint iunit);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezquickdef)(wordint *ni, wordint *nj, char *grtyp,
          wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4, wordint *iunit, F2Cl lengrtyp);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint c_ezquickdef(wordint ni, wordint nj, char *grtyp,
         wordint ig1, wordint ig2, wordint ig3, wordint ig4, wordint iunit);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdrls)(wordint *gdin);
wordint c_gdrls(wordint gdin);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezsetopt)(char *option, char *value, F2Cl lenoption, F2Cl lenvalue);
wordint c_ezsetopt(char *option, char *value);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezsetival)(char *option, wordint *ivalue, F2Cl lenoption);
wordint c_ezsetival(char *option, wordint ivalue);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezsetval)(char *option, ftnfloat *fvalue, F2Cl lenoption);
wordint c_ezsetval(char *option, ftnfloat fvalue);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezsint)(ftnfloat *zout, ftnfloat *zin);
wordint c_ezsint(ftnfloat *zout, ftnfloat *zin);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint c_ez_find_gdin(int gdin, int gdout);
wordint find_gdin_in_gset(wordint gdin, wordint gdout);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

wordint f77name(ezuvint)(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin);
wordint c_ezuvint(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezwdint)(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin);
wordint c_ezwdint(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint ftnstrclean(char *str, wordint lenstr);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdgaxes)(wordint *gdid, ftnfloat *ax, ftnfloat *ay);
wordint c_gdgaxes(wordint gdid, ftnfloat *ax, ftnfloat *ay);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdgxpndaxes)(wordint *gdid, ftnfloat *ax, ftnfloat *ay);
wordint c_gdgxpndaxes(wordint gdid, ftnfloat *ax, ftnfloat *ay);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdllfxy)(wordint *gdid, ftnfloat *lat, ftnfloat *lon, ftnfloat *x, ftnfloat *y, wordint *n);
wordint c_gdllfxy(wordint gdid, ftnfloat *lat, ftnfloat *lon, ftnfloat *x, ftnfloat *y, wordint n);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdllfxyz)(wordint *gdid, ftnfloat *lat, ftnfloat *lon, ftnfloat *x, ftnfloat *y, wordint *n);
wordint c_gdllfxyz(wordint gdid, ftnfloat *lat, ftnfloat *lon, ftnfloat *x, ftnfloat *y, wordint n);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdllsval)(wordint *gdid, ftnfloat *zout, ftnfloat *zin, ftnfloat *lat, ftnfloat *lon, wordint *n);
wordint c_gdllsval(wordint gdid, ftnfloat *zout, ftnfloat *zin, ftnfloat *lat, ftnfloat *lon, wordint n);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdllvval)(wordint *gdid, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin,
                      ftnfloat *lat, ftnfloat *lon, wordint *n);
wordint c_gdllvval(wordint gdid, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin,
               ftnfloat *lat, ftnfloat *lon, wordint n);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdllwdval)(wordint *gdid, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin,
                      ftnfloat *lat, ftnfloat *lon, wordint *n);
wordint c_gdllwdval(wordint gdid, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin,
               ftnfloat *lat, ftnfloat *lon, wordint n);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdxpncf)(wordint *gdin, wordint *i1, wordint *i2, wordint *j1, wordint *j2);
wordint c_gdxpncf(wordint gdin, wordint *i1, wordint *i2, wordint *j1, wordint *j2);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdxysval)(wordint *gdin, ftnfloat *zout, ftnfloat *zin, ftnfloat *x, ftnfloat *y, wordint *n);
wordint c_gdxysval(wordint gdin, ftnfloat *zout, ftnfloat *zin, ftnfloat *x, ftnfloat *y, wordint n);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdxywdval)(wordint *gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, ftnfloat *x, ftnfloat *y, wordint *n);
wordint c_gdxywdval(wordint gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, ftnfloat *x, ftnfloat *y, wordint n);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdxyvval)(wordint *gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, ftnfloat *x, ftnfloat *y, wordint *n);
wordint c_gdxyvval(wordint gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, ftnfloat *x, ftnfloat *y, wordint n);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gduvfwd)(wordint *gdid, ftnfloat *uugdout, ftnfloat *vvgdout,
                     ftnfloat *uullin, ftnfloat *vvllin, ftnfloat *latin, ftnfloat *lonin, wordint *npts);
wordint c_gduvfwd(wordint gdid,  ftnfloat *uugdout, ftnfloat *vvgdout, ftnfloat *uullin, ftnfloat *vvllin,
              ftnfloat *latin, ftnfloat *lonin, wordint npts);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdwdfuv)(wordint *gdid, ftnfloat *uullout, ftnfloat *vvllout, ftnfloat *uuin, ftnfloat *vvin,
              ftnfloat *latin, ftnfloat *lonin, wordint *npts);
wordint c_gdwdfuv(wordint gdid, ftnfloat *uullout, ftnfloat *vvllout, ftnfloat *uuin, ftnfloat *vvin,
              ftnfloat *latin, ftnfloat *lonin, wordint npts);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdxpngd)(wordint *gdin, ftnfloat *zxpnded, ftnfloat *zin);
wordint c_gdxpngd(wordint gdin, ftnfloat *zxpnded, ftnfloat *zin);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdxyfll)(wordint *gdid, ftnfloat *x, ftnfloat *y, ftnfloat *lat, ftnfloat *lon, wordint *n);
wordint c_gdxyfll(wordint gdid, ftnfloat *x, ftnfloat *y, ftnfloat *lat, ftnfloat *lon, wordint n);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdxyzfll)(wordint *gdid, ftnfloat *x, ftnfloat *y, ftnfloat *lat, ftnfloat *lon, wordint *n);
wordint c_gdxyzfll(wordint gdid, ftnfloat *x, ftnfloat *y, ftnfloat *lat, ftnfloat *lon, wordint n);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint c_ezgetgdin();
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint c_ezgetgdout();
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(guval)(wordint *gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin,  ftnfloat *vvin,
                   ftnfloat *x, ftnfloat *y, wordint *n);
wordint c_guval(wordint gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin,  ftnfloat *vvin, ftnfloat *x, ftnfloat *y, wordint n);
/*****************************************************************************/
void c_ezgfllfxy(ftnfloat *lonp, ftnfloat *latp,
                 ftnfloat *lon, ftnfloat *lat,
                 ftnfloat *r, ftnfloat *ri, wordint *npts,
                 ftnfloat *xlat1, ftnfloat *xlon1, ftnfloat *xlat2, ftnfloat *xlon2);
/*****************************************************************************/
void c_ezgfxyfll(ftnfloat *lonp, ftnfloat *latp,
                 ftnfloat *lon, ftnfloat *lat,
                 ftnfloat *r, ftnfloat *ri, wordint *npts,
                 ftnfloat *xlat1, ftnfloat *xlon1, ftnfloat *xlat2, ftnfloat *xlon2);
/*****************************************************************************/
void c_ezgfwfllw(ftnfloat *uullout, ftnfloat *vvllout, ftnfloat *latin, ftnfloat *lonin,
                  ftnfloat *xlatingf, ftnfloat *xloningf,
                  wordint *ni, wordint *nj,
                  char *grtyp, wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4);
/*****************************************************************************/
void  c_ezllwfgfw(ftnfloat *uullout, ftnfloat *vvllout, ftnfloat *latin, ftnfloat *lonin,
                  ftnfloat *xlatingf, ftnfloat *xloningf,
                 wordint *ni,wordint *nj,
                  char *grtyp,wordint *ig1,wordint *ig2,wordint *ig3,wordint *ig4);
/*****************************************************************************/
void c_ez_manageGrillesMemory();
int c_ez_refgrid(int grid_index);

/*****************************************************************************/

void c_ezdefxg(wordint gdid);
void c_ezdefaxes(wordint gdid, ftnfloat *ax, ftnfloat *ay);
wordint c_gdinterp(ftnfloat *zout, ftnfloat *zin, wordint gdin, ftnfloat *x, ftnfloat *y, wordint npts);
//void c_gdkey2rowcol(wordint key, wordint *row, wordint *col);
//void c_gdrowcol2key(wordint *key, wordint row, wordint col);

int f77name(gdsetmask)(int *gdid, int *mask);
int f77name(gdgetmask)(int *gdid, int *mask);
int f77name(ezsint_m)(float *zout, float *zin);
int f77name(ezuvint_m)(float *uuout, float *vvout, float *uuin, float *vvin);
int f77name(ezsint_mdm)(float *zout, int *mask_out, float *zin, int *mask_in);
int f77name(ezuvint_mdm)(float *uuout, float *vvout, int *mask_out, float *uuin, float *vvin, int *mask_in);
int f77name(ezsint_mask)(int *mask_out, int *mask_in);

int c_gdsetmask(int gdid, int *mask);
int c_gdgetmask(int gdid, int *mask);
int c_ezsint_m(float *zout, float *zin);
int c_ezuvint_m(float *uuout, float *vvout, float *uuin, float *vvin);
int c_ezsint_mdm(float *zout, int *mask_out, float *zin, int *mask_in);
int c_ezuvint_mdm(float *uuout, float *vvout, int *mask_out, float *uuin, float *vvin, int *mask_in);
int c_ezsint_mask(int *mask_out, int *mask_in);

#endif
#define _ezfuncdef
