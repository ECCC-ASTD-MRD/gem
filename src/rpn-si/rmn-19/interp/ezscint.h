#ifndef _EZSCINT

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <rpnmacros.h>
#include <string.h>
#include <pthread.h>

#define NMAXGRIDS 32
#define NMAXSETS  NMAXGRIDS * (NMAXGRIDS - 1)
#define NMAXSUBGRIDS 20

#define LAT                         1
#define AX                          2
#define XXX                         4
#define SYMMETRIQUE                64
#define TYPE_EXPANSION            128
#define NEWTON                    256
#define LATLON_PRIME_OK           512
#define SINLATLON_OK             1024
#define ZONES                    2048

#define GRID_GRAPE               1024
#define GRID_CHUNK               1024

#define GRID                      0
#define CLOUD                     1

#define SCALAIRE                  0
#define VECTEUR                   1

#define GLOBALE                   0
#define LOCALE                    1

#define NON                       0
#define OUI                       1

#define VOISIN                    0
#define NEAREST                   0
#define LINEAIRE                  1
#define LINEAR                    1
#define CUBIQUE                   3
#define DISTANCE                  4
#define TRIANGLE                  5
#define LINEAR_AND_NEAREST        6

#define EZ_EXTRAP                 1
#define EZ_NO_EXTRAP              0
#define RIEN                     -1

#define MAXIMUM                   4
#define MINIMUM                   5
#define VALEUR                    6
#define ABORT                     13

#define FICHIER         1
#define MEMOIRE         2

#define NZONES             5
#define DEHORS             0
#define AU_NORD            1
#define AU_SUD             2
#define POLE_NORD          3
#define POLE_SUD           4

#define GLOBAL          0
#define NORD            1
#define SUD             2

#define SYM             1
#define ANTISYM         0

#define CONSERVATIVE    0
#define LIBERAL         1

#define UNDEFINED       -1

#define C_TO_FTN(i,j,ni)  (wordint)((ni) * (j) + i)

#define OUT_OUT         5

#define ABSOLU          0
#define RELATIF         1

#define SWLAT           0
#define SWLON           1
#define DLAT            2
#define DLON            3

#define TD60            0
#define TDGRW           1
#define CLAT            2
#define CLON            3

#define PI              0
#define PJ              1
#define D60             2
#define DGRW            3

#define IG1             0
#define IG2             1
#define IG3             2
#define IG4             3

#define XLAT1           0
#define XLON1           1
#define XLAT2           2
#define XLON2           3

#define CHUNK           128
#define LOG2_CHUNK      7
#define MAX_LOG_CHUNK   12

#define JOINT           1
#define DISJOINT        0

#define YES             1
#define NO              0

typedef struct
{
  wordint npts;               /* nombre de points */
  ftnfloat *x;                /* vecteur de coordonnees x */
  ftnfloat *y;                /* vecteur de coordonnees y */
  wordint *idx;               /* indice du point dans le champ de destination */
} _zone;

typedef struct
{
  wordint n_wts;          /* nombre de poids */
  ftnfloat *xx, *yy;
  ftnfloat *lat, *lon;
  ftnfloat *wts;              /* tableau de poids */
  wordint *mask, *idx;               /* indice du point dans le champ de destination */
} _ygrid;                     /* Grille Y */

typedef struct
{
  wordint flags;
  ftnfloat *lat_rot, *lon_rot, *lat_true, *lon_true;
  ftnfloat *sinlat_rot, *coslat_rot, *sinlon_rot, *coslon_rot;
  ftnfloat *sinlat_true, *coslat_true, *sinlon_true, *coslon_true;
  ftnfloat r[9], ri[9];
} _gemgrid;

typedef struct
{
  wordint flags,yyflags;
  wordint use_sincos_cache;
  wordint gdin;
  wordint next_gdin;
  ftnfloat valpolesud, valpolenord;
  ftnfloat *x, *y;
  wordint *mask_in, *mask_out;
  ftnfloat *yin_maskout,*yan_maskout;
  ftnfloat *yinlat,*yinlon,*yanlat,*yanlon;
  ftnfloat *yin2yin_lat,*yin2yin_lon,*yan2yin_lat,*yan2yin_lon;
  ftnfloat *yin2yan_lat,*yin2yan_lon,*yan2yan_lat,*yan2yan_lon;
  ftnfloat *yin2yin_x,*yin2yin_y,*yan2yin_x,*yan2yin_y;
  ftnfloat *yin2yan_x,*yin2yan_y,*yan2yan_x,*yan2yan_y;
  wordint yincount_yin,yancount_yin,yincount_yan,yancount_yan;
  _gemgrid gemin, gemout;
  _ygrid ygrid;
  _zone zones[NZONES];
}_gridset;

typedef struct
   {
   int child;
   int childOf;
   int parent;
   int niOffset, njOffset;
   int sister;
   int assembly;
   int *parentOf;
   int *sisterOf;
   }_sousgrille;

typedef struct
  {
  wordint  ip1, ip2, ip3;
  wordint date;
  wordint npas, deet, nbits;
  wordint hemisphere,axe_y_inverse;
  ftnfloat xg[16], xgref[16];
  wordint  ig[16], igref[16];
  char fst_grtyp[4],fst_grref[4];
  wordint key_ax, key_ay;
  char nomvarx[8];
  char nomvary[8];
  char typvarx[4];
  char typvary[4];
  char etiketx[16];
  char etikety[16];
  } _fstinfo;

typedef struct
{
  wordint index;
  wordint grid_index;
  wordint flags;
  wordint i1, i2, j1, j2;
  wordint ni,nj;
  wordint nig, nxg;
  wordint ni_ax, nj_ay;
  wordint extension;
  wordint needs_expansion;
  wordint access_count;
  wordint structured;
  wordint next_gd;
  wordint n_gdin, next_gdin, idx_last_gdin, n_gdin_for;
  wordint log_chunk_gdin, log_chunk_gdin_for;
  wordint *gdin_for, *mask;
  wordint nsubgrids,mymaskgrid;
  wordint mymaskgridi0,mymaskgridi1;
  wordint mymaskgridj0,mymaskgridj1;
  wordint *subgrid;
  ftnfloat *lat, *lon;
  ftnfloat *ax, *ay;
  ftnfloat *ncx, *ncy;
  char grtyp[4], grref[4];
  _fstinfo fst;
  _gridset *gset;
} _Grille;


typedef struct
{
  wordint  damage_control;
  wordint  degre_interp;
  wordint  degre_extrap;
  wordint  use_1subgrid;
  wordint  valeur_1subgrid;
  wordint  symmetrie;
  wordint  vecteur;
  wordint  verbose;
  wordint  memory_use;
  wordint  polar_correction;
  wordint  wgt_num;
  wordint  msg_pt_tol;
  wordint  cld_interp_alg;
  wordint  msg_interp_alg;
  ftnfloat msg_gridpt_dist;
  ftnfloat msg_dist_thresh;
  ftnfloat valeur_extrap;
}_groptions;


#endif


#define _EZSCINT
