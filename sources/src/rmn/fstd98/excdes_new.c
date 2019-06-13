/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2014  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/*****************************************************************************
 *                        E X C D E S . C                                    *
 *                                                                           *
 *Auteur                                                                     *
 *   Mario Lépine  -  Mai  2004                                              *
 *   Michel Valin  -  Mars 2014                                              *
 *                                                                           *
 *Objet                                                                      *
 *   Definir pour le logiciel des fichiers standards des criteres            * 
 *   supplementaires de selection en plus de ceux utilises par fstinf.       *
 *                                                                           *
 *   Les criteres de selection a la desire/exclure de editfst peuvent        *
 *   etres definis directement avec des appels aux fonctions ou a l'aide     *
 *   de directives provenant d'un fichier identifie par la variable          *
 *   d'environnement FST_FILTER_FILE                                         *
 *                                                                           *
 *****************************************************************************/
#include <rpnmacros.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <convert_ip.h>

#define MAX_Nlist 40
#define MAX_requetes 20
#define Abs(x) (((x) < 0) ? -(x) : x)
#define Min(x,y) ((x < y) ? x : y)
#if defined (DEBUG)
  #define dbprint fprintf
#else
  #define dbprint ;
#endif

char **fill_string_array(char **string_array, char *farray, int nc, int ns, int rmblanks);
char **allocate_string_array(int ns);
void C_requetes_init(char *f1,char *f2);
void RequetesInit()
{
  C_requetes_init(NULL,NULL);
}
#pragma weak f_requetes_init__=f_requetes_init
#pragma weak f_requetes_init_=f_requetes_init
void f_requetes_init__();
void f_requetes_init_();
void f_requetes_init()
{
  RequetesInit();
}
#ifdef NOTUSED
enum cquoica {unused,entier, reel, deb_fin_entier, deb_fin_reel, deb_fin_entier_delta,
              deb_fin_reel_delta} parametre;
#endif
// enum onveut {desire, exclure} call_mode;

#define DESIRE 1
#define EXCLURE -1
              
#define READLX_DELTA -3
#define READLX_RANGE -2

#define UNUSED 0
#define VALUE 1
#define RANGE 2
#define DELTA 3
#define USED 255

static char *in_use[] = { "unused","value ","range ","delta " };
static int first_R = 0;
static int last_R = MAX_requetes-1;
static FILE *stddebug;
static int bundle_nb = -1;
static int desire_exclure = 1;
typedef struct {
  int in_use;
  int nelm;
  int delta;
  int data[MAX_Nlist];
} DE_int;

typedef struct {
  int in_use;
  int nelm;
  char pdata[MAX_Nlist][13];
} DE_char;
typedef struct {
  int hit;   /* hit count for this request (number of times it was satisfied) */
  int in_use;
  int in_use_supp;
  int exdes;
  DE_char etiquettes;
  DE_char nomvars;
  DE_char typvars;
  DE_int dates;
  DE_int ip1s;
  DE_int ip2s;
  DE_int ip3s;
  int nis;    /* the next 7 items are integer "supplementary" criteria */
  int njs;    /* -1 means it is not to be used */
  int nks;
  int ig1s;
  int ig2s;
  int ig3s;
  int ig4s;
  char grdtyps;    /* "supplementary" criterion grid type, ' ' means do not use */
} DesireExclure;

static DesireExclure  Requests[MAX_requetes];
static int package_not_initialized = 1;
static int DeactivateAllFilters = 0;

#pragma weak fst_deactivate_filters__=fst_deactivate_filters
#pragma weak fst_deactivate_filters_=fst_deactivate_filters
int fst_deactivate_filters__();
int fst_deactivate_filters_();
int fst_deactivate_filters(){
  int old=DeactivateAllFilters;
  DeactivateAllFilters = 1;
  return(old);
}

#pragma weak fst_reactivate_filters__=fst_reactivate_filters
#pragma weak fst_reactivate_filters_=fst_reactivate_filters
int fst_reactivate_filters__();
int fst_reactivate_filters_();
int fst_reactivate_filters(){
  int old=DeactivateAllFilters;
  DeactivateAllFilters = 0;
  return(old);
}

/* put the contents of the request table in text format
 * 
 * use_header != 0   : print separators, otherwise produce "csv" type output
 * filename != NULL  : write output into that file, without separators (ignore use_header)
 * filename == NULL  : output to stdout
 */
void WriteRequestTable(int use_header, char *filename)
{
  int i, j;
  char *sep="\n      ";
  FILE *outfile = NULL;

  if(package_not_initialized) RequetesInit();
  if(filename) {
    outfile = fopen(filename,"w");
    use_header = 0;  /* disregard use_header if filename is not NULL */
  }

  if(outfile == NULL) outfile=stdout;  /* filename NULL or error opening file */
  if(use_header)sep=",";
  for (i=first_R ; i<=last_R ; i++) {
    if(Requests[i].in_use) {
      if(use_header) fprintf(outfile,"=================== Request no %d ===================\n",i);
      if(Requests[i].ip1s.in_use){
        fprintf(outfile,"%2d, '%c', ",i,Requests[i].exdes == DESIRE ? 'D' : 'E');
        fprintf(outfile,"'IP1       ', '%6s', %2d%s %d", in_use[Requests[i].ip1s.in_use],Requests[i].ip1s.nelm,sep,Requests[i].ip1s.data[0]);
        for(j = 1 ; j < Requests[i].ip1s.nelm ; j++) fprintf(outfile,", %d",Requests[i].ip1s.data[j]);
        fprintf(outfile,"\n");
      }
      if(Requests[i].ip2s.in_use){
        fprintf(outfile,"%2d, '%c', ",i,Requests[i].exdes == DESIRE ? 'D' : 'E');
        fprintf(outfile,"'IP2       ', '%6s', %2d%s %d", in_use[Requests[i].ip2s.in_use],Requests[i].ip2s.nelm,sep,Requests[i].ip2s.data[0]);
        for(j = 1 ; j < Requests[i].ip2s.nelm ; j++) fprintf(outfile,", %d",Requests[i].ip2s.data[j]);
        fprintf(outfile,"\n");
      }
      if(Requests[i].ip3s.in_use){
        fprintf(outfile,"%2d, '%c', ",i,Requests[i].exdes == DESIRE ? 'D' : 'E');
        fprintf(outfile,"'IP3       ', '%6s', %2d%s %d", in_use[Requests[i].ip3s.in_use],Requests[i].ip3s.nelm,sep,Requests[i].ip3s.data[0]);
        for(j = 1 ; j < Requests[i].ip3s.nelm ; j++) fprintf(outfile,", %d",Requests[i].ip3s.data[j]);
        fprintf(outfile,"\n");
      }
      if(Requests[i].dates.in_use){
        fprintf(outfile,"%2d, '%c', ",i,Requests[i].exdes == DESIRE ? 'D' : 'E');
        fprintf(outfile,"'Dates     ', '%6s', %2d%s %d", in_use[Requests[i].dates.in_use],Requests[i].dates.nelm,sep,Requests[i].dates.data[0]);
        for(j = 1 ; j < Requests[i].dates.nelm ; j++) fprintf(outfile,", %d",Requests[i].dates.data[j]);
        fprintf(outfile,"\n");
      }
      if(Requests[i].nomvars.in_use){
        fprintf(outfile,"%2d, '%c', ",i,Requests[i].exdes == DESIRE ? 'D' : 'E');
        fprintf(outfile,"'Nomvar    ', '%6s', %2d%s '%-4s'", in_use[Requests[i].nomvars.in_use],Requests[i].nomvars.nelm,sep,Requests[i].nomvars.pdata[0]);
        for(j = 1 ; j < Requests[i].nomvars.nelm ; j++) fprintf(outfile,", '%-4s'",Requests[i].nomvars.pdata[j]);
        fprintf(outfile,"\n");
      }
      if(Requests[i].typvars.in_use){
        fprintf(outfile,"%2d, '%c', ",i,Requests[i].exdes == DESIRE ? 'D' : 'E');
        fprintf(outfile,"'Typvar    ', '%6s', %2d%s '%-2s'", in_use[Requests[i].typvars.in_use],Requests[i].typvars.nelm,sep,Requests[i].typvars.pdata[0]);
        for(j = 1 ; j < Requests[i].typvars.nelm ; j++) fprintf(outfile,", '%-2s'",Requests[i].typvars.pdata[j]);
        fprintf(outfile,"\n");
      }
      if(Requests[i].etiquettes.in_use){
        fprintf(outfile,"%2d, '%c', ",i,Requests[i].exdes == DESIRE ? 'D' : 'E');
        fprintf(outfile,"'Etiket    ', '%6s', %2d%s '%-12s'", in_use[Requests[i].etiquettes.in_use],Requests[i].etiquettes.nelm,sep,Requests[i].etiquettes.pdata[0]);
        for(j = 1 ; j < Requests[i].etiquettes.nelm ; j++) fprintf(outfile,", '%-12s'",Requests[i].etiquettes.pdata[j]);
        fprintf(outfile,"\n");
      }
      if(Requests[i].in_use_supp){
        fprintf(outfile,"%2d, '%c', ",i,Requests[i].exdes == DESIRE ? 'D' : 'E');
        fprintf(outfile,"'Xtra      ', 'value ',  8%s %d, %d, %d, %d, %d, %d, %d, '%c'\n",sep, 
                Requests[i].nis,Requests[i].njs,Requests[i].nks,
                Requests[i].ig1s,Requests[i].ig2s,Requests[i].ig3s,Requests[i].ig4s,Requests[i].grdtyps);
      }
    }  /* if */
  }  /* for */
  if(! use_header) fprintf(outfile," 0\n");
  if(outfile != stdout) fclose(outfile);
}
void DumpRequestTable()
{
  WriteRequestTable(1,NULL);
}
/*
 * basic validation of requests
 * do we exceed the max number of request sets ?
 * do we have a value list that is too long ?
 * check for special case where nelm is negative
 * are we consistent (desire/exclure) for this request set ?
 */
static int ValidateRequestForSet(int set_nb, int des_exc, int nelm, int nelm_lt, char *msg)
{
  if(package_not_initialized) RequetesInit();

  if (set_nb > MAX_requetes-1) {
    fprintf(stderr,"ERROR: (C_select_%s) set_nb=%d > MAX_requetes-1=%d\n",msg,set_nb,MAX_requetes-1);
    return(-1);
  }

  if (nelm > MAX_Nlist) {
    fprintf(stderr,"ERROR: (C_select_%s) nelm=%d > limit=%d\n",msg,nelm,MAX_Nlist);
    return(-2);
  }

  if (nelm <= 0) {
    if (nelm != nelm_lt) {   /* if nelm <= 0, it must be equal to nelm_lt */
      fprintf(stderr,"ERROR: (C_select_%s) nelm invalid = %d\n",msg,nelm);
      return(-3);
    }
  }

  if ((Requests[set_nb].in_use) && (((des_exc == 1) ? DESIRE : EXCLURE) != Requests[set_nb].exdes)) {
    fprintf(stderr,"ERROR: (C_select_%s) des_exc value differs from previous call for set number=%d\n",msg,set_nb);
    fprintf(stderr,"        expected %s, got %s \n",(Requests[set_nb].exdes==DESIRE) ? "desire":"exclure",(des_exc == 1) ? "desire":"exclure");
    return(-4);
  }

  return(0);
}

/*****************************************************************************
 *                    X C _ S E L E C T _ I P 1                              *
 *                                                                           *
 *Objet                                                                      *
 *   Definir la liste des IP1 desires                                        *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  set_nb  numero associe a un groupe d'elements desire/exclure         *
 *  IN  des_exc ==0 : exclure,  !=0 : desire                                 *
 *  IN  iplist  liste de 1 ou plusieurs IP1 a rechercher ou exclure          *
 *      iplist(1) = -1                                     toute valeur de x *
 *      iplist(1) = -2 , iplist(2) = v2  (nelm = 2)              x <= v2     *
 *      iplist(1) = v1 , iplist(2) = -2  (nelm = 2)        v1 <= x           *
 *      iplist(1) = v1 , iplist(2) = -2 , iplist(3) = v2   v1 <= x <= v2     *
 *  IN  nelm    nombre d'elements de la liste                                *
 *                                                                           *
 *****************************************************************************/
int Xc_Select_ip1(int set_nb, int des_exc, void *iplist, int nelm)
{
  int i;
  int *ip_entier=(int *)iplist;
  int valid;

  valid = ValidateRequestForSet(set_nb, des_exc, nelm, 1, "ip1");
  if(valid < 0) goto error ;

  if (ip_entier[0] == -1) nelm = 1;        /* universal value, rest of values if any is irrelevant */
  Requests[set_nb].in_use = USED;          /* set is in use */
  Requests[set_nb].ip1s.in_use = VALUE;    /* item in set is in use */
  Requests[set_nb].ip1s.delta = 0;         /* delta not supported */
  Requests[set_nb].exdes = (des_exc == 1) ? DESIRE : EXCLURE;
  Requests[set_nb].ip1s.nelm = nelm;     /* if range, the value of nelm does not matter, the first 2 values are used */
  Requests[set_nb].ip1s.data[0] = ip_entier[0];  /* first value from list */

  if (nelm == 1 ) return(0);               /* one value, cannot be a range */

  if(ip_entier[1] == READLX_RANGE && ip_entier[3] == READLX_DELTA && nelm == 5) {
    Requests[set_nb].ip1s.data[1] = ip_entier[2];
    Requests[set_nb].ip1s.data[2] = ip_entier[4];
    Requests[set_nb].ip1s.in_use = DELTA;
// printf("RANGE+DELTA detected %d %d %d\n",Requests[set_nb].ip1s.data[0],Requests[set_nb].ip1s.data[1],Requests[set_nb].ip1s.data[2]);
    return(0);
  }
  Requests[set_nb].ip1s.data[2] = READLX_RANGE ;   /* open interval by default for case value @*/
  for (i=1; i<nelm; i++) {                       /* rest of values */
    Requests[set_nb].ip1s.data[i] = ip_entier[i];
  }

  if (ip_entier[0] == READLX_RANGE || ip_entier[1] == READLX_RANGE) { Requests[set_nb].ip1s.in_use = RANGE; Requests[set_nb].ip1s.nelm = 2 ; }
  if (ip_entier[1] == READLX_RANGE) Requests[set_nb].ip1s.data[1] = Requests[set_nb].ip1s.data[2];   /* value @ value  or value @  or @ @ */
  return(0);
error:
  Requests[set_nb].dates.in_use = UNUSED;
  return(-1);
}

/*****************************************************************************
 *                    X C _ S E L E C T _ I P 2                              *
 *                                                                           *
 *Objet                                                                      *
 *   Definir la liste des IP2 desires                                        *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  set_nb  numero associe a un groupe d'elements desire/exclure         *
 *  IN  des_exc 0=exclure 1=desire                                           *
 *  IN  iplist  liste de 1 ou plusieurs IP2 a rechercher ou exclure          *
 *      iplist(1) = -1                                     toute valeur de x *
 *      iplist(1) = -2 , iplist(2) = v2  (nelm = 2)              x <= v2     *
 *      iplist(1) = v1 , iplist(2) = -2  (nelm = 2)        v1 <= x           *
 *      iplist(1) = v1 , iplist(2) = -2 , iplist(3) = v2   v1 <= x <= v2     *
 *  IN  nelm    nombre d'elements de la liste                                *
 *                                                                           *
 *****************************************************************************/
int Xc_Select_ip2(int set_nb, int des_exc, void *iplist, int nelm)
{
  int i;
  int *ip_entier=(int *)iplist;
  int valid;

  valid = ValidateRequestForSet(set_nb, des_exc, nelm, 1, "ip2");
  if(valid < 0) goto error ;

  if (ip_entier[0] == -1) nelm = 1;        /* universal value, rest of values if any is irrelevant */
  Requests[set_nb].in_use = USED;          /* set is in use */
  Requests[set_nb].ip2s.in_use = VALUE;    /* item in set is in use */
  Requests[set_nb].ip2s.delta = 0;         /* delta not supported */
  Requests[set_nb].exdes = (des_exc == 1) ? DESIRE : EXCLURE;
  Requests[set_nb].ip2s.nelm = nelm;     /* if range, the value of nelm does not matter, the first 2 values are used */
  Requests[set_nb].ip2s.data[0] = ip_entier[0];  /* first value from list */

  if (nelm == 1 ) return(0);               /* one value, cannot be a range */

  if(ip_entier[1] == READLX_RANGE && ip_entier[3] == READLX_DELTA && nelm == 5) {
    Requests[set_nb].ip2s.data[1] = ip_entier[2];
    Requests[set_nb].ip2s.data[2] = ip_entier[4];
    Requests[set_nb].ip2s.in_use = DELTA;
// printf("RANGE+DELTA detected %d %d %d\n",Requests[set_nb].ip2s.data[0],Requests[set_nb].ip2s.data[1],Requests[set_nb].ip2s.data[2]);
    return(0);
  }
  Requests[set_nb].ip2s.data[2] = READLX_RANGE ;   /* open interval by default for case value @*/
  for (i=1; i<nelm; i++) {                       /* rest of values */
    Requests[set_nb].ip2s.data[i] = ip_entier[i];
  }

  if (ip_entier[0] == READLX_RANGE || ip_entier[1] == READLX_RANGE) { Requests[set_nb].ip2s.in_use = RANGE;; Requests[set_nb].ip2s.nelm = 2 ; }
  if (ip_entier[1] == READLX_RANGE) Requests[set_nb].ip2s.data[1] = Requests[set_nb].ip2s.data[2];   /* value @ value  or value @  or @ @ */
  return(0);
error:
  Requests[set_nb].dates.in_use = UNUSED;
  return(-1);
}

/*****************************************************************************
 *                    X C _ S E L E C T _ I P 3                              *
 *                                                                           *
 *Objet                                                                      *
 *   Definir la liste des IP3 desires                                        *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  set_nb  numero associe a un groupe d'elements desire/exclure         *
 *  IN  des_exc 0=exclure 1=desire                                           *
 *  IN  iplist  liste de 1 ou plusieurs IP3 a rechercher ou exclure          *
 *      iplist(1) = -1                                     toute valeur de x *
 *      iplist(1) = -2 , iplist(2) = v2  (nelm = 2)              x <= v2     *
 *      iplist(1) = v1 , iplist(2) = -2  (nelm = 2)        v1 <= x           *
 *      iplist(1) = v1 , iplist(2) = -2 , iplist(3) = v2   v1 <= x <= v2     *
 *  IN  nelm    nombre d'elements de la liste                                *
 *                                                                           *
 *****************************************************************************/
int Xc_Select_ip3(int set_nb, int des_exc, void *iplist, int nelm)
{
  int i;
  int *ip_entier=(int *)iplist;
  int valid;

  valid = ValidateRequestForSet(set_nb, des_exc, nelm, 1, "ip3");
  if(valid < 0) goto error ;

  if (ip_entier[0] == -1) nelm = 1;        /* universal value, rest of values if any is irrelevant */
  Requests[set_nb].in_use = USED;          /* set is in use */
  Requests[set_nb].ip3s.in_use = VALUE;    /* item in set is in use */
  Requests[set_nb].ip3s.delta = 0;         /* delta not supported */
  Requests[set_nb].exdes = (des_exc == 1) ? DESIRE : EXCLURE;
  Requests[set_nb].ip3s.nelm = nelm;     /* if range, the value of nelm does not matter, the first 2 values are used */
  Requests[set_nb].ip3s.data[0] = ip_entier[0];  /* first value from list */

  if (nelm == 1 ) return(0);               /* one value, cannot be a range */

  if(ip_entier[1] == READLX_RANGE && ip_entier[3] == READLX_DELTA && nelm == 5) {
    Requests[set_nb].ip3s.data[1] = ip_entier[2];
    Requests[set_nb].ip3s.data[2] = ip_entier[4];
    Requests[set_nb].ip3s.in_use = DELTA;
// printf("RANGE+DELTA detected %d %d %d\n",Requests[set_nb].ip3s.data[0],Requests[set_nb].ip3s.data[1],Requests[set_nb].ip3s.data[2]);
    return(0);
  }
  Requests[set_nb].ip3s.data[2] = READLX_RANGE ;   /* open interval by default for case value @*/
  for (i=1; i<nelm; i++) {                       /* rest of values */
    Requests[set_nb].ip3s.data[i] = ip_entier[i];
  }

  if (ip_entier[0] == READLX_RANGE || ip_entier[1] == READLX_RANGE) { Requests[set_nb].ip3s.in_use = RANGE;; Requests[set_nb].ip3s.nelm = 2 ; }
  if (ip_entier[1] == READLX_RANGE) Requests[set_nb].ip3s.data[1] = Requests[set_nb].ip3s.data[2];   /* value @ value  or value @  or @ @ */
  return(0);
error:
  Requests[set_nb].dates.in_use = UNUSED;
  return(-1);
}

/*****************************************************************************
 *                    X C _ S E L E C T _ D A T E                            *
 *                                                                           *
 *Objet                                                                      *
 *   Definir la liste des DATE desirees                                      *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  set_nb     numero associe a un groupe d'elements desire/exclure      *
 *  IN  des_exc    0=exclure 1=desire                                        *
 *  IN  date_list  liste de 1 ou plusieurs DATE a rechercher ou exclure      *
 *      date_list = [-1]                 (nelm = 1)  toute valeur de x       *
 *      date_list = [-2, v2]             (nelm = 2)  x  <= v2                *
 *      date_list = [v1, -2]             (nelm = 2)  v1 <= x                 *
 *      date_list = [v1, -2, v2]         (nelm = 3)  v1 <= x <= v2           *
 *      date_list = [v1, -2, v2, -3, v3] (nelm = 5)  v1 <= x <= v2 every v3  *
 *  IN  nelm       nombre d'elements de la liste                             *
 *                                                                           *
 *****************************************************************************/
int Xc_Select_date(int set_nb, int des_exc, int *date_list, int nelm)
{
  int i, range=0;
  int valid, delta=0;
  union {
    int i;
    float f;
  }i_or_f;

  valid = ValidateRequestForSet(set_nb, des_exc, nelm, 1, "date");
  if(valid < 0) goto error ;

  if (date_list[0] == -1) nelm = 1;        /* universal value, rest of values if any is irrelevant */
  Requests[set_nb].in_use = USED;          /* set is in use */
  Requests[set_nb].dates.in_use = VALUE;    /* item in set is in use */
  Requests[set_nb].dates.delta = 0;         /* delta not supported */
  Requests[set_nb].exdes = (des_exc == 1) ? DESIRE : EXCLURE;
  Requests[set_nb].dates.nelm = nelm;     /* if range, the value of nelm does not matter, the first 2 values are used */
  Requests[set_nb].dates.data[0] = date_list[0];  /* first value from list */

  if (nelm == 1 ) return(0);               /* one value, cannot be a range */

  Requests[set_nb].dates.nelm = nelm;      /* irrelevant if we have a range of dates */
  for (i=0; i<nelm; i++) {
    Requests[set_nb].dates.data[i] = date_list[i];
    if(date_list[i] == READLX_DELTA) {                 /* DELTA keyword */
      delta++;
      if(delta > 1 || range == 0 ) goto error ;           /* more than one delta keyword or delta encountered before @  */
      if(i >= nelm-1) goto error ;                        /* no value follows delta */
      if(i != 3) goto error ;                             /* must be date1 @ date2 delta hours */
      Requests[set_nb].dates.in_use = DELTA;
      Requests[set_nb].dates.nelm = 3;
      i_or_f.i = date_list[i+1];
      if(i_or_f.i & 0x7F800000) i_or_f.i = i_or_f.f+0.5;  /* real value is in seconds and gets rounded */
      else i_or_f.i = i_or_f.i * 3600;                    /* integer value converted to seconds */
      Requests[set_nb].dates.data[2] = i_or_f.i;
      Requests[set_nb].dates.delta = i_or_f.i;
    }
    if(date_list[i] == READLX_RANGE) {   /* @  keyword   */
      range++;
      if(range>1 || i>1) goto error ;           /* more than one @ keyword or @ keyword too far in line */
//      if(i > 1) goto error ;                    /*  @ cannot be found beyond position 2 */
      Requests[set_nb].dates.in_use = RANGE ;
//      if(i==0) Requests[set_nb].dates.data[0] = 0;              /* @ date .... */
      if(i==1 && nelm>2) {                                            /* date @ date  or date @ delta */
        Requests[set_nb].dates.data[1] = date_list[2] ;
      }
      Requests[set_nb].dates.nelm = 2;
    }
    dbprint(stddebug,"Debug Requests[%d].dates.data[%d] = %d\n",set_nb,i,
           Requests[set_nb].dates.data[i]);
  }
  return(0);
error:
  Requests[set_nb].dates.in_use = UNUSED;
  return(-1);
}

/*****************************************************************************
 *                    X C _ S E L E C T _ E T I Q U E T T E                  *
 *                                                                           *
 *Objet                                                                      *
 *   Definir la liste des ETIQUETTE desirees                                 *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  set_nb     numero associe a un groupe d'elements desire/exclure      *
 *  IN  des_exc    0=exclure 1=desire                                        *
 *  IN  etiq_list  liste de 1 ou plusieurs ETIQUETTE a rechercher ou exclure *
 *  IN  nelm       nombre d'elements de la liste                             *
 *                                                                           *
 *****************************************************************************/
int Xc_Select_etiquette(int set_nb, int des_exc, char *etiq_list[], int nelm)
{
  int i;
  int valid;

  valid = ValidateRequestForSet(set_nb, des_exc, nelm, 1, "etiquette");
  if(valid < 0) goto error ;

  Requests[set_nb].in_use = 1;
  Requests[set_nb].etiquettes.in_use = 1;
  Requests[set_nb].exdes = (des_exc == 1) ? DESIRE : EXCLURE;
  Requests[set_nb].etiquettes.nelm = nelm;
  for (i=0; i<nelm; i++) {
    strncpy(Requests[set_nb].etiquettes.pdata[i],etiq_list[i],13);
    dbprint(stddebug,"Debug Requests[%i].etiquettes.pdata[%i]=%s\n",
           set_nb,i,Requests[set_nb].etiquettes.pdata[i]);
  }
  return(0);
error:
  Requests[set_nb].dates.in_use = UNUSED;
  return(-1);
}

/*****************************************************************************
 *                    X C _ S E L E C T _ N O M V A R                        *
 *                                                                           *
 *Objet                                                                      *
 *   Definir la liste des NOMVAR desirees                                    *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  set_nb     numero associe a un groupe d'elements desire/exclure      *
 *  IN  des_exc    0=exclure 1=desire                                        *
 *  IN  etiq_list  liste de 1 ou plusieurs NOMVAR a rechercher ou exclure    *
 *  IN  nelm       nombre d'elements de la liste                             *
 *                                                                           *
 *****************************************************************************/
int Xc_Select_nomvar(int set_nb, int des_exc, char *nomv_list[], int nelm)
{
  int i;
  int valid;

  valid = ValidateRequestForSet(set_nb, des_exc, nelm, 1, "nomvar");
  if(valid < 0) goto error ;

  Requests[set_nb].in_use = 1;
  Requests[set_nb].nomvars.in_use = 1;
  Requests[set_nb].exdes = (des_exc == 1) ? DESIRE : EXCLURE;
  Requests[set_nb].nomvars.nelm = nelm;
  for (i=0; i<nelm; i++) {
    strncpy(Requests[set_nb].nomvars.pdata[i],nomv_list[i],5);
    dbprint(stddebug,"Debug Requests[%i].nomvars.pdata[%i]=%s\n",
           set_nb,i,Requests[set_nb].nomvars.pdata[i]);
  }
  return(0);
error:
  Requests[set_nb].dates.in_use = UNUSED;
  return(-1);
}

/*****************************************************************************
 *                    X C _ S E L E C T _ S U P P L                          *
 *                                                                           *
 *Objet                                                                      *
 *   Definir la liste des criteres "supplementaires"                         *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  set_nb     numero associe a un groupe d'elements desire/exclure      *
 *  IN  des_exc    0=exclure 1=desire                                        *
 *  IN      ni,nj,nk,ig1,ig2,ig3,ig4,grtyp                                   *
 *          dimension et grille horizontale                                  *
 *          -1 : pas de selection   (ni,nj,nk,ig1,ig2,ig3,ig4)               *
 *          ' ': pas de selection (grtyp)                                    *
 *                                                                           *
 *****************************************************************************/
int Xc_Select_suppl(int set_nb, int des_exc, int ni, int nj, int nk, int ig1, int ig2, int ig3, int ig4, char gtyp)
{
  int i;
  int valid;

  valid = ValidateRequestForSet(set_nb, des_exc, 1, 1, "suppl");
  if(valid < 0) goto error ;

  Requests[set_nb].in_use = 1;
  Requests[set_nb].in_use_supp = 1;
  Requests[set_nb].exdes = (des_exc == 1) ? DESIRE : EXCLURE;
  Requests[set_nb].nis     = ni;
  Requests[set_nb].njs     = nj;
  Requests[set_nb].nks     = nk;
  Requests[set_nb].ig1s    = ig1;
  Requests[set_nb].ig2s    = ig2;
  Requests[set_nb].ig3s    = ig3;
  Requests[set_nb].ig4s    = ig4;
  Requests[set_nb].grdtyps = gtyp;
//  fprintf(stderr,"CRITSUP: ni=%d nj=%d nk=%d ig1=%d ig2=%d ig3=%d ig4=%d gtyp='%c'\n",ni,nj,nj,ig1,ig2,ig2,ig4,gtyp);
  return(0);
error:
  Requests[set_nb].dates.in_use = UNUSED;
  return(-1);
}

/*****************************************************************************
 *                    X C _ S E L E C T _ T Y P V A R                        *
 *                                                                           *
 *Objet                                                                      *
 *   Definir la liste des TYPVAR desirees                                    *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  set_nb     numero associe a un groupe d'elements desire/exclure      *
 *  IN  des_exc    0=exclure 1=desire                                        *
 *  IN  etiq_list  liste de 1 ou plusieurs TYPVAR a rechercher ou exclure    *
 *  IN  nelm       nombre d'elements de la liste                             *
 *                                                                           *
 *****************************************************************************/
int Xc_Select_typvar(int set_nb, int des_exc, char *typv_list[], int nelm)
{
  int i;
  int valid;

  valid = ValidateRequestForSet(set_nb, des_exc, nelm, 1, "typvar");
  if(valid < 0) goto error ;

  Requests[set_nb].in_use = 1;
  Requests[set_nb].typvars.in_use = 1;
  Requests[set_nb].exdes = (des_exc == 1) ? DESIRE : EXCLURE;
  Requests[set_nb].typvars.nelm = nelm;
  for (i=0; i<nelm; i++) {
    strncpy(Requests[set_nb].typvars.pdata[i],typv_list[i],3);
    dbprint(stddebug,"Debug Requests[%i].typvars.pdata[%i]=%s\n",
           set_nb,i,Requests[set_nb].typvars.pdata[i]);
  }
  return(0);
error:
  Requests[set_nb].dates.in_use = UNUSED;
  return(-1);
}

int ReadRequestTable(char *filename)
{
  char line[4096];
  FILE *input=NULL;
  char *cptr;
  int dirset;
  char s1[16], s2[16], s3[16];
  char sar[40][13];
  char *sarp[41];
  int nvalues;
  int a[100];
  int i, j;
  int nx[7];
  char gtyp;
  int rvd; /* range, value, delta */
  int dex ; /* desire / exclure */
  int status;
  union {
    int i;
    float f;
  }i_or_f;

  status = 0;
  if(filename != NULL) input=fopen(filename,"r");
  if(input == NULL) {
    fprintf(stderr,"ERROR: (ReadRequestTable) cannot open directive file '%s'\n",filename);
    return(-1);
  }
  readnext:
  cptr=fgets(line,sizeof(line),input);
//  fprintf(stderr,"header=%s",line);
  dirset=0;
  s1[0]='\000';
  s2[0]='\000';
  s3[0]='\000';
  nvalues=0;
  
  sscanf(line,"%d",&dirset);
  if(dirset==0) {fclose(input) ; return(0); }
  cptr=line;
  while(*cptr != '\'' ) cptr++ ; cptr++ ;
  sscanf(cptr,"%s",s1);
  while(*cptr != ',' ) cptr++ ; while(*cptr != '\'' ) cptr++ ; cptr++ ;
  sscanf(cptr,"%s",s2);
  while(*cptr != ',' ) cptr++ ; while(*cptr != '\'' ) cptr++ ; cptr++ ;
  sscanf(cptr,"%s",s3);
  while(*cptr != ',' ) cptr++ ; cptr++ ;
  sscanf(cptr,"%d",&nvalues);
//  fprintf(stderr,"%d,'%s','%s','%s',%d\n",dirset,s1,s2,s3,nvalues);
  rvd = 0 ; /* not used */
  if(s3[0] == 'v') rvd = VALUE;
  if(s3[0] == 'r') rvd = RANGE;
  if(s3[0] == 'd') rvd = DELTA;
  dex = s1[0] == 'D' ? DESIRE : EXCLURE ;
  
  cptr=fgets(line,sizeof(line),input);
  cptr = line;
//  fprintf(stderr,"reading data\n");

  if(s2[0]=='D') {                            /*  Date */
    sscanf(cptr,"%d",a);
    for (i=1;i<nvalues;i++) {
      while(*cptr != ',' ) cptr++ ; cptr++ ;
      sscanf(cptr,"%d",a+i);
    }
    if(rvd == RANGE && a[0] >= 0 && a[1] >= 0) {
      nvalues=3;
      a[2]=a[1];
      a[1]=READLX_RANGE;
    }
    if(rvd == DELTA) {
      nvalues=5;
      i_or_f.f = a[2];    /* already in seconds, send as a float */
      a[4]=i_or_f.i;
      a[3]=READLX_DELTA;
      a[2]=a[1];
      a[1]=READLX_RANGE;
    }
    status = Xc_Select_date(dirset,dex,a,nvalues);
    if(status != 0) {
      fprintf(stderr,"ERROR: (ReadRequestTable) bad value(s) in date");
      for (i=0;i<nvalues;i++) {
        fprintf(stderr," %d",a[i]);
      }
      fprintf(stderr,"\n");
    }
  }else if(s2[0]=='I') {                       /* IP1/2/3 */
    sscanf(cptr,"%d",a);
    for (i=1;i<nvalues;i++) {
      while(*cptr != ',' ) cptr++ ; cptr++ ;
      sscanf(cptr,"%d",a+i);
    }
    if(rvd == RANGE) {
      nvalues=3;
      a[2]=a[1];
      a[1]=READLX_RANGE;
    }
    if(rvd == DELTA) nvalues = 0;
    if(s2[2]=='1') status = Xc_Select_ip1(dirset,dex,a,nvalues);
    if(s2[2]=='2') status = Xc_Select_ip2(dirset,dex,a,nvalues);
    if(s2[2]=='3') status = Xc_Select_ip3(dirset,dex,a,nvalues);
    if(status != 0) {
      fprintf(stderr,"ERROR: (ReadRequestTable) bad value(s) in ip1/ip2/ip3");
      for (i=0;i<nvalues;i++) {
        fprintf(stderr," %d",a[i]);
      }
      fprintf(stderr,"\n");
    }
  }else if(s2[0]=='X'){                           /* Xtra  */
    sscanf(cptr,"%d",a);
    for (i=1;i<8;i++) {
      while(*cptr != ',' ) cptr++ ; cptr++ ;
      sscanf(cptr,"%d",a+i);
    }
    while(*cptr != '\'' ) cptr++ ; cptr++ ;
    gtyp=*cptr;
    status = Xc_Select_suppl(dirset,dex,a[0],a[1],a[2],a[3],a[4],a[5],a[6],gtyp);
    if(status != 0) {
      fprintf(stderr,"ERROR: (ReadRequestTable) bad value(s) in supplementary criteria");
      for (i=0;i<8;i++) {
        fprintf(stderr," %d",a[i]);
      }
      fprintf(stderr," %c\n",gtyp);
  }
  }else if(s2[0]=='N' || s2[0]=='T' || s2[0]=='E'){   /* Nomvar, Typvar or Etiket */
    for (i=0;i<nvalues;i++) {
      while(*cptr != '\'' ) cptr++ ; cptr++ ;
      j=0;
      while(*cptr != '\'' ) sar[i][j++]=*cptr++; cptr++ ;
      sar[i][j]='\000';
      sarp[i] = sar[i];
    }
    sarp[nvalues] = NULL;
    if(s2[0] == 'N') status = Xc_Select_nomvar(dirset,dex,sarp,nvalues) ;
    if(s2[0] == 'T') status = Xc_Select_typvar(dirset,dex,sarp,nvalues) ;
    if(s2[0] == 'E') status = Xc_Select_etiquette(dirset,dex,sarp,nvalues) ;
    if(status != 0) {
      fprintf(stderr,"ERROR: (ReadRequestTable) bad value(s) in nomvar/typvar/etiekt");
      for (i=0;i<nvalues;i++) {
        fprintf(stderr," '%s'",sar[i]);
      }
      fprintf(stderr,"\n");
    }
  }else{
    fprintf(stderr,"ERROR: (ReadRequestTable) unrecognized type s2='%s' in directive file\n",s2);
    status = -1;
  }
  if(status != 0) {
    fclose(input);
    fprintf(stderr,"ERROR: (ReadRequestTable) status=%d \n",status);
    return -1;
  }
  goto readnext;
  return 0;
}


/*****************************************************************************
 *                      C _ S E L E C T _ G R O U P S E T                    *
 *                                                                           *
 *Objet                                                                      *
 *   Definir un groupe de directives desire/exclure                          *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  first_set_nb   numero de groupe du premier desire/exclure            *
 *  IN  last_set_nb    numero de groupe du dernier desire/exclure            *
 *                                                                           *
 *****************************************************************************/
int C_select_groupset(int first_set_nb, int last_set_nb)
{
  if(package_not_initialized) RequetesInit();
  if ((first_set_nb > MAX_requetes-1) || (last_set_nb > MAX_requetes-1) || (first_set_nb > last_set_nb)) {
    fprintf(stderr,"ERROR: (C_select_groupset) first_set_nb=%d, last_set_nb=%d, MAX allowed=%d\n",
            first_set_nb,last_set_nb,MAX_requetes-1);
    return(-1);
  }
  first_R = first_set_nb;
  last_R = last_set_nb;
return 0; /*CHC/NRC*/
}

/*****************************************************************************
 *                      C _ F I L T R E _ D E S I R E                        *
 *                                                                           *
 *Objet                                                                      *
 *   Change la variable globale pour le mode desire                          *
 *   augmenter bundle_nb de 1                                                *
 *                                                                           *
 *****************************************************************************/
int C_filtre_desire()
{
  
  if(package_not_initialized) RequetesInit();
  bundle_nb++;
  desire_exclure = 1;
  if (bundle_nb > MAX_requetes-1) {
    fprintf(stderr,"ERROR: C_filtre_desire nb=%d > MAX desire/exclure =%d\n",bundle_nb,MAX_requetes-1);
    return(-1);
  }
  printf("desire bundle_nb = %d, desire_exclure = %d\n",bundle_nb,desire_exclure);
return 0; /*CHC/NRC*/
}

/*****************************************************************************
 *                      C _ F I L T R E _ E X C L U R E                      *
 *                                                                           *
 *Objet                                                                      *
 *   Change la variable globale pour le mode exclure                         *
 *   augmenter bundle_nb de 1                                                *
 *                                                                           *
 *****************************************************************************/
int C_filtre_exclure()
{
  
  if(package_not_initialized) RequetesInit();
  bundle_nb++;
  desire_exclure = 0;
  if (bundle_nb > MAX_requetes-1) {
    fprintf(stderr,"ERROR: C_filtre_exclure nb=%d > MAX desire/exclure =%d\n",bundle_nb,MAX_requetes-1);
    return(-1);
  }
  printf("exclure bundle_nb = %d, desire_exclure = %d\n",bundle_nb,desire_exclure);
return 0; /*CHC/NRC*/
}

/*****************************************************************************
 * match_ip                                                                  *
 * is ip0 equivalent to one of the values in data (nelm elements) ?          *
 *                                                                           *
 * return 1 if yes, return 0 if no                                           *
 *****************************************************************************/
static int match_ip(int in_use, int nelm, int *data, int ip1, int translatable)
{
  int amatch=0;
  int mode = -1;   /* IP to P,KIND conversion */
  int i, ip, kind1, kind2, kind3;
  float p0, p1, p2, p3, delta, modulo, error;
  float *fdata = (float *)data;

  if( ! in_use ) return 0;
  if( in_use == RANGE || in_use == DELTA ) {
//fprintf(stderr,"range matching %d\n",ip1);
    if(! translatable) return(0) ;      /* name is not translatable we are done */
    ip = ip1;
    ConvertIp(&ip, &p1, &kind1, mode);  /* convert candidate value */
    ip = data[0];
    if(ip >= 0) {
      ConvertIp(&ip, &p2, &kind2, mode);  /* convert bottom value   */
    }else{                                /* open bottom  @ value   */
      p2 = p1;
      kind2 = kind1;
    }
    ip = data[1];
    if(ip >= 0) {
      ConvertIp(&ip, &p3, &kind3, mode);  /* convert top value   */
    }else{                                /* open top   value @  */
      p3 = p1;
      kind3 = kind1;
    }
    if(p2 > p3){                          // make sure that p2 < p3
      p0 = p3 ; p3 = p2 ; p2 = p0 ;       // swap p2 and p3 (no pint in swapping kinds as they MUST match
    }
    if(kind1 != kind2 || kind1 != kind3) return 0 ;  /* not same kind, no match */
    if(p1 < p2 || p1 > p3) return 0;       /* out of value range, no match */
    if(in_use == RANGE) return 1 ;         /* we have a match if RANGE */
    if(in_use == DELTA) {                  /* we are in range, check delta if there is one */
      if(data[2] <= 0) return 0 ;          /* delta <= 0 not acceptable */
      if(p1 == p2) return 1 ;              /* we have a match, modulo will be zero */
      if(data[2] < 0xFFFFF){               /* max 18 bits for integer values */
        delta = data[2];                   /* integer interval */
      }else{
        delta = fdata[2];                  /* float interval */
      }
      if(delta <= 0) return 0 ; /* delta <= 0 not acceptable */
      modulo = fmodf( (p1 -p2) , delta ) ;
      if(modulo < 0) modulo = -modulo;    /* abs(modulo)  */
// printf("p1=%f p2=%f p3=%f delta=%f modulo=%f ratio=%f\n",p1,p2,p3,delta,modulo,modulo/delta);
      error = modulo/delta;
      if( error < 0.00001 || error > 0.99999 ) return 1 ; /* we have a match */
    }
    return 0;   /* we fell through, we have no match */
  }

//fprintf(stderr,"value matching %d, list of %d elements\n",ip1,nelm);
  if(in_use == VALUE){
    for (i=0 ; i<nelm ; i++) {
//fprintf(stderr,"list[%d]=%d\n",i,data[i]);
      if(ip1 == data[i] || data[i] == -1) return 1;  /* we have a match with raw ip values */
    }
    if(! translatable) return(0) ;/* name is not translatable we are done */
    ip = ip1;
    ConvertIp(&ip, &p1, &kind1, mode);  /* convert candidate value */
    for (i=0 ; i<nelm ; i++) {
      ip = data[i];
      ConvertIp(&ip, &p2, &kind2, mode); /* convert match reference */
      if(kind1 != kind2) continue;       /* not same kind, no match */
      if(p2 == 0 && p1 != p2) continue ; /* if one is 0, both must be, no match */
      delta = 1.0 - p1/p2 ;
      if(delta < 0) delta = -delta;    /* abs(relative error)  */
      if( delta < 0.000001) return 1 ; /* we have a match */
    }
    return 0 ;   /* if we fell through, we have no match */
  }else{
    return 0 ;   /* not a value, no match */
  }
}

/*****************************************************************************
 *                      C _ F S T _ M A T C H _ P A R M                      *
 *                                                                           *
 *Objet                                                                      *
 *   Verifier si l'enregistrement courant correspond aux criteres demandes   *
 *   Retourne 1 si l'enregistrement correspond, sinon 0                      *
 *   0 est aussi retourne si l'enregistrement correspond a un a exclure      *
 *                                                                           *
 *   en bref:  on retourne 1 si on accepte cet enregistrement, 0 sinon       *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  handle     handle de l'enregistrement courant                        *
 * IN   datevalid,...,ig4                                                    *
 *      ce qui a ete obtenu de fstprm pour handle                            *
 *                                                                           *
 *****************************************************************************/
int C_fstmatch_parm(int handle, int datevalid, int ni, int nj, int nk,
                     int ip1, int ip2, int ip3, char *typvar, char *nomvar, char *etiket,
                     char *grtyp, int ig1, int ig2, int ig3, int ig4)
{
  int i, j, set_nb, last_in_use;
  int mode=-1,flag=0;
  int ip_kind, amatch=0;
  int debut, fin, date;
  double diff_deb, diff_fin, delta8, remainder;
  char *desire_exclure;
  int translatable;
  int supp_ok;       /* supplementary parameters match */
  int nb_desire=0;

  if(package_not_initialized) {
    fprintf(stderr,"INFO: C_fstmatch_parm, initializing request tables \n");
    RequetesInit();
  }
//  if (! Requests[first_R].in_use) return(1);        /* aucune requete desire ou exclure */

#ifdef NOTUSED
  dbprint(stddebug,"Debug C_fst_match_req fstprm date=%d ip1=%d ip2=%d ip3=%d nomvar-->%s<-- typvar-->%s<-- etiket-->%s<--\n",
                   date,ip1,ip2,ip3,nomvar,typvar,etiket);
#endif
  translatable = FstCanTranslateName(nomvar) ;
  date = datevalid;
//  for (set_nb=first_R; set_nb <= last_R ; set_nb++) fprintf(stderr," %d",Requests[set_nb].in_use) ; fprintf(stderr,"\n");

  for (set_nb=first_R; (set_nb <= last_R) ; set_nb++) {
    if(Requests[set_nb].in_use == 0)continue ;

    /* process supplementary parameters if any right here */
    amatch = 0;
    if (Requests[set_nb].exdes == DESIRE) nb_desire++;
    desire_exclure = (Requests[set_nb].exdes == DESIRE)  ? "desire" : "exclure";
//fprintf(stderr,"matching request set %d\n",set_nb);
    Supplements:
      if (Requests[set_nb].in_use_supp) {   /* les criteres supplementires sont globaux */
        if( (Requests[set_nb].ig1s != ig1 && Requests[set_nb].ig1s != -1) 
        ||  (Requests[set_nb].ig2s != ig2 && Requests[set_nb].ig2s != -1)
        ||  (Requests[set_nb].ig3s != ig3 && Requests[set_nb].ig3s != -1)
        ||  (Requests[set_nb].ig4s != ig2 && Requests[set_nb].ig4s != -1)
        ||  (Requests[set_nb].nis  !=  ni &&  Requests[set_nb].nis != -1)
        ||  (Requests[set_nb].njs  !=  nj &&  Requests[set_nb].njs != -1)
        ||  (Requests[set_nb].nks  !=  nk &&  Requests[set_nb].nks != -1)
        ||  (Requests[set_nb].grdtyps  !=  grtyp[0] && Requests[set_nb].grdtyps != ' ') ) continue;  /* requete non satisfaite pour criteres supplementaires */
      }
      amatch = 1;   /* requete satisfaite jusqu'ici si criteres supplementaires actifs et OK */
    Etiquettes:
      if (Requests[set_nb].etiquettes.in_use) {
        amatch = 0;
        if(Requests[set_nb].etiquettes.pdata[0][0] == ' ') amatch = 1;
        dbprint(stddebug,"Debug C_fst_match_req verifie etiquettes du fichier=%s set_nb=%d\n",etiket,set_nb);
        for (i=0; i < Requests[set_nb].etiquettes.nelm && amatch == 0; i++) {
          if (strncmp(Requests[set_nb].etiquettes.pdata[i],etiket,Min(12,strlen(Requests[set_nb].etiquettes.pdata[i]))) == 0) {
            amatch = 1;       /* requete satisfaite jusqu'ici */
            dbprint(stddebug,"Debug C_fst_match_req match %s\n",desire_exclure);
          }
        }
//fprintf(stderr,"matching etiket, amatch=%d\n",amatch);
        if (amatch == 0) continue;  /* requete non satisfaite pour etiquettes */
      }

    Nomvars:
      if (Requests[set_nb].nomvars.in_use) {
        amatch = 0;
        if(Requests[set_nb].nomvars.pdata[0][0] == ' ') amatch = 1;
        dbprint(stddebug,"Debug C_fst_match_req verifie nomvars du fichier=%s set_nb=%d\n",nomvar,set_nb);
        for (i=0; i < Requests[set_nb].nomvars.nelm && amatch == 0; i++)
          if (strncmp(Requests[set_nb].nomvars.pdata[i],nomvar,Min(4,strlen(Requests[set_nb].nomvars.pdata[i]))) == 0) {
            amatch = 1;       /* requete satisfaite jusqu'ici */
            dbprint(stddebug,"Debug C_fst_match_req match %s\n",desire_exclure);
          }
//fprintf(stderr,"matching nomvar, amatch=%d\n",amatch);
        if (amatch == 0) continue;  /* requete non satisfaite pour nomvars */
      }

    Typvars:
      if (Requests[set_nb].typvars.in_use) {
        amatch = 0;
        if(Requests[set_nb].typvars.pdata[0][0] == ' ') amatch = 1;
        dbprint(stddebug,"Debug C_fst_match_req verifie typvars set_nb=%d\n",set_nb);
        for (i=0; i < Requests[set_nb].typvars.nelm && amatch == 0; i++)
          if (strncmp(Requests[set_nb].typvars.pdata[i],typvar,Min(2,strlen(Requests[set_nb].typvars.pdata[i]))) == 0) {
            amatch = 1;       /* requete satisfaite jusqu'ici */
            dbprint(stddebug,"Debug C_fst_match_req match %s\n",desire_exclure);
          }
//fprintf(stderr,"matching typvar, amatch=%d\n",amatch);
        if (amatch == 0) continue;  /* requete non satisfaite pour typvars */
      }

    Dates:
      if (Requests[set_nb].dates.in_use) {
          amatch = 0;
          switch (Requests[set_nb].dates.in_use) {

          case VALUE:
            dbprint(stddebug,"Debug C_fst_match_req verifie dates entier set_nb=%d\n",set_nb);
            if(Requests[set_nb].dates.data[0] == -1) amatch = 1;
            for (i=0; i < Requests[set_nb].dates.nelm && amatch == 0; i++)
              if (Requests[set_nb].dates.data[i] == date) {
                amatch = 1;
                dbprint(stddebug,"Debug C_fst_match_req match %s\n",desire_exclure);
              }
//fprintf(stderr,"matching date value, amatch=%d\n",amatch);
            if(amatch == 0) continue;  /* rien trouve qui satisfasse la requete pour dates */
            break;   /* requete satisfaite jusqu'ici */

          case RANGE:
            if (Requests[set_nb].dates.data[0] <= 0)
              debut = date;
            else {
              debut = Requests[set_nb].dates.data[0];
              }
            if (Requests[set_nb].dates.data[1] <= 0)
              fin = date;
            else
              fin = Requests[set_nb].dates.data[1];
            dbprint(stddebug,"Debug C_fst_match_req verifie dates debut=%d fin=%d set_nb=%d\n",debut,fin,set_nb);
            f77name(difdatr)(&date,&debut,&diff_deb);
            f77name(difdatr)(&date,&fin,&diff_fin);
            dbprint(stddebug,"Debug diff_deb=%f diff_fin=%f\n",diff_deb,diff_fin);
            if ((diff_deb >= 0.) && (diff_fin <= 0.)) {
              amatch = 1;
              dbprint(stddebug,"Debug C_fst_match_req match %s\n",desire_exclure);
            }
//fprintf(stderr,"matching date range, amatch=%d\n",amatch);
            if(amatch == 0) continue;  /* rien trouve qui satisfasse la requete desire/exclure */
            break;   /* requete satisfaite jusqu'ici */

          case DELTA:
            dbprint(stddebug,"Debug C_fst_match_req verifie dates debut fin delta set_nb=%d\n",set_nb);
            if (Requests[set_nb].dates.data[0] <= 0)
              debut = date;
            else {
              debut = Requests[set_nb].dates.data[0];
              }
            if (Requests[set_nb].dates.data[1] <= 0)
              fin = date;
            else
              fin = Requests[set_nb].dates.data[1];
            delta8 = Requests[set_nb].dates.delta;
            delta8 /= 3600.0;  /* put delta in hours */
            dbprint(stddebug,"Debug C_fst_match_req verifie dates debut=%d fin=%d delta=%f\n",debut,fin,delta8);
            f77name(difdatr)(&date,&debut,&diff_deb);
            f77name(difdatr)(&date,&fin,&diff_fin);
            remainder = fmod(diff_deb,delta8);
            dbprint(stddebug,"Debug diff_deb=%f diff_fin=%f modulo=%f\n",diff_deb,diff_fin,remainder);
            if ((diff_deb >= 0.) && (diff_fin <= 0.) && (remainder <= (5.0/3600.))) {
              amatch = 1;
              dbprint(stddebug,"Debug C_fst_match_req match %s\n",desire_exclure);
            }
//fprintf(stderr,"matching date delta, amatch=%d\n",amatch);
            if(amatch == 0) continue;  /* rien trouve qui satisfait la requete desire/exclure */
            break;   /* requete satisfaite jusqu'ici */

          default:
            fprintf(stderr,"ERROR: (C_fst_match_req) invalid Requests[%d].dates.in_use=%d\n",
                  set_nb,Requests[set_nb].dates.in_use);
            return(0);
            break;     /* this never gets executed */

          } /* end switch */
      }

    Ip1s:
      if (Requests[set_nb].ip1s.in_use) {
        amatch = (Requests[set_nb].ip1s.data[0] == -1) ? 1 : 0 ;
        dbprint(stddebug,"Debug C_fst_match_req verifie ip1s set_nb=%d\n",set_nb);
        if( amatch == 0)
          amatch = match_ip(Requests[set_nb].ip1s.in_use, Requests[set_nb].ip1s.nelm, Requests[set_nb].ip1s.data, ip1, translatable);
//fprintf(stderr,"%s matching ip1=%d, amatch=%d\n",in_use[Requests[set_nb].ip1s.in_use],ip1,amatch);
        if(amatch == 0) continue ;  /* requete non satisfaite pour ip1 */
      }

    Ip2s:
      if (Requests[set_nb].ip2s.in_use) {
        amatch = (Requests[set_nb].ip2s.data[0] == -1) ? 1 : 0 ;
        dbprint(stddebug,"Debug C_fst_match_req verifie ip2s set_nb=%d\n",set_nb);
        if( amatch == 0)
          amatch = match_ip(Requests[set_nb].ip2s.in_use, Requests[set_nb].ip2s.nelm, Requests[set_nb].ip2s.data, ip2, translatable);
//printf(stderr,"%s matching ip2=%d, amatch=%d\n",in_use[Requests[set_nb].ip2s.in_use],ip2,amatch);
        if(amatch == 0) continue ;  /* requete non satisfaite pour ip2 */
      }

    Ip3s:
      if (Requests[set_nb].ip3s.in_use) {
        amatch = (Requests[set_nb].ip3s.data[0] == -1) ? 1 : 0 ;
        dbprint(stddebug,"Debug C_fst_match_req verifie ip3s set_nb=%d\n",set_nb);
        if( amatch == 0)
          amatch = match_ip(Requests[set_nb].ip3s.in_use, Requests[set_nb].ip3s.nelm, Requests[set_nb].ip3s.data, ip3, translatable);
//fprintf(stderr,"%s matching ip3=%d, amatch=%d\n",in_use[Requests[set_nb].ip3s.in_use],ip3,amatch);
        if(amatch == 0) continue ;  /* requete non satisfaite pour ip3 */
      }

    Fin:
      if (amatch == 1) {
        dbprint(stddebug,"Debug C_fst_match_req fin requete %s satisfaite, handle=%d \n",desire_exclure,handle);
        Requests[set_nb].hit++ ;                                    /* add one to this request's hit count */
        return((Requests[set_nb].exdes == DESIRE) ? set_nb+1 : 0 );  /* requete desire satisfaite */
      }      
    /* Next:                  verifier le prochain desire/exclure */

  } /* end for */

  last_in_use = last_R;
  while ((last_in_use > first_R) && (! Requests[last_in_use].in_use))
    last_in_use--;
/*  fprintf(stderr,"Debug+ last_in_use=%d\n",last_in_use); */
/*  rien trouve qui satisfasse les requetes */
  if ((Requests[last_in_use].exdes == DESIRE) && (Requests[last_in_use].in_use))
    return(0);  /* ne satisfait pas la requete desire */
  else
    return(nb_desire==0 ? 1:0);  /* rien a exclure */

} /* end fst_match_req */

/*****************************************************************************
 *                      C _ F S T _ M A T C H _ R E Q                        *
 *                                                                           *
 *Objet                                                                      *
 *   Verifier si l'enregistrement courant correspond aux criteres demandes   *
 *   Retourne 1 si l'enregistrement correspond, sinon 0                      *
 *   0 est aussi retourne si l'enregistrement correspond a un a exclure      *
 *                                                                           *
 *   en bref:  on retourne 1 si on accepte cet enregistrement, 0 sinon       *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  handle     handle de l'enregistrement a valider                      *
 *                                                                           *
 *****************************************************************************/
int C_fstmatch_req(int handle)
{
  int ier;
  int ni, nj, nk, dateo, deet, npas, ip1, ip2, ip3, ig1, ig2, ig3, ig4;
  int nbits, swa, ubc, lng, dltf, datevalid, xtra1, xtra2, xtra3, datyp;
  char etiket[13]={' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ','\0'};
  char typvar[3]={' ',' ','\0'};
  char nomvar[5]={' ',' ',' ',' ','\0'};
  char grtyp[2]={' ','\0'};
  int status;
//  return(1);
  if(DeactivateAllFilters) return(1);  /* filtering deactivated */
  if(package_not_initialized) RequetesInit();
  ier = c_fstprm(handle,&dateo,&deet,&npas,&ni,&nj,&nk,
                     &nbits,&datyp,&ip1,&ip2,&ip3,typvar,
                     nomvar,etiket,grtyp,&ig1,&ig2,&ig3,&ig4,&swa,&lng,
                     &dltf,&ubc,&datevalid,&xtra2,&xtra3);

  if (ier < 0) return(0);
  status = C_fstmatch_parm(handle, datevalid, ni, nj, nk, ip1, ip2, ip3,
                            typvar, nomvar, etiket, grtyp, ig1, ig2, ig3, ig4) ;
//  fprintf(stderr,"C_fstmatch_req/fstd status=%d\n",status);
  return status ;
}

int C_fst_match_req(int handle)
{
  int ier;
  int ni, nj, nk, dateo, deet, npas, ip1, ip2, ip3, ig1, ig2, ig3, ig4;
  int nbits, swa, ubc, lng, dltf, datevalid, xtra1, xtra2, xtra3, datyp;
  char etiket[13]={' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ','\0'};
  char typvar[3]={' ',' ','\0'};
  char nomvar[5]={' ',' ',' ',' ','\0'};
  char grtyp[2]={' ','\0'};
  int status;

  if(DeactivateAllFilters) return(1);  /* filtering deactivated */
  if(package_not_initialized) RequetesInit();
  ier = c_fstprm(handle,&dateo,&deet,&npas,&ni,&nj,&nk,
                     &nbits,&datyp,&ip1,&ip2,&ip3,typvar,
                     nomvar,etiket,grtyp,&ig1,&ig2,&ig3,&ig4,&swa,&lng,
                     &dltf,&ubc,&datevalid,&xtra2,&xtra3);

  if (ier < 0) return(0);
  status = C_fstmatch_parm(handle, datevalid, ni, nj, nk, ip1, ip2, ip3,
                            typvar, nomvar, etiket, grtyp, ig1, ig2, ig3, ig4) ;
//  fprintf(stderr,"C_fst_match_req status=%d\n",status);
  return status ;
}

/*****************************************************************************
 *                     L E S   I N T E R F A C E S                           *
 *                                                                           *
 *Objet                                                                      *
 *   Les interfaces Fortran a C ainsi que les interfaces aux                 *
 *   fonctions C a liste d'arguments raccourcie aux fonctions                *
 *   a liste d'arguments etendue                                             *
 *                                                                           *
 *****************************************************************************/

int C_select_ip1(void *iplist, int nelm)
{
/*CHC/NRC*/
	return Xc_Select_ip1(bundle_nb, desire_exclure, iplist, nelm);
}

int C_select_ip2(void *iplist, int nelm)
{
/*CHC/NRC*/
	return Xc_Select_ip2(bundle_nb, desire_exclure, iplist, nelm);
}

int C_select_ip3(void *iplist, int nelm)
{
/*CHC/NRC*/
	return Xc_Select_ip3(bundle_nb, desire_exclure, iplist, nelm);
}

int C_select_date(int set_nb, int des_exc, int *date_list, int nelm)
{
/*CHC/NRC*/
	return Xc_Select_date(bundle_nb, desire_exclure, date_list, nelm);
}

int C_select_etiquette(char *etiq_list[], int nelm)
{
/*CHC/NRC*/
	return Xc_Select_etiquette(bundle_nb, desire_exclure, etiq_list, nelm);
}

int C_select_nomvar(char *nomv_list[], int nelm)
{
/*CHC/NRC*/
	return Xc_Select_nomvar(bundle_nb, desire_exclure, nomv_list, nelm);
}

int C_select_typvar(char *typv_list[], int nelm)
{
/*CHC/NRC*/
    return Xc_Select_typvar(bundle_nb, desire_exclure, typv_list, nelm);
}

int C_select_suppl(int ni, int nj, int nk, int ig1, int ig2, int ig3, int ig4, char gtyp)
{
/*CHC/NRC*/
    return Xc_Select_suppl(bundle_nb, desire_exclure, ni, nj, nk, ig1, ig2, ig3, ig4, gtyp);
}

int f77name(f_select_ip1)(void *iplist, int *nelm)
{
  return(Xc_Select_ip1(bundle_nb, desire_exclure,iplist,*nelm));
}

int f77name(f_select_ip2)(void *iplist, int *nelm)
{
  return(Xc_Select_ip2(bundle_nb, desire_exclure,iplist,*nelm));
}

int f77name(f_select_ip3)(void *iplist, int *nelm)
{
  return(Xc_Select_ip3(bundle_nb, desire_exclure,iplist,*nelm));
}

int f77name(f_select_date)(int *date_list, int *nelm)
{
  return(Xc_Select_date(bundle_nb, desire_exclure,date_list,*nelm));
}

int f77name(f_select_suppl)(int *ni, int *nj, int *nk, int *ig1, int *ig2, int *ig3, int *ig4, char *gtyp, F2Cl flng)
{
//fprintf(stderr,"f_select_suppl: ni=%d nj=%d nk=%d ig1=%d ig2=%d ig3=%d ig4=%d gtyp='%c'\n",*ni,*nj,*nj,*ig1,*ig2,*ig3,*ig4,*gtyp);
    return Xc_Select_suppl(bundle_nb, desire_exclure, *ni, *nj, *nk, *ig1, *ig2, *ig3, *ig4, *gtyp);
}

int Xf_Select_etiquette(int set_nb, int des_exc, char *etiq_list, int nelm, int flng)
{
  char **string_array;
  int i,ier;
  dbprint(stddebug,"Debug desire_etiquette lng=%d nelm=%d\n",flng,nelm);
  string_array = fill_string_array(allocate_string_array(nelm),etiq_list,flng,nelm,0);
  for (i=0; i < nelm; i++)
      dbprint(stddebug,"Debug string_array[%d]-->%s<--\n",i,string_array[i]);
  ier = Xc_Select_etiquette(set_nb, des_exc,string_array,nelm);
  free_string_array(string_array);
  return(ier);
}

int f77name(f_select_etiquette)(char *etiq_list, int *nelm, F2Cl flng)
{
  char **string_array;
  int i,ier,lng=flng;
  dbprint(stddebug,"Debug desire_etiquette lng=%d nelm=%d\n",lng,*nelm);
  string_array = fill_string_array(allocate_string_array(*nelm),etiq_list,lng,*nelm,0);
  for (i=0; i < *nelm; i++)
      dbprint(stddebug,"Debug string_array[%d]-->%s<--\n",i,string_array[i]);
  ier = Xc_Select_etiquette(bundle_nb, desire_exclure,string_array,*nelm);
  free_string_array(string_array);
  return(ier);
}

int Xf_Select_nomvar(int set_nb, int des_exc, char *nomv_list, int nelm, int flng)
{
  char **string_array;
  int i,ier;
  dbprint(stddebug,"Debug desire_etiquette lng=%d nelm=%d\n",flng,nelm);
  string_array = fill_string_array(allocate_string_array(nelm),nomv_list,flng,nelm,0);
  for (i=0; i < nelm; i++)
      dbprint(stddebug,"Debug string_array[%d]-->%s<--\n",i,string_array[i]);
  ier = Xc_Select_nomvar(set_nb, des_exc,string_array,nelm);
  free_string_array(string_array);
  return(ier);
}

int f77name(f_select_nomvar)(char *nomv_list, int *nelm, F2Cl flng)
{
  char **string_array;
  int i,ier,lng=flng;
  dbprint(stddebug,"Debug desire_nomvar lng=%d nelm=%d\n",lng,*nelm);
  string_array = fill_string_array(allocate_string_array(*nelm),nomv_list,lng,*nelm,0);
  for (i=0; i < *nelm; i++)
      dbprint(stddebug,"Debug string_array[%d]-->%s<--\n",i,string_array[i]);
  ier = Xc_Select_nomvar(bundle_nb, desire_exclure,string_array,*nelm);
  free_string_array(string_array);
  return(ier);
}

int Xf_Select_typvar(int set_nb, int des_exc, char *typv_list, int nelm, int flng)
{
  char **string_array;
  int i,ier;
  dbprint(stddebug,"Debug desire_etiquette lng=%d nelm=%d\n",flng,nelm);
  string_array = fill_string_array(allocate_string_array(nelm),typv_list,flng,nelm,0);
  for (i=0; i < nelm; i++)
      dbprint(stddebug,"Debug string_array[%d]-->%s<--\n",i,string_array[i]);
  ier = Xc_Select_typvar(set_nb, des_exc,string_array,nelm);
  free_string_array(string_array);
  return(ier);
}

int f77name(f_select_typvar)(char *typv_list, int *nelm, F2Cl flng)
{
  char **string_array;
  int i,ier,lng=flng;
  dbprint(stddebug,"Debug desire_typvar lng=%d nelm=%d\n",lng,*nelm);
  string_array = fill_string_array(allocate_string_array(*nelm),typv_list,lng,*nelm,0);
  for (i=0; i < *nelm; i++)
      dbprint(stddebug,"Debug string_array[%d]-->%s<--\n",i,string_array[i]);
  ier = Xc_Select_typvar(bundle_nb, desire_exclure,string_array,*nelm);
  free_string_array(string_array);
  return(ier);
}

int f77name(f_filtre_desire)()
{
/*CHC/NRC*/
	return C_filtre_desire();
}

int f77name(f_filtre_exclure)()
{
/*CHC/NRC*/
	return C_filtre_exclure();
}
  
int f77name(f_fst_match_req)(int *handle)
{
  return(C_fstmatch_req(*handle));
}

/*****************************************************************************
 *                      C _ R E Q U E T E S _ R E S E T                      *
 *                                                                           *
 *Objet                                                                      *
 *   Reinitialiser en tout ou en partie une requete                          *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  set_nb     numero associe a un groupe requete                        *
 *  IN  nomvars    0=initialise a zero   -1=garder tel quel                  *
 *  IN  typvars    0=initialise a zero   -1=garder tel quel                  *
 *  IN  etikets    0=initialise a zero   -1=garder tel quel                  *
 *  IN  dates      0=initialise a zero   -1=garder tel quel                  *
 *  IN  ip1s       0=initialise a zero   -1=garder tel quel                  *
 *  IN  ip2s       0=initialise a zero   -1=garder tel quel                  *
 *  IN  ip3s       0=initialise a zero   -1=garder tel quel                  *
 *                                                                           *
 *****************************************************************************/
int C_requetes_reset(int set_nb, int nomvars, int typvars, int etikets, int dates,
                    int ip1s, int ip2s, int ip3s)
{
  int j;

  if (set_nb > MAX_requetes-1) {
    fprintf(stderr,"ERROR: (C_requetes_reset) set_nb=%d > MAX allowed=%d\n",set_nb,MAX_requetes-1);
    return(-1);
  }

  Requests[set_nb].in_use = UNUSED;
  Requests[set_nb].hit = 0;
  Requests[set_nb].in_use_supp = UNUSED;
  Requests[set_nb].exdes = -1;
  if (nomvars != -1) {
    Requests[set_nb].nomvars.in_use = UNUSED;
    Requests[set_nb].nomvars.nelm = 0;
    for (j=0; j < MAX_Nlist; j++)
      strcpy(Requests[set_nb].nomvars.pdata[j],"    ");
  }
  if (typvars != -1) {
    Requests[set_nb].typvars.in_use = UNUSED;
    Requests[set_nb].typvars.nelm = 0;
    for (j=0; j < MAX_Nlist; j++)
      strcpy(Requests[set_nb].typvars.pdata[j],"  ");
  }
  if (etikets != -1) {
    Requests[set_nb].etiquettes.in_use = UNUSED;
    Requests[set_nb].etiquettes.nelm = 0;
    for (j=0; j < MAX_Nlist; j++)
      strcpy(Requests[set_nb].etiquettes.pdata[j],"            ");
  }
  if (dates != -1) {
    Requests[set_nb].dates.in_use = UNUSED;
    Requests[set_nb].dates.nelm = 0;
    for (j=0; j < MAX_Nlist; j++)
      Requests[set_nb].dates.data[j]=0;
  }
  if (ip1s != -1) {
    Requests[set_nb].ip1s.in_use = UNUSED;
    Requests[set_nb].ip1s.nelm = 0;
    for (j=0; j < MAX_Nlist; j++)
       Requests[set_nb].ip1s.data[j]=0;
  }
  if (ip2s != -1) {
    Requests[set_nb].ip2s.in_use = UNUSED;
    Requests[set_nb].ip2s.nelm = 0;
    for (j=0; j < MAX_Nlist; j++)
       Requests[set_nb].ip2s.data[j]=0;
  }
  if (ip3s != -1) {
    Requests[set_nb].ip3s.in_use = UNUSED;
    Requests[set_nb].ip3s.nelm = 0;
    for (j=0; j < MAX_Nlist; j++)
       Requests[set_nb].ip3s.data[j]=0;
  }
  Requests[set_nb].in_use_supp = UNUSED;
  Requests[set_nb].nis = 0;
  Requests[set_nb].njs = 0;
  Requests[set_nb].nks = 0;
  Requests[set_nb].ig1s = 0;
  Requests[set_nb].ig2s = 0;
  Requests[set_nb].ig3s = 0;
  Requests[set_nb].ig4s = 0;
  Requests[set_nb].grdtyps = ' ';

  return(0);
}
/*****************************************************************************
 *                      C _ R E Q U E T E S _ I N I T                        *
 *                                                                           *
 *Objet                                                                      *
 *   Initialiser le package de requetes                                      *
 *                                                                           *
 *****************************************************************************/
void C_requetes_init(char *requetes_filename, char *debug_filename)
{
  int i,j,ier;

/*  debug_filename = getenv("DEBUGFILE"); */
  if (debug_filename != NULL)
    stddebug = fopen(debug_filename,"w");
  else
    stddebug = fopen("/dev/null","w");

  first_R = 0;
  last_R = MAX_requetes-1;
  bundle_nb = -1;
  desire_exclure = 1;
  for (i=0; i<MAX_requetes; i++) {
    C_requetes_reset(i,0,0,0,0,0,0,0);
  }
#ifdef NO_LONGER_USED
  /* requetes_filename = getenv("FST_FILTER_FILE"); */
  if (requetes_filename != NULL)
    ier = C_requetes_read_file(requetes_filename);
#endif
  dbprint(stderr,"INFO: (C_requetes_init) request table initialized \n");
  package_not_initialized = 0;
}

/*****************************************   start of deprecated code *****************************************/
#ifdef NO_LONGER_USED

/*****************************************************************************
 *                      D I R E C T I V E S _ I P 1 2 3                      *
 *                                                                           *
 *Objet                                                                      *
 *   Traiter les directives desires pour ip1 ip2 et ip3                      *
 *                                                                           *
 *****************************************************************************/
int Directive_ip123(int argc , char **argv,char cmd_strt,  void *func(), char *Private_Data_2)
{
  int i, set_nb, des_exc, nelm, kind, offset, n, range;
  int list_entier[100];
  float list_reel[100];

/*  for (i=0;i<=argc;i++)
    printf("argument # %d-->%s<-- \n",i,argv[i]); */

  set_nb = bundle_nb;    
	des_exc = desire_exclure;  
  if ((strncasecmp(argv[1],"asis",3) == 0) || (strncasecmp(argv[1],"telquel",3) == 0))
  	kind = -1;    
  else if ((strncasecmp(argv[1],"height",3) == 0) || (strncasecmp(argv[1],"hauteur",3) ==0))
  	kind = 0;
  else if (strncasecmp(argv[1],"sigma",3) == 0)
  	kind = 1;
  else if (strncasecmp(argv[1],"pression",3) == 0)
  	kind = 2;
  else if (strncasecmp(argv[1],"arbitraire",3) == 0)
  	kind = 3;
  else if ((strncasecmp(argv[1],"ground",3) == 0) || (strncasecmp(argv[1],"sol",3) ==0))
  	kind = 4;    
  else if (strncasecmp(argv[1],"hybride",3) == 0)
  	kind = 5;
  else if (strncasecmp(argv[1],"theta",3) == 0)
  	kind = 6;   
  else if ((strncasecmp(argv[1],"temps",3) == 0) || (strncasecmp(argv[1],"time",3) == 0))
  	kind = 10;
  else if (strncasecmp(argv[1],"galchen",3) == 0)
  	kind = 21;
  else if (strncasecmp(argv[1],"kind_",5) == 0) 
    kind = atoi(argv[1]+5);
  else {                            
  	printf("Directive_ip123 error: unknown kind =%s\n",argv[1]);
    kind = -1;
    }

  if (argc > 2) {
    sscanf(argv[2],"[%d]",&nelm);
    offset=3;
    }
  else {
    nelm = 1;
    offset=2;
    }

  printf("ip123 set_nb=%d, des_exc=%d, nelm=%d, kind=%d\n\n",set_nb,des_exc,nelm,kind);

  range = 0;
  n = 0;
  if ((kind == -1) || (kind == -2)) {
    for (i=0; i<nelm; i++) {
      if ((i == 1) && (strncasecmp(argv[i+offset],"to",2) == 0))
        range = 1;
      else {
        list_entier[n] = atoi(argv[i+offset]);
        printf("list_entier[%d]=%d\n",n,list_entier[n]);
        n++;
        }
      }
    if (range) nelm = -2;
    func(set_nb,des_exc,list_entier,nelm,kind);
  }
  else {
    for (i=0; i<nelm; i++) {
      if ((i == 1) && (strncasecmp(argv[i+offset],"to",2) == 0))
        range = 1;
      else {
        sscanf(argv[i+offset],"%f",list_reel+n);
        printf("list_reel[%d]=%f\n",n,list_reel[n]);
        n++;
        }
      }
    if (range) nelm = -2;  
    func(set_nb,des_exc,list_reel,nelm,kind);
  }

  return(0);
}

/*****************************************************************************
 *                      D I R E C T I V E S _ D A T E S                      *
 *                                                                           *
 *Objet                                                                      *
 *   Traiter les directives desires pour la date (en format datestamp)       *
 *                                                                           *
 *****************************************************************************/
int Directive_dates(int argc , char **argv,char cmd_strt,  void *func(int set_nb, int des_exc, int *date_list, int nelm, float delta), char *Private_Data_2)
{
  int i, set_nb, des_exc, nelm, offset, n, range;
  float delta;
  int list_entier[100];
  int ier;

  for (i=0;i<=argc;i++)
    printf("argument # %d-->%s<-- \n",i,argv[i]);
	
  set_nb = bundle_nb;
  des_exc = desire_exclure;
  delta = 0.0;

  if (argc > 1) {
    sscanf(argv[1],"[%d]",&nelm);
    offset=2;
  }
  else {
    nelm = 1;
    offset=1;
  }

  n = 0;
  range = 0;
  for (i=0; i<nelm; i++) {
    if ((i == 1) && (strncasecmp(argv[i+offset],"to",2) == 0))
      range = 1;
    else {
      if ((i == 3) && (range)) {
        sscanf(argv[i+offset],"%f",&delta);
        printf("delta=%f\n",delta);
        break;
        }
      else { 
        list_entier[n] = atoi(argv[i+offset]);
        printf("list_entier[%d]=%d\n",n,list_entier[n]);
        n++;
        }
      }
    }

  if (range) nelm = -2;
  /* ier = func(set_nb,des_exc,list_entier,nelm,delta); */
  Xc_Select_date(set_nb,des_exc,list_entier,nelm);
  return(0);
}

/*****************************************************************************
 *                      D I R E C T I V E S _ D A T E V                      *
 *                                                                           *
 *Objet                                                                      *
 *   Traiter les directives desires pour la date (en format visuel)          *
 *                                                                           *
 *****************************************************************************/
int Directive_datev(int argc , char **argv,char cmd_strt,  void *func(), char *Private_Data_2)
{
  int i, j, set_nb, des_exc, nelm, offset;
  int yyyymmdd, hhmmss, stamp;
  int mode=3, range;
  float delta;
  int list_entier[100];

  for (i=0;i<=argc;i++)
    printf("argument # %d-->%s<-- \n",i,argv[i]);

/*  set_nb = atoi(argv[1]); */
/*  des_exc = atoi(argv[2]); */
  set_nb = bundle_nb;
  des_exc = desire_exclure;
  delta = 0.0;

  if (argc > 1) {
    sscanf(argv[1],"[%d]",&nelm);
    if (((nelm % 2) != 0) && (nelm != 5)) {
    	printf("Directive_datev error: visual date format list must contain an even number of date elements\n");
      return(-1);
      }
    offset=2;
  }
  else {
  	printf("Directive_datev error: visual date format list must contain an even number of date elements\n");
    return(-1);  
  }

  range = 0;
  j = 0;
  for (i=0; i<nelm; i++) {
    if ((i == 2) && (strncasecmp(argv[i+offset],"to",2) == 0))
      range = 1;
    else {  
      if ((i == 5) && (range)) {
        sscanf(argv[i+offset],"%f",&delta);
        printf("delta=%f\n",delta);
        break;
        }
      else { 
        yyyymmdd = atoi(argv[i+offset]);
        hhmmss =   atoi(argv[i+1+offset]);
        f77name(newdate)(&stamp,&yyyymmdd,&hhmmss,&mode);
        list_entier[j] = stamp;
        printf("list_entier[%d]=%d\n",j,list_entier[j]);
        j++;
        i++;
        }
      }
    }

  if (range)
    nelm = -2;
  else
    nelm = j;  
  /* func(set_nb,des_exc,list_entier,nelm,delta); */
  Xc_Select_date(set_nb,des_exc,list_entier,nelm);
  return(0);
}

/*****************************************************************************
 *                      D I R E C T I V E S _ C H A R V A R                  *
 *                                                                           *
 *Objet                                                                      *
 *   Traiter les directives desires pour nomvar, typvar et etiquette         *
 *                                                                           *
 *****************************************************************************/
int Directive_charvar(int argc , char **argv,char cmd_strt,  void *func(), char *Private_Data_2)
{
  int i, set_nb, des_exc, nelm, offset;
  char **string_array;

/*  for (i=0;i<=argc;i++)
    printf("argument # %d-->%s<-- \n",i,argv[i]); */

	set_nb = bundle_nb;
  des_exc = desire_exclure;

  if (argc > 1) {
    sscanf(argv[1],"[%d]",&nelm);
    offset=2;
  }
  else {
    nelm = 1;
    offset=1;
  }

  printf("\ncharvar set_nb=%d, des_exc=%d, nelm=%d\n",set_nb,des_exc,nelm);

  string_array = allocate_string_array(nelm);
  for (i=0; i<nelm; i++) {
    if ((*argv[i+offset] == '\'') || (*argv[i+offset] == '"')) {
    	*(argv[i+offset]+strlen(argv[i+offset])-1) = '\0';
      argv[i+offset]++;
    }
    string_array[i] = malloc(12 * sizeof(char));
    sscanf(argv[i+offset],"%s",string_array[i]);
    printf("string_array[%d]=%s\n",i,string_array[i]);
  }

  /* printf("Debug+ Directive_charvar appel la fonction avec set_nb=%d nelm=%d\n",set_nb,nelm); */
  func(set_nb,des_exc,string_array,nelm);
  free_string_array(string_array);
  return(0);
}

/*****************************************************************************
 *                      D I R E C T I V E S _ D E S I R E                    *
 *                                                                           *
 *Objet                                                                      *
 *   Change la variable globale pour le mode desire                          *
 *                                                                           *
 *****************************************************************************/
int Directive_desire(int argc , char **argv,char cmd_strt,  void *func(), char *Private_Data_2)
{
	func();
/*CHC/NRC*/
	return 0;
}

/*****************************************************************************
 *                      D I R E C T I V E S _ E X C L U R E                  *
 *                                                                           *
 *Objet                                                                      *
 *   Change la variable globale pour le mode exclure                         *
 *                                                                           *
 *****************************************************************************/
int Directive_exclure(int argc , char **argv,char cmd_strt,  void *func(), char *Private_Data_2)
{
	func();
/*CHC/NRC*/
	return 0;
}

/*****************************************************************************
 *                C _ R E Q U E T E S _ R E A D _ F I L E                    *
 *                                                                           *
 *Objet                                                                      *
 *   Definir un ensemble de requetes a partir d'un fichier                   *
 *                                                                           *
 *****************************************************************************/
int C_requetes_read_file(char *requetes_file)
{
  FILE *fp;
  int *pv2 = 0;

fprintf(stderr,"C _ R E Q U E T E S _ R E A D _ F I L E : FEATURE NOT ACTIVE\n");
return 0;

  fp = fopen(requetes_file,"r");
  if (fp == (FILE *) NULL) {
    fprintf(stderr,"C_requetes_read_file error opening file %s\n",requetes_file);
    fprintf(stderr,"filter file ignored\n");
    return(-1);
    }
  else
    fprintf(stderr,"warning: FSTD software uses filter file %s for scans\n",requetes_file);

  rpn_c_callback("ip1",Directive_ip123,"",Xc_Select_ip1,pv2);
  rpn_c_callback("ip2",Directive_ip123,"",Xc_Select_ip2,pv2);
  rpn_c_callback("ip3",Directive_ip123,"",Xc_Select_ip3,pv2);
  rpn_c_callback("datestamp",Directive_dates,"",Xc_Select_date,pv2);
  rpn_c_callback("datevisual",Directive_datev,"",Xc_Select_date,pv2);
  rpn_c_callback("nomvar",Directive_charvar,"",Xc_Select_nomvar,pv2);
  rpn_c_callback("typvar",Directive_charvar,"",Xc_Select_typvar,pv2);
  rpn_c_callback("etiquette",Directive_charvar,"",Xc_Select_etiquette,pv2);
  rpn_c_callback("desire",Directive_desire,"",C_filtre_desire,pv2);
  rpn_c_callback("exclure",Directive_exclure,"",C_filtre_exclure,pv2);
  process_c_callback(requetes_file);

	return 0; /*CHC/NRC*/
}

#endif
/*****************************************   end of deprecated code *****************************************/
#if defined (TEST)
//void f77name(c_main)(int *handle)
//c_main(int argc, char **argv)
c_main(int argc, char **argv)
{

  union {
    int i;
    float f;
  }i_or_f;
  int i,j;
  float heures;
  int ip1s_i[] = { 400,500,600,750 };
  float ip1s_r[] = { .3840, .4440, 0.6110 };
  float ip1s_r_range[] = { .3280, 0.8000 };
  int ip1s_range[] = {400, READLX_RANGE, 750};
  int ip1s_range2[] = {READLX_RANGE, 750};
  int ip1s_range3[] = { 750, READLX_RANGE};
  int ip2s[] = {0, 12, 24};
  int ip3s[] = {0, 80, 90};
  int dates[] = {313290800, 313301600};
  int dates_range[] = {313280000, READLX_RANGE, 313290800};
  int dates_range1[] = {313280000, READLX_RANGE, 313290800, READLX_DELTA, 12};
  int dates_range2[] = {313280000, READLX_RANGE};
  int dates_range3[] = {READLX_RANGE, 313290800};
  int dates_range4[] = {313280000, READLX_RANGE, 313290800, READLX_DELTA, 12};

  char *testeti[] = { "Etiquette #1", "R2428V4N", "Etiquet #3" };
  char *testnom[2] = { "TT", "GZ" };
  char *testtyp[4] = { "P", "V1", "V2", "V3" };
  int narg = 2;

  i_or_f.f = 5399.6;
  dates_range1[4] = i_or_f.i;

  dbprint(stddebug,"Debug debut \n");
  dbprint(stddebug,"Debug testeti=%s %s %s \n",testeti[0],testeti[1],testeti[2]);
  dbprint(stddebug,"Debug testeti=%s %s %s \n",testeti[0],testeti[1],testeti[2]);
  dbprint(stddebug,"Debug testnom=%s %s \n",testnom[0],testnom[1]);
  dbprint(stddebug,"Debug testtyp=%s %s %s %s \n",testtyp[0],testtyp[1],testtyp[2],testtyp[3]);
  RequetesInit();
fprintf(stderr,"DEBUG: %d\n",1);
  i = Xc_Select_etiquette(5,1,testeti,3);
  i = Xc_Select_nomvar(5,1,testnom,2);
  i = Xc_Select_typvar(5,1,testtyp,4);
  i = Xc_Select_ip1(5,1,ip1s_range,3);
//  i = Xc_Select_ip1(1,1,ip1s_i,4);
  i = Xc_Select_ip2(5,1,ip2s,3);
  i = Xc_Select_date(5,1,dates_range4,5);

  i = Xc_Select_etiquette(2,1,testeti,3);
  i = Xc_Select_nomvar(2,1,testnom,2);
  i = Xc_Select_typvar(2,1,testtyp,4);
  i = Xc_Select_ip1(2,1,ip1s_range,3);
//  i = Xc_Select_ip1(2,1,ip1s_i,4);
  i = Xc_Select_ip2(2,1,ip2s,3);
  i = Xc_Select_date(2,1,dates_range1,5);

  i = Xc_Select_ip1(3,1,ip1s_range2,2);
  i = Xc_Select_ip3(3,1,ip3s,3);
  i = Xc_Select_ip2(3,1,ip1s_range3,2);
  i = Xc_Select_date(3,1,dates_range2,2);

  i = Xc_Select_ip1(4,-1,ip1s_range3,2);
  i = Xc_Select_ip2(4,-1,ip1s_range2,2);
  i = Xc_Select_date(4,-1,dates_range3,2);
  i = C_select_groupset(2,5);
  if(argc > 1)
    WriteRequestTable(atoi(argv[1]),NULL);
  else
    WriteRequestTable(0,NULL);

//i=C_fstmatch_parm(handle,  datevalid, ni, nj, nk,       ip1,       ip2,       ip3,typvar,nomvar,      etiket,grtyp,ig1,ig2,ig3,ig4)
  i=C_fstmatch_parm(-1    ,  313280000, 50, 51, 52,       500,        12,         0,   "P",  "TT",  "R2428V4N","X"  ,  0,  0,  0,  0);
  fprintf(stderr,"match result = %d \n",i);
  i=C_fstmatch_parm(-1    ,  313280000, 50, 51, 52,  41394464,        12,         0,   "P",  "TT",  "R2428V4N","X"  ,  0,  0,  0,  0);
  fprintf(stderr,"match result = %d \n",i);
  i=C_fstmatch_parm(-1    ,  313290800, 50, 51, 52,       750,       750,     77777,   "P",  "TT",  "R2428V4N","X"  ,  0,  0,  0,  0);
  fprintf(stderr,"match result = %d \n",i);
  if(argc>narg){
    WriteRequestTable(0,argv[narg]);
    RequetesInit();
    ReadRequestTable(argv[narg]);
    fprintf(stdout,"=========== Reading Back table===========\n");
    WriteRequestTable(1,NULL);
  }
  narg++ ;
  fprintf(stdout,"\n=========================================\n\n");

  RequetesInit();
  i = Xc_Select_ip1(1,1,ip1s_i,4);
//  i = Xc_Select_suppl(1,1,-1, -1, -1, -1, -1, -1, -1, 'Z');
  i = Xc_Select_date(1,1,dates_range,3);
  i = C_select_groupset(0,1);
  if(argc > 1)
    WriteRequestTable(atoi(argv[1]),NULL);
  else
    WriteRequestTable(0,NULL);
  if(argc>narg){
    WriteRequestTable(0,argv[narg]);
    RequetesInit();
    ReadRequestTable(argv[narg]);
    fprintf(stdout,"=========== Reading Back table===========\n");
    WriteRequestTable(1,NULL);
  }
  narg++ ;
//  DumpRequestTable(0,NULL);
//  i = Xc_Select_ip1(0,1,ip1s_r,3);
/*  i = Xc_Select_ip1(1,1,ip1s_range,-2);*/
//  i = Xc_Select_ip1(1,1,ip1s_r_range,-2);
//  heures=12.0;
//  i = Xc_Select_date(1,1,dates_range,-2);
/*  heures=0.0;
  i = Xc_Select_date(1,1,dates_range,-2,heures);*/
 /* j = C_fstmatch_req(1,*handle);*/
//  j = C_fstmatch_req(*handle);

}
#endif

