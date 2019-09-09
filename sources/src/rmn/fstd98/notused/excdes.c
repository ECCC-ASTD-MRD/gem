/*****************************************************************************
 *                        E X C D E S . C                                    *
 *                                                                           *
 *Auteur                                                                     *
 *   Mario Lépine  -  Mai 2004                                               *
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
#include "qstdir.h"
#include "ip_kind.h"		/*CHC/NRC*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define MAX_Nlist 40
#define MAX_requetes 20
#define Abs(x) (((x) < 0) ? -(x) : x)
/* #if defined (DEBUG) */
#define dbprint fprintf
/* #else
   #define dbprint ;
   #endif */

char **fill_string_array(char **string_array, char *farray, int nc, int ns, int rmblanks);
char **allocate_string_array(int ns);

enum cquoica {entier, reel, deb_fin_entier, deb_fin_reel, deb_fin_entier_delta,
              deb_fin_reel_delta} parametre;
enum onveut {desire, exclure} call_mode;

static int first_R = 0;
static int last_R = MAX_requetes;
static FILE *stddebug;
static int bundle_nb = -1;
static int desire_exclure = 1;

typedef struct {
  int in_use;
  enum cquoica arg_type;
  int ip_kind;
  int nelm;
  float delta;
  union {
    int tab_elem[MAX_Nlist];
    float tab_rval[MAX_Nlist];
  } data;
} DE_int_float;

typedef struct {
  int in_use;
  int nelm;
  char *sdata[MAX_Nlist];
} DE_char;

typedef struct {
  int in_use;
  enum onveut exdes;
  DE_char etiquettes;
  DE_char nomvars;
  DE_char typvars;
  DE_int_float dates;
  DE_int_float ip1s;
  DE_int_float ip2s;
  DE_int_float ip3s;
} Desire_Exclure;

Desire_Exclure Requetes[MAX_requetes];

/*****************************************************************************
 *                    X C _ S E L E C T _ I P 1                              *
 *                                                                           *
 *Objet                                                                      *
 *   Definir la liste des IP1 desires                                        *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  set_nb  numero associe a un groupe d'elements desire/exclure         *
 *  IN  des_exc 0=exclure 1=desire                                           *
 *  IN  iplist  liste de 1 ou plusieurs IP1 a rechercher ou exclure          *
 *  IN  nelm    nombre d'elements de la liste                                *
 *              -2=interval [debut, fin]                                     *
 *                           debut=-1, indique valeur min du fichier         *
 *                           fin=-1, indique valeur max du fichier           *
 *  IN  kind    code kind de convip_plus,                                    *
 *              -1=IP1 entier deja encode                                    *
 *              -2=IP1 entier a garder tel quel (pas de convip_plus)         *
 *                                                                           *
 *****************************************************************************/
int xc_select_ip1(int set_nb, int des_exc, void *iplist, int nelm, int kind)
{
  int i, range=0;
  int *ip_entier=(int *)iplist;
  float *ip_reel=(float *)iplist;
  float p1;
  int mode=-1,flag=0;
  char string[30];

  if (set_nb > MAX_requetes) {
    fprintf(stderr,"error: c_select_ip1 set_nb=%d > MAX_requetes=%d\n",set_nb,MAX_requetes);
    return(-1);
  }

  if (nelm > MAX_Nlist) {
    fprintf(stderr,"error: c_select_ip1 nelm=%d > limits=%d\n",nelm,MAX_Nlist);
    return(-2);
  }

  if (nelm < 0) {
    if (nelm == -2) {
      range = 1;
      nelm = 2;
    }
    else {
      fprintf(stderr,"error: c_select_ip1 nelm invalid = %d\n",nelm);
      return(-3);
    }
  }
  if ((Requetes[set_nb].in_use) && (((des_exc == 1) ? desire : exclure) != Requetes[set_nb].exdes)) {
    fprintf(stderr,"c_select_ip1 error: des_exc value differs from previous call for set number=%d\n",set_nb);
    return(-4);
  }
  Requetes[set_nb].in_use = 1;
  Requetes[set_nb].ip1s.in_use = 1;
  Requetes[set_nb].exdes = (des_exc == 1) ? desire : exclure;
  Requetes[set_nb].ip1s.nelm = nelm;
  Requetes[set_nb].ip1s.ip_kind = kind;
  if (kind == -1) {
    Requetes[set_nb].ip1s.arg_type = (range ==1) ? deb_fin_reel : reel;
    for (i=0; i<nelm; i++) {
      if (ip_entier[i] == -1)
        p1 = -1.0;
      else
        f77name(convip_plus)(&ip_entier[i],&p1,&Requetes[set_nb].ip1s.ip_kind,&mode,string,&flag);
      Requetes[set_nb].ip1s.data.tab_rval[i] = p1;
      dbprint(stddebug,"Debug ip_entier=%d Requetes[%d].ip1s.data.tab_rval[%d] = %f\n",ip_entier[i],set_nb,i,
              Requetes[set_nb].ip1s.data.tab_rval[i]);
    }
  }
  else if (kind == -2) {
    Requetes[set_nb].ip1s.arg_type = (range ==1) ? deb_fin_reel : reel;
    for (i=0; i<nelm; i++) {
      Requetes[set_nb].ip1s.data.tab_rval[i] = (float) ip_entier[i];
      dbprint(stddebug,"Debug Requetes[%d].ip1s.data.tab_rval[%d] = %f\n",set_nb,i,
              Requetes[set_nb].ip1s.data.tab_rval[i]);
    }
  }
  else {
    Requetes[set_nb].ip1s.arg_type = (range ==1) ? deb_fin_reel : reel;
    for (i=0; i<nelm; i++) {
      Requetes[set_nb].ip1s.data.tab_rval[i] = ip_reel[i];
      dbprint(stddebug,"Debug Requetes[%d].ip1s.data.tab_rval[%d] = %f\n",set_nb,i,
              Requetes[set_nb].ip1s.data.tab_rval[i]);
    }
  }
  return(0);
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
 *  IN  nelm    nombre d'elements de la liste                                *
 *              -2=interval [debut, fin]                                     *
 *                           debut=-1, indique valeur min du fichier         *
 *                           fin=-1, indique valeur max du fichier           *
 *  IN  kind    code kind de convip_plus, (ex: kind=10, temp en heure)       *
 *              -1=IP2 entier deja encode (valeur classique standard)        *
 *                                                                           *
 *****************************************************************************/
int xc_select_ip2(int set_nb, int des_exc, void *iplist, int nelm, int kind)
{
  int i, range=0;
  int *ip_entier=(int *)iplist;
  float *ip_reel=(float *)iplist;

  if (set_nb > MAX_requetes) {
    fprintf(stderr,"error: c_select_ip2 set_nb=%d > MAX_requetes=%d\n",set_nb,MAX_requetes);
    return(-1);
  }

  if (nelm > MAX_Nlist) {
    fprintf(stderr,"error: c_select_ip2 nelm=%d > limits=%d\n",nelm,MAX_Nlist);
    return(-2);
  }

  if (nelm < 0) {
    if (nelm == -2) {
      range = 1;
      nelm = 2;
    }
    else {
      fprintf(stderr,"error: c_select_ip2 nelm invalid = %d\n",nelm);
      return(-3);
    }
  }
  if ((Requetes[set_nb].in_use) && (((des_exc == 1) ? desire : exclure) != Requetes[set_nb].exdes)) {
    fprintf(stderr,"c_select_ip2 error: des_exc value differs from previous call for set number=%d\n",set_nb);
    return(-4);
  }
  Requetes[set_nb].in_use = 1;
  Requetes[set_nb].ip2s.in_use = 1;
  Requetes[set_nb].exdes = (des_exc == 1) ? desire : exclure;
  Requetes[set_nb].ip2s.nelm = nelm;
  Requetes[set_nb].ip2s.ip_kind = kind;
  if (kind == -1) {
    Requetes[set_nb].ip2s.arg_type = (range ==1) ? deb_fin_entier : entier;
    for (i=0; i<nelm; i++) {
      Requetes[set_nb].ip2s.data.tab_elem[i] = ip_entier[i];
      dbprint(stddebug,"Debug Requetes[%d].ip2s.data.tab_elem[%d] = %d\n",set_nb,i,
              Requetes[set_nb].ip2s.data.tab_elem[i]);
    }
  }
  else {
    Requetes[set_nb].ip2s.arg_type = (range ==1) ? deb_fin_reel : reel;
    for (i=0; i<nelm; i++) {
      Requetes[set_nb].ip2s.data.tab_rval[i] = ip_reel[i];
      dbprint(stddebug,"Debug Requetes[%d].ip2s.data.tab_rval[%d] = %f\n",set_nb,i,
              Requetes[set_nb].ip2s.data.tab_rval[i]);
    }
  }
  return(0);
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
 *  IN  nelm    nombre d'elements de la liste                                *
 *              -2=interval [debut, fin] (pour expansion future)             *
 *                           debut=-1, indique valeur min du fichier         *
 *                           fin=-1, indique valeur max du fichier           *
 *  IN  kind    code kind de convip_plus, (pour expansion future, ignore)    *
 *              -1=IP3 entier (valeur classique standard)                    *
 *                                                                           *
 *****************************************************************************/
int xc_select_ip3(int set_nb, int des_exc, void *iplist, int nelm, int kind)
{
  int i;
  int *ip_entier=(int *)iplist;
  float *ip_reel=(float *)iplist;

  if (set_nb > MAX_requetes) {
    fprintf(stderr,"error: c_select_ip3 set_nb=%d > MAX_requetes=%d\n",set_nb,MAX_requetes);
    return(-1);
  }

  if (nelm > MAX_Nlist) {
    fprintf(stderr,"error: c_select_ip3 nelm=%d > limits=%d\n",nelm,MAX_Nlist);
    return(-2);
  }

  if (nelm < 0) {
    if (nelm == -2) {
      fprintf(stderr,"error: c_select_ip3, interval non supporte pour le moment\n");
      return(-3);
    }
    else {
      fprintf(stderr,"error: c_select_ip3 nelm invalid = %d\n",nelm);
      return(-4);
    }
  }
  if ((Requetes[set_nb].in_use) && (((des_exc == 1) ? desire : exclure) != Requetes[set_nb].exdes)) {
    fprintf(stderr,"c_select_ip3 error: des_exc value differs from previous call for set number=%d\n",set_nb);
    return(-4);
  }
  Requetes[set_nb].in_use = 1;
  Requetes[set_nb].ip3s.in_use = 1;
  Requetes[set_nb].exdes = (des_exc == 1) ? desire : exclure;
  Requetes[set_nb].ip3s.nelm = nelm;
  Requetes[set_nb].ip3s.ip_kind = kind;
  Requetes[set_nb].ip3s.arg_type = entier;
  for (i=0; i<nelm; i++) {
    Requetes[set_nb].ip3s.data.tab_elem[i] = ip_entier[i];
    dbprint(stddebug,"Debug Requetes[%d].ip3s.data.tab_elem[%d] = %d\n",set_nb,i,
           Requetes[set_nb].ip3s.data.tab_elem[i]);
  }
  return(0);
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
 *  IN  nelm       nombre d'elements de la liste                             *
 *                 -2=interval [debut, fin] ou [debut, fin] avec delta       *
 *                              debut=-1, indique valeur min du fichier      *
 *                              fin=-1, indique valeur max du fichier        *
 *  IN  delta      delta en heure ou fraction d'heure                        *
 *                 0=aucun delta desire                                      *
 *                                                                           *
 *****************************************************************************/
int xc_select_date(int set_nb, int des_exc, int *date_list, int nelm, float delta)
{
  int i, range=0;
  int *pt_delta = (int *) &delta;
  float zero = 0.0;
  int *pt_zero = (int *) &zero;
  
  if (set_nb > MAX_requetes) {
    fprintf(stderr,"error: c_select_date set_nb=%d > MAX_requetes=%d\n",set_nb,MAX_requetes);
    return(-1);
  }

  if (nelm > MAX_Nlist) {
    fprintf(stderr,"error: c_select_date nelm=%d > limits=%d\n",nelm,MAX_Nlist);
    return(-2);
  }

  if (nelm < 0) {
    if (nelm == -2) {
      range = 1;
      nelm = 2;
    }
    else {
      fprintf(stderr,"error: c_select_date nelm invalid = %d\n",nelm);
      return(-3);
    }
  }
  if ((Requetes[set_nb].in_use) && (((des_exc == 1) ? desire : exclure) != Requetes[set_nb].exdes)) {
    fprintf(stderr,"c_select_date error: des_exc value differs from previous call for set number=%d\n",set_nb);
    return(-4);
  }
  Requetes[set_nb].in_use = 1;
  Requetes[set_nb].dates.in_use = 1;
  Requetes[set_nb].exdes = (des_exc == 1) ? desire : exclure;
  Requetes[set_nb].dates.nelm = nelm;
  Requetes[set_nb].dates.delta= delta;
  if (delta != 0.) {
    Requetes[set_nb].dates.arg_type = deb_fin_entier_delta;
    dbprint(stddebug,"Debug Requetes[%i].dates.arg_type=%d delta=%f\n",
           set_nb,Requetes[set_nb].dates.arg_type,Requetes[set_nb].dates.delta);
  }
  else
    Requetes[set_nb].dates.arg_type = (range ==1) ? deb_fin_entier : entier;
  for (i=0; i<nelm; i++) {
    Requetes[set_nb].dates.data.tab_elem[i] = date_list[i];
    dbprint(stddebug,"Debug Requetes[%d].dates.data.tab_elem[%d] = %d\n",set_nb,i,
           Requetes[set_nb].dates.data.tab_elem[i]);
  }
  return(0);
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
int xc_select_etiquette(int set_nb, int des_exc, char *etiq_list[], int nelm)
{
  int i;

  if (set_nb > MAX_requetes) {
    fprintf(stderr,"error: c_select_etiquette set_nb=%d > MAX_requetes=%d\n",set_nb,MAX_requetes);
    return(-1);
  }

  if (nelm > MAX_Nlist) {
    fprintf(stderr,"error: c_select_etiquette nelm=%d > limits=%d\n",nelm,MAX_Nlist);
    return(-2);
  }

  if (nelm <= 0) {
    fprintf(stderr,"error: c_select_etiquette nelm invalid = %d\n",nelm);
    return(-3);
  }
  if ((Requetes[set_nb].in_use) && (((des_exc == 1) ? desire : exclure) != Requetes[set_nb].exdes)) {
    fprintf(stderr,"c_select_etiquette error: des_exc value differs from previous call for set number=%d\n",set_nb);
    return(-4);
  }
  Requetes[set_nb].in_use = 1;
  Requetes[set_nb].etiquettes.in_use = 1;
  Requetes[set_nb].exdes = (des_exc == 1) ? desire : exclure;
  Requetes[set_nb].etiquettes.nelm = nelm;
  for (i=0; i<nelm; i++) {
    strncpy(Requetes[set_nb].etiquettes.sdata[i],etiq_list[i],13);
    dbprint(stddebug,"Debug Requetes[%i].etiquettes.sdata[%i]=%s\n",
           set_nb,i,Requetes[set_nb].etiquettes.sdata[i]);
  }
  return(0);
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
int xc_select_nomvar(int set_nb, int des_exc, char *nomv_list[], int nelm)
{
  int i;
  
  if (set_nb > MAX_requetes) {
    fprintf(stderr,"error: c_select_nomvar set_nb=%d > MAX_requetes=%d\n",set_nb,MAX_requetes);
    return(-1);
  }

  if (nelm > MAX_Nlist) {
    fprintf(stderr,"error: c_select_nomvar nelm=%d > limits=%d\n",nelm,MAX_Nlist);
    return(-2);
  }

  if (nelm <= 0) {
    fprintf(stderr,"error: c_select_nomvar nelm invalid = %d\n",nelm);
    return(-3);
  }
  if ((Requetes[set_nb].in_use) && (((des_exc == 1) ? desire : exclure) != Requetes[set_nb].exdes)) {
    fprintf(stderr,"c_select_nomvar error: des_exc value differs from previous call for set number=%d\n",set_nb);
    return(-4);
  }
  Requetes[set_nb].in_use = 1;
  Requetes[set_nb].nomvars.in_use = 1;
  Requetes[set_nb].exdes = (des_exc == 1) ? desire : exclure;
  Requetes[set_nb].nomvars.nelm = nelm;
  for (i=0; i<nelm; i++) {
    strncpy(Requetes[set_nb].nomvars.sdata[i],nomv_list[i],5);
    dbprint(stddebug,"Debug Requetes[%i].nomvars.sdata[%i]=%s\n",
           set_nb,i,Requetes[set_nb].nomvars.sdata[i]);
  }
  return(0);
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
int xc_select_typvar(int set_nb, int des_exc, char *typv_list[], int nelm)
{
  int i;

  if (set_nb > MAX_requetes) {
    fprintf(stderr,"error: c_select_typvar set_nb=%d > MAX_requetes=%d\n",set_nb,MAX_requetes);
    return(-1);
  }

  if (nelm > MAX_Nlist) {
    fprintf(stderr,"error: c_select_typvar nelm=%d > limits=%d\n",nelm,MAX_Nlist);
    return(-2);
  }

  if (nelm <= 0) {
    fprintf(stderr,"error: c_select_typvar nelm invalid = %d\n",nelm);
    return(-3);
  }
  if ((Requetes[set_nb].in_use) && (((des_exc == 1) ? desire : exclure) != Requetes[set_nb].exdes)) {
    fprintf(stderr,"c_select_typvar error: des_exc value differs from previous call for set number=%d\n",set_nb);
    return(-4);
  }
  Requetes[set_nb].in_use = 1;
  Requetes[set_nb].typvars.in_use = 1;
  Requetes[set_nb].exdes = (des_exc == 1) ? desire : exclure;
  Requetes[set_nb].typvars.nelm = nelm;
  for (i=0; i<nelm; i++) {
    strncpy(Requetes[set_nb].typvars.sdata[i],typv_list[i],3);
    dbprint(stddebug,"Debug Requetes[%i].typvars.sdata[%i]=%s\n",
           set_nb,i,Requetes[set_nb].typvars.sdata[i]);
  }
  return(0);
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
int c_select_groupset(int first_set_nb, int last_set_nb)
{
  if ((first_set_nb > MAX_requetes) || (last_set_nb > MAX_requetes) || (first_set_nb > last_set_nb)) {
    fprintf(stderr,"error: c_select_groupset first_set_nb=%d, last_set_nb=%d, MAX_requetes=%d\n",
            first_set_nb,last_set_nb,MAX_requetes);
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
 *                                                                           *
 *****************************************************************************/
int c_filtre_desire()
{
  
  bundle_nb++;
  desire_exclure = 1;
  if (bundle_nb > MAX_requetes) {
    fprintf(stderr,"error: c_filtre_desire nb=%d > MAX desire/exclure =%d\n",bundle_nb,MAX_requetes);
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
 *                                                                           *
 *****************************************************************************/
int c_filtre_exclure()
{
  
  bundle_nb++;
  desire_exclure = 0;
  if (bundle_nb > MAX_requetes) {
    fprintf(stderr,"error: c_filtre_exclure nb=%d > MAX desire/exclure =%d\n",bundle_nb,MAX_requetes);
    return(-1);
  }
  printf("exclure bundle_nb = %d, desire_exclure = %d\n",bundle_nb,desire_exclure);
return 0; /*CHC/NRC*/
}

/*****************************************************************************
 *                      C _ F S T _ M A T C H _ R E Q                        *
 *                                                                           *
 *Objet                                                                      *
 *   Verifier si l'enregistrement courant correspond aux requetes demandees  *
 *   Retourne 1 si l'enregistrement correspond, sinon 0                      *
 *   0 est aussi retourne si l'enregistrement correspond a un a exclure      *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  set_nb     numero associe a un groupe d'elements desire/exclure      *
 *  IN  handle     handle de l'enregistrement courant                        *
 *                                                                           *
 *****************************************************************************/
int c_fst_match_req(int handle)
{
  int ier, i, j, set_nb, last_in_use;
  int ni, nj, nk, date, dateo, ip1, ip2, ip3, ig1, ig2, ig3, ig4;
  int nbits, swa, ubc, lng, dltf, datevalid, xtra2, xtra3, deet, npas, datyp;
  char etiket[13]={' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ','\0'};
  char typvar[3]={' ',' ','\0'};
  char nomvar[5]={' ',' ',' ',' ','\0'};
  char grtyp[2]={' ','\0'};
  int mode=-1,flag=0;
  int ip_kind, amatch=0;
  char string[30], *stripblanks;
  int debut, fin;
  float err_tolerance = 1.0/pow(2.0,17);
  float p1,r_debut,r_fin;
  double diff_deb, diff_fin, delta8, date8, remainder;

  if (! Requetes[first_R].in_use) return(1);        /* aucune requete desire ou exclure */

  ier = c_fstprm(handle,&dateo,&deet,&npas,&ni,&nj,&nk,&nbits,&datyp,&ip1,
                     &ip2,&ip3,typvar,nomvar,etiket,grtyp,&ig1,&ig2,
                     &ig3,&ig4,&swa,&lng,&dltf,&ubc,&datevalid,&xtra2,&xtra3);
/*  dbprint(stddebug,"Debug c_fst_match_req fstprm date=%d ip1=%d ip2=%d ip3=%d nomvar-->%s<-- typvar-->%s<-- etiket-->%s<--\n",
         date,ip1,ip2,ip3,nomvar,typvar,etiket);*/
  if (ier < 0) return(0);
  date = datevalid;
/*  
  stripblanks = strtok(nomvar," ");
  stripblanks = strtok(typvar," ");
  stripblanks = strtok(etiket," ");
*/
  for (set_nb=first_R; (set_nb <= last_R) && Requetes[set_nb].in_use; set_nb++) {

    amatch = 0;

    Etiquettes:
      if (Requetes[set_nb].etiquettes.in_use) {
      dbprint(stddebug,"Debug c_fst_match_req verifie etiquettes du fichier=%s set_nb=%d\n",etiket,set_nb);
        if (Requetes[set_nb].exdes == desire) {
          for (i=0; i < Requetes[set_nb].etiquettes.nelm; i++)
            if (strncmp(Requetes[set_nb].etiquettes.sdata[i],etiket,Min(12,strlen(Requetes[set_nb].etiquettes.sdata[i]))) == 0) {
              amatch = 1;
              dbprint(stddebug,"Debug c_fst_match_req match desire\n");
              goto Nomvars; /* satisfait la requete */
            }
          continue;  /* rien trouve qui satisfait la requete pour etiquettes */
        }
        else {
          for (i=0; i < Requetes[set_nb].etiquettes.nelm; i++)
            if (strncmp(Requetes[set_nb].etiquettes.sdata[i],etiket,Min(12,strlen(Requetes[set_nb].etiquettes.sdata[i])))== 0) {
              amatch = -1;  /* enregistrement a exclure */
              dbprint(stddebug,"Debug c_fst_match_req match exclure\n");
              goto Nomvars;
            }
          continue;   /* rien trouve a exclure */
        }
      }

    Nomvars:
      if (Requetes[set_nb].nomvars.in_use) {
      dbprint(stddebug,"Debug c_fst_match_req verifie nomvars du fichier=%s set_nb=%d\n",nomvar,set_nb);
        if (Requetes[set_nb].exdes == desire) {
          for (i=0; i < Requetes[set_nb].nomvars.nelm; i++)
            if (strncmp(Requetes[set_nb].nomvars.sdata[i],nomvar,Min(4,strlen(Requetes[set_nb].nomvars.sdata[i]))) == 0) {
              amatch = 1;
              dbprint(stddebug,"Debug c_fst_match_req match desire\n");
              goto Typvars; /* satisfait la requete */
            }
          continue;  /* rien trouve qui satisfait la requete pour nomvars */
        }
        else {
          for (i=0; i < Requetes[set_nb].nomvars.nelm; i++) {
            if (strncmp(Requetes[set_nb].nomvars.sdata[i],nomvar,Min(4,strlen(Requetes[set_nb].nomvars.sdata[i]))) == 0) {
              amatch = -1; /* enregistrement a exclure */
              dbprint(stddebug,"Debug c_fst_match_req match exclure\n");
              goto Typvars;
            }
          }
          if (amatch == -1) amatch = 0;   /* exclusion additive non satisfaite */
          continue;   /* rien trouve a exclure */
        }
      }

    Typvars:
      if (Requetes[set_nb].typvars.in_use) {
      dbprint(stddebug,"Debug c_fst_match_req verifie typvars set_nb=%d\n",set_nb);
        if (Requetes[set_nb].exdes == desire) {
          for (i=0; i < Requetes[set_nb].typvars.nelm; i++)
            if (strncmp(Requetes[set_nb].typvars.sdata[i],typvar,Min(2,strlen(Requetes[set_nb].typvars.sdata[i]))) == 0) {
              amatch = 1;
              dbprint(stddebug,"Debug c_fst_match_req match desire\n");
              goto Dates; /* satisfait la requete */
            }
          continue;  /* rien trouve qui satisfait la requete pour typvars */
        }
        else {
          for (i=0; i < Requetes[set_nb].typvars.nelm; i++)
            if (strncmp(Requetes[set_nb].typvars.sdata[i],typvar,Min(2,strlen(Requetes[set_nb].typvars.sdata[i])))== 0) {
              amatch = -1;  /* enregistrement a exclure */
              dbprint(stddebug,"Debug c_fst_match_req match exclure\n");
              goto Dates;
            }
          if (amatch == -1) amatch = 0;   /* exclusion additive non satisfaite */
          continue;   /* rien trouve a exclure */
        }
      }

    Dates:
      if (Requetes[set_nb].dates.in_use) {
          switch (Requetes[set_nb].dates.arg_type) {

          case entier:
            dbprint(stddebug,"Debug c_fst_match_req verifie dates entier set_nb=%d\n",set_nb);
            if (Requetes[set_nb].exdes == desire) {
              for (i=0; i < Requetes[set_nb].dates.nelm; i++)
                if (Requetes[set_nb].dates.data.tab_elem[i] == date) {
                  amatch = 1;
                  dbprint(stddebug,"Debug c_fst_match_req match desire\n");
                  goto Ip1s; /* satisfait la requete */
                }
              continue;  /* rien trouve qui satisfait la requete pour dates */
            }
            else {
              for (i=0; i < Requetes[set_nb].dates.nelm; i++)
                if (Requetes[set_nb].dates.data.tab_elem[i] == date) {
                  amatch = -1;  /* enregistrement a exclure */
                  dbprint(stddebug,"Debug c_fst_match_req match exclure\n");
                  goto Ip1s;
                }
              if (amatch == -1) amatch = 0;   /* exclusion additive non satisfaite */  
              continue;  /* rien trouve a exclure */
            }
            break;

          case deb_fin_entier:
            if (Requetes[set_nb].dates.data.tab_elem[0] == -1)
              debut = date;
            else {
              debut = Requetes[set_nb].dates.data.tab_elem[0];
              }
            if (Requetes[set_nb].dates.data.tab_elem[1] == -1)
              fin = date;
            else
              fin = Requetes[set_nb].dates.data.tab_elem[1];
            dbprint(stddebug,"Debug c_fst_match_req verifie dates debut=%d fin=%d set_nb=%d\n",debut,fin,set_nb);
            f77name(difdatr)(&date,&debut,&diff_deb);
            f77name(difdatr)(&date,&fin,&diff_fin);
            dbprint(stddebug,"Debug diff_deb=%f diff_fin=%f\n",diff_deb,diff_fin);
            if ((diff_deb >= 0.) && (diff_fin <= 0.))
              if (Requetes[set_nb].exdes == desire) {
                amatch = 1;
                dbprint(stddebug,"Debug c_fst_match_req match desire\n");
                goto Ip1s;  /* satisfait la requete */
                }
              else {
                amatch = -1;  /* enregistrement a exclure */
                dbprint(stddebug,"Debug c_fst_match_req match exclure\n");
                goto Ip1s;
              }
            if (amatch == -1) amatch = 0;   /* exclusion additive non satisfaite */  
            continue;  /* rien trouve qui satisfait la requete desire/exclure */
            break;

          case deb_fin_entier_delta:
            dbprint(stddebug,"Debug c_fst_match_req verifie dates debut fin delta set_nb=%d\n",set_nb);
            if (Requetes[set_nb].dates.data.tab_elem[0] == -1)
              debut = date;
            else {
              debut = Requetes[set_nb].dates.data.tab_elem[0];
              }
            if (Requetes[set_nb].dates.data.tab_elem[1] == -1)
              fin = date;
            else
              fin = Requetes[set_nb].dates.data.tab_elem[1];
            delta8 = Requetes[set_nb].dates.delta;
            dbprint(stddebug,"Debug c_fst_match_req verifie dates debut=%d fin=%d delta=%f\n",debut,fin,delta8);
            f77name(difdatr)(&date,&debut,&diff_deb);
            f77name(difdatr)(&date,&fin,&diff_fin);
            remainder = fmod(diff_deb,delta8);
            dbprint(stddebug,"Debug diff_deb=%f diff_fin=%f modulo=%f\n",diff_deb,diff_fin,remainder);
            if ((diff_deb >= 0.) && (diff_fin <= 0.) && (remainder <= (5.0/3600.)))
              if (Requetes[set_nb].exdes == desire) {
                amatch = 1;
                dbprint(stddebug,"Debug c_fst_match_req match desire\n");
                goto Ip1s;  /* satisfait la requete */
              }
              else {
                amatch = -1;  /* enregistrement a exclure */
                dbprint(stddebug,"Debug c_fst_match_req match exclure\n");
                goto Ip1s;
              }
            if (amatch == -1) amatch = 0;   /* exclusion additive non satisfaite */  
            continue;  /* rien trouve qui satisfait la requete desire/exclure */
            break;

          default:
            fprintf(stderr,"error: c_fst_match_req invalid Requetes[%d].dates.arg_type=%d\n",
                  set_nb,Requetes[set_nb].dates.arg_type);
            return(0);
            break;

          } /* end switch */
      }

    Ip1s:
      if (Requetes[set_nb].ip1s.in_use) {
        switch (Requetes[set_nb].ip1s.arg_type) {

        case reel:
          f77name(convip_plus)(&ip1,&p1,&ip_kind,&mode,string,&flag);
          dbprint(stddebug,"Debug c_fst_match_req verifie ip1s du fichier=%d reel=%f set_nb=%d\n",ip1,p1,set_nb);
          if (Requetes[set_nb].exdes == desire) {
            for (i=0; i < Requetes[set_nb].ip1s.nelm; i++)
              if (p1 == 0.) {
                if (Requetes[set_nb].ip1s.data.tab_rval[i] == 0.) {
                  amatch = 1;
                  dbprint(stddebug,"Debug c_fst_match_req match desire\n");
                  goto Ip2s; /* satisfait la requete */
                }
              }    
              else if ((Abs(1.0 - (Requetes[set_nb].ip1s.data.tab_rval[i]/p1)) < err_tolerance) &&
                       (ip_kind == Requetes[set_nb].ip1s.ip_kind)) {
                amatch = 1;
                dbprint(stddebug,"Debug c_fst_match_req match desire\n");
                goto Ip2s; /* satisfait la requete */
              }
            continue;  /* rien trouve qui satisfait la requete pour ip1s */
          }
          else {
            for (i=0; i < Requetes[set_nb].ip1s.nelm; i++)
              if (p1 == 0.) {
                if (Requetes[set_nb].ip1s.data.tab_rval[i] == 0.) {  
                  amatch = -1;  /* enregistrement a exclure */
                  dbprint(stddebug,"Debug c_fst_match_req match exclure\n");
                  goto Ip2s;
                }
              }  
              else if ((Abs(1.0 - (Requetes[set_nb].ip1s.data.tab_rval[i]/p1)) < err_tolerance) &&
                       (ip_kind == Requetes[set_nb].ip1s.ip_kind)) {
                amatch = -1;  /* enregistrement a exclure */
                dbprint(stddebug,"Debug c_fst_match_req match exclure\n");
                goto Ip2s;
              }
            if (amatch == -1) amatch = 0;   /* exclusion additive non satisfaite */  
            continue;   /* rien trouve a exclure */
          }
          break;

        case deb_fin_reel:
          f77name(convip_plus)(&ip1,&p1,&ip_kind,&mode,string,&flag);
          if (Requetes[set_nb].ip1s.data.tab_rval[0] == -1.0)
            r_debut = p1;
          else
            r_debut = Requetes[set_nb].ip1s.data.tab_rval[0];
          if (Requetes[set_nb].ip1s.data.tab_rval[1] == -1.0)
            r_fin = p1;
          else
            r_fin = Requetes[set_nb].ip1s.data.tab_rval[1];
          dbprint(stddebug,"Debug c_fst_match_req verifie ip1s r_debut=%f r_fin=%f\n",r_debut,r_fin);
          dbprint(stddebug,"Debug+ c_fst_match_req p1=%f err_tolerance=%f\n",p1,err_tolerance);
          if ((p1 >= (r_debut-err_tolerance)) && (p1 <= (r_fin+err_tolerance)) &&
              (ip_kind == Requetes[set_nb].ip1s.ip_kind)) {
            if (Requetes[set_nb].exdes == desire) {
              amatch = 1;
              dbprint(stddebug,"Debug c_fst_match_req match desire\n");
              goto Ip2s;  /* satisfait la requete */
            }
            else {
              amatch = -1;  /* enregistrement a exclure */
              dbprint(stddebug,"Debug c_fst_match_req match exclure\n");
              goto Ip2s;
            }
          }
          else {
            if (amatch == -1) amatch = 0;   /* exclusion additive non satisfaite */
            continue;  /* rien trouve qui satisfait la requete pour ip1s */
            }
          break;

        default:
          fprintf(stderr,"error: c_fst_match_req invalid Requetes[%d].ip1s.arg_type=%d\n",
                set_nb,Requetes[set_nb].ip1s.arg_type);
          return(0);
          break;


        } /* end switch */
      }

    Ip2s:
      if (Requetes[set_nb].ip2s.in_use) {
      dbprint(stddebug,"Debug c_fst_match_req verifie ip2s set_nb=%d\n",set_nb);
        switch (Requetes[set_nb].ip2s.arg_type) {

        case entier:
          if (Requetes[set_nb].exdes == desire) {
            for (i=0; i < Requetes[set_nb].ip2s.nelm; i++)
              if (Requetes[set_nb].ip2s.data.tab_elem[i] == ip2) {
                amatch = 1;
                dbprint(stddebug,"Debug c_fst_match_req match desire\n");
                goto Ip3s; /* satisfait la requete */
              }
            continue;  /* rien trouve qui satisfait la requete pour ip2s */
          }
          else {
            for (i=0; i < Requetes[set_nb].ip2s.nelm; i++)
              if (Requetes[set_nb].ip2s.data.tab_elem[i] == ip2) {
                amatch = -1;  /* enregistrement a exclure */
                dbprint(stddebug,"Debug c_fst_match_req match exclure\n");
                goto Ip3s;
              }
            if (amatch == -1) amatch = 0;   /* exclusion additive non satisfaite */  
            continue;   /* rien trouve a exclure */
          }
          break;

        default:
          fprintf(stderr,"error: c_fst_match_req invalid Requetes[%d].ip2s.arg_type=%d\n",
                set_nb,Requetes[set_nb].ip2s.arg_type);
          return(0);
          break;

        } /* end switch */
      }

    Ip3s:
      if (Requetes[set_nb].ip3s.in_use) {
      dbprint(stddebug,"Debug c_fst_match_req verifie ip3s set_nb=%d\n",set_nb);
        switch (Requetes[set_nb].ip3s.arg_type) {

        case entier:
          if (Requetes[set_nb].exdes == desire) {
            for (i=0; i < Requetes[set_nb].ip3s.nelm; i++)
              if (Requetes[set_nb].ip3s.data.tab_elem[i] == ip3) {
                amatch = 1;
                dbprint(stddebug,"Debug c_fst_match_req match desire\n");
                goto Fin; /* satisfait la requete */
              }
            continue;  /* rien trouve qui satisfait la requete pour ip3s */
          }
          else {
            for (i=0; i < Requetes[set_nb].ip3s.nelm; i++)
              if (Requetes[set_nb].ip3s.data.tab_elem[i] == ip3) {
                amatch = -1;  /* enregistrement a exclure */
                dbprint(stddebug,"Debug c_fst_match_req match exclure\n");
                goto Fin;
              }
            if (amatch == -1) amatch = 0;   /* exclusion additive non satisfaite */  
            continue;   /* rien trouve a exclure */
          }
          break;

        default:
          fprintf(stderr,"error: c_fst_match_req invalid Requetes[%d].ip3s.arg_type=%d\n",
                set_nb,Requetes[set_nb].ip3s.arg_type);
          return(0);
          break;

        } /* end switch */
      }

    Fin:
      if (amatch == 1) {
        dbprint(stddebug,"Debug c_fst_match_req fin requete desire satisfaite, handle=%d \n",handle);
        return(1);            /* requete desire satisfaite */
      }
      else if (amatch == -1) {
        dbprint(stddebug,"Debug c_fst_match_req fin requete exclure satisfaite handle=%d \n",handle);
        return(0);            /* requete exclure satisfaite */
      }

    /* Next:                  verifier le prochain desire/exclure */

  } /* end for */

  last_in_use = last_R;
  while ((last_in_use > first_R) && (! Requetes[last_in_use].in_use))
    last_in_use--;
/*  fprintf(stderr,"Debug+ last_in_use=%d\n",last_in_use); */
/* rien trouve qui satisfait les requetes */
  if ((Requetes[last_in_use].exdes == desire) && (Requetes[last_in_use].in_use))
    return(0);  /* ne satisfait pas la requete desire */
  else
    return(1);  /* rien a exclure */

} /* end fst_match_req */

/*****************************************************************************
 *                     L E S   I N T E R F A C E S                           *
 *                                                                           *
 *Objet                                                                      *
 *   Les interfaces Fortran a C ainsi que les interfaces aux                 *
 *   fonctions C a liste d'arguments raccourcie aux fonctions                *
 *   a liste d'arguments etendue                                             *
 *                                                                           *
 *****************************************************************************/

int c_select_ip1(void *iplist, int nelm, int kind)
{
/*CHC/NRC*/
	return xc_select_ip1(bundle_nb, desire_exclure, iplist, nelm, kind);
}

int c_select_ip2(void *iplist, int nelm, int kind)
{
/*CHC/NRC*/
	return xc_select_ip2(bundle_nb, desire_exclure, iplist, nelm, kind);
}

int c_select_ip3(void *iplist, int nelm, int kind)
{
/*CHC/NRC*/
	return xc_select_ip3(bundle_nb, desire_exclure, iplist, nelm, kind);
}

int c_select_date(int set_nb, int des_exc, int *date_list, int nelm, float delta)
{
/*CHC/NRC*/
	return xc_select_date(bundle_nb, desire_exclure, date_list, nelm, delta);
}

int c_select_etiquette(char *etiq_list[], int nelm)
{
/*CHC/NRC*/
	return xc_select_etiquette(bundle_nb, desire_exclure, etiq_list, nelm);
}

int c_select_nomvar(char *nomv_list[], int nelm)
{
/*CHC/NRC*/
	return xc_select_nomvar(bundle_nb, desire_exclure, nomv_list, nelm);
}

int c_select_typvar(char *typv_list[], int nelm)
{
/*CHC/NRC*/
	return xc_select_typvar(bundle_nb, desire_exclure, typv_list, nelm);
}

int f77name(select_ip1)(void *iplist, int *nelm, int *kind)
{
  return(xc_select_ip1(bundle_nb, desire_exclure,iplist,*nelm,*kind));
}

int f77name(select_ip2)(void *iplist, int *nelm, int *kind)
{
  return(xc_select_ip2(bundle_nb, desire_exclure,iplist,*nelm,*kind));
}

int f77name(select_ip3)(void *iplist, int *nelm, int *kind)
{
  return(xc_select_ip3(bundle_nb, desire_exclure,iplist,*nelm,*kind));
}

int f77name(select_date)(int *date_list, int *nelm, float *delta)
{
  return(xc_select_date(bundle_nb, desire_exclure,date_list,*nelm,*delta));
}

int f77name(select_etiquette)(char *etiq_list[], int *nelm, F2Cl flng)
{
  char **string_array;
  int i,ier,lng=flng;
  dbprint(stddebug,"Debug desire_etiquette lng=%d nelm=%d\n",lng,*nelm);
  string_array = fill_string_array(allocate_string_array(*nelm),etiq_list,lng,*nelm,0);
  for (i=0; i < *nelm; i++)
      dbprint(stddebug,"Debug string_array[%d]-->%s<--\n",i,string_array[i]);
  ier = xc_select_etiquette(bundle_nb, desire_exclure,string_array,*nelm);
  free_string_array(string_array);
  return(ier);
}

int f77name(select_nomvar)(char *nomv_list[], int *nelm, F2Cl flng)
{
  char **string_array;
  int i,ier,lng=flng;
  dbprint(stddebug,"Debug desire_nomvar lng=%d nelm=%d\n",lng,*nelm);
  string_array = fill_string_array(allocate_string_array(*nelm),nomv_list,lng,*nelm,0);
  for (i=0; i < *nelm; i++)
      dbprint(stddebug,"Debug string_array[%d]-->%s<--\n",i,string_array[i]);
  ier = xc_select_nomvar(bundle_nb, desire_exclure,string_array,*nelm);
  free_string_array(string_array);
  return(ier);
}

int f77name(select_typvar)(char *typv_list[], int *nelm, F2Cl flng)
{
  char **string_array;
  int i,ier,lng=flng;
  dbprint(stddebug,"Debug desire_typvar lng=%d nelm=%d\n",lng,*nelm);
  string_array = fill_string_array(allocate_string_array(*nelm),typv_list,lng,*nelm,0);
  for (i=0; i < *nelm; i++)
      dbprint(stddebug,"Debug string_array[%d]-->%s<--\n",i,string_array[i]);
  ier = xc_select_typvar(bundle_nb, desire_exclure,string_array,*nelm);
  free_string_array(string_array);
  return(ier);
}

int f77name(filtre_desire)()
{
/*CHC/NRC*/
	return c_filtre_desire();
}

int f77name(filtre_exclure)()
{
/*CHC/NRC*/
	return c_filtre_exclure();
}
  
int f77name(fst_match_req)(int *handle)
{
  return(c_fst_match_req(*handle));
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
int c_requetes_reset(int set_nb, int nomvars, int typvars, int etikets, int dates,
                    int ip1s, int ip2s, int ip3s)
{
  int j;

  if (set_nb > MAX_requetes) {
    fprintf(stderr,"error: c_requetes_reset set_nb=%d > MAX_requetes=%d\n",set_nb,MAX_requetes);
    return(-1);
  }

  Requetes[set_nb].in_use = 0;
  Requetes[set_nb].exdes = -1;
  if (nomvars != -1) {
    Requetes[set_nb].nomvars.in_use = 0;
    Requetes[set_nb].nomvars.nelm = 0;
    for (j=0; j < MAX_Nlist; j++)
      strcpy(Requetes[set_nb].nomvars.sdata[j],"    ");
  }
  if (typvars != -1) {
    Requetes[set_nb].typvars.in_use = 0;
    Requetes[set_nb].typvars.nelm = 0;
    for (j=0; j < MAX_Nlist; j++)
      strcpy(Requetes[set_nb].typvars.sdata[j],"  ");
  }
  if (etikets != -1) {
    Requetes[set_nb].etiquettes.in_use = 0;
    Requetes[set_nb].etiquettes.nelm = 0;
    for (j=0; j < MAX_Nlist; j++)
      strcpy(Requetes[set_nb].etiquettes.sdata[j],"            ");
  }
  if (dates != -1) {
    Requetes[set_nb].dates.in_use = 0;
    Requetes[set_nb].dates.nelm = 0;
    for (j=0; j < MAX_Nlist; j++)
      Requetes[set_nb].dates.data.tab_elem[j]=0;
  }
  if (ip1s != -1) {
    Requetes[set_nb].ip1s.in_use = 0;
    Requetes[set_nb].ip1s.nelm = 0;
    for (j=0; j < MAX_Nlist; j++)
       Requetes[set_nb].ip1s.data.tab_elem[j]=0;
  }
  if (ip2s != -1) {
    Requetes[set_nb].ip2s.in_use = 0;
    Requetes[set_nb].ip2s.nelm = 0;
    for (j=0; j < MAX_Nlist; j++)
      Requetes[set_nb].ip2s.data.tab_elem[j]=0;
  }
  if (ip3s != -1) {
    Requetes[set_nb].ip3s.in_use = 0;
    Requetes[set_nb].ip3s.nelm = 0;
    for (j=0; j < MAX_Nlist; j++)
       Requetes[set_nb].ip3s.data.tab_elem[j]=0;
  }

  return(0);
}
/*****************************************************************************
 *                      C _ R E Q U E T E S _ I N I T                        *
 *                                                                           *
 *Objet                                                                      *
 *   Initialiser le package de requetes                                      *
 *                                                                           *
 *****************************************************************************/
void c_requetes_init(char *requetes_filename, char *debug_filename)
{
  int i,j,ier;

/*  debug_filename = getenv("DEBUGFILE"); */
  if (debug_filename != NULL)
    stddebug = fopen(debug_filename,"w");
  else
    stddebug = fopen("/dev/null","w");

  for (i=0; i<MAX_requetes; i++) {
    for (j=0; j<MAX_Nlist; j++) {
      Requetes[i].etiquettes.sdata[j] = (char *) malloc(13);
      Requetes[i].nomvars.sdata[j] = (char *) malloc(5);
      Requetes[i].typvars.sdata[j] = (char *) malloc(3);
    }
    c_requetes_reset(i,1,1,1,1,1,1,1);
  }
  
  /* requetes_filename = getenv("FST_FILTER_FILE"); */
  if (requetes_filename != NULL)
    ier = c_requetes_read_file(requetes_filename);
}

/*****************************************************************************
 *                      D I R E C T I V E S _ I P 1 2 3                      *
 *                                                                           *
 *Objet                                                                      *
 *   Traiter les directives desires pour ip1 ip2 et ip3                      *
 *                                                                           *
 *****************************************************************************/
int directive_ip123(int argc , char **argv,char cmd_strt,  int *func(), char *Private_Data_2)
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
  	printf("directive_ip123 error: unknown kind =%s\n",argv[1]);
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
int directive_dates(int argc , char **argv,char cmd_strt,  int *func(int set_nb, int des_exc, int *date_list, int nelm, float delta), char *Private_Data_2)
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
  xc_select_date(set_nb,des_exc,list_entier,nelm,delta);
  return(0);
}

/*****************************************************************************
 *                      D I R E C T I V E S _ D A T E V                      *
 *                                                                           *
 *Objet                                                                      *
 *   Traiter les directives desires pour la date (en format visuel)          *
 *                                                                           *
 *****************************************************************************/
int directive_datev(int argc , char **argv,char cmd_strt,  int *func(), char *Private_Data_2)
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
    	printf("directive_datev error: visual date format list must contain an even number of date elements\n");
      return(-1);
      }
    offset=2;
  }
  else {
  	printf("directive_datev error: visual date format list must contain an even number of date elements\n");
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
  xc_select_date(set_nb,des_exc,list_entier,nelm,delta);
  return(0);
}

/*****************************************************************************
 *                      D I R E C T I V E S _ C H A R V A R                  *
 *                                                                           *
 *Objet                                                                      *
 *   Traiter les directives desires pour nomvar, typvar et etiquette         *
 *                                                                           *
 *****************************************************************************/
int directive_charvar(int argc , char **argv,char cmd_strt,  int *func(), char *Private_Data_2)
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

  /* printf("Debug+ directive_charvar appel la fonction avec set_nb=%d nelm=%d\n",set_nb,nelm); */
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
int directive_desire(int argc , char **argv,char cmd_strt,  int *func(), char *Private_Data_2)
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
int directive_exclure(int argc , char **argv,char cmd_strt,  int *func(), char *Private_Data_2)
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
int c_requetes_read_file(char *requetes_file)
{
  FILE *fp;
  int *pv2 = 0;

  fp = fopen(requetes_file,"r");
  if (fp == (FILE *) NULL) {
    fprintf(stderr,"c_requetes_read_file error opening file %s\n",requetes_file);
    fprintf(stderr,"filter file ignored\n");
    return(-1);
    }
  else
    fprintf(stderr,"warning: FSTD software uses filter file %s for scans\n",requetes_file);

  rpn_c_callback("ip1",directive_ip123,"",xc_select_ip1,pv2);
  rpn_c_callback("ip2",directive_ip123,"",xc_select_ip2,pv2);
  rpn_c_callback("ip3",directive_ip123,"",xc_select_ip3,pv2);
  rpn_c_callback("datestamp",directive_dates,"",xc_select_date,pv2);
  rpn_c_callback("datevisual",directive_datev,"",xc_select_date,pv2);
  rpn_c_callback("nomvar",directive_charvar,"",xc_select_nomvar,pv2);
  rpn_c_callback("typvar",directive_charvar,"",xc_select_typvar,pv2);
  rpn_c_callback("etiquette",directive_charvar,"",xc_select_etiquette,pv2);
  rpn_c_callback("desire",directive_desire,"",c_filtre_desire,pv2);
  rpn_c_callback("exclure",directive_exclure,"",c_filtre_exclure,pv2);
  process_c_callback(requetes_file);

	return 0; /*CHC/NRC*/
}
#if defined (TEST)
void f77name(c_main)(int *handle)
{

  int i,j;
  float heures;
  int ip1s_i[] = { 400,500,600,750 };
  float ip1s_r[] = { .3840, .4440, 0.6110 };
  float ip1s_r_range[] = { .3280, 0.8000 };
  int ip1s_range[] = {400, 750};
  int ip2s[] = {0, 12, 24};
  int ip3s[] = {0, 80, 90};
  int dates[] = {313290800, 313301600};
  int dates_range[] = {313280000, 313290800};

  char *testeti[] = { "Etiquette #1", "R2428V4N", "Etiquet #3" };
  char *testnom[2] = { "TT", "GZ" };
  char *testtyp[4] = { "P", "V1", "V2", "V3" };

  dbprint(stddebug,"Debug debut \n");
  dbprint(stddebug,"Debug testeti=%s %s %s \n",testeti[0],testeti[1],testeti[2]);
  dbprint(stddebug,"Debug testeti=%s %s %s \n",testeti[0],testeti[1],testeti[2]);
  dbprint(stddebug,"Debug testnom=%s %s \n",testnom[0],testnom[1]);
  dbprint(stddebug,"Debug testtyp=%s %s %s %s \n",testtyp[0],testtyp[1],testtyp[2],testtyp[3]);

  c_requetes_init();

  i = xc_select_ip1(2,1,ip1s_i,4,-1);
  i = xc_select_ip1(0,1,ip1s_r,3,1);
/*  i = xc_select_ip1(1,1,ip1s_range,-2,-1);*/
  i = xc_select_ip1(1,1,ip1s_r_range,-2,1);
  i = xc_select_ip2(2,1,ip2s,3,-1);
  i = xc_select_ip3(2,1,ip3s,3,-1);
  heures=12.0;
  i = xc_select_date(1,1,dates_range,-2,heures);
/*  heures=0.0;
  i = xc_select_date(1,1,dates_range,-2,heures);*/
  i = xc_select_etiquette(2,1,testeti,3);
  i = xc_select_nomvar(2,1,testnom,2);
  i = xc_select_typvar(2,1,testtyp,4);
 /* j = c_fst_match_req(1,*handle);*/
  j = c_fst_match_req(*handle);

}
#endif

