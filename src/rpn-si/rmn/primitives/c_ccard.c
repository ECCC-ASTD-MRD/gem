/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
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


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define NCARMAX 256
#define NKLEMAX 100
#define MAJUS 0
#define MINUS 1
#define PAREIL 2
#define FAUX 0
#define VRAI 1
/*************************************************
 * definition de la structure des clefs          *
 *************************************************/
struct c_jfc_definition
{
    char *clenom;
    char *cledef;
    char *cleval;
    char *clefin;
    int cletype;
};

/*******************************************************************
 *    recuperation des parametres d'appel a un programme           *
 *******************************************************************/
void c_ccard(char **argv, int argc, char **cle, char val[][NCARMAX],
             char **def, int n, int *npos)

/*
 *   James Caveen - RPN
 *
 *   Revision: J. Caveen - fevrier 1993 : modification permettant
 *                         de passer des parametres positionels
 *                         avant les clefs ainsi qu'apres les signes --
 *
 *   pour chacune des clefs, on reserve une structure dans laquelle
 *   on conserve son nom et son type (i.e., convertir a MAJUScule,
 *   MINUScule ou pas-touche), les valeurs de defauts val et def
 *   ainsi que la valeur finale attribuee a la clef.  C'est le
 *   contenu de clefin qui est retourne.
 *
 *   Si a l'entree, npos est plus grand que 0, a la sortie, il
 *   contiendra le nombre de parametres positionels.
 *   Si a l'entree, npos = -111, la premiere erreur rencontree
 *   provoquera la fin prematuree du programme appelant.
 */
{
  void free(), exit();

  int c_jfc_cherche_la_clef();
  void c_jfc_tradup(), c_jfc_traduire(),  c_jfc_les_valeurs();
  void sequence_appel();
  char **c_jfc_positionel();

  int c_jfc_majmin();
  struct c_jfc_definition les_clefs[NKLEMAX];
  int plante = FAUX;
  int retourne_npos =  FAUX;
  char *keyname, *pointeur,*prognom;

  int   pos, posc, ii, i, posmoin, npos_interne;
  int erreur = 0;

  /*  reserver la memoire et initialiser les n structures de clefs  */

  for (i=0; i<n; i++)
  {
    les_clefs[i].clenom = calloc(strlen(cle[i])+1,sizeof(char));
    les_clefs[i].cledef = calloc(strlen(def[i])+1,sizeof(char));
    les_clefs[i].cleval = calloc(NCARMAX+1,sizeof(char));
    les_clefs[i].clefin = calloc(NCARMAX+1,sizeof(char));
  }

  if (*npos == -111) {
    plante = VRAI;
  }

  if (*npos > 0)
  {
    *npos = 0;
    retourne_npos = VRAI;
  }

  for (i=0; i<n; i++)
  {
    c_jfc_tradup (cle[i],les_clefs[i].clenom);
    les_clefs[i].cletype = c_jfc_majmin(les_clefs[i].clenom);
    strcpy(les_clefs[i].cledef, def[i]);
    strcpy(les_clefs[i].cleval, val[i]);
    strcpy(les_clefs[i].clefin, val[i]);
  }

  pos = posc = -1;

  /*  obtenir le nom du programme */
  prognom = *argv;
  argv++;

  if (argc == 2)
  {
    if ((strcmp(*argv,"-h") == 0) || (strcmp(*argv,"-H")==0))
    {
      sequence_appel(les_clefs,prognom,n);
      exit(0);
    }
  }

  /*
   *  verifier si on peut recuperer en mode positionel
   *  on cherche la premiere clef '-' dans la liste
   */
  posmoin = c_jfc_cherche_la_clef("-",les_clefs,n);

  /*
   *  recuperer les premiers arguments positionels
   */
  if (posmoin > -1)
  {
    argv = c_jfc_positionel(argv,les_clefs,n,posmoin,&npos_interne,1,
                            &erreur);
    if (erreur)
    {
      fprintf(stderr,"\n ***ERREUR - TROP D'ARGUMENTS POSITIONELS \n");
      sequence_appel(les_clefs,prognom,n);
      if (plante)
        exit (-1);
    }
  }

  while (*argv)
  {
    if (**argv == '-')            /* un nom de clef */
    {
      if (*((*argv)+1) == '-')/* fin des clefs  */
      {
        /*on recupere le reste des arguments positionels */
        argv++;
        argv = c_jfc_positionel(argv,les_clefs,n,posmoin,
                                &npos_interne,0,&erreur);
        if(erreur)
        {
          fprintf(stderr,"\n ***ERREUR - TROP D'ARGUMENTS POSITIONELS \n");
          sequence_appel(les_clefs,prognom,n);
          if( plante)
            exit (-1);
        }
        break;
      }

      pointeur = (strchr(*argv,'='));
      if(pointeur != (char *) NULL)
      {
        *pointeur = '\0';
      }
      keyname = (*argv)+1;
      pos = c_jfc_cherche_la_clef(keyname,les_clefs,n);
      if (pos > -1)
      {
        posc = pos;
        if (pointeur != (char *) NULL)
        {
          *argv = pointeur+1;
          c_jfc_les_valeurs(les_clefs,argv, pos, &posc);
        }
        else
        {
          strcpy( les_clefs[posc].clefin , les_clefs[posc].cledef);
        }
      }
      else
      {
        fprintf(stderr," ***ERREUR CLEF=%s INVALIDE***\n",keyname);
        sequence_appel(les_clefs,prognom,n);
        if (plante)
          exit (-1);
      }
    }
    else
    {
      if ((posc != -1) && (strcmp(*(cle+posc),*(cle+pos)) == 0))
      {
        c_jfc_les_valeurs(les_clefs, argv, pos, &posc);
      }
      else
      {
        fprintf(stderr,"\n ***DEBORDEMENT DE LISTE \n");
        sequence_appel(les_clefs,prognom,n);
        if( plante)
          exit (-1);
      }
    }

    argv++;
  } /* while (*argv) */

  /*    recopier les valeurs finales des clefs dans val
        et liberer la memoire                           */

  for (i=0;i<n;i++)
  {
    c_jfc_traduire(les_clefs[i].clefin,les_clefs[i].cletype);
    strcpy(val[i],les_clefs[i].clefin);
    free(les_clefs[i].clenom);
    free(les_clefs[i].cledef);
    free(les_clefs[i].cleval);
    free(les_clefs[i].clefin);
  }

  /*
   *    on retourne le nombre d'arguments positionels si
   *    a l'entree, npos > 0
   */
  if (retourne_npos)
  {
    (*npos) = npos_interne;
  }
}

/********************************************************
 *  chercher une clefs dans la liste de noms de clefs   *
 *  et retourner sa position dans la liste.             *
 ********************************************************/

int c_jfc_cherche_la_clef(char *nom,struct c_jfc_definition cle[],int n)
{
  int toupper();
  int i=0;

  char tmpstr[256];

  strcpy(tmpstr, nom);
  while (*(tmpstr+i))
  {
    *(tmpstr+i) = toupper(*(tmpstr+i));
    i++;
  }

  for (i=0; i<n; i++)
  {
    if (strcmp(tmpstr,cle[i].clenom) == 0)
    {
      return(i);
    }
  }
  return(-1);
}


/****************************************************************
 *  extraire les elements d'une liste                           *
 ****************************************************************/

void c_jfc_les_valeurs(struct c_jfc_definition clefs[], char **argv,
                       int pos, int *posc)

     /*
      *   pour l'element contenu dans argv, on fait les operations
      *   suivantes:
      *              on elimine le signe = si present,
      *              on elimine les delimiteurs ':' si presents
      *              et on attribue les valeurs ainsi obtenues
      *              a la clef a traiter.
      *
      *          exemple:   *argv =  =-12:34          donnera
      *                               -12 34     qui seront associes
      *                      a deux clefs contigues ayant le meme nom dans la liste des clefs.
      */
{

  char *deux_pt;
  int nombre;

  for (nombre = 0; ;)
  {
    if ((*posc > -1) && (strcmp(clefs[*posc].clenom,clefs[pos].clenom)) == 0)
    {
      deux_pt = strchr(*argv,':');
      if (deux_pt != (char *) NULL)
      {
        *deux_pt = '\0';
        nombre = strlen(*argv);
        if (**argv == '=')
        {
          strcpy( clefs[*posc].clefin, (*argv)+1);
        }
        else
        {
          strcpy( clefs[*posc].clefin, *argv);
        }
        *argv =deux_pt + 1;
        (*posc)++;
      }
      else
      {
        if (**argv == '=')
        {
          strcpy( clefs[*posc].clefin,(*argv)+1);
        }
        else
        {
          strcpy( clefs[*posc].clefin,*argv);
        }
        (*posc)++;
        break;
      }
    }
    else
    {
      fprintf(stderr,"\n***ERREUR DEBORDEMENT DE LISTE  OU MODE POSITIONNEL\n");
    }

  }
}

/******************************************************
 *  traduire un nom a des MAJUScules                  *
 ******************************************************/

void c_jfc_tradup(char *nom, char *nommaj)

{
  int toupper();

  /*   nettoyer nommaj  */
  while(*nommaj)
  {
    *nommaj = '\0';
    nommaj++;
  }

  while(*nom)
  {
    *nommaj = toupper(*nom);
    nommaj++;  nom++;
  }
}

/**************************************************************
 *    determiner le type de clef MAJUScule/MINUScule/tel quel *
 **************************************************************/


int c_jfc_majmin(char *arg)

     /*     fonction servant a determiner si la valeur a donner a une
      *      clef est de type MAJUScule, MINUScule ou si elle reste telle
      *      que decrite lors de l'appel.  la fonction retourne le valeur
      *      du type de clef (MAJUS MINUS ou PAREIL) et met le dernier
      *      element de arg egal a nul si necessaire.
      */
{

  int lng;
  lng = strlen(arg) - 1;

  if (*(arg+lng) == '.')
  {
    *(arg+lng) = '\0';
    return(PAREIL);
  }
  if (*(arg+lng) == '_')
  {
    *(arg+lng) = '\0';
    return(MINUS);
  }
  return(MAJUS);

}

/*******************************************************************
 * traduire une valeur a MAJUScule/MINUScule selon le type de clef *
 *******************************************************************/
void c_jfc_traduire(char *cle,int cletype)

     /* pour chaque valeur de clef, on verifie si on doit faire une conversion
      * a des MAJUScules, MINUScules ou si elle doit rester telle quelle
      */
{

  if (cletype == MAJUS)
  {
    while(*cle)
    {
      *cle = toupper(*cle);
      cle++;
    }
  }
  else if (cletype == MINUS)
  {
    while(*cle)
    {
      *cle = tolower(*cle);
      cle++;
    }
  }
}


/******************************************************
 *       imprimer sur stderr la sequence d'appel      *
 ******************************************************/

void sequence_appel(struct c_jfc_definition defo[],char *scriptnom,int n)

     /*   fonction servant a imprimer sur le fichier stderr la liste des
      *    noms de clefs et leurs valeurs de defaut
      *
      *  argument
      *
      *         defo    -    structure contenant les noms de clefs ainsi que leurs
      *                      valeurs de defaut et le type.
      *
      *  cette fonction est appelle en cas d'erreur ou lorsque l'usager invoque
      *  le programme appelant avec la clef -h.
      *               prognom -h
      *
      */

{

  int i = 0;

  fprintf(stderr,"\n *** SEQUENCE D'APPEL ***\n\n");

  fprintf(stderr,"%s \n",scriptnom);

  for (i=0; i<n; i++)
  {
    fprintf(stderr,"          -%s [%s:%s]\n",defo[i].clenom,defo[i].cledef,defo[i].cleval);
  }
  fprintf(stderr,"\n");
}


/*****************************************************
 *   recuperation des parametres positionels         *
 *   lors de l'appel a un programme.                 *
 *   Chaque parametre est copie tel quel dans        *
 *   l'element clefin de la structure de clefs       *
 *****************************************************/

char **c_jfc_positionel(char **argv, struct c_jfc_definition les_clefs[],
                        int n, int pos, int *npos, int debut, int *erreur)

     /*
      * cette fonction recupere les parametres positionels (c.a.d. ceux
      * qui sont associes a une clef ayant pour nom '-'.
      *
      * Lorsque debut est VRAI, on recupere les parametres jusqu'a la rencontre
      * d'un nom de clef et on retourne la position du prochain argument contenu dans argv.
      *
      * Lorsque debut est FAUX, on recupere les parametres positionels qui
      * suivent la sequence '--'. On recupere alors jusqu'a epuisement de argv.
      *
      * La fonction retourne dans npos le total cumulatif des arguments positionels traites
      * et retourne dans erreur le cumul des erreurs rencontrees
      */
{

  static int posc;

  if (debut) {
    posc = pos;
  }

  while(*argv)
  {
    if((**argv == '-') && (debut))
    {
      return(argv);
    }
    if(posc >= n)
    {
      argv++;
      (*erreur)++;
      return(argv);
    }
    if(strcmp(les_clefs[posc].clenom,les_clefs[pos].clenom) == 0)
    {
      strcpy(les_clefs[posc].clefin,*argv);
    }
    else
    {
      (*erreur)++;
    }
    posc++;
    (*npos)++;
    argv++;
  }
  return(argv);
}
