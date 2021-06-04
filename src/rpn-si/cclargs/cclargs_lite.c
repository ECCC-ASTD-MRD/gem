#include <stdio.h>
#include <stdlib.h>

#ifndef _XOPEN_SOURCE_EXTENDED
#define _XOPEN_SOURCE_EXTENDED 
#endif
#ifdef NOUI
vmenu()
{
printf("cclargs ERROR: Interactive mode NOT SUPPORTED\n");
exit(1);
return 1;
}
#endif
/*************************************************************
*    recuperation des parametres d'appel a un script         *
*    et transformation en variables                          *
*                                                            *
*    refonte du programme cclargs ecrit par m. valin         *
*    effectuee par James Caveen                              *
*                                                            *
*    modifications subsequentes par m.valin :                *
*    support de variables de sortie,                         *
*    sortie pour perl, python                                *
*    possibilite de changer le delimiteur multi arguments    *
*    introduction de la version sans curses                  *
*    exec a la version curses si necessaire                  *
*                                                            *
**************************************************************/
#include <stdio.h>
#include <string.h>
#include <malloc.h>

#define NKLEMAX 1024

static char OUTBUF[40960];  /* buffer containing output, if exec to curses version
                               is necessary, buffer is discarded  */
static char *OUTBUFPTR=&OUTBUF[0];

static char delimiter=':';
enum interpreter {perl,python,shell} ;
static enum interpreter interp=shell;
static char *separateur="";
static char *virgule=",";
static char UI=1;

enum typecle {majus,minus,pareil};
static int NKLES=0;
int on_affiche = 0;
struct definition
{
    char *kle_nom;        /* nom de la clef */
    char *kle_def1;       /* premiere valeur de defaut */
    char *kle_def2;       /* deuxieme valeur de defaut */
    char *kle_val;        /* valeur finale donnee a la clef */
    char *kle_desc;       /* descripteur pour aide interactive */
    enum typecle type;
};


main(argc, argv)
int argc;
char **argv;
{
  /* save original values of argc and argv for an eventual exec to the GUI version */
  char **ARGV0=argv;

  char **ini_list(), **positionel(), **equivalence();
  void apres_moin_moin(), sequence_appel();
  void imprime(), interactif(), getnom();
  int obligatoire();
  char **pointeur;
  static char scriptnom[50], help_general[182];
  int i, lng;
  char *temp;

  int status = 0;  /* nombre d'erreurs rencontrees */

  struct definition defo[NKLEMAX];


/*initialisation des structures de definition  */

  for(i = 0; i < NKLEMAX ; i++)
  {
   defo[i].kle_nom=defo[i].kle_def1=defo[i].kle_def2=defo[i].kle_val=defo[i].kle_desc=NULL;
  }

  argv++;  /* on saute le nom du programme */
  if( ! strcmp(*argv,"-NOUI") )
  {
     argv++;
     UI=0;
  }
  if( ! strcmp(*argv,"-D") )
  {
     argv++;
     delimiter=**argv;
     fprintf(stderr,"Changing multi-value argument delimiter from : to %c\n",delimiter);
     argv++;
  }
  if( ! strcmp(*argv,"-python") )
  {
     interp=python;
     argv++;
  }
  if( ! strcmp(*argv,"-perl") )
  {
     interp=perl;
     argv++;
  }
  if(**argv != '-')
  {
      getnom(scriptnom,*argv,49);
      argv++;
  }
  else
  {
        scriptnom[0] ='\0';
  }

/*
 *verifier si le mode ecran est obligatoire
*/
  if(! strcmp(*argv,"-+"))
  {
     on_affiche = 1;
     argv++;
  }

/* 
 *verifier s'il y a un help general pour le script
*/
  if(**argv == '[')
  {
      lng = strlen(*argv)-1;
      *((*argv)+lng) = '\0';
      strcpy(help_general,(*argv)+1);
      argv++;
  }
       
  /*  dresser la liste des arguments et leurs valeurs de defaut */
  pointeur = ini_list(argv, defo, &status);

  if(status !=0)
  {
          printf(" exit %d ",status);
          exit (status);
  }

/*  emettre la sequence d'appel si demandee   */
    if(*pointeur)
    {
      if(!strcmp(*pointeur,"--help"))
      {
         sequence_appel(defo,scriptnom);
         printf(" exit 0 ;");
         exit (0);
      }
      if(!strcmp(*pointeur,"-h"))
      {
         sequence_appel(defo,scriptnom);
         printf(" exit 0 ;");
         exit (0);
      }
    }
if(interp== shell) {
   /*
     print list of OUTPUT keys
   */
   OUTBUFPTR+=sprintf(OUTBUFPTR,"CCLARGS_OUT_KEYS='");
   for(i=0;i<=NKLES;i++){
     if(*defo[i].kle_nom == '_') {
   /*
       defo[i].kle_nom++;
   */
       OUTBUFPTR+=sprintf(OUTBUFPTR," %s",defo[i].kle_nom);
     }
   }
   OUTBUFPTR+=sprintf(OUTBUFPTR,"';");
   /*
     print list of ALL keys
   */
   temp = defo[0].kle_nom; while ( *temp == '-' ) temp++;
   OUTBUFPTR+=sprintf(OUTBUFPTR,"CCLARGS_KEYS='%s",temp);
   for(i=1;i<=NKLES;i++) {
     temp = defo[i].kle_nom; while ( *temp == '-' ) temp++;
     OUTBUFPTR+=sprintf(OUTBUFPTR," %s",temp);
   }
   OUTBUFPTR+=sprintf(OUTBUFPTR,"';");
}

  /* recuperation des arguments en mode positionnel */
  pointeur = positionel(pointeur);

  /* recuperation des arguments en mode equivalence */
  pointeur = equivalence(pointeur, defo, scriptnom, &status);

  if(status !=0)
  {
          printf(" exit %d ",status);
          sequence_appel(defo,scriptnom);
          exit (status);
   }

/*
 * verifier si toutes les clefs obligatoires on une valeur
*/

    on_affiche = on_affiche ? on_affiche : obligatoire(defo);



/*  passer en mode ecran si demande    */
    if(on_affiche)
    {
#ifndef NOUI
         interactif(defo,scriptnom,help_general);
#else
         /*sprintf(OUTBUF,"cclargs_lite_curses_%3.3d",VERSION);*/
         sprintf(OUTBUF,"cclargs_lite_curses");
/*         fprintf(stderr,"Exec to %s\n",OUTBUF);    */
         if(UI){
           execvp(OUTBUF,ARGV0);
         } else {
           fprintf(stderr,"cclargs_lite ERROR: curses version not available\n");
           exit(1);
         }
         fprintf(stderr,"cclargs_lite ERROR: exec failed to %s\n",OUTBUF);
         exit(1);
#endif
    }

  /* recuperation des parametres positionnels apres le signe -- */
  if(*pointeur)
  {
       apres_moin_moin(pointeur);
  }
 
  if(interp == python) OUTBUFPTR+=sprintf(OUTBUFPTR,"]");
  if(interp == perl) OUTBUFPTR+=sprintf(OUTBUFPTR,"]");
  imprime(defo);
  printf("%s",OUTBUF);
/*   fprintf(stderr,"Number of characters on stdout=%d\n",strlen(OUTBUF));    */
  return(0);
}



/*********************************************************************
 *   initialiser les noms de clefs et les valeurs de defaut          *
 *********************************************************************/

char **ini_list(argv, defo, status)
char **argv;
struct definition defo[];
int *status;

/*   adaptation : j. caveen  janvier 1991 a partir du programme cclargs ecrit
*                                         par m. valin.
*
*    fonction servant a dresser la liste des noms de clefs ainsi que leurs
*    valeurs par defaut.  la fonction retourne un pointeur a l'element
*    de argv ou debute la portion de l'appel de l'usager.
*
*    les noms de clefs et leurs valeurs sont imprimes sur stdout.
*
*    arguments
*
*       argv        entree      pointeur au vecteur d'arguments
*       defo        sortie      tableau contenant le nom et les valeurs
*                               de defaut des clefs
*       status         "        compteur d'erreurs
*
*/

{
   void pas_de_deux(), converti();
   enum typecle majmin();
   int compte, i, ldesc;
   char *egal_pointeur ;

   i = -1;
   while(*argv)
   {
     if(**argv == '+')
     {
       if(*((*argv)+1) == '+')
       {
         argv++;
         break;
       }
      }
 
      else if(**argv == '-')
      {
        i++;
	NKLES=i;
        compte = 0;
        egal_pointeur = strchr((*argv)+1,'='); /* verifier pour format -kle= */
        if(egal_pointeur != (char *) NULL)
        {
           *egal_pointeur = '\0';
           if((*argv)+1 == '\0')
           {
               (*status)++;
               fprintf(stderr, "mauvais nom de clef\n");
               return NULL;
           }
           defo[i].type = majmin((*argv)+1);
           defo[i].kle_nom = (*argv)+1;
           pas_de_deux(egal_pointeur+1);
           converti(egal_pointeur+1,defo[i].type);
           defo[i].kle_def1 = defo[i].kle_val = egal_pointeur+1;
           compte++;
        }
        else
        {
           if((*argv)+1 == '\0')
           {
               fprintf(stderr, "mauvais nom de clef\n");
               (*status)++;
               return NULL;
           }
           defo[i].type = majmin((*argv)+1);
           defo[i].kle_nom = (*argv)+1;
        }
      }
      else if(**argv == '=')
      {
         if(compte == 0)
         {
           pas_de_deux((*argv)+1);
           converti((*argv)+1,defo[i].type);
           defo[i].kle_def1 = defo[i].kle_val = (*argv)+1;
           compte++;
         }
         else
         {
            pas_de_deux((*argv)+1);
            converti((*argv)+1,defo[i].type);
            defo[i].kle_def2 = (*argv)+1;
         }
            
      }
      else
      {
        if(**argv == '[')   /*  descripteur pour help   */
        {
           /* add extra blank to the end so there is no char trunc */
           ldesc = strlen(*argv) - 1 ;
           defo[i].kle_desc = malloc((ldesc+1)*sizeof(char));
           *((*argv)+ldesc) = ' ';
           ldesc++;  
           *((*argv)+ldesc) = '\0';
           strcpy(defo[i].kle_desc,(*argv)+1);
        }
        else
        {
             if(compte == 0)
             {
                converti(*argv,defo[i].type);
                defo[i].kle_def1 = defo[i].kle_val =  (*argv);
                compte++;
             }
             else
             {
                converti(*argv,defo[i].type);
                defo[i].kle_def2 = (*argv);
             }
        }
       
      }
      argv++;
    }
   i++;
   defo[i].kle_nom = '\0';
   return(argv); 
}




/**********************************************************
 *    traitement des arguments en mode positionnel        *
 **********************************************************/

char **positionel(argv)
char **argv;

/*  adaptation:  j. caveen janvier 1991 a partir du programme cclargs ecrit 
*                                       par m. valin
*
*   fonction servant a faire la recuperation des arguments en mode positionnel
*   la fonction retourne un pointeur a la position ou debute la recuperation
*   en mode equivalence
*
*   argument 
*
*       argv   entree  pointeur au mot suivant le ++ dans le vecteur d'arguments
*/

{
    int count = 0;
    
    if(interp == perl)
      OUTBUFPTR+=sprintf(OUTBUFPTR,"%cCclArgs=('--'=>[",'%');
    if(interp == python)
      OUTBUFPTR+=sprintf(OUTBUFPTR,"CclArgs={'--':[");
    if(interp == shell)
      OUTBUFPTR+=sprintf(OUTBUFPTR," set nil ; shift ; ");

    while(*argv)
    {
      if(**argv == '-')
      {
         break;
      }
      else
      {
         if (count == 0 && interp == shell)
         {
            OUTBUFPTR+=sprintf(OUTBUFPTR,"set ");
         }
         if(interp == shell) OUTBUFPTR+=sprintf(OUTBUFPTR," %s ",*argv);
         if(interp == python) OUTBUFPTR+=sprintf(OUTBUFPTR,"%s'%s'",separateur,*argv);
         if(interp == perl) OUTBUFPTR+=sprintf(OUTBUFPTR,"%s'%s'",separateur,*argv);
         separateur=virgule;
      }

      argv++;  count++;
     }

     if(count > 0) {
       if(interp == shell) OUTBUFPTR+=sprintf(OUTBUFPTR," ; ");
     }

     return(argv);
}



/*********************************************************************
 *       traitement des arguments en mode equivalence                *
 *********************************************************************/

char **equivalence(argv, defo, scriptnom,  status)
char **argv, *scriptnom;
struct definition defo[];
int *status;
/*
*   adaptation: j. caveen janvier 1991 a partir du programme cclargs 
*                                      ecrit par m. valin
*   fonction servant a la recuperation en mode equivalence des arguments 
*   de l'appel usager.  on verifie la validite du nom de chaque clef en 
*   consultant la liste de defaut defo[i].kle_nom.
*
*   arguments
*
*      argv      entree   pointeur a la premiere clef donnee par l'usager
*      defo      entree   liste des noms de clefs et de leur valeur par defaut
*      status    sortie   compteur d'erreurs
*
*/

{
   void pas_de_deux(), converti(), sequence_appel();
   int valide_kle();
   char **argu_list();
   int index;
   char *keyname;
   char *egal_pointeur;
   int nom_ecrit;
   char arg_val_buf[65538];

      while(*argv)
      {
        if(!strcmp(*argv,"--"))
        {
             argv++;
             return(argv);
        }

        if(!strcmp(*argv,"-+"))
        {
           on_affiche=1;
        }
        else
        {
        if(**argv == '-')
        {
           nom_ecrit = 0;
           egal_pointeur = strchr((*argv)+1,'=');
           if(egal_pointeur != (char *) NULL)
             *egal_pointeur = '\0';

           keyname = (*argv)+1;
           index = valide_kle(keyname,defo);
           if(index < 0)
           { 
              (*status)++;
              fprintf(stderr, "\n nom de clef pas defini, clef=%s \n",keyname);
              return(argv);
           }
           else
           {
              if(egal_pointeur != (char *) NULL)
              {
                 pas_de_deux(egal_pointeur+1);
                 converti(egal_pointeur+1,defo[index].type);
                 defo[index].kle_val = egal_pointeur+1; 
                 nom_ecrit = 1;
              }
              else
                 defo[index].kle_val = defo[index].kle_def2;
           }
        }
        else
        {
           if(! nom_ecrit)
           {
/*
           defo[index].kle_val = malloc(1000);
*/
           defo[index].kle_val = &arg_val_buf[0];
           arg_val_buf[0]='\0';
           argv = argu_list(argv,defo[index].kle_val,defo[index].type);
           defo[index].kle_val = malloc(strlen(arg_val_buf)+10);
           strcpy(defo[index].kle_val,arg_val_buf);
           }
        }
        }
        argv++;
      }

      return(argv);
}



/*********************************************************
 *             valider un nom de clef                    *
 *********************************************************/

int valide_kle(keyname, defo)
char *keyname;
struct definition defo[];

/*  auteur:  james caveen, janvier 1991
*
*   fonction s'assurant que le nom de clef utilise par l'usager est bel et
*   bien un nom valide (c.a.d. qui se retrouve dans la liste defo[i].kle_nom.)
*   la fonction retourne l'index i du nom de la clef dans la liste de defaut 
*   si le nom est valide, sinon la fonction retoure -1.
*
*    arguments
*
*        keyname        entree   nom de la clef a valider
*        defo             "      liste des nom de clefs possibles
*
*/
{
   
   int i = 0, index = -1;
  
   while(defo[i].kle_nom != '\0')
   {
      if(*defo[i].kle_nom == '_' ) {
        if(strcmp(keyname,defo[i].kle_nom+1) == 0)
             index = i; 
      } else {
        if(strcmp(keyname,defo[i].kle_nom) == 0)
             index = i; 
      }
      i++;
   }
   
   return(index);
}



/**************************************************************
 *     dresser une liste de valeurs de clefs                  *
 **************************************************************/

char **argu_list(argv,valeur,type)
char **argv;
char *valeur;
enum typecle type;

/*  dresser la liste de toutes les valeurs attribuees a une clef 
*   la fonction retourne un pointeur au dernier mot de la liste
*   des valeurs attribuees.
*
*   auteur: james caveen, janvier 1991
*
*   arguments:
*    
*       argv      "          pointeur aux valeurs a donner a kle
*
*/
{
    void converti(), pas_de_deux();


   if(**argv == '-')
   {
        argv--;
        return(argv);
   }

   converti(*argv,type);
   pas_de_deux(*argv);
   if(**argv == '=')
       valeur = strcat(valeur,(*argv)+1);
   else
       valeur = strcat(valeur,*argv);

    strcat( valeur ,"\0");
    argv++; 

    while(*argv)
    {
       if(**argv == '-')
            break;

          strcat( valeur ," \0");
          converti(*argv,type);
          pas_de_deux(*argv);
          if(**argv == '=')
              valeur = strcat(valeur,(*argv)+1);
          else
              valeur = strcat(valeur,*argv);


       strcat( valeur ,"\0");
       argv++; 
    }

    argv--;
    return(argv);
} 



/*****************************************************
 *       enlever les deux points d'une chaine        *
 *****************************************************/

void pas_de_deux(arg)
char *arg;

{
        while(*arg)
        {
           if(*arg == delimiter)
                *arg = ' ';

           arg++; 
        }
}



/********************************************************
 *          determiner le type de clef                  *
 ********************************************************/

enum typecle majmin(arg)
char *arg;

/*     fonction servant a determiner si la valeur a donner a une
*      clef est de type majuscule, minuscule ou si elle reste telle
*      que decrite lors de l'appel.  la fonction retourne le valeur
*      du type de clef (majus minus ou pareil) et met le dernier
*      element de arg egal a nul si necessaire.
*/
{

    int lng;
    lng = strlen(arg) - 1;

    if(*(arg+lng) == '+')
    { 
           *(arg+lng) = '\0';
           return(majus);
    }
    if(*(arg+lng) == '_')
    {
           *(arg+lng) = '\0';
           return(minus);
    }
    return(pareil);

}



/**********************************************************
 *    convertir une valeur de clef selon le type          *
 **********************************************************/

void converti(arg,type)
enum typecle type;
char *arg;

{
/*       fonction servant a faire la conversion d'une valeur de clef
*        si defo[index].type = majus, on force la valeur a des majuscules
*        si defo[index].type = minus, on force la valeur a des minuscules
*        si defo[index].type = pareil, on retourne la clef tel quel
*/

      if(type == pareil)
           ;

      else if(type == majus)
      {
          while(*arg)
          {
              *arg = toupper(*arg);
              arg++;
          }
      }
      else
      {
          while(*arg)
          {
              *arg = tolower(*arg);
              arg++;
          }
      }
}



/******************************************************
 *       imprimer sur stderr la sequence d'appel      *
 ******************************************************/

void sequence_appel(defo,scriptnom)
struct definition defo[];
char *scriptnom;

/*   fonction servant a imprimer sur le fichier stderr la liste des
*    noms de clefs et leurs valeurs de defaut
*
*  argument
*      
*         defo    -    structure contenant les noms de clefs ainsi que leurs
*                      valeurs de defaut et le type.
*/

{

     char *desc;
     int i = 0;

     
     fprintf(stderr,"\n *** SEQUENCE D'APPEL ***\n\n");

     fprintf(stderr,"%s [positionnels]\n",scriptnom);
   
     while(defo[i].kle_nom != '\0')
     {
       desc = defo[i].kle_desc ? defo[i].kle_desc : "";
       if(*defo[i].kle_nom == '_') /* supprimer le _ au debut des cles de sortie */
        fprintf(stderr," IN/OUT   -%s [%s:%s] %s\n",defo[i].kle_nom+1,defo[i].kle_def1,defo[i].kle_def2, desc);
       else
        fprintf(stderr," IN       -%s [%s:%s] %s\n",defo[i].kle_nom,defo[i].kle_def1,defo[i].kle_def2, desc);
         i++;
     }

     fprintf(stderr,"          [-- positionnels]\n");
     fprintf(stderr,"\n");
}



/*************************************************************
 *      imprimer les valeurs finales des clefs               *
 *************************************************************/

void imprime(defo)
struct definition defo[];

{
   int index=0;
   char *temp;

   if(interp == shell) {
     while(defo[index].kle_nom)
     {
	 temp = defo[index].kle_nom;
	 while(*temp == '-') temp++;
         if(strcmp(defo[index].kle_val,"%%%%"))
	    OUTBUFPTR+=sprintf(OUTBUFPTR," %s=\"%s\" ;",temp,defo[index].kle_val);
/*
	    OUTBUFPTR+=sprintf(OUTBUFPTR," %s=\"%s\" ;",defo[index].kle_nom,defo[index].kle_val);
*/
         index++;
     }
   }
   if(interp == perl) {
     while(defo[index].kle_nom)
	 temp = defo[index].kle_nom;
	 if(*temp == '-') temp++;
     {
         if(strcmp(defo[index].kle_val,"%%%%"))
	    OUTBUFPTR+=sprintf(OUTBUFPTR,",'%s'=>'%s'",temp,defo[index].kle_val);
/*
	    OUTBUFPTR+=sprintf(OUTBUFPTR,",'%s'=>'%s'",defo[index].kle_nom,defo[index].kle_val);
*/
         index++;
     }
     OUTBUFPTR+=sprintf(OUTBUFPTR,")");
   }
   if(interp == python) {
     while(defo[index].kle_nom)
     {
	 temp = defo[index].kle_nom;
	 if(*temp == '-') temp++;
         if(strcmp(defo[index].kle_val,"%%%%"))
	    OUTBUFPTR+=sprintf(OUTBUFPTR,",'%s':'%s'",temp,defo[index].kle_val);
/*
	    OUTBUFPTR+=sprintf(OUTBUFPTR,",'%s':'%s'",defo[index].kle_nom,defo[index].kle_val);
*/
         index++;
     }
     OUTBUFPTR+=sprintf(OUTBUFPTR,"}");
   }

}



/**********************************************************
 *    traitement des arguments  apres le signe --         *
 **********************************************************/

void apres_moin_moin(argv)
char **argv;

/*  adaptation:  j. caveen janvier 1991 a partir du programme cclargs ecrit 
*                                       par m. valin
*
*   fonction servant a faire la recuperation des arguments en mode positionnel
*   que l'on a ajoute apres le symbole -- lors de l'appel.
*
*   argument 
*
*     argv   entree  pointeur au mot suivant le ++ dans le vecteur d'arguments
*/

{
    
    if(interp == shell) OUTBUFPTR+=sprintf(OUTBUFPTR," set -- $* ");

    while(*argv)
    {
      if(interp == shell) OUTBUFPTR+=sprintf(OUTBUFPTR," %s ",*argv);
      if(interp == python) OUTBUFPTR+=sprintf(OUTBUFPTR,"%s'%s'",separateur,*argv);
      if(interp == perl) OUTBUFPTR+=sprintf(OUTBUFPTR,"%s'%s'",separateur,*argv);
      separateur=virgule;
      argv++;
    }

     if(interp == shell) OUTBUFPTR+=sprintf(OUTBUFPTR," ; ");

}


/*********************************************************************
 *   fonction servant a passer en mode interactif via vmenu          *
 *********************************************************************/
void interactif(defo,scriptnom,help_general)
struct definition defo[];
char *scriptnom, *help_general;
{

     void converti();
     int vmenu(), obligatoire();
     char **cle, **val, **aide;
     int i,nbliste,status, nfois = 0;

     /*   trouver le nombre de clefs total     */

     for(nbliste=0;(nbliste <= NKLEMAX && defo[nbliste].kle_nom); nbliste++)
                       ;

    nbliste++;
    /* allouer la memoire et initialiser les vecteurs pour vmenu  */

    cle = (char **) malloc(sizeof(char *)*nbliste);
/*    cle[0] = (char *) malloc(20*sizeof(char)); */
    for(i=1;i<nbliste;i++)
    {
      cle[i] = (char *) malloc(20*sizeof(char));
      strcpy(cle[i],defo[i-1].kle_nom);
    }

    val = (char **) malloc(sizeof(char *)*nbliste);
/*    val[0] = (char *) malloc(256*sizeof(char)); */
    for(i=1;i<nbliste;i++)
    {
      val[i] = (char *) malloc(256*sizeof(char));
      defo[i-1].kle_val && strcpy(val[i],defo[i-1].kle_val);
    }

    aide = (char **) malloc(sizeof(char *)*nbliste);
    aide[0] = (char *) malloc(182*sizeof(char));
    help_general && strcpy(aide[0],help_general);
    for(i=1;i<nbliste;i++)
    {
      aide[i] = (char *) malloc(182*sizeof(char));
      defo[i-1].kle_desc && strcpy(aide[i],defo[i-1].kle_desc);
    }

    do
    {
         status =  vmenu(scriptnom,cle,val,aide,nbliste,nfois);

         if(status == -1)
         {
                 printf(" exit %d ",status);
                 exit (status);
    
         }
/*
 *       injecter les valeurs finales dans la structure  
 */

         for(i=1;i<nbliste;i++)
         {
            converti(val[i],defo[i-1].type);
            defo[i-1].kle_val = val[i];
         }
         nfois++;
      }
      while(obligatoire(defo));   
}


/*******************************************************************
 *   obtenir le nom du script faisant l'appel a cclargs            *
 *******************************************************************/
void getnom(scriptnom,argv,size)
char *scriptnom, *argv;
int size;
{
     char  tmp[50], *pttmp = tmp;
     int lng, i, ncar=0;

     lng = strlen(argv);

/*
     extraction du nom du script a partir de la fin jusqu'a un nom 
     de chemin (i.e.  de "/bin/unscript", on extrait "unscript". 
*/
     for (i=lng-1; i>=0 && *(argv+i) != '/'; i--)
     {
        *pttmp++ = *(argv+i);
        ncar++;
     }
/*
   recopier a l'en-droit le nom du script dans scriptnom
*/
   
   if(ncar <= size)
   {
      *pttmp--;
      while(ncar)
      {
         *scriptnom++ = *pttmp--;
         ncar--;
      }
      *scriptnom = '\0';
   }
}


/************************************************************
 *    fonction servant a verifier si toutes les clefs       *
 *           obligatoires on obtenu une valeur              *
 ************************************************************/

int obligatoire(defo)
struct definition defo[];
{

      int result = 0, i = 0;
      static char *virgul = {",,,,\0"};

      while(defo[i].kle_nom)
      {
         strcmp(defo[i].kle_val,virgul) ? result : result++;    
         i++;
      }
 
      return (result);
}

                
