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

/*
 * CoMpIlAtIoN_OpTiOnS ::SX4=-O overlap::SX5=-O overlap::
 * a l'oppose de c_baseio.c
 */
#include <rpnmacros.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#ifdef WIN32    /*CHC/NRC*/
#include <fcntl.h>
#define S_IRUSR _S_IREAD
#define S_IWUSR _S_IWRITE
#define S_IRGRP _S_IREAD
#define S_IWGRP _S_IWRITE
#define S_IROTH _S_IREAD
#define S_IWOTH _S_IWRITE
typedef long int pid_t;
#else
#include <unistd.h>
#include <sys/types.h>
#if defined (__AIX__)
#include <fcntl.h>
#else
#include <sys/fcntl.h>
#endif
#endif
#include <sys/stat.h>

#if defined (NEC)
#endif
#include <stdlib.h>


/*  Pour le test testvm.f, on force la continuation 
 *  de l'execution malgre les erreurs rencontrees
 */

#if defined (TEST_VMM)
#     define SORTIR(X) /**/
#     define fd_err fdout
#else
#     define fd_err stderr
#     define SORTIR(X) f77name(tracebck)(); exit(X);
#endif

#if defined (NEC)
#     define memint INT_64
#else
#     define memint wordint
#endif

#define PRIVATE /**/

/*
 * Internal prototypes
 */
PRIVATE void ouvre_ou_ferme_controle(int , int , char *);
PRIVATE int obtient_environ();

#if defined (TEST_VMM)
#define NAMES Names
#define SLICES Slices
#define BLOCKS Blocks
#else
#define NAMES VmM__NaMeS
#define SLICES VmM__SlIcEs
#define BLOCKS VmM__BlOcKs
#endif
/*  definition du marqueur de blocs et de la macro d'initialisation INTBAR */
#define BARVAL 0xFFFA5A5A
#if defined (ALL64)
#    define INTBAR(x) *(BLOCKS[x].memadr+BLOCKS[x].size-1)=BARVAL
#else
#    define INTBAR(x) *(BLOCKS[x].memadr+BLOCKS[x].size-2)=(*(BLOCKS[x].memadr+BLOCKS[x].size-1))=BARVAL
#endif

#define SEGMENT_MASK 65535 /* Masques pour le tableau de blocs pour */
#define SEGMENT_SHIFT 16  /* le chargement par segment */
#define WEIGHT_MASK 31
#define WEIGHT_SHIFT 11
#define BKNO_MASK 2047

#define MAXSLICES 16384
#define MAXBLOCKS  2048
#define MAXNAMES    512
#define NCLASSE       9
#define NBITKEYMIN   12   /* nombre de bits pour la portion mineure de la clef */
#define NBITKEYMAJ   20   /* nombre de bits pour la portion majeure de la clef */
#define NBITKEYPAD   32
#define NCARATTR    100   /*nombre de caracteres maximal pour attributs */
#define NCARMAX     256   /*nombre de caracteres maximal pour chaine environnement */
#define NCLEMAX 1000
#define ENVDEBUG      1   /*Valeur de debug_mode dans VMM_CONFIG*/
#define ENVCHECKSUM   2   /*Valeur de checksum_mode dans VMM_CONFIG*/
#define ENVPARANOID  10   /*Valeur de mode paranoiaque dans VMM_CONFIG*/

#define UNKNOWNVAR             100  /*nom de variable pas defini */
#define BADINKEY               101  /*mauvaise valeur de KEY*/
#define NOT_IN_CORE            102 
#define ALREADY_LOCKED         103
#define NO_SPACE_LEFT_FOR_LOAD 104
#define NO_CALL_TO_VMMALLC     105
#define VMMALLC_ALREADY_CALLED 106
#define CONTROL_FILE_ERROR     107 /* erreur d'ouverture de VMM_01 ... */
#define PWD_ALREADY_SET        108
#define BAD_PASSWORD           109
#define PASSWORD_IS_SET        110
#define VAR_MUST_EXIST         111
#define CONTROLE_DAMAGE        112
#define WAS_ALTERED_RELEASE    113
#define BAD_CKSUM_MODE         114
#define NOT_ENOUGH_MEMORY      115
#define TOO_MANY_KEYS          116
#define BLOCK_DAMAGE           117
#define ATTRIBUTS_MODIFIES     118
#define SLICE_TOO_BIG          119
#define CHECKSUM_ERROR         120
#define CHECKSUM_READ_ERROR    121
#define KEEP_IN_CORE_CPK       122
#define KEEP_IN_CORE_CKMX      123
#define SPACE_STILL_ALLOCATED  124

/* callc macro servant a verifier si vmmallc a ete appele */
#define callc(X) (called_vmmallc == 0) ? vmmerr(X,NO_CALL_TO_VMMALLC) : 1

#if defined (ALL64)
typedef struct {
   unsigned int key_pading : NBITKEYPAD; /* remplissage 32 bits de gauche*/
   unsigned int key_major : NBITKEYMAJ;  /*portion majeure de la clef ( meme valeur que dans names*/
   unsigned int key_minor : NBITKEYMIN;  /*portion mineure de la clef ( ajoute a la portion majeure */
   } key_in_pieces ;
#else
#if !defined(Little_Endian)
typedef struct {
   unsigned int key_major : NBITKEYMAJ;   /*portion majeure de la clef ( meme valeur que dans names*/
   unsigned int key_minor : NBITKEYMIN;   /*portion mineure de la clef ( ajoute a la portion majeure */
   } key_in_pieces ;
#else
typedef struct {
   unsigned int key_minor : NBITKEYMIN;   /*portion mineure de la clef ( ajoute a la portion majeure */
   unsigned int key_major : NBITKEYMAJ;   /*portion majeure de la clef ( meme valeur que dans names*/
   } key_in_pieces ;
#endif
#endif
typedef union
   {
     wordint clef;
     key_in_pieces key;
   } complete_key;

#if !defined(Little_Endian)
typedef struct {
   unsigned int keep_in_core : 1;
   unsigned int is_in_core : 1;
   unsigned int in_used : 1;
   unsigned int locked : 1;
   unsigned int save : 1;
   unsigned int altered : 1;
   unsigned int was_altered : 1;
   unsigned int traced : 1;
   unsigned int hpa_alloc : 1;
   unsigned int disk_image : 1;
   unsigned int size8  : 1;
   unsigned int must_exist  : 1;
   unsigned int class : 4;
   unsigned int init   : 2;    /* initialisation: 0=aucune, 1= a zero, 2= au plus gros nombre machine */
   unsigned int weight : 4;
   unsigned int do_checksum : 1;
   unsigned int filling1 : 1;
   unsigned int filling2: 8;
   } BITFLAGS;
#else
typedef struct {
   unsigned int filling2: 8;
   unsigned int filling1 : 1;
   unsigned int do_checksum : 1;
   unsigned int weight : 4;
   unsigned int init   : 2;    /* initialisation: 0=aucune, 1= a zero, 2= au plus gros nombre machine */
   unsigned int class : 4;
   unsigned int must_exist  : 1;
   unsigned int size8  : 1;
   unsigned int disk_image : 1;
   unsigned int hpa_alloc : 1;
   unsigned int traced : 1;
   unsigned int was_altered : 1;
   unsigned int altered : 1;
   unsigned int save : 1;
   unsigned int locked : 1;
   unsigned int in_used : 1;
   unsigned int is_in_core : 1;
   unsigned int keep_in_core : 1;
   } BITFLAGS;
#endif

struct slice_table {
   union {
      unsigned int attributs;
      BITFLAGS flags;
      } info;
   int block_table_index;
   int name_table_index;
   int checksum;
   };

struct block_table {
   wordint *memadr;
   union {
      unsigned int attributs;
      BITFLAGS flags;
      } info;
   int slice_table_index;
   int file_adr;
   int size;
   int prev_fb;
   int next_fb;
   };

struct name_table {
   int base_file_adr;
   int lslice;
   int nslice;
   int major_key;
   int class;
   char nom[9];
   };

struct slice_table SLICES[MAXSLICES];
struct block_table BLOCKS[MAXBLOCKS];
struct name_table NAMES[MAXNAMES];

static int maxmem = 0, free_space = 0, nbslices = 0, nbvar = 0, nbblocks = 0;
static wordint mot_de_passe = 0, pwd_set = 0;
static int called_vmmallc = 0; /*called_vmmallc mis a  1 lors du 1er appel a vmmallc */
#if defined (_FLOAT1)
static wordfloat MAXVAL= 1.0e75;
static wordfloat zero=0.0;
static double MAXVAL8= 1.0e75;
static double zero8=0.0;
#else
static float MAXVAL= 1.0e38;
static float zero=0.0;
static double MAXVAL8= 1.0e308;
static double zero8=0.0;
#endif
 
/* variables globales pour les unites logiques des fichiers Vmm_classe
   et Vmm_controle
 */
static int fclass[NCLASSE] = { 0, 0, 0, 0, 0, 0, 0, 0, 0} ;
static char *fclass_names[NCLASSE] , *fcontrole_name ;
static wordint fcontrole=0;   /* nom du fichier: Vmm_controle */
static int wp_Vmm[NCLASSE] = { 0, 0, 0, 0, 0, 0, 0, 0, 0} ; /*longueur des Vmm_0n*/
static int fichiers_ouverts = 0; 
static int champs_bloques = 0; 
static int reprise = 0;   /* mis a 1 si le fichier Vmm_controle est non vide */
static FILE *fdout;
static int debug_mode = 0;
static int checksum_mode = 0;
static int espace_requis_max = 0; 
static int champs_bloques_max = 0; 
static int nbblocks_max = 0; 
static int nb_appels_no_lock = 0;
static int nb_appels_lock = 0;
static int nb_ecritures = 0;
static int nb_lectures = 0;
static int first_free_bloc = 0;
static int tableau_eject[MAXBLOCKS]; /* tableau ordonne des blocs du segment a
                                        utiliser pour un load */
static char cd_repertoire[128]= "./";       /* repertoire de travail pour les fichiers de controle */
      


/***s/p calc_checksum
*
*objet(calc_checksum)
*            Calculer le checksum d'un bloc memoire
*
*auteur J. Caveen - mai 1994
*
*arguments
*         in  bkno - numero du block memoire
*
**/
int
calc_checksum(int bkno)
{
     extern wordint f77name(qvmcks)();

     wordint  *adresse, longueur;
     wordint mode = 1;
     int checks;

     adresse = (wordint *) BLOCKS[bkno].memadr;
     longueur = (wordint) BLOCKS[bkno].size;

     checks=(int) (f77name(qvmcks)(adresse,&longueur,&mode));

     fprintf(fdout,"Checksum block numero %d, variable %s, tranche %d = %d\n",
               bkno,NAMES[SLICES[BLOCKS[bkno].slice_table_index].name_table_index].nom,
               BLOCKS[bkno].slice_table_index -
               NAMES[SLICES[BLOCKS[bkno].slice_table_index].name_table_index].major_key + 1,
               checks);

     return(checks);
}



/***s/p collapse_blocks
*
*objet(collapse_blocks)
*     Regrouper deux blocs vides adjacents
*
*auteur M. Lepine  -  juillet 1993
*
*arguments
*     in   i      index du premier bloc
*     in   inext  index du deuxieme bloc
*
**/
PRIVATE void 
collapse_blocks(int i, int inext)
{
   void imprime();
   int verbar();
   int j,ier;
   
   BLOCKS[i].size += BLOCKS[inext].size;
   BLOCKS[i].next_fb =
        (BLOCKS[inext].next_fb == -1 ? -1: BLOCKS[inext].next_fb-1);
   
   for (j=0; j < nbslices; j++)
      if (SLICES[j].block_table_index > i)
	 SLICES[j].block_table_index --;
   for (j=inext; j < nbblocks-1; j++) {
     ier = verbar(j);
     ier = verbar(j+1);
     BLOCKS[j].info.attributs = BLOCKS[j+1].info.attributs;
     BLOCKS[j].slice_table_index = BLOCKS[j+1].slice_table_index;
     BLOCKS[j].memadr = BLOCKS[j+1].memadr;
     BLOCKS[j].file_adr = BLOCKS[j+1].file_adr;
     BLOCKS[j].size = BLOCKS[j+1].size < 0 ? 0 :BLOCKS[j+1].size;
     BLOCKS[j].prev_fb = BLOCKS[j+1].prev_fb == -1 ? -1:BLOCKS[j+1].prev_fb -1;
     BLOCKS[j].next_fb = BLOCKS[j+1].next_fb == -1 ? -1:BLOCKS[j+1].next_fb -1;
     }
   nbblocks--;

/*
   fprintf(fdout," COLLAPSE_BLOC\n");
   imprime();
*/
   }
 


/***s/p ecrit_bloc
*
*objet(ecrit_bloc)
*     Ecrire un block sur disque ou dans le XMU
*
*auteur J. Caveen  -  juillet 1993
*
*arguments
*     in   bkno      numero du bloc a ecrire
*     in   classe      classe de la variable (1 a 9)
*     in   memadresse  adresse memoire de la tranche
*     in   fileadresse adresse sur le fichier de la tranche
*     in   nmots    longueur en mots de la tranche
*
**/
PRIVATE void 
ecrit_bloc(int bkno,int classe,wordint *memadresse,
                        int fileadresse,int nmots)
{
      int verbar(),calc_checksum();
/*
      void  ouvre_ou_ferme_controle();
*/

      wordint iun, lfileadresse, lnmots, ier;

      ier = verbar(bkno);
      if( ! fichiers_ouverts) ouvre_ou_ferme_controle(1,0,"ecrit_bloc");

      iun = (wordint) fclass[classe - 1];
      lfileadresse = (wordint) fileadresse;
      lnmots = (wordint) nmots;
/*
#if ! defined (ALL64)
      lnmots *=
       (SLICES[BLOCKS[bkno].slice_table_index].info.flags.size8 == 1) ? 2:1;
#endif
*/
      if(SLICES[BLOCKS[bkno].slice_table_index].info.flags.do_checksum)
            SLICES[BLOCKS[bkno].slice_table_index].checksum =
                                             calc_checksum(bkno);

      f77name(wawrit)(&iun,memadresse,&lfileadresse,&lnmots);
      BLOCKS[bkno].info.flags.altered = 0;
      BLOCKS[bkno].info.flags.was_altered = 0;
      BLOCKS[bkno].info.flags.disk_image = 1;
      SLICES[BLOCKS[bkno].slice_table_index].info.flags.altered = 0;
      SLICES[BLOCKS[bkno].slice_table_index].info.flags.was_altered = 0;
      SLICES[BLOCKS[bkno].slice_table_index].info.flags.disk_image = 1;

      if ((SLICES[BLOCKS[bkno].slice_table_index].info.flags.traced) ||
                                                             debug_mode)
          fprintf(fdout,
          "VMM trace: ecriture dans le fichier Vmm_0%d de la variable %s tranche %d\n",classe,
           NAMES[SLICES[BLOCKS[bkno].slice_table_index].name_table_index].nom,
           BLOCKS[bkno].slice_table_index -
           NAMES[SLICES[BLOCKS[bkno].slice_table_index].name_table_index].major_key + 1);


      nb_ecritures++;
} 

/***s/p ecrit_vmm_controle
 *
 *objet(ecrit_vmm_controle)
*      ecrire les structures du systeme de gestion de memoire virtuelle
*      dans le fichier Vmm_controle.
*
*      le fichier est structure comme suit:
*         [nbvar,tous les names,nbslices,toutes les slices]
*
*auteur: J. Caveen- juin 1993
*
*argument
**/
PRIVATE void 
ecrit_vmm_controle()
{
/*
      void  ouvre_ou_ferme_controle();
*/
   
/* ETG ofset of type off_t */
   int nmots,i;
   off_t ofset;

   if( ! fichiers_ouverts) ouvre_ou_ferme_controle(1,0,"ecrit_vmm_controle");

   ofset = 0;
   ofset = lseek(fcontrole, ofset, SEEK_SET);
   nmots = write(fcontrole,&nbvar,4);
   
   nmots = sizeof(struct name_table)*nbvar;
   nmots = write(fcontrole,&NAMES[0],nmots);
   /*
    *    ecriture des tranches
    *    on met tous les marqueurs appropries a nul
    *    (keep_in_core, is_in_core, etc
    *    apres l'ecriture, ont leur redonne leur valeur originale
    */
   
   for (i = 0; i < nbblocks; i++)
   {
      if(BLOCKS[i].info.flags.in_used)
      {
	 SLICES[BLOCKS[i].slice_table_index].block_table_index = -1;
	 SLICES[BLOCKS[i].slice_table_index].info.flags.keep_in_core = 0;
	 SLICES[BLOCKS[i].slice_table_index].info.flags.is_in_core = 0;
         SLICES[BLOCKS[i].slice_table_index].info.flags.in_used = 0;
	 SLICES[BLOCKS[i].slice_table_index].info.flags.locked = 0;
	 SLICES[BLOCKS[i].slice_table_index].info.flags.altered = 0;
	 SLICES[BLOCKS[i].slice_table_index].info.flags.was_altered = 0;
	 SLICES[BLOCKS[i].slice_table_index].info.flags.traced = 0;
      }
   }
   nmots = write(fcontrole,&nbslices,4);
   nmots = sizeof(struct slice_table)*nbslices;
   nmots = write(fcontrole,&SLICES[0],nmots);

/*
 *    on redonne aux slices les attributs des blocks
 */
   for (i = 0; i < nbblocks; i++)
   {
      if(BLOCKS[i].info.flags.in_used)
      {
	 SLICES[BLOCKS[i].slice_table_index].block_table_index = i;
	 SLICES[BLOCKS[i].slice_table_index].info.flags.keep_in_core = 
	 BLOCKS[i].info.flags.keep_in_core;
	 SLICES[BLOCKS[i].slice_table_index].info.flags.is_in_core = 
	 BLOCKS[i].info.flags.is_in_core;
	 SLICES[BLOCKS[i].slice_table_index].info.flags.in_used = 
	 BLOCKS[i].info.flags.in_used;
	 SLICES[BLOCKS[i].slice_table_index].info.flags.locked = 
	 BLOCKS[i].info.flags.locked;
	 SLICES[BLOCKS[i].slice_table_index].info.flags.altered = 
	 BLOCKS[i].info.flags.altered;
	 SLICES[BLOCKS[i].slice_table_index].info.flags.was_altered = 
	 BLOCKS[i].info.flags.was_altered;
	 SLICES[BLOCKS[i].slice_table_index].info.flags.traced = 
	 BLOCKS[i].info.flags.traced;
      }
   }

}

/***s/p eject_block
*
*objet(eject_block)
*     Ejecter un bloc de la memoire pour recuperer l'espace
*
*auteur M. Lepine  -  juillet 1993
*
*arguments
*     in   bkno     numero du bloc a expulser
*     in   save     indique si on doit sauver le bloc sur disque selon
*                   les attributs du bloc (save, altered ...)
*
**/
PRIVATE int 
eject_block(int bkno,int save,int fait_checksum)

{
   void imprime();
   void ecrit_bloc(), reserve_disk_space();
   int verbar(), calc_checksum();
   int ind,ier,cks,indx;

    if ( ! BLOCKS[bkno].info.flags.in_used)
         return(0);
   ind = BLOCKS[bkno].slice_table_index;

   if(ind != -1)
   {
      if(fait_checksum && (SLICES[ind].info.flags.do_checksum || checksum_mode))
      {
         cks = calc_checksum(bkno);
         if(cks != SLICES[ind].checksum)
            return(vmmerr("EJECT_BLOCK",CHECKSUM_ERROR));

         SLICES[ind].checksum = 0;
      }
   }
   ier = verbar(bkno);
   if ((BLOCKS[bkno].info.flags.traced) || debug_mode)
      fprintf(fdout,"VMM trace: ejection du bloc %d variable %s tranche %d\n",
      bkno,NAMES[SLICES[ind].name_table_index].nom,
      ind - NAMES[SLICES[ind].name_table_index].major_key + 1);


   if (save) {
     if ((BLOCKS[bkno].info.flags.save) && ((BLOCKS[bkno].info.flags.altered)
      || (BLOCKS[bkno].info.flags.was_altered))) {
        if(NAMES[SLICES[ind].name_table_index].base_file_adr == -1)
             reserve_disk_space(bkno);
	ecrit_bloc(bkno,BLOCKS[bkno].info.flags.class,BLOCKS[bkno].memadr,
		    BLOCKS[bkno].file_adr,BLOCKS[bkno].size);
	}
     }
   BLOCKS[bkno].info.attributs = 0;
   BLOCKS[bkno].slice_table_index= -1;
   if(ind != -1)
   {
      SLICES[ind].info.flags.is_in_core = 0;
      SLICES[ind].info.flags.keep_in_core = 0;
      SLICES[ind].info.flags.in_used = 0;
      SLICES[ind].info.flags.locked = 0;
      SLICES[ind].info.flags.altered = 0;
      SLICES[ind].info.flags.was_altered = 0;
      SLICES[ind].block_table_index = -1;
   }
   /* rechainage des pointeurs de blocs libres */
   if(first_free_bloc > bkno)
   {
      BLOCKS[bkno].prev_fb = -1;
      BLOCKS[bkno].next_fb = first_free_bloc;
      BLOCKS[first_free_bloc].prev_fb = bkno;
      first_free_bloc = bkno;
   } 
   else
   {
       indx = first_free_bloc;
       while(1)
       {
           if((BLOCKS[indx].next_fb == -1) || 
              (BLOCKS[indx].next_fb > bkno))
                      break;
  
           indx = BLOCKS[indx].next_fb;
       }
       BLOCKS[bkno].next_fb = BLOCKS[indx].next_fb;
       BLOCKS[bkno].prev_fb = indx;
       BLOCKS[BLOCKS[indx].next_fb].prev_fb = bkno;
       BLOCKS[indx].next_fb = bkno;
    }

/*
   imprime();
*/
   return(BLOCKS[bkno].size);
   }


/***s/p eject_from_tableau - ejecter les blocs d'un segment
 *
 * Auteur James Caveen - mai 1994
 *
 * Objet (eject_from_tableau) - ejecter les blocs memoires d'un segment
 *                         a partir de la liste contenue dans le tableau
 *                         tableau_eject.  L'ejection se fait en tenant
 *                         compte du poids des blocs a ejecter.  On arrete
 *                         l'ejection des qu'il y a assez d'espace dans le 
 *                         segment 
 *
 * Arguments : in size                - espace total requis
 *                tableau_eject_index - indice de depart dans le tableau
 *
 **/
PRIVATE int 
 eject_from_tableau(int size,int tableau_eject_index)
{

   int i , espace_libre, bloc_index;

     i = tableau_eject_index;
     espace_libre = 0;
     while(espace_libre < size)
     {
          bloc_index = tableau_eject[i] & BKNO_MASK;
          if(! BLOCKS[bloc_index].info.flags.in_used) 
            espace_libre += BLOCKS[bloc_index].size;
          else
            espace_libre += eject_block(bloc_index,1,1);

	  i++;
     } /* end while */

   return(espace_libre);
}



/*** fichier_vide - verifier si un fichier contient quelque chose
 *
 * Auteur : J. Caveen - juin 1994
 *
 * Objet (fichier_vide) Verifier si un fichier contient des donnees.  On ouvre
 *                      le fichier en mode lecture et on lit un caractere.
 *                      Si on a atteint la fin du fichier, le fichier est vide.
 *                      La fonction retourne 1 si le fichier est vide sinon on
 *                      retourne 0.
 *
 * Argument   in   nom -  nom du fichier a verifier
 *
 **/
PRIVATE int
fichier_vide(char *nom)
{
   FILE *fd;

   int retour, iii;
   if((fd = fopen(nom,"r")) == (FILE *) NULL)
      return(1);


   iii = getc(fd);

   retour = feof(fd);

   fclose(fd);
   return(retour);
}


/***s/p imprime 
 *
 * Auteur James Caveen - mai 1994
 *
 * Objet (imprime) - imprimer les listes de blocs libres et utilises
 *
 **/
void imprime()
{
     int i;

     printf(" Nombre de BLOCKSs = %d\n",nbblocks);
     printf(" Premier BLOCKS libre = %d\n", first_free_bloc);

     printf(" Liste des BLOCKSs libres\n");

     i = first_free_bloc;
     while(i != -1)
     {
          printf(" BLOCKS[%d].prev=%d,BLOCKS[%d].next=%d,BLOCKS[%d].size=%d\n",
                 i, BLOCKS[i].prev_fb,i, BLOCKS[i].next_fb,i,BLOCKS[i].size);
          i = BLOCKS[i].next_fb;
     }

     printf(" Liste des BLOCKSs utilises\n");
     for (i=0; i < nbblocks; i++)
     {
         if(BLOCKS[i].info.flags.in_used)
           printf(" BLOCKS[%d].next=%d,BLOCKS[%d].prev=%d,BLOCKS[%d].size=%d\n",
                 i, BLOCKS[i].prev_fb,i, BLOCKS[i].next_fb,i,BLOCKS[i].size);
     }
}


/***s/p impval - imprimer les valeurs d'un champ autour du marqueur de bloc
 *
 * Auteur: J. Caveen - avril 1994
 * 
 * Revision: james caveen - mars 1995
 *            change la nature de l'argument de float a wordint
 *            afin d'eliminer les incompatibilites de types lors de l'appel.
 *
 *objet (impval)  - fonction servant a imprimer les valeurs reelles du
 *                  contenu des adresses memoires autour d'un delimiteur
 *                  de bloc.
 *argument    entree adresse - adresse memoire a partir de laquelle on imprime
 *
 *      
 **/
void
    impval(wordint *adresse)
{

   void imp_bar();
    int i;
    int intval;
    float *fadresse;

    fadresse = (float *) adresse;

    intval = BARVAL;

    for (i=0; i<5; i++)
    {
        fprintf(fdout,"%f ",*fadresse);
        fadresse++;
    }

    imp_bar(&intval);
    
}
/***s/p imp_bar - imprimer la valeur float du delimiteur de blocs
 *
 *Auteur: james caveen - avril 1994
 *
* Revision: james caveen - mars 1995
 *            change la nature de l'argument de float a int
 *            afin d'eliminer les incompatibilites de types lors de l'appel.
 *
 *objet (imp_bar) imprimer en format float la valeur des delimiteurs
 *                de blocs memoires
 *
 *argument  entree - valeur : valeur entiere a imprimer en format float
 *
 **/
void
  imp_bar(int *valeur)
{
   float *fvaleur;
   fvaleur = (float *) valeur;
    fprintf(fdout,"\nTRUE FLOAT VALUE OF BLOCK DELIMITOR: %f\n", *fvaleur);
}


/***s/p imprime_structures
 *
 *objet(imprime_structures)
*      imprimer le contenu de toutes les structures allouees jusqu'ici
*auteur: J. Caveen- juin 1993
*
*argument
*          in mode  - indique quel type de structures imprimer
*                     mode = 0  - imprime block_table
*                     mode = 1  - imprime slice_table
*                     mode = 2  - imprime name_table
*
*/
PRIVATE void 
imprime_structures(int mode)
{
       int indice;
/*
 *     impression du contenu de la table names
 */
       switch(mode)
       {
          case 2:
              printf("\nContenu de names\n");
              for (indice=0; indice<nbvar; indice++)
              {
                printf("  Indice de la variable: %d\n",indice);
                printf("     nom          : %s\n",NAMES[indice].nom);
                printf("     base_file_adr: %d\n",NAMES[indice].base_file_adr);
                printf("     lslice       : %d\n",NAMES[indice].lslice);
                printf("     nslice       : %d\n",NAMES[indice].nslice);
                printf("     major_key    : %d\n",NAMES[indice].major_key);
                printf("     class        : %d\n",NAMES[indice].class);
              }              
              break;
          case 1:
              printf("\nContenu de slices\n");
              for (indice=0; indice<nbslices; indice++)
              {
                printf("  Indice de la slice: %d\n",indice);
                printf("     keep_in_core       : %d\n",SLICES[indice].info.flags.keep_in_core);
                printf("     is_in_core         : %d\n",SLICES[indice].info.flags.is_in_core);
                printf("     in_used            : %d\n",SLICES[indice].info.flags.in_used);
                printf("     locked             : %d\n",SLICES[indice].info.flags.locked);
                printf("     save               : %d\n",SLICES[indice].info.flags.save);
                printf("     altered            : %d\n",SLICES[indice].info.flags.altered);
                printf("     was_altered        : %d\n",SLICES[indice].info.flags.was_altered);
                printf("     traced             : %d\n",SLICES[indice].info.flags.traced);
                printf("     hpa_alloc          : %d\n",SLICES[indice].info.flags.hpa_alloc);
                printf("     disk_image         : %d\n",SLICES[indice].info.flags.disk_image);
                printf("     size8              : %d\n",SLICES[indice].info.flags.size8);
                printf("     must_exist         : %d\n",SLICES[indice].info.flags.must_exist);
                printf("     class              : %d\n",SLICES[indice].info.flags.class);
                printf("     weight             : %d\n",SLICES[indice].info.flags.weight);
                printf("     do_checksum        : %d\n",SLICES[indice].info.flags.do_checksum);
                printf("     init               : %d\n",SLICES[indice].info.flags.init);
                printf("     block_table_index  : %d\n",SLICES[indice].block_table_index);
                printf("     name_table_index   : %d\n",SLICES[indice].name_table_index);
                printf("     checksum           : %d\n",SLICES[indice].checksum);
              }
              break; 
          case 0:
              printf("\nContenu de blocks\n");
              for (indice=0; indice<nbblocks; indice++)
              {
                printf("  Indice du bloc: %d\n",indice);
                printf("     keep_in_core       : %d\n",BLOCKS[indice].info.flags.keep_in_core);
                printf("     is_in_core         : %d\n",BLOCKS[indice].info.flags.is_in_core);
                printf("     in_used            : %d\n",BLOCKS[indice].info.flags.in_used);
                printf("     locked             : %d\n",BLOCKS[indice].info.flags.locked);
                printf("     save               : %d\n",BLOCKS[indice].info.flags.save);
                printf("     altered            : %d\n",BLOCKS[indice].info.flags.altered);
                printf("     was_altered        : %d\n",BLOCKS[indice].info.flags.was_altered);
                printf("     traced             : %d\n",BLOCKS[indice].info.flags.traced);
                printf("     hpa_alloc          : %d\n",BLOCKS[indice].info.flags.hpa_alloc);
                printf("     disk_image         : %d\n",BLOCKS[indice].info.flags.disk_image);
                printf("     size8              : %d\n",BLOCKS[indice].info.flags.size8);
                printf("     must_exist         : %d\n",BLOCKS[indice].info.flags.must_exist);
                printf("     class              : %d\n",BLOCKS[indice].info.flags.class);
                printf("     weight             : %d\n",BLOCKS[indice].info.flags.weight);
                printf("     do_checksum        : %d\n",BLOCKS[indice].info.flags.do_checksum);
                printf("     init               : %d\n",BLOCKS[indice].info.flags.init);
                printf("     slice_table_index  : %d\n",BLOCKS[indice].slice_table_index);
                printf("     file_adr           : %d\n",BLOCKS[indice].file_adr);
                printf("     memadr             : %x\n",BLOCKS[indice].memadr);
                printf("     size               : %d\n",BLOCKS[indice].size);
                printf("     prev_fb            : %d\n",BLOCKS[indice].prev_fb);
                printf("     next_fb            : %d\n",BLOCKS[indice].next_fb);
              }
              break;
       }
}

/***s/p lit_bloc
*
*objet(lit_bloc)
*     Lire  un block sur disque ou dans le XMU
*
*auteur J. Caveen  -  juillet 1993
*
*arguments
*     in   bkno        numero du bloc a lire
*     in   classe      classe de la variable (1 a 9)
*     in   memadresse  adresse memoire de la tranche
*     in   fileadresse adresse sur le fichier de la tranche
*     in   nmots       longueur en mots de la tranche
*
**/
PRIVATE void 
lit_bloc(int bkno,unsigned int classe,wordint *memadresse,
                        int fileadresse,int nmots)
{
/*
      void  ouvre_ou_ferme_controle();
*/
      int calc_checksum();
      wordint iun, lfileadresse, lnmots;

      int cks;

      if( ! fichiers_ouverts) ouvre_ou_ferme_controle(1,0,"lit_bloc");

      iun = (wordint) fclass[classe - 1];
      lfileadresse = (wordint) fileadresse;
      lnmots = (wordint) nmots;
/*
#if ! defined (ALL64)
      lnmots *=
       (SLICES[BLOCKS[bkno].slice_table_index].info.flags.size8 == 1) ? 2:1;
#endif
*/
      f77name(waread)(&iun,memadresse,&lfileadresse,&lnmots);

      if ((SLICES[BLOCKS[bkno].slice_table_index].info.flags.traced) ||
                                                             debug_mode)
          fprintf(fdout,
          "VMM trace: lecture dans le fichier Vmm_0%d de la variable %s tranche %d\n",
           classe,
           NAMES[SLICES[BLOCKS[bkno].slice_table_index].name_table_index].nom,
           BLOCKS[bkno].slice_table_index -
           NAMES[SLICES[BLOCKS[bkno].slice_table_index].name_table_index].major_key + 1);
      
/*
 *     faire le checksum si slice marquee comme tel
 */
       if((SLICES[BLOCKS[bkno].slice_table_index].checksum != 0) &&
          (SLICES[BLOCKS[bkno].slice_table_index].info.flags.do_checksum ||
           checksum_mode))
       {
             cks = calc_checksum(bkno);
             if(cks != SLICES[BLOCKS[bkno].slice_table_index].checksum)
                  vmmerr("LIT_BLOC",CHECKSUM_READ_ERROR);
       }
       else if (SLICES[BLOCKS[bkno].slice_table_index].info.flags.do_checksum ||
           checksum_mode)
       {
                 SLICES[BLOCKS[bkno].slice_table_index].checksum = 
                    calc_checksum(bkno);
       }

       nb_lectures++;
} 

/***s/p lit_vmm_controle
*objet(lit_vmm_controle) relire les structures du fichier Vmm_controle
*
*auteur:  J. Caveen - juillet 1993
*
*revision: J. Caveen - juin 1994
*                      Verifier si les fichiers de donnes sont disponibles 
*                      pour les variables de type MUSTEXIST qui ont ete 
*                      utilisees lors de tranches precedentes. (Une variable
*                      MUSTEXIST qui n'a jamais ete utilisee n'a pas besoin
*                      d'avoir une image disque).
*
*                      Ajout d'un appel a une fonction servant a verifier 
*                      si le fichier a lire est vide ou non.  Avant, on 
*                      utilisait la "propriete" des fichiers WA de 
*                      retourner un nombre de mots lu de zero lorsque le
*                      fichier est vide.
*                      
*
*
**/
PRIVATE void 
lit_vmm_controle()
{

   int f77name(vmmint)(), fichier_vide(), vmmerr();

/*ETG ofset of type off_t */
   int i,j,k,ier,nmots, erreur;
   off_t ofset;
   wordint lpos,liun,lun;

   ofset = 0;
   erreur = 0;
   lun = 1;

   if( fichier_vide(fcontrole_name))
   {
        nbvar = 0;
        nbslices = 0;
   }
   else
   {
       ofset = lseek(fcontrole, ofset, SEEK_SET);
       read(fcontrole,&nbvar,4);


       if((nbvar > 0) && (nbvar <= MAXNAMES))
       {
            nmots = sizeof(struct name_table)*nbvar;
            read(fcontrole,&NAMES[0],nmots);
       }
        
       read(fcontrole,&nbslices,4);
       if((nbslices > 0) && (nbslices <= MAXSLICES))
       {
            nmots = sizeof(struct slice_table)*nbslices;
            read(fcontrole,&SLICES[0],nmots);
       }

   }

   if(nbslices > 0 && nbvar > 0)
   { 
      reprise = 1 ;
/*
 *    s'assurer que les structures lues sont bonnes
 */
      ier = f77name(vmmint)();
/*
 *    initialiser les adresses du premier mot ou ecrire a la fin des fichiers
 */
      for (i = 0; i < NCLASSE; i++)
      {
          liun = (wordint) fclass[i];
          lpos = 0;
          if(fichier_vide(fclass_names[i]))
              wp_Vmm[i] = 0;
          else
          {
              f77name(waread)(&liun,&lpos,&lun,&lun);
              wp_Vmm[i] = (int) lpos;
          }
      }
/*
*     verifier si les variables de type mustexist qui ont deja
*     ete utilisees ont bel et bien un fichier contenant leur
*     image disque.
*/
      for ( j = 0; j < nbvar; j++)
      {
          i = NAMES[j].class - 1;
          if(wp_Vmm[i] <= 2 )
          {
             if(NAMES[j].base_file_adr != -1 &&
                       SLICES[NAMES[j].major_key].info.flags.must_exist)
              {
                 fprintf(fd_err," Variable %s must exist for a restart\n and file Vmm_0%d is absent\n",NAMES[j].nom,NAMES[j].class);
                 erreur++;
              }
              NAMES[j].base_file_adr = -1;
 
              for(k = 0; k < NAMES[j].nslice; k++)
                      SLICES[NAMES[j].major_key+k].info.flags.disk_image=0;
           }
       }
     }

     if(erreur)
         erreur = vmmerr("lit_vmm_controle",VAR_MUST_EXIST);
/*
 *   initialiser la position d'ecriture pour les fichiers de controle
 */
     for (i = 0; i < NCLASSE; i++)
           wp_Vmm[i] = (wp_Vmm[i] < 2 ? 2 : wp_Vmm[i]);
}

/***s/p obtient_environ
*
*objet(obtient_environ)
*           Obtenir la valeur de VMM_CONFIG et initialiser les variables
*           debug_mode,checksum_mode ainsi que le tableau fclass_names[] et fcontrole
*
*auteur J. Caveen  -  aout 1993
*
**/
/**
*
*   revision:  mai 1994 - ajout de checksum_mode
*       revision: E. Gondet - 19 avril 2002
*                 Repertoire pour decouplage en parallele des fichiers
*		  Rangement de tous les fichiers de travail dans ce repertoire.
*   003 - Avril 2003 - M. Lepine, initialisation dbg_cks 
*
*
**/
PRIVATE int obtient_environ()
{

      char repertoire[NCARMAX], rslt_fic[NCARMAX];
      char spid[NCARMAX];      /* the pid under char* form */
      char *ptenv;
      int i, longueur ,dbg_cks=0;
      int il_err;
      pid_t pid;

      static char *noms[] = {"Vmm_01","Vmm_02","Vmm_03","Vmm_04","Vmm_05",
                        "Vmm_06","Vmm_07","Vmm_08","Vmm_09","Vmm_controle"} ;

/*
 *    initialiser repertoire a des nuls
 */
      for ( i = 0; i < NCARMAX; i++)
      {
        repertoire[i] = '\0';
        rslt_fic[i] = '\0';
      }
/*ETG 17/04/2002 rajout
*/
      il_err = sprintf(repertoire,"%s", cd_repertoire);

      ptenv = (char *) getenv("VMM_CONFIG");
     
      if(ptenv != (char *) NULL)
      {
        sscanf(ptenv,"%s %d %s",&repertoire[0],&dbg_cks,&rslt_fic[0]);
      }
/*
 *    initialiser debug_mode et/ou checksum_mode selon la valeur de dbg_cks
 */

      switch(dbg_cks)
      {
         case ENVDEBUG:
              debug_mode = 1;
              break;
         case ENVCHECKSUM:
              checksum_mode = 1;
              break;
         case ENVPARANOID:
              debug_mode = 1;
              checksum_mode = 1;
              break;
       }

/* 
 *    initialiser les noms de fichiers Vmm_0n et Vmm_controle
 */

      longueur = strlen(repertoire);
      if (longueur > 0)
      {
          if(repertoire[longueur-1] != '/')
          {
               repertoire[longueur] = '/';
               longueur++;
          }
      }


/*ETG 10/04/01 on attrappe le pid                                               
      pid = getpid();
      sprintf(spid,"%i", pid);
*/


      for (i = 0; i < 9; i++)
      {
           fclass_names[i] = (char *) calloc(longueur+7,sizeof(char));
           strcpy(fclass_names[i],repertoire);
           strcat(fclass_names[i],noms[i]);                                     
/*ETG  fclass_names[i] need strlen(spid) and +1 for the . more characters       */
/*
           fclass_names[i] = (char *) calloc(longueur+7+strlen(spid)+1,sizeof(char));
	   sprintf(fclass_names[i],"%s%s.%s",repertoire, noms[i], spid);
*/
      }

      fcontrole_name = (char *) calloc(longueur+13,sizeof(char));   
      strcpy(fcontrole_name,repertoire);
      strcat(fcontrole_name,noms[9]);                                          
/*ETG  fclass_names[i] need strlen(spid) and +1 for the . more characters       */
/*
      fcontrole_name = (char *) calloc(longueur+13+strlen(spid)+1,sizeof(char));
      sprintf(fcontrole_name,"%s%s.%s",repertoire, noms[9], spid);
*/

/*    
 *    ouvrir le fichier pour imprimer les resultats d'execution
 */
      if(strlen(rslt_fic) > 0)
      {
          if( ! strncmp(rslt_fic,"stdout",6))
             fdout = stdout;
          else if( ! strncmp(rslt_fic,"fd_err",6))
             fdout = fd_err;
          else
          {
            if((fdout = fopen(rslt_fic,"w")) == (FILE *) NULL)
            {
               fprintf(fd_err," WARNING - CANNOT OPEN OUTPUT FILE %s\n",rslt_fic);
               fprintf(fd_err,"           USING STDOUT  INSTEAD\n");
               fdout = stdout ;
            }
          }
       }
       else
       {
             fdout = stdout;
       }

      if(debug_mode)
      {
         fprintf(fdout," VMM_CONFIG=%s\n",ptenv);
         fprintf(fdout," Repertoire pour fichiers de controle=%s\n",repertoire);
         fprintf(fdout," Fichier de sortie=%s\n",rslt_fic);
      }

      return(longueur);
}

/***s/p ouvre_ou_ferme_controle
*
*objet(ouvre_ou_ferme_controle)
*     Ouvrir les fichiers necessaire au systeme de gestion de memoire virtuelle
*
*auteur J. Caveen  -  juillet 1993
*       revision: E. Gondet - 19 avril 2002
*		  Porting on FUJITSU (__uxpv__)
*
*
*argument 
*           in ouvre           si ouvre != 0, on ouvre
*                              si ouvre  = 0, on ferme
*              premiere_fois   si premiere_fois, on appelle fnom pour faire
*                              l'association unite logique-nom de fichier
*              fonction        nom de la fonction ayant fait l'appel
**/
PRIVATE void 
ouvre_ou_ferme_controle(int ouvre, int premiere_fois, char *fonction)
{

         extern void c_waopen(), c_waclos();
         int  vmmerr();
	 /*
         int  obtient_environ();
	 */
         
         
         int i, lng,  mzero = 0;
         wordint  iun, ier;
         char *ouvmod = "RND+R/W";
/*
 *   au premier appel, on fait l'association fichiers-numero d'unite
 */
     if(premiere_fois)
     {
/*
 *        on obtient la valeur de VMM_CONFIG
 *        et on compose les noms de fichiers a ouvrir
 */
/*ETG obtient_environ renvoit la taille du path uniquement                           */
          lng = obtient_environ();
          ier = 0;
          for ( i = 0; i < NCLASSE; i++)
           {
              iun = 0;
              ier += f77name(fnom)(&iun,fclass_names[i],ouvmod,&mzero,lng+6,7); 
/*ETG BUG There is a bug with the following line chksum become false ETG*/
/*ETG              ier += f77name(fnom)(&iun,fclass_names[i],ouvmod,&mzero,lng+strlen(fclass_names[i]),7); */
              fclass[i] = (int) iun;
           }
          if(ier != 0)
         vmmerr(fonction,CONTROL_FILE_ERROR);
     }

     if(ouvre)
     {
         for (i = 0; i < NCLASSE; i++)
                  c_waopen(fclass[i]);

         fcontrole = open(fcontrole_name,O_RDWR | O_CREAT,
             S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
         fichiers_ouverts = 1;
     }
     else
     {
         for (i = 0; i < NCLASSE; i++)
                  c_waclos(fclass[i]);

         close(fcontrole);
         fichiers_ouverts = 0;
     }
}



/***s/p pack_blocks
*
*objet(pack_blocks)
*     Reorganisation des blocs memoire, recuperation des espaces vides
*     et retourne la dimension du plus gros bloc libre
*
*auteur M. Lepine  -  juillet 1993
*
*argument
*   out   biggest_free_block_index  indice du plus gros bloc libre
*
**/
PRIVATE int 
pack_blocks(int *biggest_free_block_index)
{
   void swap_blocks(),collapse_blocks();
   int i=0, big=0, ind = -1;
   
   while ((i < nbblocks-1) && (! BLOCKS[i].info.flags.hpa_alloc)) {
     if (BLOCKS[i].info.flags.in_used)
	i++;
     else {
        while ((i<nbblocks-1) && (! BLOCKS[i+1].info.flags.in_used))
	   collapse_blocks(i,i+1);
	if (i < nbblocks)
	   if ((BLOCKS[i].info.flags.locked) || (BLOCKS[i+1].info.flags.locked))
	      i++;
	   else
	      while ((i<nbblocks-1) && (BLOCKS[i+1].info.flags.in_used) &&
		     (! BLOCKS[i+1].info.flags.locked)) {
	         swap_blocks(i,i+1);
		 i++;
		 }
	} /* end else */
     } /* end while */
   for (i=0; i < nbblocks; i++)
      if (! BLOCKS[i].info.flags.in_used)
	 if (BLOCKS[i].size > big) {
	    big = BLOCKS[i].size;
	    ind = i;
	    }
   *biggest_free_block_index = ind;
   return(big);
   }
            

/***s/p pack_segment
*
*objet(pack_segment)
*     Reorganisation des blocs memoire, recuperation des espaces vides
*     et retourne la dimension du plus gros bloc libre du segment debutant
*     au bloc bkno
*
*auteur M. Lepine  -  juillet 1993
*Modification: J.Caveen - mai 1994 
*               adaptation du pack_bloc de M.Lepine pour ne
*               traiter que les blocs a l'interieur d'un segment
*
*argument
*   in    bkno numero du bloc marquant le debut du segment
*   out   biggest_free_block_index  indice du plus gros bloc libre
*
**/
PRIVATE int 
pack_segment(int bkno, int *biggest_free_block_index)
{
   void swap_blocks(),collapse_blocks();
   int  big=0, ind = -1, index;
   int i;
   /*
    * On avance au premier bloc non bloque
    */
   i = bkno;
   while ((i < nbblocks-1) && (BLOCKS[i].info.flags.locked))
              i++;

   index = i;

   while ((i < nbblocks-1) && ((! BLOCKS[i].info.flags.hpa_alloc)||
                               (! BLOCKS[i].info.flags.locked))) 
   {
     if (BLOCKS[i].info.flags.in_used)
	i++;
     else 
     {
        while ((i<nbblocks-1) && (! BLOCKS[i+1].info.flags.in_used))
	   collapse_blocks(i,i+1);

	if (i < nbblocks)
        {
	   if ((BLOCKS[i].info.flags.locked) || (BLOCKS[i+1].info.flags.locked))
	       break;
	   else
	      while ((i<nbblocks-1) && (BLOCKS[i+1].info.flags.in_used) &&
		     (! BLOCKS[i+1].info.flags.locked)) 
              {
	         swap_blocks(i,i+1);
		 i++;
              }
        }
     } /* end else */
   } /* end while */

    while ((! BLOCKS[index].info.flags.locked) && ( index < nbblocks))
    {
      if (! BLOCKS[index].info.flags.in_used)
      {
	 if (BLOCKS[index].size > big) 
         {
	    big = BLOCKS[index].size;
	    ind = index;
	 }
      }
      index++;
    }

   *biggest_free_block_index = ind;
   return(big);
   }


/***s/p qvmindex_from_key
*
*objet(qvmindex_from_key)
*    Valider et retourner l'indice dans la table slices a partir d'une clef de l'usager
*
*auteur J. Caveen  -  juin 1993
*
*argument
*     in   inkey   clef usager pointant a une des slices de la variable
*
**/
PRIVATE int 
qvmindex_from_key(complete_key inkey)

{

    int  laclef ;


    laclef = inkey.key.key_major + (inkey.key.key_minor == 0 ? 0: inkey.key.key_minor-1); 

    if(laclef > nbslices || laclef < 0)
         return(-BADINKEY);

    return ( inkey.key.key_major == NAMES[SLICES[laclef].name_table_index].major_key ? laclef : -BADINKEY);
}
                 

/***s/p qvmlod
*
*objet(qvmlod)
*      Charger un champ en memoire
*
*auteur M. Lepine  -  juillet 1993
*
*revisions J.Caveen -  octobre 1993
*                  vmmlod devient qvmlod qui est appele par le
*                  nouveau vmmlod. Si il y a des champs bloques
*                  les champs sont charges un a la fois, sinon on
*                  les charge d'un coup.
*
*       J. Caveen  -  mars 1994
*                  changement de l'algorithme de chargement:  pour charger
*                  un bloc, on fait les operations dans l'ordre suivant,
*                  si necessaire:
*                        1 - cherche un best fit, 
*                        2 - eject jusqu'a espace_requis <= espace_libre
*                            et pack_blocks
*                        3 - eject tous les blocs ejectables et pack_blocks
*                        4 - si espace_requis > espace_libre abort
*
*                  Cette modif vise a eliminer un probleme qui survenait
*                  lorsque l'espace disponible etait suffisant des le debut
*                  du chargement.  On ne faisait alors jamais d'ejection par
*                  la suite.
*
*        
*       J. Caveen - juin,juillet 1994
*                  Introduction de la methode de chargement par segments
*                  et de l'utilisation de listes liees de blocs vides.
*                  Un segment est une suite de champs non bloques compris
*                  entre deux champs bloques.
*                  Nouvelle methode de chargement:
*
*                    1 - parcourir la liste des blocs vides pour trouver
*                        un best_fit ( le plus petit bloc assez grand pour
*                        contenir le champ a charger).
*                    2 - si il n'y a pas de best-fit, trouver le segment 
*                        le moins couteux a utiliser ( c.a.d. celui duquel
*                        on ejectera le moins de champs possible du plus
*                        petit poids possible), ejecter des champs du segment
*                        en tenant compte des poids, jusqu'a l'obtention de
*                        l'espace et comprimer le segment.
*
*                    La fonction trouve_best_free() se charge de gerer tout le
*                    travail et de retourner l'index du bloc qui sera utilise
*                    pour le chargement du champ.
*              
*arguments
*         in inlkey  -  liste de clefs des champs a charger
*         in nkey    -  nombre de clefs dans inlkey
*
**/
wordint
qvmlod(complete_key inlkey[], wordint *nkey)
{
   
   void lit_bloc();
   int verbar(),calc_checksum(), trouve_best_free();

   int i, j, espace_requis=0, *p, best_fit;
   int slice_lng, *indices, indic[NCLEMAX];
   int t_espace_requis_max;
#if !defined (ALL64)
   int tempor;
#endif
   int qvmindex_from_key(), vmmerr();

   wordint lslice_lng, *lbest_fit;
   double *dlbest_fit;
       
   if(callc("VMMLOD"))   ;
       


   if(pwd_set)
        return(vmmerr("VMMLOD",PASSWORD_IS_SET));

   

   indices = &indic[0];

/*
 *   boucle sur toutes les clefs
 *   pour cette version, on ne peut depasser NCLEMAX clefs par appel
 */

   if(*nkey > NCLEMAX)
       return(vmmerr("VMMLOD",TOO_MANY_KEYS));

   for (i=0,p=indices; i < *nkey; i++,p++)
   {
      if ((*p = qvmindex_from_key(inlkey[i])) < 0)
         return vmmerr("VMMLOD",*p);

   }

/*
 *   trouver l'espace requis pour contenir toutes les tranches
 *   si une tranche est deja en memoire, on ne calcule pas son 
 *   espace mais on met le marqueur keep_in_core a 1 afin d'eviter
 *   son ejection par pack_blocks()
 */

   for (i=0,p=indices; i < *nkey; i++,p++) {
      if (! SLICES[*p].info.flags.is_in_core) {
#if !defined (ALL64)
         tempor = NAMES[SLICES[*p].name_table_index].lslice;
         espace_requis += (tempor *= (SLICES[*p].info.flags.size8 == 1) ? 2 : 1);
#else
         espace_requis += (NAMES[SLICES[*p].name_table_index].lslice);
#endif
	 }
         else
         {
            SLICES[*p].info.flags.keep_in_core = 1;
            BLOCKS[SLICES[*p].block_table_index].info.flags.keep_in_core = 1;
            if(debug_mode)
                 fprintf(fdout," VMM-trace champ deja en memoire:%s, tranche%d\n",
                   NAMES[SLICES[*p].name_table_index].nom,*p-NAMES[SLICES[*p].name_table_index].major_key+1);
         }
      }
      if(debug_mode)
      {
         fprintf(fdout," VMMLOD-Espace requis = %d\n",espace_requis);
         if(champs_bloques)
            fprintf(fdout," VMMLOD-Nombre de champs deja bloques: %d\n",
                           champs_bloques);
      }

/*
 *    calculer l'espace requis maximal pour
 *    fin diagnostique
 */
      t_espace_requis_max = 0;
      for (i = 0; i < nbblocks; i++)
      {
           if(BLOCKS[i].info.flags.keep_in_core)
                t_espace_requis_max+= BLOCKS[i].size;
      }

      t_espace_requis_max+= espace_requis;
      espace_requis_max = espace_requis_max > t_espace_requis_max ?
                          espace_requis_max : t_espace_requis_max ;
/*
 *   chercher un "best fit"
 */

   best_fit = trouve_best_free(espace_requis);
   if (best_fit == -1) 
   {
            fprintf(fdout," VMMLOD-Espace requis minimal = %d\n",
                            espace_requis_max);
	    return(vmmerr("VMMLOD",NO_SPACE_LEFT_FOR_LOAD));
   }
   
   for (i=0, p=indices; i < *nkey; i++,p++) {
      if (! SLICES[*p].info.flags.is_in_core) {
	 slice_lng = NAMES[SLICES[*p].name_table_index].lslice;

#if ! defined (ALL64)
         slice_lng *=
            (SLICES[*p].info.flags.size8 ==1) ? 2:1;
#endif
	 if (BLOCKS[best_fit].size != slice_lng) 
         { 
/*
 *         fragmentation du bloc libre
 */
           nbblocks++;
           nbblocks_max = nbblocks_max > nbblocks ?
                          nbblocks_max : nbblocks ;
	   for (j=nbblocks-1; j > best_fit+1; j--)
           {
	      BLOCKS[j].info.attributs = BLOCKS[j-1].info.attributs;
	      BLOCKS[j].slice_table_index = BLOCKS[j-1].slice_table_index;
       	      BLOCKS[j].memadr = BLOCKS[j-1].memadr;
       	      BLOCKS[j].file_adr = BLOCKS[j-1].file_adr;
	      BLOCKS[j].size = BLOCKS[j-1].size;
              BLOCKS[j].prev_fb = (BLOCKS[j-1].prev_fb >= best_fit?
                      BLOCKS[j-1].prev_fb+1:BLOCKS[j-1].prev_fb);
              BLOCKS[j].next_fb = (BLOCKS[j-1].next_fb != -1 ?
                      BLOCKS[j-1].next_fb+1 : -1);
	   }
	   for (j=0; j < nbslices; j++)
	      if (SLICES[j].block_table_index > best_fit)
	         SLICES[j].block_table_index++;


	   BLOCKS[best_fit+1].info.attributs = 0;
	   BLOCKS[best_fit+1].slice_table_index = -1;
	   BLOCKS[best_fit+1].file_adr = -1;
           BLOCKS[best_fit+1].size = BLOCKS[best_fit].size - slice_lng;
	   BLOCKS[best_fit+1].memadr = BLOCKS[best_fit].memadr + slice_lng;
           BLOCKS[best_fit+1].prev_fb = BLOCKS[best_fit].prev_fb;
           BLOCKS[best_fit+1].next_fb = (BLOCKS[best_fit].next_fb != -1?
                         BLOCKS[best_fit].next_fb+1: -1);

           if(BLOCKS[best_fit+1].prev_fb != -1)
             BLOCKS[BLOCKS[best_fit+1].prev_fb].next_fb=best_fit+1;

           if(BLOCKS[best_fit+1].next_fb != -1)
             BLOCKS[BLOCKS[best_fit+1].next_fb].prev_fb=best_fit+1;

           if(first_free_bloc == best_fit)
                 first_free_bloc++;
         }
         else
         {
             /* rechainage des blocs libres */
             if (BLOCKS[best_fit].prev_fb != -1)
             {
                 BLOCKS[BLOCKS[best_fit].prev_fb].next_fb=BLOCKS[best_fit].next_fb;
                 if(first_free_bloc == best_fit)
                        first_free_bloc = BLOCKS[best_fit].prev_fb;
             }
             if(BLOCKS[best_fit].next_fb != -1)
             {
                 BLOCKS[BLOCKS[best_fit].next_fb].prev_fb = BLOCKS[best_fit].prev_fb; 
                 if(first_free_bloc == best_fit)
                         first_free_bloc = BLOCKS[best_fit].next_fb;
             }
         }
	 SLICES[*p].info.flags.is_in_core = 1;
	 SLICES[*p].info.flags.keep_in_core = 1;
	 SLICES[*p].info.flags.in_used = 1;
	 SLICES[*p].block_table_index = best_fit;
         BLOCKS[best_fit].next_fb = BLOCKS[best_fit].prev_fb =  -1;
	 BLOCKS[best_fit].slice_table_index = *p;
	 BLOCKS[best_fit].info.attributs = SLICES[*p].info.attributs;
	 BLOCKS[best_fit].size = slice_lng;
	 BLOCKS[best_fit].file_adr = 
	    (NAMES[SLICES[*p].name_table_index].base_file_adr == -1) ? -1 :
	   NAMES[SLICES[*p].name_table_index].base_file_adr + (slice_lng *
           (inlkey[i].key.key_minor == 0 ? 0: inlkey[i].key.key_minor-1));
         if ((BLOCKS[best_fit].info.flags.traced) || debug_mode)
             fprintf(fdout,"VMM trace: chargement memoire  de variable %s tranche %d en position %d\n",
              NAMES[SLICES[*p].name_table_index].nom,
              *p - NAMES[SLICES[*p].name_table_index].major_key + 1,best_fit);
	 if (SLICES[*p].info.flags.disk_image) {
/*
 *    lire la tranche sur disque
 */
	    lit_bloc(best_fit,BLOCKS[best_fit].info.flags.class,
		       BLOCKS[best_fit].memadr,BLOCKS[best_fit].file_adr,
		       BLOCKS[best_fit].size);
         }
	 else if (SLICES[*p].info.flags.init) {
/*
 *    initialiser le champ
 */
            lslice_lng = (wordint) slice_lng;
            if (SLICES[*p].info.flags.size8)
            {
#if ! defined (ALL64)
              lslice_lng /= 2 ;
#endif
              dlbest_fit = (double *) BLOCKS[best_fit].memadr;
	      switch (SLICES[*p].info.flags.init)
              {
	      
	         case 1:
	            f77name(afix8)(dlbest_fit,&zero8,&lslice_lng);
	            break;
                 case 2:
                    f77name(afix8)(dlbest_fit,&MAXVAL8,&lslice_lng);
                    break;
                 default:
                    fprintf(fd_err,"vmmlod error: bad init mode, init =%d",
			  SLICES[*p].info.flags.init);
                    break;
              }
            }
            else
            {
              lbest_fit = (wordint *) BLOCKS[best_fit].memadr;
              switch (SLICES[*p].info.flags.init)
              {
	      
                 case 1:
                    f77name(afix)(lbest_fit,&zero,&lslice_lng);
                    break;
                 case 2:
                    f77name(afix)(lbest_fit,&MAXVAL,&lslice_lng);
                    break;
                 default:
                    fprintf(fd_err,"vmmlod error: bad init mode, init =%d",
			  SLICES[*p].info.flags.init);
                    break;
               }
             }
/*
 *           Ecrire le marqueur de bloc a la fin
 */
             INTBAR(best_fit);
/*
 *           Calculer le checksum si requis
 */
             if(SLICES[*p].info.flags.do_checksum || checksum_mode)
                  SLICES[*p].checksum = calc_checksum(best_fit);

	 }
         else
         {
/*
 *           Ecrire le marqueur de bloc a la fin
 */
             INTBAR(best_fit);
/*
 *           Calculer le checksum si requis
 */
             if(SLICES[*p].info.flags.do_checksum || checksum_mode)
                  SLICES[*p].checksum = calc_checksum(best_fit);
         }

         verbar(best_fit);

	 best_fit++;
	 }
      } /* end for */


/*
 *       fprintf(fdout," espace_requis_max %d\n",espace_requis_max);
 */
   return 0;
}



/***s/p reserve_disk_space 
*    objet(reserve_disk_space) - reserver la memoire disque pour
*                                l'ecriture des variables.  L'espace disque
*                                n'est reserve qu'au moment opportun.
*                                On reserve l'espace pour toutes les tranches
*                                d'une variable.
*
*auteur J. Caveen - mai 1994
*
*arguments
*         in  bkno - numero du block memoire
*
*/
PRIVATE void
reserve_disk_space(int bkno)
{

   int ind,cl,i,slice_lng;
   wordint lun, lpos, liun, *bidon ,nmots;

   ind = SLICES[BLOCKS[bkno].slice_table_index].name_table_index;

   if(debug_mode)
   {
      fprintf(fdout," RESERVE_DISK_SPACE-Allocation d'espace disque, variable=%s,lslice=%d,nslice=%d\n",
      NAMES[ind].nom,NAMES[ind].lslice,NAMES[ind].nslice);
   }

   cl  = NAMES[ind].class;
   nmots =  (wordint) NAMES[ind].lslice;
#if !defined (ALL64)
   nmots *= (BLOCKS[bkno].info.flags.size8 == 1) ? 2 :1;
#endif

/*
 * init. des adresses d'ecriture
 */
   NAMES[ind].base_file_adr = wp_Vmm[cl - 1];
   lpos = (wordint) wp_Vmm[cl - 1];
   bidon = BLOCKS[0].memadr;
   liun = (wordint) fclass[cl - 1];
/*
 *  on reserve l'espace pour toutes les tranches de la variable
 */
   for( i = 0; i < NAMES[ind].nslice; i++)
   {
       f77name(wawrit)(&liun,bidon, &lpos,&nmots);
       lpos+=nmots;
   }
/*
 * Mise a jour de la position d'ecriture dans le premier mot du fichier
 */
   wp_Vmm[cl - 1] =  (int) lpos;
   lun = 1;
   f77name(wawrit)(&liun,&lpos,&lun,&lun);

/*
 * faire la mise a jour de tous les blocs memoires associes
 * a la meme variable: on calcule la bonne adresse fichier 
 * pour chacun de ces blocs
 */
   
   slice_lng = NAMES[ind].lslice;
#if ! defined (ALL64)
         slice_lng *= (SLICES[NAMES[ind].major_key].info.flags.size8 ==1) ? 2:1;
#endif

   for (i = 0; i < NAMES[ind].nslice ; i++)
   {
     if(SLICES[NAMES[ind].major_key + i].block_table_index != -1)
        BLOCKS[SLICES[NAMES[ind].major_key + i].block_table_index].file_adr = 
           NAMES[ind].base_file_adr + (slice_lng * i);
   }

}



/***s/p strfind
*
*objet(strfind)
*     Recherche d'une sous-chaine dans une chaine de caractere
*     et retourne la position trouvee (-1 si pas trouvee)
*
*auteur Y. Chartier
*
*revision 001 - M. Lepine - comparaison insensible aux maj./min. +
*                           documentation
*argument
*     in  SousChaine   chaine de caractere recherchee
*     in  Chaine       chaine complete
*
**/
PRIVATE int 
strfind(char *SousChaine, char *Chaine)
{
int i,j, LongueurChaine, LongueurSousChaine, PositionTrouvee;

PositionTrouvee = -1;
LongueurChaine = strlen (Chaine);
LongueurSousChaine = strlen (SousChaine);

if (LongueurChaine < LongueurSousChaine)
   return(PositionTrouvee);

i=0;
while (i < LongueurChaine)
      {
      if (tolower(SousChaine[0]) == tolower(Chaine[i]))
	 {
	 j = i + LongueurSousChaine - 1;

	 if (j > LongueurChaine)
	    return(PositionTrouvee);

	 while (j > i && tolower(SousChaine[j-i]) == tolower(Chaine[j]))
	       j--;

         if (j==i)
	    PositionTrouvee = i;
         }
      i++;
      }
return(PositionTrouvee); 
}



/***s/p swap_blocks
*
*objet(swap_blocks)
*     Intervertir un bloc vide avec un bloc non vide deplacable
*
*auteur M. Lepine  -  juillet 1993
*       revision: E. Gondet - 19 avril 2002
*		  Porting on FUJITSU (__uxpv__) + appel libmp pour copie de blocs: _MmCopy
*
*arguments
*     in   i      index du bloc vide
*     in   inext  index du bloc non vide deplacable
*
**/
PRIVATE void 
swap_blocks(int i,int inext)

{

    void imprime();
   int t_attributs, t_slice_table_index, t_size;
   int  t_file_adr;
   wordint *dest_adr, *src_adr;
   wordint nmots, ind;
   extern void f77name(movlev8)();

   ind = BLOCKS[inext].slice_table_index;
   if ((BLOCKS[inext].info.flags.traced) || debug_mode)
      fprintf(fdout,"VMM trace: deplacement du bloc %d variable %s tranche %d en position %d\n",
      inext,NAMES[SLICES[ind].name_table_index].nom,
      ind - NAMES[SLICES[ind].name_table_index].major_key + 1,i);
   dest_adr =  BLOCKS[i].memadr;
   src_adr =  BLOCKS[inext].memadr;
   nmots = (wordint) BLOCKS[inext].size;
#if ! defined (ALL64)
   nmots *= BLOCKS[inext].info.flags.size8 == 1 ? 2 : 1;
#endif
   t_attributs = BLOCKS[i].info.attributs;
   t_slice_table_index = BLOCKS[i].slice_table_index;
   t_size = BLOCKS[i].size;
   t_file_adr = BLOCKS[i].file_adr;
   BLOCKS[i].info.attributs = BLOCKS[inext].info.attributs;
   BLOCKS[i].slice_table_index = BLOCKS[inext].slice_table_index;
   BLOCKS[i].size = BLOCKS[inext].size;
   BLOCKS[i].file_adr = BLOCKS[inext].file_adr;
   BLOCKS[inext].info.attributs = t_attributs;
   BLOCKS[inext].slice_table_index = t_slice_table_index;
   BLOCKS[inext].memadr = BLOCKS[i].memadr + (wordint) BLOCKS[i].size;
   BLOCKS[inext].size = t_size;
   BLOCKS[inext].file_adr = t_file_adr;
 /* rechainage des blocs libres */
   if(BLOCKS[i].next_fb != -1)
        BLOCKS[BLOCKS[i].next_fb].prev_fb = inext;
   if(BLOCKS[i].prev_fb != -1)
        BLOCKS[BLOCKS[i].prev_fb].next_fb = inext;
   if(first_free_bloc == i)
        first_free_bloc = inext;

   BLOCKS[inext].next_fb = BLOCKS[i].next_fb;
   BLOCKS[inext].prev_fb = BLOCKS[i].prev_fb;
   BLOCKS[i].prev_fb = BLOCKS[i].next_fb = -1;
   SLICES[BLOCKS[i].slice_table_index].block_table_index = i;

#if defined(__uxpv__)
   _MmCopy(dest_adr,src_adr,nmots*sizeof(wordint));
#else 
#  if defined(NEC)
   f77name(movlev8)(src_adr,dest_adr,&nmots);      /* assure la vectorisation */
#  else
   memcpy(dest_adr,src_adr,nmots*sizeof(wordint));
#  endif
#endif


/*
   imprime();
*/
   }


/***s/p trie_le_tableau - trier un tableau par ordre croissant
 *
 * Auteur : J. Caveen - mai 1994
 *
 * Objet(trie_le_tableau) _ faire un bubble sort sur un tableau d'entiers
 *                          le tableau est trie par ordre croissant
 *
 * Arguments
 *           in-out  - table    - tableau a trier
 *           in      - longueur - nombre d'elements dans le tableau
 **/
PRIVATE void
trie_le_tableau(int *table,int longueur)
{
   int i, j;
   
   
   for(i = 0; i< longueur-1; i++)
   {
      for(j=longueur-1; j > i; j--)
      {
	 if(*(table+j) < *(table+j-1))
	 {
	    *(table+j) = *(table+j)^*(table+j-1);
	    *(table+j-1) = *(table+j)^*(table+j-1);
	    *(table+j) = *(table+j)^*(table+j-1);
	 }
      }
   }
}


/***s/p trouve_best_fit
 *
 *  Auteur : James Caveen -  mai 1994
 *
 *  Objet (trouve_best_fit) - trouver le plus petit bloc memoire 
 *                            assez gros pour contenir un champ de
 *                            dimension size.
 *                            La recherche d'un bloc se fait a partir
 *                            du premier bloc libre en parcourant la chaine
 *                            des blocs libres jusqu'au bout.
 *
 * Argument:   in size - grosseur du champ a charger en memoire
 **/
int trouve_best_fit(int size)
{

     int i,best_fit, min_diff;

     best_fit = -1;
     if (first_free_bloc == -1)
          return(best_fit);

     i = first_free_bloc;
     min_diff =9999999;
     while( i != -1)
     {
        if(BLOCKS[i].size >= size)
        {
             if(min_diff > (BLOCKS[i].size - size))
             {
                 min_diff = BLOCKS[i].size - size;
                 best_fit = i;
             }
        }
         if(min_diff == 0) 
               return(best_fit);

         i = BLOCKS[i].next_fb;
     }
          
     return(best_fit);
}


/***s/p trouve_best_free - trouver le meilleur bloc memoire pour 
 *                         charger un champ.
 * Auteur James Caveen - mai 1994
 *
 * Objet (trouve_best_free) - fonction servant a trouver le meileur bloc 
 *                            memoire pouvant loger un champ de dimension
 *                            size.  
 *                            Methode : 1 - chercher un best_fit parmis
 *                                          les blocs vides.
 *                                      2 - Chercher le segment pouvant fournir
 *                                          l'espace a moindre cout, ejecter
 *                                          des champs en fonctions des poids
 *                                          si necessaire et comprimer
 *                                          le segment.
 *                            La fonction retourne l'indice du bloc ayant
 *                            la meilleure grosseur. Si aucun bloc n'est
 *                            trouve, on retourne -1.
 *
 * Argument : in size - grosseur du champ a charger en memoire
 *
 **/
int trouve_best_free(int size)
{
     int trouve_best_fit(), pack_segment(), eject_from_tableau();
     int trouve_best_segment();
     int best_fit, bkno,grosseur;
     int tableau_eject_index;

     best_fit=trouve_best_fit(size);

     if(best_fit > -1) return (best_fit);
/*
 *   Il n'y a pas de bloc libre assez gros
 *   on cherche un segment appropie;
 *   bkno pointe au debut du segment
 */
     bkno = trouve_best_segment(size,&tableau_eject_index);

     if(bkno == -1)  /* Aucun segment approprie */
        return (bkno);

/*
 *   on ejecte les champs du segment jusqu'a
 *   l'obtention de l'espace requis
 *   Cette operation doit se faire avant
 *   l'appel a pack_segment puisque la liste
 *   de blocs a ejecter est contenu dans une
 *   liste ordonnee par poids dans le tableau
 *   tableau_eject
 *   L'operation pack_segment rend tableau_eject inutilisable
 */
     grosseur = eject_from_tableau(size,tableau_eject_index);
     
/*
 *   On comprime le segment
 */
     grosseur = pack_segment(bkno,&best_fit);

     return (best_fit);


}










/***s/p trouve_best_segment - trouver le meilleur segment pouvant 
 *                           contenir un bloc de grosseur size
 *                           Un segment est une suite de blocs compris
 *                           entre deux blocs bloques non-contigues.
 *
 * Auteur: J. Caveen - juin 1994
 *
 *Objet (trouve_best_segment)- fonction servant a determiner quel
 *                             segment memoire est le plus approprie
 *                             pour recevoir un bloc de grosseur size.
 *                             La fonction retourne le numero du premier
 *                             bloc du segment.  
 *
 * Methode:       Parcourir chaque bloc de chaque segment et 
 *                dresser une liste de tous les blocs.  Chaque entree de la 
 *                liste est un entier qui contient le numero du segment,
 *                le poids du bloc et le numero du bloc:
 *
 *                  17 bits     4 bits   11 bits
 *                | no segment | poids | no de bloc |
 *
 *                Note: un bloc libre a un poids de zero
 *                Cette liste triee se trouve dans le tableau tableau_eject.
 *                La fonction retourne l'indice du premier bloc
 *                du segment pouvant fournir l'espace necessaire 
 *                au moindre cout.
 *
 *                  
 *
 *
 * Argument in size - grosseur du champ a charger en memoire
 *          out tableau_eject_index - indice dans le tableau trie tableau_eject
 *                                    a partir d'ou on fait l'ejection.
 **/

PRIVATE int
trouve_best_segment(int size, int *tableau_eject_index)
{

   void trie_le_tableau();
   int nseg, i, index, depart, bkno;
   int nblocs;

   int poid_max , npoid_max ;
   int t_poid_max, t_npoid_max, t_nblocs, t_bkno;
   int t_size, t_tablo_index;
   int l_poid;
   index = nseg = 0;

   bkno = -1;

   nblocs = 0;
   poid_max =  npoid_max =  999999;
   for (i = 0; i < nbblocks; i++)
   {
        /* avancer au premier bloc non- bloque */
        while(BLOCKS[i].info.flags.locked && i < nbblocks)
               i++;

        if(i >= nbblocks)
          break;

        nseg++;
        depart = index;
        while(! BLOCKS[i].info.flags.locked && i < nbblocks)
        {
            if(! BLOCKS[i].info.flags.keep_in_core)
            {
	         if(! BLOCKS[i].info.flags.in_used )
	            l_poid = 0;
	         else
	            l_poid = (BLOCKS[i].info.flags.must_exist == 1 ?
			 BLOCKS[i].info.flags.weight :
			 BLOCKS[i].info.flags.weight + 9);

                 tableau_eject[index] = (nseg & SEGMENT_MASK)<< SEGMENT_SHIFT |
                                (l_poid & WEIGHT_MASK)<< WEIGHT_SHIFT |
                                 i & BKNO_MASK;
                 index++;
            }
            i++;
        }
   }

   /* on trie le segment par ordre croissant */
   trie_le_tableau(&tableau_eject[0],index);


/*
 *   trouver le segment le moins couteux a utiliser
 */

   i = 0;
   *tableau_eject_index = 0;
   while(i < index)
   {
     nseg = tableau_eject[i]>>SEGMENT_SHIFT & SEGMENT_MASK;
     t_tablo_index = i;
     t_nblocs = t_poid_max =  t_npoid_max = t_size = 0;
     t_bkno = tableau_eject[i] & BKNO_MASK;
     while((((tableau_eject[i]>>SEGMENT_SHIFT) & SEGMENT_MASK) == nseg) && (i < index) )
     {
        if(((tableau_eject[i]>>WEIGHT_SHIFT) & WEIGHT_MASK) != t_poid_max)
        {
             t_poid_max = (tableau_eject[i]>>WEIGHT_SHIFT) & WEIGHT_MASK;
             t_npoid_max = 0;
        }

        t_bkno = (tableau_eject[i] & BKNO_MASK) < t_bkno ? (tableau_eject[i] & BKNO_MASK): t_bkno;
        t_nblocs++;
        t_npoid_max++;
        t_size += BLOCKS[tableau_eject[i] & BKNO_MASK].size;
        i++;
        if(t_size >= size)
        {
            while((((tableau_eject[i]>>SEGMENT_SHIFT) & SEGMENT_MASK) == nseg) && ( i < index))
                i++; 

            break;
        }
     }
        

 /*
  *  Si on a assez d'espace, on evalue son cout d'utilisation
  *  relatif aux autres segments utilisables
  */

     if(t_size >= size)
     {
       if (((t_poid_max < poid_max) ||
             ((t_poid_max == poid_max) && (t_npoid_max < npoid_max))))
       {
             bkno =  t_bkno ;
             poid_max =  t_poid_max ;
             npoid_max = t_npoid_max ;
             nblocs  = t_nblocs;
             *tableau_eject_index = t_tablo_index;
       } 
     }
   }
         
    if(debug_mode)
    {
       if(bkno > -1)
       {
          printf(" On utilise:bkno %d poid_max %d npoid_max %d nblocs %d size %d\n",
               bkno,poid_max,npoid_max,nblocs,t_size);
       }
       else
       {
	  printf(" Aucun segment assez gros pour chargement\n");
       }
    } 

    return (bkno);

}

/***s/p verbar
*
*objet(verbar)
*          Fonction servant a verifier si les marqueurs de debut et de fin
*          d'un bloc memoire sont intacts
*
*auteur J. Caveen  -  avril 1994
*
*arguments
*       bkno - entree - numero du bloc pour lequel on verifie les marqueurs
*
**/

PRIVATE int
 verbar(int bkno)
{

      int vmmerr();
      void imprime_structures();
      void impval();
      extern void f77name(tracebak)();

      int erreur;

/*  si les marqueurs sont endommages, on imprime un message d'erreur, et on met 
 *  fin a l'execution du programme.
 *  Sinon, on retourne la valeur 0
 */ 

    erreur = 0;
    if(! BLOCKS[bkno].info.flags.in_used)
         return(erreur);

    if(*(BLOCKS[bkno].memadr -1)^BARVAL)
    {
        erreur++;
        fprintf(fdout," ERROR - BEGINNING BLOCK DELIMITOR FOR BLOCK %d IS DAMAGED\n",bkno);

        if((bkno > 0) && (BLOCKS[bkno-1].info.flags.in_used))
        {
           fprintf(fdout,"       - POSSIBLE MEMORY OVERLAP: VARIABLE %s, SLICE %d\n",
                             NAMES[SLICES[BLOCKS[bkno-1].slice_table_index].name_table_index].nom,
                             BLOCKS[bkno-1].slice_table_index -
                             NAMES[SLICES[BLOCKS[bkno-1].slice_table_index].name_table_index].major_key + 1);
           fprintf(fdout,"                              AND VARIABLE %s, SLICE %d\n",
                             NAMES[SLICES[BLOCKS[bkno].slice_table_index].name_table_index].nom,
                             BLOCKS[bkno].slice_table_index -
                             NAMES[SLICES[BLOCKS[bkno].slice_table_index].name_table_index].major_key + 1);

           fprintf(fdout,"BLOCK DELIMITOR ADDRESS +- 2 WORDS\n");
           impval(BLOCKS[bkno].memadr-2);
        }
        else
        {
           fprintf(fdout,"       - POSSIBLE ADDRESSING ERROR: VARIABLE %s, SLICE %d\n",
                             NAMES[SLICES[BLOCKS[bkno].slice_table_index].name_table_index].nom,
                             BLOCKS[bkno].slice_table_index -
                             NAMES[SLICES[BLOCKS[bkno].slice_table_index].name_table_index].major_key + 1);

           fprintf(fdout,"BLOCK DELIMITOR ADDRESS + 4  WORDS\n");
            impval(BLOCKS[bkno].memadr-1);
        }

    }

    if(*(BLOCKS[bkno].memadr+BLOCKS[bkno].size-1)^BARVAL)
    {
        erreur++;
        fprintf(fdout," ERROR - END BLOCK DELIMITOR FOR BLOCK %d IS DAMAGED\n",bkno);
        if((bkno <  (nbblocks-1)) && (BLOCKS[bkno+1].info.flags.in_used))
        {
           fprintf(fdout,"       - POSSIBLE MEMORY OVERLAP: VARIABLE %s, SLICE %d\n",
                             NAMES[SLICES[BLOCKS[bkno].slice_table_index].name_table_index].nom,
                             BLOCKS[bkno].slice_table_index -
                             NAMES[SLICES[BLOCKS[bkno].slice_table_index].name_table_index].major_key + 1);
           fprintf(fdout,"         AND VARIABLE %s, SLICE %d\n",
                             NAMES[SLICES[BLOCKS[bkno+1].slice_table_index].name_table_index].nom,
                             BLOCKS[bkno+1].slice_table_index -
                             NAMES[SLICES[BLOCKS[bkno+1].slice_table_index].name_table_index].major_key + 1);

           fprintf(fdout,"BLOCK DELIMITOR ADDRESS +- 2 WORDS\n");
            impval(BLOCKS[bkno].memadr+BLOCKS[bkno].size-2);
        }
        else 
        {
           fprintf(fdout,"       - POSSIBLE ADDRESSING ERROR: VARIABLE %s, SLICE %d\n",
                             NAMES[SLICES[BLOCKS[bkno].slice_table_index].name_table_index].nom,
                             BLOCKS[bkno].slice_table_index -
                             NAMES[SLICES[BLOCKS[bkno].slice_table_index].name_table_index].major_key + 1);

           fprintf(fdout,"BLOCK DELIMITOR ADDRESS - 4 WORDS\n");
           impval(BLOCKS[bkno].memadr+BLOCKS[bkno].size-5);
        }
    }

    if(erreur) 
    {
      if(debug_mode)
      {
         imprime_structures(2);
         imprime_structures(1);
         imprime_structures(0);
      }
      return (vmmerr("VERBAR",BLOCK_DAMAGE));
    }


    return(0);
}


/***s/p vmmallc2
*
*objet(vmmallc2)
*     Reserver l'espace memoire necessaire au systeme de gestion de
*     memoire virtuelle. 
*     vmmallc doit ete la premiere fonction du systeme de gestion de
*     memoire virtuelle appelee par l'usager.
*
*     vmmallc fait l'ouverture de tous les fichiers associes au systeme
*     de gestion de memoire virtuelle (VMM_01,... et VMM_controle)
*     et initialise les wp_Vmm (mot ou ecrire a la fin des fichiers)
*
*     vmmallc relit le fichier de controle si il existe et est non vide
*
*auteur M. Lepine  -  juin 1993
*       J. Caveen
*
*       revision: J. Caveen - 1 avril 1994
*                 ajout d'un marqueur de bloc 64 bits au debut et
*                 decalage de l'adresse de depart du bloc[0] en
*                 consequence.
*       revision: E. Gondet - 19 avril 2002
*                 Repertoire en argument pour decouplage en parallele des fichiers
*		  Porting on FUJITSU (__uxpv__)
*argument
*     in  memry    quantite de memoire a reserver (en mots)
*     in  cd_repertoire Nom du repertoire de travail pour le processus vmm
*
**/
wordint
f77name(vmmallc2)(wordint *memry, char *cd_rep, F2Cl lng)
{
   int vmmerr();
   void  lit_vmm_controle();
/*
      void  ouvre_ou_ferme_controle();
*/
 
   int i;
  
#if defined (NEC) && defined (_FLOAT0) && !defined(__uxpv__)
   void *f77name(malloc2)();
   long long noctets; 
#else
   size_t noctets;
#endif
/*
 * s'assurer que c'est bien le premier appel a vmmallc
 */

   if(called_vmmallc)
       return vmmerr("VMMALLC",VMMALLC_ALREADY_CALLED);

/* 
 * initialiser les structures
 */
   for (i = 0; i < MAXSLICES; i++)
          SLICES[i].info.attributs = 0;
   for (i = 0; i < MAXBLOCKS; i++)
   {
          BLOCKS[i].info.attributs = 0;
          BLOCKS[i].prev_fb = BLOCKS[i].next_fb = -1;
   }
/*
 * ouverture de tous les fichiers du systeme de gestion 
 */
   strncpy(cd_repertoire,cd_rep,((lng > 256) ? 256 : lng));
   cd_repertoire[lng] = '\0';
   ouvre_ou_ferme_controle(1,1,"VMMALLC");

/*
 * relire toutes les variables et toutes les slices du fichier de controle
 * en cas de reprise
 */
   lit_vmm_controle();
       

   maxmem = (int) *memry;
   /*   fprintf(stdout,"Debug vmmallc maxmem=%d \n",maxmem); */
/* Rendre maxmem pair */
   maxmem = maxmem + (maxmem & 1);
/*
 * reserver la memoire demandee
 */

   noctets = maxmem;
/*ETG noctets=noctets*4 si wordint sur 4 par 8 si wordint sur 8 */
   noctets <<= (sizeof(wordint)==4?2:3) ;
   free_space = maxmem;
#if defined (NEC) && defined (_FLOAT0) && !defined(__uxpv__)
   /*   fprintf(stdout,"Debug VMMALLC appel a malloc2\n"); */
   noctets += 8;
   BLOCKS[0].memadr = (wordint *) f77name(malloc2)(&noctets); 
   /*   fprintf(stdout,"Debug vmmallc memry=%d, maxmem=%d, noctets=%ld\n",*memry,maxmem,noctets); */
   /*   fprintf(stdout,"Debug vmmallc BLOCKS[0].memadr=%d\n",BLOCKS[0].memadr); */
#else
   BLOCKS[0].memadr = (wordint *) malloc(noctets + 8);
#endif
   if(BLOCKS[0].memadr == (wordint *) NULL)
        return(vmmerr("VMMALLC",NOT_ENOUGH_MEMORY));

/*
 *  Ecriture du marqueur de bloc dans les premiers 64 bits et
 *  decalage de l'adresse de depart du bloc[0]
 */
#if defined (ALL64)
   *(BLOCKS[0].memadr) = BARVAL; 
#else
   *(BLOCKS[0].memadr) = *(BLOCKS[0].memadr+1) = BARVAL;
#endif
   BLOCKS[0].memadr+=(8/sizeof(wordint));
   BLOCKS[0].slice_table_index = -1;
   BLOCKS[0].size = maxmem ;
   nbblocks++;
   first_free_bloc = 0;
   called_vmmallc = 1;

#if defined CALL_SEQ
     fprintf(fdout,"CALL- vmmallc(%d)\n",*memry);
#endif
   if(debug_mode)
         fprintf(fdout," VMMALLC-allocation memoire de %d mots\n",maxmem);

   return(0);
   }

/***s/p vmmallc
*
*objet(vmmallc)
*     Reserver l'espace memoire necessaire au systeme de gestion de
*     memoire virtuelle. 
*     vmmallc doit ete la premiere fonction du systeme de gestion de
*     memoire virtuelle appelee par l'usager.
*
*     vmmallc fait l'ouverture de tous les fichiers associes au systeme
*     de gestion de memoire virtuelle (VMM_01,... et VMM_controle)
*     et initialise les wp_Vmm (mot ou ecrire a la fin des fichiers)
*
*     vmmallc relit le fichier de controle si il existe et est non vide
*
*auteur M. Lepine  -  juin 1993
*       J. Caveen
*
*       revision: J. Caveen - 1 avril 1994
*                 ajout d'un marqueur de bloc 64 bits au debut et
*                 decalage de l'adresse de depart du bloc[0] en
*                 consequence.
*       revision: E. Gondet - 19 avril 2002
*                 Repertoire en argument pour decouplage en parallele des fichiers
*		  Porting on FUJITSU (__uxpv__)
*argument
*     in  memry    quantite de memoire a reserver (en mots)
*     in  cd_repertoire Nom du repertoire de travail pour le processus vmm
*
**/
wordint
f77name(vmmallc)(wordint *memry)
{
  wordint ier;
  char current_dir[2] = "./";
  
  ier = f77name(vmmallc2)(memry,current_dir,2);
  return(ier);
}


/***s/p vmmatt
*
*objet(vmmatt)
*     Retourner les caracteristiques d'une variable
*
*auteur J. Caveen  -  juin 1993
*
*argument
*     in   namevar  nom de la variable (max 8 car)
*     out  lpiece   dimension de chacun des morceaux qui composent le champ
*     out  npiece   nombre de morceaux qui composent le champ
*     out  attr     chaine de caracteres decrivant les attributs du champ
*
**/
wordint
f77name(vmmatt)(char *namevar,wordint *lpiece,wordint *npiece,char *attr,F2Cl l1,F2Cl l2)
{
   int vmmerr();

   int ind, cl,  wt, nsauv, ninit, nsize8,mexist, slice_num;
   char cnsauv, cninit ;
   int cnsize8;
   complete_key keyout;

   int i;
   char innamevar[9];

   if(callc("VMMATT"))  ;
       


   if(pwd_set)
        return(vmmerr("VMMATT",PASSWORD_IS_SET));

#if defined CALL_SEQ
     fprintf(fdout,
        "CALL- vmmatt(%s,%d,%d,%s)\n",namevar,*lpiece,*npiece,attr);
#endif
  strncpy(innamevar,namevar,l1);

/*
 * on remplie le reste de innamevar avec des blancs
 */
   for (i = l1; i < 8; i++)
          innamevar[i] = ' ';

   innamevar[8] = '\0';


/*
 * trouver la variable
 */
   ind = 0;
   while ((strncmp(innamevar,NAMES[ind].nom,8) != 0) && (ind < MAXNAMES))
      ind++;
   
   if(ind == MAXNAMES)
           return vmmerr("VMMATT", UNKNOWNVAR);

   *lpiece = (wordint) NAMES[ind].lslice - 8/sizeof(wordint); 
   *npiece = (wordint) NAMES[ind].nslice; 
   slice_num = NAMES[ind].major_key;
/*
 * recuperation des attributs
 */
   cl = SLICES[slice_num].info.flags.class;
   wt = SLICES[slice_num].info.flags.weight;
   nsauv = SLICES[slice_num].info.flags.save;
   ninit = SLICES[slice_num].info.flags.init;
   nsize8 = SLICES[slice_num].info.flags.size8;
   mexist = SLICES[slice_num].info.flags.must_exist;
 
   switch (nsauv)
   {
        case 1:
             cnsauv = 'Y';
                break;
        default:
             cnsauv = 'N';
   }
 
   switch (ninit)
   {
        case 1: 
            cninit  = '0';
            break;
        case 2:
            cninit  = 'R';
            break;
        default:
            cninit  = '-';
   }

   switch (nsize8)
   {
        case 1:
            cnsize8 = 8;
            break;
        default:
            cnsize8 = 0;
   }
        
/*
 * recomposer la liste d'attributs
 */
   sprintf(attr,"SAVE=%c,CL=%d,W=%d,INIT=%c,SIZE=%d%c",cnsauv,cl,wt,cninit,cnsize8,'\0');
   if(mexist)
       strcat(attr,",MUSTEXIST");

   keyout.clef = 0;
   keyout.key.key_major = slice_num;
   keyout.key.key_minor = 0;
   return (keyout.clef);
}

/***s/p vmmcks
*
*objet(vmmcks)
*          Calculer et retourner le chek sum du bloc
*          memoire associe a la clef inkey
*
*auteur J. Caveen  -  novembre 1993
*
*arguments
*         in  inkey  -  clef pour laquelle on fait le check-sum
*         in  mode   -  mode de calcul 
*                       pour cette version, mode=1 seulement => xor du champ
*
**/
wordint
f77name(vmmcks)(complete_key *inkey, wordint *mode)
{

     extern wordint f77name(qvmcks)();
     int qvmindex_from_key(), vmmerr();

     int slice_ind, i ;
     wordint check_sum, *adresse_du_bloc;
     wordint nbelem;

     if(callc("VMMCKS")) ;
       


#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmcks(%d,%d)\n",(wordint) inkey->clef, *mode);
#endif
     if(*mode != 1)
        return(vmmerr("VMMCKS",BAD_CKSUM_MODE)); 
            
     slice_ind = qvmindex_from_key(*inkey);

     if(slice_ind <0)
        return(vmmerr("VMMCKS",slice_ind));


     if(! SLICES[slice_ind].info.flags.is_in_core)
        return(vmmerr("VMMCKS",NOT_IN_CORE));


     adresse_du_bloc = BLOCKS[SLICES[slice_ind].block_table_index].memadr; 
     nbelem = (wordint) BLOCKS[SLICES[slice_ind].block_table_index].size;


/*
 *   on fait le check-sum mode 1
 */
     check_sum = f77name(qvmcks)(adresse_du_bloc,&nbelem,mode);
     return (check_sum);
}


/***s/p vmmcpk
 *
 *objet(vmmcpk)
 *    Faire un check-point ( sauver toutes les variables dans les fichiers Vmm_0n
 *    et sauver les tables de controle dans Vmm_controle
 *    On ferme les fichiers pour vider tous les tampons
 *
 *auteur J. Caveen  -  juillet 1993
 *
 *revision  james caveen - mars 1995
 *          forcer l'ejection des blocs qui n'ont pas l'attribut save
 *          
 **/
wordint f77name(vmmcpk)()
{
   
   void ecrit_vmm_controle(), ecrit_bloc();
   void reserve_disk_space();
   int vmmerr(), eject_block();
/*
      void  ouvre_ou_ferme_controle();
*/
   
   int i;
   
   if(callc("VMMCPK"))  ;
   
   
   
   if(pwd_set)
   return(vmmerr("VMMCPK",PASSWORD_IS_SET));
   
#if defined CALL_SEQ
   fprintf(fdout,"CALL- vmmcpk()\n");
#endif
   for (i = 0; i< nbblocks; i++)
   {
      if(BLOCKS[i].info.flags.in_used)
      {
	 if(BLOCKS[i].info.flags.keep_in_core)
	 return(vmmerr("VMMCPK",KEEP_IN_CORE_CPK));

	 if(BLOCKS[i].info.flags.save && (BLOCKS[i].info.flags.altered ||
					  BLOCKS[i].info.flags.was_altered))
	 {
	    if(BLOCKS[i].file_adr == -1)
	    reserve_disk_space(i);
	    ecrit_bloc(i,BLOCKS[i].info.flags.class,BLOCKS[i].memadr,
		       BLOCKS[i].file_adr,BLOCKS[i].size);
	 }
	 else
	    eject_block(i,0,0);
      }
   }
   
   ecrit_vmm_controle();
   ouvre_ou_ferme_controle(0,0,"VMMCPK");
   return 0;
}

/***s/p vmmcre
*
*objet(vmmcre)
*     Definir une variable et ses caracteristiques
*
*auteur M. Lepine  -  juin 1993
*       J. Caveen
*
*modification J. Cavven - mai 1994
*      - enlever la portion de code servant a reserver l'espace disque:
*        on reserve l'espace au moment opportun via reserve_disk_space()
*
*      - ajout de verification suplementaires quant a la longueur des
*        tranches et a l'integrite des caracteristiques des variables.
*
*argument
*     in   namevar  nom de la variable (max 8 car)
*     in   lpiece   dimension de chacun des morceaux qui composent le champ
*     in   npiece   nombre de morceaux qui composent le champ
*     in   inattr   chaine de caracteres decrivant les attributs du champ
*
*     in   l1       longueur de namevar
*     in   l2       longueur de inattr
**/
wordint
f77name(vmmcre)(char innamevar[],wordint *lpiece,wordint *npiece,
                char *inattr,F2Cl l1,F2Cl l2)
{
   int vmmerr();

   int i, pos, ind, cl, j, wt, nsauv, ninit, nsize8, mexist,slice_num;
   int innpiece,inlpiece;
   
   char junk[20], csauv, cinit, csize8, attr[NCARATTR], namevar[9];
   complete_key keyout;

   wordint nmots;

   if(callc("VMMCRE")) ;
       

   
   if(pwd_set)
        return(vmmerr("VMMCRE",PASSWORD_IS_SET));

/*
 * on remplie le reste de namevar avec des blancs
 */
   strncpy(namevar,innamevar,l1);
   for (i = l1; i < 8; i++)
          namevar[i] = ' ';

   namevar[8] = '\0';

   strncpy(attr,inattr,l2);
   attr[l2] = '\0';
#if defined CALL_SEQ
     fprintf(fdout,
        "CALL- vmmcre(%s,%d,%d,%s)\n",namevar,*lpiece,*npiece,attr);
     fprintf(fdout,"longueur de innamevar et attr: %d, %d\n",l1,l2);
#endif

   ind = 0;
   while ((strncmp(namevar,NAMES[ind].nom,8) != 0) && (ind < MAXNAMES))
      ind++;
/*
 * On cree une nouvelle variable
 */
   if (ind == MAXNAMES)
   {
/*
 *    si il s'agit d'une reprise, on s'assure que la variable 
 *    n'a pas l'attribut MUSTEXIST
 *    if(reprise && ((pos = strfind("MUSTEXIST",attr))  != -1))
 *          return vmmerr("VMMCRE",VAR_MUST_EXIST);
 *
 */
/*
 *    initialisation des attributs
 */
      ninit =  nsize8 = wt = mexist =  nsauv = 0;
      cl = 1;
      ind = nbvar++;
      strcpy(NAMES[ind].nom,namevar);
      NAMES[ind].base_file_adr = -1;
      NAMES[ind].nslice = (int) *npiece;
      NAMES[ind].major_key = nbslices;
/*      NAMES[ind].lslice = (int) *lpiece; */
       NAMES[ind].lslice = (int) ((((sizeof(wordint)*(*lpiece) +7)/8)*8)/sizeof(wordint));
      NAMES[ind].lslice+=(8/sizeof(wordint));
      if((pos = strfind("CL=",attr)) != -1)
         sscanf(&attr[pos],"%3s%d",junk,&cl);

      if((pos = strfind("W=",attr))  != -1)
         sscanf(&attr[pos],"%2s%d",junk,&wt);
      
      if((pos = strfind("SAVE=",attr))  != -1)
      {
         sscanf(&attr[pos],"%5s%1c",junk,&csauv);
         if ( (char) tolower(csauv) == 'y' )
              nsauv = 1;
      }

      if((pos = strfind("INIT=",attr))  != -1)
      {
         sscanf(&attr[pos],"%5s%1c",junk,&cinit);
         if ( (char) tolower(cinit) == '0' )
              ninit = 1;
         else if ( (char) tolower(cinit) == 'r' )
              ninit = 2;
      }

      if((pos = strfind("SIZE=",attr))  != -1)
      {
         sscanf(&attr[pos],"%5s%1c",junk,&csize8);
         if ( csize8 == '8')  
              nsize8 = 1;
      }
      

      if((pos = strfind("MUSTEXIST",attr))  != -1)
            mexist = 1;

      NAMES[ind].class = cl;
      for (j=0; j < *npiece; j++)
      {
	    SLICES[nbslices+j].info.flags.save = nsauv;
	    SLICES[nbslices+j].info.flags.size8 = nsize8;
	    SLICES[nbslices+j].info.flags.must_exist = mexist;
	    SLICES[nbslices+j].info.flags.class = cl;
	    SLICES[nbslices+j].info.flags.weight = wt;
	    SLICES[nbslices+j].info.flags.init = ninit;
	    SLICES[nbslices+j].name_table_index = ind;
	    SLICES[nbslices+j].block_table_index = -1;
      }
        

/*
 *    composer la clef de sortie
 */
      keyout.clef = 0;
      keyout.key.key_major = nbslices;
      keyout.key.key_minor = 0;
      nbslices += *npiece;
   }

/*
 * on modifie les attributs d'une variable
 * on ne peut modifier que save= et w=
 */
   else
   {
      slice_num = NAMES[ind].major_key;
      innpiece    = NAMES[ind].nslice;
      nsauv = SLICES[slice_num].info.flags.save;
      wt = SLICES[slice_num].info.flags.weight;
      ninit = SLICES[slice_num].info.flags.init;
      nsize8 = 0;
      
      if((pos = strfind("W=",attr))  != -1)
         sscanf(&attr[pos],"%2s%d",junk,&wt);
      
      if((pos = strfind("SAVE=",attr))  != -1)
      {
         sscanf(&attr[pos],"%5s%1c",junk,&csauv);
         if ( (char) tolower(csauv) == 'y' )
              nsauv = 1;
         else
              nsauv = 0;
      }
      if((pos = strfind("INIT=",attr))  != -1)
      {
         sscanf(&attr[pos],"%5s%1c",junk,&cinit);
         if ( (char) tolower(cinit) == '0' )
              ninit = 1;
         else if ( (char) tolower(cinit) == 'r' )
              ninit = 2;
         else
              ninit = 0;
      }
      if((pos = strfind("CL=",attr)) != -1)
         sscanf(&attr[pos],"%3s%d",junk,&cl);

      if((pos = strfind("SIZE=",attr))  != -1)
      {
         sscanf(&attr[pos],"%5s%1c",junk,&csize8);
         if ( csize8 == '8')
              nsize8 = 1;
      }

      if((pos = strfind("MUSTEXIST",attr))  != -1)
            mexist = 1;
      else
            mexist = 0;

/*
 *    si il s'agit d'une reprise, on s'assure que le fichier
 *    d'une variable MUSTEXIST est bien la si la variable a
 *    une adresse disque.

***Mis en commentaire: la verification se fait dans vmmallc
      if(reprise)
      {
       if (mexist &&
            ((NAMES[ind].base_file_adr != -1) && (wp_Vmm[cl - 1] <=2)))
       {
            fprintf(fdout," Variable %s must exist for a restart\n and file Vmm_0%d is absent\n",namevar,cl);
            return vmmerr("VMMCRE",VAR_MUST_EXIST);
       }
      }
*/


/* 
 *    on s'assure que les caracteristiques inchangeables
 *    sont bel et bien inchangees avant de faire les modifications
 */
      inlpiece = (int) ((((sizeof(wordint)*(*lpiece) +7)/8)*8)/sizeof(wordint));
      inlpiece +=(8/sizeof(wordint));


      if((*npiece != NAMES[ind].nslice) ||
         (inlpiece != NAMES[ind].lslice) ||
         (SLICES[slice_num].info.flags.class != cl) ||
         (SLICES[slice_num].info.flags.size8 != nsize8))
      {
             fprintf(fd_err,
                     " VMMCRE - THE CARACTERISTICS OF VARIABLE %s HAVE BEEN ALTERED\n",
                     NAMES[ind].nom); 
             fprintf(fd_err,
                 "npiece=%d, NAMES.nlsice=%d\nlpiece=%d, NAMES.lsice=%d\ncl=%d, NAMES.class=%d\nsize8=%d, NAMES.size8=%d\n",
      *npiece,NAMES[ind].nslice,inlpiece,NAMES[ind].lslice,
      cl,SLICES[slice_num].info.flags.class,
      nsize8,SLICES[slice_num].info.flags.size8);
                   
                 return (vmmerr("VMMCRE",ATTRIBUTS_MODIFIES));
      }

      for (j=0; j < innpiece; j++)
      {
	    SLICES[slice_num+j].info.flags.save = nsauv;
	    SLICES[slice_num+j].info.flags.weight = wt;
	    SLICES[slice_num+j].info.flags.init = ninit;
	    SLICES[slice_num+j].info.flags.must_exist = mexist;
      }
/*
 *    composer la clef de sortie
 */
      keyout.clef = 0;
      keyout.key.key_major = slice_num;
      keyout.key.key_minor = 0;
    
   }

/*
 *  s'assurer que la longueur d'une tranche n'est pas plus longue que
 *  l'espace demande dans vmmallc
 */
    nmots =  (wordint) NAMES[ind].lslice;
#if !defined (ALL64)
    nmots *= (nsize8 == 1) ? 2 :1;
#endif
    if (nmots > maxmem)
    {
         fprintf(fd_err," SLICE(S) OF VARIABLE %s CANNOT FIT IN THE ALLOCATED MEMORY\n",
                 NAMES[ind].nom);
         return(vmmerr("VMMCRE",SLICE_TOO_BIG));
    }

   return (keyout.clef);
}


/***s/p vmmdbg
*
*objet(vmmdbg)
*      Controle l'aide a la mise au point selon le contenu de "comand"
*      specifiant les options a activer.
*
*auteur M. Lepine  -  juillet 1993
*
*modification J.Caveen - mai 1994
*
*         - ajout de l'option CHECKSUM
*
*
*arguments
*         in command -  option a activer
*         in inlkey  -  liste de clefs des champs a charger
*         in nkey    -  nombre de clefs dans inlkey
*
**/
wordint
f77name(vmmdbg)(char command[],complete_key inlkey[], wordint *nkey,F2Cl l1)
{

  char cmd[NCARATTR], junk[20], diag_file[80], msg[80];
  int pos, ind, i, nvar, iii;


#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmdbg(%s,[",command);
        for (iii = 0; iii < *nkey; iii++)
            fprintf(fdout,"%d ",inlkey[iii]);

     fprintf(fdout,"],%d)\n",*nkey);
#endif
  strncpy(cmd,command,l1);
  cmd[l1] = '\0';

  if ((pos = strfind("OUTFILE=",cmd)) != -1) {
     sscanf(&cmd[pos],"%8s%s",junk,diag_file);
     fdout = fopen(diag_file,"w");
     }

  if ((pos = strfind("MSG=",cmd)) != -1) {
     strncpy(msg,&cmd[pos+4],l1-4);
     msg[l1-4] = '\0';
     fprintf(fdout,"%s\n",msg);
     }

  nvar = (int) ((inlkey[0].clef == -1) ? nbslices : *nkey);

  if ((pos = strfind("TRACE",cmd)) != -1) {
     for (i=0; i < nvar; i++) {
        ind = (inlkey[0].clef == -1) ? i : qvmindex_from_key(inlkey[i]);
	SLICES[ind].info.flags.traced = 1;
	}
     }

  if ((pos = strfind("CHECKSUM",cmd)) != -1) {
     for (i=0; i < nvar; i++) {
        ind = (inlkey[0].clef == -1) ? i : qvmindex_from_key(inlkey[i]);
	SLICES[ind].info.flags.do_checksum = 1;
	}
     }

  if ((pos = strfind("MEMDMP",cmd)) != -1) {
     for (i=0; i < nvar; i++) {
        ind = (inlkey[0].clef == -1) ? i : qvmindex_from_key(inlkey[i]);

                fprintf(fdout,"  Variable %s , tranche %d slice_table_index %d block_table_index %d\n",
                              NAMES[SLICES[ind].name_table_index].nom,
                              ind - NAMES[SLICES[ind].name_table_index].major_key + 1,ind,SLICES[ind].block_table_index);
                fprintf(fdout,"     keep_in_core       : %d\n",SLICES[ind].info.flags.keep_in_core);
                fprintf(fdout,"     is_in_core         : %d\n",SLICES[ind].info.flags.is_in_core);
                fprintf(fdout,"     in_used            : %d\n",SLICES[ind].info.flags.in_used);
                fprintf(fdout,"     locked             : %d\n",SLICES[ind].info.flags.locked);
                fprintf(fdout,"     save               : %d\n",SLICES[ind].info.flags.save);
                fprintf(fdout,"     altered            : %d\n",SLICES[ind].info.flags.altered);
                fprintf(fdout,"     was_altered        : %d\n",SLICES[ind].info.flags.was_altered);
                fprintf(fdout,"     traced             : %d\n",SLICES[ind].info.flags.traced);
                fprintf(fdout,"     hpa_alloc          : %d\n",SLICES[ind].info.flags.hpa_alloc);
                fprintf(fdout,"     disk_image         : %d\n",SLICES[ind].info.flags.disk_image);
                fprintf(fdout,"     size8              : %d\n",SLICES[ind].info.flags.size8);
                fprintf(fdout,"     must_exist         : %d\n",SLICES[ind].info.flags.must_exist);
                fprintf(fdout,"     class              : %d\n",SLICES[ind].info.flags.class);
                fprintf(fdout,"     weight             : %d\n",SLICES[ind].info.flags.weight);
                fprintf(fdout,"     do_checksum        : %d\n",SLICES[ind].info.flags.do_checksum);
                fprintf(fdout,"     init               : %d\n",SLICES[ind].info.flags.init);
	}
     }

  return(0);
  }


/***s/p vmmdiag
*
*objet(vmmdiag)
*          Imprimer sur fdout certaines statistiques concernant
*          l'utilisation faite de vmm
*
*auteur J. Caveen  -  novembre 1993
*
*arguments
*                aucun
*
**/
wordint
f77name(vmmdiag)()
{
   int vmmerr();


   if(callc("VMMDIAG")) ;
       
#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmdiag()\n");
#endif
   fprintf(fdout," VMMDIAG-MINIMUM MEMORY REQUIRED : %d WORDS\n",espace_requis_max);
   fprintf(fdout," VMMDIAG-MAXIMUM NUMBER OF MEMORY BLOCKS : %d\n",nbblocks_max);
   fprintf(fdout," VMMDIAG-MAXIMUM NUMBER OF SIMULTANEOUSLY LOCKED FIELDS : %d\n",champs_bloques_max);
   fprintf(fdout," VMMDIAG-NUMBER OF CALLS TO VMMLOD WITH NO LOCKED FIELDS : %d\n",nb_appels_no_lock);
   fprintf(fdout," VMMDIAG-NUMBER OF CALLS TO VMMLOD WITH LOCKED FIELDS : %d\n",nb_appels_lock);
   fprintf(fdout," VMMDIAG-NUMBER OF DISK READS : %d\n",nb_lectures);
   fprintf(fdout," VMMDIAG-NUMBER OF DISK WRITES : %d\n",nb_ecritures);
   return (0) ;
}

/***s/p vmmdmp
*
*objet(vmmdmp)
*
*auteur J. Caveen - juillet 1993
*
*
*argument
*          in mode   mode d'impression
*                    mode = 001   imprime block_table
*                    mode = 010   imprme slice_table
*                    mode = 100   imprime name_table
*
*
**/
wordint 
f77name(vmmdmp)(unsigned wordint *mode)
{
         void imprime_structures();
         int i;
#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmdmp(%d)\n",*mode);
#endif
         for(i = 0; i < 3; i++)
         {
              if((*mode >> i) & 1)   
                       imprime_structures(i);
         }
	 return(0);
}

/***s/p vmmerr
*
*objet(vmmerr)
*     Emettre les messages d'erreur et mettre fin a l'execution si necessaire
*
*auteur J. Caveen  -  juillet 1993
*revision james caveen - mars 1995 
*       changement a la definition de la macro sortir afin de lui
*       passer un argument.
*       Ajout des erreurs pour appel a vmmcpk ou vmmckmx avec champs
*       ayant l'attribut keep_in_core.
*       Elimine l'erreur associee a vmmckmx: CKP_MUSTEXIST_DONE
*
*arguments
*     in   valeur      code de l'erreur
*     in   fonction    nom de la fonction appelante
*
**/
PRIVATE int 
vmmerr(char *fonction,wordint valeur)
{
      switch (valeur)
      {
           case NO_CALL_TO_VMMALLC :                               /*FATAL*/
                    fprintf(fd_err,
                      "ERROR - %s - NO PREVIOUS CALL TO VMMALLC\n",fonction);
                    SORTIR(valeur)
                    break;
           case VMMALLC_ALREADY_CALLED :                           /*FATAL*/
                    fprintf(fd_err,
                      "ERROR - %s - VMMALLC ALREADY CALLED\n",fonction);
                    SORTIR(valeur)
                    break;
           case NO_SPACE_LEFT_FOR_LOAD:                           /*FATAL*/
                    fprintf(fd_err,
                      "ERROR - %s - NO SPACE LEFT FOR LOAD\n",fonction);
                    SORTIR(valeur)
                    break;
           case CONTROL_FILE_ERROR:                               /*FATAL*/
                    fprintf(fd_err,
                      "ERROR - %s - CANNOT OPEN CONTROL FILES\n",fonction);
                    SORTIR(valeur)
                    break;
           case PWD_ALREADY_SET:                                  /*FATAL*/
                    fprintf(fd_err,
                      "ERROR - %s - PASSWORD IS ALREADY SET\n",fonction);
                    SORTIR(valeur)
                    break;
           case BAD_PASSWORD:                                     /*FATAL*/
                    fprintf(fd_err,
                      "ERROR - %s - WRONG PASSWORD\n",fonction);
                    SORTIR(valeur)
                    break;
           case PASSWORD_IS_SET:                                     /*FATAL*/
                    fprintf(fd_err,
                      "ERROR - %s - PASSWORD IS SET\n",fonction);
                    SORTIR(valeur)
                    break;
           case VAR_MUST_EXIST:                                     /*FATAL*/
                    fprintf(fd_err,
                      "ERROR - %s - VARIABLE MUST EXIST FOR A RESTART\n",fonction);
                    SORTIR(valeur)
                    break;
           case CONTROLE_DAMAGE:                                 /*FATAL*/
                    fprintf(fd_err,
                      "ERROR - %s - NAMES-SLICES OR BLOCK-SLICES INCONSISTENCIES\n",
                       fonction);
                    SORTIR(valeur)
                    break;
           case NOT_ENOUGH_MEMORY:                                 /*FATAL*/
                    fprintf(fd_err,
                      "ERROR - %s - CANNOT ALLOCATE MEMORY REQUESTED\n",fonction);
                    SORTIR(valeur)
                    break;
           case TOO_MANY_KEYS:                                 /*FATAL*/
                    fprintf(fd_err,
                      "ERROR - %s - NKEYS > %d, LIMIT EXCEEDED\n",fonction,NCLEMAX);
                    SORTIR(valeur)
                    break;
           case BAD_CKSUM_MODE:                                 /*FATAL*/
                    fprintf(fd_err,
                      "ERROR - %s - BAD MODE FOR CHECK SUM\n",fonction);
                    SORTIR(valeur)
                    break;
           case UNKNOWNVAR:
                    fprintf(fd_err,
                      "ERROR - %s - UNKNOWN VARIABLE\n",fonction);
                    SORTIR(valeur)
                    break;
           case -BADINKEY:
                    fprintf(fd_err,
                      "ERROR - %s - BAD KEY\n",fonction);
                    SORTIR(-valeur)
		    valeur = -valeur;   /* cas du test vmm pour retourner une valeur negative */
                    break;
           case NOT_IN_CORE:
                    fprintf(fd_err,
                      "ERROR - %s - SLICE NOT IN CORE\n",fonction);
                    SORTIR(valeur)
                    break;
           case ALREADY_LOCKED:
                    fprintf(fd_err,
                      "ERROR - %s - SLICE ALREADY LOCKED\n",fonction);
                    SORTIR(valeur)
                    break;
           case BLOCK_DAMAGE:
                    fprintf(fd_err,
                      "ERROR - %s - MEMORY BLOCK DAMAGE\n",fonction);
                    SORTIR(valeur)
                    break;
           case ATTRIBUTS_MODIFIES:
                    fprintf(fd_err,
                      "ERROR - %s - MODIFICATION TO UNCHANGEABLE ATTRIBUTES OF A VARIABLE\n",fonction);
                    SORTIR(valeur)
                    break;
           case SLICE_TOO_BIG:
                    fprintf(fd_err,
                      "ERROR - %s - SLICE LONGER THAN TOTAL MEMORY REQUESTED\n",fonction);
                    SORTIR(valeur)
                    break;
           case CHECKSUM_ERROR:
                    fprintf(fd_err,
                      "ERROR - %s - CHECKSUM MODIFIED FOR AN UNLOCKED FIELD\n",fonction);
                    SORTIR(valeur)
                    break;
           case CHECKSUM_READ_ERROR:
                    fprintf(fd_err,
                      "ERROR - %s - CHECKSUM ERROR WHILE READING FROM FILE\n",fonction);
                    SORTIR(valeur)
                    break;
           case KEEP_IN_CORE_CPK:
                    fprintf(fd_err,
                      "ERROR - %s - CALL TO VMMCPK WITH KEEP IN CORE FIELDS\n",fonction);
                    SORTIR(valeur)
                    break;
           case KEEP_IN_CORE_CKMX:
                    fprintf(fd_err,
                      "ERROR - %s - CALL TO VMMCKMX WITH KEEP IN CORE FIELDS\n",fonction);
                    SORTIR(valeur)
                    break;
           case WAS_ALTERED_RELEASE:
                    fprintf(fd_err,
                      "WARNING - %s - RELEASING A POSSIBLY MODIFIED FIELD\n",fonction);
                    break;
      }   

      return( (int) (-valeur));
}

/***s/p vmmfgt
*
*objet(vmmfgt)
*     Rendre un champ non defini en memoire et sur disque
*
*auteur J. Caveen  -  aout 1993
*
*arguments
*     in   inlkey      liste des clefs de tranches a oublier
*     in   nkey        nombre de clefs dans inlkey
*
**/
f77name(vmmfgt)(complete_key inlkey[], wordint *nkey)
{
       int qvmindex_from_key(), vmmerr();
       int indice, i, iii, bloc_indice;


       if(callc("VMMFGT"))   ;
       


   if(pwd_set)
        return(vmmerr("VMMFGT",PASSWORD_IS_SET));

#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmfgt(");
        for (iii = 0; iii < *nkey; iii++)
            fprintf(fdout,"%d ",inlkey[iii]);

     fprintf(fdout,"],%d)\n",*nkey);
#endif

        for (i=0; i< *nkey; i++)
        {
              indice = qvmindex_from_key(inlkey[i]);
              if (indice < 0)
                    return vmmerr("VMMFGT",indice);
            
              bloc_indice = SLICES[indice].block_table_index ;
              SLICES[indice].info.flags.keep_in_core = 0;
              SLICES[indice].info.flags.is_in_core = 0;
              SLICES[indice].info.flags.in_used = 0;
              SLICES[indice].info.flags.locked = 0;
              SLICES[indice].info.flags.altered = 0;
              SLICES[indice].info.flags.was_altered = 0;
              SLICES[indice].info.flags.disk_image = 0;
              SLICES[indice].checksum = 0;
              SLICES[indice].block_table_index = -1;
              if(bloc_indice != -1)
              {
                 BLOCKS[bloc_indice].info.attributs = 0;
                 BLOCKS[bloc_indice].slice_table_index = -1;
                 BLOCKS[bloc_indice].file_adr = -1;
              }
       }

       return 0;
}


/***s/p vmmget
*
*objet(vmmget)
*         Obtenir le pointeur (l'adresse) d une tranche en memoire
*         Si le champ n'est pas bloque, le bloquer
*
*auteur J. Caveen  -  juillet 1993
*
*revision james caveen - mars 1995
*       bug fix: ajout de la verification de l'attribut keep_in_core
*                avant on pouvait faire un get d'un champ ayant subit un 
*                unload mais etant toujours resident en memoire
*
*
*arguments
*         in  inkey     -  clef du champ desire
*         out pointeur  -  pointeur au champ desire (pointeur fortran)
*         out tablo     -  champ a obtenir (inutile en mode dynamique)
*
**/
wordint 
#if defined (_FLOAT1)
f77name(vmmget)(complete_key  *inkey, wordint *pointeur,wordint *tablo)
#else
f77name(vmmget)(complete_key  *inkey, void **pointeur,wordint *tablo)
#endif
{

          int qvmindex_from_key(), vmmerr(), verbar(),calc_checksum();

          int indice,ier,cks2;
          memint intptr;

          if(callc("VMMGET"))  ;
       


#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmget(%d,%d,%d)\n",(wordint) inkey->clef , pointeur, tablo);
#endif
          indice = qvmindex_from_key(*inkey);

 
          if(indice < 0)
              return vmmerr("VMMGET",indice);

          if( (! SLICES[indice].info.flags.is_in_core)) 
/*          if( (! SLICES[indice].info.flags.is_in_core) || (! SLICES[indice].info.flags.keep_in_core)) **/
          {
              if(debug_mode)
              {
                  fprintf(fdout,"VMM-trace : VARIABLE %s, SLICE %d NOT IN CORE\n",
                         NAMES[SLICES[indice].name_table_index].nom,
                         indice-NAMES[SLICES[indice].name_table_index].major_key+1);
              }
                         
              return vmmerr("VMMGET",NOT_IN_CORE);
          }

          if(! SLICES[indice].info.flags.locked)
          {
                 champs_bloques++;
/*
 *            faire le checksum si requis et le comparer a celui deja dans 
 *            SLICE[indice]
 */
                     if(SLICES[indice].info.flags.do_checksum || checksum_mode)
                     {
                         cks2 = calc_checksum(SLICES[indice].block_table_index);
                         if(cks2 != SLICES[indice].checksum)
                            return(vmmerr("VMMGET",CHECKSUM_ERROR));
                     }
          }

          ier = verbar(SLICES[indice].block_table_index);
          SLICES[indice].info.flags.locked = 1;
          BLOCKS[SLICES[indice].block_table_index].info.flags.locked = 1;
          SLICES[indice].info.flags.was_altered =
               SLICES[indice].info.flags.was_altered ||
               SLICES[indice].info.flags.altered;
          SLICES[indice].info.flags.altered = SLICES[indice].info.flags.save; 
          BLOCKS[SLICES[indice].block_table_index].info.flags.was_altered = 
               SLICES[indice].info.flags.was_altered;
          BLOCKS[SLICES[indice].block_table_index].info.flags.altered = 
	       SLICES[indice].info.flags.save;


#if defined (_FLOAT1)
          intptr = (memint) BLOCKS[SLICES[indice].block_table_index].memadr;
#if defined (ALL64)
          *pointeur = (wordint) (intptr >> 3);
#else
          *pointeur = (wordint)  (intptr >> 2);
#endif
#else
          *pointeur = (void *) BLOCKS[SLICES[indice].block_table_index].memadr;
#endif
          
          champs_bloques_max = champs_bloques_max > champs_bloques ?
                               champs_bloques_max : champs_bloques ;
          return 0;
}

/***s/p vmmhpa
*
*objet(vmmhpa)
*     Allocation d'espace temporaire a partir du bloc de memoire gere par
*     le systeme VMM
*
*auteur M. Lepine  -  juin 1993
*
*argument
*     out  ptr    pointeur au bloc de memoire demande
*     in   memry  quantite d'espace demande (en mots)
*     in   mode   si mode=8, forcer une allocation en mots de 8 bytes
*
**/
wordint
#if defined (_FLOAT1)
f77name(vmmhpa)(wordint *ptr,wordint *memry,wordint *mode)
#else
f77name(vmmhpa)(void **ptr,wordint *memry,wordint *mode)
#endif
{

   int vmmerr();
   int nbytes;
   memint lptr;
   wordint *pointeur;

   if( callc("VMMHPA")) ; 
       
#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmhpa(%d,%d,%d)\n",ptr,*memry,*mode);
#endif

   nbytes = sizeof(wordint);
   nbytes = (int) ((nbytes ==4 && *mode ==8) ? nbytes * *memry * 2 : nbytes * *memry);
/*
   if (free_space - nbytes/sizeof(wordint) < 0)
      fprintf(fd_err,"vmmhpa debug : allocation supplementaire (malloc)\n");
   else
      free_space -=  nbytes/sizeof(wordint);
*/
   pointeur = (wordint *) malloc(nbytes); 

   if (pointeur == (wordint *) NULL)
        return(vmmerr("VMMHPA",NOT_ENOUGH_MEMORY));
   else
   {
#if defined (_FLOAT1)
      lptr = (memint) pointeur;
#if defined (ALL64)
      *ptr = (wordint) (lptr >> 3);
#else
      *ptr = (wordint) (lptr >> 2);
#endif
#else
      *ptr = pointeur;
#endif
      return(0);
    }
   }

/***s/p vmmhpd
*
*objet(vmmhpd)
*     Desallocation d'espace temporaire a partir du bloc de memoire gere par
*     le systeme VMM
*
*auteur M. Lepine  -  aout 1993
*
*argument
*     in   ptr    pointeur au bloc de memoire
*
**/
wordint
#if defined(_FLOAT1)
f77name(vmmhpd)(wordint *ptr)
#else
f77name(vmmhpd)(void **ptr)
#endif
{

   memint lptr;
   wordint *pointeur;

   int vmmerr();

   if( callc("VMMHPD")) ; 
       


#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmhpd(%d)\n",ptr);
#endif

#if defined (_FLOAT1)
   lptr = (memint) *ptr;
#if defined (ALL64)
   pointeur = (wordint *) (lptr << 3);
#else
   pointeur = (wordint *) (lptr << 2);
#endif
#else
   pointeur =  *ptr;
#endif
   free(pointeur);
   return(0);
   }

/***s/p vmmint
*
*objet(vmmint)
*          S'assurer que les structures names, slices et blocks
*          sont en bon etat.
*
*          pour chaque bloc -  on s'assure que son slice_table_index
*                              est valable.
*          pour chaque names-  on s'assure qu'il y a le bon nombre
*                              de slices et que les name_table_index
*                              sont valables.
*
*          En cas d'erreur, vmmint arrete l'execution
*          du programme.
*
*auteur J. Caveen  -  aout 1993
*
*/
int f77name(vmmint)()
{

       int vmmerr(), verbar();
       int i,j,ier,ntranches;

#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmint()\n");
#endif
/*
 *     boucle sur les blocs
 */
       ier = 0;
       for(i = 0; i<nbblocks; i++)
       {
         if(BLOCKS[i].slice_table_index != -1)
         {
	    ier = verbar(i);
            if(SLICES[BLOCKS[i].slice_table_index].block_table_index != i)
            {
                fprintf(fdout," ERROR - INDEX MISMATCH BLOCKS[%d].slice_table_index = %d, SLICES[%d].block_table_index = %d\n",
                  i,BLOCKS[i].slice_table_index,
                  BLOCKS[i].slice_table_index,
                  SLICES[BLOCKS[i].slice_table_index].block_table_index); 
                    ier--;
            }
         }
       }

/*
 *     boucle sur les names et slices 
 */
       for(i=0; i<nbvar; i++)
       {
            ntranches = NAMES[i].nslice;
            for(j=0; j<ntranches; j++)
            {
                 if(SLICES[NAMES[i].major_key+j].name_table_index != i)
                 {
                    fprintf(fdout,
                      " ERROR - INDEX MISMATCH SLICES[%d].name_table_index = %d for NAMES[%d] (%s)\n",
                      NAMES[i].major_key+j,
                      SLICES[NAMES[i].major_key+j].name_table_index,
                      i,NAMES[i].nom);
                    ier--;
                 } 
                 if(SLICES[NAMES[i].major_key+j].info.flags.class !=
                                                        NAMES[i].class)
                 {
                    fprintf(fdout,
                       " ERROR - CLASS MISMATCH SLICES[%d].class = %d, NAMES[%d].class = %d\n",
                       NAMES[i].major_key+j,
                       SLICES[NAMES[i].major_key+j].info.flags.class,
                       i,NAMES[i].class );
                    ier--;
                 } 
                 
             }
                     
       }
       
       return (ier == 0 ? ier : vmmerr("VMMINT",CONTROLE_DAMAGE));
}

/***s/p vmmlck
*
*objet(vmmlck)
*      Bloquer un ou plusieurs champs en memoire.
*      Un champ bloque ne peut etre deplace en memoire.
*      Un champ ne peut etre bloque que s'il est resident en memoire
*      et n'est pas deja bloque.
*
*auteur J. Caveen  -  juillet 1993
*
*
*arguments
*         in inlkey  -  liste de clefs des champs a bloquer
*         in nkey    -  nombre de clefs dans inlkey
*
**/
wordint
f77name(vmmlck)(complete_key inlkey[], wordint *nkey)
{
       int qvmindex_from_key(), vmmerr(), verbar(),calc_checksum();
       int indice, i,iii,ier;

       if (callc("VMMLCK"))   ;
       


   if(pwd_set)
        return(vmmerr("VMMLCK",PASSWORD_IS_SET));


#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmlck([");
        for (iii = 0; iii < *nkey; iii++)
            fprintf(fdout,"%d ",(wordint) inlkey[iii].clef);

     fprintf(fdout,"],%d)\n",*nkey);
#endif
/*
 *     boucle sur toutes les clefs
 */

       for (i=0; i< *nkey; i++)
       {
              indice = qvmindex_from_key(inlkey[i]);
              if (indice < 0)
                    return vmmerr("VMMLCK",indice);
            
              if( ! SLICES[indice].info.flags.is_in_core)
                 return vmmerr("VMMLCK",NOT_IN_CORE); 

              if(SLICES[indice].info.flags.locked)
                 return vmmerr("VMMLCK",ALREADY_LOCKED);
              ier = verbar(SLICES[indice].block_table_index);   
              SLICES[indice].info.flags.locked = 1;
              BLOCKS[SLICES[indice].block_table_index].info.flags.locked = 1;
              champs_bloques++;
              if ((SLICES[indice].info.flags.traced) || debug_mode)
                fprintf(fdout,"VMM trace: blocage de %s tranche %d\n",
                NAMES[SLICES[indice].name_table_index].nom,
                indice - NAMES[SLICES[indice].name_table_index].major_key + 1);

              if (SLICES[indice].info.flags.do_checksum || checksum_mode)
                    SLICES[indice].checksum =
                       calc_checksum(SLICES[indice].block_table_index);
       }
     
       champs_bloques_max = champs_bloques_max > champs_bloques ?
                            champs_bloques_max : champs_bloques ;
       return 0;
}

/***s/p vmmlod
*
*objet(vmmlod)
*      Interface usager pour faire les chargements de
*      champs en memoire.  Vmmlod appelle qvmlod.
*
*auteur J. Caveen - octobre 1993
*
*revision : J. Caveen - juillet 1994
*                       On met l'attribut keep_in_core a 1
*                       pour les champs deja en memoire dans le
*                       cas ou il y a  des champs bloques.
*                       Ceci evite l'ejection inutile de blocs. 
*                       
*                       
*                       
*arguments
*         in inlkey  -  liste de clefs des champs a charger
*         in nkey    -  nombre de clefs dans inlkey
*
**/
wordint
f77name(vmmlod)(complete_key inlkey[], wordint *nkey)
{
        wordint qvmlod();
	int qvm_index_frim_key();
        wordint ier, i,iii, un = 1;
	int clef;

         

#if defined CALL_SEQ
     if(*nkey > NCLEMAX)
          fprintf(fdout,"CALL- vmmlod(nkey > %d)\n",NCLEMAX);
     else
     {
         fprintf(fdout,
          "CALL- vmmlod([");
            for (iii = 0; iii < *nkey; iii++)
                fprintf(fdout,"%d ",(wordint) inlkey[iii].clef);

         fprintf(fdout,"],%d)\n",*nkey);
         fprintf(fdout," adresse de inlkey = %x\n",&inlkey[0]);
     }
#endif

        if(champs_bloques)
        {
             nb_appels_lock++;
	     
	     for(i = 0;  i < *nkey; i++)
	     {
		clef = qvmindex_from_key(inlkey[i]);
		if(SLICES[clef].info.flags.is_in_core)
		{
		   SLICES[clef].info.flags.keep_in_core = 1;
		   BLOCKS[SLICES[clef].block_table_index].info.flags.keep_in_core = 1;
		}
	     }
	       
             for( i = 0; i < *nkey; i++)
             {
                  ier = qvmlod(&inlkey[i],&un);
/*ETG                  if(ier < 0)
                         return(ier); */

             }
        }
        else
        {
             nb_appels_no_lock++;
/*ETG             return( ier = qvmlod(inlkey,nkey)); */
             ier = qvmlod(inlkey,nkey);
        }
        return(ier);
}
                  


/***s/p vmmlse
*
*objet(vmmlse)
*     Retourne a l'usager la quantite de memoire encore disponible dans le
*     systeme de gestion de memoire virtuelle
*
*auteur M. Lepine  -  juin 1993
*
*
**/
wordint
f77name(vmmlse)()
{
   int vmmerr();
   int biggest_free_block_size, i;

   if( callc("VMMLSE")) ;
       


   if(pwd_set)
        return(vmmerr("VMMLSE",PASSWORD_IS_SET));

#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmlse()\n");
#endif

   biggest_free_block_size = 0;
   for (i=0; i < nbblocks; i++)
      if (! BLOCKS[i].info.flags.in_used)
	 if (BLOCKS[i].size > biggest_free_block_size)
	    biggest_free_block_size = BLOCKS[i].size;
   return(biggest_free_block_size);
   }


/***s/p vmmpak
*
*objet(vmmpak)
*     Reoptimisation (compression) de l'espace memoire du systeme de gestion de
*     memoire virtuelle
*
*auteur M. Lepine  -  juin 1993
*
**/
wordint
f77name(vmmpak)()
{
   int vmmerr(), pack_blocks();
   int biggest, ind;

   if( callc("VMMPAK")) ;
       


   if(pwd_set)
        return(vmmerr("VMMPAK",PASSWORD_IS_SET));

#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmpak()\n");
#endif
   
   ind = pack_blocks(&biggest);

   return(0);
   }

      
/***s/p vmmpwd
*
*objet(vmmpwd)
*       Fonction servant a verrouiller ou a deverrouiller le systeme
*       de gestion de memoire virtuelle au moyen d'un mot de passe.
*       Lorsque le mot de passe est cree, la seule fonction
*       du systeme de gestion de memoire virtuelle que l'usager
*       peut appeler est vmmget.
*       On deverrouille le systeme en faisant un autre appel a vmmpwd
*       en lui passant le mot de passe ayant servi au verrouillage
*
*auteur J. Caveen  -  juillet 1993
*
*
*arguments
*         in mot_passe     mot de passe servant a verrouiller ou deverrouiller
*         in mode          si mode = 0, on verrouille
*                             mode = 1, on deverrouille
*
**/
wordint
f77name(vmmpwd)(wordint *mot_passe, wordint *mode)
{


         int vmmerr();

         if(callc("VMMPWD"))   ;
       


        if(*mode == 0)
        {
           if(pwd_set)
               vmmerr("VMMPWD",PWD_ALREADY_SET);

#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmpwd(%d,%d)\n",*mot_passe, *mode);
#endif
           pwd_set = 1;
           mot_de_passe = *mot_passe; 
          
        }
        else
        {
           if(*mot_passe != mot_de_passe)
               vmmerr("VMMPWD",BAD_PASSWORD);

           pwd_set = mot_de_passe = 0;
         }

         return (0);
}

/***s/p vmmrls
*
*objet(vmmrls)
*     ejection d'un champ sans sauvegarde 
*
*auteur J. Caveen  -  juillet 1993
*
*revision j. caveen - fevrier 1995
*                     enleve le mode inlkey[0] = -1
*                     pour faire un release de toutes les clefs
*                     Cette methode a des effets de bords
*                     pour le moins dangeureux
*
*arguments
*     in   inlkey      liste des clefs de tranches a expulser
*     in   nkey        nombre de clefs dans inlkey
*
**/
f77name(vmmrls)(complete_key inlkey[], wordint *nkey)
{
       int qvmindex_from_key(), vmmerr(),eject_block();
       int indice, i,iii, bloc_indice, ier;


       if(callc("VMMRLS"))   ;
       


   if(pwd_set)
        return(vmmerr("VMMRLS",PASSWORD_IS_SET));

#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmrls([");
        for (iii = 0; iii < *nkey; iii++)
            fprintf(fdout,"%d ",(wordint) inlkey[iii].clef);

       fprintf(fdout,"],%d)\n",*nkey);
       fprintf(fdout," adresse de inlkey = %x\n",&inlkey[0]);
       
#endif
       for (i=0; i< *nkey; i++)
       {
	  indice = qvmindex_from_key(inlkey[i]);
	  if (indice < 0)
	  return vmmerr("VMMRLS",indice);
	  
	  if( SLICES[indice].info.flags.was_altered)
	  ier = vmmerr("VMMRLS",WAS_ALTERED_RELEASE); 
	  
	  bloc_indice = SLICES[indice].block_table_index ;
	  if(bloc_indice != -1)
	  {
	     SLICES[indice].checksum = 0;
	     if(BLOCKS[bloc_indice].info.flags.locked)
	     champs_bloques--;
	     ier = eject_block(bloc_indice,0,0);
	  }
       }
       
       return 0;
}

/***s/p vmmrnm
*
*objet(vmmrnm)
*    Changer le nom d'une variable dans les tables du systeme de 
*    gestion de memoire 
*
*auteur J. Caveen  -  juin 1993
*
*argument
*     in   oldkey   clef pointant a une des slices de la variable
*     in   newname  nouveau nom a remiser dans la table names
*     in   l1       longueur de newname
*
**/
wordint
f77name(vmmrnm)(complete_key *oldkey,char *newname, F2Cl l1)
{
     int vmmerr();
     int qvmindex_from_key();
     int slice_index,name_index;
     char innewname[9];
     int i;

     if (callc("VMMRNM"))   ;
       


   if(pwd_set)
        return(vmmerr("VMMRNM",PASSWORD_IS_SET));

   strncpy(innewname,newname,l1);

/*
 * on remplie le reste de innewname avec des blancs
 */
   for (i = l1; i < 8; i++)
          innewname[i] = ' ';

   innewname[8] = '\0';



#if defined CALL_SEQ
     fprintf(fdout," longueur de newname: %d\n",l1);
     fprintf(fdout,
      "CALL- vmmrnm(%d,%s)\n",(wordint) oldkey->clef,innewname);
#endif
/*
 *   trouver la slice 
 */
     slice_index = qvmindex_from_key(*oldkey);

     if(slice_index < 0)
            return vmmerr("VMMRNM",slice_index);
/*
 *   trouver l'entree dans la table name et changer le nom
 */
     name_index = SLICES[slice_index].name_table_index;
     strcpy(NAMES[name_index].nom,innewname);

     
     return 0;
}


/***s/p vmmsav
*
*objet(vmmsav)
*     Sauver une ou plusieurs variables sur disque 
*
*auteur J. Caveen  -  juillet 1993
*
*arguments
*     in   inlkey      liste des clefs de tranches a expulser
*     in   nkey        nombre de clefs dans inlkey
*
**/
wordint 
f77name(vmmsav)(complete_key inlkey[], wordint *nkey)
{
         void ecrit_bloc(),reserve_disk_space();
         int qvmindex_from_key(), vmmerr();
         int indice,  i,iii;

         if(callc("VMMSAV"))  ;
       


         if(pwd_set)
             return(vmmerr("VMMSAV",PASSWORD_IS_SET));

#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmsav([");
        for (iii = 0; iii < *nkey; iii++)
            fprintf(fdout,"%d ",(wordint) inlkey[iii].clef);

     fprintf(fdout,"],%d)\n",*nkey);
     fprintf(fdout," adresse de inlkey = %x\n",&inlkey[0]);

#endif

         if(inlkey[0].clef == -1)
         {
             for (i = 0; i< nbblocks; i++)
             {
                if(BLOCKS[i].info.flags.in_used)
                {
                   if(BLOCKS[i].info.flags.save && (BLOCKS[i].info.flags.altered || 
                      BLOCKS[i].info.flags.was_altered))
                   {
                      if(BLOCKS[i].file_adr == -1)
                          reserve_disk_space(i);
                            ecrit_bloc(i,BLOCKS[i].info.flags.class,BLOCKS[i].memadr,
                                       BLOCKS[i].file_adr,BLOCKS[i].size);
                   }
                }
             }
         }
         else
         {
             for (i = 0; i < *nkey; i++)
             {
                indice = qvmindex_from_key(inlkey[i]);
                if(indice < 0)
                      return vmmerr("VMMSAV",indice);

                if(BLOCKS[SLICES[indice].block_table_index].info.flags.in_used)
                {
                   if(BLOCKS[SLICES[indice].block_table_index].info.flags.save &&
                      (BLOCKS[SLICES[indice].block_table_index].info.flags.altered ||
                      BLOCKS[SLICES[indice].block_table_index].info.flags.was_altered))
                   {
                       if(BLOCKS[SLICES[indice].block_table_index].file_adr == -1)
                            reserve_disk_space(SLICES[indice].block_table_index);
                       ecrit_bloc(SLICES[indice].block_table_index,
                             BLOCKS[SLICES[indice].block_table_index].info.flags.class,
                             BLOCKS[SLICES[indice].block_table_index].memadr,
                             BLOCKS[SLICES[indice].block_table_index].file_adr,
                             BLOCKS[SLICES[indice].block_table_index].size);
                   }
                }
              }
          }

          return 0;
}




                

/***s/p vmmuld
*
*objet(vmmuld)
*     Rendre un champ ejectable apres sauvegarde 
*     vmmuld ne fait ni la sauvegarde ni l'ejection
*     sauf si le champ est nosave
*
*auteur J. Caveen  -  juillet 1993
*
*arguments
*     in   inlkey      liste des clefs de tranches a expulser
*     in   nkey        nombre de clefs dans inlkey
*
**/
wordint
f77name(vmmuld)(complete_key inlkey[], wordint *nkey)
{
       int qvmindex_from_key(), vmmerr(), eject_block(), verbar();
       int calc_checksum();
       int indice, i,iii, bloc_indice,ier;


       if(callc("VMMULD"))   ;
       


   if(pwd_set)
        return(vmmerr("VMMULD",PASSWORD_IS_SET));

#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmuld([");
        for (iii = 0; iii < *nkey; iii++)
            fprintf(fdout,"%d ",(wordint) inlkey[iii].clef);

     fprintf(fdout,"],%d)\n",*nkey);
     fprintf(fdout," adresse de inlkey = %x\n",&inlkey[0]);

#endif

/*  
 *     si inlkey[0] = -1, on ejecte toutes les slices
 */
 
/*
 *          boucle sur toutes les clefs
 */
       if ( inlkey[0].clef == -1)
       {
            for (i=0; i< nbblocks; i++)
            {
              if(BLOCKS[i].info.flags.in_used)
              {
                if ((BLOCKS[i].info.flags.traced) || debug_mode)
                    fprintf(fdout,"VMM trace: vmmuld du bloc %d variable %s tranche %d\n",
                          i,NAMES[SLICES[BLOCKS[i].slice_table_index].name_table_index].nom,
                           BLOCKS[i].slice_table_index - NAMES[SLICES[BLOCKS[i].slice_table_index].name_table_index].major_key + 1 );
                if(BLOCKS[i].slice_table_index != -1)
                {
                   SLICES[ BLOCKS[i].slice_table_index].info.flags.keep_in_core = 0;
                   SLICES[ BLOCKS[i].slice_table_index].info.flags.locked = 0;
                   if(SLICES[ BLOCKS[i].slice_table_index].info.flags.do_checksum || checksum_mode)
                      SLICES[ BLOCKS[i].slice_table_index].checksum = 
                         calc_checksum(i);
                }
                if (BLOCKS[i].info.flags.save)
                {
                   ier = verbar(i);

                   BLOCKS[i].info.flags.keep_in_core = 0;
                   BLOCKS[i].info.flags.locked = 0;
                }
                else
                   ier = eject_block(i,0,0);

              }
            }
            champs_bloques = 0;
       }
       else
       {
            for (i=0; i< *nkey; i++)
            {
              indice = qvmindex_from_key(inlkey[i]);
              if (indice < 0)
                    return vmmerr("VMMULD",indice);
            
              SLICES[indice].info.flags.keep_in_core = 0;
              SLICES[indice].info.flags.locked = 0;
              bloc_indice = SLICES[indice].block_table_index ;

              if(bloc_indice != -1)
              {
                 if(SLICES[indice].info.flags.do_checksum || checksum_mode)
                       SLICES[indice].checksum = calc_checksum(bloc_indice);

                 if(SLICES[indice].info.flags.save)
                 {
                     ier = verbar(bloc_indice);
                     if ((BLOCKS[bloc_indice].info.flags.traced) || debug_mode)
                         fprintf(fdout,"VMM trace: vmmuld du bloc %d variable %s tranche %d\n",
                     bloc_indice,NAMES[SLICES[indice].name_table_index].nom,
                     indice - NAMES[SLICES[indice].name_table_index].major_key + 1);

                     BLOCKS[bloc_indice].info.flags.keep_in_core = 0;
                     if(BLOCKS[bloc_indice].info.flags.locked)
                             champs_bloques--;
                     BLOCKS[bloc_indice].info.flags.locked = 0;
                  } 
                  else
                  {
                      if(BLOCKS[bloc_indice].info.flags.locked)
                            champs_bloques--;
     
                      ier = eject_block(bloc_indice,0,0);
                  }

              }
          }
       }
       return 0;
}

/***s/p vmmulk
*
*objet(vmmulk)
*      Debloquer un ou plusieurs champs en memoire.
*
*auteur J. Caveen  -  juillet 1993
*
*
*arguments
*         in inlkey  -  liste de clefs des champs a bloquer
*                       si inlkey = -1, on debloque toutes les variables
*         in nkey    -  nombre de clefs dans inlkey
*
**/
wordint
f77name(vmmulk)(complete_key inlkey[], wordint *nkey)
{
       int qvmindex_from_key(), vmmerr(),verbar(),calc_checksum();
       int indice, i,iii,ier;


       if(callc("VMMULK"))   ;
       


   if(pwd_set)
        return(vmmerr("VMMULK",PASSWORD_IS_SET));

#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmulk([");
        for (iii = 0; iii < *nkey; iii++)
            fprintf(fdout,"%d ",(wordint) inlkey[iii].clef);

     fprintf(fdout,"],%d)\n",*nkey);
#endif

/*  
 *     si inlkey[0] = -1, on boucle sur les blocs pour mettre locked = 0 
 *     sinon, on boucle sur les slices
 */
 
       if ( inlkey[0].clef == -1)
       {
            for (i=0; i< nbblocks; i++)
            {
              ier = verbar(i);
	      indice = BLOCKS[i].slice_table_index;
              BLOCKS[i].info.flags.locked = 0;
              if(indice != -1)
              {
                  SLICES[indice].info.flags.locked = 0;

                  if ((SLICES[indice].info.flags.traced) || debug_mode)
                       fprintf(fdout,"VMM trace: deblocage de %s tranche %d\n",
                          NAMES[SLICES[indice].name_table_index].nom,
                          indice - NAMES[SLICES[indice].name_table_index].major_key + 1);

                  if (SLICES[indice].info.flags.do_checksum || checksum_mode)
                        SLICES[indice].checksum = calc_checksum(i);
              }
           }
           champs_bloques = 0;
       }
       else
       {
            for (i = 0; i < *nkey; i++)
            {
              indice = qvmindex_from_key(inlkey[i]);
              if (indice < 0)
                    return vmmerr("VMMULK",indice);
            
              SLICES[indice].info.flags.locked = 0;
              if(SLICES[indice].block_table_index != -1)
              {
                ier = verbar(SLICES[indice].block_table_index);
                if(BLOCKS[SLICES[indice].block_table_index].info.flags.locked)
                    champs_bloques--;

                 BLOCKS[SLICES[indice].block_table_index].info.flags.locked = 0;

                 if(SLICES[indice].info.flags.do_checksum || checksum_mode)
                       SLICES[indice].checksum = calc_checksum(SLICES[indice].block_table_index);
              }

              if ((SLICES[indice].info.flags.traced) || debug_mode)
                   fprintf(fdout,"VMM trace: deblocage de %s tranche %d\n",
                   NAMES[SLICES[indice].name_table_index].nom,
                   indice - NAMES[SLICES[indice].name_table_index].major_key + 1);

            }
        }

       return 0;
}

/***s/p vmmuln
*
*objet(vmmuln)
*     Rendre un champ ejectable sans sauvegarde 
*     vmmuln ne fait pas l'ejection
*     sauf si le bloc est nosave
*
*auteur J. Caveen  -  juillet 1993
*
*revision: j. caveen - fevrier 1995 
*                      enlever la possibilite de passer
*                      inlk[0] = -1 pour faire un uln
*                      tous les blocs memoires; (cette
*                      methode a des effets de bord
*                      pour le moins dangeureux)
*                      
*arguments
*     in   inlkey      liste des clefs de tranches a expulser
*     in   nkey        nombre de clefs dans inlkey
*
**/
f77name(vmmuln)(complete_key inlkey[], wordint *nkey)
{
       int qvmindex_from_key(), vmmerr(), eject_block(),verbar();
       int calc_checksum();
       int indice, i,iii, bloc_indice,ier;


       if(callc("VMMULN"))   ;
       


   if(pwd_set)
        return(vmmerr("VMMULN",PASSWORD_IS_SET));

#if defined CALL_SEQ
     fprintf(fdout,
      "CALL- vmmuln([");
        for (iii = 0; iii < *nkey; iii++)
            fprintf(fdout,"%d ",(wordint) inlkey[iii].clef);

     fprintf(fdout,"],%d)\n",*nkey);
#endif
       
       for (i=0; i< *nkey; i++)
       {
	  indice = qvmindex_from_key(inlkey[i]);
	  if (indice < 0)
	  return vmmerr("VMMULN",indice);
	  
	  SLICES[indice].info.flags.keep_in_core = 0;
	  SLICES[indice].info.flags.locked = 0;
	  SLICES[indice].info.flags.altered =
	         SLICES[indice].info.flags.was_altered;
	  bloc_indice = SLICES[indice].block_table_index ;
	  if(bloc_indice != -1)
	  {
	     if(SLICES[indice].info.flags.do_checksum || checksum_mode)
	     SLICES[indice].checksum = calc_checksum(bloc_indice);
	     
	     if(SLICES[indice].info.flags.save) 
	     {
		ier = verbar(bloc_indice);
		if(BLOCKS[bloc_indice].info.flags.locked)
		champs_bloques--;
		BLOCKS[bloc_indice].info.flags.keep_in_core = 0;
		BLOCKS[bloc_indice].info.flags.locked = 0;
		BLOCKS[bloc_indice].info.flags.altered =
		BLOCKS[bloc_indice].info.flags.was_altered;
		
		if ((BLOCKS[bloc_indice].info.flags.traced) || debug_mode)
		{
		   fprintf(fdout,"VMM trace: vmmuln du bloc %d variable %s tranche %d\n",
			   bloc_indice,NAMES[SLICES[indice].name_table_index].nom,
			   indice - NAMES[SLICES[indice].name_table_index].major_key + 1);
		   if(BLOCKS[bloc_indice].info.flags.altered)
		   fprintf(fdout,"           Block will be saved upon ejection\n");
		   else
                   fprintf(fdout,"           Block will not be saved upon ejection\n");
		}
	     }
	     else
	     {
		if(BLOCKS[bloc_indice].info.flags.locked)
		champs_bloques--;
		eject_block(bloc_indice,0,0);
	     }
	  }
	  else
	  {
	     fprintf(fdout,"VMM trace: vmmuln  variable %s tranche %d pas en memoire\n",
		     NAMES[SLICES[indice].name_table_index].nom,
		     indice - NAMES[SLICES[indice].name_table_index].major_key + 1);
	  }
       }
       return 0;
}



/***s/p vmmwho
*
*objet(vmmwho)
*          Determiner a quelle variable correspond l'adresse
*          recherchee
*
*auteur M. Lepine  -  avril 1994
*
*arguments
*         in  adr  -  adresse recherchee
*
**/
wordint
f77name(vmmwho)(wordint *adr)
{
   int i,sind,nind;

   if ((adr < BLOCKS[0].memadr) || (adr > BLOCKS[nbblocks-1].memadr)) {
     fprintf(fdout,"VMMWHO address %#x not in block table\n",adr);
     return(-1);
     }
   for (i=0; i < nbblocks; i++) {
      if ((adr >= BLOCKS[i].memadr) && 
          (adr <= BLOCKS[i].memadr+BLOCKS[i].size-1)) {
        sind = BLOCKS[i].slice_table_index;
	nind = SLICES[sind].name_table_index;
	fprintf(fdout,
                "  VMMWHO\n\t address in block #%d, variable %s, slice %d\n\n",
		i,NAMES[nind].nom,sind - NAMES[nind].major_key + 1);
	fprintf(fdout,"\t          address = %#x\n",adr);
	fprintf(fdout,"\t BLOCKS[%d].memadr = %#x\n",i,BLOCKS[i].memadr);
	fprintf(fdout,"\t          indice  = %d\n\n",1+(adr-BLOCKS[i].memadr));
        fprintf(fdout,"\t keep_in_core       : %d\n",SLICES[sind].info.flags.keep_in_core);
        fprintf(fdout,"\t is_in_core         : %d\n",SLICES[sind].info.flags.is_in_core);
        fprintf(fdout,"\t in_used            : %d\n",SLICES[sind].info.flags.in_used);
        fprintf(fdout,"\t locked             : %d\n",SLICES[sind].info.flags.locked);

        fprintf(fdout,"\t save               : %d\n",SLICES[sind].info.flags.save);
        fprintf(fdout,"\t altered            : %d\n",SLICES[sind].info.flags.altered);
        fprintf(fdout,"\t was_altered        : %d\n",SLICES[sind].info.flags.was_altered);
        fprintf(fdout,"\t traced             : %d\n",SLICES[sind].info.flags.traced);
        fprintf(fdout,"\t hpa_alloc          : %d\n",SLICES[sind].info.flags.hpa_alloc);
        fprintf(fdout,"\t disk_image         : %d\n",SLICES[sind].info.flags.disk_image);
        fprintf(fdout,"\t size8              : %d\n",SLICES[sind].info.flags.size8);
        fprintf(fdout,"\t must_exist         : %d\n",SLICES[sind].info.flags.must_exist);
        fprintf(fdout,"\t class              : %d\n",SLICES[sind].info.flags.class);
        fprintf(fdout,"\t weight             : %d\n",SLICES[sind].info.flags.weight);
        fprintf(fdout,"\t init               : %d\n",SLICES[sind].info.flags.init);
        fprintf(fdout,"\t block_table_index  : %d\n",SLICES[sind].block_table_index);
        fprintf(fdout,"\t name_table_index   : %d\n",SLICES[sind].name_table_index);
	return(0);
	}
      }
  return 0; /*CHC/NRC*/
}

/***s/p vmmckmx
*
*objet(vmmckmx)
*    Faire un check-point des variables MUSTEXIST (sauver les variables dans les fichiers Vmm_0n
*    et sauver les tables de controle dans Vmm_controle
*    On ferme les fichiers pour vider tous les tampons
*
*auteur J. Caveen  -  juillet 1993
*
*revision james caveen - mars 1995
*         forcer l'ejection des blocs qui n'ont pas l'attribut mustexist
*         et indiquer que toutes les variables qui ne sont pas mustexist
*         n'ont pas d'image disque.
*         elimination du marqueur ckp_mustexist_done qui servait a prevenir
*         la continuation de l'execution apres vmmckmx.  On peut maintenant
*         faire plusieurs appels a vmmckmx au cours d'une meme execution.
*

wordint f77name(vmmckmx)()
{

   void ecrit_vmm_controle(), ecrit_bloc();
   void reserve_disk_space();
   int vmmerr(), strfind(), eject_block();
   void  ouvre_ou_ferme_controle();
   
   int i,pos;
   
   if(callc("VMMCKMX"))  ;
   
   
   
   if(pwd_set)
   return(vmmerr("VMMCKMX",PASSWORD_IS_SET));
   
#if defined CALL_SEQ
   fprintf(fdout,"CALL- vmmckmx()\n");
#endif
**/   
/*
 * On ecrit tous les blocs memoires qui ont l'attribut mustexist   
 * et on ejecte tous les autres sans sauvegarde
   for (i = 0; i< nbblocks; i++)
   {
      if(BLOCKS[i].info.flags.in_used)
      {
	 if(BLOCKS[i].info.flags.keep_in_core)
	 return(vmmerr("VMMCKMX",KEEP_IN_CORE_CKMX));

	 if(BLOCKS[i].info.flags.must_exist && (BLOCKS[i].info.flags.altered ||
					     BLOCKS[i].info.flags.was_altered))
	 {
	    if(BLOCKS[i].file_adr == -1)
	               reserve_disk_space(i);

	    ecrit_bloc(i,BLOCKS[i].info.flags.class,BLOCKS[i].memadr,
		       BLOCKS[i].file_adr,BLOCKS[i].size);
	 }
	 else
	    eject_block(i,0,0);
      }
   }

 * On indique que toutes les tranches qui ne 
 * sont pas mustexist n'ont pas d'image disque

   for (i = 0; i < nbslices; i++)
      if(! SLICES[i].info.flags.must_exist)
	 SLICES[i].info.flags.disk_image = 0;


   
   ecrit_vmm_controle();
   ouvre_ou_ferme_controle(0,0,"VMMCKMX");
   return 0;
}

**/
/***s/p vmmckmx
*
*objet(vmmckmx)
*    Faire un check-point des variables MUSTEXIST (sauver les variables dans les fichiers Vmm_0n
*    et sauver les tables de controle dans Vmm_controle
*    On ferme les fichiers pour vider tous les tampons
*
*auteur J. Caveen  -  juillet 1993
*
**/
wordint f77name(vmmckmx)()
{

   void ecrit_vmm_controle(), ecrit_bloc();
 /*   
    ouvre_ou_ferme_controle(); 
*/
   void reserve_disk_space();
   int vmmerr(), strfind(); 
   
   int i,pos;
   
   if(callc("VMMCKMX"))  ;
   
   
   if(pwd_set)
   return(vmmerr("VMMCKMX",PASSWORD_IS_SET));
   
#if defined CALL_SEQ
   fprintf(fdout,"CALL- vmmckmx()\n");
#endif
   
   
   for (i = 0; i< nbblocks; i++)
   {
      if(BLOCKS[i].info.flags.in_used)
      {
	 if(BLOCKS[i].info.flags.must_exist && (BLOCKS[i].info.flags.altered ||
						BLOCKS[i].info.flags.was_altered))
	 {
	    if(BLOCKS[i].file_adr == -1)
	               reserve_disk_space(i);

	    ecrit_bloc(i,BLOCKS[i].info.flags.class,BLOCKS[i].memadr,
		       BLOCKS[i].file_adr,BLOCKS[i].size);
	 }
	 else if(BLOCKS[i].info.flags.save && !(BLOCKS[i].info.flags.must_exist))
	 {
	    BLOCKS[i].info.flags.altered =
	    BLOCKS[i].info.flags.was_altered =
	    BLOCKS[i].info.flags.disk_image = 0;
	    SLICES[BLOCKS[i].slice_table_index].info.flags.altered =
	    SLICES[BLOCKS[i].slice_table_index].info.flags.was_altered =
	    SLICES[BLOCKS[i].slice_table_index].info.flags.disk_image = 0;
	 }
      }
   }
   
/**   ckp_mustexist_done = 1; **/
   
   ecrit_vmm_controle();
   ouvre_ou_ferme_controle(0,0,"VMMCKMX");
   return 0;
}

/***s/p vmmdel
*
*objet(vmmdel)
*     vmmdel detruit le fichier de controle si il existe ainsi que les fichiers
*     de sauvegarde
*
*auteur E. Gondet  -  April 2001
*
*argument
*    yes if 1 delete 10 files Vmm_* 
*
**/
wordint f77name(vmmdel)(int yes) {
   int vmmerr();
   int i, ier;

   if (yes == 1) 
   {
     for ( i = 0; i < NCLASSE; i++)
     {
       printf("Fichier i = %s\n", fclass_names[i]); 
       ier = unlink(fclass_names[i]);
       printf("\n unlink rend : %i\n", ier);
     }
     if (ier != 0)
       vmmerr("vmmdel",CONTROL_FILE_ERROR);

     ier = unlink(fcontrole_name);
       printf("\n unlink rend pour le fichier de controle : %i\n", ier);
     if (ier != 0)
       vmmerr("vmmdel",CONTROL_FILE_ERROR);

     fichiers_ouverts = 0;
   }

   return (0);
}

/***s/p vmmend
*
*objet(vmmend)
*     Liberer l'espace memoire du systeme de gestion de memoire virtuelle. 
*     Remettre a zero les tableaux de structures statiques controlant 
*     l'espace memoire.
*     vmmend doit etre la derniere fonction du systeme de gestion de
*     memoire virtuelle appelee par l'usager.
*
*     vmmend fait l'ouverture de tous les fichiers associes au systeme
*     de gestion de memoire virtuelle (VMM_01,... et VMM_controle)
*     et initialise les wp_Vmm (mot ou ecrire a la fin des fichiers)
*
*     vmmend detruit le fichier de controle si il existe.
*
*auteur E. Gondet  -  April 2001
*
*argument
*     None
*
**/
wordint
f77name(vmmend)()
{
   int vmmerr();
   void  lit_vmm_controle();
 /*   
    ouvre_ou_ferme_controle(); 
*/
 
   int i;
  
/*
 * s'assurer que vmmallc a ete appele avant vmmend
 */

   if(!called_vmmallc)
       return vmmerr("VMMEND",NO_CALL_TO_VMMALLC);

/* 
 * Reinitialiser les structures NAMES, SLICES, BLOCKS
 */
/*ETG OLD from vmmallc   for (i = 0; i < MAXSLICES; i++)
          SLICES[i].info.attributs = 0;
   for (i = 0; i < MAXBLOCKS; i++)
   {
          BLOCKS[i].info.attributs = 0;
          BLOCKS[i].prev_fb = BLOCKS[i].next_fb = -1;
   }                                                     ETG*/

   memset(BLOCKS, 0, MAXBLOCKS*sizeof(struct block_table) );
   memset(NAMES , 0, MAXNAMES*sizeof(struct name_table) );
   memset(SLICES, 0, MAXSLICES*sizeof(struct slice_table) );

/*
 * fermerture de tous les fichiers du systeme de gestion 
 */
   ouvre_ou_ferme_controle(0,0,"VMMEND");

/*
 * Destruction des fichiers Vmm_*
 */
   f77name(vmmdel)(1);

/*
 * Liberer la memoire allouee par VMM
 */

   free(BLOCKS[0].memadr); 
   fprintf(stdout,"Debug vmmend BLOCKS[0].memadr=%d\n",BLOCKS[0].memadr);

   if(BLOCKS[0].memadr != (wordint *) NULL)
        return(vmmerr("VMMALLC",SPACE_STILL_ALLOCATED));


   BLOCKS[0].slice_table_index = -1;
   BLOCKS[0].size = 0 ;

   if(debug_mode)
         fprintf(fdout," VMMEND-deallocation complete de l espace memoire de VMM\n");

 maxmem = 0, free_space = 0, nbslices = 0, nbvar = 0, nbblocks = 0;
 mot_de_passe = 0, pwd_set = 0;
 called_vmmallc = 0; /*called_vmmallc mis a  1 lors du 1er appel a vmmallc */
 
/* variables globales pour les unites logiques des fichiers Vmm_classe
   et Vmm_controle
 */
   for (i=0; i<NCLASSE; i++)
   {
     fclass[i] = 0;
     wp_Vmm[i] = 0; ; /*longueur des Vmm_0n*/
   }
   fcontrole=0;   /* nom du fichier: Vmm_controle */
   fichiers_ouverts = 0; 
   champs_bloques = 0; 
   reprise = 0;       /* mis a 1 si le fichier Vmm_controle est non vide */
   debug_mode = 0;
   checksum_mode = 0;
   espace_requis_max = 0; 
   champs_bloques_max = 0; 
   nbblocks_max = 0; 
   nb_appels_no_lock = 0;
   nb_appels_lock = 0;
   nb_ecritures = 0;
   nb_lectures = 0;
   first_free_bloc = 0;

   return(0);
   }

