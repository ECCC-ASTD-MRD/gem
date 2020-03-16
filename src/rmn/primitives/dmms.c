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

/************************************************************************
 *    D y n a m i c _ M e m o r y _ M a n a g e m e n t _ S y s t e m   *
 ************************************************************************/

#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <malloc.h>
#if defined (NEC)
#include <sys/time.h>
#include <sys/resource.h>
#endif
#include <rpnmacros.h>
#define ERROR 0
#define STACK 0
#define HEAP  1

struct blocmem {

  struct blocmem *fwd;
  struct blocmem *bwd;
  int *data[4];
  };

static struct blocmem stack_first, stack_last, heap_first, heap_last;
static struct blocmem *badptr;

static int init = 0, initmem = 0;

static wordint con;

static int ptrsize, *pointer, debug_mode=0, dejala=0, dmms_noabort=0;

#define single() {\
  if (dejala) {\
    fprintf(stderr," * * * ERROR * * *: more than one task in dmms\n");\
    f77name(tracebck)();\
    exit(50);\
    }\
  else dejala=1;\
}

#define sortie() {\
  dejala=0;\
}

#if defined (NEC)
/************************************************************************
 *                        c _ m e m u s e                               *
 ************************************************************************/
void c_memuse()
{
        INT_64 use, max;
	struct rlimit lims;    
	struct jresourcemem jusage;

	if (getresourcej( 0, CURR_UMEM, &jusage) == -1) {
		perror("getresource failed");
		exit(1);
	}
	use = jusage.jr_umem.mv_used / 1024 ; 

	if (getrlimitj(0,RLIMIT_UMEM, &lims) == -1) {
		perror("getrlimit failed");
		exit(1);
	}
	max = lims.rlim_max / 1024; 
	fprintf(stderr, "memory usage of jid (%d) is: %d kbytes\n", 
		getjid(0), use );
	fprintf(stderr, "max memory usage of jid (%d) is: %d kbytes\n", 
		getjid(0), max );
}
#endif


/************************************************************************
 *                      b l o c _ a l l o c                             *
 ************************************************************************/

/**
 *
 *auteur   M. Lepine - fev 92
 *
 *objet(bloc_alloc)
 *
 *     allocation d'un bloc de memoire (heap ou stack ) avec
 *     pointeurs avant et arriere et double pointeurs en boucle comme
 *     mecanisme de validation
 *
 *         ___________________   __________________>
 *        |  ________         | |  ________
 *        | |        |        | | |        |
 *     ----------------      ----------------
 *     |*|*|*| data |*|      |*|*|*| data |*|
 *     ----------------      ----------------
 *  <___|   |________|        |   |________|
 *                            |
 *      |_____________________| 
 *
 **/

struct blocmem *bloc_alloc(nbytes,mode)
int nbytes,mode;

{
  struct blocmem *ptbloc;
  unsigned int errptr;
  int lng, nitem, n;
  char *value, *getenv();

  single();
  ptrsize = sizeof(pointer);

  nitem = (nbytes + ptrsize - 1) / ptrsize;
  lng = sizeof(struct blocmem) + (nitem * ptrsize);

  ptbloc = (struct blocmem *) malloc(lng);
  if (ptbloc == NULL) {
     if (dmms_noabort)
       return((struct blocmem *) NULL);
     else {
       perror("bloc_alloc error can't allocate");
       fprintf(stderr,"bloc_alloc trying to allocate lng=%d bytes\n",lng);
#if defined(NEC)
       c_memuse();
#endif
       f77name(tracebck)();
       exit(7);
     }
  }

  if (! init) {
    /*    mallopt(M_DEBUG,1); */

    heap_first.bwd = (struct blocmem *) NULL;
    heap_first.fwd = &heap_last;
    heap_last.bwd = &heap_first;
    heap_last.fwd = (struct blocmem *) NULL;

    stack_first.bwd = (struct blocmem *) NULL;
    stack_first.fwd = &stack_last;
    stack_last.bwd = &stack_first;
    stack_last.fwd = (struct blocmem *) NULL;

    heap_first.data[0] = (int *) &(heap_first.data[1]);
    heap_first.data[1] = (int *) &(heap_first.data[0]);
    heap_last.data[0] = (int *) &(heap_last.data[1]);
    heap_last.data[1] = (int *) &(heap_last.data[0]);

    stack_first.data[0] = (int *) &(stack_first.data[1]);
    stack_first.data[1] = (int *) &(stack_first.data[0]);
    stack_last.data[0] = (int *) &(stack_last.data[1]);
    stack_last.data[1] = (int *) &(stack_last.data[0]);

    value = getenv("BAD_POINTER");
    if (value != NULL) {
      n = sscanf(value,"%x",&errptr);
      badptr = (struct blocmem *) errptr;
      fprintf(stderr,"Debug bad_pointer to look for is %#x\n",badptr);
      }
    else
      badptr = (struct blocmem *) 0;

    initmem = 0;
    value = getenv("INITMEM");
    if (value != NULL) {
      initmem = 1;
      if (strcmp(value,"ON") == 0)
         con = 0xFFFA5A5A;
      else
	 n = sscanf(value,"%x",&con);
      }

    value = getenv("DEBUG_MODE");
    debug_mode = ((value != NULL) && (strcmp(value,"OFF") != 0)
                   && (strcmp(value,"0") != 0));    
    init = 1;
    if (debug_mode) {
       fprintf(stdout,"DEBUG_MODE %s\n",value);
       fprintf(stdout,"Debug &heap_first =%#x\n",&heap_first);
       fprintf(stdout,"Debug &heap_last =%#x\n",&heap_last);
       fprintf(stdout,"Debug &stack_first =%#x\n",&stack_first);
       fprintf(stdout,"Debug &stack_last =%#x\n",&stack_last);
       }
    }

  if (badptr != (struct blocmem *) 0)
     if (badptr ==  ptbloc) {
        fprintf(stderr,"bloc_alloc bad_pointer %#x\n",ptbloc);
        f77name(tracebck)();
	exit(10);
	}

  if (mode == HEAP) {
    ptbloc->bwd = heap_last.bwd;
    ptbloc->fwd = &heap_last;
    heap_last.bwd = ptbloc;
    (ptbloc->bwd)->fwd = ptbloc;
    }
  else {
    ptbloc->bwd = stack_last.bwd;
    ptbloc->fwd = &stack_last;
    stack_last.bwd = ptbloc;
    (ptbloc->bwd)->fwd = ptbloc;
    }

  ptbloc->data[0] = (int *) &(ptbloc->data[nitem+1]);
  ptbloc->data[nitem+1] = (int *) &(ptbloc->data[0]);
  if (debug_mode) {
     fprintf(stdout,"\n");
     fprintf(stdout,"Debug alloc_bloc nitem = %d\n",nitem);
     fprintf(stdout,"Debug alloc_bloc lng = %d\n",lng);
     fprintf(stdout,"Debug alloc_bloc ptbloc =%#x\n",ptbloc);
     fprintf(stdout,"Debug alloc_bloc ptbloc->bwd =%#x\n",ptbloc->bwd);
     fprintf(stdout,"Debug alloc_bloc ptbloc->fwd =%#x\n",ptbloc->fwd);
     fprintf(stdout,"Debug alloc_bloc ptbloc->data[0] =%#x\n",ptbloc->data[0]);
     fprintf(stdout,"Debug alloc_bloc ptbloc->data[nitem+1] =%#x\n",ptbloc->data[nitem+1]);
     }

  if (initmem) {
     lng = (nitem-2) * ptrsize / sizeof(wordint);
     f77name(afix)(&(ptbloc->data[2]),&con,&lng);
/*     while (ptint < &(ptbloc->data[nitem+1])) {
       *ptint = con;
       ptint ++;
       }
*/
     }

  sortie();
  return(ptbloc);
  }

/************************************************************************
 *                      b l o c _ d e a l l o c                         *
 ************************************************************************/

/**
 *
 *auteur   M. Lepine - fev 92
 *
 *objet(bloc_dealloc)
 *
 *     desallocation d'un bloc de memoire (heap ou stack ) avec
 *     validation des pointeurs du bloc
 **/

int bloc_dealloc(ptbloc,mode)
struct blocmem *ptbloc;
int mode;
{
   int err;
   struct blocmem *pt;

   single();
   if (debug_mode) {
      fprintf(stdout,"\n");
      fprintf(stdout,"Debug bloc_dealloc ptbloc =%#x\n",ptbloc);
      }
   if (mode == HEAP) {
     err = bloc_check(ptbloc,0);
     if (err < 0) {
        f77name(tracebck)();
        exit(12);
	}
     (ptbloc->bwd)->fwd = ptbloc->fwd;
     (ptbloc->fwd)->bwd = ptbloc->bwd;
     free(ptbloc);
     sortie();
     return(0);
     }
   else {
     pt = ptbloc;
     stack_last.bwd = ptbloc->bwd;
     (ptbloc->bwd)->fwd = &stack_last;
     while (pt != &stack_last) {
       err = bloc_check(pt,0);
       if (err < 0) {
          f77name(tracebck)();
          exit(14);
	  }
       pt = ptbloc->fwd;
       free(ptbloc);
       ptbloc = pt;
       }
     }
   sortie();
   return(0);
   
   }

/************************************************************************
 *                      b l o c _ c h e c k                             *
 ************************************************************************/

/**
 *
 *auteur   M. Lepine - fev 92
 *
 *objet(bloc_check)
 *
 *     validation des differents pointeurs d'un bloc de memoire
 **/

int bloc_check(ptbloc,msg_level)
struct blocmem *ptbloc;
int msg_level;

{
   int **pt;

   if ((debug_mode) || (msg_level > 1)) {
      fprintf(stdout,"\n");
      fprintf(stdout,"Debug check ptbloc =%#x\n",ptbloc);
      fprintf(stdout,"Debug check ptbloc->bwd =%#x\n",ptbloc->bwd);
      fprintf(stdout,"Debug check ptbloc->fwd =%#x\n",ptbloc->fwd);
      }

   if (ptbloc->bwd == NULL) {
     fprintf(stderr,
     "block_check error: NULL backward pointer ptbloc=%#x\n",ptbloc);
     return(-1);
     }

   if (ptbloc->fwd == NULL) {
     fprintf(stderr,
     "block_check error: NULL forward pointer ptbloc=%#x\n",ptbloc);
     return(-2);
     }

   pt = (int **) ptbloc->data[0];

   if ((debug_mode) || (msg_level > 1)) {
      fprintf(stdout,"Debug check ptbloc->data[0] =%#x\n",ptbloc->data[0]);
      fprintf(stdout,"Debug check ptbloc->data[nitem+1] =%#x\n",*pt);
      }

   if (*pt != (int *) &(ptbloc->data[0])) {
     fprintf(stderr,
 "block_check error: internal pointers destroyed ptbloc=%#x\n",ptbloc);
     return(-3);
     }

   if (msg_level > ERROR)
       fprintf(stderr,"block_check OK \n");

   return(0);
   }

/************************************************************************
 *                        m e m _ c h e c k                             *
 ************************************************************************/

/**
 *
 *auteur   M. Lepine - fev 92
 *
 *objet(mem_check)
 *
 *     parcourir la liste des blocs de memoire (heap ou stack )
 *     pour validation des differents pointeurs
 **/

int mem_check(mode,msg_level)
int mode, msg_level;

{
   struct blocmem *ptbloc;
   int err;
   if (mode == HEAP)
      ptbloc = heap_first.fwd;
   else
      ptbloc = stack_first.fwd;
   while (ptbloc->fwd != (struct blocmem *) NULL) {
     if ((err = bloc_check(ptbloc,msg_level)) < 0)
        return(err);
     ptbloc = ptbloc->fwd;
     }

   return(0);
   }


/************************************************************************
 *                            m e m o i r c                             *
 ************************************************************************/

/**
 *
 *auteur   M. Lepine - mars 92
 *
 *objet(memoirc)
 *
 *     parcourir la liste des blocs de memoire (heap et stack )
 *     pour verification d'integrite des differents blocs de memoire
 **/

int f77name(memoirc)(msg_level)
wordint *msg_level;

{
   int errs, errh;

   if (! init)
      return(0);
   if (stack_first.fwd != &stack_last)
      fprintf(stderr,"memoirc warning: stack not empty \n");
   if (*msg_level > 1) {
       fprintf(stdout,"Debug &heap_first =%#x\n",&heap_first);
       fprintf(stdout,"Debug &heap_last =%#x\n",&heap_last);
       fprintf(stdout,"Debug &stack_first =%#x\n",&stack_first);
       fprintf(stdout,"Debug &stack_last =%#x\n",&stack_last);
       }
   errh = mem_check(HEAP,*msg_level);
   errs = mem_check(STACK,*msg_level);
   return((errh !=0) ? errh : (errs != 0) ? errs : 0);
   }





/************************************************************************
 *                            d m m s d b g                             *
 ************************************************************************/

/**
 *
 *auteur   M. Lepine - janvier 94
 *
 *objet(dmmsdbg)
 *
 *     Activer ou desactiver les messages de debug
 **/

void f77name(dmmsdbg)(dbgr)
wordint *dbgr;
{
   debug_mode = (*dbgr == 1) ? 1 : 0;
   }

/************************************************************************
 *                            d m m s n a b t                           *
 ************************************************************************/

/**
 *
 *auteur   M. Lepine - mars 2000
 *
 *objet(dmmsnabt)
 *
 *     Activer ou desactiver le mode abort pour manque de memoire
 **/

void f77name(dmmsnabt)(abort)
wordint *abort;
{
   dmms_noabort = (*abort == 1) ? 1 : 0;
   }


#if !defined (_FLOAT1)
/************************************************************************
 *                           h p a l l o c                              *
 ************************************************************************/

void
f77name(hpalloc)( addr, length, errcode, abort )
wordint *length, *errcode, *abort;
void **addr;
{
   struct blocmem *ptbloc;
   if (*length == 0) {
      fprintf(stderr,"HPALLOC error: 0 length\n");
      f77name(tracebck)();
      exit(13);
      }
   ptbloc = bloc_alloc(8 + *length * sizeof(wordint) * ((*abort==8) ? 2 : 1),HEAP);
   *addr =  (void *) &(ptbloc->data[2]);
   *errcode = (ptbloc == (struct blocmem *) NULL) ? 1 : 0;
   }

/************************************************************************
 *                          h p d e a l l c                             *
 ************************************************************************/

void 
f77name(hpdeallc)(addr, errcode, abort)
wordint *errcode, *abort;
char **addr;
{
   int offset=4*sizeof(addr);
   *errcode = bloc_dealloc((*addr)-offset,HEAP);
   }

/************************************************************************
 *                          c a _ a l l o c                             *
 *                                                                      *
 *   cached align allocation                                            *
 *                                                                      *
 ************************************************************************/

void
f77name(ca_alloc)(void **addr, wordint *length, wordint *errcode, wordint *abort, wordint *fpw2 )
{
   struct blocmem *ptbloc;
   int **pt_data1, **pt_aligned;
#if defined (AIX)
   int alignment[3] = {128, 128, 512};
#else
   int alignment[3] = {32, 32, 32};
#endif
   int nbytes, ajout, i;
   int pw2 = *fpw2;

   if (*length == 0) {
      fprintf(stderr,"CA_ALLOC error: 0 length\n");
      f77name(tracebck)();
      exit(13);
      }
   if (pw2 < 0) {
     pw2 = -pw2;
     if ((pw2 < 1) || (pw2 > 3)) {
       fprintf(stderr,"ca_alloc wrong value for alignment:%d\n",-pw2);
       exit(33);
     }
     nbytes = alignment[pw2-1];
   }
   else
     nbytes = 1 << pw2;
   ptbloc = bloc_alloc(nbytes + 8 + *length * sizeof(wordint) * ((*abort==8) ? 2 : 1),HEAP);
   pt_data1 = &(ptbloc->data[1]);
   pt_aligned = &(ptbloc->data[2]) + (nbytes / sizeof(pt_aligned));
   pt_aligned = (void *) (((PTR_AS_INT) pt_aligned) >> pw2);
   pt_aligned = (void *) (((PTR_AS_INT) pt_aligned) << pw2);
   ajout =  pt_aligned -  pt_data1;
   ptbloc->data[1] = (int *) ptbloc;
   for (i=0; i <= ajout; i++)
     ptbloc->data[1+i] = (int *) ptbloc;
   *addr =  (void *) pt_aligned;
   *errcode = (ptbloc == (struct blocmem *) NULL) ? 1 : 0;
   }

/************************************************************************
 *                          c a _ d e a l l c                           *
 ************************************************************************/

void 
f77name(ca_deallc)(addr, errcode, abort)
wordint *errcode, *abort;
int **addr;
{
   int **ptr;
   ptr = *addr - sizeof(addr);
   *errcode = bloc_dealloc(*ptr,HEAP);
   }
   
#if defined (OLD_32BIT_CODE)
/************************************************************************
 *                          m e m o i r h                               *
 ************************************************************************/
void
f77name(memoirh)(buf,ind,nw)
wordint buf[], *ind, *nw;
{
   int errcode, **ptr;
   struct blocmem *ptbloc;
   wordint *adr1;

#if !defined (OLD_32BIT_CODE)
	fprintf(stderr, "****************************************************\n");
	fprintf(stderr, "* ERROR: MEMOIRH                                   *\n");
	fprintf(stderr, "* This code is obsolete, will not work on a 64     *\n");
	fprintf(stderr, "* bit architecture and should not be used.         *\n");
	fprintf(stderr, "* allocate shoud be used instead                   *\n");
	fprintf(stderr, "* EXITING                                          *\n");
	fprintf(stderr, "****************************************************\n");
	exit(33);
#else
	
   if (*nw > 0) {
      ptbloc = bloc_alloc(8 + *nw * sizeof(wordint),HEAP);
      adr1 = (wordint *) &(ptbloc->data[2]);
      *ind = (adr1 - buf) + 1;
      }
   else {
      ptr = (int **) &buf[*ind - 1];
      errcode = bloc_dealloc(ptr-4,HEAP);
      }
#endif
   }

/************************************************************************
 *                           m e m o i r                                *
 ************************************************************************/
void
f77name(memoir)(buf,ind,nw)
wordint buf[], *ind, *nw;
{
   int  errcode, **ptr;
   struct blocmem *ptbloc;
   wordint *adr1;
#if !defined (OLD_32BIT_CODE)
	fprintf(stderr, "****************************************************\n");
	fprintf(stderr, "* ERROR: MEMOIR                                    *\n");
	fprintf(stderr, "* This code is obsolete, will not work on a 64     *\n");
	fprintf(stderr, "* bit architecture and should not be used.         *\n");
	fprintf(stderr, "* allocate shoud be used instead                   *\n");
	fprintf(stderr, "* EXITING                                          *\n");
	fprintf(stderr, "****************************************************\n");
	exit(44);
#else

   if (*nw > 0) {
      ptbloc = bloc_alloc(8 + *nw * sizeof(wordint),STACK);
      adr1 = (wordint *) &(ptbloc->data[2]);
      *ind = (adr1 - buf) + 1;
      }
   else {
      ptr = (int **) &buf[*ind - 1];
      errcode = bloc_dealloc(ptr-4,STACK);
      }
#endif
}
#endif

/************************************************************************
 *                          b k c h e c k                               *
 ************************************************************************/

void
f77name(bkcheck)(addr, errcode)
wordint **addr, *errcode;
{
   *errcode = bloc_check((*addr)-4,1);
   }

#endif

/************************************************************************
 *                          h p c h e c k                               *
 ************************************************************************/
void
f77name(hpcheck)(errcode)
wordint *errcode;
{
   if (*errcode == 0) 
      *errcode = mem_check(HEAP,0);
   else
      *errcode = mem_check(HEAP,1);
   }

/************************************************************************
 *                           m c h e c k                                *
 ************************************************************************/
void
f77name(mcheck)(errcode)
int *errcode;
{
   if (*errcode == 0) 
      *errcode = mem_check(STACK,0);
   else
      *errcode = mem_check(STACK,1);
   }



#if defined (_FLOAT1)
/************************************************************************
 *                          b k c h e c k                               *
 ************************************************************************/

void
f77name(bkcheck)(addr, errcode)
wordint *addr, *errcode;
{
   INT_64 temp;
   int **p;
   temp = (INT_64) *addr;
#if defined (ALL64)
   temp = temp << 3;
#else
   temp = temp << 2;
#endif
   p = (int **) temp;
   *errcode = bloc_check(p-3,1);
   }

/************************************************************************
 *                           h p a l l o c                              *
 ************************************************************************/
void
f77name(hpalloc)( addr, length, errcode, abort )
wordint *addr, *length, *errcode, *abort;
{
   struct blocmem *ptbloc;
   INT_64 laddr;
   if (*length == 0) {
      fprintf(stderr,"HPALLOC error: 0 length\n");
      f77name(tracebck)();
      exit(13);
      }
   ptbloc = bloc_alloc(*length * sizeof(wordint) * ((*abort==8) ? 2 : 1),HEAP);
   laddr = (INT_64) &(ptbloc->data[1]);
#if defined (ALL64)
   *addr = (wordint) (laddr >> 3);
#else
   *addr = (wordint) (laddr >> 2);
#endif
   *errcode = (*addr == NULL) ? 1 : 0;
   }

/************************************************************************
 *                          h p d e a l l c                             *
 ************************************************************************/
void 
f77name(hpdeallc)(addr, errcode, abort)
wordint *addr, *errcode, *abort;
{
   INT_64 temp;
   int **p;
   temp = (INT_64) *addr;
#if defined (ALL64)
   temp = temp << 3;
#else
   temp = temp << 2;
#endif
   p = (int **) temp;
   *errcode = bloc_dealloc(p-3,HEAP);
   }

/************************************************************************
 *                          m e m o i r h                               *
 ************************************************************************/
void
f77name(memoirh)(buf,ind,nw)
wordint buf[], *ind, *nw;
{
   int errcode, **ptr;
   INT_64 adr1, adr2;
   struct blocmem *ptbloc;

   if (*nw > 0) {
      ptbloc = bloc_alloc(*nw * sizeof(wordint),HEAP);
      adr1 = (INT_64) &(ptbloc->data[1]);
      adr2 = (INT_64) buf;
      *ind = (wordint) ((adr1 - adr2) / sizeof(wordint)) + 1;
      }
   else {
      ptr = (int **) &buf[*ind - 1];
      errcode = bloc_dealloc(ptr-3,HEAP);
      }
   }

/************************************************************************
 *                           m e m o i r                                *
 ************************************************************************/
void
f77name(memoir)(buf,ind,nw)
wordint buf[], *ind, *nw;
{
   int errcode, **ptr;
   INT_64 adr1, adr2;
   struct blocmem *ptbloc;

   if (*nw > 0) {
      ptbloc = bloc_alloc(*nw * sizeof(wordint),STACK);
      adr1 = (INT_64) &(ptbloc->data[1]);
      adr2 = (INT_64) buf;
      *ind = (wordint) ((adr1 - adr2) / sizeof(wordint)) + 1;
      }
   else {
      ptr = (int **) &buf[*ind - 1];
      errcode = bloc_dealloc(ptr-3,STACK);
      }
   }
#endif
