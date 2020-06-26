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

/*splitpoint a000000 */
/*****************************************************************************
 *                               X D F                                       *
 *                                                                           *
 *  EXtensible Directory File: a data bank system software                   *
 *                                                                           *
 *Authors                                                                    *
 *  M. Lepine - M. Valin                                                     *
 *  Revised version of original fortran code by J. Caveen and M. Lepine      *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include <stdio.h>
#ifndef WIN32    /*CHC/NRC*/
#include <unistd.h>
#endif
#include <stdlib.h>
#include <string.h>
#define XDF_OWNER
#include "qstdir.h"

static int endian_int=1;
static char *little_endian=(char *)&endian_int;
static int req_no=0;

static int init_package_done=0;
static key_descriptor truc[] = {
   { 'STI1', 7,7,33,0},
   { 'STI2',15,7,33,0},
   { 'STI3',23,7,33,0},
   { 'STI4',31,7,33,0},
   { 'STI5',39,7,33,0},
   { 'STI6',47,7,33,0},
   { 'STI7',55,7,33,0},
   { 'STI8',63,7,33,0},
   { 'STI9',71,7,33,0},
   { 'FLGS',95,23,0,0},
   { 'LATI',111,15,0,0},
   { 'LONG',127,15,0,0},
   { 'DATE',147,19,0,0},
   { 'DX  ',159,11,0,0},
   { 'IDTP',167,7,0,0},
   { 'DY  ',179,11,0,0},
   { 'HEUR',185,5,0,0},
   { 'MIN ',191,5,0,0}
};

static key_descriptor auxkey[] = {
   { 'NBLK',15,15,0,0},
   { 'OARS',31,15,0,0},
   { 'ELEV',44,12,0,0},
   { 'DRCV',55,10,0,0},
   { 'RUNN',63,7,0,0},
};

static key_descriptor fstcles[] = {
  { 'DEET', 23,23, 0,0},
  { 'NBIT', 31, 7, 0,0},
  { 'NI  ', 55,23, 0,0},
  { 'GTYP', 63, 7, 0,0},
  { 'NJ  ', 87,23, 0,0},
  { 'DTYP', 95, 7, 1,0},
  { 'NK  ',115,19, 0,0},
  { 'UBC ',127,11, 0,0},
  { 'NPAS',153,25, 0,0},
  { 'PAD7',159, 5, 0,0},
  { 'IG4 ',183,23, 0,0},
  { 'IG2A',191, 7, 0,0},
  { 'IG1 ',215,23, 0,0},
  { 'IG2B',223, 7, 0,0},
  { 'IG3 ',247,23, 0,0},
  { 'IG2C',255, 7, 0,0},
  { 'ET15',285,29,21,0},
  { 'PAD1',287, 1, 0,0},
  { 'ET6A',317,29,21,0},
  { 'PAD2',319, 1, 0,0},
  { 'ETBC',331,11,18,0},
  { 'TVAR',343,11,18,0},
  { 'PAD3',351, 7, 0,0},
  { 'NVAR',375,23,20,0},
  { 'PAD4',383, 7, 0,0},
  { 'IP1 ',411,27, 0,0},
  { 'LTYP',415, 3, 0,0},
  { 'IP2 ',443,27, 0,0},
  { 'PAD5',447, 3, 0,0},
  { 'IP3 ',475,27, 0,0},
  { 'PAD6',479, 3, 0,0},
  { 'DATE',511,31, 0,0}
};
/* prototypes declaration */

static int get_free_index();
static void init_file(int i);
static void init_package();
static INT_32 scan_random(int file_index);
static INT_32 add_dir_page(int file_index,int wflag);
static INT_32 rewind_file(int file_index, int handle);
static int create_new_xdf(int index, int iun, word_2 *pri, int npri,
                          word_2 *aux, int naux, char *appl);
static word next_match(int file_index);
static void build_gen_prim_keys(word *buf, word *keys, word *mask,
                                word *mskkeys, int index, int mode);
static void build_gen_info_keys(word *buf, word *keys, int index,
                                int mode);
int C_fst_match_req(int set_nb, int handle);
#include "proto.h"

/*splitpoint add_dir_page */
/*****************************************************************************
 *                   A D D _ D I R E C T O R Y _ P A G E                     *
 *                                                                           *
 *Object                                                                     *
 *  add a directory page to file file_index                                  *
 *                                                                           *
 *Arguments                                                                  *
 *  IN   file_index  file index in file table                                *
 *  IN   wflag       write flag. In case of a new file, or a file extenxion  *
 *                   the directory page has to be written on disk.           *
 *                                                                           *
 *****************************************************************************/

/* add a directory page to file file_index */
static INT_32 add_dir_page(int file_index,int wflag)
{
   register file_table_entry *f;
   int page_size, i, wdlng;
   full_dir_page *page;

   /* check if file exists */
   if( (f = file_table[file_index]) == NULL) return (ERR_NO_FILE);

   /* check if we can add a directory page, and if memory is available */
   if(f->npages >= MAX_DIR_PAGES) {
      sprintf(errmsg,"Too many records, no more directory pages available");
      return(error_msg("c_qdfrstr",ERR_DIR_FULL,ERROR));
   }
   page_size = sizeof(full_dir_page) +
                8*(f->primary_len)*ENTRIES_PER_PAGE;

   /* allocate directory page , link it into the directory page chain */
   if( (page=(page_ptr) calloc(1,page_size)) == NULL ) return (ERR_MEM_FULL);
   f->dir_page[f->npages] = page;
   (f->dir_page[f->npages])->next_page = NULL;
   if(f->npages == 0)        /* first page has no predecessor */
      (f->dir_page[f->npages])->prev_page = NULL;
   else{                /* update succesor to preceding page */
      (f->dir_page[f->npages])->prev_page = f->dir_page[f->npages-1];
      (f->dir_page[f->npages-1])->next_page = f->dir_page[f->npages];
   }

   /* initialize directory header and directory page entries */
   page->modified=0;
   page->true_file_index = file_index;
   page->dir.idtyp=0;
   wdlng=ENTRIES_PER_PAGE*f->primary_len +4;
   page->dir.lng=wdlng;
   page->dir.addr=WDTO64(f->nxtadr -1) +1;
   page->dir.reserved1=0;
   page->dir.reserved2=0;
   page->dir.nxt_addr=0;
   page->dir.nent=0;
   page->dir.chksum=0;
   for (i=0 ; i<=f->primary_len*ENTRIES_PER_PAGE ; i++ ) page->dir.entry[i]=0;
   if (wflag) {                        /* write directory page to file */
      if(f->npages != 0){   /* first page has no predecessor */
         (f->dir_page[f->npages-1])->dir.nxt_addr=page->dir.addr;
         }
      c_wawrit(f->iun,&page->dir,f->nxtadr,W64TOWD(wdlng));
      f->nxtadr += W64TOWD(wdlng);
      f->header->fsiz = WDTO64(f->nxtadr -1);
      f->header->nbd++;
      f->header->plst = page->dir.addr;
      page->modified=1;                /* checksum has to be computed */
      }
   f->npages++;
   return(0);
}

/*splitpoint address_from_handle */
/***************************************************************************** 
 *                 A D D R E S S _ F R O M _ H A N D L E                     * 
 *                                                                           *
 *Object                                                                     *
 *  calculates an address from an handle for a sequential file               *
 *                                                                           *
 *Arguments                                                                  *
 *  IN   handle      handle (cluster:2,address:22,file_index:8)              *
 *  IN   f           pointer to xdf file information structure               *
 *                                                                           *
 *****************************************************************************/

static INT_32 address_from_handle(int handle, file_table_entry *f)
{
  int addr;
  
  addr = (ADDRESS_FROM_HNDL(handle) << (2 * CLUSTER_FROM_HANDLE(handle)));
  if (f->fstd_vintage_89)
    addr = (addr * 15);
  addr = W64TOwd(addr) + 1;
  return(addr);
}

/*splitpoint build_fstd_info_keys */
/***************************************************************************** 
 *                  B U I L D _ F S T D _ I N F O _ K E Y S                  *
 *                                                                           *
 *Object                                                                     *
 *  Pack fstd info keys into buffer or get info keys from buffer             *
 *  depending on mode argument.                                              *
 *                                                                           *
 *Arguments                                                                  *
 *  IN/OUT  buf     buffer to contain the keys                               *
 *  IN/OUT  keys    info keys                                                *
 *  IN      index   file index in file table                                 *
 *  IN      mode    if mode = WMODE,write to buffer otherwise get keys from  *
 *                  buffer.                                                  *
 *                                                                           * 
 *****************************************************************************/

void build_fstd_info_keys(word *buf, word *keys, int index, int mode)
{
}

/*splitpoint build_fstd_prim_keys */
/***************************************************************************** 
 *                  B U I L D _ F S T D _ P R I M _ K E Y S                  *
 *                                                                           *
 *Object                                                                     *
 *  Pack fstd primary keys into buffer or get primary keys from buffer       *
 *  depending on mode argument.                                              *
 *                                                                           *
 *Arguments                                                                  *
 *  IN/OUT  buf     buffer to contain the keys                               *
 *  IN/OUT  keys    primary keys                                             *
 *  OUT     mask    search mask                                              *
 *  IN      mskkeys unpacked masks                                           *
 *  IN      index   file index in file table                                 *
 *  IN      mode    if mode = WMODE,write to buffer otherwise get keys from  *
 *                  buffer.                                                  *
 *                                                                           * 
 *****************************************************************************/

void build_fstd_prim_keys(word *buf, word *keys, word *mask, word *mskkeys,
                                int index, int mode)
{
   file_header *fh;
   int i, wi, sc, rmask, key, wfirst, wlast;

   buf += W64TOWD(1);         /* skip first 64 bit header */

   mask[0] = 0;              /* first 64 bits not part of the search mask */
   if (W64TOWD(1) > 1) mask[1] = 0;

   mask += W64TOWD(1);        /* skip first 64 bit header */

   fh = file_table[index]->header;

   if (mode == WMODE) {        /* write keys to buffer */

     for (i=0; i < W64TOWD(fh->lprm -1); i++) {
       buf[i] = keys[i];
       mask[i] = mskkeys[i];
     }
   }
   else
     for (i=0; i < W64TOWD(fh->lprm -1); i++)
       keys[i] = buf[i];
}

/*splitpoint build_gen_info_keys */
/***************************************************************************** 
 *                  B U I L D _ G E N _ I N F O _ K E Y S                    *
 *                                                                           *
 *Object                                                                     *
 *  Pack generic info keys into buffer or get info keys from buffer          *
 *  depending on mode argument.                                              *
 *                                                                           *
 *Arguments                                                                  *
 *  IN/OUT  buf     buffer to contain the keys                               *
 *  IN/OUT  keys    info keys                                                *
 *  IN      index   file index in file table                                 *
 *  IN      mode    if mode = WMODE,write to buffer otherwise get keys from  *
 *                  buffer.                                                  *
 *                                                                           * 
 *****************************************************************************/

static void build_gen_info_keys(word *buf, word *keys, int index,
                                int mode)
{
   file_header *fh;
   int i, wi, sc, rmask, key;
   int bitmot=32;


   fh = file_table[index]->header;

   if (mode == WMODE) {        /* write keys to buffer */

      for (i=0; i < fh->naux; i++)
         if (keys[i] != -1) {
            wi = fh->keys[i+fh->nprm].bit1 / bitmot;
            sc = (bitmot-1) - (fh->keys[i+fh->nprm].bit1 % bitmot);
            rmask = -1 << (fh->keys[i+fh->nprm].lcle);
            rmask = ~(rmask << 1);
 /* equivalent of << lcle+1 and covers 32 bit case */   
            key = keys[i];
            if ((fh->keys[i+fh->nprm].tcle /32) > 0)
               key = key & (~((key & 0x40404040) >> 1));
            buf[wi] = buf[wi] & (~(rmask << sc));       /* clear bits */
            buf[wi] = buf[wi] | ((key & rmask) << sc);
            }
      }
   else
      for (i=0; i < fh->naux; i++) {
         wi = fh->keys[i+fh->nprm].bit1 / bitmot;
         sc = (bitmot-1) - (fh->keys[i+fh->nprm].bit1 % bitmot);
         rmask = -1 << (fh->keys[i+fh->nprm].lcle);
         rmask = ~(rmask << 1);
         keys[i] = (buf[wi] >> sc) & rmask;
         }
}

/*splitpoint build_gen_prim_keys */
/***************************************************************************** 
 *                  B U I L D _ G E N _ P R I M _ K E Y S                    *
 *                                                                           *
 *Object                                                                     *
 *  Pack generic primary keys into buffer or get primary keys from buffer    *
 *  depending on mode argument.                                              *
 *                                                                           *
 *Arguments                                                                  *
 *  IN/OUT  buf     buffer to contain the keys                               *
 *  IN/OUT  keys    primary keys                                             *
 *  OUT     mask    search mask                                              *
 *  IN      mskkeys unpacked masks                                           *
 *  IN      index   file index in file table                                 *
 *  IN      mode    if mode = WMODE,write to buffer otherwise get keys from  *
 *                  buffer.                                                  *
 *                                                                           * 
 *****************************************************************************/

static void build_gen_prim_keys(word *buf, word *keys,
                                word *mask, word *mskkeys,
                                int index, int mode)
{
   file_header *fh;
   int i, wi, sc, rmask, key, wfirst, wlast;
   int bitmot=32;

   buf += 2;                  /* skip first 64 bit header */

   mask[0] = 0;              /* first 64 bits not part of the search mask */
   mask[1] = 0;

   mask += 2;                 /* skip first 64 bit header */

   fh = file_table[index]->header;

   if (mode == WMODE) {        /* write keys to buffer */

      wfirst = fh->keys[0].bit1 / bitmot;
      wlast = fh->keys[fh->nprm-1].bit1 / bitmot;
      for (i=wfirst; i <= wlast; i++) {
         mask[i] = 0;
         }

      for (i=0; i < fh->nprm; i++)
         if (keys[i] != -1) {
            wi = fh->keys[i].bit1 / bitmot;
            sc = (bitmot-1) - (fh->keys[i].bit1 % bitmot);
            rmask = -1 << (fh->keys[i].lcle);
            rmask = ~(rmask << 1);
 /* equivalent of << lcle+1 and covers 32 bit case */   
            key = keys[i];
            if ((fh->keys[i].tcle /32) > 0)
               key = key & (~((key & 0x40404040) >> 1));
            buf[wi] = buf[wi] & (~(rmask << sc));       /* clear bits */
            buf[wi] = buf[wi] | ((key & rmask) << sc);
            mask[wi] = mask[wi] & (~(rmask << sc));     /* clear bits */
            mask[wi] = mask[wi] | (rmask << sc);
            }
      }
   else
      for (i=0; i < fh->nprm; i++) {
         wi = fh->keys[i].bit1 / bitmot;
         sc = (bitmot-1) - (fh->keys[i].bit1 % bitmot);
         rmask = -1 << (fh->keys[i].lcle);
         rmask = ~(rmask << 1);
         keys[i] = (buf[wi] >> sc) & rmask;
         }
}

/*splitpoint c_qdfdiag */
/***************************************************************************** 
 *                             C _ Q D F D I A G                             *
 *                                                                           * 
 *Object                                                                     * 
 *   Establishes diagnosis and statistics for a xdf file                     *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    iun      unit number associated to the file                        * 
 *                                                                           * 
 *****************************************************************************/

int c_qdfdiag(int iun)
{
#define swap_4(mot) { register unsigned INT_32 tmp =(unsigned INT_32)mot; \
   mot = (tmp>>24) | (tmp<<24) | ((tmp>>8)&0xFF00) | ((tmp&0xFF00)<<8); }

   int index, index_fnom, ier, wasopen=0, addr, nw;
   int nrec_tot=0, nrec_act=0, nrec_eff=0, ndirect=0, leplusgros=0;
   int readpos, eofile=0, thesame;
   char appl[5], vers[5];
   file_header *fh;
   file_record header64;
   xdf_record_header header;

   index_fnom = fnom_index(iun);
   if (index_fnom == -1) {
      sprintf(errmsg,"file is not connected with fnom");
      return(error_msg("c_qdfrstr",ERR_NO_FNOM,WARNING));
      }

   if ((index = file_index(iun)) == ERR_NO_FILE) {
/*
 * open file and read file header
 */
     c_waopen(iun);
     c_waread(iun,&header64,1,W64TOWD(2));
     if (header64.data[0] != 'XDF0' && header64.data[0] !='xdf0') {
     /*if (strncmp(&header64.data[0], "XDF0", 4) != 0 && strncmp(&header64.data[0], "xdf0", 4) != 0) {*/
        sprintf(errmsg,"file is not XDF type\n");
        return(error_msg("c_qdfdiag",ERR_NOT_XDF,ERRFATAL));
        }
     if ((fh = calloc(1,header64.lng*8)) == NULL) {
        sprintf(errmsg,"memory is full\n");
        return(error_msg("c_qdfdiag",ERR_MEM_FULL,ERRFATAL));
        }
     c_waread(iun,fh,1,W64TOWD(header64.lng));
     }

   else {
/*
 * file is open, file_table[index]->header already contains the required info
 */
      wasopen = 1;
      fh = file_table[index]->header;
      }
   nw = c_wasize(iun);
   if(*little_endian) {  
     int  ct=fh->vrsn ;
     swap_4(ct);
     strncpy(vers,(char *)&(ct),4);
   }
   else strncpy(vers,(char *)&(fh->vrsn),4);
   if(*little_endian)
   {
     int  ct=fh->sign ;
     swap_4(ct);
     strncpy(appl,(char *)&(ct),4);
   }
   else strncpy(appl,(char *)&(fh->sign),4);

   vers[4] = '\0';
   appl[4] = '\0';
   readpos = 1 + W64TOwd(header64.lng);
   eofile = 0;
   while (! eofile) {
     c_waread(iun,&header,readpos,W64TOwd(1));
     addr = W64TOwd(header.addr-1) + 1;
/*     printf("Debug+ qdfdiag, readpos=%d, addr=%d \n",readpos,addr); */
     if (addr == readpos) {
       if (header.lng < W64TOwd(1)) {
         sprintf(errmsg,"Invalid record length=%d, addr=%d\n",header.lng,addr);
         return(error_msg("c_qdfdiag",ERR_BAD_LEN,ERRFATAL));
       }
       if (header.idtyp == 0)
         ndirect++;
       else {
         if (header.idtyp == 255)
           nrec_eff++;
         else
           nrec_act++;
         nrec_tot++;
         leplusgros = (header.lng > leplusgros) ? header.lng : leplusgros;
       }
       readpos += W64TOwd(header.lng);
       if (readpos > nw)
         eofile = 1;
     }
     else
       eofile = 1;
   } /* end while */
   thesame = ((fh->nxtn == nrec_tot) && (fh->nrec == nrec_act) &&
              ((fh->nxtn - fh->nrec) == nrec_eff) &&
              (fh->nbig == leplusgros) && (fh->nbd == ndirect));
   fprintf(stdout,"\nStatistics from file header for %s\n",
           FGFDT[index_fnom].file_name);
   fprintf(stdout,"\t file size (64 bit units)        %d\n",fh->fsiz); 
   fprintf(stdout,"\t number of rewrites              %d\n",fh->nrwr); 
   fprintf(stdout,"\t number of extensions            %d\n",fh->nxtn); 
   fprintf(stdout,"\t number of directory pages       %d\n",fh->nbd); 
   fprintf(stdout,"\t last directory page address     %d\n",fh->plst); 
   fprintf(stdout,"\t size of biggest record          %d\n",fh->nbig); 
   fprintf(stdout,"\t number erasures                 %d\n",fh->neff); 
   fprintf(stdout,"\t number of valid records         %d\n",fh->nrec); 
   fprintf(stdout,"\t XDF version                     %s\n",vers); 
   fprintf(stdout,"\t application signature           %s\n",appl);

   if (! thesame) {
      fprintf(stdout,"\n **** This file has been damaged ****\n"); 
      fprintf(stdout,"\nStatistics from file scan\n");
      fprintf(stdout,"\t number of extensions            %d\n",nrec_tot); 
      fprintf(stdout,"\t number of directory pages       %d\n",ndirect); 
      fprintf(stdout,"\t size of biggest record          %d\n",leplusgros); 
      fprintf(stdout,"\t number erasures                 %d\n",nrec_eff); 
      fprintf(stdout,"\t number of valid records         %d\n",nrec_act); 
   }
   else
      fprintf(stdout,"\n **** This file is OK ****\n");

   if (! wasopen) c_waclos(iun);
   free(fh);
   if (! thesame)
     return(ERR_DAMAGED);
   else 
     return(0);
}

/*splitpoint c_qdfmsig */
/***************************************************************************** 
 *                           C _ Q D F M S I G                               *
 *                                                                           * 
 *Object                                                                     * 
 *   Change the application signature of a file. This is used for example    *
 *   in conjuncture with burp files when the table used is not the           *
 *   official one.                                                           *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    iun      unit number associated to the file                        *
 *  IN    newappl  new application signature                                 *
 *                                                                           *
 *****************************************************************************/
int c_qdfmsig(int iun, char* newappl)
{
   file_header *fh;
   int index, index_fnom;

   index_fnom = fnom_index(iun);
   if (index_fnom == -1) {
      sprintf(errmsg,"file is not connected with fnom");
      return(error_msg("c_qdfmsig",ERR_NO_FNOM,ERROR));
      }

   if ((index = file_index(iun)) == ERR_NO_FILE) {
     sprintf(errmsg,"file is not open");
     return(error_msg("c_qdfmsig",ERR_NO_FILE,ERROR));
     }

   fh = file_table[index]->header;
   /*   sprintf(fh->sign,"%-4.4s",newappl); */
   fh->sign = newappl[0] << 24 | newappl[1] << 16 | newappl[2] << 8 | newappl[3];
   return(0);
}

/*splitpoint c_qdfput */
/***************************************************************************** 
 *                             C _ Q D F P U T                               * 
 *                                                                           * 
 *Object                                                                     * 
 *   Add to the buffer element elem of length nbits ending a bit dernit.     *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN/OUT  buf      buffer to contain the element to add                    * 
 *  IN      elem     element to add                                          *
 *  IN      derbit   position of last bit in buf                             *
 *  IN      nbits    length in bit of element to add                         *
 *                                                                           *
 *****************************************************************************/

int c_qdfput(word *buf,int elem,int derbit,int nbits)
{
   int wi,sc,msk,bitmot=32;

   wi = derbit / bitmot;
   sc = (bitmot-1) - (derbit % bitmot);
   msk = ~(-1 << nbits);
   buf[wi] = buf[wi] & (~(msk << sc));
   buf[wi] = buf[wi] | ((elem & msk) << sc);
   return(0);
}

/*splitpoint c_qdfrstr */
/***************************************************************************** 
 *                             C _ Q D F R S T R                             *
 *                                                                           * 
 *Object                                                                     * 
 *   Restores a file to its original length (the one before a task that      *
 *   has prematurely terminated while doing appends)                         *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    inp      input unit number associated to the file                  * 
 *  IN    outp     output unit number associated to the file                 * 
 *                                                                           * 
 *****************************************************************************/

int c_qdfrstr(int inp, int outp)
{
   int index_fnom, ier, nw, rwpos, lng, i, alire;
#define Buflen 8192
   int buffer[Buflen];
   file_header fh;
   file_record header64;

   index_fnom = fnom_index(inp);
   if (index_fnom == -1) {
      sprintf(errmsg,"file (unit=%d) is not connected with fnom",inp);
      return(error_msg("c_qdfrstr",ERR_NO_FNOM,WARNING));
      }

   index_fnom = fnom_index(outp);
   if (index_fnom == -1) {
      sprintf(errmsg,"file (unit=%d) is not connected with fnom",outp);
      return(error_msg("c_qdfrstr",ERR_NO_FNOM,WARNING));
      }

/*
 * open file and read file header
 */
   c_waopen(inp);
   c_waread(inp,&header64,1,W64TOWD(2));
   if (header64.data[0] != 'XDF0' && header64.data[0] !='xdf0') {
   /*if (strncmp(&header64.data[0], "XDF0", 4) != 0 && strncmp(&header64.data[0], "xdf0", 4) != 0) {*/
     sprintf(errmsg,"file is not XDF type\n");
     return(error_msg("c_qdfrstr",ERR_NOT_XDF,ERRFATAL));
   }
   lng = W64TOWD(header64.lng);
   nw = c_wasize(inp);

   if (lng > nw) {
     sprintf(errmsg,"Invalid header file length=%d\n",header64.lng);
     return(error_msg("c_qdfrstr",ERR_BAD_LEN,ERRFATAL));
   }
   c_waread(inp,&fh,1,W64TOWD(header64.lng));
   c_waopen(outp);
/*
 * reset read/write flag to zero
 */
   fh.rwflg = 0;
   c_wawrit(outp,&fh,1,W64TOwd(header64.lng));
   rwpos = 1 + W64TOwd(header64.lng);

   lng = W64TOwd(fh.fsiz - header64.lng);
   alire = (lng < Buflen) ? lng : Buflen;
   while (alire > 0) {
     for (i=0; i < Buflen; i++)
       buffer[i] = 0;
     c_waread(inp,&buffer,rwpos,alire);
     c_wawrit(outp,&buffer,rwpos,alire);
     lng -= alire;
     rwpos += alire;
     alire = (lng < Buflen) ? lng : Buflen;
   }
   
   c_waclos(inp);
   c_waclos(outp);

   return(0);
}

/*splitpoint c_xdfadd */
/***************************************************************************** 
 *                             C _ X D F A D D                               * 
 *                                                                           * 
 *Object                                                                     * 
 *   Add to the end of the record contain into buf, nelm*nbits bits from     *
 *   donnees. The number of bits must be a multiple of 64.                   *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN/OUT  buf      buffer to contain the record                            * 
 *  IN      donnees  data bits to add to buffer                              *
 *  IN      nelm     number of elements                                      *
 *  IN      nbits    number of bits per element                              *
 *  IN      datyp    data type of elements to add                            *
 *                                                                           *
 *****************************************************************************/

int c_xdfadd(word *buffer, word *donnees, int nelm, int nbits, int datyp)
{
   int index_word, nbwords, mode, i;
   buffer_interface_ptr buf = (buffer_interface_ptr) buffer;
   int ier;
   

   if (((datyp == 3) || (datyp == 5)) && (nbits != 8)) {
      sprintf(errmsg,"nbits must be 8 for datyp %d",datyp);
      return(error_msg("c_xdfadd",ERR_BAD_DATYP,ERRFATAL));
      }

   nbwords = (nelm * nbits + 63) / 64;
   nbwords = W64TOWD(nbwords);

   index_word = buf->nbits / (sizeof(word) * 8);

   if ((index_word + nbwords - 1) > buf->nwords) {
      sprintf(errmsg,"buffer not big enough for insertion");
      return(error_msg("c_xdfadd",ERR_BAD_DIM,ERROR));
      }

   switch (datyp)
      {
      case 0:         /* transparent mode */
      case 3:

         for (i=0; i < nbwords; i++)
            buf->data[index_word+i] = donnees[i];
         break;

#if !defined(NEC)
      case 6:
      case 8:

         for (i=0; i < nbwords; i++)
            buf->data[index_word+i] = donnees[i];
         break;

      case 7:
      case 9:

         if (*little_endian) 
           for (i=0; i < nbwords; i+=2) {
             buf->data[index_word+i] = donnees[i+1];
             buf->data[index_word+i+1] = donnees[i];
           }
         else
           for (i=0; i < nbwords; i++)
             buf->data[index_word+i] = donnees[i];
         break;
#endif
      
      case 5:        /* upper char only */
         for (i=0; i < nbwords; i++)
            buf->data[index_word+i] = upper_case_word(donnees[i]);
         break;
      
      case 2: mode = 1;
         ier = compact_integer(donnees,(void *) NULL,&(buf->data[index_word]),nelm,
                               nbits,0,xdf_stride,mode);
         break;

      case 4: mode = 3;
         ier = compact_integer(donnees,(void *) NULL,&(buf->data[index_word]),nelm,
                               nbits,0,xdf_stride,mode);
         break;

      default: sprintf(errmsg,"invalid datyp=%d\n",datyp);
         return(error_msg("c_xdfadd",ERR_BAD_DATYP,ERROR));

      } /* end switch */

   buf->nbits += nbwords * sizeof(word) * 8;
   return(0);
}

/*splitpoint c_xdfcle */
/***************************************************************************** 
 *                             C _ X D F C L E                               *
 *                                                                           * 
 *Object                                                                     * 
 *   Pack key descriptors into 2 different 32 bit wide elements.             *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    keyname  name of the key (max 4 char)                              * 
 *  IN    bit1     last right bit of key in record                           *
 *  IN    lkey     key length in bits                                        *
 *  IN    tkey     type of key                                               *
 *                    0 : unsigned int                                       *
 *                    1 : ascii char                                         *
 *                    2 : ascii char                                         *
 *                    3 : ascii char                                         *
 *                    4 : ascii char                                         *
 * OUT    desc1    first element to contain the keyname                      *
 * OUT    desc2    second element to contain the bit position length and type*
 *                                                                           *
 *****************************************************************************/
int c_xdfcle(char *keyname,int bit1,int lkey,int tkey,int *desc1,int *desc2)
{
   int i, bitpos, rmask, shift_count;
   int bitmot=32;
   
   *desc1=0; *desc2=0;
   
   for(i=0; (i<4 && *keyname); i++,*keyname++)
     *desc1 = (*desc1 <<8) | (*keyname & 0xff);

   while (i++ < 4)
     *desc1 = (*desc1 <<8) | ' ';

   bitpos = BPBIT1 - (64-bitmot);
   shift_count = (bitmot-1) - bitpos & (bitmot-1);
   rmask = -1 >> (bitmot -LBIT1);
   *desc2 = *desc2 | ((bit1 & rmask) << shift_count);

   bitpos = BPLCLE - (64-bitmot);
   shift_count = (bitmot-1) - bitpos & (bitmot-1);
   rmask = -1 >> (bitmot -LLCLE);
   *desc2 = *desc2 | (((lkey -1) & rmask) << shift_count);

   bitpos = BPTCLE - (64-bitmot);
   shift_count = (bitmot-1) - bitpos & (bitmot-1);
   rmask = -1 >> (bitmot -LTCLE);
   *desc2 = *desc2 | ((tkey & rmask) << shift_count);

   return(0);
}

/*splitpoint c_xdfcls */
/***************************************************************************** 
 *                            C _ X D F C L S                                *
 *                                                                           * 
 *Object                                                                     * 
 *   Closes the XDF file. Rewrites file header, computes directory checksum  *
 *   and rewrites directory pages if modified.                               *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *                                                                           * 
 *****************************************************************************/

int c_xdfcls(int iun)
{
   int index,index_fnom,i,j,lng64,width,open_mode;
   file_table_entry *f;
   xdf_dir_page * curpage;
   word32 * check32, checksum;
   word *entry;
   xdf_record_header *rec;
 
   index_fnom = fnom_index(iun);
   if (index_fnom == -1) {
      sprintf(errmsg,"file is not connected with fnom");
      return(error_msg("c_xdfcls",ERR_NO_FNOM,WARNING));
      }

   if ((index = file_index(iun)) == ERR_NO_FILE) {
     sprintf(errmsg,"file is not open");
     return(error_msg("c_xdfcls",ERR_NO_FILE,WARNING));
     }

   f = file_table[index];

   if ((f->header->rwflg != RDMODE) && (!FGFDT[index_fnom].attr.read_only))
      /* rewrite file header */
     c_wawrit(iun,f->header,1,W64TOWD(f->header->lng)); 

   if (f->modified) {  /* file has been modified rewrite dir. pages */
      for (i=0; i < f->header->nbd; i++) {
         if (f->dir_page[i]->modified) {
/*
 * reset idtyp entries to original value if modified by scan_dir_page
 */
            width = W64TOWD(f->primary_len);
            entry = (f->dir_page[i])->dir.entry;
            for (j=0; j < (f->dir_page[i])->dir.nent; j++) {
               rec = (xdf_record_header *) entry;
               if ((rec->idtyp | 0x80) == 254) {
                  rec->idtyp = 255;
                  c_wawrit(iun,rec,W64TOWD(rec->addr-1)+1,W64TOWD(1));
                  }
               rec->idtyp = 
                    ((rec->idtyp | 0x80) == 255) ? 255 : (rec->idtyp & 0x7f);
               entry += width;
               }
/*
 * compute checksum and rewrite page
 */
            lng64 = f->primary_len * ENTRIES_PER_PAGE + 4;
            curpage = &((f->dir_page[i])->dir);
            checksum = curpage->chksum;
            check32 = (word32 *) curpage;
            for (j=4; j < W64TOWD(lng64); j++)
               checksum ^= check32[j];
            curpage->chksum = checksum;
            c_wawrit(iun,curpage,W64TOWD(curpage->addr-1)+1,W64TOWD(lng64));
            f->dir_page[i]->modified = 0;
            } /* end if page modified */
         } /* end for i */
      if (f->xdf_seq) {
        { int trunc_to;
        trunc_to = FGFDT[index_fnom].file_size * sizeof(word);
        c_secateur(FGFDT[index_fnom].file_name,trunc_to);
        }
      }
      f->modified = 0;
      } /* end if file modified */

   if (!xdf_checkpoint) {
     if ((f->header->rwflg != RDMODE) && (!FGFDT[index_fnom].attr.read_only)) {
      /* rewrite file header */
       f->header->rwflg = 0;
       c_wawrit(iun,f->header,1,W64TOWD(f->header->lng)); 
     }

     c_waclos(iun);
/*
 * free allocated pages
 */
     for (i=0; i < f->npages; i++)
       free(f->dir_page[i]);

/* reset file informations */
     init_file(index);
   }
   else
     xdf_checkpoint = 0;
   return(0);
}

/*splitpoint c_xdfcut */
/***************************************************************************** 
 *                             C _ X D F C U T                               *
 *                                                                           * 
 *Object                                                                     * 
 *   Retrieve nelm elements from buf starting at bit position bitpos.        *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    nelm     number of elements to cut from buf                        * 
 *  IN    bitpos   starting bit position for cutting                         *
 *  IN    nbit     number of bits kept per element                           *
 *  IN    datyp    data type                                                 *
 * IN/OUT buf      buffer to contain the modified record                     *
 *                                                                           *
 *****************************************************************************/

int c_xdfcut(void *buffer, int bitpos, int nelm, int nbits, int datyp)
{
   int nbwords, index_word, last_ind, i;
   buffer_interface_ptr buf = (buffer_interface_ptr) buffer;

   if ((bitpos % 64) != 0) {
      sprintf(errmsg,"bitpos must be a multiple of 64");
      return(error_msg("c_xdfcut",ERR_BAD_ADDR,ERRFATAL));
      }

   if ((datyp == 3) || (datyp == 5) && (nbits != 8)) {
      sprintf(errmsg,"nbits must be 8 for datyp %d",datyp);
      return(error_msg("c_xdfcut",ERR_BAD_DATYP,ERRFATAL));
      }

   nbwords = (nelm * nbits + 63) / 64;
   nbwords = W64TOWD(nbwords);

   index_word = buf->data_index + (bitpos / (sizeof(word) * 8));

   last_ind = buf->record_index + (buf->nbits / (sizeof(word) *8));

   /* move buffer content nbwords to the left */
   if (last_ind != index_word)
      for (i=index_word; i <= last_ind; i++)
         buf->data[i] = buf->data[i+nbwords];

   buf->nbits -= nbwords * sizeof(word) * 8;
   return(0);
}

/*splitpoint c_xdfdel */
/***************************************************************************** 
 *                            C _ X D F D E L                                * 
 *                                                                           * 
 *Object                                                                     * 
 *   Delete record referenced by handle.                                     *
 *   Deleted record are marked as idtyp = X111111X (X = don't care bits)     *
 *   and will be marked as idtyp=255 upon closing of the file.
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    handle  file index, page number and record number to record        * 
 *                                                                           * 
 *****************************************************************************/

int c_xdfdel(int handle)
{
   int index, page_number, record_number, idtyp, i, addr;
   file_table_entry *f;
   file_record *record;
   word *rec;
   page_ptr target_page;
   xdf_record_header header;

   index =         INDEX_FROM_HANDLE(handle);
   page_number =   PAGENO_FROM_HANDLE(handle);
   record_number = RECORD_FROM_HANDLE(handle);

   /* validate index, page number and record number */

   if ((index >= MAX_XDF_FILES) || (file_table[index] == NULL) ||
       (file_table[index]->iun < 0)) {
      sprintf(errmsg,"invalid handle, invalid file index\n");
      return(error_msg("c_xdfdel",ERR_BAD_HNDL,ERROR));
      }

   f = file_table[index];

   if ((f->header->rwflg == RDMODE) || (f->header->rwflg == APPEND)) {
      sprintf(errmsg,"file is open in read or append mode only\n");
      return(error_msg("c_xdfdel",ERR_RDONLY,ERROR));
      }

   if (f->cur_info->attr.read_only) {
      sprintf(errmsg,"file is read only\n");
      return(error_msg("c_xdfdel",ERR_RDONLY,ERROR));
   }

   if (! f->xdf_seq) {
     if (page_number < f->npages) {   /* page is in current file */
       target_page = f->dir_page[page_number];
     }
     else {                           /* page is in a link file */
       if (f->link == -1) {
         sprintf(errmsg,"page number=%d > last page=%d and file not linked\n",
                 page_number,f->npages-1);
           return(error_msg("c_xdfdel",ERR_BAD_PAGENO,ERROR));
       }
       target_page = f->dir_page[f->npages-1];
       for (i=0; (i<= (page_number - f->npages)) && target_page; i++)
         target_page = (target_page)->next_page;
       if (target_page == NULL) {
         sprintf(errmsg,"invalid handle, invalid page number\n");
         return(error_msg("c_xdfdel",ERR_BAD_PAGENO,ERROR));
       }
     }

     if (record_number > target_page->dir.nent) {
       sprintf(errmsg,"invalid handle, invalid record number\n");
       return(error_msg("c_xdfdel",ERR_BAD_HNDL,ERROR));
     }

     rec = target_page->dir.entry + record_number * W64TOWD(f->primary_len);
     record = (file_record *) rec;
   
     idtyp = record->idtyp;
   }
   else { /* xdf sequential */
     addr = address_from_handle(handle,f);
     c_waread(f->iun,&header,addr,W64TOWD(1));
     idtyp = header.idtyp;
   }

   if (idtyp == 0) {
      sprintf(errmsg,"special record idtyp=0\n");
      return(error_msg("c_xdfdel",ERR_SPECIAL,ERROR));
      }

   if ((idtyp & 0x7E) == 0x7E) {
      sprintf(errmsg,"record already deleted\n");
      return(error_msg("c_xdfdel",ERR_DELETED,WARNING));
      }
   
   if (! f->xdf_seq) {
     /* update directory entry */
     record->idtyp = 0xFE;               /* 254 */
     target_page->modified =1;
   }
   else { /* xdf sequential */
     header.idtyp = 255;      /* deleted */
     c_wawrit(f->iun,&header,addr,W64TOWD(1));
   }

   /* update file header */
   f->header->neff++;
   f->header->nrec--;

   f->modified = 1;
   return(0);
}

/*splitpoint c_xdfget */
/***************************************************************************** 
 *                            C _ X D F G E T                                * 
 *                                                                           * 
 *Object                                                                     * 
 *   Obtain record referenced by handle in buf                               * 
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    handle  file index, page number and record number to record        * 
 *  OUT   buf     buffer to contain record                                   * 
 *                                                                           * 
 *****************************************************************************/

int c_xdfget(int handle, buffer_interface_ptr buf)
{
   int *aux_keys = NULL;

   return(c_xdfget2(handle,buf,aux_keys));
}

/*splitpoint c_xdfget2 */
/***************************************************************************** 
 *                            C _ X D F G E T 2                              * 
 *                                                                           * 
 *Object                                                                     * 
 *   Obtain record referenced by handle in buf                               * 
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    handle  file index, page number and record number to record        * 
 *  OUT   buf     buffer to contain record                                   * 
 *                                                                           * 
 *****************************************************************************/

int c_xdfget2(int handle, buffer_interface_ptr buf, int *aux_ptr)
{
   int index, record_number, page_number, i, idtyp, addr, lng, lngw;
   int offset, nw, nread;
   file_table_entry *f;
   file_record *record;
   word *rec;
   max_dir_keys argument_not_used;
   page_ptr target_page;

   index =         INDEX_FROM_HANDLE(handle);
   page_number =   PAGENO_FROM_HANDLE(handle);
   record_number = RECORD_FROM_HANDLE(handle);

   /* validate index, page number and record number */

   if ((index >= MAX_XDF_FILES) || (file_table[index] == NULL) ||
       (file_table[index]->iun < 0)) {
      sprintf(errmsg,"invalid handle, invalid file index\n");
      return(error_msg("c_xdfget",ERR_BAD_HNDL,ERROR));
      }

   f = file_table[index];

   if (! f->xdf_seq) {
     if (page_number < f->npages) {   /* page is in current file */
       target_page = f->dir_page[page_number];
     }
     else {                           /* page is in a link file */
       if (f->link == -1) {
         sprintf(errmsg,"page number=%d > last page=%d and file not linked\n",
                 page_number,f->npages-1);
           return(error_msg("c_xdfget",ERR_BAD_PAGENO,ERROR));
       }
       target_page = f->dir_page[f->npages-1];
       for (i=0; (i<= (page_number - f->npages)) && target_page; i++)
         target_page = (target_page)->next_page;
       if (target_page == NULL) {
         sprintf(errmsg,"invalid handle, invalid page number\n");
         return(error_msg("c_xdfget",ERR_BAD_PAGENO,ERROR));
       }
       f = file_table[target_page->true_file_index];
     }
     
     if (record_number > target_page->dir.nent) {
       sprintf(errmsg,"invalid handle, invalid record number\n");
       return(error_msg("c_xdfget",ERR_BAD_HNDL,ERROR));
     }

     rec = target_page->dir.entry + record_number * W64TOWD(f->primary_len);
     record = (file_record *) rec;
   }
   else {
     if (! f->valid_pos) {
       sprintf(errmsg,"no valid file position for sequential file\n");
       return(error_msg("c_xdfget",ERR_NO_POS,ERROR));
     }
     record = (file_record *) f->head_keys;
     if (address_from_handle(handle,f) != W64TOWD(record->addr-1)+1) {
       sprintf(errmsg,"invalid handle, invalid address=%d record address=%d\n",
               address_from_handle(handle,f),W64TOWD(record->addr-1)+1);
       return(error_msg("c_xdfget",ERR_BAD_HNDL,ERROR));
     }
   }
   
   idtyp = record->idtyp;
   addr = record->addr;
   lng = record->lng;
   lngw = W64TOWD(lng);

   if (idtyp == 0) {
      sprintf(errmsg,"special record idtyp=0\n");
      return(error_msg("c_xdfget",ERR_SPECIAL,ERROR));
      }

   if ((idtyp & 0x7E) == 0x7E) {
      sprintf(errmsg,"deleted record\n");
      return(error_msg("c_xdfget",ERR_DELETED,ERROR));
      }

   nw = buf->nwords;
   offset = 0;
   if (nw < 0) {
     if (buf->nbits != -1) {
       sprintf(errmsg,"dimension of buf is invalid = %d\n",nw);
       return(error_msg("c_xdfget",ERR_BAD_DIM,ERROR));
     }
     nw = -nw;          /* data only, no directory entry in record */
     if (! f->fstd_vintage_89) 
       offset = W64TOWD(f->primary_len + f->info_len); 
     else
       if (f->xdf_seq)  /* old standard sequential */
         offset = 30;
   }
   if (lngw > (nw - RECADDR +1)) {
       sprintf(errmsg,"dimension of buf (%d) < record size (%d)\n",
              nw,lngw);
      return(error_msg("c_xdfget",ERR_BAD_DIM,ERROR));
      }

   buf->nbits = lngw * 8 * sizeof(word);
   buf->record_index = RECADDR;
   buf->data_index = buf->record_index + W64TOWD(f->primary_len + f->info_len);
   buf->iun = f->iun;
   buf->aux_index = buf->record_index + W64TOWD(f->primary_len);

   if (aux_ptr != NULL) {
     *aux_ptr = 0;
     *(aux_ptr+1) = 0;
   }
   if ( (aux_ptr != NULL) && (!f->fstd_vintage_89) && (!f->xdf_seq) ) {
     c_waread(buf->iun,aux_ptr,W64TOWD(addr-1)+1+W64TOWD(f->primary_len),W64TOWD(f->info_len));
   }

   for(i=0; i < lngw; i++)
      buf->data[i] = 0;

   nread = c_waread2(buf->iun,&(buf->data),W64TOWD(addr-1)+1+offset,lngw-offset);
   if (nread != lngw-offset) {
       sprintf(errmsg,"short read, truncated record, asking for %d, got %d\n",lngw-offset,nread);
       return(error_msg("c_xdfget",ERR_SHORT_READ,ERROR));
   }
   else
     return(0);

}

/*splitpoint c_xdfgop */
/***************************************************************************** 
 *                            C _ X D F G O P                                * 
 *                                                                           * 
 *Object                                                                     * 
 *   Get different options settings values.                                  *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    optname  name of option to get                                     * 
 *  OUT   optc     value of option if type is character                      *
 *  OUT   optv     value of option if type is integer                        *
 *                                                                           * 
 *****************************************************************************/

int c_xdfgop(char *optname, char *optc, int *optv)
{
   if (strstr(optname,"ERRTOLR") || strstr(optname,"errtolr")) {

      if (xdf_toler == TRIVIAL)
         strcpy(optc,"TRIVIAL");
      else if (xdf_toler == INFORM)
         strcpy(optc,"INFORM");
      else if (xdf_toler == WARNING)
         strcpy(optc,"WARNING");
      else if (xdf_toler == ERROR)
         strcpy(optc,"ERROR");
      else if (xdf_toler == ERRFATAL)
         strcpy(optc,"ERRFATAL");
      else if (xdf_toler == SYSTEM)
         strcpy(optc,"SYSTEM");
      }
   else if (strstr(optname,"MSGLVL") || strstr(optname,"msglvl")) {

      if (msg_level == TRIVIAL)
         strcpy(optc,"TRIVIAL");
      else if (msg_level == INFORM)
         strcpy(optc,"INFORM");
      else if (msg_level == WARNING)
         strcpy(optc,"WARNING");
      else if (msg_level == ERROR)
         strcpy(optc,"ERROR");
      else if (msg_level == ERRFATAL)
         strcpy(optc,"ERRFATAL");
      else if (msg_level == SYSTEM)
         strcpy(optc,"SYSTEM");
      }
   else {
      sprintf(errmsg,"invalid option name: %s",optname);
      return(error_msg("c_xdfgop",ERR_BAD_OPT,ERROR));
      }
   return(0);
}

/*splitpoint c_xdfhdr */
/***************************************************************************** 
 *                            C _ X D F H D R                                *
 *                                                                           * 
 *Object                                                                     * 
 *   Get the descriptive parameters of the record contain into buf           *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    buf      buffer that contains the record                           * 
 *  OUT   addr     record starting address                                   *
 *  OUT   lng      length of the record in 64 bit units                      *
 *  OUT   idtyp    record id type                                            *
 *  OUT   primk    primary keys                                              *
 *  IN    nprim    number of primary keys                                    *
 *  OUT   info     info keys                                                 *
 *  IN    ninfo    number of info keys                                       *
 *                                                                           * 
 *****************************************************************************/

int c_xdfhdr(buffer_interface_ptr buf ,int *addr,int *lng,int *idtyp,
             word *primk,int nprim,word *info,int ninfo)
{
   int index, record_number, page_number, i;
   file_table_entry *f;
   file_record *record;
   max_dir_keys argument_not_used;
   word *mskkeys = NULL;

   record = (file_record *) buf->data;
  
   *idtyp = record->idtyp;
   *addr = record->addr;
   *lng = record->lng;
   
   if ((index = file_index(buf->iun)) == ERR_NO_FILE) {
      sprintf(errmsg,"file is not open");
      return(error_msg("c_xdfhdr",ERR_NO_FILE,WARNING));
      }

   f = file_table[index];

   f->build_primary(buf->data,primk,argument_not_used,mskkeys,
                    index,RDMODE);
   if (ninfo > 0) f->build_info(buf->data + W64TOWD(f->primary_len),
                                info,index,RDMODE);

   return(0);

}

/*splitpoint c_xdfimp */
/***************************************************************************** 
 *                            C _ X D F I M P                                *
 *                                                                           * 
 *Object                                                                     * 
 *   Prints the statistics associated to a file as given by xdfsta           *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    iun      unit number associated to the file                        * 
 *  IN    stats    statistics of the file                                    *
 *  IN    nstat    number of statistics to print                             *
 *  IN    pri      primary keys                                              *
 *  IN    aux      auxiliary keys                                            *
 *  IN    vers     software version                                          *
 *  IN    appl     application signature                                     *
 *                                                                           *
 *****************************************************************************/

int c_xdfimp(int iun,word *stat,int nstat,word_2 *pri,word_2 *aux,
                    char *vers,char *appl)
{
   int i,ind,temp;
   key_descriptor * kdp;
   word wtemp[2];
   char nomcle[5];

   ind = fnom_index(iun);
   if (ind == -1) {
      sprintf(errmsg,"file is not connected with fnom");
      return(error_msg("c_xdfimp",ERR_NO_FNOM,WARNING));
      }

   fprintf(stdout,"\n  Statistiques pour le fichier %d, Nom %s\n",iun,
                  FGFDT[ind].file_name);
   fprintf(stdout,"  Taille du fichier \t\t\t %d\n",stat[0]);
   fprintf(stdout,"  Nombre de reecriture \t\t\t %d\n",stat[1]);
   fprintf(stdout,"  Nombre d'extensions \t\t\t %d\n",stat[2]);
   fprintf(stdout,"  Nombre de pages de Repertoires \t %d\n",stat[3]);
   fprintf(stdout,"  Addr. derniere page de repertoire \t %d\n",stat[4]);
   fprintf(stdout,"  Taille du plus gros enregistrement \t %d\n",stat[5]);
   fprintf(stdout,"  Nombre d'enregistrements effaces \t %d\n",stat[10]);
   fprintf(stdout,"  Nombre d'enregistrements valides \t %d\n",stat[11]);
   fprintf(stdout,"  Version du progiciel XDF \t\t %s\n",vers);
   fprintf(stdout,"  Nom de l'application \t\t\t %s\n",appl);

   fprintf(stdout,"\n  Definition des cles primaires \n\n");
   fprintf(stdout," \tNom \tBit1 \tLongueur \tType \n\n");

   for (i=0; i < stat[6]; i++) {
      wtemp[0] = pri[i].wd1;
      wtemp[1] = pri[i].wd2;
      kdp = (key_descriptor *) &wtemp;
      temp = kdp->ncle;
      nomcle[0] = temp >> 24 & 0xff;
      nomcle[1] = temp >> 16 & 0xff;
      nomcle[2] = temp >>  8 & 0xff;
      nomcle[3] = temp & 0xff;
      nomcle[4] = '\0';
      fprintf(stdout," \t%s \t%d \t%d \t\t%d \n",
                       nomcle,kdp->bit1,kdp->lcle+1,kdp->tcle);
      }
   
   if (stat[8] > 0) {
      fprintf(stdout,"\n  Definition des cles auxiliaires \n\n");
      fprintf(stdout," \tNom \tBit1 \tLongueur \tType \n\n");
      for (i=0; i < stat[8]; i++) {
         wtemp[0] = aux[i].wd1;
         wtemp[1] = aux[i].wd2;
         kdp = (key_descriptor *) &wtemp;
         temp = kdp->ncle;
         nomcle[0] = temp >> 24 & 0xff;
         nomcle[1] = temp >> 16 & 0xff;
         nomcle[2] = temp >>  8 & 0xff;
         nomcle[3] = temp & 0xff;
         nomcle[4] = '\0';
         fprintf(stdout," \t%s \t%d \t%d \t\t%d \n",
                          nomcle,kdp->bit1,kdp->lcle+1,kdp->tcle);
         }
      } /* end if */

   return(0);
}

/*splitpoint c_xdfini */
/***************************************************************************** 
 *                          C _ X D F I N I                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Initialize the keys contain in buffer.                                  *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN      iun    unit number associated to the file                        * 
 *  IN/OUT  buf    buffer to contain the record (buf(0) contains dim of buf) *
 *  IN      idtyp  id type                                                   *
 *  IN      keys   primary keys                                              *
 *  IN      nkeys  number of primary keys                                    *
 *  IN      info   info keys                                                 *
 *  IN      ninfo  number of info keys                                       *
 *                                                                           *
 *****************************************************************************/

int c_xdfini(int iun,buffer_interface_ptr buf,int idtyp,
             word *keys,int nkeys,word *info,int ninfo)
{
   int lngbuf, i, index, index_fnom;
   file_record *record;
   file_table_entry *f;
   max_dir_keys mask;
   word *mskkeys = NULL, *buffer = (word *) buf;

   index_fnom = fnom_index(iun);
   if (index_fnom == -1) {
      sprintf(errmsg,"file is not connected with fnom");
      return(error_msg("c_xdfini",ERR_NO_FNOM,WARNING));
      }

   if ((index = file_index(iun)) == ERR_NO_FILE) {
      sprintf(errmsg,"file is not open");
      return(error_msg("c_xdfini",ERR_NO_FILE,WARNING));
      }

   lngbuf = buf->nwords;
   for (i=1; i < lngbuf; i++)
      buffer[i] = 0;
   
   if ((idtyp < 1) || (idtyp > 126)) {
      sprintf(errmsg,"invalid idtyp=%d, must be between 1 and 126",idtyp);
      return(error_msg("c_xdfini",ERR_BAD_DATYP,ERROR));
      }
   
   buf->record_index = RECADDR;
   buf->iun = iun;
   
   record = (file_record *) buf->data;
   record->idtyp = idtyp;
   
   f = file_table[index];
   buf->data_index = buf->record_index + W64TOWD(f->primary_len + f->info_len);
   buf->nbits = (f->primary_len + f->info_len) * 64;
   buf->aux_index = buf->record_index + W64TOWD(f->primary_len);

   if (nkeys > 0) f->build_primary(buf->data,keys,mask,mskkeys,index,WMODE);
   if (ninfo > 0) f->build_info(buf->data + W64TOWD(f->primary_len),
                                info,index,WMODE);

   return(0); 
}

/*splitpoint c_xdfins */
/***************************************************************************** 
 *                             C _ X D F I N S                               *
 *                                                                           * 
 *Object                                                                     * 
 *   Insert content of donnees into buf. nelm elements are inserted into     *
 *   buf starting at position bitpos.                                        *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    nelm     number of elements to add into buf                        * 
 *  IN    bitpos   bit position of insertion into buf                        *
 *  IN    donnees  data bits to add                                          *
 *  IN    nbit     number of bits kept per element                           *
 *  IN    datyp    data type                                                 *
 * IN/OUT buf      buffer to contain the modified record                     *
 *                                                                           *
 *****************************************************************************/

int c_xdfins(word *buffer, word *donnees, int bitpos,
             int nelm, int nbits, int datyp)
{
   int nbwords, index_word, last_ind, i, mode;
   buffer_interface_ptr buf = (buffer_interface_ptr) buffer;
   int  ier;

   if ((bitpos % 64) != 0) {
      sprintf(errmsg,"bitpos must be a multiple of 64");
      return(error_msg("c_xdfins",ERR_BAD_ADDR,ERRFATAL));
      }

   if (((datyp == 3) || (datyp == 5)) && (nbits != 8)) {
      sprintf(errmsg,"nbits must be 8 for datyp %d",datyp);
      return(error_msg("c_xdfins",ERR_BAD_DATYP,ERRFATAL));
      }

   nbwords = (nelm * nbits + 63) / 64;
   nbwords = W64TOWD(nbwords);

   index_word = buf->data_index + (bitpos / (sizeof(word) * 8));

   last_ind = buf->record_index + (buf->nbits / (sizeof(word) *8));

   if ((last_ind + nbwords - 1) > buf->nwords) {
      sprintf(errmsg,"buffer not big enough for insertion");
      return(error_msg("c_xdfins",ERR_BAD_DIM,ERROR));
      }

   /* move buffer content nbwords to the right for insertion */
   if (last_ind != index_word)
      for (i=last_ind; i >= index_word; i--)
         buf->data[i+nbwords] = buf->data[i];

   /* insert data */
   switch (datyp)
      {
      case 0:         /* transparent mode */
      case 3:

         for (i=0; i < nbwords; i++)
            buf->data[index_word+i] = donnees[i];
         break;

#if !defined(NEC)
      case 6:
      case 8:

         for (i=0; i < nbwords; i++)
            buf->data[index_word+i] = donnees[i];
         break;

      case 7:
      case 9:

         if (*little_endian) 
           for (i=0; i < nbwords; i+=2) {
             buf->data[index_word+i] = donnees[i+1];
             buf->data[index_word+i+1] = donnees[i];
           }
         else
           for (i=0; i < nbwords; i++)
             buf->data[index_word+i] = donnees[i];
         break;
#endif
      
      case 5:        /* upper char only */
         for (i=0; i < nbwords; i++)
            buf->data[index_word+i] = upper_case_word(donnees[i]);
         break;
      
      case 2: mode = 1;
         ier = compact_integer(donnees,(void *) NULL,&(buf->data[index_word]),nelm,
                               nbits,0,xdf_stride,mode);
         break;

      case 4: mode = 3;
         ier = compact_integer(donnees,(void *) NULL,&(buf->data[index_word]),nelm,
                               nbits,0,xdf_stride,mode);
         break;

      default: sprintf(errmsg,"invalid datyp=%d",datyp);
         return(error_msg("c_xdfins",ERR_BAD_DATYP,ERROR));

      } /* end switch */

   buf->nbits += nbwords * sizeof(word) * 8;
   return(0);
}

/*splitpoint c_xdflnk */
/***************************************************************************** 
 *                             C _ X D F L N K                               * 
 *                                                                           * 
 *Object                                                                     * 
 *   Links the list of random files together for record search purpose.      *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    liste    unit number associated to the file                        * 
 *  IN    n        number of files to be linked                              *
 *                                                                           *
 *****************************************************************************/

int c_xdflnk(word *liste, int n)
{
  int index, indnext, i, index_fnom;
  file_table_entry *f, *fnext;

  index_fnom = fnom_index(liste[0]);
  if (index_fnom == -1) {
    sprintf(errmsg,"file is not connected with fnom");
    return(error_msg("c_xdflnk",ERR_NO_FNOM,ERROR));
  }
  
  if ((index = file_index(liste[0])) == ERR_NO_FILE) {
    sprintf(errmsg,"file is not open");
    return(error_msg("c_xdflnk",ERR_NO_FILE,ERROR));
  }

  f = file_table[index];
  for (i=1; i < n; i++) {
    if ((index_fnom = fnom_index(liste[i])) == -1) {
      sprintf(errmsg,"file is not connected with fnom");
      return(error_msg("c_xdflnk",ERR_NO_FNOM,ERROR));
    }
    if ((indnext = file_index(liste[i])) == ERR_NO_FILE) {
      sprintf(errmsg,"file is not open");
      return(error_msg("c_xdflnk",ERR_NO_FILE,ERROR));
    }
    if (msg_level <= TRIVIAL)
      fprintf(stdout,"Debug xdflink %d avec %d\n",liste[i-1],liste[i]);
    fnext = file_table[indnext];
    f->link = indnext;
    (f->dir_page[f->npages-1])->next_page = fnext->dir_page[0];
    index = indnext;
    f = file_table[index];
  }
   
  return(0);
}

/*splitpoint c_xdfloc */
/***************************************************************************** 
 *                             C _ X D F L O C                               * 
 *                                                                           * 
 *Object                                                                     * 
 *   Find the position of the record as described by the given primary keys. *
 *   The search begins from the record pointed by handle.  If handle is 0,   *
 *   the search is from beginning of file. If handle is -1, the search       *
 *   begins from the current position. If a specific key has the value of    *
 *   -1, this key will not used as a selection criteria.                     *
 *   If nprim is -1, the last selection criterias will be used for the       *
 *   search. Upon completion a "pointer" to the record (handle) is returned. *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    iun      unit number associated to the file                        * 
 *  IN    handle   handle to the starting search position                    *
 *  IN    primk    primary search keys                                       *
 *  IN    nprim    number search primary search keys                         *
 *                                                                           *
 *****************************************************************************/

int c_xdfloc(int iun, int handle, word *primk,int nprim)
{
  word *mskkeys = NULL;
  int ier;
  
  ier = c_xdfloc2(iun,handle,primk,nprim,mskkeys);
  return(ier);
}

/*splitpoint c_xdfloc2 */
/***************************************************************************** 
 *                             C _ X D F L O C 2                             * 
 *                                                                           * 
 *Object                                                                     * 
 *   Find the position of the record as described by the given primary keys. *
 *   The search begins from the record pointed by handle.  If handle is 0,   *
 *   the search is from beginning of file. If handle is -1, the search       *
 *   begins from the current position. If a specific key has the value of    *
 *   -1, this key will not used as a selection criteria.                     *
 *   If nprim is -1, the last selection criterias will be used for the       *
 *   search. Upon completion a "pointer" to the record (handle) is returned. *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    iun      unit number associated to the file                        * 
 *  IN    handle   handle to the starting search position                    *
 *  IN    primk    primary search keys                                       *
 *  IN    nprim    number search primary search keys                         *
 *  IN    mskkeys  search mask                                               *
 *                                                                           *
 *****************************************************************************/

int c_xdfloc2(int iun, int handle, word *primk,int nprim, word *mskkeys)
{
   int i,record,pageno;
   int index, index_h, was_allocated, new_handle;
   word *pmask, *psmask;
   file_table_entry *f;
   xdf_record_header header;
   seq_dir_keys seq_entry;

   if ((index = file_index(iun)) == ERR_NO_FILE) {
     sprintf(errmsg,"file is not open,iun=%d\n",iun);
     return(error_msg("c_xdfloc",ERR_NO_FILE,ERROR));
   }
   f = file_table[index];

   was_allocated = 0;
   if ((mskkeys == NULL) && (nprim != 0)) {
     mskkeys = malloc(nprim *sizeof(word));
     was_allocated = 1;
     for (i = 0; i < nprim; i++) {
       if (*primk == -1)
         *mskkeys = 0;
       else
         *mskkeys = -1;
     }
   }

   if (handle > 0) {          /* search begins from handle given position */
      index_h = INDEX_FROM_HANDLE(handle);
      record = RECORD_FROM_HANDLE(handle)+1;   /* search from next record */
      pageno = PAGENO_FROM_HANDLE(handle);
      if (index_h != index) {
         sprintf(errmsg,"invalid handle=%d or iun=%d\n",handle,iun);
         return(error_msg("c_xdfloc",ERR_BAD_HNDL,ERROR));
         }
      if (f->xdf_seq) {
        f->cur_addr = address_from_handle(handle,f);
        if (f->fstd_vintage_89) {
          c_waread(iun,&seq_entry,f->cur_addr,sizeof(seq_entry)/bytesperword);
          header.lng = ((seq_entry.lng + 3) >> 2) + 15; 
         }
        else
          c_waread(iun,&header,f->cur_addr,W64TOWD(1));
        f->cur_addr += W64TOWD(header.lng);
      }
   }
   else if (handle == 0) {    /* search from beginning of file */
     record = 0;
     pageno = 0;
     f->cur_pageno = -1;
     if (f->xdf_seq) 
       f->cur_addr = f->seq_bof;
   }
   else if (handle == -1) {   /* search from current position */
     if (((f->cur_entry == NULL) || (f->cur_pageno == -1))
         && (! f->xdf_seq)) {
       sprintf(errmsg,"current file position is invalid\n");
       return(error_msg("c_xdfloc",ERR_NO_POS,ERROR));
     }
   }
   else {
     sprintf(errmsg,"invalid handle\n");
     return(error_msg("c_xdfloc",ERR_BAD_HNDL,ERROR));
   }
   
   if (nprim) {       /* if nprim == 0 keep same search target */
     f->build_primary(f->target,primk,f->cur_mask,mskkeys,index,WMODE);
     pmask = (word *) f->cur_mask;
     psmask = (word *) f->srch_mask;
     for (i=0; i < W64TOWD(f->primary_len); i++, pmask++, psmask++)
       *pmask &= *psmask;
     f->valid_target = 1;
   }

   if ((handle != -1) && (! f->xdf_seq)) {
     if (pageno != f->cur_pageno) {
       if (pageno < f->npages) {   /* page is in current file */
         f->cur_dir_page = f->dir_page[pageno];
         f->cur_pageno = pageno;
       }
       else {                      /* page is in a link file */
         if (f->link == -1) {
          sprintf(errmsg,"page number=%d > last page=%d and file not linked\n",
                       pageno,f->npages-1);
          f->cur_entry = NULL;
          return(error_msg("c_xdfloc",ERR_BAD_PAGENO,ERROR));
         }
         f->cur_dir_page = f->dir_page[f->npages-1];
         f->cur_pageno = f->npages-1;
         for (i=0; (i<= pageno - f->cur_pageno) && (f->cur_dir_page); i++) {
           f->cur_dir_page = (f->cur_dir_page)->next_page;
           f->cur_pageno++;
         }
         if (f->cur_dir_page == NULL) {
           sprintf(errmsg,"invalid handle, invalid page number\n");
           f->cur_entry = NULL;
           return(error_msg("c_xdfloc",ERR_BAD_PAGENO,ERROR));
         }
       }
     }
     else
       f->cur_dir_page = f->dir_page[f->cur_pageno]; /* just to make sure */
     f->cur_entry = (f->cur_dir_page)->dir.entry + 
       record * W64TOWD(f->primary_len);
     f->page_record = record;
   }

   if (((f->cur_entry == NULL) || (f->cur_pageno == -1)) && (! f->xdf_seq)) {
     sprintf(errmsg,"no valid current file position\n");
     return(error_msg("c_xdfloc",ERR_NO_POS,ERROR));
   }
   
   if (! f->valid_target) {
     sprintf(errmsg,"no valid current search target\n");
     return(error_msg("c_xdfloc",ERR_NO_TARGET,ERROR));
   }
   
   new_handle = next_match(index);   
   if (was_allocated) free(mskkeys);
   return(new_handle);
}

/*splitpoint c_xdfopn */
/*****************************************************************************
 *                          C _ X D F O P N                                  *
 *                                                                           *
 *Object                                                                     *
 *   Open an XDF file.                                                       *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN    iun    unit number associated to the file                          *
 *  IN    mode   open mode (READ, WRITE, R-W, CREATE, APPEND)                *
 *  IN    keys   primary keys                                                *
 *  IN    nkeys  number of primary keys                                      *
 *  IN    info   info keys                                                   *
 *  IN    ninfo  number of info keys                                         *
 *  IN    appl   application signature                                       *
 *                                                                           *
 *****************************************************************************/

int c_xdfopn(int iun,char *mode,word_2 *pri,int npri,
             word_2 *aux,int naux,char *appl)
{
  int index,index_fnom,ier,i,j,nrec=0;
  file_table_entry *f;
  ftnword f_datev;
  double nhours;
  int deet,npas,i_nhours,run,datexx;
  word STDR_sign = 'S' << 24 | 'T' << 16 | 'D' << 8 | 'R';
  word STDS_sign = 'S' << 24 | 'T' << 16 | 'D' << 8 | 'S';

  if (!init_package_done) {
     init_package();
     init_package_done = 1;
     }

  if ((iun <= 0) || (iun > 999)) {
     sprintf(errmsg,"invalid unit number=%d",iun);
     return(error_msg("c_xdfopn",ERR_BAD_UNIT,ERROR));
     }

  if (file_index(iun) != ERR_NO_FILE) {
     sprintf(errmsg,"file (unit=%d) is already open",iun);
     return(error_msg("c_xdfopn",ERR_FILE_OPN,WARNING));
     }

  index = get_free_index();
  file_table[index]->iun = iun;
  file_table[index]->file_index = index;

  index_fnom = fnom_index(iun);
  if (index_fnom == -1) {
     sprintf(errmsg,"file (unit=%d) is not connected with fnom",iun);
     return(error_msg("c_xdfopn",ERR_NO_FNOM,ERROR));
     }

  f = file_table[index];
  f->cur_info = &Fnom_General_File_Desc_Table[index_fnom];

  if (! f->cur_info->attr.rnd) {
     sprintf(errmsg,"file must be random\n file in error: %s\n",FGFDT[index_fnom].file_name);
     return(error_msg("c_xdfopn",-1,ERROR));
     }

  if (strstr(appl,"BRP0")) f->cur_info->attr.burp = 1;
  if (strstr(appl,"STD"))  f->cur_info->attr.std = 1;

  if (f->cur_info->attr.burp) {
     f->build_primary = (fn_b_p *) build_burp_prim_keys;
     f->build_info = (fn_ptr *) build_burp_info_keys;
     }
  else {
    if (f->cur_info->attr.std) {
      f->build_primary = (fn_b_p *) build_fstd_prim_keys;
      f->build_info = (fn_ptr *) build_fstd_info_keys;
      f->file_filter = (fn_ptr *) C_fst_match_req;
    }
    else {
    f->build_primary = (fn_b_p *) build_gen_prim_keys;
    f->build_info = (fn_ptr *) build_gen_info_keys;
    }
  }

  if ((strstr(f->cur_info->file_type,"SEQ")) || (strstr(f->cur_info->file_type,"seq"))) {
    f->xdf_seq = 1;
    STDSEQ_opened = 1;             /* at least one seq file is opened, limit number of xdf files is now 128 */
    if (index > 127) {
      sprintf(errmsg,"while opening std/seq file, limit of 128 opened file reached");
      return(error_msg("c_xdfopn",-1,ERROR));
    }
  }


  if (msg_level <= TRIVIAL)
    fprintf(stdout,"Debug c_xdfopn f->xdf_seq=%d\n",f->xdf_seq);

  if (strstr(mode,"CREATE") || strstr(mode,"create")) {
/*
 *  create new xdf file
 */
     c_waopen(iun);
     ier = create_new_xdf(index,iun,pri,npri,aux,naux,appl);
     if (! f->xdf_seq)
       add_dir_page(index,WMODE);
     else {
       f->cur_addr = f->nxtadr;
       f->header->fsiz = WDTO64(f->nxtadr -1);
     }
     f->header->rwflg = CREATE;
     f->modified = 1;
     return(ier);
     }
/*
 * file exist, read directory pages
 */
  if (! FGFDT[index_fnom].open_flag)
     c_waopen(iun);

  { int unit=iun, wdaddress=1 , wdlng_header, lng64; file_record header64;
    word32 *check32, checksum;
    int header_seq[30];
    int lng;
    stdf_struct_RND header_rnd;
    rnd_dir_keys *directory;
    stdf_dir_keys *stdf_entry;
    xdf_dir_page * curpage;

    c_waread(unit,&header64,wdaddress,W64TOWD(2));
    if (header64.data[0] == 'XDF0' || header64.data[0] == 'xdf0') {
    /*if (strncmp(&header64.data[0], "XDF0", 4) == 0 || strncmp(&header64.data[0], "xdf0", 4) == 0) {*/
      if ((f->header = calloc(1,header64.lng*8)) == NULL) {
        sprintf(errmsg,"memory is full");
        return(error_msg("c_xdfopn",ERR_MEM_FULL,ERRFATAL));
      }
/*
 * read file header
 */
      wdlng_header = W64TOWD(header64.lng);
      c_waread(unit,f->header,wdaddress,wdlng_header);
      f->primary_len = f->header->lprm;
      f->info_len = f->header->laux;
      f->nxtadr = W64TOWD(f->header->fsiz) + 1;    /* nxtadr = fsiz +1 */
      wdaddress += wdlng_header;
      if (f->cur_info->attr.std)
        if ((f->header->sign != STDR_sign) && (f->header->sign != STDS_sign)) {
          sprintf(errmsg,"%s is not a standard file",FGFDT[index_fnom].file_name);
          return(error_msg("c_xdfopn",ERR_WRONG_FTYPE,ERRFATAL));
          }
      if (strstr(mode,"READ") || strstr(mode,"read")) {
        f->header->rwflg = RDMODE;
      }
      else {
        if (f->header->rwflg != RDMODE) {
          sprintf(errmsg,
          "file (unit=%d) currently used by another application in write mode",
                  unit);
          return(error_msg("c_xdfopn",ERR_STILL_OPN,ERRFATAL));
        }
        if (strstr(mode,"WRITE") || strstr(mode,"write")) {
          f->header->rwflg = WMODE;
        }
        else
          if (strstr(mode,"R-W") || strstr(mode,"r-w")) {
            f->header->rwflg = RWMODE;
          }
          else
            if (strstr(mode,"APPEND") || strstr(mode,"append")) {
              f->header->rwflg = APPEND;
            }
      }

      if (f->header->nbd ==0) {
        if ( (f->cur_info->attr.std) && (header64.data[1] == 'STDR' || header64.data[1] == 'stdr') ) {
          sprintf(errmsg,"File probably damaged\n file in error: %s\n",FGFDT[index_fnom].file_name);
          return(error_msg("c_xdfopn",ERR_BAD_DIR,ERROR));
          }
        else
          f->xdf_seq = 1;
        }
      else
        if (f->xdf_seq == 1)    /* random file opened in seqential mode */
          f->header->rwflg = RDMODE;

      if (msg_level <= TRIVIAL)
        fprintf(stdout,"Debug c_xdfopn fichier existe f->xdf_seq=%d\n",
                f->xdf_seq);

      if (! f->xdf_seq) {
/*
 * read directory pages and compute checksum
 */
        for (i=0; i < f->header->nbd; i++) {
          add_dir_page(index,RDMODE);
          lng64 = f->primary_len * ENTRIES_PER_PAGE + 4;
          curpage = &((f->dir_page[f->npages-1])->dir);
          c_waread(unit,curpage,wdaddress,W64TOWD(lng64));
          checksum = 0;
          check32 = (word32 *) curpage;
          for (j=4; j < W64TOWD(lng64); j++)
            checksum ^= check32[j];
          if (checksum != 0) {
            sprintf(errmsg,"incorrect checksum in page %d, directory is probably damaged\n file in error: %s\n",i,FGFDT[index_fnom].file_name);
          /*          return(error_msg("c_xdfopn",ERR_BAD_CHKS,ERROR)); */
            error_msg("c_xdfopn",ERR_BAD_CHKS,ERROR);
          }
          wdaddress = W64TOWD(curpage->nxt_addr - 1) +1;
          if (((wdaddress == 0) && (i != f->header->nbd-1)) || (wdaddress > FGFDT[index_fnom].file_size)) {
            sprintf(errmsg,"number of directory pages is incorrect\n file in error: %s\n",FGFDT[index_fnom].file_name);
            return(error_msg("c_xdfopn",ERR_BAD_DIR,ERROR));
          }
          nrec += curpage->nent;
        } /* end for */

        f->nrecords = nrec;
      }
      else {  /* file is xdf sequential, position address to first record */
        f->cur_addr = wdaddress;
        f->seq_bof = wdaddress;
      }
    }
    else {                              /* signature != XDF0 */     
      check32 = (word32 *) &header64;
      if (*check32 == STDF_RND_SIGN) {        /* old random standard file */
        f->cur_info->attr.read_only = 1;
        f->fstd_vintage_89 = 1;
        lng = sizeof(header_rnd) / sizeof(word);
        c_waread(iun,&header_rnd,wdaddress,lng);
        wdaddress += lng;
        if ((directory =
             calloc(header_rnd.nutil,
                    sizeof(word)*sizeof(rnd_dir_keys))) == NULL) {
          sprintf(errmsg,"memory is full");
          return(error_msg("c_xdfopn",ERR_MEM_FULL,ERRFATAL));
        }
        lng = header_rnd.nutil * sizeof(rnd_dir_keys) / sizeof(word); 
        c_waread(iun,directory,wdaddress,lng);
        create_new_xdf(index,iun,(word_2 *)&stdfkeys,16,aux,0,"STDF");
        add_dir_page(index,RDMODE);
        f->cur_dir_page = f->dir_page[f->npages-1];
        for (i=0; i < header_rnd.nutil; i++) {
          if (f->cur_dir_page->dir.nent >= ENTRIES_PER_PAGE) {
            f->nrecords += f->page_nrecords;
            add_dir_page(index,RDMODE);
            f->cur_dir_page = f->dir_page[f->npages-1];
          }
          f->cur_entry = f->cur_dir_page->dir.entry +
                         f->cur_dir_page->dir.nent * W64TOWD(f->primary_len);
          stdf_entry = (stdf_dir_keys *) f->cur_entry;
          f->page_nrecords = ++f->cur_dir_page->dir.nent;
          if (! directory[i].dltf) {
            stdf_entry->deleted = 0;
            stdf_entry->select = 1;
            stdf_entry->lng = (directory[i].lng + 3) >> 2;
            stdf_entry->addr = (directory[i].swa >> 2) +1;
            stdf_entry->deet = directory[i].deet;
            stdf_entry->nbits = directory[i].nbits;
            stdf_entry->ni = directory[i].ni;
            stdf_entry->gtyp = directory[i].grtyp;
            stdf_entry->nj = directory[i].nj;
            stdf_entry->datyp = directory[i].datyp;
            stdf_entry->nk = directory[i].nk;
            stdf_entry->ubc = 0;
            stdf_entry->npas = (directory[i].npas2 << 16) | 
                               directory[i].npas1;
            stdf_entry->pad7 = 0;
            stdf_entry->ig4 = directory[i].ig4;
            stdf_entry->ig2a = 0;
            stdf_entry->ig1 = directory[i].ig1;
            stdf_entry->ig2b = directory[i].ig2 >> 8;
            stdf_entry->ig3 = directory[i].ig3;
            stdf_entry->ig2c = directory[i].ig2 & 0xff;
            stdf_entry->etik15 = 
              (ascii6(directory[i].etiq14 >> 24) << 24) |
              (ascii6((directory[i].etiq14 >> 16) & 0xff) << 18) |
              (ascii6((directory[i].etiq14 >>  8) & 0xff) << 12) |
              (ascii6((directory[i].etiq14      ) & 0xff) <<  6) |
               (ascii6((directory[i].etiq56 >>  8) & 0xff));
            stdf_entry->pad1 = 0;
            stdf_entry->etik6a = 
              (ascii6((directory[i].etiq56      ) & 0xff) << 24) |
              (ascii6((directory[i].etiq78 >>  8) & 0xff) << 18) |
              (ascii6((directory[i].etiq78      ) & 0xff) << 12);
            stdf_entry->pad2 = 0;
            stdf_entry->etikbc = 0;
            stdf_entry->typvar = ascii6(directory[i].typvar) << 6;
            stdf_entry->pad3 = 0;
            stdf_entry->nomvar = 
               (ascii6((directory[i].nomvar >>  8) & 0xff) << 18) |
               (ascii6((directory[i].nomvar      ) & 0xff) << 12);
            stdf_entry->pad4 = 0;
            stdf_entry->ip1 = directory[i].ip1;
            stdf_entry->levtyp = 0;
            stdf_entry->ip2 = directory[i].ip2;
            stdf_entry->pad5 = 0;
            stdf_entry->ip3 = directory[i].ip3;
            stdf_entry->pad6 = 0;
            stdf_entry->date_stamp = directory[i].date;
            deet = stdf_entry->deet;
            npas = stdf_entry->npas;
            if (((deet*npas) % 3600) != 0) {
              /*
               *  recompute datev to take care of rounding used with 1989 version
               *  de-octalise the date_stamp
               */
              run = stdf_entry->date_stamp & 0x7;
              datexx = (stdf_entry->date_stamp >> 3) * 10 + run;
              
              f_datev = (ftnword) datexx;
              i_nhours = (deet*npas - ((deet*npas+1800)/3600)*3600);
              nhours = (double) (i_nhours / 3600.0);
              f77name(incdatr)(&f_datev,&f_datev,&nhours);
              datexx = (int) f_datev;
              /*
               *  re-octalise the date_stamp
               */
              stdf_entry->date_stamp = 8 * (datexx/10) + (datexx % 10);
            }
          }
          else {
            stdf_entry->deleted = 1;
            /*            fprintf(stdout,"Debug c_xdfopn i=%d is deleted\n",i); */
          }
        } /* end for */
        f->nrecords += f->page_nrecords;
        free(directory);
      } /* end if random */
      else {
        c_waread(iun,&header_seq,1,30);
        if (header_seq[896/32] == STDF_SEQ_SIGN) { /* old sequential stdf */
          f->cur_info->attr.read_only = 1;
          f->fstd_vintage_89 = 1;
          f->xdf_seq = 1;
          create_new_xdf(index,iun,(word_2 *)&stdfkeys,16,aux,0,"STDF");
          f->cur_addr = 1;
          f->seq_bof = 1;
          return(0);
        }
        else {
          sprintf(errmsg,"file is not XDF type or old standard random type\n");
          return(error_msg("c_xdfopn",ERR_NOT_XDF,ERRFATAL));
        }
      }
    }
  }
  return(f->header->nrec);
}


/*splitpoint c_xdfopt */
/***************************************************************************** 
 *                            C _ X D F O P T                                * 
 *                                                                           * 
 *Object                                                                     * 
 *   Set different options in xdf.                                           *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    optname  name of option to set                                     * 
 *  IN    optc     value of option if type is character                      *
 *  IN    optv     value of option if type is integer                        *
 *                                                                           * 
 *****************************************************************************/

int c_xdfopt(char *optname, char *optc, int optv)
{
   if (strstr(optname,"ERRTOLR") || strstr(optname,"errtolr")) {

      if (strstr(optc,"TRIVIAL") || strstr(optc,"trivial"))
         xdf_toler = TRIVIAL;
      else if (strstr(optc,"INFORMATIF") || strstr(optc,"informatif"))
         xdf_toler = INFORM;
      else if (strstr(optc,"WARNING") || strstr(optc,"warning"))
         xdf_toler = WARNING;
      else if (strstr(optc,"ERROR") || strstr(optc,"error"))
         xdf_toler = ERROR;
      else if (strstr(optc,"FATAL") || strstr(optc,"fatal"))
         xdf_toler = ERRFATAL;
      else if (strstr(optc,"SYSTEM") || strstr(optc,"SYSTEM")) 
         xdf_toler = SYSTEM;
      else {
         sprintf(errmsg,"invalid option value: %s",optc);
         return(error_msg("c_xdfopt",ERR_BAD_OPT,ERROR));
         }
      
      }
   else if (strstr(optname,"MSGLVL") || strstr(optname,"msglvl")) {

      if (strstr(optc,"TRIVIAL") || strstr(optc,"trivial"))
         msg_level = TRIVIAL;
      else if (strstr(optc,"INFORMATIF") || strstr(optc,"informatif"))
         msg_level = INFORM;
      else if (strstr(optc,"WARNING") || strstr(optc,"warning"))
         msg_level = WARNING;
      else if (strstr(optc,"ERROR") || strstr(optc,"error"))
         msg_level = ERROR;
      else if (strstr(optc,"FATAL") || strstr(optc,"fatal"))
         msg_level = ERRFATAL;
      else if (strstr(optc,"SYSTEM") || strstr(optc,"SYSTEM"))
         msg_level = SYSTEM;
      else {
         sprintf(errmsg,"invalid option value: %s",optc);
         return(error_msg("c_xdfopt",ERR_BAD_OPT,ERROR));
         }
      }

   else if (strstr(optname,"STRIPING") || strstr(optname,"striping")) {
      xdf_nsplit = optv;
      }
   else {
      sprintf(errmsg,"invalid option name: %s",optname);
      return(error_msg("c_xdfopt",ERR_BAD_OPT,ERROR));
      }
   return(0);

}

/*splitpoint c_xdfprm */
/*****************************************************************************
 *                            C _ X D F P R M                                *
 *                                                                           *
 *Object                                                                     *
 *   Get the descriptive parameters of the record pointed by handle          *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    handle   buffer that contains the record                           * 
 *  OUT   addr     record starting address                                   *
 *  OUT   lng      length of the record in 64 bit units                      *
 *  OUT   idtyp    record id type                                            *
 *  OUT   primk    primary keys                                              *
 *  IN    nprim    number of primary keys                                    *
 *                                                                           * 
 *****************************************************************************/

int c_xdfprm(int handle,int *addr,int *lng,int *idtyp,word *primk,int nprim)
{
   int index, record_number, page_number, i;
   file_table_entry *f;
   file_record *record;
   word *rec;
   max_dir_keys argument_not_used;
   word *mskkeys = NULL;
   page_ptr target_page;

   index =         INDEX_FROM_HANDLE(handle);
   page_number =   PAGENO_FROM_HANDLE(handle);
   record_number = RECORD_FROM_HANDLE(handle);

/*   printf("Debug+ c_xdfprm index=%d page_number=%d record_number=%d\n",index,page_number,record_number); */
   /* validate index, page number and record number */

   if ((file_table[index] == NULL) || (file_table[index]->iun < 0)) {
      sprintf(errmsg,"invalid handle, invalid file index\n");
      return(error_msg("c_xdfprm",ERR_BAD_HNDL,ERROR));
      }

   f = file_table[index];


   if (! f->xdf_seq) {
     if (page_number < f->npages) {   /* page is in current file */
       target_page = f->dir_page[page_number];
     }
     else {                           /* page is in a link file */
       if (f->link == -1) {
         sprintf(errmsg,"page number=%d > last page=%d and file not linked\n",
                 page_number,f->npages-1);
           return(error_msg("c_xdfprm",ERR_BAD_PAGENO,ERROR));
       }
       target_page = f->dir_page[f->npages-1];
       for (i=0; (i<= (page_number - f->npages)) && target_page; i++)
         target_page = (target_page)->next_page;
       if (target_page == NULL) {
         sprintf(errmsg,"invalid handle, invalid page number\n");
         return(error_msg("c_xdfprm",ERR_BAD_PAGENO,ERROR));
       }
     }
     
     if (record_number > target_page->dir.nent) {
      sprintf(errmsg,"invalid handle, invalid record number\n");
      return(error_msg("c_xdfprm",ERR_BAD_HNDL,ERROR));
     }
     
     rec = target_page->dir.entry + record_number * W64TOWD(f->primary_len);
     record = (file_record *) rec;
   }
   else {
     if (! f->valid_pos) {
       sprintf(errmsg,"no valid file position for sequential file\n");
       return(error_msg("c_xdfprm",ERR_NO_POS,ERROR));
     }
     record = (file_record *) f->head_keys;
     if (address_from_handle(handle,f) != W64TOWD(record->addr-1)+1) {
       sprintf(errmsg,
               "invalid handle=%d, invalid address=%d record address=%d\n",
               handle,address_from_handle(handle,f),W64TOWD(record->addr-1)+1);
       return(error_msg("c_xdfprm",ERR_BAD_HNDL,ERROR));
     }
   }
  
   *idtyp = record->idtyp;
   *addr = record->addr;
   *lng = record->lng;

   f->build_primary((word *) record,primk,argument_not_used,mskkeys,index,RDMODE);

   return(0);

}

/*splitpoint c_xdfput */
/***************************************************************************** 
 *                            C _ X D F P U T                                * 
 *                                                                           * 
 *Object                                                                     * 
 *   Write record (from buf) to xdf file. If handle is not 0, rewrite        * 
 *   record referenced by handle. If handle is negative, the record will be  *
 *   append to end of file.                                                  *
 *   If handle is 0, the record is added at the end of file.                 *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           *
 *  IN    iun     unit number of the file to be written                      *
 *  IN    handle  file index, page number and record number to record        * 
 *  IN    buf     buffer to contain record                                   * 
 *                                                                           * 
 *****************************************************************************/

int c_xdfput(int iun, int handle, buffer_interface_ptr buf)
{
   int index, record_number, page_number, i, idtyp, addr, lng, lngw;
   int index_from_buf, index_from_iun, write_to_end=0, nwords;
   int write_addr, err, index_fnom, cluster_size, next_cluster_addr;
   file_table_entry *f, *f_buf;
   file_record *record, *bufrec;
   max_dir_keys primk, argument_not_used, mskkeys;
   xdf_record_header header64;
   burp_record *debug_record;
   postfix_seq postfix;

   index_fnom = fnom_index(iun);
   if ((index_from_buf = file_index(buf->iun)) == ERR_NO_FILE) {
      sprintf(errmsg,"record not properly initialized\n");
      return(error_msg("c_xdfput",ERR_BAD_INIT,ERROR));
      }

   if ((index_from_iun = file_index(iun)) == ERR_NO_FILE) {
      sprintf(errmsg,"invalid iun (%d)\n",iun);
      return(error_msg("c_xdfput",ERR_BAD_UNIT,ERROR));
      }

   if ((buf->nbits & 0x3f) != 0) {
      sprintf(errmsg,"buf->nbits is not a multiple of 64 bits\n");
      return(error_msg("c_xdfput",ERR_BAD_ADDR,SYSTEM));
      }

   nwords = buf->nbits / (8 * sizeof(word));
   index = index_from_iun;
   f = file_table[index_from_iun];

   if ((f->header->rwflg == RDMODE) || (FGFDT[index_fnom].attr.read_only)) {  /* mode read only */
      sprintf(errmsg,"file is open in read only mode or no write permission\n");
      return(error_msg("c_xdfput",ERR_NO_WRITE,ERROR));
      }

   if ((handle != 0) && (f->header->rwflg == APPEND)) {
      sprintf(errmsg,"file is open in append mode only\n");
      return(error_msg("c_xdfput",ERR_NO_WRITE,ERROR));
      }

   if (handle <= 0)
      write_to_end = 1;

   if (handle != 0) {
      if (handle < 0 ) 
         handle = -handle;

      index =         INDEX_FROM_HANDLE(handle);
      page_number =   PAGENO_FROM_HANDLE(handle);
      record_number = RECORD_FROM_HANDLE(handle);

      /* validate index, page number and record number */

      if (index != index_from_iun) {
         sprintf(errmsg,"iun and handle do not match\n");
         return(error_msg("c_xdfput",ERR_BAD_HNDL,ERROR));
         }

      /* make sure that files are the same type */

      if (index_from_iun != index_from_buf) {
         f_buf = file_table[index_from_buf];
         if ((f->header->nprm != f_buf->header->nprm) ||
             (f->header->lprm != f_buf->header->lprm) ||
             (f->header->naux != f_buf->header->naux) ||
             (f->header->laux != f_buf->header->laux)) {
            sprintf(errmsg,"source and destination files are different type of file\n");
            return(error_msg("c_xdfput",ERR_NOT_COMP,ERROR));
            }
         }
      if (! f->xdf_seq) {
        if (page_number < f->npages) {   /* page is in current file */
          f->cur_dir_page = f->dir_page[page_number];
        }
        else {                           /* page is in a link file */
          if (f->link == -1) {
          sprintf(errmsg,"page number=%d > last page=%d and file not linked\n",
                    page_number,f->npages-1);
            return(error_msg("c_xdfput",ERR_BAD_PAGENO,ERROR));
          }
          f->cur_dir_page = f->dir_page[f->npages-1];
          for (i=0; ((i<= page_number - f->npages) && f->cur_dir_page); i++)
            f->cur_dir_page = (f->cur_dir_page)->next_page;
          if (f->cur_dir_page == NULL) {
            sprintf(errmsg,"invalid handle, invalid page number\n");
            return(error_msg("c_xdfput",ERR_BAD_PAGENO,ERROR));
          }
        }

        if (record_number > f->cur_dir_page->dir.nent) {
          sprintf(errmsg,"invalid handle, invalid record number\n");
          return(error_msg("c_xdfput",ERR_BAD_HNDL,ERROR));
        }
        
        f->cur_entry =
          f->cur_dir_page->dir.entry + record_number * W64TOWD(f->primary_len);
   
        record = (file_record *) f->cur_entry;
      }
      else {   /* file is xdf sequential */
        if (handle > 0) {   /* rewrite record */
          addr = address_from_handle(handle,f);
          c_waread(iun,f->head_keys,addr,MAX_PRIMARY_LNG);
          record = (file_record *) f->head_keys;
        }
      }
      idtyp = record->idtyp;
      addr = record->addr;
      lng = record->lng;
      lngw = W64TOWD(lng);

      if (idtyp == 0) {
         sprintf(errmsg,"special record idtyp=0\n");
         return(error_msg("c_xdfput",ERR_SPECIAL,ERROR));
         }

      if ((idtyp & 0x7E) == 0x7E) {
         sprintf(errmsg,"deleted record\n");
         return(error_msg("c_xdfput",ERR_DELETED,ERROR));
         }

      if (lngw != nwords)      /* enforce rewrite to end of file */
         write_to_end = 1;
      else
         write_addr = W64TOWD(addr-1)+1;
      
      } /* end if handle != 0 */
           
   if ((write_to_end) && (! f->xdf_seq)) {
      f->cur_dir_page = f->dir_page[f->npages-1];
      f->cur_pageno = f->npages-1;
      if (f->cur_dir_page->dir.nent >= ENTRIES_PER_PAGE) {
         f->cur_dir_page->modified = 1;
         err = add_dir_page(index,WMODE);
         if (err < 0) return(err);
         f->cur_dir_page = f->dir_page[f->npages-1];
         f->cur_entry = f->cur_dir_page->dir.entry;
         }
      else {
         f->cur_entry = f->cur_dir_page->dir.entry +
                        f->cur_dir_page->dir.nent * W64TOWD(f->primary_len);
         }
      f->page_nrecords = f->cur_dir_page->dir.nent++;
      write_addr = f->nxtadr;
      }

   if ((write_to_end) && (f->xdf_seq))
     write_addr = f->cur_addr;

   if (handle != 0) {      /* rewrite, delete old record */
      err = c_xdfdel(handle);
      if (err < 0) return(err);
      }


   /* update record header */
   bufrec = (file_record *) buf->data;
   bufrec->addr = WDTO64(write_addr - 1) + 1;
   if (f->xdf_seq) {
     next_cluster_addr = f->cur_addr -1 + nwords + W64TOwd(2);
     if ((next_cluster_addr >> 18) >= 512) {
       cluster_size = 128;
     }
     else if ((next_cluster_addr >> 18) >= 128) {
       cluster_size = 32;
     }
     else if ((next_cluster_addr >> 18) >= 32) {
       cluster_size = 8;
     }
     else
       cluster_size = 2;
     next_cluster_addr = ((next_cluster_addr + cluster_size-1) / cluster_size)
                         * cluster_size;
     nwords = next_cluster_addr - f->cur_addr - W64TOwd(2) +1;
   }
   bufrec->lng = WDTO64(nwords);

   /* write record to file */
   c_wawrit(iun,&(buf->data),write_addr,nwords);
   if (f->xdf_seq)
     f->cur_addr += nwords;

   if (! f->xdf_seq) {
     /* update directory entry */
     debug_record = (burp_record *) f->cur_entry;
     f->build_primary(buf->data,primk,argument_not_used,
                      mskkeys,index,RDMODE);              /* get keys */
     f->build_primary(f->cur_entry,primk,argument_not_used,
                      mskkeys,index,WMODE); 
     record = (file_record *) f->cur_entry;
     record->idtyp = bufrec->idtyp;
     record->addr = WDTO64(write_addr - 1) + 1;
     record->lng = WDTO64(nwords);
   }

   /* update file header */
   f->header->nrec++;
   if (write_to_end) {
      f->header->nxtn++;
      f->header->fsiz += WDTO64(nwords);
      f->nxtadr = W64TOWD(f->header->fsiz) + 1;  /* nxtadr = fsiz +1 */
      f->header->nbig = 
         (WDTO64(nwords) > f->header->nbig) ? WDTO64(nwords) : f->header->nbig;
      if (f->xdf_seq) {  /* add postfix and eof marker */
        postfix.idtyp = 0;
        postfix.lng = 2;
        postfix.addr = -1;
        postfix.prev_idtyp = bufrec->idtyp;
        postfix.prev_lng = bufrec->lng;
        postfix.prev_addr = WDTO64(write_addr -1) + 1;
        c_wawrit(iun,&postfix,f->cur_addr,W64TOWD(2));
        f->cur_addr += W64TOWD(2);
        header64.idtyp = 127;         /* 112 + 15, 15 means EOF */
        header64.lng = 1;
        header64.addr = WDTO64(f->cur_addr -1) + 1;
        f->nxtadr = f->cur_addr;
        c_wawrit(iun,&header64,f->cur_addr,W64TOWD(1));
        FGFDT[index_fnom].file_size = f->nxtadr+W64TOWD(1)-1;
      }
   }
   if (handle != 0)
      f->header->nrwr++;

   f->modified = 1;
   if (! f->xdf_seq)
     f->cur_dir_page->modified = 1;
   return(0);

}

/*splitpoint c_xdfrep */
/***************************************************************************** 
 *                             C _ X D F R E P                               *
 *                                                                           * 
 *Object                                                                     * 
 *   Replace in buf the content of data. nelm elements are replaced in       *
 *   buf starting at position bitpos.                                        *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    nelm     number of elements to replace in buf                      * 
 *  IN    bitpos   bit position of replacement in buf                        *
 *  IN    donnees  replacement data bits                                     *
 *  IN    nbit     number of bits kept per element                           *
 *  IN    datyp    data type                                                 *
 * IN/OUT buf      buffer to contain the modified record                     *
 *                                                                           *
 *****************************************************************************/

int c_xdfrep(word *buffer, word *donnees, int bitpos,
             int nelm, int nbits, int datyp)
{
   int nbwords, index_word, last_ind, i, mode;
   buffer_interface_ptr buf = (buffer_interface_ptr) buffer;
   int ier;

   if ((bitpos % 64) != 0) {
      sprintf(errmsg,"bitpos must be a multiple of 64");
      return(error_msg("c_xdfrep",ERR_BAD_ADDR,ERRFATAL));
      }

   if ((datyp == 3) || (datyp == 5) && (nbits != 8)) {
      sprintf(errmsg,"nbits must be 8 for datyp %d",datyp);
      return(error_msg("c_xdfrep",ERR_BAD_DATYP,ERRFATAL));
      }

   nbwords = (nelm * nbits + 63) / 64;
   nbwords = W64TOWD(nbwords);

   index_word = buf->data_index + (bitpos / (sizeof(word) * 8));

   last_ind = buf->record_index + (buf->nbits / (sizeof(word) *8));

   if ((index_word + nbwords - 1) > buf->nwords) {
      sprintf(errmsg,"buffer not big enough for replacement");
      return(error_msg("c_xdfrep",ERR_BAD_DIM,ERROR));
      }

   /* initialize region to 0 */
   for (i=0; i < nbwords; i++)
      buf->data[index_word+i] = 0;

   /* insert data */
   switch (datyp)
      {
      case 0:         /* transparent mode */
      case 3:

         for (i=0; i < nbwords; i++)
            buf->data[index_word+i] = donnees[i];
         break;
      
#if !defined(NEC)
      case 6:
      case 8:

         for (i=0; i < nbwords; i++)
            buf->data[index_word+i] = donnees[i];
         break;

      case 7:
      case 9:

         if (*little_endian) 
           for (i=0; i < nbwords; i+=2) {
             buf->data[index_word+i] = donnees[i+1];
             buf->data[index_word+i+1] = donnees[i];
           }
         else
           for (i=0; i < nbwords; i++)
             buf->data[index_word+i] = donnees[i];
         break;
#endif

      case 5:        /* upper char only */
         for (i=0; i < nbwords; i++)
            buf->data[index_word+i] = upper_case_word(donnees[i]);
         break;
      
      case 2: mode = 1;
         ier = compact_integer(donnees,(void *) NULL,&(buf->data[index_word]),nelm,
                               nbits,0,xdf_stride,mode);
         break;

      case 4: mode = 3;
         ier = compact_integer(donnees,(void *) NULL,&(buf->data[index_word]),nelm,
                               nbits,0,xdf_stride,mode);
         break;

      default: sprintf(errmsg,"invalid datyp=%d",datyp);
         return(error_msg("c_xdfrep",ERR_BAD_DATYP,ERROR));

      } /* end switch */

   return(0);
}

/*splitpoint c_xdfsta */
/***************************************************************************** 
 *                             C _ X D F S T A                               *
 *                                                                           * 
 *Object                                                                     * 
 *   Get the statistics associated to a file.                                *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    iun      unit number associated to the file                        * 
 *  OUT   stats    statistics of the file                                    *
 *  IN    nstat    number of statistics to print                             *
 *  OUT   pri      primary keys                                              *
 *  IN    npri     number of primary keys                                    *
 *  OUT   aux      auxiliary keys                                            *
 *  IN    naux     number of auxiliary keys                                  *
 *  OUT   vers     software version                                          *
 *  OUT   appl     application signature                                     *
 *                                                                           * 
 *****************************************************************************/

int c_xdfsta(int iun,word *stat,int nstat,
                    word_2 *pri,int npri,word_2 *aux,int naux,
                    char *vers,char *appl)
{
   int index, index_fnom, ier, wasopen=0, i, nn;
   file_header *fh;
   file_record header64;

   index_fnom = fnom_index(iun);
   if (index_fnom == -1) {
      sprintf(errmsg,"file is not connected with fnom");
      return(error_msg("c_xdfsta",ERR_NO_FNOM,WARNING));
      }

   if ((index = file_index(iun)) == ERR_NO_FILE) {
/*
 * open file and read file header
 */
     c_waopen(iun);
     c_waread(iun,&header64,1,W64TOWD(2));
     if (header64.data[0] != 'XDF0' && header64.data[0] !='xdf0') {
     /*if (strncmp(&header64.data[0], "XDF0", 4) != 0 && strncmp(&header64.data[0], "xdf0", 4) != 0) {*/
        sprintf(errmsg,"file is not XDF type\n");
        return(error_msg("c_xdfsta",ERR_NOT_XDF,ERRFATAL));
        }
     if ((fh = calloc(1,header64.lng*8)) == NULL) {
        sprintf(errmsg,"memory is full\n");
        return(error_msg("c_xdfsta",ERR_MEM_FULL,ERRFATAL));
        }
     c_waread(iun,fh,1,W64TOWD(header64.lng));
     }

   else {
/*
 * file is open, file_table[index]->header already contains the required info
 */
      wasopen = 1;
      fh = file_table[index]->header;
      }
   
   switch (nstat)
      { 
      case 12: stat[11] = fh->nrec;
      case 11: stat[10] = fh->neff;
      case 10: stat[9] = fh->laux;
      case  9: stat[8] = fh->naux;
      case  8: stat[7] = fh->lprm;
      case  7: stat[6] = fh->nprm;
      case  6: stat[5] = fh->nbig;
      case  5: stat[4] = fh->plst;
      case  4: stat[3] = fh->nbd;
      case  3: stat[2] = fh->nxtn;
      case  2: stat[1] = fh->nrwr;
      case  1: stat[0] = fh->fsiz;
      case  0: vers[0] = fh->vrsn >> 24 & 0xff;
               vers[1] = fh->vrsn >> 16 & 0xff;
               vers[2] = fh->vrsn >>  8 & 0xff;
               vers[3] = fh->vrsn & 0xff;
               vers[4] = '\0';
               appl[0] = fh->sign >> 24 & 0xff;
               appl[1] = fh->sign >> 16 & 0xff;
               appl[2] = fh->sign >>  8 & 0xff;
               appl[3] = fh->sign & 0xff;
               appl[4] = '\0';
               break;
      default: sprintf(errmsg,"wrong number of stat nstat=%d\n",nstat);
               return(error_msg("c_xdfsta",ERR_BAD_NSTAT,ERROR));
      }

/* 
 * primary keys description
 */
   nn = 0;
   npri = (npri < fh->nprm) ? npri : fh->nprm;
   for (i=0; i < npri; i++,nn++) {
      pri[i].wd1 = fh->keys[nn].ncle;
      pri[i].wd2 = (fh->keys[nn].bit1 << 19) | (fh->keys[nn].lcle << 14) |
                   (fh->keys[nn].tcle << 8) | (fh->keys[nn].reserved);
      }
/*
 * auxiliary keys description
 */
   naux = (naux < fh->naux) ? naux : fh->naux;
   for (i=0; i < naux; i++,nn++) {
      aux[i].wd1 = fh->keys[nn].ncle;
      aux[i].wd2 = (fh->keys[nn].bit1 << 19) | (fh->keys[nn].lcle << 14) |
                   (fh->keys[nn].tcle << 8) | (fh->keys[nn].reserved);
      }
   if (! wasopen)
      c_waclos(iun);
   return(0);
}

/*splitpoint c_xdfunl */
/***************************************************************************** 
 *                             C _ X D F U N L                               * 
 *                                                                           * 
 *Object                                                                     * 
 *   Unlinks the list of random files previously linked by c_xdflnk.         *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    liste    unit number associated to the file                        * 
 *  IN    n        number of files to be unlinked                            *
 *                                                                           *
 *****************************************************************************/

int c_xdfunl(word *liste, int n)
{
  int index, index_fnom, i;
  file_table_entry *f;

  for (i=0; i < n; i++) {
    if ((index_fnom = fnom_index(liste[i])) == -1) {
      sprintf(errmsg,"file is not connected with fnom");
      return(error_msg("c_xdfunl",ERR_NO_FNOM,ERROR));
    }
    if ((index = file_index(liste[i])) == ERR_NO_FILE) {
      sprintf(errmsg,"file is not open");
      return(error_msg("c_xdfunl",ERR_NO_FILE,ERROR));
    }
    f = file_table[index];
    f->link = -1;
    (f->dir_page[f->npages-1])->next_page = NULL;
  }
   
  return(0);
}


/*splitpoint c_xdfupd */
/***************************************************************************** 
 *                          C _ X D F U P D                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Update primary keys and info keys in buffer. If a key value is -1,      *
 *   there is no update of this element.                                     *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    iun      unit number associated to the file                        * 
 * IN/OUT buf      buffer to contain the modified record                     *
 *  IN    idtyp    record type                                               *
 *  IN    keys     list of primary keys                                      *
 *  IN    nkeys    number of primary keys                                    *
 *  IN    info     list of secondary keys                                    *
 *  IN    ninfo    number of secondary keys                                  *
 *                                                                           * 
 *****************************************************************************/

int c_xdfupd(int iun,buffer_interface_ptr buf,int idtyp,
             word *keys,int nkeys,word *info,int ninfo)
{
   int i, index, index_fnom;
   file_record *record;
   file_table_entry *f;
   max_dir_keys argument_not_used;
   word *mskkeys = NULL;

   index_fnom = fnom_index(iun);
   if (index_fnom == -1) {
      sprintf(errmsg,"file is not connected with fnom");
      return(error_msg("c_xdfupd",ERR_NO_FNOM,ERROR));
      }

   if ((index = file_index(iun)) == ERR_NO_FILE) {
      sprintf(errmsg,"file is not open");
      return(error_msg("c_xdfupd",ERR_NO_FILE,ERROR));
      }

   if ((idtyp < 1) && (idtyp != -1) || (idtyp > 126)) {
      sprintf(errmsg,
              "invalid idtyp=%d, must be between 1 and 126 or -1",idtyp);
      return(error_msg("c_xdfupd",ERR_BAD_DATYP,ERROR));
      }
   
   buf->iun = iun;
   
   record = (file_record *) buf->data;
   if (idtyp > -1)
      record->idtyp = idtyp;
   
   f = file_table[index];

   if (nkeys > 0) f->build_primary(buf->data,keys,argument_not_used,
                                   mskkeys,index,WMODE);
   if (ninfo > 0) f->build_info(buf->data + W64TOWD(f->primary_len),
                                info,index,WMODE);

   return(0);
}

/*splitpoint c_xdfuse */
/*****************************************************************************
 *                          C _ X D F U S E                                  *
 *                                                                           *
 *Object                                                                     *
 *   This function creates or updates a destination file containing all the  *
 *   valid records of the source file. If one of the two files is not opened *
 *   it will be open and close at the end except for the destination file    *
 *   which has to be close by the user.                                      *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN    src_unit   unit number associated to the source file               *
 *  IN    dest_unit  unit number associated to the destination file          *
 *                                                                           *
 *****************************************************************************/

int c_xdfuse(int src_unit, int dest_unit)
{
   int i, index_src, index_fnom_src, index_dest, index_fnom_dest;
   int old_xdf_toler = 0, old_msg_level = 0, err, nprim, ninfo;
   int close_src = 0, match;
   file_table_entry *f;
   max_dir_keys argument_not_used;
   static int stat[MAX_STAT];
   static word_2 primk[MAX_KEYS], info[MAX_KEYS];
   word_2 src_primk[MAX_KEYS], src_info[MAX_KEYS];
   int src_stat[MAX_STAT];
   char vers[5],appl[5];

   if (xdf_toler > ERROR) {
      old_xdf_toler = xdf_toler;
      xdf_toler = ERROR;
      }

   if (msg_level > ERROR) {
      old_msg_level = msg_level;
      msg_level = ERROR;
      }

   index_fnom_src = fnom_index(src_unit);
   if (index_fnom_src == -1) {
      sprintf(errmsg,"source file is not connected with fnom");
      return(error_msg("c_xdfuse",ERR_NO_FNOM,ERROR));
      }

   index_fnom_dest = fnom_index(dest_unit);
   if (index_fnom_dest == -1) {
      sprintf(errmsg,"destination file is not connected with fnom");
      return(error_msg("c_xdfuse",ERR_NO_FNOM,ERROR));
      }

/*
 * process destination file
 */
   if ((index_dest = file_index(dest_unit)) == ERR_NO_FILE) {
      if (FGFDT[index_fnom_dest].file_size > 0) {   /* file exist */
         err = c_xdfsta(dest_unit,(word *)&stat,MAX_STAT,primk,MAX_KEYS,
                        info,MAX_KEYS,vers,appl);
         if (err < 0) return(err);
         err = c_xdfopn(dest_unit,"APPEND",primk,stat[6],info,stat[8],appl);
         if (err < 0) return(err);
         }
      else {
         err = c_xdfsta(src_unit,(word *)&stat,MAX_STAT,primk,MAX_KEYS,info,
                        MAX_KEYS,vers,appl);
         if (err < 0) return(err);
         for (i=0; i < xdf_nsplit; i++) {
           err = c_xdfopn(dest_unit+i,"CREATE",primk,stat[6],info,stat[8],appl);
           if (err < 0) return(err);
           }
         }
      }
   else {      /* file is open */
      err = c_xdfsta(dest_unit,(word *)&stat,MAX_STAT,primk,MAX_KEYS,info,
                     MAX_KEYS,vers,appl);
      if (err < 0) return(err);
      }
   nprim = stat[6];
   ninfo = stat[8];

/*
 * process source file
 */
   if ((index_src = file_index(src_unit)) == ERR_NO_FILE) {
      close_src = 1;
      err = c_xdfopn(src_unit,"READ",src_primk,0,src_info,0,appl);
      if (err < 0) return(err);
      }

   err = c_xdfsta(src_unit,(word *)&src_stat,MAX_STAT,src_primk,MAX_KEYS,
                  src_info,MAX_KEYS,vers,appl);
   if (err < 0) return(err);

/*
 * make sure that source file and destination file are the same type
 */
   match = 1;
   for (i=6; i <= 9; i++)
      match = match && (stat[i] == src_stat[i]);
   for (i=0; i < nprim; i++)
      match = match && (primk[i].wd1 == src_primk[i].wd1)
                    && (primk[i].wd2 == src_primk[i].wd2);
   for (i=0; i < ninfo; i++)
      match = match && (info[i].wd1 == src_info[i].wd1)
                    && (info[i].wd2 == src_info[i].wd2);

   if (! match) {
      if (close_src)
         err = c_xdfcls(src_unit);
      sprintf(errmsg,
              "incompatible source (unit %d) and destination (unit %d) files",
              src_unit,dest_unit);
      return(error_msg("c_xdfuse",ERR_NOT_COMP,ERROR));
      }
   if (msg_level <= INFORM)
      err = c_xdfimp(src_unit,(word *)&src_stat,src_stat[6],src_primk,
                     src_info,vers,appl);

/*
 * copy source file to destination
 */
   { xdf_record_header header;
     buffer_interface_ptr buf;
     int readpos=1, nbrec=0, idtyp, lng, addr, nbwords, file_size, nomore;
     int maxmem, initial_mem=5000, last_word;

     if ((file_size = FGFDT[index_fnom_src].file_size) > 0)
        nomore = 0;
     else
        nomore = 1;

     maxmem = WDTO64(initial_mem + RECADDR);

     if ((buf = malloc(maxmem*8)) == NULL) {
        sprintf(errmsg,"malloc can't allocate, no more memory available");
        return(error_msg("c_xdfuse",ERR_MEM_FULL,ERROR));
        }

     buf->nwords = initial_mem + RECADDR;
     buf->data_index = RECADDR;
     buf->iun = src_unit;

     while (! nomore) {
        c_waread(src_unit,&header,readpos,W64TOWD(1));
        addr = (header.addr == 0) ? 1 : header.addr;
        idtyp = header.idtyp;
        lng = header.lng;
        if (W64TOWD(addr-1)+1 == readpos) {     /* valid address */
           if ((idtyp > 0) && (idtyp < 127)) {    /* valid regular record */
              if (lng+100+RECADDR > maxmem) {
                /* Ajout de 100 elements pour couvrir le debut du buffer buf */
                 maxmem = lng+100+RECADDR;
                 if ((buf = realloc(buf,maxmem*8)) == NULL) {
                    sprintf(errmsg,
                            "malloc can't allocate, no more memory available");
                    return(error_msg("c_xdfuse",ERR_MEM_FULL,ERROR));
                    }
                 buf->nwords = W64TOWD(maxmem) + RECADDR;
                 }  /* end if lng > maxmem */
              nbwords = W64TOWD(lng);
              last_word = readpos + nbwords -1;
              if (last_word <= file_size) {
                 buf->nbits = nbwords * sizeof(word) * 8;
                 for (i=0; i < nbwords; i++)
                    buf->data[i] = 0;
                 c_waread(src_unit,&(buf->data),readpos,nbwords);
                 if (xdf_nsplit > 1)
                   err = c_xdfput(dest_unit+(nbrec%xdf_nsplit),0,buf);
                 else
                   err = c_xdfput(dest_unit,0,buf);
                 nbrec++;
                 }
              } /* end valid regular record */
           readpos += W64TOWD(lng);
           if (readpos >= file_size)
              nomore = 1;
           } /* end valid address */
        else
           nomore = 1;
        } /* end while */
     free(buf);
     fprintf(stdout,
             ">>> c_xdfuse: copy of %d records from unit %d to unit %d <<<\n",
             nbrec,src_unit,dest_unit);
     } /* end copy files block */

   if (close_src)
      err = c_xdfcls(src_unit);

   if (old_xdf_toler)
      xdf_toler = old_xdf_toler;

   if (old_msg_level)
      msg_level = old_msg_level;

   return(0);
}

/*splitpoint c_xdfxtr */
/***************************************************************************** 
 *                             C _ X D F X T R                               *
 *                                                                           * 
 *Object                                                                     * 
 *   Extract a portion of the record contained in buf and copy it into       *
 *   vector donnees.                                                         *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    nelm     number of elements to extract into buf                    * 
 *  IN    bitpos   bit position of starting extraction                       *
 *  OUT   donnees  data bits to get                                          *
 *  IN    nbit     number of bits kept per element                           *
 *  IN    datyp    data type                                                 *
 *  IN    buf      buffer to contain the modified record                     *
 *                                                                           *
 *****************************************************************************/

int c_xdfxtr(word *buffer, word *donnees, int bitpos,
             int nelm, int nbits, int datyp)
{
   int nbwords, index_word, last_ind, i, mode;
   buffer_interface_ptr buf = (buffer_interface_ptr) buffer;
   int ier;

   if ((bitpos % 64) != 0) {
      sprintf(errmsg,"bitpos must be a multiple of 64");
      return(error_msg("c_xdfxtr",ERR_BAD_ADDR,ERRFATAL));
      }

   if (((datyp == 3) || (datyp == 5)) && (nbits != 8)) {
      sprintf(errmsg,"nbits must be 8 for datyp %d",datyp);
      return(error_msg("c_xdfxtr",ERR_BAD_DATYP,ERRFATAL));
      }

   nbwords = (nelm * nbits + 63) / 64;
   nbwords = W64TOWD(nbwords);

   index_word = buf->data_index + (bitpos / (sizeof(word) * 8));

   /* extract data */
   switch (datyp)
      {
      case 0:
      case 3:
      case 5:

         for (i=0; i < nbwords; i++)
            donnees[i] = buf->data[index_word+i];
         break;

#if !defined(NEC)
      case 6:
      case 8:
         for (i=0; i < nbwords; i++)
            donnees[i] = buf->data[index_word+i];
         break;

      case 7:
      case 9:
         if (*little_endian) 
           for (i=0; i < nbwords; i+=2) {
             donnees[i+1] = buf->data[index_word+i];
             donnees[i] = buf->data[index_word+i+1];
           }
         else
           for (i=0; i < nbwords; i++)
             donnees[i] = buf->data[index_word+i];
         break;
#endif

      case 2: mode = 2;
         ier = compact_integer(donnees,(void *) NULL,&(buf->data[index_word]),nelm,
                               nbits,0,xdf_stride,mode);
         break;

      case 4: mode = 4;
         ier = compact_integer(donnees,(void *) NULL,&(buf->data[index_word]),nelm,
                               nbits,0,xdf_stride,mode);
         break;

      default: sprintf(errmsg,"invalid datyp=%d",datyp);
         return(error_msg("c_xdfxtr",ERR_BAD_DATYP,ERROR));

      } /* end switch */

   return(0);
}

/*splitpoint create_new_xdf */
/***************************************************************************** 
 *                      C R E A T E _ N E W _ X D F                          *
 *                                                                           * 
 *Object                                                                     * 
 *   Create a new XDF file.                                                  *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    index    file index in table                                       * 
 *  IN    iun      unit number associated to the file                        *
 *  IN    pri      primary keys                                              *
 *  IN    npri     number of primary keys                                    *
 *  IN    aux      auxiliary keys                                            *
 *  IN    naux     number of auxiliary keys                                  *
 *  IN    appl     application signature                                     *
 *                                                                           * 
 *****************************************************************************/

static int create_new_xdf(int index, int iun, word_2 *pri, int npri,
                          word_2 *aux, int naux, char *appl)
{
   file_header *file;
   int ikle=0, lprm=0, laux=0, lng_header=naux+npri+512/64;

   if ((file_table[index]->header = malloc(lng_header*8)) == NULL) {
      sprintf(errmsg,"memory is full\n");
      return(error_msg("create_new_xdf",ERR_MEM_FULL,ERRFATAL));
      }
   file = file_table[index]->header;
   /*   sprintf(file->vrsn,"%s","XDF0"); */
   file->vrsn = 'X' << 24 | 'D' << 16 | 'F' << 8 | '0';
   /*   sprintf(file->sign,"%-4.4s",appl); */
   file->sign = appl[0] << 24 | appl[1] << 16 | appl[2] << 8 | appl[3];
   file->idtyp=0;
   file->lng=lng_header;           /* keys + fixed part of 512 bits */
   file->addr=0;
   file->fsiz=lng_header;          /* keys + fixed part of 512 bits */
   file->nrwr=0;
   file->nxtn=0;
   file->nbd=0;
   file->plst=0;
   file->nbig=0;
   file->nprm=npri;
   file->naux=naux;
   file->neff=0;
   file->nrec=0;
   file->rwflg=0;
   file->reserved=0;

   { int i=0 , bit1=0, lcle=0, tcle=0;
   while (npri--) {
      file->keys[ikle].ncle=pri[i].wd1;
      bit1=pri[i].wd2 >> 19;
      lcle=0x1F & (pri[i].wd2 >> 14);
      tcle=0x3F  & (pri[i].wd2 >> 8);
      file->keys[ikle].bit1=bit1;
      file->keys[ikle].lcle=lcle;
      file->keys[ikle].tcle=tcle;
      file->keys[ikle].reserved=0;
      lprm += lcle;
      i++; ikle++;
      }
   }
   lprm = (lprm+63)/64 + 1 ;        /* primary keys + 64 bit header */
   file->lprm=lprm;
   file_table[index]->primary_len=lprm;
   { int i=0 , bit1=0, lcle=0, tcle=0;
   while (naux--) {
      file->keys[ikle].ncle=aux[i].wd1;
      bit1=aux[i].wd2 >> 19;
      lcle=0x1F & (aux[i].wd2 >> 14);
      tcle=0x3F  & (aux[i].wd2 >> 8);
      file->keys[ikle].bit1=bit1;
      file->keys[ikle].lcle=lcle;
      file->keys[ikle].tcle=tcle;
      file->keys[ikle].reserved=0;
      laux += lcle;
      i++; ikle++;
      }
   }
   laux = (laux+63)/64;
   file->laux=laux;
   file_table[index]->info_len=laux;
   if (! file_table[index]->cur_info->attr.read_only) {
     { int unit=iun, waddress=1 , nwords = file->fsiz*8/sizeof(word);
     c_wawrit(unit,file_table[index]->header,waddress,nwords);
     file_table[index]->nxtadr += nwords;
     }
   }
   return(0);
}

/*splitpoint error_msg */
/***************************************************************************** 
 *                           E R R O R _ M S G                               * 
 *                                                                           * 
 *Object                                                                     * 
 *   Write an error message to stderr and conditionally exit if the level    *
 *   of error is greater than the tolerance level                            *
 *                                                                           *
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    function_name name of the calling function                         * 
 *  IN    errcode       error code                                           *
 *  IN    errlevel      error level                                          *
 *                                                                           *
 *****************************************************************************/

int error_msg(char *function_name, int errcode, int errlevel)
{
   static char *errtab[7] = {"TRIVIAL","INFORMATION","WARNING","ERROR",
                             "FATAL ERROR","SYSTEM SOFTWARE ERROR",
                             "CATASTROPHIC"};

        if (errlevel >= msg_level)
           fprintf(stderr,"*** %s #%d from module %s: %s\n",
                   errtab[errlevel],-errcode,function_name,errmsg);
        if (errlevel > xdf_toler) 
           exit(-errcode);
        return(errcode);
}

/*splitpoint file_index */
/***************************************************************************** 
 *                            F I L E _ I N D E X                            * 
 *                                                                           * 
 *Object                                                                     * 
 *   Find position of file iun in file table.                                *
 *                                                                           *
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    iun   unit number associated to the file                           * 
 *                                                                           *
 *****************************************************************************/

int file_index(int iun)
{
int i;
   for (i=0 ; i<MAX_XDF_FILES ; i++)
      if (file_table[i] != NULL)
         if (file_table[i]->iun == iun)
            return(i);
   return(ERR_NO_FILE);
}

/*splitpoint fnom_index */
/***************************************************************************** 
 *                            F N O M _ I N D E X                            *
 *                                                                           * 
 *Object                                                                     * 
 *   Find index position in master file table (fnom file table).             *
 *                                                                           *
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    iun   unit number associated to the file                           * 
 *                                                                           * 
 *****************************************************************************/

int fnom_index(int iun)
{
   int i;
   for (i=0; i < MAXFILES; i++)           /* find file index in file table */
      if (FGFDT[i].iun == iun) return(i);
   return(-1);
}

/*splitpoint get_free_index */
/***************************************************************************** 
 *                          G E T _ F R E E _ I N D E X                      *
 *                                                                           * 
 *Object                                                                     * 
 *   Find a free position in file table and initialize file attributes.      *
 *                                                                           *
 *                                                                           * 
 *****************************************************************************/

static int get_free_index()
{
   int i, nlimite;
   
   if (STDSEQ_opened == 1) 
     nlimite = 128;
   else
     nlimite = MAX_XDF_FILES;
   for (i=0; i < nlimite; i++) {
      if (file_table[i] == NULL) {
         if ((file_table[i] = (file_table_entry_ptr) calloc(1,sizeof(file_table_entry))) == NULL) {
            sprintf(errmsg,"can't alocate file_table_entry\n");
            return(error_msg("get_free_index",ERR_MEM_FULL,ERRFATAL));
            }
         file_table[i]->file_index = -1;    /* assure first time use of index i */
         init_file(i);
         return(i);
         }
      else {
         if (file_table[i]->iun == -1)
            return(i);
         }
   }
   sprintf(errmsg,"xdf file table is full\n");
   return(error_msg("get_free_index",ERR_FTAB_FULL,ERRFATAL));
}

/*splitpoint init_file */
/***************************************************************************** 
 *                            I N I T _ F I L E                              *
 *Object                                                                     * 
 *  Initialize a file table entry.                                           *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    i   index in the file table                                        * 
 *                                                                           * 
 *****************************************************************************/

static void init_file(int i)
{
   int j;

   for ( j=1 ; j<MAX_DIR_PAGES ; j++ )
      file_table[i]->dir_page[j] = NULL;
   file_table[i]->cur_dir_page = NULL;
   file_table[i]->build_primary = NULL;
   file_table[i]->build_info = NULL;

   file_table[i]->scan_file = NULL;
   file_table[i]->file_filter = NULL;
   file_table[i]->cur_entry = NULL;
   if ((file_table[i]->file_index == i) && (file_table[i]->header != NULL)) {   /* reuse file */
     /* fprintf(stderr,"Debug freeing pointer file_index=%d i=%d\n",file_table[i]->file_index,i); */
     free(file_table[i]->header);
   }
   file_table[i]->header = NULL;
   file_table[i]->nxtadr = 1;
   file_table[i]->primary_len = 0;
   file_table[i]->info_len = 0;
   file_table[i]->link = -1;
   file_table[i]->iun = -1;
   file_table[i]->file_index = i;
   file_table[i]->modified = 0;
   file_table[i]->npages = 0;
   file_table[i]->nrecords = 0;
   file_table[i]->cur_pageno = -1;
   file_table[i]->page_record = 0;
   file_table[i]->page_nrecords = 0;
   file_table[i]->file_version = 0;
   file_table[i]->valid_target = 0;
   file_table[i]->xdf_seq = 0;
   file_table[i]->valid_pos = 0;
   file_table[i]->cur_addr = -1;
   file_table[i]->seq_bof = 1;
   file_table[i]->fstd_vintage_89 = 0;
   for ( j=0 ; j< MAX_SECONDARY_LNG ; j++) {
      file_table[i]->info_keys[j] = 0;
      }
   for ( j=0 ; j< MAX_PRIMARY_LNG ; j++) {
      file_table[i]->head_keys[j] = 0;
      file_table[i]->cur_keys[j] = 0;
      file_table[i]->target[j] = 0;
      file_table[i]->srch_mask[j] = -1;
      file_table[i]->cur_mask[j] = -1;
   }
}


/*splitpoint init_package */
/*****************************************************************************
 *                         I N I T _ P A C K A G E                           *
 *Object                                                                     *
 *  Initialize all file table entries                                        *
 *                                                                           *
 *****************************************************************************/
static void init_package()
{
   int i,ind;

   for (i=0 ; i<MAX_XDF_FILES ; i++)
      file_table[i]=NULL;
   ind = get_free_index(); /* init first entry to start with file index=1 */
   file_table[ind]->iun = 1234567;
   init_package_done = 1;
}

/*splitpoint make_seq_handle */
/***************************************************************************** 
 *                       M A K E _ S E Q _ H A N D L E                       * 
 *                                                                           *
 *Object                                                                     *
 *  calculates an handle for a sequential file from address and index        *
 *                                                                           *
 *Arguments                                                                  *
 *  IN   address     address (in 32bit unit)                                 *
 *  IN   file_index  file index in table                                     *
 *  IN   f           pointer to xdf file information structure               *
 *                                                                           *
 *****************************************************************************/

static INT_32 make_seq_handle(int address, int file_index, file_table_entry *f)
{
  int cluster, addr;
  static int MB512 = 0x4000000, MB128 = 0x1000000, MB32 = 0x400000;
  /* MB512 , MB128 and MB32 are represented in 64 bit unit */
  
  address = (address-1) >> 1;     /* 32 bit word to 64 bit unit address */
  if (f->fstd_vintage_89) {
    cluster = 0;
    address = address / 15;
    addr = address;
  }
  else {

    if (address >= MB512) {
      cluster = 3;
      addr = address >> 6;
    }
    else
      if (address >= MB128) {
        cluster = 2;
        addr = address >> 4;
      }
      else
        if (address >= MB32) {
          cluster = 1;
          addr = address >> 2;
        }
        else {
          cluster = 0;
          addr = address;
        }
  }

  return(MAKE_SEQ_HANDLE(cluster,addr,file_index));
}


/*splitpoint next_match */
/*****************************************************************************
 *                           N E X T _ M A T C H                             *
 *Object                                                                     *
 *  Find the next record that matches the current search criterias.          *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN    file_index   index of file in the file table                       *
 *                                                                           *
 *****************************************************************************/

static word next_match(int file_index)
{
   register file_table_entry *f;
   int found, i, record_in_page, page_no, j, width, match;
   int end_of_file, nw, iun, addr_match;
   /*   xdf_record_header *rec; */
   word *entry, *search, *mask, handle;
   stdf_dir_keys *stds, *stdm, *stde;
   seq_dir_keys *seq_entry;
   xdf_record_header *header;
   ftnword f_datev;
   double nhours;
   int deet,npas,i_nhours,run,datexx;

   /* check if file exists */
   if( (f = file_table[file_index]) == NULL) return (ERR_NO_FILE);

   iun = file_table[file_index]->iun;
   width = W64TOWD(f->primary_len);
   found = 0;
   f->valid_pos = 0;

   if (! f->xdf_seq) {
     /* check if there is a valid current page */
     if(f->cur_dir_page == NULL) return (ERR_NO_POS);
     if(f->cur_entry == NULL) return (ERR_NO_POS);
     if(f->cur_entry - (f->cur_dir_page)->dir.entry != f->page_record *
        W64TOWD(f->primary_len)) return (ERR_NO_POS);

     f->page_nrecords = (f->cur_dir_page)->dir.nent;
     while (! found) {
       if(f->page_record >= f->page_nrecords ){
             /* no more records in page */
         f->page_record = 0;
         f->page_nrecords = 0;
             /* position to next dir page */
         f->cur_dir_page = (f->cur_dir_page)->next_page ;
         if(f->cur_dir_page == NULL) break;  /* no more pages, end of file */
         f->cur_pageno++;
         f->page_nrecords = (f->cur_dir_page)->dir.nent;
         f->cur_entry = (f->cur_dir_page)->dir.entry;
         }
      else{
         for (j = f->page_record; (j < f->page_nrecords) && !found; j++) {
           entry = f->cur_entry;
           header = (xdf_record_header *) entry;
           if (header->idtyp < 127) {
             search = (word *) f->target;
             mask = (word *) f->cur_mask;
             stde = (stdf_dir_keys *) entry;
             stdm = (stdf_dir_keys *) mask;
             stds = (stdf_dir_keys *) search;
             match = 0;
             for (i = 0; i < width; i++, mask++, search++, entry++)
               match |= (((*entry) ^ (*search)) & (*mask));
             found = (match == 0);
             if ( (f->file_filter != NULL) && found ) {   /* no need to call filter if not 'found' */
               handle= MAKE_RND_HANDLE( f->cur_pageno, f->page_record, f->file_index );
               found = found && f->file_filter(handle);
/*               if (found) printf("Debug+ found handle=%d\n",handle);*/
               }
           }
           f->page_record++;  /* position to next record for next search */
           f->cur_entry += width;
         }
      }
     } /* end while */
   }
   else {
     while (! found) {
       nw = c_waread2(iun,f->head_keys,f->cur_addr,width);
       header = (xdf_record_header *) f->head_keys;
       if ((header->idtyp >= 112) || (nw < width)) {
         if ((header->idtyp >= 112) && (header->idtyp < 127)) {
            f->cur_addr += W64TOWD(1);
         }
         end_of_file = 1;
         break;
       }
       if (f->fstd_vintage_89) {   /* old sequential standard */
         if ((stde = calloc(1,sizeof(stdf_dir_keys))) == NULL) {
           sprintf(errmsg,"memory is full");
           return(error_msg("next_match",ERR_MEM_FULL,ERRFATAL));
         }
         seq_entry = (seq_dir_keys *) f->head_keys;
         if (seq_entry->dltf) {
           f->cur_addr += W64TOWD(((seq_entry->lng + 3) >> 2)+15);
           continue;
         }
         if (seq_entry->eof > 0) {
           header->idtyp = 112 + seq_entry->eof;
           header->lng = 1;
           end_of_file = 1;           
           break;
         }
         stde->deleted = 0;
         stde->select = 1;
         stde->lng = ((seq_entry->lng + 3) >> 2) + 15;
         stde->addr = (seq_entry->swa >> 2) +1;
         stde->deet = seq_entry->deet;
         stde->nbits = seq_entry->nbits;
         stde->ni = seq_entry->ni;
         stde->gtyp = seq_entry->grtyp;
         stde->nj = seq_entry->nj;
         stde->datyp = seq_entry->datyp;
         stde->nk = seq_entry->nk;
         stde->ubc = 0;
         stde->npas = (seq_entry->npas2 << 16) | 
           seq_entry->npas1;
         stde->pad7 = 0;
         stde->ig4 = seq_entry->ig4;
         stde->ig2a = 0;
         stde->ig1 = seq_entry->ig1;
         stde->ig2b = seq_entry->ig2 >> 8;
         stde->ig3 = seq_entry->ig3;
         stde->ig2c = seq_entry->ig2 & 0xff;
         stde->etik15 = 
           (ascii6(seq_entry->etiq14 >> 24) << 24) |
           (ascii6((seq_entry->etiq14 >> 16) & 0xff) << 18) |
           (ascii6((seq_entry->etiq14 >>  8) & 0xff) << 12) |
           (ascii6((seq_entry->etiq14      ) & 0xff) <<  6) |
           (ascii6((seq_entry->etiq56 >>  8) & 0xff));
         stde->pad1 = 0;
         stde->etik6a = 
           (ascii6((seq_entry->etiq56      ) & 0xff) << 24) |
           (ascii6((seq_entry->etiq78 >>  8) & 0xff) << 18) |
           (ascii6((seq_entry->etiq78      ) & 0xff) << 12);
         stde->pad2 = 0;
         stde->etikbc = 0;
         stde->typvar = ascii6(seq_entry->typvar) << 6;
         stde->pad3 = 0;
         stde->nomvar = 
           (ascii6((seq_entry->nomvar >>  8) & 0xff) << 18) |
           (ascii6((seq_entry->nomvar      ) & 0xff) << 12);
         stde->pad4 = 0;
         stde->ip1 = seq_entry->ip1;
         stde->levtyp = 0;
         stde->ip2 = seq_entry->ip2;
         stde->pad5 = 0;
         stde->ip3 = seq_entry->ip3;
         stde->pad6 = 0;
         stde->date_stamp = seq_entry->date;
         deet = stde->deet;
         npas = stde->npas;
         if (((deet*npas) % 3600) != 0) {
           /*
            *  recompute datev to take care of rounding used with 1989 version
            *  de-octalise the date_stamp
            */
           run = stde->date_stamp & 0x7;
           datexx = (stde->date_stamp >> 3) * 10 + run;
           
           f_datev = (ftnword) datexx;
           i_nhours = (deet*npas - ((deet*npas+1800)/3600)*3600);
           nhours = (double) (i_nhours / 3600.0);
           f77name(incdatr)(&f_datev,&f_datev,&nhours);
           /*
            *  re-octalise the date_stamp
            */
           datexx = (int) f_datev;
           stde->date_stamp = 8 * (datexx/10) + (datexx % 10);
         }
         
         entry = (word *) stde;
         search = (word *) f->head_keys;
         for (i = 0; i < width; i++, entry++, search++)
           *search = *entry;
         free(stde);
       } /* end if fstd_vintage_89 */
       if ((header->idtyp < 1) || (header->idtyp > 127)) {
         f->cur_addr += W64TOWD(header->lng);
         continue;
       }
       entry = (word *) f->head_keys;
       search = (word *) f->target;
       mask = (word *) f->cur_mask;
       match = 0;
       for (i = 0; i < width; i++, mask++, search++, entry++) {
         match |= (((*entry) ^ (*search)) & (*mask));
       }
       found = (match == 0);

       if (found) {
         f->valid_pos = 1;
         addr_match = f->cur_addr;
       }
       f->cur_addr += W64TOWD(header->lng);   /* position to next record */
       if (! f->fstd_vintage_89)
         f->cur_addr += W64TOWD(2);           /* skip postfix */
     }
   }

   if (! found) return(ERR_NOT_FOUND);

   if (! f->xdf_seq) {
     if (msg_level <= TRIVIAL) 
       fprintf(stdout,"Record found at page# %d, record# %d\n",
               f->cur_pageno,f->page_record-1);
     handle= MAKE_RND_HANDLE( f->cur_pageno, f->page_record-1, f->file_index );
   }
   else {
     if (msg_level <= TRIVIAL)
       fprintf(stdout,"Record found at address %d,\n",addr_match);
     stde = (stdf_dir_keys *) f->head_keys; 
     handle = make_seq_handle(addr_match,f->file_index,f);
   }
   return(handle);
}
#include "if_xdf98.hc"
