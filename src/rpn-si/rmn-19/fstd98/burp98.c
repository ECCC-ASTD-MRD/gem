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
#ifdef WIN32    /*CHC/NRC*/
#include <string.h>
#else
#include <unistd.h>
#endif
#include <stdlib.h>
#include "qstdir.h"
#include "proto.h"
int BurP_nele;
int BurP_ntot;
/*splitpoint c_buf89a0 */
/***************************************************************************** 
 *                            C _ B U F 8 9 A 0                              * 
 *****************************************************************************/

void c_buf89a0(word *buffer)
{
  buffer_interface_ptr buf = (buffer_interface_ptr) buffer;

  buf->buf78.buf8 = 0;
  buf->buf9 = 0;
}
/*splitpoint c_getbuf8 */
/***************************************************************************** 
 *                            C _ G E T B U F 8                              * 
 *****************************************************************************/

int c_getbuf8(word *buffer)
{
  buffer_interface_ptr buf = (buffer_interface_ptr) buffer;

  return(buf->buf78.buf8);
}

/*splitpoint build_burp_info_keys */
/***************************************************************************** 
 *                  B U I L D _ B U R P _ I N F O _ K E Y S                  * 
 *                                                                           *
 *Object                                                                     *
 *  Pack burp info keys into buffer or get info keys from buffer depending   *
 *  on mode argument.                                                        *
 *                                                                           *
 *Arguments                                                                  *
 *  IN/OUT  buf     buffer to contain the keys                               *
 *  IN/OUT  keys    info keys                                                *
 *  IN      index   file index in file table                                 *
 *  IN      mode    if mode = WMODE,write to buffer otherwise get keys from  *
 *                  buffer.                                                  *
 *                                                                           *
 *****************************************************************************/

void build_burp_info_keys(word *buf, word *keys, int index, int mode)
{
   burp_dir_info *info;
   
   info = (burp_dir_info *) buf;

   if (mode == WMODE) {     /* write info keys to burp record */
      if (keys[0] != -1) info->nblks = keys[0];
      if (keys[1] != -1) info->oars  = keys[1];
      if (keys[2] != -1) info->elev  = keys[2];
      if (keys[3] != -1) info->drcv  = keys[3];
      if (keys[4] != -1) info->runn  = keys[4];
      }
   else {                   /* read info keys from burp record */
      keys[0] = info->nblks;
      keys[1] = info->oars;
      keys[2] = info->elev;
      keys[3] = info->drcv;
      keys[4] = info->runn;
      }
}

/*splitpoint build_burp_prim_keys */
/***************************************************************************** 
 *                  B U I L D _ B U R P _ P R I M _ K E Y S                  *
 *                                                                           *
 *Object                                                                     *
 *  Pack burp primary keys into buffer or get primary keys from buffer       *
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

void build_burp_prim_keys(burp_record *brpk, word *keys,
                                 burp_record *mask, word *mskkeys,
				 int index, int mode)
{
   int *mskptr;
   int i;

   mskptr = (int *) mask;

   if (mode == WMODE) {     /* write keys to burp record */
/*
 * initialize compress keys to 0 and bit masks to all ones
 */

      for (i=0; i < sizeof(burp_dir_keys) / sizeof(int); i++) {
	 *mskptr = -1;
	 mskptr++;
	 }

      mask->keys.idtyp = 0;
      mask->keys.lng = 0;
      mask->keys.addr = 0;

      if (keys[0] == -1)  mask->keys.sti1 = 0;
      else
	 brpk->keys.sti1 = upper_case(keys[0]);
      if (keys[1] == -1)  mask->keys.sti2 = 0;
      else
	 brpk->keys.sti2 = upper_case(keys[1]);
      if (keys[2] == -1)  mask->keys.sti3 = 0;
      else
	 brpk->keys.sti3 = upper_case(keys[2]);
      if (keys[3] == -1)  mask->keys.sti4 = 0;
      else
	 brpk->keys.sti4 = upper_case(keys[3]);
      if (keys[4] == -1)  mask->keys.sti5 = 0;
      else
	 brpk->keys.sti5 = upper_case(keys[4]);
      if (keys[5] == -1)  mask->keys.sti6 = 0;
      else
	 brpk->keys.sti6 = upper_case(keys[5]);
      if (keys[6] == -1)  mask->keys.sti7 = 0;
      else
	 brpk->keys.sti7 = upper_case(keys[6]);
      if (keys[7] == -1)  mask->keys.sti8 = 0;
      else
	 brpk->keys.sti8 = upper_case(keys[7]);
      if (keys[8] == -1)  mask->keys.sti9 = 0;
      else
	 brpk->keys.sti9 = upper_case(keys[8]);
      if (keys[9] == -1)  mask->keys.flgs = 0;
      else
	 brpk->keys.flgs = keys[9];
      if (keys[10] == -1) mask->keys.lati = 0;
      else
	 brpk->keys.lati = keys[10];
      if (keys[11] == -1) mask->keys.lon  = 0;
      else
	 brpk->keys.lon  = keys[11];
      if (keys[12] == -1) mask->keys.date = 0;
      else
	 brpk->keys.date = keys[12];
      if (keys[13] == -1) mask->keys.dx   = 0;
      else
	 brpk->keys.dx   = keys[13];
      if (keys[14] == -1) mask->keys.idtp = 0;
      else
	 brpk->keys.idtp = keys[14];
      if (keys[15] == -1) mask->keys.dy   = 0;
      else 
	 brpk->keys.dy   = keys[15];
      if (keys[16] == -1) mask->keys.heur = 0;
      else
	 brpk->keys.heur = keys[16];
      if (keys[17] == -1) mask->keys.min  = 0;
      else
	 brpk->keys.min  = keys[17];
      }

   else {             /* read keys from burp record */ 
   
      keys[0] = brpk->keys.sti1;
      keys[1] = brpk->keys.sti2;
      keys[2] = brpk->keys.sti3;
      keys[3] = brpk->keys.sti4;
      keys[4] = brpk->keys.sti5;
      keys[5] = brpk->keys.sti6;
      keys[6] = brpk->keys.sti7;
      keys[7] = brpk->keys.sti8;
      keys[8] = brpk->keys.sti9;
      keys[9] = brpk->keys.flgs;
      keys[10] = brpk->keys.lati;
      keys[11] = brpk->keys.lon;
      keys[12] = brpk->keys.date;
      keys[13] = brpk->keys.dx;
      keys[14] = brpk->keys.idtp;
      keys[15] = brpk->keys.dy;
      keys[16] = brpk->keys.heur;
      keys[17] = brpk->keys.min;
      }
}

/**----------------------------------------------------------------------------
 * @brief  Check burp file for corruption
 * @date   October 2020
 *    @param[in]  Filename Path to the file
 *
 *    @return     Valid(0) or not(-1)
*/
int c_burpcheck(const char* filename)
{
   return(c_xdfcheck(filename));
}
/*splitpoint burp_nbit_datyp */
/***************************************************************************** 
 *                  B U R P _ N B I T _ D A T Y P                            *
 *                                                                           *
 *Object                                                                     *
 *  Rewritten version of original fortran QRBNBDT by James Caveen.           *
 *                                                                           *
 *  Determine the number of bits and the proper datyp needed to insert       *
 *  tblval into a data block.                                                *
 *                                                                           *
 *Arguments                                                                  *
 *  IN/OUT  nbit    number of bits                                           *
 *  IN/OUT  datyp   data type of tblval                                      *
 *  IN      tblval  array of values to be inserted subsequently              *
 *  IN      tbldim  tblval dimension                                         *
 *  IN      stride  stride to use within tblval                              *
 *                                                                           *
 *****************************************************************************/
static int burp_nbit_datyp(int *nbits, int *datyp, word *tval,int tbldim,
			   int stride)
{
  int i, needed, inbits=*nbits;
  INT_32 *tblval = (INT_32 *) tval, temp, tblmax, tblmin;
  int v_nbits = inbits, v_datyp = *datyp;

  /* if transparent mode, nomodifications */
  if (((*datyp==2) && (*nbits == 32)) || (datyp == 0))
    return(0);

#if defined(NEC)
      /* check if datyp is greater than 5 --> not allowed on the NEC */
  if (*datyp > 5) { 
    sprintf(errmsg,"datyp 6,7,8 or 9 is not allowed on the NEC");
    return(error_msg("burp_nbit_datyp",ERR_BAD_DATYP,ERRFATAL));
  }
#endif

  /* if datyp > 5 (real, real*8, complex, complex*8) nbits is set to 32 */
  if (*datyp > 5) {
    *nbits = 32;
    return(0);
  }
  
  /* if datyp = character type, nbits is set to 8 */
  if ((*datyp == 3) || (*datyp == 5)) {
    *nbits = 8;
    return(0);
  }
  
  /* find min and max */
  tblmin = tblval[0]; tblmax = tblval[0];
  for (i=0; i < tbldim*stride; i+=stride) {
    tblmax = Max(tblval[i],tblmax);
    tblmin = Min(tblval[i],tblmin);
  }

  if ((tblmin == tblmax) && (tblmin == -1)) return(0);

  /* check if signed integer is required and set max to ABS value */
  if (tblmin < -1) {
    tblmax = (tblmax > (-tblmin)) ? tblmax : (-tblmin);
    *datyp = 4;
  }
  
  /* determine the number of bits needed to represent the maximum value. */
  /* If there are missing values (-1) one additional bit is needed */
  temp = tblmax;
  needed = 0;
  while ((temp) && (needed < 32)) {
    needed++;
    temp >>=1;
  }
  if (tblmax == (~(-1 << needed))) needed++;
  
  *nbits = needed;
  if (*datyp == 4) {
    (*nbits)++;
    if (*nbits > 31) {
      *nbits = 32;
      *datyp = 2;
      sprintf(errmsg,"encoding values < 0 with nbit=32 and datyp=2");
      return(error_msg("burp_nbit_datyp",BURP_ERR_CMPR,WARNING));
    }
  }
  *nbits = (*nbits > 32) ? 32 : *nbits;
  *nbits = (inbits > *nbits) ? inbits : *nbits;

  return(0);
}

/*splitpoint burp_valid789 */
/***************************************************************************** 
 *                          B U R P _ V A L I D 7 8 9                        *
 *                                                                           *
 *Object                                                                     *
 *  Rewritten version of original fortran QRV7A9 by James Caveen.            *
 *     objet (qrv7a9) 
 *     Traitement de la liste d'elements pour les data de type 
 *     7 = real*8
 *     8 = complex
 *     9 = complex*8
 *     
 *     Selon le datyp, on s'assure que le nombre d'elements est bon
 *     et que la valeur des elements est adequate.  Si la valeur d'un
 *     element est zero, on lui attribut la bonne valeur.  Les elements 
 *     verifies sont les suivants:
 *     
 *     055204  32 BITS INFERIEURS D'UN REEL*8 OCTETS      
 *     055205  IMAGINAIRE D'UN COMPLEXE*4 OCTETS          
 *     055206  32 BITS SUPERIEURS, IMAGINAIRE COMPLEXE*8  
 *     055207  32 BITS INFERIEURS, IMAGINAIRE COMPLEXE*8  
 *     
 *     Exemple:  temperature du point de rose en reel*8  (datyp = 7)
 *     code burp 012003  = temperature du point de rose
 *     055204  = 32 bits du bas d'un reel*8
 *     
 *     Si la liste d'elements contient 2 elements soit:
 *     elem1     elem2
 *     012003    055204           --> O.K.
 *     012003         0           --> remplace 0 par 055204
 *     012003    !=055204 | !=0   --> erreur 
 *     
 *     Si mod(nele,2) != 0        --> erreur
 *     
 *                                                                           *
 *Arguments                                                                  *
 *  IN/OUT  lstele  list of elements                                         *
 *  IN      nele    number of elements in list                               *
 *  IN      datyp   data type                                                *
 *                                                                           *
 *****************************************************************************/
static int burp_valid789(word *lstele, int nele, int datyp)
{

  int i, j, codval[3];
  
  switch (datyp) {

  case 7:
  case 8:
    if ((nele &1) != 0) {
      sprintf(errmsg,"datyp=%d, nele must be even, nele=%d",datyp,nele);
      return(error_msg("burp_valid789",ERR_BAD_DATYP,ERRFATAL));
    }
    if (datyp == 7)
      codval[0] = MRBCOV(0,55,204);
    else
      codval[0] = MRBCOV(0,55,205);
    for (i=1; i < nele; i+=2) {
      if (lstele[i] != codval[0]) {
	if (lstele[i] == 0)
	  lstele[i] = codval[0];
	else {
	  sprintf(errmsg,"invalid code for datyp %d",datyp);
	  return(error_msg("burp_valid789",BURP_ERR_CODE,ERRFATAL));
	}
      }
    }
    break;

  case 9:
    if ((nele & 3) != 0) {
      sprintf(errmsg,"datyp=%d, nele must be a multiple of 4, nele=%d",
              datyp,nele);
      return(error_msg("burp_valid789",ERR_BAD_DATYP,ERRFATAL));
    }
    codval[0] = MRBCOV(0,55,204);
    codval[1] = MRBCOV(0,55,206);
    codval[2] = MRBCOV(0,55,207);
    for (i=1; i < nele; i+=4) {
      for (j=0; j < 3; j++) {
        if (lstele[i+j] != codval[j]) {
          if (lstele[i+j] == 0)
            lstele[i+j] = codval[j];
          else {
            sprintf(errmsg,"invalid code for datyp %d",datyp);
            return(error_msg("burp_valid789",BURP_ERR_CODE,ERRFATAL));
          }
        }
      }
    }
    break;
  }
  return(0);
}
    

/*splitpoint c_mrbadd */
/***************************************************************************** 
 *                             C _ M R B A D D                               *
 *                                                                           * 
 *Object                                                                     * 
 *   Add a data block at the end of the report.                              *
 *                                                                           *
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN/OUT buffer vector to contain the report                               * 
 *    OUT  bkno   number of blocks in buf                                    *
 *    IN   nele   number of meteorogical elements in block                   *
 *    IN   nval   number of data per elements                                *
 *    IN   nt     number of group of nele*nval values in block               *
 *    IN   bfam   block family (12 bits, bdesc no more used)                 *
 *    IN   bdesc  kept for backward compatibility                            *
 *    IN   btyp   block type                                                 *
 *    IN   nbit   number of bit to keep per values                           *
 *    OUT  bit0   position of first bit of the report                        *
 *    IN   datyp  data type for packing                                      *
 *    IN   lstele list of nele meteorogical elements                         *
 *    IN   tblval array of values to write (nele*nval*nt)                    *
 *                                                                           * 
 *****************************************************************************/

int c_mrbadd(void *buffer, int *bkno, int nele, int nval, int nt, int bfam,
	     int bdesc, int btyp, int nbit, int *bit0, int datyp,
	     word *lstele, word *tblval)
{
   burp_block_header entete;
   buffer_interface_ptr buf = (buffer_interface_ptr) buffer;
   word *pos;
   int inbit,idatyp,nombre,err,temp,left,i,bits_added,bfamho;
   int done,indx;
   word r_nele,r_nval,r_nt,r_bfam,r_bdesc;
   word r_btyp,r_nbit,r_bit0,r_datyp;

   /* initialize block header to 0 */
   pos = (word *) &entete;
   for (i=0; i < (sizeof(entete) / sizeof(word)); i++)
      *pos = 0;

#if defined(NEC)
   if (datyp > 5) {
      sprintf(errmsg,"datyp 6,7,8 or 9 not allowed on the NEC");
      return(error_msg("c_mrbadd",ERR_BAD_DATYP,ERRFATAL));
      }
#endif
   
   if (((datyp == 3) || (datyp == 5)) && (nbit != 8)) {
      sprintf(errmsg,"nbits must be 8 for datyp 3 or 5");
      return(error_msg("c_mrbadd",ERR_BAD_DATYP,ERROR));
      }
   
   /* validate (possibly modify) nbit and datyp */
   inbit = nbit;
   idatyp = datyp;
   nombre = nele*nval*nt;
   err = burp_nbit_datyp(&inbit,&idatyp,tblval,nombre,xdf_stride);
   if (err < 0) return(err);

   /* validates lstele for datyp 7,8 or 9 */
   if (datyp >= 7) {
      err = burp_valid789(lstele,nele,datyp);
      if (err < 0) return(err);
      }
      
   if (bdesc != 0) {
     bfamho = (bfam >> 6) & 0x3f;
     if ((bfamho != 0) && (bfamho != bdesc)) {
        sprintf(errmsg,"illegal use of bdesc");
	return(error_msg("c_mrbadd",BURP_ERR_BDESC,ERRFATAL));
	}
     entete.bfamdesc = (bfam & 0x3f) << 6;
     entete.bfamdesc |= (bdesc & 0x3f);
     }
   else {
      entete.bfamdesc = (bfam & 0x3f) << 6;
      entete.bfamdesc |= ((bfam >> 6) & 0x3f);
      }
   
   entete.btyp = btyp;
   entete.nbit = inbit-1;
   entete.datyp = idatyp;
   *bit0 = (buf->nbits - NBENTR - buf->buf9) / 64;
   *bit0 = (*bit0 < 0) ? 0 : *bit0;
   entete.bit0 = *bit0;
   if ((nele >= GROSNELE) || (nval >= GROSDIM) || (nt >= GROSDIM)) {
      entete.flag = 1;
      entete.nele = GROSNELE;
      entete.elem1 = nele;
      entete.elem2 = nval;
      entete.elem3 = nt;
      done = 0;
      }
   else {
      entete.flag = 0;
      entete.nt = nt;
      entete.nele = nele;
      entete.nval = nval;
      entete.elem1 = lstele[0];
      entete.elem2 = lstele[1];
      entete.elem3 = lstele[2];
      done = 3;
      }
   indx = (buf->nbits) / (8 * sizeof(word));
   pos = (word *) &(buf->data[indx]);
   *pos = 0;
   temp = 0;
   left = 8 * sizeof(word);
   bits_added = 0;
/*
   for (i=done; i < nele; i++) {
      if (left == 0) {
         *pos = temp;
	 pos ++;
	 temp = 0;
	 left = 8 * sizeof(word);
	 bits_added += left;
	 }
      temp = (temp << 16) | (lstele[i] & 0xffff);
      left -= 16;
      }
   if (temp != 0) {
      temp <<= left;
      *pos = temp;
      bits_added += 8*sizeof(word);
      }
*/
   i = done;
   while (i < nele) {
     *pos = (*pos << (32-left)) | (lstele[i] & 0xffff);
     bits_added += 16;
     left -= 16;
     if (left == 0) {
       pos++;
       *pos = 0;
       left = 8*sizeof(word);
     }
     i++;
   }
   if (left != 8*sizeof(word)) {
     *pos <<= 16;
     bits_added += 16;
   }
  
   buf->nbits += (bits_added +63) / 64 * 64;
   err = c_xdfins((word *)buf,(word *)&entete,buf->buf9,DIMENT,
		  8*sizeof(word),0);
   if (err < 0) return(err);
   
   err = c_xdfadd((word *)buf,tblval,nombre,inbit,idatyp);
   if (err < 0) return(err);
   
   buf->buf78.buf8++;
   buf->buf9 += NBENTB;
   *bkno = buf->buf78.buf8;

   if (msg_level <= INFORM) {
      err = c_mrbprm((word *)buf,*bkno,&r_nele,&r_nval,&r_nt,&r_bfam,&r_bdesc,
		     &r_btyp,&r_nbit,&r_bit0,&r_datyp);
      fprintf(stdout,"MRBADD - write block #%5d NELE=%5d NVAL=%5d NT=%5d BFAM=%4d BTYP=%4d NBITS=%2d BIT0=%8d DATYP=%1d\n",
	      *bkno,r_nele,r_nval,r_nt,r_bfam,r_btyp,r_nbit,r_bit0,r_datyp);
      }
   return(0);
}

/*splitpoint c_mrbdel */
/***************************************************************************** 
 *                             C _ M R B D E L                               *
 *                                                                           * 
 *Object                                                                     * 
 *   Delete a particular block of the report.                                *
 *                                                                           *
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN/OUT buffer   vector to contain the report                             * 
 *    IN   number   block number to be deleted                               *
 *                                                                           * 
 *****************************************************************************/

int c_mrbdel(void *buffer, int number)
{
   buffer_interface_ptr buf = ( buffer_interface_ptr) buffer;
   int err, bitpos, bit0, prebit, difbit, nele, nval, nt, nbit;
   int nelements, datyp, nombre, i;
   burp_block_header entete, *block;
   burp_record *burprec;
   word *pos;

   if ((number < 1) || (number > buf->buf78.buf8)) {
      sprintf(errmsg,"invalid block number");
      return(error_msg("c_mrbdel",BURP_ERR_BNUM,ERROR));
      }

   bitpos = NBENTB * (number -1);
   err = c_xdfxtr((word *)buf,(word *)&entete,bitpos,DIMENT,8*sizeof(word),0);
   if (err < 0) return(err);

   if (entete.flag != 0) {
      nele = entete.elem1;
      nval = entete.elem2;
      nt = entete.elem3;
      nelements = nele;
      }
   else {
      nele = entete.nele;
      nval = entete.nval;
      nt = entete.nt;
      nelements = ((nele-3) > 0) ? nele-3 : 0;
      }
   nbit = entete.nbit + 1;
   datyp = entete.datyp;
   bit0 = entete.bit0;

   prebit = bit0;
   bit0 = bit0 * 64 + buf->buf9;

/* cut the meteorogical element list */
   if (nelements > 0) err = c_xdfcut(buf,bit0,nelements,16,2);

/* cut the value array */
   nombre = nele*nval*nt;
   err = c_xdfcut(buf,bit0,nombre,nbit,datyp);

/* cut block header */
   err = c_xdfcut(buf,bitpos,DIMENT,8*sizeof(word),0);

/* update number of blocks */
   buf->buf78.buf8--;
   
/* update bit position of block */
   buf->buf9 -= NBENTB;

/* update bit0 for each block */
   burprec = (burp_record *) buf->data;
   block = (burp_block_header *) burprec->data;
   bitpos = block[number-1].bit0;
   difbit = bitpos - prebit;

   for (i=number-1; i < buf->buf78.buf8; i++)
      block[i].bit0 -= difbit;

   return(0);
}

/*splitpoint c_mrbhdr */
/***************************************************************************** 
 *                             C _ M R B H D R                               *
 *                                                                           * 
 *Object                                                                     * 
 *  Return the description parameters of the data block of order bkno.       *
 *                                                                           *
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    buf    vector to contain the report                                * 
 *  IN    bkno   block number                                                *
 *  OUT   nele   number of elements                                          *
 *  OUT   nval   number of values per element                                *
 *  OUT   nt     number of NELE*NVAL values                                  *
 *  OUT   bfam   block family type (12 bits)                                 *
 *  OUT   bdesc  block descriptor (set to zero)                              *
 *  OUT   btyp   block type                                                  *
 *  OUT   nbit   number of bits kept per value                               *
 *  OUT   bit0   first bit of array values                                   *
 *  OUT   datyp  data compaction type                                        *
 *                                                                           *
 *****************************************************************************/
int c_mrbhdr(word *buf, int *temps, int *flgs, char *stnid, int *idtyp,
             int *lati, int *lon, int *dx, int *dy, int *elev,
	     int *drcv, int *date, int *oars, int *run, int *nblk,
	     word *sup, int nsup, word *xaux, int nxaux)
{

  buffer_interface_ptr buffer = (buffer_interface_ptr) buf;
  burp_record *burprec;
  int AA,MM,JJ,annee,mois;

  if (nsup > NPRISUP) {
    sprintf(errmsg,"there is too many supplementary prim keys");
    error_msg("c_mrbhdr",BURP_ERR_CLEF,WARNING);
  }
  if (nxaux > NAUXSUP) {
    sprintf(errmsg,"there is too many supplementary aux keys");
    error_msg("c_mrbhdr",BURP_ERR_CLEF,WARNING);
  }
  
  burprec = (burp_record *) buffer->data;
  stnid[0] = burprec->keys.sti1;
  stnid[1] = burprec->keys.sti2;
  stnid[2] = burprec->keys.sti3;
  stnid[3] = burprec->keys.sti4;
  stnid[4] = burprec->keys.sti5;
  stnid[5] = burprec->keys.sti6;
  stnid[6] = burprec->keys.sti7;
  stnid[7] = burprec->keys.sti8;
  stnid[8] = burprec->keys.sti9;
  *temps = burprec->keys.heur * 100 + burprec->keys.min;
  *flgs = burprec->keys.flgs;
  *idtyp = burprec->keys.idtp;
  *lati = burprec->keys.lati;
  *lon = burprec->keys.lon;
  *dx = burprec->keys.dx;
  *dy = burprec->keys.dy;
  *elev = burprec->info.elev;
  *drcv = burprec->info.drcv;
  *oars = burprec->info.oars;
  *run = burprec->info.runn;
  *nblk = buffer->buf78.buf8;
  *date = burprec->keys.date;
  if ((((*date/100) % 100) > 12) || (xdf_enforc8)) {
/*
 *     retourner la date en format AAAAMMJJ
 */
    AA = (*date/10000) % 100;
    MM = (*date/100) % 100;
    JJ = *date % 100;
    annee = 1900 + AA + (((MM-1)/12)*100);
    mois = 1 + ((MM-1) % 12);
    *date = (annee * 10000) + (mois * 100) + JJ;
  }

  /* CLES SUPPLEMENTAIRES */
  return(0);
}

/*splitpoint c_mrblen */
/***************************************************************************** 
 *                             C _ M R B L E N                               *
 *                                                                           * 
 *Object                                                                     * 
 *   Return the number of bits used in buf and the number of bits left.      *
 *                                                                           *
 *Arguments                                                                  * 
 *                                                                           * 
 *    IN   buffer   vector that contains the report                          * 
 *   OUT   lbits    number of bits used                                      *
 *   OUT   left     number of bits left                                      *
 *                                                                           * 
 *****************************************************************************/

int c_mrblen(void *buffer, int *lbits, int *left)
{

   buffer_interface_ptr buf = (buffer_interface_ptr) buffer;

   *lbits = buf->nbits;
   *left = ((buf->nwords - 9) * 8*sizeof(word)) - *lbits;
   return(0);
}

/*splitpoint c_mrbloc */
/***************************************************************************** 
 *                             C _ M R B L O C                               *
 *                                                                           * 
 *Object                                                                     * 
 *   Search for a specific block in the buffer. Search starts at block       *
 *   blkno. If blkno = 0 search starts from the beginning.                   *
 *                                                                           *
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN   buffer vector to contain the report                                 * 
 *  IN   bfam   block family (12 bits, bdesc no more used)                   *
 *  IN   bdesc  kept for backward compatibility                              *
 *  IN   btyp   block type                                                   *
 *  in   bkno   number of blocks in buf                                      *
 *                                                                           * 
 *****************************************************************************/

int c_mrbloc(void *buffer, int bfam, int bdesc, int btyp, int blkno)
{
   buffer_interface_ptr buf = (buffer_interface_ptr) buffer;
   int i,bfamho,bfamdesc,err;
   int match, mskbtyp, mskfamdesc;
   int r_nele,r_nval,r_nt,r_bfam,r_bdesc,r_btyp,r_nbit,r_bit0,r_datyp;
   int bno;
   burp_block_header *block;
   burp_record *burprec;


   if (bfam == -1) {
      if ((bdesc == 0) || (bdesc == -1))
	 bfamdesc = -1;
      else {
         sprintf(errmsg,"illegal use of bdesc");
	 return(error_msg("c_mrbloc",BURP_ERR_BDESC,ERRFATAL));
	 }
      }
   else if ((bdesc != 0) && (bdesc != -1)) {
     bfamho = (bfam >> 6) & 0x3f;
     if ((bfamho != 0) && (bfamho != bdesc)) {
        sprintf(errmsg,"illegal use of bdesc");
	return(error_msg("c_mrbloc",BURP_ERR_BDESC,ERRFATAL));
	}
     bfamdesc = (bfam & 0x3f) << 6;
     bfamdesc |= (bdesc & 0x3f);
     }
   else {
      bfamdesc = (bfam & 0x3f) << 6;
      bfamdesc |= ((bfam >> 6) & 0x3f);
      }

   if (bfamdesc == -1)
      mskfamdesc = 0;
   else
      mskfamdesc = -1;

   if (btyp == -1)
      mskbtyp = 0;
   else {
      mskbtyp = 0;
      mskbtyp |= (0x10000000 & btyp) ? 0 : 0x000f;
      mskbtyp |= (0x20000000 & btyp) ? 0 : 0x07f0;
      mskbtyp |= (0x40000000 & btyp) ? 0 : 0x7800;
      }

   burprec = (burp_record *) buf->data;
   block = (burp_block_header *) burprec->data;
   match = 0;
   for (i=blkno; i < buf->buf78.buf8; i++) {
      match = ((((bfamdesc ^ block[i].bfamdesc) & mskfamdesc) == 0) &&
	       (((btyp ^ block[i].btyp) & mskbtyp) == 0));
      if (match) {
         if (msg_level <= INFORM) {
	    bno = i+1;
            err = c_mrbprm(buf,bno,&r_nele,&r_nval,&r_nt,
				  &r_bfam,&r_bdesc,&r_btyp,&r_nbit,
				  &r_bit0,&r_datyp);
	    fprintf(stdout,"MRBLOC - find block #%5d NELE=%5d NVAL=%5d NT=%5d BFAM=%4d BTYP=%4d NBITS=%2d BIT0=%8d DATYP=%1d\n",
		    i+1,r_nele,r_nval,r_nt,r_bfam,r_btyp,r_nbit,r_bit0,r_datyp);
	    }
	 return(i+1);
	 } /* end match */
      } /* end for */
   if (msg_level <= INFORM) {
      fprintf(stdout,"MRBLOC - block not found bfam=%d, bdesc=%d, btyp=%d\n",
	      bfam,bdesc,btyp);
      }
   return(-1);
}

/*splitpoint c_mrbprm */
/***************************************************************************** 
 *                             C _ M R B P R M                               *
 *                                                                           * 
 *Object                                                                     * 
 *  Return the description parameters of the data block of order bkno.       *
 *                                                                           *
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN    buf    vector to contain the report                                * 
 *  IN    bkno   block number                                                *
 *  OUT   nele   number of elements                                          *
 *  OUT   nval   number of values per element                                *
 *  OUT   nt     number of NELE*NVAL values                                  *
 *  OUT   bfam   block family type (12 bits)                                 *
 *  OUT   bdesc  block descriptor (set to zero)                              *
 *  OUT   btyp   block type                                                  *
 *  OUT   nbit   number of bits kept per value                               *
 *  OUT   bit0   first bit of array values                                   *
 *  OUT   datyp  data compaction type                                        *
 *                                                                           *
 *****************************************************************************/
int c_mrbprm(word *buf,int  bkno, int *nele, int *nval, int *nt, int *bfam,
             int *bdesc, int *btyp, int *nbit, int *bit0, int *datyp)
{
  burp_block_header header;
  int ier, bitpos;

  bitpos = (bkno-1) * NBENTB;
  ier = c_xdfxtr(buf,(word *)&header,bitpos,DIMENT,32,0);
  
  *btyp = header.btyp;
  *nbit = header.nbit +1;
  *bit0 = header.bit0;
  *datyp = header.datyp;
  *nele = header.nele;

  if (header.flag) {
    *nele = header.elem1;
    *nval = header.elem2;
    *nt = header.elem3;
  }
  else {
    *nval = header.nval;
    *nt = header.nt;
  }

  *bfam = (header.bfamdesc & 0x3f) << 6;
  *bfam |= ((header.bfamdesc >> 6) & 0x3f);
  *bdesc = 0;
  return(ier);
}

/*splitpoint c_mrbrep */
/***************************************************************************** 
 *                             C _ M R B R E P                               *
 *                                                                           * 
 *Object                                                                     * 
 *   Replace a data block by an other one with the same variables and        *
 *   dimensions.                                                             *
 *                                                                           *
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN/OUT buffer vector that contains the report                            * 
 *    IN   bkno   block number to be replaced                                *
 *    IN   tblval array of values to write (nele*nval*nt)                    *
 *                                                                           * 
 *****************************************************************************/

int c_mrbrep(void *buffer, int blkno, word *tblval)
{
   burp_block_header *block;
   buffer_interface_ptr buf = (buffer_interface_ptr) buffer;
   int nele,nval,nt,nelements;
   int btyp,nbit,bit0,datyp,bitpos,err,done,i,rmask,nnbits,allones;
   int nmots,new_nmots,diff_nmots;
   int new_nbit,new_datyp,nombre;
   burp_record *burprec;

   if ((blkno < 1) || (blkno > buf->buf78.buf8)) {
      sprintf(errmsg,"invalid block number");
      return(error_msg("c_mrbrep",BURP_ERR_BNUM,ERROR));
      }

   burprec = (burp_record *) buf->data;
   block = (burp_block_header *) burprec->data;

   if (block[blkno-1].flag != 0) {
      nele = block[blkno-1].elem1;
      nval = block[blkno-1].elem2;
      nt = block[blkno-1].elem3;
      nelements = nele;
      done = 0;
      }
   else {
      nele = block[blkno-1].nele;
      nval = block[blkno-1].nval;
      nt = block[blkno-1].nt;
      nelements = ((nele-3) > 0) ? nele-3 : 0;
      done = (nele < 3) ? nele : 3;
      }
   nbit = block[blkno-1].nbit + 1;
   datyp = block[blkno-1].datyp;
   bit0 = block[blkno-1].bit0;

#if defined(NEC)
   if (datyp > 6) {
      sprintf(errmsg,"datyp 6,7,8 or 9 not allowed on the NEC");
      return(error_msg("c_mrbrep",ERR_BAD_DATYP,ERRFATAL));
      }
#endif

   bitpos = bit0 * 64 + buf->buf9;

   nnbits = (((nele-done)*16 + 63) / 64) * 64;
   bitpos = bitpos + nnbits;
   nombre = nele * nval * nt;
   new_nbit = nbit;
   new_datyp = datyp;

/* check if precision should be ajusted */
   err = burp_nbit_datyp(&new_nbit,&new_datyp,tblval,nombre,xdf_stride);
   if (err < 0) return(err);

   if (new_datyp != datyp)
      block[blkno-1].datyp = new_datyp;

   if (new_nbit != nbit) {
      block[blkno-1].nbit = new_nbit -1;
      nmots = (nbit * nombre +63) / 64;
      new_nmots = (new_nbit * nombre + 63) /64;
      diff_nmots = new_nmots - nmots;
      new_nmots = diff_nmots * 2;    /* insertion by 32 bit slices */
      if (new_nmots > 0) {
         err = c_xdfins((word *)buf,tblval,bitpos,new_nmots,32,0);
	 for (i=blkno; i < buf->buf78.buf8; i++)
	    block[i].bit0 += diff_nmots;
	 }
      }
   
/* replace block */
   err = c_xdfrep((word *)buf,tblval,bitpos,nombre,new_nbit,new_datyp);

/* set all bits of missing values to 1 
   if ((datyp == 2) || (datyp == 4) || (datyp == 6)) {
      allones = ~(0);
      rmask = ~(-1 << nbit);
      for (i=0; i < nombre; i++)
	 if (tblval[i] == rmask) tblval[i] = allones;
      }
*/
   return(err);
}

/*splitpoint c_mrbxtr */
/***************************************************************************** 
 *                             C _ M R B X T R                               *
 *                                                                           * 
 *Object                                                                     * 
 *   Extract list of element and values from buffer.                         *
 *                                                                           *
 *Arguments                                                                  * 
 *                                                                           * 
 *    IN   buffer vector that contains the report                            * 
 *    IN   bkno   number of blocks in buf                                    *
 *   OUT   lstele list of nele meteorogical elements                         *
 *   OUT   tblval array of values to write (nele*nval*nt)                    *
 *                                                                           * 
 *****************************************************************************/

int c_mrbxtr(void *buffer, int bkno, word *lstele, word *tblval)
{
   burp_block_header entete;
   buffer_interface_ptr buf = (buffer_interface_ptr) buffer;
   int nele,nval,nt,nelements;
   int btyp,nbit,bit0,datyp,bitpos,err,done,i,rmask,nnbits,nombre,allones;
   int in_header;
   word *pos;

   err = 0;
   if ((bkno < 1) || (bkno > buf->buf78.buf8)) {
      sprintf(errmsg,"invalid block number");
      return(error_msg("c_mrbxtr",BURP_ERR_BNUM,ERROR));
      }

   bitpos = NBENTB * (bkno -1);
   err = c_xdfxtr((word *)buf,(word *)&entete,bitpos,DIMENT,8*sizeof(word),0);
   if (err < 0) return(err);

   if (entete.flag != 0) {
      nele = entete.elem1;
      nval = entete.elem2;
      nt = entete.elem3;
      nelements = nele;
      done = 0;
      }
   else {
      nele = entete.nele;
      nval = entete.nval;
      nt = entete.nt;
      in_header = (nele < 3) ? nele : 3;
      nelements = ((nele-in_header) > 0) ? nele-in_header : 0;
      if (in_header > 0) lstele[0] = entete.elem1;
      if (in_header > 1) lstele[1] = entete.elem2;
      if (in_header > 2) lstele[2] = entete.elem3;
      lstele += in_header;
      done = (nele < in_header) ? nele : in_header;
      }
   BurP_nele = nele;
   BurP_ntot = nele*nval*nt;
   nbit = entete.nbit + 1;
   datyp = entete.datyp;
   bit0 = entete.bit0;

#if defined(NEC)
   if (datyp > 6) {
      sprintf(errmsg,"datyp 6,7,8 or 9 not allowed on the NEC");
      return(error_msg("c_mrbxtr",ERR_BAD_DATYP,ERRFATAL));
      }
#endif

   bitpos = bit0 * 64 + buf->buf9;

/* extract list of elements */
   if (nelements > 0) 
      err = c_xdfxtr((word *)buf,lstele,bitpos,nelements,16,2);

/* extract array of values */
   nnbits = (((nele-done)*16 + 63) / 64) * 64;
   bitpos = bitpos + nnbits;
   nombre = nele * nval * nt;
   err = c_xdfxtr((word *)buf,tblval,bitpos,nombre,nbit,datyp);

   if ((datyp == 2) || (datyp == 4) || (datyp == 6)) {
/* set all bits of missing values to 1 */
      allones = ~(0);
      rmask = (nbit >= 32) ? allones : ~(-1 << nbit);
/*      rmask = ~(-1 << nbit); */
      for (i=0; i < nombre; i++) {
	 if (tblval[i] == rmask) tblval[i] = allones;
      }
   }
   return(err);
}

/*splitpoint c_mrfapp */
/***************************************************************************** 
 *                            C _ M R F A P P                                *
 *                                                                           * 
 *Object                                                                     * 
 *   Position at the end of a sequential file for an append.                 *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *                                                                           * 
 *****************************************************************************/

int c_mrfapp(int iun)
{
  int index_fnom, index, width, nw, end_of_file;
  file_table_entry *f;
  xdf_record_header *header;
  seq_dir_keys * seq_entry;

  index_fnom = fnom_index(iun);
  if (index_fnom == -1) {
    sprintf(errmsg,"file (unit=%d) is not connected with fnom",iun);
    return(error_msg("c_mrfapp",ERR_NO_FNOM,ERROR));
  }

  if ((index = file_index(iun)) == ERR_NO_FILE) {
    sprintf(errmsg,"file (unit=%d) is not open",iun);
    return(error_msg("c_mrfapp",ERR_NO_FILE,ERROR));
  }

  f = file_table[index];

  if (!f->xdf_seq) {
    sprintf(errmsg,"file (unit=%d) is not sequential",iun);
    return(error_msg("c_mrfapp",ERR_BAD_FTYPE,WARNING));
  }

  end_of_file = 0;
  width = W64TOWD(f->primary_len);

  while (!end_of_file) {
    nw = c_waread2(iun,f->head_keys,f->cur_addr,width);
    header = (xdf_record_header *) f->head_keys;
    if (nw < W64TOWD(1)) {
      end_of_file = 1;
      header->idtyp = 127;
      header->lng = 1;
      break;
    }
    if ((header->idtyp >= 112) && (header->idtyp <= 127)) {
      end_of_file = 1;
      break;
    }
    f->cur_addr += W64TOWD(header->lng);
  }
  f->nxtadr = f->cur_addr;
  return(0);
}

/*splitpoint c_mrfapp */
/***************************************************************************** 
 *                            C _ M R F B F L                                *
 *                                                                           * 
 *Object                                                                     * 
 *   Return the length of the longer report in the file                      *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *                                                                           * 
 *****************************************************************************/

int c_mrfbfl(int iun)
{
  int bfl;
  
  bfl = f77name(mrfbfl)(&iun);
  return(bfl);
}

  /*splitpoint c_mrfget */
/***************************************************************************** 
 *                             C _ M R F G E T                               *
 *                                                                           * 
 *Object                                                                     * 
 *   Read the report referenced by handle from the file.                     *
 *                                                                           *
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN   handle logical pointer to the record                                *
 *  OUT  buffer vector that contains the report                              * 
 *                                                                           * 
 *****************************************************************************/
int c_mrfget(int handle, void *buffer)
{
   buffer_interface_ptr buf = (buffer_interface_ptr) buffer;
   burp_record *burprec;
   int err;

   err = c_xdfget(handle,buf);
   if (err < 0) {
     return(error_msg("c_mrfget",err,ERROR));
   }
   if (msg_level <= INFORM) {
     fprintf(stdout,"RECORD READ\n");
   }

   burprec = (burp_record *) buf->data;
   buf->buf78.buf8 =  burprec->info.nblks;
   buf->buf9 =  burprec->info.nblks * NBENTB;
   
   return(0);
}

/*splitpoint c_mrfput */
/***************************************************************************** 
 *                             C _ M R F P U T                               *
 *                                                                           * 
 *Object                                                                     * 
 *   Write a report to the file. If handle is not 0, record referenced by    *
 *   by handle is written at end of file. If handle is 0, a new record is    *
 *   written.  If hanlde is > 0, it will be forced to be negative to write   *
 *   at end of file.                                                         *
 *                                                                           *
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN   iun    unit number associated to the file                           *
 *  IN   handle logical pointer to the record                                *
 *  IN   buffer vector that contains the report                              * 
 *                                                                           * 
 *****************************************************************************/
int c_mrfput(int iun, int handle, void *buffer)
{
   buffer_interface_ptr buf = (buffer_interface_ptr) buffer;
   burp_record *burprec;
   int err,new_handle;
   int temps,flgs,idtyp,lat,lon,dx,dy,elev,drnd,date;
   int oars,runn,nblk,sup,nsup,xaux,nxaux;
   char stnid[10];

   burprec = (burp_record *) buf->data;
   burprec->info.nblks = buf->buf78.buf8;

   new_handle = (handle > 0) ? -handle : handle;
   err = c_xdfput(iun,new_handle,buf);
   if (msg_level <= INFORM) {
      nsup = 0;
      nxaux = 0;
      err = c_mrbhdr((word *)buf,&temps,&flgs,stnid,&idtyp,&lat,&lon,
		     &dx,&dy,&elev,&drnd,&date,&oars,&runn,&nblk,
		     (word *)&sup,nsup,(word *)&xaux,nxaux);
      stnid[9] = '\0';
      fprintf(stdout,"MRFPUT - WRITE - STNID=%s IDTYP=%3d LAT=%5d LON=%5d DX=%4d DY=%4d DATE=%8d TEMPS=%4d, FLGS=%8d\n",
	      stnid,idtyp,lat,lon,dx,dy,date,temps,flgs);
      }
   return(0);
}

/*splitpoint c_mrfrwd */
/***************************************************************************** 
 *                           C _ M R F R W D                                 *
 *                                                                           * 
 *Object                                                                     * 
 *   Rewinds a BURP sequential file.                                         *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *                                                                           * 
 *****************************************************************************/
int c_mrfrwd(int iun)
{
  int index, index_fnom;
  file_table_entry *f;

  index_fnom = fnom_index(iun);
  if (index_fnom == -1) {
    sprintf(errmsg,"file (unit=%d) is not connected with fnom",iun);
    return(error_msg("c_mrfrwd",ERR_NO_FNOM,ERROR));
  }

  if ((index = file_index(iun)) == ERR_NO_FILE) {
    sprintf(errmsg,"file (unit=%d) is not open",iun);
    return(error_msg("c_mrfrwd",ERR_NO_FILE,ERROR));
  }

  f = file_table[index];

  if (! f->cur_info->attr.burp) {
    sprintf(errmsg,"file (unit=%d) is not a BURP file",iun);
    return(error_msg("c_fstrwd",ERR_NO_FILE,ERROR));
  }

  if (! f->xdf_seq) {
    sprintf(errmsg,"file (unit=%d) is not sequential",iun);
    return(error_msg("c_mrfrwd",ERR_BAD_FTYPE,WARNING));
  }

  f->cur_addr = f->seq_bof;
  f->valid_pos = 0;
  return(0);
}

#include "if_burp98.hc"
