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

/*splitpoint buf89a0 */
/***************************************************************************** 
 *                            B U F 8 9 A 0                                  * 
 *****************************************************************************/

void f77name(buf89a0)(word *buf)
{
#if defined(NEC64)
  BUF_C;
  c_buf89a0(buf+1);
  BUF_F;
#else
  c_buf89a0(buf);
#endif
}

/*splitpoint getbuf8 */
/***************************************************************************** 
 *                            G E T B U F 8                                  * 
 *****************************************************************************/

ftnword f77name(getbuf8)(word *buf)
{
  int n;
#if defined(NEC64)
  BUF_C;
  n = c_getbuf8(buf+1);
  BUF_F;
#else
  n = c_getbuf8(buf);
#endif
  return((ftnword) n);
}

/*splitpoint genvdt8 */
/***************************************************************************** 
 *                            G E N V D T 8                                  * 
 *****************************************************************************/

void f77name(genvdt8)(ftnword *val)
{
  char *enforc8;
  
  enforc8 = getenv("ENFORC8");
  if (enforc8) {
    *val = 1;
    xdf_enforc8 = 1;
  }
  else {
    *val = 0;
    xdf_enforc8 = 0;
  }
}

/*splitpoint mrbadd */
/***************************************************************************** 
 *                              M R B A D D                                  * 
 *****************************************************************************/

ftnword f77name(mrbadd)(word *buf, ftnword *f_bkno, ftnword *f_nele,
			ftnword *f_nval, ftnword *f_nt, ftnword *f_bfam,
			ftnword *f_bdesc, ftnword *f_btyp, ftnword *f_nbit,
			ftnword *f_bit0, ftnword *f_datyp,
			ftnword *lstele, word *tblval)
{
   int bkno=*f_bkno, nele=*f_nele, nval=*f_nval, nt=*f_nt, bfam=*f_bfam;
   int bdesc=*f_bdesc, btyp=*f_btyp, nbit=*f_nbit, bit0=*f_bit0;
   int datyp=*f_datyp, err;

#if defined(NEC64)
   INT_32 listele[1024], *pliste;
   int was_allocated = 0, i;

   BUF_C;
   xdf_stride = 2;
   if (nele > 1024) {
     pliste = calloc(nele,sizeof(INT_32));
     was_allocated = 1;
   }
   else 
     pliste = listele;
   for (i=0; i < nele; i++)
     pliste[i] = lstele[i];
   err = c_mrbadd(buf+1,&bkno,nele,nval,nt,bfam,bdesc,btyp,nbit,&bit0,
		  datyp,pliste,tblval+1);
   if (was_allocated) free(pliste);
   xdf_stride = 1;
   BUF_F;
#else
   err = c_mrbadd(buf,&bkno,nele,nval,nt,bfam,bdesc,btyp,nbit,&bit0,
		  datyp,(word *)lstele,tblval);
#endif
   *f_bit0 = (ftnword) bit0;
   *f_bkno = (ftnword) bkno;
   return(err);
}

/*splitpoint mrbdel */
/***************************************************************************** 
 *                              M R B D E L                                  * 
 *****************************************************************************/

ftnword f77name(mrbdel)(word *buf, ftnword *f_number)
{
   int number = *f_number;
   int ier;
#if defined(NEC64)
   BUF_C;
   ier = c_mrbdel(buf+1,number);
   BUF_F;
#else
   ier = c_mrbdel(buf,number);
#endif
   return((ftnword) ier);
   
}

/*splitpoint mrbhdr */
/***************************************************************************** 
 *                              M R B H D R                                  *
 *****************************************************************************/
ftnword f77name(mrbhdr)(word *buf, ftnword *f_temps,
			ftnword *f_flgs, char *f_stnid,
			ftnword *f_idtyp, ftnword *f_lati, ftnword *f_lon,
			ftnword *f_dx, ftnword *f_dy, ftnword *f_elev,
			ftnword *f_drcv, ftnword *f_date, ftnword *f_oars,
			ftnword *f_run, ftnword *f_nblk,
			ftnword *f_sup, ftnword *f_nsup,
			ftnword *f_xaux, ftnword *f_nxaux, int ll1)
{
int temps, flgs, idtyp, lati, lon;
int dx, dy,elev, drcv, date, oars, run, nblk;
int nsup=*f_nsup, nxaux=*f_nxaux;
char stnid[11] = {' ',' ',' ',' ',' ',' ',' ',' ',' ',' ','\0'};
int ier, l1;

#if defined(NEC64)
   BUF_C;
   ier = c_mrbhdr(buf+1, &temps, &flgs, stnid, &idtyp,
		  &lati, &lon, &dx, &dy, &elev,
		  &drcv, &date, &oars, &run, &nblk,
		  (word *)f_sup, nsup, (word *)f_xaux, nxaux);
   BUF_F;
#else
   ier = c_mrbhdr(buf, &temps, &flgs, stnid, &idtyp,
		  &lati, &lon, &dx, &dy, &elev,
		  &drcv, &date, &oars, &run, &nblk,
		  (word *)f_sup, nsup, (word *)f_xaux, nxaux);
#endif
   *f_temps = (ftnword) temps;
   *f_flgs = (ftnword) flgs;
   *f_idtyp = (ftnword) idtyp;
   *f_lati  = (ftnword) lati;
   *f_lon = (ftnword) lon;
   *f_dx = (ftnword) dx;
   *f_dy = (ftnword) dy;
   *f_elev = (ftnword) elev;
   *f_drcv = (ftnword) drcv;
   *f_date = (ftnword) date;
   *f_oars = (ftnword) oars;
   *f_run = (ftnword) run;
   *f_nblk = (ftnword) nblk;
   l1 = (ll1 > 11) ? 11: ll1;
   string_copy(f_stnid,stnid,l1);
   return((ftnword) ier);
}

/*splitpoint mrblen */
/***************************************************************************** 
 *                              M R B L E N                                  * 
 *****************************************************************************/

ftnword f77name(mrblen)(word *buf, ftnword *f_lbits, ftnword *f_left)
{
   int lbits, left, err;
#if defined(NEC64)
   BUF_C;
   err = c_mrblen(buf+1,&lbits,&left);
   BUF_F;   
#else
   err = c_mrblen(buf,&lbits,&left);
#endif
   *f_lbits = (ftnword) lbits;
   *f_left = (ftnword) left;
   return((ftnword) err);
   
}

/*splitpoint mrbloc */
/***************************************************************************** 
 *                              M R B L O C                                  * 
 *****************************************************************************/

ftnword f77name(mrbloc)(word *buf, ftnword *f_bfam, ftnword *f_bdesc,
			ftnword *f_btyp, ftnword *f_blkno)
{
   int blkno=*f_blkno, bfam=*f_bfam;
   int bdesc=*f_bdesc, btyp=*f_btyp;
   int ier;

#if defined(NEC64)
   BUF_C;
   ier = c_mrbloc(buf+1,bfam,bdesc,btyp,blkno);
   BUF_F;   
#else
   ier = c_mrbloc(buf,bfam,bdesc,btyp,blkno);
#endif
   return ((ftnword) ier);
}

/*splitpoint mrbrep */
/***************************************************************************** 
 *                              M R B R E P                                  * 
 *****************************************************************************/

ftnword f77name(mrbrep)(word *buf, ftnword *f_blkno, word *tblval)
{
   int blkno=*f_blkno;
   int ier;
#if defined(NEC64)
   xdf_stride = 2;
   BUF_C;
   ier = c_mrbrep(buf+1,blkno,tblval+1);
   xdf_stride = 1;
   BUF_F;   
#else
   ier = c_mrbrep(buf,blkno,tblval);
#endif
   return((ftnword) ier);
}

/*splitpoint mrbxtr */
/***************************************************************************** 
 *                              M R B X T R                                  * 
 *****************************************************************************/

ftnword f77name(mrbxtr)(word *buf, ftnword *f_bkno, ftnword *lstele,
			ftnword *tblval)
{
   int bkno=*f_bkno, ier, i;
#if defined(NEC64)
   INT_32 *plong;
   BUF_C;
   ier = c_mrbxtr(buf+1,bkno,(word *)lstele,(word *)tblval);
   BUF_F;   
   plong = (INT_32 *) lstele;
   for (i=BurP_nele-1; i >= 0; i--)
     lstele[i] = plong[i];
   plong = (INT_32 *) tblval;
   for (i=BurP_ntot-1; i >= 0; i--)
     tblval[i] = plong[i];
#else
   ier = c_mrbxtr(buf,bkno,(word *) lstele,(word *) tblval);
#endif
   return((ftnword) ier);
}
/*splitpoint fstapp */
/***************************************************************************** 
 *                             M R F A P P                                   *
 *****************************************************************************/

ftnword f77name(mrfapp)(ftnword *f_iun)
{
  int ier, iun = *f_iun;
  
  ier = c_mrfapp(iun);
  return((ftnword) ier);
}



/*splitpoint mrfget */
/***************************************************************************** 
 *                              M R F G E T                                  * 
 *****************************************************************************/

ftnword f77name(mrfget)(ftnword *f_handle, word *buf)
{
   int handle=*f_handle, ier;
#if defined(NEC64)
   BUF_C;
   ier = c_mrfget(handle,buf+1);
   BUF_F;
#else
   ier = c_mrfget(handle,buf);
#endif
   return((ftnword) ier);
}

/*splitpoint mrfput */
/***************************************************************************** 
 *                              M R F P U T                                  * 
 *****************************************************************************/

ftnword f77name(mrfput)(ftnword *f_iun, ftnword *f_handle, word *buf)
{
   int iun=*f_iun, handle=*f_handle,ier;
   
#if defined(NEC64)
   BUF_C;
   ier = c_mrfput(iun,handle,buf+1);
   BUF_F;
#else
   ier = c_mrfput(iun,handle,buf);
#endif
   return((ftnword) ier);
}

/*splitpoint mrfrwd */
/***************************************************************************** 
 *                               M R F R W D                                 *
 *****************************************************************************/
ftnword f77name(mrfrwd)(ftnword *f_iun)
{
  int err, iun = *f_iun;

  err = c_mrfrwd(iun);
  return((ftnword) err);
}
