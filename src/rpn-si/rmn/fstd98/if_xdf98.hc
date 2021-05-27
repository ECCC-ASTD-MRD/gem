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

/*splitpoint qdfdiag */
/***************************************************************************** 
 *                             Q D F D I A G                                 *
 *****************************************************************************/

ftnword  f77name(qdfdiag)(ftnword *f_iun)
{
  int iun = *f_iun, ier;

  ier = c_qdfdiag(iun);
  return((ftnword) ier);
}
/*splitpoint qdferr */
/***************************************************************************** 
 *                              Q D F E R R                                  * 
 *****************************************************************************/

ftnword f77name(qdferr)(char *subname,char *msg,ftnword *ferrlevl,
			ftnword *ferrcode,F2Cl l1,F2Cl l2)
{
   int errlevl = *ferrlevl, errcode = *ferrcode, lng;
   char c_subname[128];

   errcode = (errcode > 0) ? -errcode : errcode;
   lng = (l1 < 128) ? l1 : 127;
   strncpy(c_subname,subname,lng);
   c_subname[lng] = '\0';

   lng = (l2 < 1024) ? l2 : 1023;
   strncpy(errmsg,msg,lng);

   return((ftnword) error_msg(c_subname,errcode,errlevl));
}

/*splitpoint qdfind */
/***************************************************************************** 
 *                              Q D F I N D                                  * 
 *****************************************************************************/

ftnword f77name(qdfind)(ftnword *iun)
{
   int ind;
   ind = file_index(*iun);
   ind = (ind != ERR_NO_FILE) ? ind : 9999;
   
   return((ftnword) ind);
}

/*splitpoint qdfmsig */
/***************************************************************************** 
 *                            Q D F M S I G                                  * 
 *****************************************************************************/

ftnword f77name(qdfmsig)(ftnword *fiun,char *appl,F2Cl l1)
{
   int iun = *fiun, lng;
   char c_appl[257];

   lng = (l1 <= 256) ? l1 : 256;
   strncpy(c_appl,appl,lng);
   c_appl[lng] = '\0';

   return(c_qdfmsig(iun,c_appl));
}

/*splitpoint qdfput */
/***************************************************************************** 
 *                              Q D F P U T                                  * 
 *****************************************************************************/

ftnword f77name(qdfput)(word *buf, ftnword *felem, ftnword *fderbit,
			ftnword *fnbits)
{
   int elem = *felem, nbits = *fnbits, derbit = *fderbit;
   int ier;
#if defined(NEC64)
   BUF_C;
   ier = c_qdfput(buf+1,elem,derbit,nbits);
   BUF_F;
#else
   ier = c_qdfput(buf,elem,derbit,nbits);
#endif
   return((ftnword) ier);
}
/*splitpoint qdfrstr */
/***************************************************************************** 
 *                             Q D F R S T R                                 *
 *****************************************************************************/

ftnword  f77name(qdfrstr)(ftnword *f_inp, ftnword *f_outp)
{
  int inp = *f_inp, outp = *f_outp, ier;

  ier = c_qdfrstr(inp,outp);
  return((ftnword) ier);
}

/*splitpoint rewind_file */
/***************************************************************************** 
 *                          R E W I N D _ F I L E                            * 
 *****************************************************************************/

/* set file position, if handle=-1, rewind file */
static INT_32 rewind_file(int file_index, int handle)
{
   register file_table_entry *f, *f2;
   int linked=0, file_index2;

   /* check if file exists */
   if( (f = file_table[file_index]) == NULL) return (ERR_NO_FILE);

   if (handle==-1){
      f->cur_dir_page = f->dir_page[0];
      f->cur_entry = (f->cur_dir_page)->dir.entry;
      f->page_nrecords = (f->cur_dir_page)->dir.nent;
      f->page_record = 0;
   }
   else{
      file_index2=INDEX_FROM_HANDLE(handle);

      if(file_index!=file_index2){
         /* check if file2 exists and is linked to file */
         if( (f2 = file_table[file_index2]) == NULL) return (ERR_NO_FILE);
         f2=f;
         linked=f2->link;
         while((linked>0) && (linked!=file_index2)){
            f2=file_table[linked];
            linked=f2->link;
         };
         if(linked!=file_index2) return (ERR_BAD_LINK);
      }
      f2=file_table[file_index2];
      f->cur_dir_page = f2->dir_page[PAGENO_FROM_HANDLE(handle)];
      f->page_record = RECORD_FROM_HANDLE(handle);
      f->page_nrecords = (f->cur_dir_page)->dir.nent;
      f->cur_entry = (f->cur_dir_page)->dir.entry + (f->page_record)*(f->primary_len);
   }
   return(0);
}

/*splitpoint xdfadd */
/***************************************************************************** 
 *                              X D F A D D                                  * 
 *****************************************************************************/

ftnword f77name(xdfadd)(word *buf,word *donnees,
                        ftnword *fnelm, ftnword *fnbits, ftnword *fdatyp)
{
   int nelm = *fnelm, nbits = *fnbits, datyp = *fdatyp;
   int ier;

#if defined(NEC64)
   BUF_C;
   if ((datyp == 2) || (datyp == 4)) {
     xdf_stride = 2;
     ier = c_xdfadd(buf+1,donnees+1,nelm,nbits,datyp);
     xdf_stride = 1;
   }
   else
     ier = c_xdfadd(buf+1,donnees,nelm,nbits,datyp);
   BUF_F;
#else
   ier = c_xdfadd(buf,donnees,nelm,nbits,datyp);
#endif
   return((ftnword) ier);

}


/*splitpoint xdfcle */
/***************************************************************************** 
 *                                X D F C L E                                *
 *****************************************************************************/
ftnword f77name(xdfcle)(char *fkeyname,ftnword *fbit1,ftnword *flkey,
			ftnword *ftkey,ftnword *fdesc1,ftnword *fdesc2,F2Cl l1)
{
   char keyname[5]={' ',' ',' ',' ','\0'};
   int lkey = *flkey, tkey = *ftkey, bit1 = *fbit1;
   int desc1, desc2;
   int lng, err;

   lng = (l1 <= 4) ? l1 : 4;
   strncpy(keyname,fkeyname,lng);
   
   err = c_xdfcle(keyname,bit1,lkey,tkey,&desc1,&desc2);
   
   *fdesc1 = (ftnword) desc1;
   *fdesc2 = (ftnword) desc2;
   
   return(err);
}

/*splitpoint xdfcls */
/***************************************************************************** 
 *                              X D F C L S                                  * 
 *****************************************************************************/

ftnword f77name(xdfcls)(ftnword *fiun)
{
   int iun = *fiun;

   return(c_xdfcls(iun));
}

/*splitpoint xdfcut */
/***************************************************************************** 
 *                              X D F C U T                                  * 
 *****************************************************************************/

ftnword f77name(xdfcut)(word *buf,
			ftnword *fbitpos, ftnword *fnelm,
			ftnword *fnbits, ftnword *fdatyp)
{
   int nelm = *fnelm, nbits = *fnbits, datyp = *fdatyp, bitpos = *fbitpos;
   int ier;

#if defined(NEC64)
   BUF_C;
   ier = c_xdfcut(buf+1,bitpos,nelm,nbits,datyp);
   BUF_F;
#else
   ier = c_xdfcut(buf,bitpos,nelm,nbits,datyp);
#endif
   return((ftnword) ier);
}

/*splitpoint xdfdel */
/***************************************************************************** 
 *                              X D F D E L                                  * 
 *****************************************************************************/

ftnword f77name(xdfdel)(ftnword *fhandle)
{
   int handle = *fhandle;

   return((ftnword) c_xdfdel(handle));
}

/*splitpoint xdfget */
/***************************************************************************** 
 *                              X D F G E T                                  * 
 *****************************************************************************/

ftnword f77name(xdfget)(ftnword *fhandle, word *buf)
{
   int handle = *fhandle, ier;

#if defined(NEC64)
   BUF_C;
   ier = c_xdfget(handle,buf+1);
   BUF_F;
#else
   ier = c_xdfget(handle,(buffer_interface_ptr) buf);
#endif
   return((ftnword) ier);
}

/*splitpoint xdfgop */
/***************************************************************************** 
 *                              X D F G O P                                  * 
 *****************************************************************************/

ftnword f77name(xdfgop)(char *foptname, char *foptc, ftnword *foptv,
			F2Cl ll1, F2Cl ll2)
{
   int optv, err, l1=ll1, l2=ll2;
   char optname[257], optc[257];

   l1 = (l1 <= 256) ? l1 : 256; 
   strncpy(optname,foptname,l1);
   optname[l1] = '\0';

   err = c_xdfgop(optname,optc,&optv);

   l2 = (l2 <= 256) ? l2 : 256;
   strncpy(foptc,optc,l2);

   *foptv = (ftnword) optv;
   return((ftnword) err);

}

/*splitpoint xdfhdr */
/***************************************************************************** 
 *                              X D F H D R                                  * 
 *****************************************************************************/

ftnword f77name(xdfhdr)(word *buf,ftnword *addr,ftnword *lng,
                        ftnword *idtyp,ftnword *primk,ftnword *fnprim,
			ftnword *info, ftnword *fninfo)
{
   int nprim = *fnprim, ninfo = *fninfo, ier, i;
   int l_addr, l_lng, l_idtyp;
   word l_primk[MAX_KEYS];
   word l_info[MAX_KEYS];

#if defined(NEC64)
   BUF_C;
   ier = c_xdfhdr(buf+1,&l_addr,&l_lng,&l_idtyp,l_primk,nprim,l_info,ninfo);
   BUF_F;
#else
   ier = c_xdfhdr((buffer_interface_ptr)buf,&l_addr,&l_lng,&l_idtyp,l_primk,nprim,l_info,ninfo);
#endif

   *addr =  (ftnword) l_addr;
   *lng =   (ftnword) l_lng;
   *idtyp = (ftnword) l_idtyp;

   if ((nprim > MAX_KEYS) || (ninfo >MAX_KEYS)) {
      sprintf(errmsg,"nprim=%d or ninfo=%d > MAX_KEYS must recompile",nprim,ninfo);
      return(error_msg("xdfhdr",ERR_OUT_RANGE,SYSTEM));
      }

   for (i=0; i < nprim; i++)
      primk[i] = (ftnword) l_primk[i];

   for (i=0; i < ninfo; i++)
      info[i] = (ftnword) l_info[i];

   return ((ftnword) ier);
}

/*splitpoint xdfimp */
/***************************************************************************** 
 *                              X D F I M P                                  * 
 *****************************************************************************/

ftnword f77name(xdfimp)(ftnword *fiun,ftnword *stat,ftnword *fnstat,
                    ftnword_2 *pri,ftnword_2 *aux,
                    char *vers,char *appl,F2Cl l1,F2Cl l2)
{
   int iun = *fiun, nstat = *fnstat,lng, ier, i, nkeys, ninfo;
   char c_vers[257], c_appl[257];
   word_2 primk[MAX_KEYS], infok[MAX_KEYS];
   word lstat[12];
   
   lng = (l1 <= 256) ? l1 : 256;
   strncpy(c_vers,vers,lng);
   c_vers[lng] = '\0';

   lng = (l2 <= 256) ? l2 : 256;
   strncpy(c_appl,appl,lng);
   c_appl[lng] = '\0';

#if defined(NEC64)
   for (i=0; i < nstat; i++)
     lstat[i] = stat[i];
   nkeys = lstat[6];
   ninfo = lstat[8];
   if ((nkeys > MAX_KEYS) || (ninfo >MAX_KEYS)) {
      sprintf(errmsg,"nkeys=%d or ninfo=%d > MAX_KEYS must recompile",nkeys,ninfo);
      return(error_msg("xdfimp",ERR_OUT_RANGE,SYSTEM));
      }
   for (i=0; i < nkeys; i++) {
     primk[i].wd1 = pri[i].wd1;
     primk[i].wd2 = pri[i].wd2;
     }
   for (i=0; i < ninfo; i++) {
     infok[i].wd1 = aux[i].wd1;
     infok[i].wd2 = aux[i].wd2;
     }
   ier = c_xdfimp(iun,lstat,nstat,primk,infok,c_vers,c_appl);
#else
   ier = c_xdfimp(iun,stat,nstat,pri,aux,c_vers,c_appl);
#endif
   return((ftnword) ier);
}

/*splitpoint xdfini */
/***************************************************************************** 
 *                              X D F I N I                                  * 
 *****************************************************************************/

ftnword f77name(xdfini)(ftnword *fiun,word *buf,ftnword *fidtyp,
			ftnword *keys,ftnword *fnkeys,ftnword *info,
			ftnword *fninfo)
{
   int iun = *fiun, idtyp = *fidtyp, nkeys = *fnkeys, ninfo = *fninfo;
   int ier, i;
#if defined(NEC64)
   word primk[MAX_KEYS], infok[MAX_KEYS];

   if ((nkeys > MAX_KEYS) || (ninfo >MAX_KEYS)) {
      sprintf(errmsg,"nkeys=%d or ninfo=%d > MAX_KEYS must recompile",nkeys,ninfo);
      return(error_msg("xdfini",ERR_OUT_RANGE,SYSTEM));
      }
   for (i=0; i < nkeys; i++) 
     primk[i] = keys[i];
   for (i=0; i < ninfo; i++) 
     infok[i] = info[i];
   BUF_C;
   ier = c_xdfini(iun,buf+1,idtyp,primk,nkeys,infok,ninfo);
   BUF_F;
#else
   ier = c_xdfini(iun,(buffer_interface_ptr)buf,idtyp,keys,nkeys,info,ninfo);
#endif
   return((ftnword) ier);
}

/*splitpoint xdfins */
/***************************************************************************** 
 *                              X D F I N S                                  * 
 *****************************************************************************/

ftnword f77name(xdfins)(word *buf,word *donnees,
                        ftnword *fbitpos, ftnword *fnelm,
			ftnword *fnbits, ftnword *fdatyp)
{
   int nelm = *fnelm, nbits = *fnbits, datyp = *fdatyp, bitpos = *fbitpos;
   int ier;

#if defined(NEC64)
   BUF_C;
   if ((datyp == 2) || (datyp == 4)) {
     xdf_stride = 2;
     ier = c_xdfins(buf+1,donnees+1,bitpos,nelm,nbits,datyp);
     xdf_stride = 1;
   }
   else
     ier = c_xdfins(buf+1,donnees,bitpos,nelm,nbits,datyp);
   BUF_F;
#else
   ier = c_xdfins(buf,donnees,bitpos,nelm,nbits,datyp);
#endif
   return((ftnword) ier);
}

/*splitpoint xdflnk */
/***************************************************************************** 
 *                              X D F L N K                                  * 
 *****************************************************************************/

ftnword f77name(xdflnk)(ftnword *liste, ftnword *fn)
{
   int n = *fn, ier;

   ier = c_xdflnk(liste,n);
   return((ftnword) ier);

}

/*splitpoint xdfloc */
/***************************************************************************** 
 *                              X D F L O C                                  * 
 *****************************************************************************/

ftnword f77name(xdfloc)(ftnword *fiun, ftnword *fhandle, ftnword *primk,
			ftnword *fnprim)
{
   int iun = *fiun, nprim = *fnprim, i;
   int handle = *fhandle;
   word l_primk[MAX_KEYS];

   if (nprim > MAX_KEYS) {
      sprintf(errmsg,"nprim=%d > MAX_KEYS must recompile",nprim);
      return(error_msg("xdfloc",ERR_OUT_RANGE,SYSTEM));
      }
   for (i=0; i<nprim; i++)
      l_primk[i] = primk[i];

   return((ftnword) c_xdfloc(iun,handle,l_primk,nprim));

}

/*splitpoint xdfopn */
/***************************************************************************** 
 *                              X D F O P N                                  * 
 *****************************************************************************/

ftnword f77name(xdfopn)(ftnword *fiun, char *mode,
			ftnword_2 *pri, ftnword *fnpri,
			ftnword_2 *aux, ftnword *fnaux,
			char *appl, F2Cl l1, F2Cl l2)
{
   int iun = *fiun, npri = *fnpri, naux = *fnaux, lng, i, ier;
   char c_mode[257], c_appl[257];
   word_2 primk[MAX_KEYS], infok[MAX_KEYS];

   lng = (l1 <= 256) ? l1 : 256;
   strncpy(c_mode,mode,lng);
   c_mode[lng] = '\0';

   lng = (l2 <= 256) ? l2 : 256;
   strncpy(c_appl,appl,lng);
   c_appl[lng] = '\0';

   if ((npri > MAX_KEYS) || (naux >MAX_KEYS)) {
      sprintf(errmsg,"npri=%d or naux=%d > MAX_KEYS must recompile",
	      npri,naux);
      return(error_msg("xdfopn",ERR_OUT_RANGE,SYSTEM));
      }
   for (i=0; i < npri; i++) {
     primk[i].wd1 = pri[i].wd1;
     primk[i].wd2 = pri[i].wd2;
     }
   for (i=0; i < naux; i++) {
     infok[i].wd1 = aux[i].wd1;
     infok[i].wd2 = aux[i].wd2;
     }
   ier = c_xdfopn(iun,c_mode,primk,npri,infok,naux,c_appl);
   return((ftnword) ier);
}

/*splitpoint xdfopt */
/***************************************************************************** 
 *                              X D F O P T                                  * 
 *****************************************************************************/

ftnword f77name(xdfopt)(char *foptname, char *foptc, ftnword *foptv,
			F2Cl ll1, F2Cl ll2)
{
   int optv = *foptv, l1=ll1, l2=ll2;
   char optname[257], optc[257];

   l1 = (l1 <= 256) ? l1 : 256; 
   strncpy(optname,foptname,l1);
   optname[l1] = '\0';

   l2 = (l2 <= 256) ? l2 : 256;
   strncpy(optc,foptc,l2);
   optc[l2] = '\0';

   return((ftnword) c_xdfopt(optname,optc,optv));
}

/*splitpoint xdfprm */
/***************************************************************************** 
 *                              X D F P R M                                  * 
 *****************************************************************************/

ftnword f77name(xdfprm)(ftnword *fhandle,ftnword *addr,ftnword *lng,
                        ftnword *idtyp,ftnword *primk,ftnword *fnprim)
{
   int nprim = *fnprim, ier, i;
   int handle = *fhandle;
   int l_addr, l_lng, l_idtyp;
   word l_primk[MAX_KEYS];

   ier = c_xdfprm(handle,&l_addr,&l_lng,&l_idtyp,l_primk,nprim);
   *addr =  (ftnword) l_addr;
   *lng =   (ftnword) l_lng;
   *idtyp = (ftnword) l_idtyp;

   for (i=0; i < nprim; i++)
      primk[i] = (ftnword) l_primk[i];

   return ((ftnword) ier);
}

/*splitpoint xdfput */
/***************************************************************************** 
 *                              X D F P U T                                  * 
 *****************************************************************************/

ftnword f77name(xdfput)(ftnword *fiun, ftnword *fhandle,
			word *buf)
{
   int handle = *fhandle;
   int iun = *fiun, ier;

#if defined(NEC64)
   BUF_C;
   ier = c_xdfput(iun,handle,buf+1);
   BUF_F;
#else
   ier = c_xdfput(iun,handle,(buffer_interface_ptr)buf);
#endif
   return((ftnword) ier);
}

/*splitpoint xdfrep */
/***************************************************************************** 
 *                              X D F R E P                                  * 
 *****************************************************************************/

ftnword f77name(xdfrep)(word *buf,word *donnees,
                        ftnword *fbitpos, ftnword *fnelm,
			ftnword *fnbits, ftnword *fdatyp)
{
   int nelm = *fnelm, nbits = *fnbits, datyp = *fdatyp, bitpos = *fbitpos;
   int ier;

#if defined(NEC64)
   BUF_C;
   if ((datyp == 2) || (datyp == 4)) {
     xdf_stride = 2;
     ier = c_xdfrep(buf+1,donnees+1,bitpos,nelm,nbits,datyp);
     xdf_stride = 1;
   }
   else
     ier = c_xdfrep(buf+1,donnees,bitpos,nelm,nbits,datyp);
   BUF_F;
#else
   ier = c_xdfrep(buf,donnees,bitpos,nelm,nbits,datyp);
#endif
   return((ftnword) ier);
}

/*splitpoint xdfsta */
/***************************************************************************** 
 *                              X D F S T A                                  * 
 *****************************************************************************/

ftnword f77name(xdfsta)(ftnword *fiun,ftnword *stat,ftnword *fnstat,
			ftnword_2 *pri,ftnword *fnpri,
			ftnword_2 *aux,ftnword *fnaux,
			char *vers,char *appl,F2Cl l1,F2Cl l2)
{
   int iun = *fiun, npri = *fnpri, naux = *fnaux, nstat = *fnstat,lng, ier;
   char c_vers[257], c_appl[257];
   int i;

#if defined(NEC64)
   word_2 primk[MAX_KEYS], infok[MAX_KEYS];
   word lstat[12];

   ier = c_xdfsta(iun,lstat,nstat,primk,npri,infok,naux,c_vers,c_appl);

   if ((npri > MAX_KEYS) || (naux >MAX_KEYS)) {
      sprintf(errmsg,"npri=%d or naux=%d > MAX_KEYS must recompile",
	      npri,naux);
      return(error_msg("xdfsta",ERR_OUT_RANGE,SYSTEM));
      }
   for (i=0; i < npri; i++) {
     pri[i].wd1 = primk[i].wd1;
     pri[i].wd2 = primk[i].wd2;
     }
   for (i=0; i < naux; i++) {
     aux[i].wd1 = infok[i].wd1;
     aux[i].wd2 = infok[i].wd2;
     }
   for (i=0; i < nstat; i++)
     stat[i] = lstat[i];
#else
   ier = c_xdfsta(iun,stat,nstat,pri,npri,aux,naux,c_vers,c_appl);
#endif
   lng = (l1 <= 256) ? l1 : 256;
   c_vers[lng] = '\0';
   strncpy(vers,c_vers,lng);

   lng = (l2 <= 256) ? l2 : 256;
   c_appl[lng] = '\0';
   strncpy(appl,c_appl,lng);

   return((ftnword) ier);
   
}

/*splitpoint xdfupd */
/***************************************************************************** 
 *                              X D F U P D                                  * 
 *****************************************************************************/

ftnword f77name(xdfupd)(ftnword *fiun,word *buf,ftnword *fidtyp,
			ftnword *keys,ftnword *fnkeys,
			ftnword *info,ftnword *fninfo)
{
   int iun = *fiun, idtyp = *fidtyp, nkeys = *fnkeys, ninfo = *fninfo;
   int ier, i;
#if defined(NEC64)
   word l_keys[MAX_KEYS],l_info[MAX_KEYS];

   BUF_C;
   for (i=0; i < nkeys; i++)
      l_keys[i] = (ftnword) keys[i];
   for (i=0; i < ninfo; i++)
      l_info[i] = (ftnword) info[i];
   ier = c_xdfupd(iun,buf+1,idtyp,l_keys,nkeys,l_info,ninfo);
   BUF_F;
#else
   ier = c_xdfupd(iun,(buffer_interface_ptr)buf,idtyp,keys,nkeys,info,ninfo);
#endif
   return((ftnword) ier);
}

/*splitpoint xdfuse */
/***************************************************************************** 
 *                              X D F U S E                                  * 
 *****************************************************************************/

ftnword f77name(xdfuse)(ftnword *fsrc_unit, ftnword *fdest_unit)
{
   int src_unit = *fsrc_unit, dest_unit = *fdest_unit;

   return((ftnword)c_xdfuse(src_unit,dest_unit));

}

/*splitpoint xdfxtr */
/***************************************************************************** 
 *                              X D F X T R                                  * 
 *****************************************************************************/

ftnword f77name(xdfxtr)(word *buf,word *donnees,
                        ftnword *fbitpos, ftnword *fnelm,
			ftnword *fnbits, ftnword *fdatyp)
{
   int nelm = *fnelm, nbits = *fnbits, datyp = *fdatyp, bitpos = *fbitpos;
   int ier;

#if defined(NEC64)
   BUF_C;
   if ((datyp == 2) || (datyp == 4)) {
     xdf_stride = 2;
     ier = c_xdfxtr(buf+1,donnees+1,bitpos,nelm,nbits,datyp);
     xdf_stride = 1;
   }
   else
     ier = c_xdfxtr(buf+1,donnees,bitpos,nelm,nbits,datyp);
   BUF_F;
#else
   ier = c_xdfxtr(buf,donnees,bitpos,nelm,nbits,datyp);
#endif
   return((ftnword) ier);
}


/***************************************************************************** 
 *                            S E C A T E U R                                *
 *                                                                           * 
 *Object                                                                     * 
 *  The file whose name is given by filename has its size truncated to       *
 *  the number of bytes given by where.                                      *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  filename file name                                                   * 
 *  IN  where    number of bytes for the file zise                           * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(secateur)(char *filename, ftnword *f_where, F2Cl l1)
{
  int ier, where = *f_where;
  
  ier = c_secateur(filename,where);
  return((ftnword) ier);
}


/***************************************************************************** 
 *                        C _ S E C A T E U R                                *
 *                                                                           * 
 *Object                                                                     * 
 *  The file whose name is given by filename has its size truncated to       *
 *  the number of bytes given by where.                                      *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  filename file name                                                   * 
 *  IN  where    number of bytes for the file zise                           * 
 *                                                                           * 
 *****************************************************************************/

int c_secateur(char *filename, int where)
{
  int ier;
  
  if (msg_level <= TRIVIAL)
    fprintf(stdout,"Truncating %s to \t %d Bytes\n",filename,where);

  ier = truncate(filename,where);
  if (ier == -1) perror("secateur");
  return(ier);
}
