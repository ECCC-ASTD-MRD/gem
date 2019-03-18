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


/*splitpoint backto64 */
/***************************************************************************** 
 *                          B A C K T O 6 4                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Copy back content of field to a Fortran 64 bit array of elements.       *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN/OUT field    Fortran 64 bit array                                     * 
 *  IN     temp     c 32 bit array (same array, different pointer type)      * 
 *  IN     nelm     number of elements                                       * 
 *                                                                           * 
 *****************************************************************************/

void backto64(ftnword *field, word *temp, int nelm)
{
#if defined(NEC64)
  int i;
  INT_64 *pll = (INT_64 *) field;
  INT_32 *pl = (INT_32 *) temp;
  unsigned INT_64 *upll = (unsigned INT_64 *) field;
  unsigned INT_32 *upl = (unsigned INT_32 *) temp;

  if (xdf_datatyp == 4) {
      for (i=nelm-1; i>=0; i--) 
        pll[i] = pl[i];
  }
  else
    if ((xdf_datatyp == 2) || (xdf_datatyp == 3)) {
      for (i=nelm-1; i>=0; i--) 
        upll[i] = upl[i];
    }
#else
    exit(1);
#endif
}

/*splitpoint fstapp */
/***************************************************************************** 
 *                             F S T A P P                                   *
 *                                                                           * 
 *Object                                                                     * 
 *   Position at the end of a sequential file for an append.                 *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *  IN  option  kept for backward compatibility (not used)                   * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstapp)(ftnword *f_iun, char *option, F2Cl lng)
{
  int ier, iun = *f_iun;
  
  ier = c_fstapp(iun,option);    /* option not used anymore by c_fstapp */
  return((ftnword) ier);
}

/*splitpoint fstckp */
/***************************************************************************** 
 *                                F S T C K P                                *
 *                                                                           * 
 *Object                                                                     * 
 *   Checkpoint. Clear buffers, rewrite headers.                             *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *                                                                           * 
 *****************************************************************************/
ftnword f77name(fstckp)(ftnword *f_iun)
{
  int iun = *f_iun;
  int ier;

  ier = c_fstckp(iun);
  return((ftnword) ier);
}

/*splitpoint fstcvt */
/***************************************************************************** 
 *                                F S T C V T                                *
 *                                                                           * 
 *Object                                                                     * 
 *         VARIABLE UTILISEE COMME HOLLERITH SERA TRANSFORME EN CARACTERE    *
 *         OU L'INVERSE POUR NOMVAR,TYPVAR,GRID TYPE, ET ETIKET LE MOT       *
 *         ETIKET AURA 2 LOCATIONS SUR UNE MACHINE A 32 BITS                 *
 *                                                                           *
 *ARGUMENTS                                                                  *
 *  IN OUT    NOM       HOLLERITH *2                       [NOMVAR]          *
 *  IN OUT    TYP       HOLLERITH *1                       [TYPVAR]          * 
 *  IN OUT    ETIK      HOLLERITH *8   2 MOTS A4 POUR SUN  [ETIKET]          *
 *  IN OUT    GRTP      HOLLERITH *1                       [GRTYP]           *
 *  OUT IN    CNOM      CARACTERE *2                                         *
 *  OUT IN    CTYP      CARACTERE *1                                         *
 *  OUT IN    CETIK     CARACTERE *8                                         *
 *  OUT IN    CGRTP     CARACTRE *1                                          *
 *  IN        HOLACAR   LOGICAL .TRUE.  HOLLERITH A CARATERE                 *
 *                      LOGICAL .FALSE. CARACTERE A HOLLERITH                *
 *                                                                           * 
 *****************************************************************************/
ftnword f77name(fstcvt)(ftnword *name, ftnword *type, ftnword *etik, ftnword *grtyp,
                        char *cname, char *ctype, char *cetik, char *cgrtyp, ftnword *holocar,
                        F2Cl l1, F2Cl l2, F2Cl l3, F2Cl l4)
{
  int ier;

  ier = f77name(fstcvt2)(name,type,etik,grtyp,cname,ctype,cetik,cgrtyp,holocar,l1,l2,l3,l4);
  return((ftnword) ier);
}

/*splitpoint fst_data_length */
/***************************************************************************** 
 *                      F S T _ D A T A _ L E N G T H                        *
 *                                                                           * 
 *Object                                                                     * 
 *   Gives information on data lenght of the elements passed to fstecr       *
 *   and fstlir (double, short integer, byte ...)                            *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  length_type     data length kind                                     * 
 *                      1: byte                                              *
 *                      2: short (16 bits)                                   *
 *                      4: regular 32 bits                                   *
 *                      8: double (64 bits)                                  *
 *                                                                           * 
 *****************************************************************************/
ftnword f77name(fst_data_length)(int *f_length_type)
{
  int ier, length_type=*f_length_type;
  
  ier = c_fst_data_length(length_type);
  return((ftnword) ier);
}

/*splitpoint fstecr */
/***************************************************************************** 
 *                              F S T E C R                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Writes record to file.                                                  *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  field   field to write to the file                                   * 
 *  IN  work    work field (kept for backward compatibility)                 * 
 *  IN  npak    number of bits kept for the elements of the field (-npak)    * 
 *  IN  iun     unit number associated to the file                           * 
 *  IN  date    date time stamp                                              * 
 *  IN  deet    length of a time step in seconds                             * 
 *  IN  npas    time step number                                             * 
 *  IN  ni      first dimension of the data field                            * 
 *  IN  nj      second dimension of the data field                           * 
 *  IN  nk      third dimension of the data field                            * 
 *  IN  ip1     vertical level                                               * 
 *  IN  ip2     forecast hour                                                * 
 *  IN  ip3     user defined identifier                                      * 
 *  IN  typvar  type of field (forecast, analysis, climatology)              * 
 *  IN  nomvar  variable name                                                * 
 *  IN  etiket  label                                                        * 
 *  IN  grtyp   type of geographical projection                              * 
 *  IN  ig1     first grid descriptor                                        * 
 *  IN  ig2     second grid descriptor                                       * 
 *  IN  ig3     third grid descriptor                                        * 
 *  IN  ig4     fourth grid descriptor                                       * 
 *  IN  datyp   data type of the elements                                    * 
 *  IN  rewrit  rewrite flag (true=rewrite existing record, false=append)    *
 *                                                                           * 
 *****************************************************************************/
ftnword f77name(fstecr)(word *field, ftnword *work, ftnword *f_npak,
                        ftnword *f_iun, ftnword *f_date,
                        ftnword *f_deet, ftnword *f_npas,
                        ftnword *f_ni, ftnword *f_nj, ftnword *f_nk,
                        ftnword *f_ip1, ftnword *f_ip2, ftnword *f_ip3,
                        char *f_typvar, char *f_nomvar, char *f_etiket,
                        char *f_grtyp, ftnword *f_ig1, ftnword *f_ig2,
                        ftnword *f_ig3, ftnword *f_ig4,
                        ftnword *f_datyp, ftnword *f_rewrit,
                        F2Cl ll1, F2Cl ll2, F2Cl ll3, F2Cl ll4)
{
  int iun=*f_iun, npak=*f_npak, date=*f_date, deet=*f_deet;
  int npas=*f_npas, ip1=*f_ip1, ip2=*f_ip2, ip3=*f_ip3;
  int ni=*f_ni, nj=*f_nj, nk=*f_nk;
  int ig1=*f_ig1, ig2=*f_ig2, ig3=*f_ig3, ig4=*f_ig4;
  int datyp=*f_datyp, rewrit=*f_rewrit;
  int ier;
  int l1=ll1, l2=ll2, l3=ll3, l4=ll4;

  char etiket[13];
  char typvar[3];
  char nomvar[5];
  char grtyp[2];

/*  l1 = (l1 < 2) ? l1 : 2;            /* typvar length */
/*  l2 = (l2 < 4) ? l2 : 4;            /* nomvar length */
/*  l3 = (l3 < 12) ? l3 :12;           /* etiket length */
/*  l4 = (l4 < 1) ? l4 : 1;            /* grtyp length */

  str_cp_init(typvar,3,f_typvar,l1);
  str_cp_init(nomvar,5,f_nomvar,l2);
  str_cp_init(etiket,13,f_etiket,l3);
  str_cp_init(grtyp,2,f_grtyp,l4);

#if defined (NEC64)
  if ((datyp == 1) || (datyp == 5)) {      /* floating point */
    xdf_double = 1;
    ier = c_fstecr(field,work,npak,iun,date,deet,npas,
                   ni,nj,nk,ip1,ip2,ip3,
                   typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,
                   datyp,rewrit);
    xdf_double = 0;
  }
  else if ((datyp != 0) && (!image_mode_copy)) {
    xdf_stride = 2;
    ier = c_fstecr(field+1,work,npak,iun,date,deet,npas,
                   ni,nj,nk,ip1,ip2,ip3,
                   typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,
                   datyp,rewrit);
    xdf_stride = 1;
  }
  else {
  ier = c_fstecr(field,work,npak,iun,date,deet,npas,
                 ni,nj,nk,ip1,ip2,ip3,
                 typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,
                 datyp,rewrit);
  }
#else
  ier = c_fstecr(field,work,npak,iun,date,deet,npas,
                 ni,nj,nk,ip1,ip2,ip3,
                 typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,
                 datyp,rewrit);
#endif
  return((ftnword) ier);
}
/*splitpoint fstecr_s */
/***************************************************************************** 
 *                              F S T E C R _ S                              *
 *                                                                           * 
 *Object                                                                     * 
 *   Writes record to file.                                                  *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  string  character string to write to the file                        * 
 *  IN  work    work field (kept for backward compatibility)                 * 
 *  IN  npak    number of bits kept for the elements of the field (-npak)    * 
 *  IN  iun     unit number associated to the file                           * 
 *  IN  date    date time stamp                                              * 
 *  IN  deet    length of a time step in seconds                             * 
 *  IN  npas    time step number                                             * 
 *  IN  ni      first dimension of the data field                            * 
 *  IN  nj      second dimension of the data field                           * 
 *  IN  nk      third dimension of the data field                            * 
 *  IN  ip1     vertical level                                               * 
 *  IN  ip2     forecast hour                                                * 
 *  IN  ip3     user defined identifier                                      * 
 *  IN  typvar  type of field (forecast, analysis, climatology)              * 
 *  IN  nomvar  variable name                                                * 
 *  IN  etiket  label                                                        * 
 *  IN  grtyp   type of geographical projection                              * 
 *  IN  ig1     first grid descriptor                                        * 
 *  IN  ig2     second grid descriptor                                       * 
 *  IN  ig3     third grid descriptor                                        * 
 *  IN  ig4     fourth grid descriptor                                       * 
 *  IN  datyp   data type of the elements                                    * 
 *  IN  rewrit  rewrite flag (true=rewrite existing record, false=append)    *
 *                                                                           * 
 *****************************************************************************/
ftnword f77name(fstecr_s)(void *string, ftnword *work, ftnword *f_npak,
                        ftnword *f_iun, ftnword *f_date,
                        ftnword *f_deet, ftnword *f_npas,
                        ftnword *f_ni, ftnword *f_nj, ftnword *f_nk,
                        ftnword *f_ip1, ftnword *f_ip2, ftnword *f_ip3,
                        char *f_typvar, char *f_nomvar, char *f_etiket,
                        char *f_grtyp, ftnword *f_ig1, ftnword *f_ig2,
                        ftnword *f_ig3, ftnword *f_ig4,
                        ftnword *f_datyp, ftnword *f_rewrit,
                        int lng_string, F2Cl ll1, F2Cl ll2, F2Cl ll3, F2Cl ll4)
{
int ninjnk;
ftnword ier=0;

ninjnk = Max(1,*f_ni) * Max(1,*f_nj) * Max(1,*f_nk);
if (ninjnk > lng_string * *f_nj) {
  sprintf(errmsg,"ni*nj*nk (%d) > string length (%d)",ninjnk,lng_string);
  return(error_msg("FSTECR_S",ERR_BAD_DIM,ERROR));
  }
else
  {
    ier = f77name(fstecr)(string,work,f_npak,f_iun,f_date,f_deet,f_npas,f_ni,f_nj,f_nk,f_ip1,f_ip2,f_ip3,
                          f_typvar,f_nomvar,f_etiket,f_grtyp,f_ig1,f_ig2,f_ig3,f_ig4,f_datyp,f_rewrit,
                          ll1,ll2,ll3,ll4);
    return(ier);
  }
}                        
/*splitpoint fstecr_h */
/***************************************************************************** 
 *                              F S T E C R _ H                              *
 *                                                                           * 
 *Object                                                                     * 
 *   Writes record to file.                                                  *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  haft_w  haft word (16 bit) array to write to the file                * 
 *  IN  work    work field (kept for backward compatibility)                 * 
 *  IN  npak    number of bits kept for the elements of the field (-npak)    * 
 *  IN  iun     unit number associated to the file                           * 
 *  IN  date    date time stamp                                              * 
 *  IN  deet    length of a time step in seconds                             * 
 *  IN  npas    time step number                                             * 
 *  IN  ni      first dimension of the data field                            * 
 *  IN  nj      second dimension of the data field                           * 
 *  IN  nk      third dimension of the data field                            * 
 *  IN  ip1     vertical level                                               * 
 *  IN  ip2     forecast hour                                                * 
 *  IN  ip3     user defined identifier                                      * 
 *  IN  typvar  type of field (forecast, analysis, climatology)              * 
 *  IN  nomvar  variable name                                                * 
 *  IN  etiket  label                                                        * 
 *  IN  grtyp   type of geographical projection                              * 
 *  IN  ig1     first grid descriptor                                        * 
 *  IN  ig2     second grid descriptor                                       * 
 *  IN  ig3     third grid descriptor                                        * 
 *  IN  ig4     fourth grid descriptor                                       * 
 *  IN  datyp   data type of the elements                                    * 
 *  IN  rewrit  rewrite flag (true=rewrite existing record, false=append)    *
 *                                                                           * 
 *****************************************************************************/
ftnword f77name(fstecr_h)(void *haft_w, ftnword *work, ftnword *f_npak,
                        ftnword *f_iun, ftnword *f_date,
                        ftnword *f_deet, ftnword *f_npas,
                        ftnword *f_ni, ftnword *f_nj, ftnword *f_nk,
                        ftnword *f_ip1, ftnword *f_ip2, ftnword *f_ip3,
                        char *f_typvar, char *f_nomvar, char *f_etiket,
                        char *f_grtyp, ftnword *f_ig1, ftnword *f_ig2,
                        ftnword *f_ig3, ftnword *f_ig4,
                        ftnword *f_datyp, ftnword *f_rewrit,
                        F2Cl ll1, F2Cl ll2, F2Cl ll3, F2Cl ll4)
{
// int ninjnk;
ftnword ier=0;
  xdf_short = 1;
  ier = f77name(fstecr)(haft_w,work,f_npak,f_iun,f_date,f_deet,f_npas,f_ni,f_nj,f_nk,f_ip1,f_ip2,f_ip3,
                        f_typvar,f_nomvar,f_etiket,f_grtyp,f_ig1,f_ig2,f_ig3,f_ig4,f_datyp,f_rewrit,
                        ll1,ll2,ll3,ll4);
  xdf_short = 0;
  return(ier);
}                        
/*splitpoint fstecr_b */
/***************************************************************************** 
 *                              F S T E C R _ B                              *
 *                                                                           * 
 *Object                                                                     * 
 *   Writes record to file.                                                  *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  bytes   byte array to write to the file                              * 
 *  IN  work    work field (kept for backward compatibility)                 * 
 *  IN  npak    number of bits kept for the elements of the field (-npak)    * 
 *  IN  iun     unit number associated to the file                           * 
 *  IN  date    date time stamp                                              * 
 *  IN  deet    length of a time step in seconds                             * 
 *  IN  npas    time step number                                             * 
 *  IN  ni      first dimension of the data field                            * 
 *  IN  nj      second dimension of the data field                           * 
 *  IN  nk      third dimension of the data field                            * 
 *  IN  ip1     vertical level                                               * 
 *  IN  ip2     forecast hour                                                * 
 *  IN  ip3     user defined identifier                                      * 
 *  IN  typvar  type of field (forecast, analysis, climatology)              * 
 *  IN  nomvar  variable name                                                * 
 *  IN  etiket  label                                                        * 
 *  IN  grtyp   type of geographical projection                              * 
 *  IN  ig1     first grid descriptor                                        * 
 *  IN  ig2     second grid descriptor                                       * 
 *  IN  ig3     third grid descriptor                                        * 
 *  IN  ig4     fourth grid descriptor                                       * 
 *  IN  datyp   data type of the elements                                    * 
 *  IN  rewrit  rewrite flag (true=rewrite existing record, false=append)    *
 *                                                                           * 
 *****************************************************************************/
ftnword f77name(fstecr_b)(void *bytes, ftnword *work, ftnword *f_npak,
                        ftnword *f_iun, ftnword *f_date,
                        ftnword *f_deet, ftnword *f_npas,
                        ftnword *f_ni, ftnword *f_nj, ftnword *f_nk,
                        ftnword *f_ip1, ftnword *f_ip2, ftnword *f_ip3,
                        char *f_typvar, char *f_nomvar, char *f_etiket,
                        char *f_grtyp, ftnword *f_ig1, ftnword *f_ig2,
                        ftnword *f_ig3, ftnword *f_ig4,
                        ftnword *f_datyp, ftnword *f_rewrit,
                        F2Cl ll1, F2Cl ll2, F2Cl ll3, F2Cl ll4)
{
//int ninjnk;
ftnword ier=0;
  xdf_byte = 1;
  ier = f77name(fstecr)(bytes,work,f_npak,f_iun,f_date,f_deet,f_npas,f_ni,f_nj,f_nk,f_ip1,f_ip2,f_ip3,
                        f_typvar,f_nomvar,f_etiket,f_grtyp,f_ig1,f_ig2,f_ig3,f_ig4,f_datyp,f_rewrit,
                        ll1,ll2,ll3,ll4);
  xdf_byte = 0;
  return(ier);
}                        

/*splitpoint c_fst_edit_dir */
/*****************************************************************************
 *                     C _ F S T _ E D I T _ D I R                           *
 *                                                                           *
 *Object                                                                     * 
 *   Edits the directory content of a RPN standard file.                     *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  handle     handle to the directory entry to edit                     * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77_name(fst_edit_dir_plus)(ftnword *f_handle,
                               ftnword *f_date, ftnword *f_deet, ftnword *f_npas,
                               ftnword *f_ni, ftnword *f_nj, ftnword *f_nk,
                               ftnword *f_ip1, ftnword *f_ip2, ftnword *f_ip3,
                               char *f_typvar, char *f_nomvar, char *f_etiket,
                               char *f_grtyp, ftnword *f_ig1, ftnword *f_ig2,
                               ftnword *f_ig3, ftnword *f_ig4, ftnword *f_datyp,
                               F2Cl l1, F2Cl l2, F2Cl l3, F2Cl l4)
{
  int ier;
  int handle=*f_handle;
  int date=*f_date, deet=*f_deet;
  int npas=*f_npas, ip1=*f_ip1, ip2=*f_ip2, ip3=*f_ip3;
  int ni=*f_ni, nj=*f_nj, nk=*f_nk;
  int ig1=*f_ig1, ig2=*f_ig2, ig3=*f_ig3, ig4=*f_ig4;
  int datyp=*f_datyp;

  char etiket[13];
  char typvar[3];
  char nomvar[5];
  char grtyp[2];

  str_cp_init(typvar,3,f_typvar,l1);
  str_cp_init(nomvar,5,f_nomvar,l2);
  str_cp_init(etiket,13,f_etiket,l3);
  str_cp_init(grtyp,2,f_grtyp,l4);

  ier = c_fst_edit_dir_plus(handle,date,deet,npas,ni,nj,nk,ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,
                       ig1,ig2,ig3,ig4,datyp);
  return((ftnword) ier);
}

ftnword f77_name(fst_edit_dir)(ftnword *f_handle,
                               ftnword *f_date, ftnword *f_deet, ftnword *f_npas,
                               ftnword *f_ni, ftnword *f_nj, ftnword *f_nk,
                               ftnword *f_ip1, ftnword *f_ip2, ftnword *f_ip3,
                               char *f_typvar, char *f_nomvar, char *f_etiket,
                               char *f_grtyp, ftnword *f_ig1, ftnword *f_ig2,
                               ftnword *f_ig3, ftnword *f_ig4, ftnword *f_datyp,
                               F2Cl l1, F2Cl l2, F2Cl l3, F2Cl l4)
{
  int ier;
  int handle=*f_handle;
  int date=*f_date, deet=*f_deet;
  int npas=*f_npas, ip1=*f_ip1, ip2=*f_ip2, ip3=*f_ip3;
  int ni=*f_ni, nj=*f_nj, nk=*f_nk;
  int ig1=*f_ig1, ig2=*f_ig2, ig3=*f_ig3, ig4=*f_ig4;
  int datyp=*f_datyp;

  char etiket[13];
  char typvar[3];
  char nomvar[5];
  char grtyp[2];

  str_cp_init(typvar,3,f_typvar,l1);
  str_cp_init(nomvar,5,f_nomvar,l2);
  str_cp_init(etiket,13,f_etiket,l3);
  str_cp_init(grtyp,2,f_grtyp,l4);

  ier = c_fst_edit_dir(handle,date,deet,npas,ni,nj,nk,ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,
                       ig1,ig2,ig3,ig4,datyp);
  return((ftnword) ier);
}

/*splitpoint fsteff */
/***************************************************************************** 
 *                             F S T E F F                                   *
 *                                                                           * 
 *Object                                                                     * 
 *   Deletes the record associated to handle.                                *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  handle  handle to the record to delete                               * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fsteff)(ftnword *f_handle)
{
  int ier,handle = *f_handle;

  ier = c_fsteff(handle);
  return((ftnword) ier);
}

/*splitpoint fsteof */
/***************************************************************************** 
 *                             F S T E O F                                   *
 *                                                                           * 
 *Object                                                                     * 
 *   Return the level of end of file for the sequential file.                *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fsteof)(ftnword *f_iun)
{
  int eof, iun = *f_iun;
  
  eof = c_fsteof(iun);
  return((ftnword) eof);
}

/*splitpoint fstfrm */
/***************************************************************************** 
 *                              F S T F R M                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Closes a RPN standard file.                                             *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstfrm)(ftnword *f_iun)
{
  int iun = *f_iun;
  return ((ftnword) c_fstfrm(iun));
}

/*splitpoint fstinf */
/***************************************************************************** 
 *                              F S T I N F                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Locate the next record that matches the research keys.                  *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *  OUT ni      dimension 1 of the data field                                * 
 *  OUT nj      dimension 2 of the data field                                * 
 *  OUT nk      dimension 3 of the data field                                * 
 *  IN  datev   valid date                                                   * 
 *  IN  etiket  label                                                        * 
 *  IN  ip1     vertical level                                               * 
 *  IN  ip2     forecast hour                                                * 
 *  IN  ip3     user defined identifier                                      * 
 *  IN  typvar  type of field                                                * 
 *  IN  nomvar  variable name                                                * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstinf)(ftnword *f_iun, ftnword *f_ni, ftnword *f_nj,
                        ftnword *f_nk, ftnword *f_datev, char *f_etiket,
                        ftnword *f_ip1, ftnword *f_ip2, ftnword *f_ip3,
                        char *f_typvar, char *f_nomvar,
                        F2Cl ll1, F2Cl ll2, F2Cl ll3)
{
  int iun = *f_iun, datev=*f_datev, ip1=*f_ip1, ip2=*f_ip2, ip3=*f_ip3;
  int ier,ni,nj,nk;
  int l1=ll1, l2=ll2, l3=ll3;

  char etiket[13];
  char typvar[3];
  char nomvar[5];

/*  l1 = (l1 < 12) ? l1 :12;           /* etiket length */
/*  l2 = (l2 < 2) ? l2 : 2;            /* typvar length */
/*  l3 = (l3 < 4) ? l3 : 4;            /* nomvar length */

  str_cp_init(etiket,13,f_etiket,l1);
  str_cp_init(typvar,3,f_typvar,l2);
  str_cp_init(nomvar,5,f_nomvar,l3);

  ier = c_fstinf(iun,&ni,&nj,&nk,datev,etiket,ip1,ip2,ip3,typvar,nomvar);
  *f_ni = (ftnword) ni;
  *f_nj = (ftnword) nj;
  *f_nk = (ftnword) nk;
  return((ftnword) ier);
}

/*splitpoint fstinfx */
/***************************************************************************** 
 *                              F S T I N F X                                *
 *                                                                           * 
 *Object                                                                     * 
 *   Locate the next record that matches the research keys.                  *
 *   The search begins at the position given by handle.                      * 
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *  OUT ni      dimension 1 of the data field                                * 
 *  OUT nj      dimension 2 of the data field                                * 
 *  OUT nk      dimension 3 of the data field                                * 
 *  IN  datev   valid date                                                   * 
 *  IN  etiket  label                                                        * 
 *  IN  ip1     vertical level                                               * 
 *  IN  ip2     forecast hour                                                * 
 *  IN  ip3     user defined identifier                                      * 
 *  IN  typvar  type of field                                                * 
 *  IN  nomvar  variable name                                                * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstinfx)(ftnword *f_handle, ftnword *f_iun,
                         ftnword *f_ni, ftnword *f_nj,
                         ftnword *f_nk, ftnword *f_datev, char *f_etiket,
                         ftnword *f_ip1, ftnword *f_ip2, ftnword *f_ip3,
                         char *f_typvar, char *f_nomvar,
                         F2Cl ll1, F2Cl ll2, F2Cl ll3)
{
  int iun = *f_iun, datev=*f_datev, ip1=*f_ip1, ip2=*f_ip2, ip3=*f_ip3;
  int handle = *f_handle;
  int ier,ni,nj,nk;
  int l1=ll1, l2=ll2, l3=ll3;

  char etiket[13];
  char typvar[3];
  char nomvar[5];

  str_cp_init(etiket,13,f_etiket,l1);
  str_cp_init(typvar,3,f_typvar,l2);
  str_cp_init(nomvar,5,f_nomvar,l3);

  ier = c_fstinfx(handle,iun,&ni,&nj,&nk,datev,etiket,
                     ip1,ip2,ip3,typvar,nomvar);
  *f_ni = (ftnword) ni;
  *f_nj = (ftnword) nj;
  *f_nk = (ftnword) nk;
  return((ftnword) ier);
}

/*splitpoint fstinl */
/***************************************************************************** 
 *                              F S T I N L                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Locates all the records that matches the research keys.                 *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *  OUT ni      dimension 1 of the data field                                * 
 *  OUT nj      dimension 2 of the data field                                * 
 *  OUT nk      dimension 3 of the data field                                * 
 *  IN  datev   valid date                                                   * 
 *  IN  etiket  label                                                        * 
 *  IN  ip1     vertical level                                               * 
 *  IN  ip2     forecast hour                                                * 
 *  IN  ip3     user defined identifier                                      * 
 *  IN  typvar  type of field                                                * 
 *  IN  nomvar  variable name                                                *
 *  OUT liste   list of handle to the records                                *
 *  OUT infon   number of elements for the list (number of records found)    *
 *  OUT nmax    dimension of list as given by caller                         *
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstinl)(ftnword *f_iun, ftnword *f_ni, ftnword *f_nj,
                        ftnword *f_nk, ftnword *f_datev, char *f_etiket,
                        ftnword *f_ip1, ftnword *f_ip2, ftnword *f_ip3,
                        char *f_typvar, char *f_nomvar,
                        ftnword *liste, ftnword *f_infon, ftnword *f_nmax,
                        F2Cl ll1, F2Cl ll2, F2Cl ll3)
{
  int iun = *f_iun, datev=*f_datev, ip1=*f_ip1, ip2=*f_ip2, ip3=*f_ip3;
  int infon, nmax = *f_nmax;
  int ier,ni,nj,nk, i;
  int l1=ll1, l2=ll2, l3=ll3;
  INT_32 *plong;
  char etiket[13];
  char typvar[3];
  char nomvar[5];

  str_cp_init(etiket,13,f_etiket,l1);
  str_cp_init(typvar,3,f_typvar,l2);
  str_cp_init(nomvar,5,f_nomvar,l3);

  ier = c_fstinl(iun,&ni,&nj,&nk,datev,etiket,ip1,ip2,ip3,typvar,nomvar,
                 (word *)liste,&infon,nmax);
  *f_ni = (ftnword) ni;
  *f_nj = (ftnword) nj;
  *f_nk = (ftnword) nk;
  *f_infon = (ftnword) infon;
  return((ftnword) ier);
}


/*splitpoint fstlic */
/***************************************************************************** 
 *                             F S T L I C                                   *
 *                                                                           * 
 *Object                                                                     * 
 *   Search for a record that matches the research keys and check that the   *
 *   remaining parmeters match the record descriptors                        *
 *                                                                           *
 *Arguments                                                                  * 
 *                                                                           * 
 *  OUT field    data field to be read                                       * 
 *  IN  iun      unit number associated to the file                          * 
 *  IN  niin     dimension 1 of the data field                               * 
 *  IN  njin     dimension 2 of the data field                               * 
 *  IN  nkin     dimension 3 of the data field                               * 
 *  IN  datein   valid date                                                  * 
 *  IN  etiketin label                                                       * 
 *  IN  ip1in    vertical level                                              * 
 *  IN  ip2in    forecast hour                                               * 
 *  IN  ip3in    user defined identifier                                     * 
 *  IN  typvarin type of field                                               * 
 *  IN  nomvarin variable name                                               * 
 *  IN  ig1      first grid descriptor                                       * 
 *  IN  ig2      second grid descriptor                                      * 
 *  IN  ig3      third grid descriptor                                       * 
 *  IN  ig4      fourth grid descriptor                                      * 
 *  IN  grtypin  type of geographical projection                             * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstlic)(word *field, ftnword *f_iun,
                         ftnword *f_ni, ftnword *f_nj,
                         ftnword *f_nk, ftnword *f_date, char *f_etiket,
                         ftnword *f_ip1, ftnword *f_ip2, ftnword *f_ip3,
                         char *f_typvar, char *f_nomvar,
                         ftnword *f_ig1, ftnword *f_ig2, ftnword *f_ig3,
                         ftnword *f_ig4, char *f_grtyp, 
                         F2Cl ll1, F2Cl ll2, F2Cl ll3, F2Cl ll4)
{

  int iun = *f_iun, ni = *f_ni, nj = *f_nj, nk = *f_nk;
  int date = *f_date, ip1 = *f_ip1, ip2 = *f_ip2, ip3 = *f_ip3;
  int ig1 = *f_ig1, ig2 = *f_ig2, ig3 = *f_ig3, ig4 = *f_ig4;
  int ier;
  int l1=ll1, l2=ll2, l3=ll3, l4=ll4;

  char etiket[13];
  char typvar[3];
  char nomvar[5];
  char grtyp[2];

  str_cp_init(etiket,13,f_etiket,l1);
  str_cp_init(typvar,3,f_typvar,l2);
  str_cp_init(nomvar,5,f_nomvar,l3);
  str_cp_init(grtyp,2,f_grtyp,l4);

  ier = c_fstlic(field,iun,ni,nj,nk,date,etiket,ip1,ip2,ip3,
                     typvar,nomvar,ig1,ig2,ig3,ig4,grtyp);
  return((ftnword) ier);
}

/*splitpoint fstlir */
/***************************************************************************** 
 *                              F S T L I R                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Reads the next record that matches the research keys.                   *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *  OUT field   data field to be read                                        * 
 *  OUT ni      dimension 1 of the data field                                * 
 *  OUT nj      dimension 2 of the data field                                * 
 *  OUT nk      dimension 3 of the data field                                * 
 *  IN  datev   valid date                                                   * 
 *  IN  etiket  label                                                        * 
 *  IN  ip1     vertical level                                               * 
 *  IN  ip2     forecast hour                                                * 
 *  IN  ip3     user defined identifier                                      * 
 *  IN  typvar  type of field                                                * 
 *  IN  nomvar  variable name                                                * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstlir)(void *field, ftnword *f_iun,
                        ftnword *f_ni, ftnword *f_nj,
                        ftnword *f_nk, ftnword *f_datev, char *f_etiket,
                        ftnword *f_ip1, ftnword *f_ip2, ftnword *f_ip3,
                        char *f_typvar, char *f_nomvar,
                        F2Cl ll1, F2Cl ll2, F2Cl ll3)
{
  int iun = *f_iun, datev=*f_datev, ip1=*f_ip1, ip2=*f_ip2, ip3=*f_ip3;
  int ier,ni,nj,nk;
  int l1=ll1, l2=ll2, l3=ll3;

  char etiket[13];
  char typvar[3];
  char nomvar[5];

  str_cp_init(etiket,13,f_etiket,l1);
  str_cp_init(typvar,3,f_typvar,l2);
  str_cp_init(nomvar,5,f_nomvar,l3);

#if defined(NEC64)
  xdf_double = 1;
  ier = c_fstlir(field,iun,&ni,&nj,&nk,datev,etiket,
                     ip1,ip2,ip3,typvar,nomvar);
  backto64(field,(word *) field,ni*nj*nk);
  xdf_double = 0;
#else
  ier = c_fstlir(field,iun,&ni,&nj,&nk,datev,etiket,
                     ip1,ip2,ip3,typvar,nomvar);
#endif
  if (ier < 0) return((ftnword) ier);
  *f_ni = (ftnword) ni;
  *f_nj = (ftnword) nj;
  *f_nk = (ftnword) nk;
  return((ftnword) ier);
}
/*splitpoint fstlir_s */
/***************************************************************************** 
 *                              F S T L I R _ S                              *
 *                                                                           * 
 *Object                                                                     * 
 *   Reads the next record that matches the research keys.                   *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *  OUT string  character string to be read                                  * 
 *  OUT ni      dimension 1 of the data field                                * 
 *  OUT nj      dimension 2 of the data field                                * 
 *  OUT nk      dimension 3 of the data field                                * 
 *  IN  datev   valid date                                                   * 
 *  IN  etiket  label                                                        * 
 *  IN  ip1     vertical level                                               * 
 *  IN  ip2     forecast hour                                                * 
 *  IN  ip3     user defined identifier                                      * 
 *  IN  typvar  type of field                                                * 
 *  IN  nomvar  variable name                                                * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstlir_s)(void *string, ftnword *f_iun,
                        ftnword *f_ni, ftnword *f_nj,
                        ftnword *f_nk, ftnword *f_datev, char *f_etiket,
                        ftnword *f_ip1, ftnword *f_ip2, ftnword *f_ip3,
                        char *f_typvar, char *f_nomvar,
                        int lng_string, F2Cl ll1, F2Cl ll2, F2Cl ll3)
{
  ftnword ier;
  int i;
  char *ptr=(char*) string;
  
  for (i=0; i < lng_string; i++)
    *ptr++ = ' ';
  ier = f77name(fstlir)(string,f_iun,f_ni,f_nj,f_nk,f_datev,f_etiket,
                     f_ip1,f_ip2,f_ip3,f_typvar,f_nomvar,ll1,ll2,ll3);
  return(ier);
}
/*splitpoint fstlir_h */
/***************************************************************************** 
 *                              F S T L I R _ H                              *
 *                                                                           * 
 *Object                                                                     * 
 *   Reads the next record that matches the research keys.                   *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *  OUT haft_w  haft word short integer (16 bit) array to be read            * 
 *  OUT ni      dimension 1 of the data field                                * 
 *  OUT nj      dimension 2 of the data field                                * 
 *  OUT nk      dimension 3 of the data field                                * 
 *  IN  datev   valid date                                                   * 
 *  IN  etiket  label                                                        * 
 *  IN  ip1     vertical level                                               * 
 *  IN  ip2     forecast hour                                                * 
 *  IN  ip3     user defined identifier                                      * 
 *  IN  typvar  type of field                                                * 
 *  IN  nomvar  variable name                                                * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstlir_h)(void *haft_w, ftnword *f_iun,
                        ftnword *f_ni, ftnword *f_nj,
                        ftnword *f_nk, ftnword *f_datev, char *f_etiket,
                        ftnword *f_ip1, ftnword *f_ip2, ftnword *f_ip3,
                        char *f_typvar, char *f_nomvar,
                        F2Cl ll1, F2Cl ll2, F2Cl ll3)
{
  ftnword ier;
  
  xdf_short = 1;
  ier = f77name(fstlir)(haft_w,f_iun,f_ni,f_nj,f_nk,f_datev,f_etiket,
                     f_ip1,f_ip2,f_ip3,f_typvar,f_nomvar,ll1,ll2,ll3);
  xdf_short = 0;
  return(ier);
}
/*splitpoint fstlir_b */
/***************************************************************************** 
 *                              F S T L I R _ B                              *
 *                                                                           * 
 *Object                                                                     * 
 *   Reads the next record that matches the research keys.                   *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *  OUT bytes   byte array to be read                                        * 
 *  OUT ni      dimension 1 of the data field                                * 
 *  OUT nj      dimension 2 of the data field                                * 
 *  OUT nk      dimension 3 of the data field                                * 
 *  IN  datev   valid date                                                   * 
 *  IN  etiket  label                                                        * 
 *  IN  ip1     vertical level                                               * 
 *  IN  ip2     forecast hour                                                * 
 *  IN  ip3     user defined identifier                                      * 
 *  IN  typvar  type of field                                                * 
 *  IN  nomvar  variable name                                                * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstlir_b)(void *bytes, ftnword *f_iun,
                        ftnword *f_ni, ftnword *f_nj,
                        ftnword *f_nk, ftnword *f_datev, char *f_etiket,
                        ftnword *f_ip1, ftnword *f_ip2, ftnword *f_ip3,
                        char *f_typvar, char *f_nomvar,
                        F2Cl ll1, F2Cl ll2, F2Cl ll3)
{
  ftnword ier;
  
  xdf_byte = 1;
  ier = f77name(fstlir)(bytes,f_iun,f_ni,f_nj,f_nk,f_datev,f_etiket,
                     f_ip1,f_ip2,f_ip3,f_typvar,f_nomvar,ll1,ll2,ll3);
  xdf_byte = 0;
  return(ier);
}

/*splitpoint fstlirx */
/***************************************************************************** 
 *                              F S T L I R X                                *
 *                                                                           * 
 *Object                                                                     * 
 *   Reads the next record that matches the research keys.                   *
 *   The search begins at the position given by handle.                      * 
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *  OUT field   data field to be read                                        * 
 *  OUT ni      dimension 1 of the data field                                * 
 *  OUT nj      dimension 2 of the data field                                * 
 *  OUT nk      dimension 3 of the data field                                * 
 *  IN  datev   valid date                                                   * 
 *  IN  etiket  label                                                        * 
 *  IN  ip1     vertical level                                               * 
 *  IN  ip2     forecast hour                                                * 
 *  IN  ip3     user defined identifier                                      * 
 *  IN  typvar  type of field                                                * 
 *  IN  nomvar  variable name                                                * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstlirx)(word *field, ftnword *f_handle, ftnword *f_iun,
                        ftnword *f_ni, ftnword *f_nj,
                        ftnword *f_nk, ftnword *f_datev, char *f_etiket,
                        ftnword *f_ip1, ftnword *f_ip2, ftnword *f_ip3,
                        char *f_typvar, char *f_nomvar,
                        F2Cl ll1, F2Cl ll2, F2Cl ll3)
{
  int iun = *f_iun, datev=*f_datev, ip1=*f_ip1, ip2=*f_ip2, ip3=*f_ip3;
  int handle = *f_handle;
  int ier,ni,nj,nk;
  int l1=ll1, l2=ll2, l3=ll3;

  char etiket[13];
  char typvar[3];
  char nomvar[5];

  str_cp_init(etiket,13,f_etiket,l1);
  str_cp_init(typvar,3,f_typvar,l2);
  str_cp_init(nomvar,5,f_nomvar,l3);

#if defined(NEC64)
  xdf_double = 1;
  ier = c_fstlirx(field,handle,iun,&ni,&nj,&nk,datev,etiket,
                     ip1,ip2,ip3,typvar,nomvar);
  backto64(field,(word *) field,ni*nj*nk);
  xdf_double = 0;
#else
  ier = c_fstlirx(field,handle,iun,&ni,&nj,&nk,datev,etiket,
                     ip1,ip2,ip3,typvar,nomvar);
#endif
  *f_ni = (ftnword) ni;
  *f_nj = (ftnword) nj;
  *f_nk = (ftnword) nk;
  return((ftnword) ier);
}

/*splitpoint fstlis */
/***************************************************************************** 
 *                            F S T L I S                                    *
 *                                                                           * 
 *Object                                                                     * 
 *   Reads the next record that matches the last search criterias            *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  OUT field   data field to be read                                        * 
 *  IN  iun     unit number associated to the file                           * 
 *  OUT ni      dimension 1 of the data field                                * 
 *  OUT nj      dimension 2 of the data field                                * 
 *  OUT nk      dimension 3 of the data field                                * 
 *                                                                           * 
 *****************************************************************************/
ftnword f77name(fstlis)(word *field, ftnword *f_iun,
                        ftnword *f_ni, ftnword *f_nj, ftnword *f_nk)
{
  int iun = *f_iun;
  int ier,ni,nj,nk;

#if defined(NEC64)
  xdf_double = 1;
  ier = c_fstlis((word *)field,iun,&ni,&nj,&nk);
  backto64(field,(word *) field,ni*nj*nk);
  xdf_double = 0;
#else
  ier = c_fstlis((word *)field,iun,&ni,&nj,&nk);
#endif
  *f_ni = (ftnword) ni;
  *f_nj = (ftnword) nj;
  *f_nk = (ftnword) nk;
  return((ftnword) ier);
}


/*splitpoint fstlnk */
/***************************************************************************** 
 *                              F S T L N K                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Links a list of files together for search purpose.                      *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  liste   list of unit numbers associated to the files                 * 
 *  IN  n       number of files to link                                      * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstlnk)(ftnword *liste, ftnword *f_n)
{
  int ier;
  int i;

  link_n = *f_n;
  for (i=0; i<link_n; i++)
    link_list[i] = liste[i];
  ier = c_xdflnk(link_list,link_n);
  return ((ftnword) ier);
}

/*splitpoint fstluk */
/***************************************************************************** 
 *                              F S T L U K                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Read the record at position given by handle.                            *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  OUT field   data field to be read                                        * 
 *  IN  handle  positioning information to the record                        * 
 *  OUT ni      dimension 1 of the data field                                * 
 *  OUT nj      dimension 2 of the data field                                * 
 *  OUT nk      dimension 3 of the data field                                * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstluk)(word *field, ftnword *f_handle,
                        ftnword *f_ni, ftnword *f_nj, ftnword *f_nk)
{
  int handle = *f_handle;
  int ier,ni,nj,nk;

#if defined(NEC64)
  xdf_double = 1;
  ier = c_fstluk(field,handle,&ni,&nj,&nk);
  backto64(field,(word *) field,ni*nj*nk);
  xdf_double = 0;
#else
  ier = c_fstluk(field,handle,&ni,&nj,&nk);
#endif
  *f_ni = (ftnword) ni;
  *f_nj = (ftnword) nj;
  *f_nk = (ftnword) nk;
  return((ftnword) ier);
}

/*splitpoint fstmsq */
/***************************************************************************** 
 *                            F S T M S Q                                    *
 *                                                                           * 
 *Object                                                                     * 
 *   Mask a portion of the research keys.                                    *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *   IN    iun     unit number associated to the file                        * 
 * IN/OUT  mip1    mask for vertical level                                   * 
 * IN/OUT  mip2    mask for forecast hour                                    * 
 * IN/OUT  mip3    mask for ip3 (user defined identifier)                    * 
 * IN/OUT  metiket mask for label                                            * 
 *   IN    getmode logical (1: getmode 0:set mode)                           * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstmsq)(ftnword *f_iun, ftnword *f_mip1, ftnword *f_mip2,
                        ftnword *f_mip3, char *f_metiket, ftnword *f_getmode,
                        F2Cl ll1)
{
  int err, iun = *f_iun, mip1 = *f_mip1, mip2 = *f_mip2, mip3 = *f_mip3;
  int getmode = *f_getmode;
  int l1=ll1;

  char metiket[13];

  str_cp_init(metiket,13,f_metiket,l1);
  err = c_fstmsq(iun,&mip1,&mip2,&mip3,metiket,getmode);

  if (getmode) {
    *f_mip1 = (ftnword) mip1;
    *f_mip2 = (ftnword) mip2;
    *f_mip3 = (ftnword) mip3;
  }
  return(err);
}

/*splitpoint fstnbr */
/***************************************************************************** 
 *                              F S T N B R                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Returns the number of records of the file associated with unit number.  *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstnbr)(ftnword *f_iun)
{
  int iun = *f_iun;
  return ((ftnword) c_fstnbr(iun));
}  

/*splitpoint fstnbrv */
/***************************************************************************** 
 *                              F S T N B R V                                *
 *                                                                           * 
 *Object                                                                     * 
 *   Returns the number of valid records (excluding deleted records) of the  *
 *   file associated with unit number.                                       *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstnbrv)(ftnword *f_iun)
{
  int iun = *f_iun;
  return ((ftnword) c_fstnbrv(iun));
}  

/*splitpoint fstopc */
/*****************************************************************************
 *                              F S T O P C                                  *
 *                                                                           *
 *Object                                                                     *
 *   Print out or set a fstd or xdf global variable option.                  *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *   IN     option   option name to be set/printed                           *
 *   IN     value    option value                                            *
 *   IN     getmode  logical (1: get option, 0: set option)                  *
 *                                                                           *
 *****************************************************************************/
ftnword f77name(fstopc)(char *f_option, char *f_value, ftnword *f_getmode,
                        F2Cl ll1, F2Cl ll2)
{
  int getmode = *f_getmode, ier;
  char option[17];
  char value[129];
  int l1=ll1, l2=ll2;

  l1 = (l1 > 16) ? 16 : l1;
  l2 = (l2 > 128) ? 128 : l2;
  strncpy(option,f_option,l1);
  option[l1] = '\0';
  l1--;
  while ((l1 > 0) && (option[l1] == ' ')) {
    option[l1]='\0';
    l1--;
  }
  strncpy(value,f_value,l2);
  value[l2] = '\0';
  l2--;
  while ((l2 > 0) && (value[l2] == ' ')) {
    value[l2]='\0';
    l2--;
  }
  
  ier = c_fstopc(option,value,getmode);
  return((ftnword) ier);
}

/*splitpoint fstopi */
/*****************************************************************************
 *                              F S T O P I                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Print out or set a fstd or xdf global variable option.                  *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *   IN     option   option name to be set/printed                           * 
 *   IN     value    option value                                            * 
 *   IN     getmode  logical (1: get option, 0: set option)                  * 
 *                                                                           * 
 *****************************************************************************/
ftnword f77name(fstopi)(char *f_option, ftnword *f_value, ftnword * f_getmode,
                        F2Cl ll1)
{
  int getmode = *f_getmode, value = *f_value, ier;
  int l1=ll1;
  char option[7] = {' ',' ',' ',' ',' ',' ','\0'};
    
  l1 = (l1 > 6) ? 6 : l1;
  strncpy(option,f_option,l1);

  ier = c_fstopi(option,value,getmode);
  return((ftnword) ier);
}

/*splitpoint fstopl */
/***************************************************************************** 
 *                              F S T O P L                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Print out or set a fstd or xdf global variable option.                  *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *   IN     option   option name to be set/printed                           * 
 *   IN     value    option value                                            * 
 *   IN     getmode  logical (1: get option, 0: set option)                  * 
 *                                                                           * 
 *****************************************************************************/
ftnword f77name(fstopl)(char *f_option, ftnword *f_value, ftnword * f_getmode,
                        F2Cl ll1)
{
  int getmode = *f_getmode, value = *f_value, ier;
  int l1=ll1;
  char option[17];
    
  l1 = (l1 > 16) ? 16 : l1;
  strncpy(option,f_option,l1);
  option[l1] = '\0';

  ier = c_fstopl(option,value,getmode);
  return((ftnword) ier);
}

/*splitpoint fstopr */
/***************************************************************************** 
 *                              F S T O P R                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Print out or set a fstd or xdf global variable option.                  *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *   IN     option   option name to be set/printed                           * 
 *   IN     value    option value                                            * 
 *   IN     getmode  logical (1: get option, 0: set option)                  * 
 *                                                                           * 
 *****************************************************************************/
ftnword f77name(fstopr)(char *f_option, ftnfloat *f_value, ftnword * f_getmode,
                        F2Cl ll1)
{
  int getmode = *f_getmode, ier;
  float value = *f_value;
  int l1=ll1;
  char option[7];
    
  l1 = (l1 > 6) ? 6 : l1;
  strncpy(option,f_option,l1);
  option[l1] = '\0';

  ier = c_fstopr(option,value,getmode);
  return((ftnword) ier);
}

/*splitpoint fstouv */
/*****************************************************************************
 *                              F S T O U V                                  *
 *                                                                           *
 *Object                                                                     *
 *   Opens a RPN standard file.                                              *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  iun     unit number associated to the file                           *
 *  IN  options random or sequential access                                  *
 *                                                                           *
 *****************************************************************************/

ftnword f77name(fstouv)(ftnword *f_iun, char *options, F2Cl lng)
{
  int iun = *f_iun, ier;
  ier = c_fstouv(iun,options);
  return ((ftnword) ier);
}

/*splitpoint fstprm */
/*****************************************************************************
 *                              F S T P R M                                  *
 *                                                                           *
 *Object                                                                     *
 *   Get all the description informations of the record.                     *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  handle  positioning information to the record                        *
 *  OUT date    date time stamp                                              *
 *  OUT deet    length of a time step in seconds                             *
 *  OUT npas    time step number                                             *
 *  OUT ni      first dimension of the data field                            *
 *  OUT nj      second dimension of the data field                           *
 *  OUT nk      third dimension of the data field                            * 
 *  OUT nbits   number of bits kept for the elements of the field            * 
 *  OUT datyp   data type of the elements                                    * 
 *  OUT ip1     vertical level                                               * 
 *  OUT ip2     forecast hour                                                * 
 *  OUT ip3     user defined identifier                                      * 
 *  OUT typvar  type of field (forecast, analysis, climatology)              * 
 *  OUT nomvar  variable name                                                * 
 *  OUT etiket  label                                                        * 
 *  OUT grtyp   type of geographical projection                              * 
 *  OUT ig1     first grid descriptor                                        * 
 *  OUT ig2     second grid descriptor                                       * 
 *  OUT ig3     third grid descriptor                                        * 
 *  OUT ig4     fourth grid descriptor                                       * 
 *  OUT swa     starting word address                                        * 
 *  OUT lng     record length                                                * 
 *  OUT dltf    delete flag                                                  * 
 *  OUT ubc     unused bit count                                             * 
 *  OUT extra1  extra parameter                                              * 
 *  OUT extra2  extra parameter                                              * 
 *  OUT extra3  extra parameter                                              * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstprm)(ftnword *f_handle,
                        ftnword *f_dateo, ftnword *f_deet, ftnword *f_npas,
                        ftnword *f_ni, ftnword *f_nj, ftnword *f_nk,
                        ftnword *f_nbits, ftnword *f_datyp, ftnword *f_ip1,
                        ftnword *f_ip2, ftnword *f_ip3, char *f_typvar,
                        char *f_nomvar, char *f_etiket, char *f_grtyp,
                        ftnword *f_ig1, ftnword *f_ig2, ftnword *f_ig3,
                        ftnword *f_ig4, ftnword *f_swa, ftnword *f_lng,
                        ftnword *f_dltf, ftnword *f_ubc, ftnword *f_extra1,
                        ftnword *f_extra2, ftnword *f_extra3, 
                        F2Cl ll1, F2Cl ll2, F2Cl ll3, F2Cl ll4)
{
  int handle = *f_handle;
  int ni,nj,nk,nbits,datyp,ip1,ip2,ip3,dateo,deet,npas;
  int ig1,ig2,ig3,ig4,swa,lng,dltf,ubc,extra1,extra2,extra3,ier;
  int l1=ll1, l2=ll2, l3=ll3, l4=ll4;
  char etiket[13]={' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '
                   ,'\0'};
  char typvar[3]={' ',' ','\0'};
  char nomvar[5]={' ',' ',' ',' ','\0'};
  char grtyp[2]={' ','\0'};

  l1 = (l1 < 2) ? l1 : 2;            /* typvar length */
  l2 = (l2 < 4) ? l2 : 4;            /* nomvar length */
  l3 = (l3 < 12) ? l3 :12;           /* etiket length */
  l4 = (l4 < 1) ? l4 : 1;            /* grtyp length */

  ier = c_fstprm(handle,&dateo,&deet,&npas,&ni,&nj,&nk,
                     &nbits,&datyp,&ip1,&ip2,&ip3,typvar,
                     nomvar,etiket,grtyp,&ig1,&ig2,&ig3,&ig4,&swa,&lng,
                     &dltf,&ubc,&extra1,&extra2,&extra3);
  *f_ni = (ftnword) ni;
  *f_nj = (ftnword) nj;
  *f_nk = (ftnword) nk;
  *f_dateo = (ftnword) dateo;
  *f_deet = (ftnword) deet;
  *f_npas = (ftnword) npas;
  *f_nbits = (ftnword) nbits;
  *f_datyp = (ftnword) datyp;
  *f_ip1 = (ftnword) ip1;
  *f_ip2 = (ftnword) ip2;
  *f_ip3 = (ftnword) ip3;
  *f_ig1 = (ftnword) ig1;
  *f_ig2 = (ftnword) ig2;
  *f_ig3 = (ftnword) ig3;
  *f_ig4 = (ftnword) ig4;
  *f_swa = (ftnword) swa;
  *f_lng = (ftnword) lng;
  *f_dltf = (ftnword) dltf;
  *f_ubc = (ftnword) ubc;
  *f_extra1 = (ftnword) extra1;
  *f_extra2 = (ftnword) extra2;
  *f_extra3 = (ftnword) extra3;
  string_copy(f_typvar,typvar,l1);
  string_copy(f_nomvar,nomvar,l2);
  string_copy(f_etiket,etiket,l3);
  string_copy(f_grtyp,grtyp,l4);
  return((ftnword) ier);
}

/*splitpoint fstrwd */
/***************************************************************************** 
 *                   F S T R E S E T _ I P _ F L A G S                       *
 *                                                                           * 
 *Object                                                                     * 
 *   Reset all the flags previously set by ip(1-3)_val                       *
 *                                                                           * 
 *****************************************************************************/
void f77name(fstreset_ip_flags)()
{
  int ier;
  ier = init_ip_vals();
}


/*splitpoint fstrwd */
/***************************************************************************** 
 *                               F S T R W D                                 *
 *                                                                           * 
 *Object                                                                     * 
 *   Rewinds a RPN standard sequential file.                                 *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *                                                                           * 
 *****************************************************************************/
ftnword f77name(fstrwd)(ftnword *f_iun)
{
  int err, iun = *f_iun;

  err = c_fstrwd(iun);
  return((ftnword) err);
}

/*splitpoint fstskp */
/***************************************************************************** 
 *                                F S T S K P                                *
 *                                                                           * 
 *Object                                                                     * 
 *   Skip nrec records forward or backward in the sequential file.           *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *  IN  nrec    number of records to skip (negative nrec means backward)     * 
 *                                                                           * 
 *****************************************************************************/
ftnword f77name(fstskp)(ftnword *f_iun, ftnword *f_nrec)
{
  int iun = *f_iun, nrec = *f_nrec;
  int ier;

  ier = c_fstskp(iun,nrec);
  return((ftnword) ier);
}

/*splitpoint fstsui */
/***************************************************************************** 
 *                            F S T S U I                                    *
 *                                                                           * 
 *Object                                                                     * 
 *   Finds the next record that matches the last search criterias            *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *  OUT ni      dimension 1 of the data field                                * 
 *  OUT nj      dimension 2 of the data field                                * 
 *  OUT nk      dimension 3 of the data field                                * 
 *                                                                           * 
 *****************************************************************************/
ftnword f77name(fstsui)(ftnword *f_iun,
                        ftnword *f_ni, ftnword *f_nj, ftnword *f_nk)
{
  int iun = *f_iun;
  int ier,ni,nj,nk;

  ier = c_fstsui(iun,&ni,&nj,&nk);
  *f_ni = (ftnword) ni;
  *f_nj = (ftnword) nj;
  *f_nk = (ftnword) nk;
  return((ftnword) ier);
}

/*splitpoint fstunl */
/***************************************************************************** 
 *                              F S T U N L                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Unlinks a list of files previously linked by fstlnk.                    *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  liste   list of unit numbers associated to the files                 * 
 *  IN  n       number of files to link                                      * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstunl)()
{
  int ier;

  ier = c_xdfunl(link_list,link_n);
  return ((ftnword) ier);
}  

/*splitpoint fst_version */
/***************************************************************************** 
 *                           F S T  _ V E R S I O N                          *
 *                                                                           * 
 *Object                                                                     * 
 *   Returns package version number.                                          *
 *                                                                           * 
 *****************************************************************************/

wordint f77name(fst_version)()
{
  return((wordint) stdf_version);
}

/*splitpoint fstvoi */
/***************************************************************************** 
 *                              F S T V O I                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Opens a RPN standard file.                                              *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *  IN  options random or sequential access                                  * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstvoi)(ftnword *f_iun,char *f_options, F2Cl ll1)
{
  int iun = *f_iun, l1=ll1;
  char options[80] = 
  {' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
   ' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
   ' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
   ' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
   ' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
   ' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
   ' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
   ' ',' ',' ',' ',' ',' ',' ',' ',' ','\0'};

  l1 = (l1 > 79) ? 79 : l1;
  strncpy(options,f_options,l1);
  options[l1] = '\0';

  return ((ftnword) c_fstvoi(iun,options));
}

/*splitpoint fstweo */
/***************************************************************************** 
 *                             F S T W E O                                   *
 *                                                                           * 
 *Object                                                                     * 
 *   Writes a logical end of file on a sequential file.                      *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  iun     unit number associated to the file                           * 
 *  IN  level   level of logical end of file                                 * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(fstweo)(ftnword *f_iun, ftnword *f_level)
{
  int ier, iun = *f_iun, level = *f_level;
  
  ier = c_fstweo(iun,level);
  return((ftnword) ier);
}

/*splitpoint ip1_all */
/***************************************************************************** 
 *                             I P 1 _ A L L                                 *
 *                                                                           * 
 *Object                                                                     * 
 *   Generates all possible coded ip1 values for a given level               *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  level          ip1 level (float value)                               * 
 *  IN  kind           level kind as defined in convip                       * 
 *                                                                           * 
 *****************************************************************************/

ftnword f77name(ip1_all)(ftnfloat *f_level, ftnword *f_kind)
{
  int kind = *f_kind, ip1;
  float level = *f_level;
  
  ip1 = c_ip1_all(level,kind);
  return((ftnword) ip1);
}

ftnword f77name(ip2_all)(ftnfloat *f_level, ftnword *f_kind)
{
  int kind = *f_kind, ip2;
  float level = *f_level;
  
  ip2 = c_ip2_all(level,kind);
  return((ftnword) ip2);
}

ftnword f77name(ip3_all)(ftnfloat *f_level, ftnword *f_kind)
{
  int kind = *f_kind, ip3;
  float level = *f_level;
  
  ip3 = c_ip3_all(level,kind);
  return((ftnword) ip3);
}

ftnword f77name(ip1_val)(ftnfloat *f_level, ftnword *f_kind)
{
  int kind = *f_kind, ip1;
  float level = *f_level;
  
  ip1 = c_ip1_val(level,kind);
  return((ftnword) ip1);
}

ftnword f77name(ip2_val)(ftnfloat *f_level, ftnword *f_kind)
{
  int kind = *f_kind, ip2;
  float level = *f_level;
  
  ip2 = c_ip2_val(level,kind);
  return((ftnword) ip2);
}

ftnword f77name(ip3_val)(ftnfloat *f_level, ftnword *f_kind)
{
  int kind = *f_kind, ip3;
  float level = *f_level;
  
  ip3 = c_ip3_val(level,kind);
  return((ftnword) ip3);
}


/*splitpoint fst_can_translate_name */
static char exception_vars[256]="~^[<>!^]";  /* by default ignore names starting with >!^ */
static int read_done=0;
static regex_t pattern;

int FstCanTranslateName(char *varname)  /* is this name NOT FOUND in do no translate table */
{
  FILE *fileref;
  static char filename[256];
  char *FST_NOIP_NAME, *BASENAME;
  int result, i;
  regmatch_t match_table;

  if (! read_done) {  /* first call, get do not translate table */
    read_done=1;
    FST_NOIP_NAME=getenv("FST_NOIP_NAME");
    ARMNLIB=getenv("ARMNLIB");
    BASENAME=ARMNLIB;
    if(FST_NOIP_NAME)  /* environment variable contains the table */
    {
      strncpy( exception_vars , FST_NOIP_NAME , sizeof(exception_vars) );
      BASENAME=NULL;
      if(exception_vars[0]=='|') BASENAME=exception_vars+1; /* FST_NOIP_NAME contains a file name */
    }
    if(BASENAME){ /* get table from $ARMNLIB/data/exception_vars file if it exists */
      if(BASENAME==ARMNLIB)
        snprintf(filename,sizeof(filename),"%s/data/exception_regex_var",ARMNLIB);
      else
        snprintf(filename,sizeof(filename),"%s",BASENAME);
      if ((fileref = fopen(filename,"r")) != NULL) {
        if(NULL == fgets(exception_vars,sizeof(exception_vars),fileref) ) exception_vars[0]='\0' ;
        fprintf(stderr,"OPENING exception file: %s\n",filename);
        fclose(fileref) ;
      }
    }
    if(exception_vars[0]=='~') 
    {
      for (i=0 ; exception_vars[i]!='\0' && exception_vars[i]!='\n' ; i++) ; exception_vars[i]='\0';
      result = regcomp(&pattern,exception_vars+1,REG_EXTENDED|REG_NOSUB);
      if (msg_level < INFORM) fprintf(stderr,"exception pattern: '%s'\n",exception_vars+1);
    }
  }
  if(exception_vars[0]=='~')  /* this is a regex pattern */
  {
    result = regexec(&pattern,varname,(size_t) 0,NULL,0)!=0;  /* name not in pattern, it can be translated */
  }else{  /* this is a straight list of 4 char tokens */
    result = strstr(exception_vars,varname)==NULL ; /* name not in list, it can be translated */
  }
  return result;
}
/*splitpoint print_std_parms */
/***************************************************************************** 
 *                      P R I N T _ S T D _ P A R M S                        *
 *                                                                           * 
 *Object                                                                     * 
 *   Prints out the standard file record descriptors                         *
 *                                                                           * 
 *Revision                                                                   * 
 * 002   B. Dugas -   Oct 2012, ip23 buffer overflow correction              * 
 * 003   M. Valin -   Feb 2013, introduction of missing data flag            * 
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  stdf_entry     directory entry that contains the descriptors         * 
 *  IN  pre            preamble string                                       * 
 *  IN  option         print option (newstyle or oldstyle)                   * 
 *                                                                           * 
 *****************************************************************************/

static void print_std_parms(stdf_dir_keys *stdf_entry, char *pre, char *option,
                            int header)
{

  stdf_special_parms cracked;
  char cdt[9]={'X','R','I','C','S','E','F','A','Z'};
  char cmsgp=' '; /* initialize for case where there are no missing value(s) in record */
  int dat2,dat3,minus3=-3;
  int iip1,kind, mode=-1, flag=1;
  int ig1, ig2, ig3, ig4;
  float level;
  char c_level[16], pg1[7], pg2[7], pg3[8], pg4[8];
  char h_dims[23], h_dateo[16], h_stampo[10], h_datev[26], h_level[16], h_ip1[10], h_grid[32];
  char v_dims[23], v_dateo[16], v_stampo[10], v_datev[26], v_level[16], v_ip1[10], v_grid[32];
  char h_decoded[39], v_decoded[39];
  char h_nomv[5], h_typv[3], h_etiq[13], h_ip23[20], h_deet[9], h_npas[9], h_dty[5]; 
  char v_nomv[5], v_typv[3], v_etiq[13], v_ip23[20], v_deet[9], v_npas[9], v_dty[5]; 
  int posc, posv;
  static char *ARMNLIB=NULL;                        /* ARMNLIB environment variable */
  static char filename[256];
  FILE *fileref;
  
  /* printf("Debug+ print_std_parms option=%s\n",option); */
  crack_std_parms(stdf_entry,&cracked);

  if (header) {                      /* build and print header line */

    if (strstr(option,"NONOMV"))
      h_nomv[0]='\0';
    else
      snprintf(h_nomv,sizeof(h_nomv),"%s","NOMV");

    if (strstr(option,"NOTYPV"))
      h_typv[0]='\0';
    else
      snprintf(h_typv,sizeof(h_typv),"%s","TV");

    if (strstr(option,"NOETIQ"))
      h_etiq[0]='\0';
    else
      snprintf(h_etiq,sizeof(h_etiq),"%s","  ETIQUETTE ");

    if (strstr(option,"NINJNK"))
      snprintf(h_dims,sizeof(h_dims),"%s","      NI      NJ    NK");
    else
      h_dims[0]='\0';

    if (strstr(option,"DATEO"))
      /*      snprintf(h_dateo,"%s","YYYYMMDD HHMMSS"); */
      snprintf(h_dateo,sizeof(h_dateo),"%s","(DATE-O  h m s)");
    else
      h_dateo[0]='\0';

    if (strstr(option,"DATESTAMPO"))
      snprintf(h_stampo,sizeof(h_stampo),"%s","  STAMP-O");
    else
      h_stampo[0]='\0';

    if (strstr(option,"DATEV"))
      /*      snprintf(h_datev,"%s","YYYYMMDD HHMMSS     DATEV"); */
      snprintf(h_datev,sizeof(h_datev),"%s","(DATE-V  h m s)   STAMP-V");
    else
      h_datev[0]='\0';

    if (strstr(option,"LEVEL"))
      snprintf(h_level,sizeof(h_level),"%s","       LEVEL   ");
    else
      h_level[0]='\0';

    if (strstr(option,"IPALL"))
      snprintf(h_decoded,sizeof(h_decoded),"%s","          DECODED IP1/IP2/IP3         ");
    else
      h_decoded[0]='\0';
    
    if (strstr(option,"IP1"))
      snprintf(h_ip1,sizeof(h_ip1),"%s","      IP1");
    else
      h_ip1[0]='\0';

    if (strstr(option,"NOIP23"))
      h_ip23[0]='\0';
    else
      snprintf(h_ip23,sizeof(h_ip23),"%s","      IP2       IP3");

    if (strstr(option,"NODEET"))
      h_deet[0]='\0';
    else
      snprintf(h_deet,sizeof(h_deet),"%s","    DEET");

    if (strstr(option,"NONPAS"))
      h_npas[0]='\0';
    else
      snprintf(h_npas,sizeof(h_npas),"%s","    NPAS");

    if (strstr(option,"NODTY"))
      h_dty[0]='\0';
    else
      snprintf(h_dty,sizeof(h_dty),"%s","DTY ");

    if (strstr(option,"GRIDINFO"))
      snprintf(h_grid,sizeof(h_grid),"%s","G    XG1    XG2     XG3     XG4");
    else
      if (strstr(option,"IG1234"))
        snprintf(h_grid,sizeof(h_grid),"%s","G   IG1   IG2   IG3   IG4");
      else
        h_grid[0]='\0';

    fprintf(stdout,"\n       %s %s %s %s %s %s %s %s %s %s %s %s %s  %s  %s\n\n",h_nomv,h_typv,h_etiq,h_dims,h_dateo,h_stampo,h_datev,h_level,h_decoded,h_ip1,h_ip23,h_deet,h_npas,h_dty,h_grid);
    /*    fprintf(stdout,"\n       NOMV TV ETIQUETTE       NI    NJ    NK %s %s %s %s %s   IP2   IP3     DEET     NPAS  DTY  %s\n\n",h_dateo,h_stampo,h_datev,h_level,h_ip1,h_grid); */
  }

  if (strstr(option,"NONOMV"))
    v_nomv[0]='\0';
  else
    snprintf(v_nomv,sizeof(v_nomv),"%4s",cracked.nomvar);
  
  if (strstr(option,"NOTYPV"))
    v_typv[0]='\0';
  else
    snprintf(v_typv,sizeof(v_typv),"%2s",cracked.typvar);

  if (strstr(option,"NOETIQ"))
    v_etiq[0]='\0';
  else
    snprintf(v_etiq,sizeof(v_etiq),"%12s",cracked.etiket);

  if (strstr(option,"NINJNK"))
    snprintf(v_dims,sizeof(v_dims)," %7d %7d %5d",stdf_entry->ni,stdf_entry->nj,stdf_entry->nk);
  else
    v_dims[0]='\0';

  if (strstr(option,"DATEO")) {
    f77name(newdate)(&cracked.date_stamp,&dat2,&dat3,&minus3);
    snprintf(v_dateo,sizeof(v_dateo),"%08d %06d",dat2,dat3/100);
  }
  else
    v_dateo[0]='\0';

  if (strstr(option,"DATESTAMPO"))
    snprintf(v_stampo,sizeof(v_stampo),"%09d",cracked.date_stamp);
  else
    v_stampo[0]='\0';

  if (strstr(option,"DATEV")) {
    f77name(newdate)(&cracked.date_valid,&dat2,&dat3,&minus3);
    if (cracked.date_valid < -1)
      snprintf(v_datev,sizeof(v_datev),"%08d %06d %10d",dat2,dat3/100,cracked.date_valid);
    else
      snprintf(v_datev,sizeof(v_datev),"%08d %06d %09d",dat2,dat3/100,cracked.date_valid);
  }
  else
    v_datev[0]='\0';

  v_level[0]='\0';
  v_decoded[0]='\0';
  if ( strstr(option,"LEVEL") || strstr(option,"IPALL") )
    {
/*
      if (! read_done) {
        ARMNLIB=getenv("ARMNLIB");
        sprintf(filename,"%s/data/exception_vars",ARMNLIB);
        if ((fileref = fopen(filename,"r")) != NULL) {
          if(NULL == fgets(exception_vars,sizeof(exception_vars),fileref) ) exception_vars[0]='\0' ;
          fclose(fileref) ;
        }
        read_done=1;
      }
*/
      iip1 = stdf_entry->ip1;
//      if (strstr(exception_vars,cracked.nomvar)) {     /* special variable, no decoding */
      if(! FstCanTranslateName(cracked.nomvar)) {
        snprintf(c_level,sizeof(c_level),"%12d   ",iip1);
        if (strstr(option,"LEVEL")) snprintf(v_level,sizeof(v_level),"%15s","     -----     ");
//        if (strstr(option,"IPALL")) snprintf(v_decoded,sizeof(v_decoded),"%16s------%16s","","");
        if (strstr(option,"IPALL")) snprintf(v_decoded,sizeof(v_decoded),"[%10d] [%10d] [%10d]",stdf_entry->ip1,stdf_entry->ip2,stdf_entry->ip3);
      }
      else     /* not a special variable  */
      {
        if (strstr(option,"LEVEL")) /* good old level option */
        {
          f77name(convip)(&iip1,&level,&kind,&mode,c_level,&flag,(F2Cl) 15);
          c_level[15] = '\0';
          snprintf(v_level,sizeof(v_level),"%s","               ");        /* blank initialisation */
          posc=14;
          posv=14;
          while ((posc >= 0) && (isspace(c_level[posc])))  /* skip blanks and right justify string */
            posc--;
          if (isdigit(c_level[posc]))
            posv -= 3;
          while ((posv >=0) && (posc >=0)) {
            v_level[posv] = c_level[posc];
            posv--;
            posc--;
          }
        }
        if (strstr(option,"IPALL")) /* full IP1/IP2/IP3 triplet decoding */
        {
          float p1,p2,p3;
          int kind1,kind2,kind3, StatusIP;
          StatusIP=ConvertIPtoPK(&p1,&kind1,&p2,&kind2,&p3,&kind3,stdf_entry->ip1,stdf_entry->ip2,stdf_entry->ip3);
          if(kind1<0 || kind2<0 || kind3<0 || (StatusIP&CONVERT_ERROR) ) { /* decode error somewhere */
            kind1 = 15; kind2 = 15 ; kind3 = 15;  /* integer code P=IP */
            p1 = stdf_entry->ip1 ; p2 = stdf_entry->ip2 ; p3 = stdf_entry->ip3;
          }
          kind1 &= 0x1F ; kind2 &= 0x1F ; kind3 &= 0x1F ;   /* force modulo 32 */
          snprintf(v_decoded,sizeof(v_decoded),"%10g%s %10g%s %10g%s",p1,kinds(kind1),p2,kinds(kind2),p3,kinds(kind3));
        }
      }     /* special variable, no decoding */
    }

  if (strstr(option,"IP1"))
    snprintf(v_ip1,sizeof(v_ip1),"%9d",stdf_entry->ip1);
  else
    v_ip1[0]='\0';

  if (strstr(option,"NOIP23"))
    v_ip23[0]='\0';
  else
    snprintf(v_ip23,sizeof(v_ip23),"%9d %9d",stdf_entry->ip2,stdf_entry->ip3);

  if (strstr(option,"NODEET"))
    v_deet[0]='\0';
  else
    snprintf(v_deet,sizeof(v_deet),"%8d",stdf_entry->deet);
  
  if (strstr(option,"NONPAS"))
    v_npas[0]='\0';
  else
    snprintf(v_npas,sizeof(v_npas),"%8d",stdf_entry->npas);
  
  if(stdf_entry->datyp & 64)cmsgp='m';  /* m will be added to data type if there are missing values in record */
  if (strstr(option,"NODTY"))
    v_dty[0]='\0';
  else 
    if (stdf_entry->datyp > 128)  /* force lower case data type code if compressed */
      snprintf(v_dty,sizeof(v_dty),"%1c%1c%2d",tolower(cdt[stdf_entry->datyp&0x3F]),cmsgp,stdf_entry->nbits);  /* suppress bits for 64 and 128 */
    else
      snprintf(v_dty,sizeof(v_dty),"%1c%1c%2d",cdt[stdf_entry->datyp&0x3F],cmsgp,stdf_entry->nbits);  /* suppress bits for 64 and 128 */

  if (strstr(option,"GRIDINFO")) {
    F2Cl lc1=1,lc2=7,lc3=7,lc4=8,lc5=8;
    ig1=stdf_entry->ig1; ig2=cracked.ig2;
    ig3=stdf_entry->ig3; ig4=stdf_entry->ig4;
    f77name(igapg)(cracked.gtyp,pg1,pg2,pg3,pg4,&ig1,&ig2,&ig3,&ig4,
                   lc1,lc2,lc3,lc4,lc5);
            /*     1,7,7,8,8);       */
    pg1[6]='\0'; pg2[6]='\0'; pg3[7]='\0'; pg4[7]='\0';
    snprintf(v_grid,sizeof(v_grid),"%1s %6s %6s %7s %7s",cracked.gtyp,pg1,pg2,pg3,pg4);
  }
  else
    if (strstr(option,"IG1234"))
      snprintf(v_grid,sizeof(v_grid),"%1s %5d %5d %5d %5d",cracked.gtyp,stdf_entry->ig1,cracked.ig2,stdf_entry->ig3,stdf_entry->ig4);
    else
      v_grid[0]='\0';
  
  fprintf(stdout,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s  %s\n",
          pre,v_nomv,v_typv,v_etiq,
          v_dims,v_dateo,v_stampo,v_datev,
          v_level,v_decoded,v_ip1,v_ip23,
          v_deet,v_npas,v_dty,
          v_grid);
}

/*splitpoint zzz_stub */
/***************************************************************************** 
 *                               S T U B                                     *
 *                                                                           * 
 *Object                                                                     * 
 *   Stubs for standard file routines not implemented yet.                   *
 *                                                                           * 
 *****************************************************************************/
ftnword f77name(fstabt)()
{
  sprintf(errmsg,"this routine is not implemented in FSTD98");
  return(error_msg("FSTABT",ERR_NOT_IMPL,ERRFATAL));
}
ftnword f77name(fstsel)()
{
  sprintf(errmsg,"this routine is not implemented in FSTD98\n \t\t fstinfx or fstlirx must be used instead");
  return(error_msg("FSTSEL",ERR_NOT_IMPL,WARNING));
}
ftnword f77name(zfstcvt)()
{
  sprintf(errmsg,"this routine is not implemented yet in FSTD98");
  return(error_msg("FSTCVT",ERR_NOT_IMPL,ERROR));
}
ftnword f77name(fstpos)()
{
  sprintf(errmsg,"this routine is not implemented in FSTD98\n \t\t fstinfx or fstlirx must be used instead");
  return(error_msg("FSTPOS",ERR_NOT_IMPL,WARNING));
}
