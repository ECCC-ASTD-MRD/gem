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

/*****************************************************************************
 *                                                                           *
 *INTERFACE C POUR BURP                                                      *
 *                                                                           *
 *AUTEUR  M. LEPINE  -  NOV 1990                                             *
 *                                                                           *
 *OBJET                                                                      *
 *     INTERFACER LES SOUS-PROGRAMMES FORTRAN DU LOGICIEL BURP AVEC LE       *
 *     LANGAGE C EN TENANT COMPTE DES PARTICULARITES DES APPELS              *
 *     INTER-LANGAGES                                                        *
 *                                                                           *
 *         EX:   PASSAGE D'ARGUMENT PAR ADRESSE EN FORTRAN                   *
 *               PASSAGE SUPPLEMENTAIRE DE LA LONGUEUR POUR UN ARGUMENT DE   *
 *               TYPE CARACTERE ...                                          *
 *                                                                           *
 *****************************************************************************/
#include <rpnmacros.h>

/*****************************************************************************
 *                              C _ M R B C O L                              *
 *****************************************************************************/
int
c_mrbcol(liste,cliste,nele)
int liste[], cliste[], nele;

   {
   int lnele;
   lnele = nele;
   return(f77name(mrbcol)(liste,cliste,&lnele));
   }
/*****************************************************************************
 *                              C _ M R B C O V                              *
 *****************************************************************************/
int 
c_mrbcov(elem)
int elem;

   {
   int lelem;
   lelem = elem;
   return(f77name(mrbcov)(&lelem));
   }




/*****************************************************************************
 *                              C _ M R B C V T                              *
 *****************************************************************************/
int
c_mrbcvt(liste,tblval,rval,nele,nval,nt,mode)
int liste[], tblval[], nele, nval, nt, mode;
float rval[];

   {
   int lnele, lnval, lnt, lmode;
   lnele = nele; lnval = nval; lnt = nt; lmode = mode;
   return(f77name(mrbcvt)(liste,tblval,rval,&lnele,&lnval,&lnt,&lmode));
   }

/*****************************************************************************
 *                              C _ M R B D C L                              *
 *****************************************************************************/
int
c_mrbdcl(cliste,liste,nele)
int liste[], cliste[], nele;

   {
   int lnele;
   lnele = nele;
   return(f77name(mrbdcl)(cliste,liste,&lnele));
   }

/*****************************************************************************
 *                              C _ M R B D C V                              *
 *****************************************************************************/
int 
c_mrbdcv(elem)
int elem;

   {
   int lelem;
   lelem = elem;
   return(f77name(mrbdcv)(&lelem));
   }


/*****************************************************************************
 *                              C _ M R B I N I                              *
 *****************************************************************************/
int
c_mrbini(iun,buf,temps,flgs,stnid,idtp,lati,longi,dx,dy,elev,drcv,date,
         oars,runn,sup,nsup,xaux,nxaux)
int buf[],temps,flgs,idtp,lati,longi,elev,drcv,date,oars,runn,sup[],nsup;
int dx, dy;
int xaux[],nxaux;
char stnid[];

  {
  int liun,ltemps,lflgs,lidtp,llati,llongi,lelev,ldrcv,ldate,loars;
  int ldx, ldy;
  int lrunn,lnsup,lnxaux;
  F2Cl l1;
  ltemps = temps; lflgs = flgs; lidtp = idtp; llati = lati; llongi = longi;
  ldx = dx; ldy = dy;
  lelev = elev; ldrcv = drcv; ldate = date; loars = oars; lrunn = runn;
  lnsup = nsup; lnxaux = nxaux; liun = iun;
  l1 = strlen(stnid);
  return(f77name(mrbini)(&liun,buf,&ltemps,&lflgs,stnid,&lidtp,&llati,&llongi,
                 &ldx,&ldy,&lelev,&ldrcv,&ldate,&loars,&lrunn,sup,&lnsup,xaux,
                 &lnxaux,l1));
  }

/*****************************************************************************
 *                              C _ M R B L O C X                            *
 *****************************************************************************/
int
c_mrblocx(buf,bfam,bdesc,bknat,bktyp,bkstp,blk0)
int buf[],bfam,bdesc,bknat,bktyp,bkstp,blk0;

  {
  int lbfam,lbdesc,lbknat,lbktyp,lbkstp,lblk0;
  lbfam = bfam; lbdesc = bdesc; lblk0 = blk0;
  lbknat = bknat; lbktyp = bktyp; lbkstp = bkstp;
  return(f77name(mrblocx)(buf,&lbfam,&lbdesc,&lbknat,&lbktyp,&lbkstp,&lblk0));
  }

/*****************************************************************************
 *                              C _ M R B P R M L                            *
 *****************************************************************************/
int
c_mrbprml(buf,bkno,tblprm,nprm,inblocs)
int buf[],tblprm[],bkno,nprm,inblocs;

  {
    int lbkno,lnprm,linblocs;

    lbkno = bkno; lnprm = nprm; linblocs = inblocs;

    return(f77name(mrbprml)(buf,&lbkno,tblprm,&lnprm,&linblocs));

  }

 
/*****************************************************************************
 *                              C _ M R B R P T                              *
 *****************************************************************************/
int 
c_mrbrpt(elem)
int elem;

   {
   int lelem;
   lelem = elem;
   return(f77name(mrbrpt)(&lelem));
   }


/*****************************************************************************
 *                              C _ M R B S C T                              *
 *****************************************************************************/
int
c_mrbsct(tablusr,neleusr)
int tablusr[], neleusr;

  {
  int lneleusr;
  lneleusr = neleusr;
  return(f77name(mrbsct)(tablusr,&lneleusr));
  }

/*****************************************************************************
 *                              C _ M R B T B L                              *
 *****************************************************************************/
int
c_mrbtbl(tablusr,nslots,neleusr)
int tablusr[], neleusr,nslots;

  {
  int lneleusr, lnslots;
  lneleusr = neleusr;
  lnslots  = nslots;
  return(f77name(mrbtbl)(tablusr,&lnslots,&lneleusr));
  }


/*****************************************************************************
 *                              C _ M R B T Y P                              *
 *****************************************************************************/
int
c_mrbtyp(hbknat,hbktyp,hbkstp,hbtyp)
int hbtyp, *hbknat, *hbktyp, *hbkstp;

{
      int lbtyp, iii;
      lbtyp = hbtyp;
      return(f77name(mrbtyp)(hbknat,hbktyp,hbkstp,&lbtyp));
}

/*****************************************************************************
 *                              C _ M R B U P D                              *
 *****************************************************************************/
int
c_mrbupd(iun,buf,temps,flgs,stnid,idtp,lati,longi,dx,dy,elev,drcv,date,
         oars,runn,sup,nsup,xaux,nxaux)
int buf[],temps,flgs,idtp,lati,longi,elev,drcv,date,oars,runn,sup[],nsup;
int dx, dy;
int xaux[],nxaux;
char stnid[];

  {
  int liun,ltemps,lflgs,lidtp,llati,llongi,lelev,ldrcv,ldate,loars;
  int ldx, ldy;
  int lrunn,lnsup,lnxaux;
  F2Cl l1;
  ltemps = temps; lflgs = flgs; lidtp = idtp; llati = lati; llongi = longi;
  ldx = dx; ldy = dy;
  lelev = elev; ldrcv = drcv; ldate = date; loars = oars; lrunn = runn;
  lnsup = nsup; lnxaux = nxaux; liun = iun;
  l1 = strlen(stnid);
  return(f77name(mrbupd)(&liun,buf,&ltemps,&lflgs,stnid,&lidtp,&llati,&llongi,
                 &ldx,&ldy,&lelev,&ldrcv,&ldate,&loars,&lrunn,sup,&lnsup,xaux,
                 &lnxaux,l1));
  }

/*****************************************************************************
 *                              C _ M R F C L S                              *
 *****************************************************************************/
int
c_mrfcls(iun)
int iun;

  {
  int liun;
  liun = iun;
  return(f77name(mrfcls)(&liun));
  }

/*****************************************************************************
 *                              C _ M R F G O C                              *
 *****************************************************************************/
int
c_mrfgoc(optnom,opvalc)
char optnom[],opvalc[9];

  {
  F2Cl l1,l2;
  int iii;
  l1 = strlen(optnom);
  l2 = strlen(opvalc);
  iii = f77name(mrfgoc)(optnom,opvalc,l1,l2);
  opvalc[8] = '\0';
  return (iii);
  }

/*****************************************************************************
 *                              C _ M R F G O R                              *
 *****************************************************************************/
int
c_mrfgor(optnom,opvalr)
char optnom[];
float *opvalr;

  {
  F2Cl l1;
  l1 = strlen(optnom);
  return(f77name(mrfgor)(optnom,opvalr,l1));
  }

/*****************************************************************************
 *                              C _ M R F L O C                              *
 *****************************************************************************/
int
c_mrfloc(iun,handle,stnid,idtyp,lat,lon,date,temps,sup,nsup)
int iun,handle,idtyp,lat,lon,date,temps,sup[],nsup;
char stnid[];

  {
  int liun,lhandle,lidtyp,llat,llon,ldate,ltemps,lnsup;
  F2Cl l1;
  l1 = strlen(stnid);
  liun = iun; lhandle = handle; lidtyp = idtyp; llat = lat; llon = lon;
  ldate = date; ltemps = temps; lnsup = nsup;
  return(f77name(mrfloc)(&liun,&lhandle,stnid,&lidtyp,&llat,&llon,&ldate,
                 &ltemps,sup,&lnsup,l1));
  }

/*****************************************************************************
*                               C _ M R F M X L                              *
******************************************************************************/
int
c_mrfmxl(iun)
int iun;
   {
   int liun;
   liun = iun;
   return(f77name(mrfmxl)(&liun));
   }


/*****************************************************************************
 *                              C _ M R F N B R                              *
 *****************************************************************************/
int
c_mrfnbr(iun)
int iun;

  {
  int liun;
  liun = iun;
  return(f77name(mrfnbr)(&liun));
  }

/*****************************************************************************
 *                              C _ M R F O P N                              *
 *****************************************************************************/
int
c_mrfopn(iun,mode)
int iun;
char mode[];

  {
  int liun;
  F2Cl l1;
  liun = iun;
  l1 = strlen(mode);
  return(f77name(mrfopn)(&liun,mode,l1));
  }

/*****************************************************************************
 *                              C _ M R F O P C                              *
 *****************************************************************************/
int
c_mrfopc(optnom,opvalc)
char optnom[],opvalc[];

  {
  F2Cl l1,l2;
  l1 = strlen(optnom);
  l2 = strlen(opvalc);
  return(f77name(mrfopc)(optnom,opvalc,l1,l2));
  }

/*****************************************************************************
 *                              C _ M R F O P R                              *
 *****************************************************************************/
int
c_mrfopr(optnom,opvalr)
char optnom[];
float opvalr;

  {
  F2Cl l1;
  float lopvalr;
  l1 = strlen(optnom);
  lopvalr = opvalr;
  return(f77name(mrfopr)(optnom,&lopvalr,l1));
  }

/*****************************************************************************
 *                              C _ M R F P R M                              *
 *****************************************************************************/
int
c_mrfprm(handle,stnid,idtyp,lat,lon,dx,dy,date,temps,flgs,sup,nsup,lng)
int handle,*idtyp,*lat,*lon,*date,*temps,*flgs,sup[],nsup,*lng;
int *dx, *dy;
char stnid[10];

  {
  int lhandle,lnsup,iii;
  F2Cl l1;
  lhandle = handle; lnsup = nsup;
  l1 = strlen(stnid);
  iii = f77name(mrfprm)(&lhandle,stnid,idtyp,lat,lon,dx,dy,date,temps,flgs,
                sup,&lnsup,lng,l1);
  stnid[9] = '\0';
  return(iii);
  }
       
/*****************************************************************************
 *                              C _ M R F V O I                              *
 *****************************************************************************/
int 
c_mrfvoi(iun)
int iun;

  {
  int liun;
  liun = iun;
  return(f77name(mrfvoi)(&liun));
  }

/******************************************************************************
*                               C _ M R F D E L                               *
*******************************************************************************/
int
c_mrfdel(handle)
int handle;

    {
    int lhandle;
    lhandle = handle;
    return(f77name(mrfdel)(&lhandle));
    }

