*/* RMNLIB - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
* *                          Environnement Canada
* *
* * This library is free software; you can redistribute it and/or
* * modify it under the terms of the GNU Lesser General Public
* * License as published by the Free Software Foundation,
* * version 2.1 of the License.
* *
* * This library is distributed in the hope that it will be useful,
* * but WITHOUT ANY WARRANTY; without even the implied warranty of
* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* * Lesser General Public License for more details.
* *
* * You should have received a copy of the GNU Lesser General Public
* * License along with this library; if not, write to the
* * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
* * Boston, MA 02111-1307, USA.
* */
*
**********************************************************************
***s/r hyb2pres - to return pressure given a unit number and 
*                 a list of fstkeys of the RPN standard file
**********************************************************************
*
      integer function hyb2pres(iun,fstkeys,NK,NI,NJ,PX,logPX_L)
*
      implicit none
*
      integer iun,NI,NJ,NK
      real    PX(NI,NJ,NK) 
      integer fstkeys(NK)
      logical logPX_L
*
*author
*     V.Lee Oct 2009 
*
*revision
* ________________________________________________________________
*  Name        I/O      Description
* ----------------------------------------------------------------
* iun           I       unit number to input RPNstd file
* fstkeys       I       list of keys tagged to records from RPNstd file
* NK            I       number of keys in fstkeys given
* NI            I       Ni dimension of field
* NJ            I       Nj dimension of field
* PX            O       if logPX_L=F, pressure field is in mb
*                       if logPX_L=T, pressure field in ln of Pascals
* logPX_L       I       True to have output PX in ln of Pascals
*                       False to have output PX in mb.
* ----------------------------------------------------------------
* the function will return 0 upon success, -1 if there is an error
**

      integer   fstinf, fstprm, fstluk, fnom,
     $          read_decode_hyb, hyb_to_pres,eta_to_pres,
     $          sigma_to_pres,etasef_to_pres
      external  fstinf, fstprm, fstluk, fnom,
     $          read_decode_hyb, hyb_to_pres,eta_to_pres,
     $          sigma_to_pres,etasef_to_pres

      integer   ip1(NK)
      real      hyb(NK),p0(NI,NJ),work(NI,NJ)
      integer   nia,nja,nka,ni1,nj1,nk1,i,j,k
      integer   e1_key,hy_key,pt_key,p0_key,xx_key
      integer   datev, dateo, deet, ipas, ip1a, ip2a, ip3a, 
     $          ig1a, ig2a, ig3a, ig4a, bit, datyp,
     $          swa, lng, dlf, ubc, ex1, ex2, ex3, kind
      real      lev,ptop,pref,rcoef,etatop
      character*1   tva, grda, blk_S
      character*4   var
      character*12  etik_S

      hyb2pres=0

      hyb2pres = fstprm (fstkeys(1), dateo, deet, ipas, nia, nja, nka, 
     $     bit, datyp, ip1a,ip2a,ip3a, tva, var, etik_S, grda, 
     $     ig1a,ig2a,ig3a,ig4a, swa,lng, dlf, ubc, ex1, ex2, ex3 )
      if (hyb2pres.lt.0) then
          print *,'hyb2pres error: fstprm on key',fstkeys(1)
          return
      else
         call convip_plus (ip1a, lev, kind,-1, blk_S, .false.)
         ip1(1)=ip1a
         hyb(1)=lev
      endif
      do k=2,NK
         hyb2pres = fstprm(fstkeys(k), dateo, deet, ipas, ni1, nj1, nk1, 
     $     bit, datyp, ip1a,ip2a,ip3a, tva, var, etik_S, grda, 
     $     ig1a,ig2a,ig3a,ig4a, swa,lng, dlf, ubc, ex1, ex2, ex3 )
           if (ni1.ne.nia.and.nj1.ne.nja.and.nk1.ne.nka.or.
     $                                               hyb2pres.lt.0) then
             print *,'hyb2pres error: fstprm on key',fstkeys(k),
     $               'dim mismatch'
             return
         endif
         call convip_plus (ip1a, lev, kind,-1, blk_S, .false.)
         ip1(k)=ip1a
         hyb(k)=lev
      enddo
      if (kind.ne.1.and.kind.ne.2.and.kind.ne.5) then
          print *,'hyb2pres error: kind = ',kind,' has to be 1,2 or 5'
          hyb2pres = -1
          return
      endif
      hy_key=fstinf (iun,ni,nj,nk,-1,' ',-1,  -1,  -1,' ','HY')
      pt_key=fstinf (iun,ni,nj,nk,-1,' ',-1,  -1,  -1,' ','PT')
      e1_key=fstinf (iun,ni,nj,nk,-1,' ',-1,  -1,  -1,' ','E1')

      call incdatr(datev,dateo,ipas*deet/3600.0d0)

      if (kind.eq.1) then
          p0_key=fstinf(iun,ni,nj,nk,datev,etik_S,-1,ip2a,ip3a,' ','P0')
          if (p0_key.lt.0) then
              print *,'hyb2pres error: No p0 found, kind = 1'
              hyb2pres = -1
              return
          else
              hyb2pres = fstluk(p0,p0_key,ni,nj,nk)
          endif
          if (pt_key.ge.0) then
              hyb2pres = fstluk(work,pt_key,ni,nj,nk)
              if (hyb2pres.lt.0) then
                   print *,'hyb2pres error: in fstluk PT'
                   return
              endif
              ptop = work(1,1)
              if (e1_key.ge.0) then
* etasef coordinate found
                  hyb2pres = fstluk(work,e1_key,ni,nj,nk)
                  if (hyb2pres.lt.0) then
                      print *,'hyb2pres error: in fstluk E1'
                      return
                  endif
                  etatop = work(1,1)
                  hyb2pres=etasef_to_pres(PX,hyb,ptop,etatop,p0,
     $                                    NI,NJ,NK)
                  if (logPX_L) PX(:,:,:)=log( 100.0*PX(:,:,:))
                  return
              else
* eta coordinate found
                  hyb2pres=eta_to_pres(PX,hyb,ptop,p0,NI,NJ,NK)
                  if (logPX_L) PX(:,:,:)=log( 100.0*PX(:,:,:))
                  return
              endif
          else if (hy_key.ge.0) then
* hybrid (normalized) coordinate found
              hyb2pres=read_decode_hyb (iun,'HY',  -1,  -1,' ',
     $                                           -1, ptop,pref,rcoef)
              if (hyb2pres.lt.0) then
                  print *,'hyb2pres error: in read_decode_hyb'
                  return
              endif
              hyb2pres=hyb_to_pres(PX,hyb,ptop,rcoef,pref,kind,p0,NI,NJ,
     $                             NK)
              if (logPX_L) PX(:,:,:)=log( 100.0*PX(:,:,:))
              return
          else 
* sigma coordinate found
              hyb2pres=sigma_to_pres(PX,hyb,p0,NI,NJ,NK)
              if (logPX_L) PX(:,:,:)=log( 100.0*PX(:,:,:))
          endif
      endif

      if (kind.eq.2) then
* pressure coordinate found
          do k=1,NK
          do j=1,nJ
          do i=1,nI
             PX(i,j,k)=hyb(k)*100.0
          enddo
          enddo
          enddo
          return
      endif

      if (kind.eq.5) then
* vstag coordinate found
         xx_key=fstinf (iun,ni, nj, nk, -1,etik_S,-1,-1,-1 ,' ','!!  ')
         if (xx_key.ge.0) then
             print *,'ERROR: coordinate is not ready'
             return
         else if (hy_key.ge.0) then
* hybrid (un-normalized) coordinate found
              hyb2pres=read_decode_hyb (iun,'HY',  -1,  -1,' ',
     $                                           -1,ptop,pref,rcoef)
              if (hyb2pres.lt.0) then
                  print *,'hyb2pres error: in read_decode_hyb'
                  return
              endif
              p0_key=fstinf (iun,ni,nj,nk,datev,etik_S,-1,ip2a,ip3a,' ',
     $                       'P0')
              if (p0_key.lt.0) then
                 print *,'hyb2pres error: HY found,No p0 found, kind=5'
                 hyb2pres = -1
                 return
              else
                 hyb2pres = fstluk(p0,p0_key,ni,nj,nk)
              endif
              hyb2pres=hyb_to_pres(PX,hyb,ptop,rcoef,pref,kind,p0,NI,NJ,
     $                             NK)
              return
         else
             print *,'hyb2pres error: kind=5 but !! nor  HY NOT FOUND'
             hyb2pres = -1
             return
         endif
      endif
      end
*
**********************************************************************
***function hybrid_to_pres  - function to convert from hybrid to pressure
*                    WARNING : only good for normalized hybrid
**********************************************************************
*
      integer function hybrid_to_pres(pressure,hybm,ptop,ps,NI,NJ,rcoef,
     %                                                 pref,hyb,NK)

      implicit none
      integer NI,NJ,NK
      real pressure(NI*NJ,NK),ptop,ps(NI*NJ),rcoef,pref,hyb(NK),hybm(NK)
*
*author
*     Vivian Lee/Michel Valin    Nov.28, 2001
*
*revision
* Lee V. - warning added to routine (check for normalized versus unnormalized
*                                    hybrid levels is absent)
*
*object
*     To derive pressure fields and model hybrid levels (levels used by
*     the model) from user-defined hybrid levels and the hybrid reference 
*     parameters(ptop,rcoef,pref)
*arguments
* ________________________________________________________________
*  Name        I/O      Description
* ----------------------------------------------------------------
* pressure      O       array of pressure levels (same units as ps)
* hybm          O       array of model hybrid levels (0.0 to 1.0)
*                       - calculated using ptop,rcoef and pref
* ptop          I       average pressure at the top (mb)
* ps            I       pressure at the surface (mb or pascals)
* rcoef         I       coefficient (1.0 to 2.0)
* pref          I       reference pressure in mb (normally = 800mb)
* hyb           I       array of user-defined hybrid levels (0.0 to 1.0)
* NI            I       Ni dimension of field
* NJ            I       Nj dimension of field
* NK            I       number of levels in hybm
* ----------------------------------------------------------------
* the function will return 0 upon success, -1 if there is an error
**

      integer i,k
      real*8 hybm_8(nk),prpref,pr1,pibb(nk),pia(nk)
      real*8 conv,fact

      hybrid_to_pres=-1
      print *,'!!WARNING!! function hybrid_to_pres will calculate only'
      print *,'a NORMALIZED (kind=1) hybrid coordinate. '
      print *,'!!RECOMMEND using hyb_to_pres function!!'

      if (rcoef.lt.1.0.or.rcoef.gt.2.0) then
          print *,
     $      'ERROR in hybrid_to_pres: rcoef must be between 1.0 and 2.0'
          return
      endif
      if (pref .lt.400 .or. pref .gt.1050) then
          print *,'ERROR in hybrid_to_pres: pref must be a value',
     $            ' between 400 and 1050'
          return
      endif
      if (ptop .lt.0 .or. ptop .gt.1200) then
          print *,'ERROR in hybrid_to_pres: ptop  must be a value',
     $            ' between 0 and 1200'
          return
      endif
      if (abs(rcoef-1.0).lt.1.0e-5) rcoef = 1.0

      fact = 1.0
*     detect if ps is in millibars or pascals
      if (ps(1).lt.40000.0) then
          conv = 100.0
          fact = fact/100.0
      else
          conv = 1.0
      endif


      do k=1,nk
         if (hyb(k).lt.0 .or. hyb(k).gt. 1.0) then
           print *,'ERROR in hybrid_to_pres: invalid value(s) in',
     $             ' hybrid coordinate array'
           return
         endif
         hybm_8(k)= hyb(k) + (1-hyb(k)) * ptop/pref
      enddo

      prpref = 100.*ptop/hybm_8(1)

      pr1 = 1./(1. - hybm_8(1))
      do k = 1,nk
         pibb(k)  = ((hybm_8(k) - hybm_8(1))*pr1 ) ** rcoef
         pia(k)  = prpref * ( hybm_8(k) - pibb(k) )
      enddo

      do k=1,nk
         pibb(k) = pibb(k)*conv
         do i=1,ni*nj
            pressure(i,k) = (pia(k)+pibb(k)*ps(i)) * fact
         enddo
         hybm(k) = hybm_8(k)
      enddo 
      print *,'WARNING: recommend to use hyb_to_pres'
      hybrid_to_pres=0
      return
      end
*
**********************************************************************
***function hybref_to_ig hybrid coordinate markers coding
**********************************************************************
*
      integer function hybref_to_ig(ig1,ig2,ig3,ig4,rcoef,pref,x1,x2)

      implicit none
      integer ig1,ig2,ig3,ig4
      real rcoef,pref,x1,x2
*
*author
*     Vivian Lee/Michel Valin    Nov.28, 2001
*
*revision
*
*object
*     To derive ig? values given hybrid reference values
*
*arguments
* ________________________________________________________________
*  Name        I/O      Description
* ----------------------------------------------------------------
* ig1           O       ig1 value
* ig2           O       ig2 value
* ig3           O       ig3 value
* ig4           O       ig4 value
* pref          I       reference pressure (normally = 800mb)
* rcoef         I       coefficient (1.0 to 2.0)
* x1            I       not used
* x2            I       not used
* ----------------------------------------------------------------
* the function will return 0 upon success, -1 if there is an error
**
      hybref_to_ig=-1
      if (pref.lt.400.0.or.pref.gt.1050.0) then
          print *,
     %    'ERROR in hybref_to_ig: pref must be between 400.0 and 1050.0'
          return
      endif
      if (rcoef.lt.1.0.or.rcoef.gt.2.0) then
          print *,
     %    'ERROR in hybref_to_ig: rcoef must be between 1.0 and 2.0'
          return
      endif
      ig1 = pref
      ig2 = rcoef*1000.0
      ig3 = 0
      ig4 = 0
      hybref_to_ig=0
      return
      end
*
**********************************************************************
***function ig_to_hybref hybrid coordinate markers decoding
**********************************************************************
*
      integer function ig_to_hybref(ig1,ig2,ig3,ig4,rcoef,pref,x1,x2)

      implicit none
      integer ig1,ig2,ig3,ig4
      real ptop,rcoef,pref,x1,x2
*
*author
*     Vivian Lee/Michel Valin    Nov.28, 2001
*
*revision
*
*object
*     To derive hybrid reference values given ig? values
*
*arguments
* ________________________________________________________________
*  Name        I/O      Description
* ----------------------------------------------------------------
* ig1           I       ig1 value
* ig2           I       ig2 value
* ig3           I       ig3 value
* ig4           I       ig4 value
* pref          O       reference pressure (normally = 800mb)
* rcoef         O       coefficient (1.0 to 2.0)
* x1            I       not used
* x2            I       not used
* ----------------------------------------------------------------
* the function will return 0 upon success, -1 if there is an error
**
      ig_to_hybref=-1
      pref = ig1
      rcoef = ig2/1000.0
      if (pref.lt.400.0.or.pref.gt.1050.0) then
          print *,
     %    'ERROR in ig_to_hybref: ig1 must be between 400 and 1050'
          return
      endif
      if (rcoef.lt.1.0.or.rcoef.gt.2.0) then
          print *,
     %    'ERROR in ig_to_hybref: rcoef(ig2/1000) must be',
     %    ' between 1.0 and 2.0'
          return
      endif
      if (ig3.ne.0) then
          print *,'ERROR in ig_to_hybref: ig3 must be 0'
          return
      endif
      if (ig4.ne.0) then
          print *,'ERROR in ig_to_hybref: ig4 must be 0'
          return
      endif
      ig_to_hybref=0
      return
      end
*
**********************************************************************
***function read_decode_hyb - to read a FSTD record and then decode the
*                                reference values for a hybrid coordinate
**********************************************************************
*
      integer function read_decode_hyb(iun,nom,ip2,ip3,etik,
     %                                            date,ptop,pref,rcoef)

      implicit none
      integer iun,ip2,ip3,date
      real rcoef,pref,ptop
      character*(*) nom
      character*(*) etik
*
*author
*     Vivian Lee/Michel Valin    Dec.19, 2001
*
*revision
*     Vivian Lee BUG in return code Dec 03 2007
*
*object
*     To derive hybrid reference values given a selected FSTD record
*
*arguments
* ________________________________________________________________
*  Name        I/O      Description
* ----------------------------------------------------------------
* iun           I       unit number for input file to read from
* nom           I       variable name for search
* ip2           I       ip2 value for search reference
* ip3           I       ip3 value for search reference
* etik          I       label for search reference
* datev         I       date of validity for search reference
* ptop          O       surface pressure at the top
* pref          O       reference pressure in mb(normally = 800mb)
* rcoef         O       coefficient (1.0 to 2.0)
* ----------------------------------------------------------------
* the function will return KEY of FSTD record upon success, a negative error code if there is an error
*
      integer  fstinf,fstprm,ig_to_hybref
      external fstinf,fstprm,ig_to_hybref


      integer  l,deet,ip1a, ip2a, ip3a, ig1a, ig2a, ig3a, ig4a, bit
      integer  idayo, dty,  swa,  lng,  dlf,  ubc,  ex1,  ex2, ex3
      integer  npas, nia, nja, i, j, k,ierr,kind
      real     x1,x2
      character*1 typ,grda,blk_S
      character*(4) var
      character*(12) labanl

*     typvar of HY must be X 
      l = fstinf(iun, i, j, k, date, etik, -1, ip2, ip3, 'X', nom)
      read_decode_hyb = l   !!! BUG FIX read_decode_hyb now properly set
      if (l.ge.0) then
          ierr= fstprm ( l, idayo, deet, npas, nia, nja, k, bit, dty,
     x                 ip1a, ip2a, ip3a, typ, var, labanl,grda,
     $                 ig1a,ig2a,ig3a,ig4a,swa,lng,dlf,ubc,ex1,ex2,ex3 )
          call convip_plus(ip1a,ptop,kind,-1,blk_S,.false.)
          if (ptop.lt.0.0.or.ptop.gt.1200..or. kind.ne.2) then
              print *, 'ERROR in read_decode_hyb: Decoding of ip1 in ',
     $                  nom
              read_decode_hyb = -1
          endif
          ierr=ig_to_hybref(ig1a,ig2a,ig3a,ig4a,rcoef,pref,x1,x2)
          if (ierr.lt.0) then
              print *, 'ERROR in read_decode_hyb: Decoding of ig?? in ',
     $                  nom
              read_decode_hyb = -1
          endif
*         read_decode_hyb = l   !!! BUG , misplaced statement
      else
          print *,'ERROR in read_decode_hyb:Record ',nom,
     $            ' of typvar X is not found'
          read_decode_hyb = -2  !!! SPECIFIC error code if record not found
      endif
      return
      end
*
**********************************************************************
***function write_encode_hyb - to encode hybrid reference values into the
*                          IG* parameters and then write out the FSTD record
**********************************************************************
*
      integer function write_encode_hyb(iun,nom,ip2,ip3,etik,date,
     %                                               ptop,pref,rcoef)

      implicit none
      integer iun,ip2,ip3,date
      real rcoef,pref,ptop
      character*(*) nom
      character*(*) etik
*
*author
*     Vivian Lee/Michel Valin    Dec.19, 2001
*
*revision
*
*object
*     To encode the given hybrid reference values into IG* parameters
*                and write out the FSTD record with given IP*,etik,datev
*
*arguments
* ________________________________________________________________
*  Name        I/O      Description
* ----------------------------------------------------------------
* iun           I       unit number for output file to write to
* nom           I       variable name to be used
* ip2           I       ip2 value to be used
* ip3           I       ip3 value to be used
* etik          I       label to be used
* datev         I       date of validity to be used
* ptop          I       surface pressure at the top(mb)
* pref          I       reference pressure (normally = 800mb)
* rcoef         I       coefficient (1.0 to 2.0)
* ----------------------------------------------------------------
* the function will return 0 upon success, -1 if there is an error
*
      integer  fstecr,hybref_to_ig
      external fstecr,hybref_to_ig


      integer  l,deet,ip1, ig1, ig2, ig3, ig4, bit
      integer  idayo, dty,  swa,  lng,  dlf,  ubc,  ex1,  ex2, ex3
      integer  npas, nia, nja, i, j, k,ierr
      real     x1,x2
      character*1 blk_S

      if (ptop.lt.0.0.or.ptop.gt.1200.) then
          print *, 'ERROR in write_encode_hyb: Encoding of ip1 in ',nom
          write_encode_hyb = -1
          return
      endif
      call convip_plus(ip1,ptop,2,+1,blk_S,.false.)
      ierr = hybref_to_ig(ig1,ig2,ig3,ig4,rcoef,pref,x1,x2)
      if (ierr.lt.0) then
          print *, 'ERROR in write_encode_hyb: Encoding of ig?? in ',nom
          write_encode_hyb = -1
          return
      endif
      x1=ptop
      ierr = fstecr(x1,x2,-32,iun,date,0,0,1,1,1,ip1,ip2,ip3,
     x              'X',nom,etik,'X',ig1,ig2,ig3,ig4,5,.true.)
      write_encode_hyb=ierr
      return
      end
*
**********************************************************************
***function write_bin_hyb- to write hybrid reference values into a
*                          binary file
**********************************************************************
*
      integer function write_bin_hyb(iun,nom,ip2,ip3,etik,datev,
     %                                               ptop,pref,rcoef)

      implicit none
      integer iun,ip2,ip3,datev
      real rcoef,pref,ptop
      character*4 nom
      character*12 etik
*
*author
*     Vivian Lee    Oct.19, 2009
*
*revision
*
*object
*     To write out the given hybrid reference values into a binary file
*
*arguments
* ________________________________________________________________
*  Name        I/O      Description
* ----------------------------------------------------------------
* iun           I       unit number for output file to write to
* nom           I       variable name to be used
* ip2           I       ip2 value to be used
* ip3           I       ip3 value to be used
* etik          I       label to be used
* datev         I       date of validity to be used
* ptop          I       surface pressure at the top(mb)
* pref          I       reference pressure (normally = 800mb)
* rcoef         I       coefficient (1.0 to 2.0)
* ----------------------------------------------------------------
* the function will return 0 upon success, -1 if there is an error
*
      write_bin_hyb=0
      if (ptop.lt.0.0.or.ptop.gt.1200.) then
          print *, 'ERROR in write_encode_hyb: ptop out of range',nom
          write_bin_hyb = -1
          return
      endif
      if (pref.lt.400.0.or.pref.gt.1050.0) then
          print *,
     %   'ERROR in write_bin_hyb: pref must be between 400.0 and 1050.0'
          write_bin_hyb = -1
          return
      endif
      if (rcoef.lt.1.0.or.rcoef.gt.2.0) then
          print *,
     %    'ERROR in write_bin_hyb: rcoef must be between 1.0 and 2.0'
          write_bin_hyb = -1
          return
      endif
      write(iun)nom,etik,ip2,ip3,datev,ptop,rcoef,pref
      return
      end
*
**********************************************************************
***function read_bin_hyb- to read hybrid reference values from a
*                          binary file
**********************************************************************
*
      integer function read_bin_hyb(iun,nom,ip2,ip3,etik,datev,
     %                                               ptop,pref,rcoef)

      implicit none
      integer iun,ip2,ip3,datev
      real rcoef,pref,ptop
      character*4 nom
      character*12 etik
*
*author
*     Vivian Lee    Oct.19, 2009
*
*revision
*
*object
*     To read the hybrid reference values from a binary file
*
*arguments
* ________________________________________________________________
*  Name        I/O      Description
* ----------------------------------------------------------------
* iun           I       unit number for output file to write to
* nom           I       variable name to be used
* ip2           I       ip2 value to be used
* ip3           I       ip3 value to be used
* etik          I       label to be used
* datev         I       date of validity to be used
* ptop          I       surface pressure at the top(mb)
* pref          I       reference pressure (normally = 800mb)
* rcoef         I       coefficient (1.0 to 2.0)
* ----------------------------------------------------------------
* the function will return 0 upon success, -1 if there is an error
*
      read_bin_hyb=0
      read(iun)nom,etik,ip2,ip3,datev,ptop,rcoef,pref
      if (ptop.lt.0.0.or.ptop.gt.1200.) then
          read_bin_hyb = -1
          print *, 'ERROR in read_bin_hyb: ptop out of range',nom
          return
      endif
      if (pref.lt.400.0.or.pref.gt.1050.0) then
          print *,
     %    'ERROR in read_bin_hyb: pref must be between 400.0 and 1050.0'
          read_bin_hyb = -1
          return
      endif
      if (rcoef.lt.1.0.or.rcoef.gt.2.0) then
          print *,
     %    'ERROR in read_bin_hyb: rcoef must be between 1.0 and 2.0'
          read_bin_hyb = -1
          return
      endif
      return
      end
*
**********************************************************************
***function hyb_to_pres  - new function to convert from hybrid to pressure(mb)
*                          which includes the kind value
**********************************************************************
*
      integer function hyb_to_pres(pressure,hyb,ptop,rcoef,pref,
     %                             kind,ps,NI,NJ,NK) 
      implicit none
      integer NI,NJ,NK,kind
      real pressure(NI*NJ,NK),ptop,ps(NI*NJ),rcoef,pref,hyb(NK)
*
*author
*     Vivian Lee    Oct.16, 2009
*
*revision
*     Jeff Blezius  Jan.13, 2010
*      pibb(k)=(hybm_8(k) - ptop/pref/(1-ptop/pref))**rcoef
*      instead of
*      pibb(k)=(hybm_8(k) - hybm_8(1)*ptop/hybm_8(1))**rcoef
*      in case hybm_8(1) is not given!
*
*object
*     To derive pressure fields from levels derived by "convip_plus", the
*     kind value from the ip1 codes and the hybrid reference 
*     parameters(ptop,rcoef,pref)
*
*arguments
* ________________________________________________________________
*  Name        I/O      Description
* ----------------------------------------------------------------
* pressure      O       array of pressure levels (mb)
* hyb           I       array of user-defined hybrid levels (0.0 to 1.0)
* ptop          I       average pressure at the top (mb)
* rcoef         I       coefficient (1.0 to 2.0)
* pref          I       reference pressure in mb (normally = 800mb)
* kind          I       1=normalized, 5=unnormalized
* ps            I       2D pressure at the surface (mb)
* NI            I       Ni dimension of field
* NJ            I       Nj dimension of field
* NK            I       number of levels in hyb
* ----------------------------------------------------------------
* the function will return 0 upon success, -1 if there is an error
**

      integer i,k
      real*8 hybm_8(nk),pr1,pibb(nk),pia(nk)

      hyb_to_pres=-1

      if (kind.eq.1) then
          do k=1,nk
             hybm_8(k)= hyb(k) + (1.-hyb(k)) * ptop/pref
          enddo
      else if (kind.eq.5) then
          do k=1,nk
             hybm_8(k)= hyb(k)
          enddo
      else
          print *,'ERROR in hyb_to_pres: kind is not 1 nor 5'
          return
      endif

      pr1 = 1./(1. - ptop/pref)
      do k = 1,nk
         pibb(k)  = (dmax1(hybm_8(k) - ptop/pref,0.0d0)*pr1 ) ** rcoef
         pia(k)  = pref * ( hybm_8(k) - pibb(k) )
      enddo

      do k=1,nk
         do i=1,ni*nj
            pressure(i,k) = pia(k)+pibb(k)*ps(i)
         enddo
      enddo 
      hyb_to_pres=0
      return
      end
*
**********************************************************************
***function eta_to_pres  - function to convert from eta to pressure(mb)
**********************************************************************
*
      integer function eta_to_pres(pressure,hybm,ptop,ps,NI,NJ,NK)

      implicit none
      integer NI,NJ,NK
      real pressure(NI*NJ,NK),ptop,ps(NI*NJ),hybm(NK)
*
*author
*     Vivian Lee    Oct.16, 2009
*
*revision
*
*object
*     To derive pressure fields from model eta levels (levels used by
*     the model) and parameter(ptop) and p0
*arguments
* ________________________________________________________________
*  Name        I/O      Description
* ----------------------------------------------------------------
* pressure      O       array of pressure levels (mb)
* hybm          I       array of model eta levels (0.0 to 1.0)
* ptop          I       average pressure at the top (mb)
* ps            I       2D pressure at the surface (mb)
* NI            I       Ni dimension of field
* NJ            I       Nj dimension of field
* NK            I       number of levels in hybm
* ----------------------------------------------------------------
* the function will return 0 upon success, -1 if there is an error
**

      integer i,k
      real*8 pibb(nk),pia(nk)

      eta_to_pres=-1

      do k=1,nk
         pibb(k)  = hybm(k)
         pia(k)  = ptop* ( 1.0d0 -hybm(k))
      enddo

      do k=1,nk
         do i=1,ni*nj
            pressure(i,k) = pia(k)+pibb(k)*ps(i)
         enddo
      enddo 
      eta_to_pres=0
      return
      end
*
**********************************************************************
***function etasef_to_pres  - function to convert from eta in
*                             model SEF to pressure(mb)
**********************************************************************
*
      integer function etasef_to_pres(pressure,hybm,ptop,etatop,ps,
     $                                NI,NJ,NK)

      implicit none
      integer NI,NJ,NK
      real pressure(NI*NJ,NK),ptop,ps(NI*NJ),hybm(NK),etatop
*
*author
*     Vivian Lee    Oct.16, 2009
*
*revision
*
*object
*     To derive pressure fields from SEF model levels (value from
*     IP1 using convip_plus), the ptop (PT), the etatop (E1) and p0
*
*arguments
* ________________________________________________________________
*  Name        I/O      Description
* ----------------------------------------------------------------
* pressure      O       array of pressure levels (mb)
* hybm          I       array of model eta levels (0.0 to 1.0)
* ptop          I       average pressure at the top (mb)
* etatop        I       eta at top of model (0.0 to 1.0)
* ps            I       2D pressure at the surface (mb)
* NI            I       Ni dimension of field
* NJ            I       Nj dimension of field
* NK            I       number of levels in hybm
* ----------------------------------------------------------------
* the function will return 0 upon success, -1 if there is an error
**

      integer i,k
      real*8 pibb(nk),pia(nk),eta1

      etasef_to_pres=-1

      eta1 = 1./(1.-etatop)
      do k=1,nk
         pibb(k)  = (hybm(k)- etatop)*eta1
         pia(k)  = ptop* ( 1.0d0 -pibb(k))
      enddo

      do k=1,nk
         do i=1,ni*nj
            pressure(i,k) = pia(k)+pibb(k)*ps(i)
         enddo
      enddo 
      etasef_to_pres=0
      return
      end
*
**********************************************************************
***function sigma_to_pres  - function to convert from sigma to pressure(mb)
**********************************************************************
*
      integer function sigma_to_pres(pressure,hybm,ps,NI,NJ,NK)

      implicit none
      integer NI,NJ,NK
      real pressure(NI*NJ,NK),ps(NI*NJ),hybm(NK)
*
*author
*     Vivian Lee    Oct.16, 2009
*
*revision
*
*object
*     To derive pressure fields from model sigma levels 
*
*arguments
* ________________________________________________________________
*  Name        I/O      Description
* ----------------------------------------------------------------
* pressure      O       array of pressure levels (mb)
* hybm          I       array of model sigma levels (0.0 to 1.0)
* ps            I       2D pressure at the surface (mb)
* NI            I       Ni dimension of field
* NJ            I       Nj dimension of field
* NK            I       number of levels in hybm
* ----------------------------------------------------------------
* the function will return 0 upon success, -1 if there is an error
**

      integer i,k

      sigma_to_pres=-1
      do k=1,nk
         do i=1,ni*nj
            pressure(i,k) = hybm(k)*ps(i)
         enddo
      enddo 
      sigma_to_pres=0
      return
      end
