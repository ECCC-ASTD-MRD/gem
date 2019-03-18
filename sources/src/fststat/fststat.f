***s/p fststat
*
      subroutine fststat(iun, NI, NJ, NK, datev,
     $     etiket, ip1, ip2, ip3,typvar, nomvar)
      implicit none
      integer iun,ni,nj,nk,ip1,ip2,ip3,datev
      character*12 etiket
      character*2 typvar
      character*4 nomvar
*
*AUTHOR   Yves Chartier                      July 1993
* 
*REVISION 001  M. Lepine - Mars 2002 - ajout du mode reduction 32 pour IEEE a 64 bits
*REVISION 002  M. Lepine - Juin 2003 - remplacement de memoirh
*
*Version 6.2   D. Bouhemhem - Nov. 2012 - Compilation avec librmn_013
*Version 6.3   M. Lepine    - Mars 2014 - Compilation avec librmn_014
*Version 6.4   M. Lepine    - Dec  2014 - Compilation avec librmn_015.1
*Version 6.5   M. Lepine    - Fev  2015 - Compilation avec librmn_015.2
*Version 6.6   M. Lepine    - Jan  2017 - Eviter le traitement des
*      enregistrements '!!' qui contiennent un melange entiers, reels et
*      caracteres 
*
*LANGUAGE:  fortran 77
*
*OBJECT (fststat)
*
*FILES
*     tape1: TSF file
*     tape10-49: RPN standard files 
*
*ARGUMENTS 
*
*IMPLICIT     
*
*MODULES
      external fstinf,fstprm,fstluk,fstsui,fstopl
      integer ier,fstprm,fstinf,fstsui,fstluk,fstopl
      character*1 grtyp
      integer i
      integer key, date0, deet, npas, nbits, datyp
      integer swa, lng, dltf, ubc
      integer ig1, ig2, ig3, ig4, extra1, extra2, extra3
      real rtemp
      integer itemp
      equivalence (rtemp,itemp)
      real, allocatable, dimension(:) :: buf
      character*12 etiket2
      key = fstinf(iun, ni, nj, nk,  datev, etiket,
     $     ip1,ip2,ip3,typvar,nomvar)
      ier = fstopl('REDUCTION32',.true.,.false.)
 10   if (key.ge.0) then
*         call memoirh(buf,ibuf,ni*nj*nk)
         allocate(buf(ni*nj*nk))
         ier = fstprm(key, date0, deet, npas, ni, nj, nk, nbits,
     *        datyp, ip1, ip2, ip3, typvar, nomvar, etiket2, grtyp,
     *        ig1, ig2, ig3, ig4, swa, lng, dltf, ubc,
     *        extra1, extra2, extra3)
         if (nomvar .eq. '!!') then
           WRITE(6,*)' **   SKIPPING RECORD "!!", CAN''T PROCESS  **'
           goto 20
         endif
         ier = fstluk(buf,key,ni,nj,nk)
         if (datyp.eq.2.or.datyp.eq.4) then
            do i=1,ni*nj*nk
               rtemp = buf(i)
               buf(i)=real(itemp)
            end do
         endif
         call statfld4(nomvar,typvar,ip1,ip2,ip3,date0,etiket2,buf,
     $        ni,nj,nk)
 20      key = fstsui(iun,ni,nj,nk)
         deallocate(buf)
*     call memoirh(buf,ibuf,0)
         goto 10
      endif
      return
      end
