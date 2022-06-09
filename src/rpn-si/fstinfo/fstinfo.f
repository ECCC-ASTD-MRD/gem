      program fstinfo
c
! Revision 006 - M. Lepine, augmenter la longueur pour les noms de fichiers a 4k
! Revision 007 - M. Lepine, reload avec librmn_014
! Revision 008 - M. Lepine, reload avec librmn_015.1
! Revision 009 - M. Lepine, reload avec librmn_015.2
      implicit none
c
      integer pnier,pnfstsrc,pnnis,pnnjs,pnnks
      integer fnom,fstouv,fstfrm,fstinf,pnkey,fstsui
c
      integer pndateo,pndatprinto,pndatev,pndatprintv,pntimeo
      integer pntimev,pndeet,pnnpas,pnnbits
      integer pndatyp, pnip1,pnip2,pnip3
      integer pnig1,pnig2,pnig3,pnig4,pnswa,pnlng,pndltf,pnubc,pnextra1
      integer pnextra2,pnextra3
c
      integer fstprm,pnseqout,ilen
      integer wkoffit,ikind
      external qqexit,ccard,wkoffit
c
      real*8 prtmp8
c
      character*1 ptgrtyp,ptdel
      character*2 pttypvar
      character*4 ptnomvar
      character*12 ptetiket
c
      character(len=8) :: cltimev
      character*8 ptcle(12)
      character*4096 ptvar(12), ptdefvar(12)
c
      data ptcle /'izfst.','datev.','vdatev.','etiket.','ip1.','ip2.'
     &     ,'ip3.','typvar.','nomvar.','otxt.','del.','champs.'/
      data ptvar /'bidon','-1','-1',' ','-1','-1','-1',' ',' ','fstlist'
     &     ,':','0'/
      data ptdefvar /'bidon','-1','-1',' ','-1','-1','-1',' ',' '
     &     ,'fstlist',':','1'/
      data pnfstsrc /11/
c
      call ccard(ptcle,ptdefvar,ptvar,12,-1)
c
c     ------ Liste de champs -----
      if (ptvar(12).eq.'1') then
c
         write(*,'(27(/,1x,a))')
     $        '1 nomvar',
     $        '2 typvar',
     $        '3 ip1',
     $        '4 ip2',
     $        '5 ip3',
     $        '6  ni',
     $        '7  nj',
     $        '8  nk',
     $        '9 etiket',
     $        '10  date d origine',
     $        '11 date de validite',
     $        '12 deet',
     $        '13  npas',
     $        '14 grtyp',
     $        '15 ig1',
     $        '16 ig2',
     $        '17 ig3',
     $        '18 ig4',
     $        '19  datyp',
     $        '20  nbits',
     $        '21 swa',
     $        '22 lng',
     $        '23 dltf',
     $        '24 ubc',
     $        '25 extra1',
     $        '26 extra2',
     $        '27 extra3'
         stop
      endif
c
      ikind = wkoffit(ptvar(1))
c      write(*,*) 'IKIND ====== ',ikind
      if(ikind .ne. 33 .and. ikind .ne. 34 .and. ikind .ne. 1) call qqex
     %it(1)
c     ------ Initialisation des clefs de recherche -----
c
      if (ptvar(2) .ne. '-1') then
        read(ptvar(2),*) pndatev
      elseif(ptvar(3) .ne. '-1') then
        read(ptvar(3)(1:8),*) pndatprintv
        read(ptvar(3)(9:),'(a8)') cltimev
        ilen = len_trim(cltimev)
        cltimev = '00000000'
        read(ptvar(3)(9:9+ilen),'(a)') cltimev(1:ilen)
        read(cltimev,*) pntimev
        call newdate(pndatev,pndatprintv,pntimev,3)
      else
        pndatev = -1
      endif
      read(ptvar(5),*) pnip1
      read(ptvar(6),*) pnip2
      read(ptvar(7),*) pnip3
c
      ptetiket=ptvar(4)
      pttypvar=ptvar(8)
      ptnomvar=ptvar(9)
      ptdel=ptvar(11)
c
c     ------ Ouverture des fichiers ------
c
      pnseqout = 21
      open(pnseqout,FILE=ptvar(10),ACCESS='SEQUENTIAL')
c
      pnier =  fnom  (pnfstsrc,ptvar(1),'RND+STD+R/O',0)
      pnier =  fstouv(pnfstsrc,'RND')
c
      pnkey=fstinf(pnfstsrc,pnnis,pnnjs,pnnks,pndatev,
     $     ptetiket, pnip1, pnip2, pnip3, pttypvar,
     $     ptnomvar)
c
 100  if (pnkey.ge.0) then
         pnier=fstprm(pnkey,pndateo,pndeet,pnnpas,pnnis,pnnjs,pnnks,pnnb
     %its,
     $        pndatyp,pnip1,pnip2,pnip3,pttypvar,ptnomvar,ptetiket,ptgrt
     %yp,pnig1,pnig2,
     $        pnig3,pnig4,pnswa,pnlng,pndltf,pnubc,pnextra1,pnextra2,pne
     %xtra3)
c
         prtmp8=pndeet
         prtmp8=prtmp8*pnnpas/3600.
c
         call incdatr(pndatev,pndateo,prtmp8)
c
         call newdate(pndateo,pndatprinto,pntimeo,-3)
         call newdate(pndatev,pndatprintv,pntimev,-3)
c
         write (pnseqout,1000) ptnomvar,ptdel,pttypvar,ptdel,pnip1,ptdel
     &        ,pnip2,ptdel,pnip3,ptdel,pnnis,ptdel,pnnjs,ptdel,pnnks
     &        ,ptdel,ptetiket,ptdel,pndatprinto,pntimeo,ptdel
     &        ,pndatprintv,pntimev,ptdel,pndeet,ptdel,pnnpas,ptdel
     &        ,ptgrtyp,ptdel,pnig1,ptdel,pnig2,ptdel,pnig3,ptdel,pnig4
     &        ,ptdel,pndatyp,ptdel,pnnbits,ptdel,pnswa,ptdel,pnlng,ptdel
     &        ,pndltf,ptdel,pnubc,ptdel,pnextra1,ptdel,pnextra2,ptdel
     &        ,pnextra3
c
         pnkey=fstsui(pnfstsrc,pnnis,pnnjs,pnnks)
         goto 100
       endif
c
      pnier =  fstfrm(pnfstsrc)
      close(pnseqout)
c
 1000 format(2(a,a),6(i10,a),1(a,a),2(i8.8,i8.8,a),2(i10,a),1(a,a)
     &     ,12(i10,a),i10)
c
      call qqexit(0)
      stop
      end
      character *128 function product_id_tag()
      product_id_tag='$Id$'
      return
      end
