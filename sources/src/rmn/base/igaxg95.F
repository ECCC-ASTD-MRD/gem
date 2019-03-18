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
      integer function igagds95(gds,gtyin,ig1,ig2,ig3,ig4)
      implicit none
      character*1 gtyin,gtyout
      character*256 gds
      integer ig1,ig2,ig3,ig4

      integer qqqgbld,qgblukig
      external qqqgbld,qgblukig

      real xglst(20)
      integer nxg

      integer iun,ier
#include "qqqcgb1.cdk"      

      nxg=20
      if (gtyin.ne.'!') then
         print 80
         print *, ' <igagds95> Error: bad grid type ', gtyin
         return
      else
         ier = qqqgbld()
         if (ier.lt.0) then
            igagds95 = -1
            return
         endif
         
         ier = qgblukig(gtyout,xglst,nxg,gtyin,ig1,ig2,ig3,ig4)
         if (ier.lt.0) then
            igagds95 = -1
            return
         else
            gds=gdsinfo(ier)
            igagds95 = nint(xglst(14))
         endif
      endif
      
 80   format(' ',79('*'))

      return
      end


      subroutine igaxg95(grtypout,xglst,nxg,grtypin,ig1,ig2,ig3,ig4)
      implicit none
      character*1 grtypout,grtypin
      integer nxg
      real xglst(nxg)
      integer ig1,ig2,ig3,ig4

      integer qqqgbld,qgblukig
      external qqqgbld,qgblukig

      integer iun,ier
      
      grtypout= '*'
      if (grtypin.ne.'!') then
         if (nxg.lt.4) then
            print 80
            print *, ' <igaxg95> Error: xglist not long enough'
            print *, 
     $           ' <igaxg95> need 4 elements, list given has only ',nxg
            return
         endif
         call cigaxg(grtypin,xglst(1),xglst(2),xglst(3),xglst(4),
     $        ig1,ig2,ig3,ig4)
         grtypout= grtypin
         return
      else
         ier = qqqgbld()
         if (ier.lt.0) then
            return
         endif
            
         ier = qgblukig(grtypout,xglst,nxg,grtypin,ig1,ig2,ig3,ig4)
         ig3 = 256*nint(xglst(13))+nint(xglst(14))
         ig4 = 0
         if (ier.lt.0) then
            return
         endif
         
      endif
      
 80   format(' ',79('*'))

      return
      end


      integer function qgblukig(gtyout,xglst,nxg,gtyin,ig1,ig2,ig3,ig4)
      implicit none
      integer nxg
      character*1 gtyin,gtyout
      real xglst(nxg)
      integer ig1,ig2,ig3,ig4


#include "qqqcgb1.cdk"
#include "qqqgb2.cdk"
      
      integer i,gribcode,ni,nj,centercode,subcentercode,projcode
      real xlat00,xlon00,xlatninj,xlonninj,dx,dy,
     $     yaxislon,latin1,latin2
      logical found

      found = .false.

      gtyout = '*'
      i = 1
 10   if (i.gt.qqqngb.and..not.found) goto 999
      read (cgbinfo(i),*,err=110) gribcode
      if (gribcode.eq.ig1) then
         if (nxg.lt.14) then
            print 80
            print *, ' <qgblukig> Error: xglist not long enough'
            print *, 
     $        ' <qgblukig> need 14 elements, list given has only', nxg
            qgblukig = -1
            return
         endif
         read (cgbinfo(i),*,err=110) gribcode,
     $        centercode,subcentercode,projcode,gtyout, 
     $        xlat00,xlon00,dx,dy,
     $        yaxislon,latin1,latin2,ni,nj,xlatninj,xlonninj,gdsinfo(i)
         
         
         xglst( 1)  = xlat00
         xglst( 2)  = xlon00
         xglst( 3)  = dx
         xglst( 4)  = dy
         xglst( 5)  = yaxislon
         xglst( 6)  = latin1
         xglst( 7)  = latin2
         xglst( 8)  = ni*1.0
         xglst( 9)  = nj*1.0
         xglst(10)  = xlatninj
         xglst(11)  = xlonninj
         xglst(12)  = centercode * 1.0
         xglst(13)  = subcentercode * 1.0
         xglst(14)  = projcode *1.0
         
         if (ig2.eq.nint(xglst(12)).or.ig2.eq.0) then
            found = .true.
         endif
      else
         i = i+1
         goto 10
      endif
      
      qgblukig = i 
      return
      
 80   format(' ',79('*'))
 256  format(256(' '))
 110  print 80
      print *,' <qgblukig> read error encountered in file "gribtable"'
      print *,
     $    ' <qgblukig> contents of file don''t match field description'
      goto 999

 999  qgblukig = -1
      print *,' <qgblukig> GRIB code ',ig1,'not found...'
      

      return
      end

      
      subroutine qqqcltab(str,strlen)
      integer strlen
      character str(strlen)

      integer i

      do i=1,strlen
         if (str(i).lt.' ') str(i) = ' '
      enddo

      return
      end

      integer function qqqgbld()
      implicit none
#include "qqqcgb1.cdk"      
#include "qqqgb2.cdk"

      character*256 filename,armnlib
      character*512 buffer
      external fnom,fclos,qqqlen
      integer fnom
      integer ier,strlen,iun,qqqlen

      logical once
      data once /.false./
      save once

      data qqqngb /0/


      if (once) then
         qqqgbld = 1
         return
      endif
      
      filename='gribtable'
      iun = 0
      ier = fnom(iun, filename, 'FTN+FMT+R/O+OLD', 0)
      if (ier.ne.0) then
        print *,' <qqqgbld> no local "gribtable" file found'
        call fclos(iun)
        call getenvc('AFSISIO',armnlib)
        strlen = qqqlen(armnlib,len(armnlib))
        print *,len(armnlib)
        filename = armnlib(1:strlen)//'/datafiles/constants/'//filename
        ier = fnom(iun, filename, 'FTN+FMT+R/O+OLD', 0)
        if (ier.ne.0) then
          call getenvc('ARMNLIB',armnlib)
          strlen = qqqlen(armnlib,len(armnlib))
          print *,len(armnlib)
          filename = armnlib(1:strlen)//'/data/'//filename
          ier = fnom(iun, filename, 'FTN+FMT+R/O+OLD', 0)
          if (ier.ne.0) then
            print 80
            print *,' <qqqgbld> no "gribtable" file found'
            print *,' <qqqgbld> no "./gribtable" and  no "'
     &           ,filename(1:qqqlen(filename,len(filename))),'"'
            call fclos(iun)
            print *,' <qqqgbld> using grib table in '
     &           ,filename(1:qqqlen(filename,len(filename)))
          endif
        endif
      endif
      
      rewind(iun)
 20   read(iun,256,err=110,end=120) buffer
      if (buffer(1:1).ne.'#') then
         call qqqcltab(buffer,len(buffer))
         qqqngb = qqqngb + 1
         if (qqqngb.gt.ncodemax) goto 111
         cgbinfo(qqqngb) = buffer
      endif
      goto 20
      
 120  print *,' <qqqgbld> end-of-file reached'
      print *,' <qqqgbld> found ',qqqngb, ' grids'
      once = .true.

      call fclos(iun)
      qqqgbld = 0
      return

 110  print 80
      print *,' <qqqgbld> read error encountered in file "gribtable"'
      call fclos(iun)
      qqqgbld = -1

 111  print 80
      print *,' <qqqgbld> internal string table too small'
      print *,' <qqqgbld> modify code to allocate more space'
      call fclos(iun)
      qqqgbld = -1

 80   format(' ',79('*'))
 256  format(a512)

      return
      end


      integer function qqqlen(str,strlen)
      integer strlen
      character str(strlen)

      integer i
      i = strlen

 20   if (str(i).gt.' '.and.i.ge.1) then
         qqqlen = i
      else
         i = i-1
         goto 20
      endif


      return
      end
