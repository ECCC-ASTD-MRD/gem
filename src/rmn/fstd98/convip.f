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
      SUBROUTINE CONVIP( ip, p, kind, mode, string, flag )
      implicit none
      integer ip, kind, mode
      real p
      character *(*) string 
      logical flag

**********************************************************************
*     Codage/Decodage P de/a IP pour IP1, IP2, IP3
*     necessaire avant de lire/ecrire un enregistrement
*     sur un fichier standard.
*
*     Etendu des valeurs encodes: 10e-5 -> 10e10
*     1024x1024-1 = 1048575    1048001 -> 1048575 non utilise
*
*     Auteurs: N. Ek et B. Dugas - Mars 1996
*     Revision 001  M. Lepine - juin 1997 convpr devient convip
*     Revision 002  M. Valin  - mai  1998 fichiers std 98
*     Revision 003  B. Dugas  - juillet 2000 code arbitraire 
*     Revision 004  M. Lepine - fevrier 2002 kind = 4, hauteur au sol +
*                               possibilite de forcer newstyle ou non avec mode=2 et mode=3
*     Revision 005  M. Lepine - avril 2002 kind = 5 (hybride), kind = 21 (GalChen)
*                               valeur min, max, zero et facteur multiplicatif
*     Revision 006  M. Lepine - Juin 2002 kind = 6 (Theta)
*     Revision 007  M. Lepine - Oct 2003 kind = 10 (temps en heure)
*     Revision 008  M. Lepine - Dec 2005 kind = 17 (indice de matrice de niveaux)
*     Revision 009  M. Valin  - Mars 2008 kind = 21 (metres pression remplacant GalChen)
*                               introduction de zero_val2 pour la conversion ip->p
*     Revision 010  M. Lepine - Mai 2010 traitement des valeurs en dehors des intervals connus
*                               comme valeurs arbitraires
*
*     Input:    MODE = -1, de IP -->  P
*               MODE =  0, forcer conversion pour ip a 31 bits
*                          (default = ip a 15 bits)
*                          (appel d'initialisation)
*               MODE = +1, de P  --> IP
*               MODE = +2, de P  --> IP en mode NEWSTYLE force a true
*               MODE = +3, de P  --> IP en mode NEWSTYLE force a false
*               FLAG = .true. , ecriture de P avec format dans string
*
*     Input/
*     Ouput:    IP  =   Valeur codee 
*               P    =   Valeur reelle
*               KIND =0, p est en hauteur (m) par rapport au niveau de la mer (-20,000 -> 100,000)
*               KIND =1, p est en sigma                                       (0.0 -> 1.0)
*               KIND =2, p est en pression (mb)                               (0 -> 1100)
*               KIND =3, p est un code arbitraire                             (-4.8e8 -> 1.0e10)
*               KIND =4, p est en hauteur (M) par rapport au niveau du sol    (-20,000 -> 100,000)
*               KIND =5, p est en coordonnee hybride                          (0.0 -> 1.0)
*               KIND =6, p est en coordonnee theta                            (1 -> 200,000)
*               KIND =10, p represente le temps en heure                      (0.0 -> 1.0e10)
*               KIND =15, reserve (entiers)                                   
*               KIND =17, p represente l'indice x de la matrice de conversion (1.0 -> 1.0e10)
*               KIND =21, p est en metres-pression  (partage avec kind=5 a cause du range exclusif)
*                                                                             (0 -> 1,000,000) fact=1e4
*               STRING = valeur de P formattee
**********************************************************************
      real *8 TEN
      parameter (TEN=10.0)
      real *8 limit1, limit2, temp
      real abs_p
      integer iexp,  offset, itemp, lstring
      character *128 var_fmt

      INTEGER, PARAMETER :: Max_Kind = 31
      integer maxkind
      logical NEWSTYLE, NEWENCODING
      real *8 exptab(0:15)
      character *2 kinds(0:Max_Kind)

      INTEGER :: i

      LOGICAL :: validkind(0:Max_Kind) =
     %  (/ (.true.,i=0,6), (.false.,i=7,9), .true., (.false.,i=11,16),
     %     .true., (.false., i=18,20), .true., (.false., i=22,31) /)

      REAL, PARAMETER, DIMENSION(0:Max_Kind) :: low_val =
     %   (/ -20 000., 0., 0.,    -4.8e+8, -20 000., 0.,
     %      1.0, (-4.8e+8,i=7,9), 0.0, (-4.8e+8,i=11,16), 
     %      1.0, (-4.8e+8,i=18,20), 0., (-4.8e+8,i=22,31) /)
      REAL, PARAMETER :: hi_val(0:Max_Kind) =
     %   (/  100 000., 1., 1100., 1.0e+10, 100 000., 1.,
     %       200 000., (1.0e+10,i=7,9), 1.0e+10, (1.0e+10,i=11,16),
     %       1.0e+10, (1.0e+10,i=18,20), 1000000., (1.0e+10,i=22,31) /)
      REAL, PARAMETER :: zero_val(0:Max_Kind) =
     %   (/ 0., 0., 0., 0., 0., 0., 1., (0.0,i=7,16),
     %      1.0, (0.0,i=18,20), 1.001e-4, (0.0,i=22,31) /)
      REAL, PARAMETER :: zero_val2(0:Max_Kind) =
     %   (/ 0., 0., 0., 0., 0., 0., 1., (0.0,i=7,16),
     %      1.0, (0.0,i=18,20), 0.0, (0.0,i=22,31) /)
      REAL, PARAMETER :: fact_val(0:Max_Kind) =
     %   (/ 1., 1., 1., 1., 1., 1., 1., (1.0,i=7,16),
     %      -1.0, (1.0,i=18,20), 1.0e+4, (1.0,i=22,31) /)

      save NEWSTYLE, exptab, kinds, maxkind

      external qqexit

      data NEWSTYLE /.false./

      data exptab /0.0001D0, 0.001D0, 0.01D0, 0.1D0, 1.0, 10.0, 100.0,
     %   1000.0, 10000.0, 100000.0, 1000000.0, 10000000.0,
     %   100000000.0, 1000000000.0, 10000000000.0, 100000000000.0 /

!      CHARACTER(len=2), PARAMETER :: kinds(0:Max_Kind) =
!     %     (/ 'm ', 'sg', 'mb', '  ', 'M ', 'hy', 'th', 25*'??' /)
      data kinds
     %     / 'm ', 'sg', 'mb', '  ', 'M ', 'hy', 'th', 14*'??',
     %       'mp', 10*'??' /

      if (mode .eq.0) then
         NEWSTYLE = .true.
         return
      endif
      NEWENCODING = NEWSTYLE
      if (mode .eq. 2) NEWENCODING = .true.
      if (mode .eq. 3) NEWENCODING = .false.
      if ((NEWENCODING) .or. (mode .eq. -1)) then
         maxkind = Max_Kind
      else
         maxkind = 3
      endif
      if(flag)lstring=len(string)
*
      if (mode.gt.0) then
*
*        Conversion P a IP
*
       if ( kind.lt.0 .or. kind.gt.maxkind ) then
	  write(6,6004) kind
	  call qqexit(1)
	  return
       elseif ( .not. validkind(kind) ) then
*          write(6,6007) kind
          ip = -999999
          return
       endif
       if (kind .eq. 2 .and. p .eq. 0.) then
          ip = 0
          return
       endif
       if(NEWENCODING)then
         if (p .lt. low_val(kind) .or. p .gt. hi_val(kind)) then
            write(6,6006) p,low_val(kind),hi_val(kind)
            ip = -999999
            return
         endif
         iexp = 4
         temp = p
         if (abs(temp) .lt. zero_val(kind)) temp = zero_val(kind)
         temp = temp * fact_val(kind)
         if ( temp .ge. 0) then
            limit1 = 1000000.0
            limit2 = 100000.0
            offset = 0
         else
            temp = -temp
            limit1 = 48000.0
            limit2 = 4800.0
            offset = 1000000
         endif
100      if ( iexp .gt. 0 .and. iexp .lt. 15 ) then
c            ip = temp
c            if ( ip .eq. temp .and. temp .le. limit1)
c     %         goto 101
            if (temp .ge. limit1 ) then
               temp = temp / TEN
               iexp = iexp -1
            else if ( temp .lt. limit2 ) then
               temp = temp * TEN
               iexp = iexp + 1
            else
               goto 101
            endif
         goto 100
         endif
101      continue
         if ( temp .gt. limit1 ) then
            ip = -1
         else
            ip = offset + temp + .5
         endif
         ip = ior (ip,ishft(iexp,20))
         ip = ior (ip,ishft(iand(15,kind),24))
       else

         if (kind.eq.0) then

*           ... de hauteur ...

            ip = max( 12001,min( 32000,nint( p/5.0 + 12001 ) ) )

         elseif (kind.eq.1) then

*           ... de sigma ...

            if ( .not. (  0.0 .le. p .and. p .le. 1.0 ) ) then
               write(6,6001) p
               ip = -999999
               return
            endif

            ip = nint( p * 10000. ) + 2000

         elseif (kind.eq.2) then

*           ... de pression ...

            if (  .not. (0.0 .le. p .and. p .lt. 1100. ) ) then
               write(6,6002) p
               ip = -999999
               return
            endif

            if (0.999999e+1 .le. p .and. p .lt. 1100. ) then
               ip = nint ( p )
            elseif ( p .lt. 0.999999e+1 ) then
               if( p .ge. 0.999999e0 ) then
                  ip = 1800 + nint(20.*p)
               elseif ( p .ge. 0.999999e-1 ) then
                  ip = 1600 + nint(200.*p)
               elseif ( p .ge. 0.999999e-2 ) then
                  ip = 1400 + nint(2000.*p)
               elseif ( p .ge. 0.999999e-3 ) then
                  ip = 1200 + nint(20000.*p)
               else
                  ip = 0
               endif
            endif

         elseif (kind.eq.3) then

*           ... ou de code arbitraire

	    ip = nint( p )
            if ( 0 .le. ip .and. ip .le. 100 ) then
               ip = 1200 - ip
            else
               write(6,6003) p
               ip = -999999
               return
            endif

         elseif (kind.eq.10) then

*           ... temps en heure
            ip = nint( p )
            if (ip > 32767) then
               write(6,6005) p
               ip = -999999
               return
            endif
 
         else

*           Valeur illegale de kind.

            write(6,6004) kind
            ip = -999999
	    return
         endif

       endif

      elseif (mode.lt.0) then

*        Conversion de ip ...

         if ( ip .gt. 32767 ) then

*           nouveau codage

            kind = iand(15,ishft(ip,-24))
            if ( kind.lt.0 .or. kind.gt.maxkind ) then
               write(6,6004) kind
               p = -999999.0
               kind = -1
               if(flag) string = 'Invalid'
               return
            elseif ( .not. validkind(kind) ) then
*              write(6,6007) kind
               p = -999999.0
               kind = -1
               if(flag) string = 'Invalid'
            endif
            iexp = iand (15,ishft(ip,-20))
            itemp = iand (1048575, ip)
            if (itemp > 1 000 000) itemp = -(itemp - 1 000 000)
 555        continue
            p = itemp / exptab(iexp)
            p = p / fact_val(kind)                          ! retirer le facteur multiplicatif
            if (kind .eq. 5) then
               if ((p .lt. low_val(kind)) .or. 
     %             (p .gt. hi_val(kind))) then
                  kind = 21         !  (en dehors du range, on essaye code partage avec kind=5)
*                  print *,'Debug+ o.b. kind 5 --> 21'
                  goto 555  
               endif
            endif
            if (kind .eq. 1) then
               if ((p .lt. low_val(kind)) .or. 
     %             (p .gt. hi_val(kind))) then
                  kind = 17         !  (en dehors du range, on essaye code partage avec kind=2)
*                  print *,'Debug+ o.b. kind 1 --> 17'
                  goto 555  
               endif
            endif
            if (p .lt. low_val(kind)) p = low_val(kind)     ! clipping a la valeur minimale
            if (p .gt. hi_val(kind)) p = hi_val(kind)       ! clipping a la valeur maximale
            if (abs(p) .lt. 1.001*zero_val(kind)) p = zero_val2(kind)
            abs_p = abs(p)
            if (flag) then
               if (len(string) .ge. 15) then
                  if(abs_p .eq. int(abs_p) .and. abs_p.lt.1000000.) then
                     write(var_fmt,'(i12,1x,a2)')int(p),kinds(kind)
                  elseif (abs_p.ge. 1000000.) then
                     write(var_fmt,'(e12.6,1x,a2)')p,kinds(kind)
                  elseif (abs_p.ge. 100000.) then
                     write(var_fmt,'(f12.0,1x,a2)')p,kinds(kind)
                  elseif (abs_p.ge. 10000.) then
                     write(var_fmt,'(f12.1,1x,a2)')p,kinds(kind)
                  elseif (abs_p.ge. 1000.) then
                     write(var_fmt,'(f12.2,1x,a2)')p,kinds(kind)
                  elseif (abs_p.ge. 100.) then
                     write(var_fmt,'(f12.3,1x,a2)')p,kinds(kind)
                  elseif (abs_p.ge. 10.) then
                     write(var_fmt,'(f12.4,1x,a2)')p,kinds(kind)
                  elseif (abs_p.ge. 1.) then
                     write(var_fmt,'(f12.5,1x,a2)')p,kinds(kind)
                  elseif (abs_p.ge. 0.1) then
                     write(var_fmt,'(f12.6,1x,a2)')p,kinds(kind)
                  elseif (abs_p.ge. 0.01) then
                     write(var_fmt,'(f12.7,1x,a2)')p,kinds(kind)
                  elseif (abs_p.ge. 0.001) then
                     write(var_fmt,'(f12.8,1x,a2)')p,kinds(kind)
                  else
                     write(var_fmt,'(e12.6,1x,a2)')p,kinds(kind)
                  endif
               else
                  if(abs_p .eq. int(abs_p) .and. abs_p.lt.1000000.) then
                     write(var_fmt,'(i6,1x,a2)')int(abs_p),kinds(kind)
                  elseif (abs_p.gt. 1000000.) then
                     write(var_fmt,'(e9.4,1x,a2)')p,kinds(kind)
                  elseif (abs_p.gt. 100.) then
                     write(var_fmt,'(f6.2,1x,a2)')p,kinds(kind)
                  elseif (abs_p.gt. 10.) then
                     write(var_fmt,'(f6.3,1x,a2)')p,kinds(kind)
                  elseif (abs_p.gt. 1.) then
                     write(var_fmt,'(f6.4,1x,a2)')p,kinds(kind)
                  elseif (abs_p.ge. 0.001) then
                     write(var_fmt,'(f6.5,1x,a2)')p,kinds(kind)
                  else
                     write(var_fmt,'(e9.4,1x,a2)')p,kinds(kind)
                  endif
               endif
               string=var_fmt
            endif

         elseif (  12000 .lt. ip .and. ip .le. 32000) then

*           ... a hauteur ...

            kind = 0
            p = 5 * ( ip -12001)
            if (flag) write(string,'(i6,1x,a1)') nint(p),'m'


         elseif (  2000 .le. ip .and. ip .le. 12000 ) then

*           ... a sigma ...

            kind = 1
            p = float (ip - 2000) / 10000.
            if (flag) write(string,'(f6.4,1x,a2)') p,'sg'

         elseif (( 0    .le. ip .and. ip .lt. 1100 )  .or.
     +           ( 1200 .lt. ip .and. ip .lt. 2000 )) then

*           ... a pression ...
     
            kind = 2
            if ( 0 .le. ip .and. ip .lt. 1100 ) then
               p = float(ip)
               if (flag) write(string,'(i6,1x,a2)') int(p),'mb'
            elseif ( ip .lt. 1400 ) then
                  p = float(ip-1200) / 20000.D0
                  if (flag) write(string,'(f6.5,1x,a2)') p,'mb'
            elseif ( ip .lt. 1600) then
                  p = float(ip-1400) / 2000.D0
                  if (flag) write(string,'(f6.4,1x,a2)') p,'mb'
            elseif ( ip .lt. 1800) then
                  p = float(ip-1600) / 200.D0
                  if (flag) write(string,'(f6.3,1x,a2)') p,'mb'
            elseif  ( ip .lt. 2000) then
                  p = float(ip-1800) / 20.D0
                  if (flag) write(string,'(f6.2,1x,a2)') p,'mb'
            endif

         elseif ( 1100 .le. ip .and. ip .le. 1200) then

*           ... ou a code arbitraire.

            kind = 3
            p = float( ip )
            p = 1200. - p
            if (flag) write(string,'(i6,3x)') nint(p)

         else

*           Valeur inderminee de ip.
            kind = 3
            p = float( ip )
            if (flag) write(string,'(i6,3x)') nint(p)

         endif

      endif
      
      return

**********************************************************************
 6001 format(' Error in convip: sigma value =',e10.5,
     %       ' returned ip is -999999')
 6002 format(' Error in convip: pressure value =',e10.5,
     %       ' returned ip is -999999')
 6003 format(' Error in convip: arbitrairy value=',e10.5,
     %       ' returned ip is -999999')
 6004 format(' Error in convip: invalid kind =',I10)
 6005 format(' Error in convip: kind=10 (oldstyle) value out of range='
     %       ,e10.5,' returned ip is -999999')
 6006 format(' Error in convip: p is out of bounds =',e10.5,' min=',
     %       e10.5,' max=',e10.5,' returned ip is -999999')
 6007 format(' Warning in convip: undetermined kind used =',I10)

      end
      
      subroutine igapg(grtyp,pg1,pg2,pg3,pg4,ig1,ig2,ig3,ig4)
      implicit none
      character *1 grtyp
      character *(*) pg1,pg2,pg3,pg4
      integer ig1,ig2,ig3,ig4
*
      real xg1,xg2,xg3,xg4
      external cigaxg
*
      if ((grtyp .eq. 'Y') .or. (grtyp .eq. 'Z')
     %    .or. (grtyp .eq. '!')) then
         write(pg1,'(i6)') ig1
         write(pg2,'(i6)') ig2
         write(pg3,'(i7)') ig3
         write(pg4,'(i7)') ig4
         return
      else if ((grtyp .eq. 'A') .or. (grtyp .eq. 'B') .or.
     %         (grtyp .eq. 'G') .or. (grtyp .eq. 'N') .or.
     %         (grtyp .eq. 'S') .or. (grtyp .eq. 'E')) then
         call cigaxg(grtyp,xg1,xg2,xg3,xg4,ig1,ig2,ig3,ig4)
      else
         write(pg1,'(i6)') ig1
         write(pg2,'(i6)') ig2
         write(pg3,'(i7)') ig3
         write(pg4,'(i7)') ig4
         return
      endif
      if ((grtyp .eq. 'N') .or. (grtyp .eq. 'S')) then
         write(pg1,'(f6.1)') xg1
         write(pg2,'(f6.1)') xg2
         write(pg3,'(f5.1,a2)') (xg3/1000.), 'Km'
         write(pg4,'(f7.3)') xg4
      else if ((grtyp .eq. 'A') .or. (grtyp .eq. 'B') .or.
     %         (grtyp .eq. 'G')) then
         write(pg1,'(i6)') int(xg1)
         write(pg2,'(i6)') int(xg2)
         write(pg3,'(i7)') int(xg3)
         write(pg4,'(i7)') int(xg4)
      else
         write(pg1,'(f6.2)') xg1
         write(pg2,'(f6.2)') xg2
         write(pg3,'(f7.3)') xg3
         write(pg4,'(f7.3)') xg4
      endif
      return
      end
