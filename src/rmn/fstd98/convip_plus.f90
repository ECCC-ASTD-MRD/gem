!> Encode/Decode IP pour IP1, IP2, IP3
!> This is required before reading or writting to a FST file
!> @date 1996
!> @author Nils Ek
!> @author Bernard Dugas
!> @author Mario Lepine
!> @author Michel Valin
SUBROUTINE CONVIP_plus( ip, p, kind, mode, string, flagv )
  use convert_ip123_int
  implicit none
  include 'convip_plus.inc'

  integer, intent(INOUT) :: ip, kind
  integer, intent(IN) :: mode
  real, intent(INOUT) :: p
  character *(*), intent(OUT) :: string 
  logical, intent(IN) :: flagv

!  Etendue des valeurs reelles encodees: 10e-5 -> 10e10
!  pseudo mantisses :
!  1024x1024-1 = 1048575
!  1048001 -> 1048575 non utilise
!  1000001 -> 1048000 utilise pour valeurs negatives

!     Revision 001  M. Lepine - juin 1997 convpr devient convip
!     Revision 002  M. Valin  - mai  1998 fichiers std 98
!     Revision 003  B. Dugas  - juillet 2000 code arbitraire 
!     Revision 004  M. Lepine - fevrier 2002 kind = 4, hauteur au sol +
!                               possibilite de forcer newstyle ou non avec mode=2 et mode=3
!     Revision 005  M. Lepine - avril 2002 kind = 5 (hybride), kind = 21 (GalChen)
!                               valeur min, max, zero et facteur multiplicatif
!     Revision 006  M. Lepine - Juin 2002 kind = 6 (Theta)
!     Revision 007  M. Lepine - Oct 2003 kind = 10 (temps en heure)
!     Revision 008  M. Lepine - Dec 2005 kind = 17 (indice de matrice de niveaux)
!     Revision 009  M. Valin  - Mars 2008 kind = 21 (metres pression remplacant GalChen)
!                               introduction de zero_val2 pour la conversion ip->p
!     Revision 010  M. Lepine - Mai 2010 traitement des valeurs en dehors des intervalles connus
!                               comme valeurs arbitraires
!     Revision 011  M. Valin  - Mai/Juin 2013 activation du code 15, ajout de la conversion groupee,
!                               menage dans le code, changement de nom, refactoring
!     Revision 012  M. Valin  - Oct/Nov 2013 bug de conversion corrige pour certains cas limites
!                               enleve une amelioration qui entrainait une non compatibilite avec convip

! INPUTS
!    MODE = -1, de IP -->  P
!    MODE =  0, forcer conversion pour ip a 31 bits
!                          (default = ip a 15 bits)
!                          (appel d'initialisation)
!    MODE = +1, de P  --> IP
!    MODE = +2, de P  --> IP en mode NEWSTYLE force a true
!    MODE = +3, de P  --> IP en mode NEWSTYLE force a false
!    FLAGV = .true. , ecriture de P avec format dans string
! INOUTS
!    IP  =   Valeur codee 
!    P    =   Valeur reelle
!    KIND =0, p est en hauteur (m) par rapport au niveau de la mer (-20,000 -> 100,000)
!    KIND =1, p est en sigma                                       (0.0 -> 1.0)
!    KIND =2, p est en pression (mb)                               (0 -> 1100)
!    KIND =3, p est un code arbitraire                             (-4.8e8 -> 1.0e10)
!    KIND =4, p est en hauteur (M) par rapport au niveau du sol    (-20,000 -> 100,000)
!    KIND =5, p est en coordonnee hybride                          (0.0 -> 1.0)
!    KIND =6, p est en coordonnee theta                            (1 -> 200,000)
!    KIND =10, p represente le temps en heure                      (0.0 -> 1.0e10)
!    KIND =15, reserve (entiers)                                   
!    KIND =17, p represente l'indice x de la matrice de conversion (1.0 -> 1.0e10)
!              (partage avec kind=1 a cause du range exclusif
!    KIND =21, p est en metres-pression                            (0 -> 1,000,000) fact=1e4
!              (partage avec kind=5 a cause du range exclusif)                                                             
! OUTPUTS
!    STRING = valeur de P formattee

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
  character (len=12) :: string2
  integer :: status

  INTEGER :: i
  logical :: flag

  LOGICAL, PARAMETER, DIMENSION(0:Max_Kind) :: validkind =                    &
  & (/ (.true.,i=0,6), (.false.,i=7,9), .true., (.false.,i=11,14),            &
  &    .true., .false.,.true.,                                                &
  &    (.false., i=18,20), .true., (.false., i=22,30),.true. /)   ! kind 31 valide

  REAL, PARAMETER, DIMENSION(0:Max_Kind) :: low_val =                         &
  &  (/ -20000., 0., 0.,    -4.8e+8, -20000., 0.,                             &
  &    1.0, (-4.8e+8,i=7,9), 0.0, (-4.8e+8,i=11,16),                          &
  &    1.0, (-4.8e+8,i=18,20), 0., (-4.8e+8,i=22,31) /)
  REAL, PARAMETER, DIMENSION(0:Max_Kind) :: hi_val =                          &
  &  (/  100000., 1., 1100., 1.0e+10, 100000., 1.,                            &
  &     200000., (1.0e+10,i=7,9), 1.0e+10, (1.0e+10,i=11,16),                 &
  &     1.0e+10, (1.0e+10,i=18,20), 1000000., (1.0e+10,i=22,31) /)
  REAL, PARAMETER, DIMENSION(0:Max_Kind) :: zero_val =                        &
  &  (/ 0., 0., 0., 0., 0., 0., 1., (0.0,i=7,16),                             &
  &    1.0, (0.0,i=18,20), 1.001e-4, (0.0,i=22,31) /)
  REAL, PARAMETER, DIMENSION(0:Max_Kind) :: zero_val2 =                       &
  &  (/ 0., 0., 0., 0., 0., 0., 1., (0.0,i=7,16),                             &
  &    1.0, (0.0,i=18,20), 0.0, (0.0,i=22,31) /)
  REAL, PARAMETER, DIMENSION(0:Max_Kind) :: fact_val =                        &
  &  (/ 1., 1., 1., 1., 1., 1., 1., (1.0,i=7,16),                             &
  &    -1.0, (1.0,i=18,20), 1.0e+4, (1.0,i=22,31) /)

  save NEWSTYLE, exptab, maxkind

  data NEWSTYLE /.false./

  data exptab /0.0001D0, 0.001D0, 0.01D0, 0.1D0, 1.0, 10.0, 100.0,            &
  &  1000.0, 10000.0, 100000.0, 1000000.0, 10000000.0,                        &
  &  100000000.0, 1000000000.0, 10000000000.0, 100000000000.0 /

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

  if (mode.gt.0) then  ! .... Conversion P,KIND a IP  ....

      if ( is_invalid_kind(kind) ) then
          write(6,6004) kind
!           call qqexit(1)    !  force excessive ?
          ip = -999999 
          return  ! erreur si kind pas valide
      endif
      if (kind .eq. 2 .and. p .eq. 0.) then  ! ou ajouter .and. .not. NEWENCODING
         ip = 0
!          if(NEWENCODING) ip =  ishft(2,24) +  ishft(1,20) ! si  newstyle (kind=2, mantissa=0, exponent=1)
         return
      endif

      if(NEWENCODING)then    ! new style excoding
          if(iand(kind,15) == 15) then  ! kind 15 and subkinds done elsewhere (pure integers)
            status = conv_kind_15(p,kind,ip,mode)
            return
          endif
          if (p .lt. low_val(kind) .or. p .gt. hi_val(kind)) then
              write(6,6006) p,low_val(kind),hi_val(kind)
              ip = -999999
              return
          endif
          iexp = 4
          temp = p
          if (abs(temp) .lt. zero_val(kind)) temp = zero_val(kind)
          temp = temp * fact_val(kind)  ! apply scaling factor before conversion to mantissa and pseudo exponent
          if ( temp .ge. 0) then
              limit1 = 1000000
              limit2 = 100000
              offset = 0
          else
              temp = -temp
              limit1 = 48000
              limit2 = 4800
              offset = 1000000
          endif
!          temp=temp*1.00000005D0
          do while ( iexp .gt. 0 .and. iexp .lt. 15 )  ! must keep pseudo exponent in range
              if (temp .ge. limit1 ) then        ! too big, divide by 10 and adjust pseudo exponent
                temp = temp / TEN
                iexp = iexp -1
              else if ( temp .lt. limit2 ) then  ! too small multiply by 10 and adjust pseudo exponent
                temp = temp * TEN
                iexp = iexp + 1
              else   ! >=limit2 and <limit1
                EXIT
              endif
          enddo
          if ( temp .gt. limit1 ) then          ! number is too big, cannot code
              ip = -1
          else
              ip = offset + nint(temp)
          endif
          ip = ior (ip,ishft(iexp,20))          ! add pseudo exponent
          ip = ior (ip,ishft(iand(15,kind),24)) ! add primary kind flag

      else ! OLD style encoding

          if (kind.eq.0) then   ! ...  hauteur ...
              ip = max( 12001,min( 32000,nint( p/5.0 + 12001 ) ) )
          elseif (kind.eq.1) then  ! ...  sigma ...
              if ( .not. (  0.0 .le. p .and. p .le. 1.0 ) ) then
                write(6,6001) p
                ip = -999999
                return
              endif
              ip = nint( p * 10000. ) + 2000
          elseif (kind.eq.2) then  ! ...  pression ...
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
          elseif (kind.eq.3) then  ! ...  code arbitraire
              ip = nint( p )
              if ( 0 .le. ip .and. ip .le. 100 ) then
                ip = 1200 - ip
              else
                write(6,6003) p
                ip = -999999
                return
              endif
          else  !  OLD encoding not valid for this kind
              write(6,6004) kind
              ip = -999999
              return
          endif
      endif  ! .not NEWENCODING, OLD style encoding

  elseif (mode.lt.0) then  ! ....  Conversion de ip a p,kind .....
      flag = flagv
      lstring=0
      if(flag) then
        lstring=len(string)
        if(lstring<9) flag=.false.   ! if less than 9 characters are available, no attempt to format value will be made
      endif
      if ( ip .gt. 32767 ) then  !   tous types, nouveau codage
          p = 0.0
          kind = iand(15,ishft(ip,-24))
          if(kind == 15) then  ! type 15 et associes traite a part
            if(conv_kind_15(p,kind,ip,mode) /= 0) goto 777  ! il y a une erreur de decodage pour ip
            if (flag) goto 666  ! impression dans string
            return   ! decodage termine, return
          endif
          if ( .not. validkind(kind) ) goto 777

          iexp = iand (15,ishft(ip,-20))
          itemp = iand (1048575, ip)
          if (itemp > 1000000) itemp = -(itemp - 1000000)
 555        continue
          p = itemp / exptab(iexp)             ! apply pseudo exponent
          p = p / fact_val(kind)               ! apply scaling factor
 
          if (p < low_val(kind) .or. p>hi_val(kind)) then ! hors limite, essayer le type associe si valide
            if(kind+16 <= Max_Kind) then
              if(validkind(kind) .and. validkind(kind+16)) then
                kind = kind+16
                goto 555         ! try new kind
              else
                goto 777  ! invalid kind
              endif
            else
              goto 777  ! invalid kind
            endif
          endif
          p = max(p,low_val(kind))     ! clipping a la valeur minimale
          p = min(p,hi_val(kind))      ! clipping a la valeur maximale
          if (abs(p) .lt. 1.001*zero_val(kind)) p = zero_val2(kind)   ! mise a "zero" si valeur absolue trop faible
666       abs_p = abs(p)
          if (flag) then  ! convert P into a formatted string with appropriate units for kind
             string2=""
             status=value_to_string(p , string2 , min(len(string2),len(string)-3) )
             string=trim(string2)//' '//kind_to_string(kind)
          endif
      elseif (  12000 .lt. ip .and. ip .le. 32000) then  !  ...  hauteur old style ...
          kind = 0
          p = 5 * ( ip -12001)
          if (flag) write(string,'(i6,1x,a1)') nint(p),'m'
      elseif (  2000 .le. ip .and. ip .le. 12000 ) then  !  ...  sigma old style ...
          kind = 1
          p = float (ip - 2000) / 10000.
          if (flag) write(string,'(f6.4,1x,a2)') p,'sg'
      elseif (( 0    .le. ip .and. ip .lt. 1100 )  .or. ( 1200 .lt. ip .and. ip .lt. 2000 )) then  !  ... pression old style ...
          kind = 2
          if ( 0 .le. ip .and. ip .lt. 1100 ) then
             p = float(ip)
             if (flag) write(string,'(i6,1x,a2)') ip,'mb'
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
      elseif ( 1100 .le. ip .and. ip .le. 1200) then  ! ...  code arbitraire old style ...
          kind = 3
          p = float( ip )
          p = 1200. - p
          if (flag) write(string,'(i6,3x)') nint(p)
      else  !  Valeur inderminee de ip  old style
          kind = 3
          p = float( ip )
      endif   !    ip .gt. 32767 elseif, elseif, 
  endif  ! ....  Conversion de xx a yy .....
      
  return

777  continue  ! invalid ip, return kind = -1
  p = -999999
  kind = -1
  if(flagv) string = 'Invalid'
  return

  6001 format(' Error in convip: sigma value =',e12.5,' returned ip is -999999')
  6002 format(' Error in convip: pressure value =',e12.5,' returned ip is -999999')
  6003 format(' Error in convip: arbitrary value=',e12.5,' returned ip is -999999')
  6004 format(' Error in convip: invalid kind =',I10)
  6005 format(' Error in convip: kind=10 (oldstyle) value out of range=',e12.5,' returned ip is -999999')
  6006 format(' Error in convip: p is out of bounds =',e12.5,' min=',e12.5,' max=',e12.5,' returned ip is -999999')
! 6007 format(' Warning in convip: undetermined kind used =',I10)

  contains
!=========================== start of private function =========================================
function conv_kind_15(p,mykind,ip,mode) result(status) ! convert kind = 15 and subkinds
  implicit none
  integer :: status
  integer, intent(INOUT) :: mykind,ip
  integer, intent(IN) :: mode ! -1, 1, 2, 3 (see convip code for meaning of mode)
  real, intent(INOUT) :: p

  type e15
    integer :: lo      ! lowest value of ipv for this sub kind
    integer :: hi      ! lhighest value of ipv for this sub kind
    integer :: base    ! offset for this sub kind
  end type
  type(e15), dimension(2), save :: t15 = & 
         (/ &
         e15(       0,  2000000,      0), &     ! values between 0 and 1 999 999    (kind 15)
         e15(16000000, 15000001,  -1000)  &     ! values between -1000 and 998 999 (kind 31) ! entries swapped to deactivate
         /)
  integer :: i, subt, ipv
  integer, parameter :: FFFFFF=Z'FFFFFF'

  status = -1               ! precondition for failure

  if(ip > 0 .and. ishft(ip,-24) == 15 .and. mode == -1) then  ! kind 15 and sub kinds ip to p conversion
    mykind = -1
    ipv = iand(ip,FFFFFF)  ! get rid of kind 15 indicator
    subt = -1
    do i=1,size(t15)   ! lookup in bounds table
      if(ipv >= t15(i)%lo .and. ipv <= t15(i)%hi ) then ! is ipv in the range of this sub kind ?
        subt = i    ! yes
        exit
      endif
    enddo
    if(subt == -1) return  ! invalid ip value for kind = 15 and associated sub kinds
    p = ipv - t15(subt)%lo + t15(subt)%base   ! convert ipv to actual value
    mykind = 15 + 16*(subt-1)    ! return proper kind type
    status = 0
  endif

  if(15 == iand(mykind,15) .and. mode > 0 .and. mode <3) then  ! newstyle p to ip conversion
    ip = -1                      ! precondition for fail
    subt = 1 + ishft(mykind,-4)
    if(subt <= 0 .or. subt > size(t15)) return   ! sub kind out of range
    ipv = nint(p) - t15(subt)%base + t15(subt)%lo
    if(ipv < t15(subt)%lo .or. ipv > t15(subt)%hi) return  ! p is out of range
    ip = ior(ipv,ishft(15,24))  ! add type 15 flag
    status = 0
  endif
! total failure if we get here
  return
end function conv_kind_15

end SUBROUTINE CONVIP_plus
!===============================================================================================
!============================= end of private function =========================================
!===============================================================================================
!****f* rmnlib/value_to_string
! SUMMARY
! write value val into string using at most maxlen characters
! SYNOPSIS
integer function value_to_string(val,string,maxlen)  
! AUTHOR
!  M.Valin 2013
!  Revision 001 :  M.Valin  Oct 2013 alignement a droite corrige pour valeurs entieres > 0
! ARGUMENTS
  integer, intent(IN) :: maxlen
  character (len=*), intent(OUT) :: string
  real *4, intent(IN) :: val
! INPUTS
!  maxlen : use at most maxlen characters to code value
!  val    : real value to be encoded (left aligned) into string
! OUTPUTS
!  string : result of encoding
! RESULT
!  if something went wrong, the return value will be <= 0
!  if an integer format (I) is used, the return value is the number of characters used
!  if a floating format (F/G) is used, return value =- 100*field_width + nb_of_significant_digits
!  ex: if F8.3 is used, the result will be 803, if I6 is used, the result will be 6
!******
  character (len=32) :: fstring
  character (len=128) :: tempstring
  real *4 :: value
  integer :: after, before
  integer :: grosint, maxc, intdig

  string=" "
  maxc=min(maxlen,len(string))
  value=abs(val)
  after=min(7,maxc-6)
  before=maxc
  value_to_string=-(100*before+after*10)
  write(fstring,11)maxc,after    ! default G format

  if(value /= 0.0) then
    if(value >= 1000000000000.0 .or. value < .0001) goto 666   ! use G format
  endif

  if(nint(value)==value) then ! exact integral value
    grosint=1
    intdig=2
    do i=1,min(9,maxc-1)
      if(nint(value) > grosint) intdig=intdig+1
      grosint=grosint*10  ! largest integer value that will fit in maxc characters
    enddo
    if(value >= grosint) goto 444   ! try something else
    if(val>0) intdig = intdig - 1   ! one less character if number is positive
    write(fstring,12)'(I',min(maxc,intdig),')'    ! use I format
    value_to_string=min(maxc,intdig)
    goto 777
  endif

  444 continue  ! real values within "civilized" range
  if (value >= 1.0) then
    before = 0
    after = 0
    do while(value>=1.0)
      before = before +1
      value = value * .1
    enddo
    if(before<6) after=min(6-before,maxc-before-2) ! we have at best 6 significant digits
!    if(before<8) after=max(0,maxc-before-2)
  else   ! value < 1.0
    after = 5
    before = 0
    do while(value<1.0)
      value = value * 10.0
      after = after  +1
    enddo
    after=min(after,maxc-2)
  endif

  after=min(9,after)  ! never more than 9 digits after the decimal point

  if(before+after+2 > maxc) goto 666  ! use G format

  before=before+after+2
  value_to_string=100*before+after*10

  write(fstring,10)before,after       ! use F format
  !print *,'=',trim(fstring)
  write(string,fstring)val            ! F format
  i=len(trim(string))
  do while(i>4 .and. string(i:i)=='0')
    string(i:i)=' '
    i=i-1
  enddo
  if(string(1:1)==' ') then
    tempstring=string(2:len(string))
    string=tempstring
  endif
  return

  666 continue
  if(maxc-6<=0) goto 888
  !print *,'=',trim(fstring)
  write(string,fstring)val            ! G format
  if(string(1:1)==' ') then
    tempstring=string(2:len(string))
    string=tempstring
  endif
  return

  777 continue
  !print *,'=',trim(fstring)
  write(string,fstring)nint(val)      ! I format
  if(string(1:1)==' ') then
    tempstring=string(2:len(string))
    string=tempstring
  endif
  return

  888 continue
  value_to_string=0
  return

10 format(2H(F,I2,1H.,I1,1H))
11 format(2H(G,I2,1H.,I1,1H))
12 format(A,I2,A)
end function value_to_string
!===============================================================================================
subroutine test_value_to_string
  implicit none
  include 'convip_plus.inc'
  character (len=8) :: stringa
  character (len=12) :: stringb
  character (len=15) :: stringc
  integer :: i
  integer :: status
  real *4 :: value

  value=1.000001
  do i=1,9
    status=value_to_string(real(nint(value)),stringa,15)
    stringc = stringa
    print 101,stringc,trim(stringa),'',status*.01
    value=value*10.0
  enddo

  value=1.000001
  do i=1,9
    status=value_to_string(real(nint(-value)),stringa,15)
    stringc = stringa
    print 101,stringc,trim(stringa),'mb',status*.01
    value=value*10.0
  enddo

  value=1.234567
  do i=1,12
    status=value_to_string(-value,stringb,15)
    stringc = stringb
    print 101,stringc,trim(stringb),'mb',status*.01
    value=value*10.0
  enddo

  value=1.23456789
  do i=1,12
    status=value_to_string(-value,stringc,12)
    print 101,stringc,trim(stringc),'mb',status*.01
    value=value*0.1
  enddo

101 format(1H|,A15,1H|,A15,1H|,1X,A2,3X,f6.2)
return
end subroutine test_value_to_string
!===============================================================================================
subroutine test_convip_plus() ! test routine for convip_plus
  use ISO_C_BINDING
  implicit none
  include 'convert_ip123.inc'
  integer :: ip1, ip2, i, j, nip, nip2, k2, nip3
  real :: p,p2
  character(len=15) :: string
  integer :: kind
!     test #1, ip1 -> p,kind -> ip2 (ip2 must be = ip1)
  nip = 0
  nip2 = 0
  nip3 = 0
  do j=0,15
  do i=0,1047999
    ip1 = i
    ip1=ior(ip1, ishft(j,20))  ! add exponent
    ip1=ior(ip1, ishft(3,24))  ! kind =3 to test full range
    call CONVIP_plus( ip1, p, kind, -1, string, .false. )  ! ip1 -> p,kind 
    call CONVIP_plus( ip2, p, kind, +2, string, .false. )  ! p,kind -> ip2
    call CONVIP_plus( ip2, p2, k2, -1, string, .false. )
    if(ip1/=ip2) nip=nip+1
    if(p/=p2 .and. abs(p2/p-1.0) < .0000002) nip2=nip2+1
    if(ip1/=ip2 .and. p /= p2) then
      if(abs(p2/p-1.0) >= .0000002) then ! not within tolerance
        print 111, j,i,ip1,ip2,iand(ip1,1048575),iand(ip2,1048575),p,p2,abs(p2/p-1.0)
        nip3 = nip3+1
!           stop
      endif
    endif
  enddo
  enddo
  print 112,'ip1<>ip2 (normalization aliases)=',nip,' p1 ~= p2 (within 2.0E-7 relative error)',nip2,' errors=',nip3
  do i=5,1005,100
     ip1 = i
     call CONVIP_plus( ip1, p, kind, -1, string, .false.)
     p=p+.751
     call CONVIP_plus( ip1, p, kind, +2, string, .false. )
     call CONVIP_plus( ip1, p, kind, -1, string, .true.)
     print 113,i,p,kind,':'//trim(string)//':'
  enddo
111 format(2I9,2Z8,2I9,3G16.8)
112 format(A,I9,A,I9,A,I9)
113 format(I9,F11.5,I3,A)
  return
end subroutine test_convip_plus
!===============================================================================================
subroutine c_kind_to_string(code,s1,s2) BIND(C,name='KindToString')  ! interface for C routines
! translate kind integer code to 2 character string, gateway to Fortran kind_to_string
  use ISO_C_BINDING
  implicit none
  include 'convip_plus.inc'

  integer(C_INT), intent(IN), value :: code
  character(C_CHAR), intent(OUT) :: s1, s2

  character(len=2) :: temp
  temp = kind_to_string(code)
  s1 = transfer(temp(1:1),s1)
  s2 = transfer(temp(2:2),s2)
  return
end subroutine c_kind_to_string

FUNCTION kind_to_string(code) RESULT(string)  ! translate ip kind into a 2 character string code
  implicit none
  integer, intent(IN) :: code
  character(len=2) :: string
  integer, parameter :: Max_Kind=31
  character(len=2), save, dimension(0:Max_Kind) :: kinds = &
    (/    ' m', 'sg', 'mb', '  ', ' M', 'hy', 'th', '??',                       &
          '??', '??', ' H', '??', '??', '??', '??', '_0',                       &
          '??', '[]', '??', '??', '??', 'mp', '??', '??',                       &
          '??', '??', '??', '??', '??', '??', '??', '_1' /)


  string = '!!'    ! precondition to fail

  if(code<0) return   ! invalid type code

  if(code<=Max_Kind) then
    string = kinds(code)  ! straight table lookup
    return
  endif

  if(mod(code,16)/=15) return  ! not a subkind of kind 15

  if(code/16>9)  return  ! not a potentially valid subkind of kind 15

  write(1,string)'I',code/16   ! 'In' code where n=code/16 (n=0,1..,9)
1 format(A,I1)

  return
end FUNCTION kind_to_string
!===============================================================================================

!> C language interface with no string option
subroutine C_CONV_IP( ip, p, kind, mode ) BIND(C,name='ConvertIp')
  use ISO_C_BINDING
  implicit none
    integer(C_INT), intent(INOUT) :: ip, kind
    integer(C_INT), intent(IN), value :: mode
    real(C_FLOAT), intent(INOUT) :: p
    character (len=1) :: string
    integer :: mode2
    mode2 = mode
    call CONVIP_plus( ip, p, kind, mode2,string,.false.)
end subroutine C_CONV_IP
