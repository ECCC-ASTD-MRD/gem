
subroutine yyencode()
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
   include "rmnlib_basics.inc"
   !@author   M.Desgagne   -   Spring 2012
   !@revision V.Lee Spring 2014 (rewrit in fstecr is not date sensitive)
!!$   integer, external :: fnom,fstouv,fstinf,fstinl,fstprm,fstluk,exdb,exfin, &
!!$        fstecr,fstfrm,fclos,fstopl,fstsel,fstlis,fstnbr

   character(len=1024) :: LISTEc(3), DEF(3), VAL(3)
   integer NPOS
   data LISTEc /'yin.'   , 'yan.' , 'o.'    /
   data VAL    /'/null'  , '/null', '/null' /
   data DEF    /'/null'  , '/null', '/null' /

   character(len=1)    :: family_uencode_S
   character(len=2)    :: typ_S, grd_S
   character(len=4)    :: var_S
   character(len=12)   :: lab_S
   character(len=1024) :: yin_S,yan_S,out_S

   integer :: dte, det, ipas, p1, p2, p3, g1, g2, g3, g4, bit, &
        dty, swa, lng, dlf, ubc, ex1, ex2, ex3,ip1,ip2,ip3
   integer :: p2a, p3a, g2a, g3a

   integer :: iun1,iun2,iun3,maxni,maxnj,i,datev,niyy,vesion_uencode
   integer :: nlis,lislon, key, ni1,nj1,nk1,ni,nj,err,sindx_yin,sindx
   integer :: yinlislon,yanlislon,taclislon,j
   integer, dimension(:), allocatable :: liste,niv,yinliste,yanliste,tacliste

   real  :: xlat1,xlon1, xlat2,xlon2
   real, dimension(:),allocatable :: champ, yy
   real(REAL64) :: nhours

   !-------------------------------------------------------------------

   err = exdb('YYENCODE','3.0', 'NON')
   NPOS = 1
   call CCARD(LISTEc,DEF,VAL,3,NPOS)
   yin_S = val(1)
   yan_S = val(2)
   out_S = val(3)

   iun1 = 0 ; iun2 = 0 ; iun3 = 0

   if (fnom(iun1,yin_S,'RND+OLD',0) >= 0) then
      if (fstouv(iun1,'RND') < 0) then
         write (6,8001) yin_S
         stop
      endif
   else
      write (6,8000) yin_S
      stop
   endif
   if (fnom(iun2,yan_S,'RND+OLD',0) >= 0) then
      if (fstouv(iun2,'RND') < 0) then
         write (6,8001) yan_S
         stop
      endif
   else
      write (6,8000) yan_S
      stop
   endif
   if (fnom(iun3,out_S,'RND',0) >= 0) then
      if (fstouv(iun3,'RND') < 0) then
         write (6,8001) out_S
         stop
      endif
   else
      write (6,8000) out_S
      stop
   endif

   nlis = fstnbr(iun1)
   allocate(liste(nlis),niv(nlis),yinliste(nlis),yanliste(nlis),tacliste(nlis))

   err= fstinl(iun1,ni1,nj1,nk1,-1,' ',-1,-1,-1,' ',' ',&
        liste,lislon,nlis)
   maxni=0 ; maxnj=0
   do i=1,lislon
      err= fstprm (liste(i), dte, det, ipas, ni, nj, nk1, bit , &
           dty, p1, p2, p3, typ_S, var_S, lab_S, grd_S, g1, &
           g2, g3, g4, swa, lng, dlf, ubc, ex1, ex2, ex3)
      maxni= max(maxni,ni)
      maxnj= max(maxni,nj)
   end do


   if (lislon < 1) then
      write(6,'(/3x,"NOTHING to DO -- QUIT"/)')
      stop
   endif

   err= fstinl(iun1,ni,nj1,nk1,-1,' ',-1,-1,-1,' ','>>',&
        yinliste,yinlislon,nlis)

   if (yinlislon == 0.or.yinliste(1) < 0) then
      write(6,'(/3x,"YIN positionnal parameters >> not available - ABORT"/)')
      stop
   endif

   err= fstinl(iun1,ni1,nj,nk1,-1,' ',-1,-1,-1,' ','^^',&
        tacliste,taclislon,nlis)
   if (taclislon == 0 .or. tacliste(1) < 0) then
      write(6,'(/3x,"YIN positionnal parameters ^^ not available - ABORT"/)')
      stop
   endif

   if (taclislon.ne.yinlislon) then
      write(6,'(/3x,"YIN positionnal parameters not all available - ABORT"/)')
      stop
   endif

   err= fstinl(iun2,ni,nj1,nk1,-1,' ',-1,-1,-1,' ','>>',&
        yanliste,yanlislon,nlis)

   if (yanlislon == 0 .or. yanliste(1) < 0) then
      write(6,'(/3x,"YAN positionnal parameters >> not available - ABORT"/)')
      stop
   endif

   err= fstinl(iun2,ni1,nj,nk1,-1,' ',-1,-1,-1,' ','^^',&
        tacliste,taclislon,nlis)
   if (taclislon == 0 .or. tacliste(1) < 0) then
      write(6,'(/3x,"YAN positionnal parameters ^^ not available - ABORT"/)')
      stop
   endif

   if (taclislon.ne.yanlislon) then
      write(6,'(/3x,"YAN positionnal parameters not all available - ABORT"/)')
      stop
   endif

   allocate (champ(maxni*2*maxni))

   DO_FLDS: do i=1,lislon
      err= fstprm(liste(i), dte, det, ipas, ni, nj, nk1, bit , &
           dty, p1, p2, p3, typ_S, var_S, lab_S, grd_S, g1, &
           g2, g3, g4, swa, lng, dlf, ubc, ex1, ex2, ex3)

      nhours = det * ipas / 3600.d0
      datev  = -1
      if (dte > 0) call incdatr(datev, dte, nhours)

      key= FSTINF(iun2, NI1, NJ1, NK1, datev, ' ', p1, p2, p3, typ_S, var_S)

      if (var_S == '!!' .or. var_S == '>>' .or. var_S == '^^' .or. var_S == 'META') then
         if (var_S == '!!') then
            err = fstluk(champ,liste(i),ni1,nj1,nk1)
            err = FSTECR(champ, champ, -bit, iun3, dte, det, ipas, ni1, nj1, &
                 nk1, p1, p2, p3, typ_S, var_S, lab_S, grd_S, &
                 g1, g2, g3, g4, dty, .false.)
         endif
         if (var_S == 'META') then
            err = fstluk(champ,liste(i),ni1,nj1,nk1)
            err = FSTECR(champ, champ, -bit, iun3, dte, det, ipas, ni1, nj1, &
                 nk1, p1, p2, p3, typ_S, var_S, lab_S, grd_S, &
                 g1, g2, g3, g4, dty, .false.)
            key = FSTINF(iun2, NI1, NJ1, NK1, datev, ' ', -1, -1, -1, typ_S, var_S)
            err = fstluk(champ,key,ni1,nj1,nk1)
            err = FSTECR(champ, champ, -bit, iun3, dte, det, ipas, ni1, nj1, &
                 nk1, p1, p2, p3, typ_S, var_S, lab_S, grd_S, &
                 g1, g2, g3, g4, dty, .false.)
         endif
      else
         if (key < 0) then
            write(6,'(/3x,"Corresponding YAN variable: ",a," NOT FOUND - ABORT"/)') var_S
            stop
         endif
         g3a  =  1                            ! points de masse
         g2a = g2
         if (trim(var_S) == 'UT1' .or. trim(var_S) == 'URT1')  then
            g3a = 2 ! points U
            g2a= g2-1
         endif
         if (trim(var_S) == 'VT1' .or. trim(var_S) == 'VRT1')  then
            g3a = 3 ! points V
            g2a= g2-2
         endif
         err = fstluk(champ,liste(i),ni1,nj1,nk1)
         err = fstluk(champ(ni1*nj1+1),key,ni1,nj1,nk1)
         err = FSTECR(champ, champ, -bit, iun3, dte, det, ipas, ni1, 2*nj1, &
              nk1, p1, p2, p3, typ_S, var_S, lab_S, 'U', &
              g1, g2a, g3a, 0, dty, .false.)
      endif

   enddo DO_FLDS

   !Read and write out new grid descriptors'

   DO_TICTACS: do j=1,yinlislon

      err= fstprm(yinliste(j), dte, det, ipas, ni1, nj1, nk1, bit, &
           dty, ip1, ip2, ip3, typ_S, var_S, lab_S, grd_S, g1, &
           g2, g3, g4, swa, lng, dlf, ubc, ex1, ex2, ex3)

      call cigaxg('E', xlat1,xlon1, xlat2,xlon2, g1,g2,g3,g4)
      do i=1,lislon
         err= fstprm(liste(i), dte, det, ipas, ni, nj, nk1, bit , &
              dty, p1, p2, p3, typ_S, var_S, lab_S, grd_S, g1, &
              g2, g3, g4, swa, lng, dlf, ubc, ex1, ex2, ex3)
         if(trim(var_S).ne.'!!'.and.trim(var_S).ne.'META'.and.&
              trim(var_S).ne.'>>'.and.trim(var_S).ne.'^^') then
            if (g1 == ip1 .and. g2 == ip2 .and. g3 == ip3) then
               if (trim(var_S) == 'UT1' .or. trim(var_S)  ==  'URT1') then
                  p2a = g2-1
                  p3a = 2
               else if (trim(var_S) == 'VT1' .or. trim(var_S) == 'VRT1') then
                  p2a = g2-2
                  p3a = 3
               else
                  p2a = g2
                  p3a = 1
               endif
               exit
            endif
         endif
      enddo

      niyy=5+2*(10+ni+nj)
      allocate(yy(niyy))

      vesion_uencode    = 1
      family_uencode_S = 'F'

      yy(1) = iachar(family_uencode_S)
      yy(2) = vesion_uencode
      yy(3) = 2 ! 2 grids (Yin & Yang)
      yy(4) = 1 ! the 2 grids have same resolution
      yy(5) = 1 ! the 2 grids have same area extension on the sphere

      !YIN
      sindx  = 6
      yy(sindx  ) = ni
      yy(sindx+1) = nj
      yy(sindx+6) = xlat1
      yy(sindx+7) = xlon1
      yy(sindx+8) = xlat2
      yy(sindx+9) = xlon2
      err = fstluk(yy(sindx+10   ),yinliste(j),ni1,nj1,nk1)
      err = fstluk(yy(sindx+10+ni),tacliste(j),ni1,nj1,nk1)
      yy(sindx+2) = yy(sindx+10      )
      yy(sindx+3) = yy(sindx+ 9+ni   )
      yy(sindx+4) = yy(sindx+10+ni   )
      yy(sindx+5) = yy(sindx+ 9+ni+nj)
      sindx_yin= sindx

      !YAN

      err= fstprm(yanliste(j), dte, det, ipas, ni1, nj1, nk1, bit , &
           dty, p1, p2, p3, typ_S, var_S, lab_S, grd_S, g1   , &
           g2, g3, g4, swa, lng, dlf, ubc, ex1, ex2, ex3)

      call cigaxg('E', xlat1,xlon1, xlat2,xlon2, g1,g2,g3,g4)

      sindx   = sindx+10+ni+nj
      yy(sindx  ) = ni
      yy(sindx+1) = nj
      yy(sindx+2) = yy(sindx_yin+10      )
      yy(sindx+3) = yy(sindx_yin+ 9+ni   )
      yy(sindx+4) = yy(sindx_yin+10+ni   )
      yy(sindx+5) = yy(sindx_yin+ 9+ni+nj)
      yy(sindx+6) = xlat1
      yy(sindx+7) = xlon1
      yy(sindx+8) = xlat2
      yy(sindx+9) = xlon2
      yy(sindx+10    :sindx+9+ni  )= yy(sindx_yin+10   :sindx_yin+9+ni   )
      yy(sindx+10+ni:sindx+9+ni+nj)= yy(sindx_yin+10+ni:sindx_yin+9+ni+nj)

      err = FSTECR(yy, yy, -32, iun3, 0, 0, 0, niyy, 1, 1  , &
           ip1, p2a,  p3a, 'X', '^>', 'YYG_UE_GEMV4', &
           family_uencode_S, vesion_uencode,0,0,0, 5, .false.)
      deallocate(yy, stat=err)
   enddo DO_TICTACS

   err = fstfrm(iun3)

   err = exfin('YYENCODE','3.0', 'OK')

8000 format (/' Unable to fnom: '  ,a/)
8001 format (/' Unable to fstouv: ',a/)

   return
end subroutine yyencode
