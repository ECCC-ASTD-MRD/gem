

subroutine yydecode()
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
   include "rmnlib_basics.inc"

   !@author   M.Desgagne   -   Spring 2012

!!$   integer, external :: fnom,fstouv,fstinf,fstinl,fstprm,fstluk,exdb,exfin, &
!!$        fstecr,fstfrm,fclos,fstopl,fstsel,fstlis,fstnbr

   character(LEN=1024) :: LISTEc(3), DEF(3), VAL(3)
   integer NPOS
   data LISTEc /'yin.'   , 'yan.' , 'i.'    /
   data VAL    /'/null'  , '/null', '/null' /
   data DEF    /'/null'  , '/null', '/null' /

   character(len=2)    :: typ_S, grd_S
   character(len=4 )   :: var_S
   character(len=12)   :: lab_S
   character(len=1024) :: yin_S,yan_S,in_S

   integer dte, det, ipas, p1, p2, p3, g1, g2, g3, g4, bit, &
        dty, swa, lng, dlf, ubc, ex1, ex2, ex3, nhgrids
   integer, dimension (:,:,:), allocatable :: igs

   integer i,j,nhours,datev,key,iun1,iun2,iun3,sindx
   integer nlis,lislon, ni1,nj1,nk1,ni,nj,err,lindex
   integer, dimension(:), allocatable :: liste,niv
   logical secondMETA_L

   real(REAL64), dimension(:),allocatable :: posz
   real  , dimension(:),allocatable :: yy,champ
   real dummy
   real(REAL64) :: nhours_8
   !
   !-------------------------------------------------------------------
   !
   err = exdb('YYDECODE','3.0', 'NON')

   NPOS = 1
   call CCARD(LISTEc,DEF,VAL,3,NPOS)
   yin_S = val(1)
   yan_S = val(2)
   in_S  = val(3)

   iun1 = 0 ; iun2 = 0 ; iun3 = 0
   secondmeta_L=.false.

   if (fnom(iun1,in_S,'RND+OLD',0) >= 0) then
      if (fstouv(iun1,'RND') < 0) then
         write(6,8001) in_S
         stop
      endif
   else
      write (6,8000) in_S
      stop
   endif

   if (fnom(iun2,yin_S,'RND',0) >= 0) then
      if (fstouv(iun2,'RND') < 0) then
         write (6,8001) yin_S
         stop
      endif
   else
      write (6,8000) yin_S
      stop
   endif
   if (fnom(iun3,yan_S,'RND',0) >= 0) then
      if (fstouv(iun3,'RND') < 0) then
         write (6,8001) yan_S
         stop
      endif
   else
      write (6,8000) yan_S
      stop
   endif

   nlis= fstnbr(iun1)
   !print *,'nlis=',nlis
   allocate(liste(nlis),niv(nlis))

   err= fstinl(iun1,ni1,nj1,nk1,-1,' ',-1,-1,-1,' ','^>',&
        liste,nhgrids,nlis)
   !print *,'err of fstinl=',err,'nhgrids=',nhgrids

   allocate(igs(3,3,max(nhgrids,1)))
   do i=1,nhgrids
      err= fstprm(liste(i), dte, det, ipas, ni, nj, nk1, bit , &
           dty, igs(1,1,i), igs(2,1,i), igs(3,1,i), typ_S , &
           var_S, lab_S, grd_S, g1, &
           g2, g3, g4, swa, lng, dlf, ubc, ex1, ex2, ex3)
      allocate(yy(ni))
      err = fstluk(yy, liste(i), ni1,nj1,nk1)
      !print *,'err of fstluk=',err,'ni for yy=',ni,' ni1,nj1,nk1=',ni1,nj1,nk1
      sindx  = 6

      ni = nint(yy(sindx  ))
      nj = nint(yy(sindx+1))
      !print *,'ni is redone to',ni,'nj=',nj
      call cxgaig('E', g1,g2,g3,g4, &
           yy(sindx+6), yy(sindx+7), yy(sindx+8), yy(sindx+9))
      !print *,'sindx+10+ni=',sindx+10+ni
      call set_igs2(igs(1,2,i), igs(2,2,i)            , &
           yy(sindx+10),yy(sindx+10+ni),ni,nj, &
           g1,g2,g3,g4, 1,ni,1,1,nj,1)
      igs(3,2,i)= igs(3,1,i)
      err = FSTECR(yy(sindx+10), yy, -32, iun2, 0, 0, 0, ni, 1, 1  , &
           igs(1,2,i),igs(2,2,i),igs(3,2,i), 'X', '>>', 'YYG_POSX', 'E', &
           g1, g2, g3, g4, dty, .true.)
      err = FSTECR ( yy(sindx+10+ni), yy, -32, iun2, 0, 0, 0, 1, nj, 1, &
           igs(1,2,i),igs(2,2,i),igs(3,2,i), 'X', '^^', 'YYG_POSY', 'E', &
           g1, g2, g3, g4, dty, .true.)
      sindx = sindx+10+ni+nj
      call cxgaig('E', g1,g2,g3,g4, &
           yy(sindx+6), yy(sindx+7), yy(sindx+8), yy(sindx+9))
      !print *,'sindx+10+ni=',sindx+10+ni, 'for yy'
      call set_igs2(igs(1,3,i), igs(2,3,i)            , &
           yy(sindx+10),yy(sindx+10+ni),ni,nj, &
           g1,g2,g3,g4, 1,ni,1,1,nj,1)
      igs(3,3,i)= igs(3,1,i)
      err = FSTECR(yy(sindx+10), yy, -32, iun3, 0, 0, 0, ni, 1, 1  , &
           igs(1,3,i),igs(2,3,i),igs(3,3,i), 'X', '>>', 'YYG_POSX', 'E', &
           g1, g2, g3, g4, dty, .true.)
      err = FSTECR(yy(sindx+10+ni), yy, -32, iun3, 0, 0, 0, 1, nj, 1, &
           igs(1,3,i),igs(2,3,i),igs(3,3,i), 'X', '^^', 'YYG_POSY', 'E', &
           g1, g2, g3, g4, dty, .true.)
      deallocate (yy)
   end do

   err= fstinl(iun1,ni1,nj1,nk1,-1,' ',-1,-1,-1,' ',' ',&
        liste,lislon,nlis)
   !print *,'viv err=',err
   do i=1,lislon
      err= fstprm(liste(i), dte, det, ipas, ni, nj, nk1, bit , &
           dty, p1, p2, p3, typ_S, var_S, lab_S, grd_S, g1, &
           g2, g3, g4, swa, lng, dlf, ubc, ex1, ex2, ex3)

      nhours = det * ipas / 3600.d0
      if (dte  >  0) then
         nhours_8 = nhours
         call incdatr(datev, dte, nhours_8)
         key= FSTINF(iun1, NI1, NJ1, NK1, datev, lab_S, p1, p2, p3, typ_S, var_S)
      else
         key= FSTINF(iun1, NI1, NJ1, NK1,    -1, lab_S, -1, -1, -1, typ_S, var_S)
      endif

      if (trim(var_S) == '^>') cycle
      if ((trim(var_S) == '^^') .or. (trim(var_S) == '>>') &
           .or. (trim(var_S) == '!!') .or. (trim(var_S) == 'META' )) then
         allocate(posz(ni*nj))
         err = fstluk(posz, liste(i),ni1,nj1,nk1)
         if (trim(var_S) == 'META') then
             if (secondmeta_L) then
                 err = FSTECR(posz, dummy, -bit, iun3, dte, det, ipas, ni1, nj1,&
                 nk1, p1, p2, p3, typ_S, var_S, lab_S, grd_S,&
                 g1, g2, g3, g4, dty, .false.)
             else
                 err = FSTECR(posz, dummy, -bit, iun2, dte, det, ipas, ni1, nj1,&
                 nk1, p1, p2, p3, typ_S, var_S, lab_S, grd_S,&
                 g1, g2, g3, g4, dty, .false.)
                 secondmeta_L=.true.
             endif
         else
              err = FSTECR(posz, dummy, -bit, iun2, dte, det, ipas, ni1, nj1,&
              nk1, p1, p2, p3, typ_S, var_S, lab_S, grd_S,&
              g1, g2, g3, g4, dty, .true.)
         endif
         if (trim(var_S) == '!!' ) &
              err = FSTECR(posz, dummy, -bit, iun3, dte, det, ipas, ni1, nj1,&
              nk1, p1, p2, p3, typ_S, var_S, lab_S, grd_S,&
              g1, g2, g3, g4, dty, .true.)
         deallocate(posz, stat=err)

      else

         lindex=0
         do j=1,nhgrids
            if ((igs(1,1,j).eq.g1) .and. &
                (igs(2,1,j).eq.g2) .and. &
                (igs(3,1,j).eq.g3)) lindex= j
         end do

         allocate(champ(ni*nj))
         err = fstluk(champ,liste(i),ni1,nj1,nk1)
         if (lindex  >  0 .and. grd_S(1:1) == 'U') then
            err = FSTECR(champ, champ, -bit, iun2, dte, det, ipas, ni1, nj1/2, nk1,&
                 p1, p2, p3, typ_S, var_S, lab_S, 'Z', &
                 igs(1,2,lindex),igs(2,2,lindex),igs(3,2,lindex),g4,dty,.false. )
            err = FSTECR(champ(ni1*nj1/2+1), champ, -bit, iun3, dte, det, ipas, &
                 ni1, nj1/2, nk1, p1, p2, p3, typ_S, var_S, lab_S, 'Z', &
                 igs(1,3,lindex),igs(2,3,lindex),igs(3,3,lindex),g4,dty,.false. )
         else
            err = FSTECR(champ, champ, -bit, iun2, dte, det, ipas, ni1, nj1, nk1,&
                 p1, p2, p3, typ_S, var_S, lab_S, grd_S(1:1), &
                 g1,g2,g3,g4,dty,.false. )
         endif
         deallocate(champ, stat=err)

      endif
   end do

99 err = fstfrm(iun1)
   err = fstfrm(iun2)
   err = fstfrm(iun3)

   err = exfin('YYDECODE','3.0', 'OK')

8000 format (/' Unable to fnom: '  ,a/)
8001 format (/' Unable to fstouv: ',a/)

   !-------------------------------------------------------------------
   return
end subroutine yydecode
