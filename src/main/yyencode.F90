
#include "modelutils_build_info.h"

program yyencode
   use, intrinsic :: iso_fortran_env, only: REAL64
   use app
   use rmn_fst24

   implicit none
   include "rmnlib_basics.inc"
   !@author   M.Desgagne   -   Spring 2012
   !@revision V.Lee Spring 2014 (rewrit in fstecr is not date sensitive)

   character(len=1024) :: LISTEc(3), DEF(3), VAL(3)
   integer NPOS
   data LISTEc /'yin.'   , 'yan.' , 'o.'    /
   data VAL    /'/null'  , '/null', '/null' /
   data DEF    /'/null'  , '/null', '/null' /

   character(len=1)    :: family_uencode_S
   character(len=4)    :: nvar
   character(len=1024) :: yin_S,yan_S,out_S

   integer :: dte, det, ipas, p1, p2, p3, g1, g2, g3, g4, bit, &
        dty, swa, lng, dlf, ubc, ex1, ex2, ex3,ip1,ip2,ip3
   integer :: p2a, p3a, g2a, g3a

   integer :: iun1,iun2,iun3,maxni,maxnj,i,datev,niyy,vesion_uencode
   integer :: nlis,lislon, key, ni1,nj1,nk1,ni,nj,err,sindx_yin,sindx
   integer :: yinlislon,yanlislon,taclislon,j
   type(fst_record), dimension(:), allocatable :: liste,niv,yinliste,yanliste,tacliste

   real  :: xlat1,xlon1, xlat2,xlon2
   real, dimension(:),allocatable, target :: champ, yy
   real(REAL64) :: nhours

   type(fst_file)   :: file_in,file_out,file_yin,file_yan
   type(fst_record) :: record
   type(fst_query)  :: query, query2
   logical          :: success
   !-------------------------------------------------------------------

   app_ptr=app_init(0,'yyencode',YYENCODE_VERSION,'',BUILD_TIMESTAMP)
   call app_start()

   NPOS = 1
   call CCARD(LISTEc,DEF,VAL,3,NPOS)
   yin_S = val(1)
   yan_S = val(2)
   out_S = val(3)

   iun1 = 0 ; iun2 = 0 ; iun3 = 0

   if (.not. file_out%open(out_S,'RND+R/W')) then
      call app_log(APP_ERROR,'Unable to open '//trim(out_S))
      call qqexit(app_end(-1))
   endif 

   if (.not. file_yin%open(yin_S,'RND')) then
      call app_log(APP_ERROR,'Unable to open '//trim(yin_S))
      call qqexit(app_end(-1))
   endif 

   if (.not. file_yan%open(yan_S,'RND')) then
      call app_log(APP_ERROR,'Unable to open '//trim(yan_S))
      call qqexit(app_end(-1))
   endif 

   nlis = file_yin%get_num_records()
   allocate(liste(nlis),niv(nlis),yinliste(nlis),yanliste(nlis),tacliste(nlis))

   query = file_yin%new_query()
   lislon=query%find_all(liste)

   if (lislon < 1) then
      call app_log(APP_INFO,'Nothing to do')
      app_status=app_end(-1)
      stop
   endif

   maxni=0 ; maxnj=0
   do i=1,lislon
      maxni= max(maxni,liste(i)%ni)
      maxnj= max(maxni,liste(i)%nj)
   end do

   query = file_yin%new_query(nomvar='>>  ')
   yinlislon=query%find_all(yinliste)
   call query%free()

   if (yinlislon == 0) then
      call app_log(APP_ERROR,'YIN positionnal parameters >> not available')
      call qqexit(app_end(-1))
   endif

   query = file_yin%new_query(nomvar='^^  ')
   taclislon=query%find_all(tacliste)
   call query%free()

   if (taclislon == 0) then
      call app_log(APP_ERROR,'YIN positionnal parameters ^^ not available')
      call qqexit(app_end(-1))
   endif

   if (taclislon.ne.yinlislon) then
      call app_log(APP_ERROR,'YIN positionnal parameters not all available')
      call qqexit(app_end(-1))
   endif

   query = file_yan%new_query(nomvar='>>  ')
   yanlislon=query%find_all(yanliste)
   call query%free()

   if (yanlislon == 0) then
      call app_log(APP_ERROR,'YAN positionnal parameters >> not available')
      call qqexit(app_end(-1))
   endif

   query = file_yan%new_query(nomvar='^^  ')
   taclislon=query%find_all(tacliste)
   call query%free()

   if (taclislon == 0) then
      call app_log(APP_ERROR,'YAN positionnal parameters ^^ not available')
      call qqexit(app_end(-1))
   endif

   if (taclislon.ne.yanlislon) then
      call app_log(APP_ERROR,'YAN positionnal parameters not all available')
      call qqexit(app_end(-1))
   endif

   allocate (champ(maxni*2*maxni))

   DO_FLDS: do i=1,lislon

      datev  = -1
      if (dte > 0) datev=liste(i)%datev

      query = file_yan%new_query(datev=datev,ip1=liste(i)%ip1,ip2=liste(i)%ip2,ip3=liste(i)%ip3,typvar=liste(i)%typvar,nomvar=liste(i)%nomvar)
      success=query%find_next(record)
      nvar=trim(liste(i)%nomvar)

      if (nvar == '!!' .or. nvar == '>>' .or. nvar == '^^' .or. nvar == 'META') then
         if (nvar == '!!') then
            success = liste(i)%read()
            success = file_out%write(liste(i))
         endif
         if (nvar == 'META') then
            success = liste(i)%read()
            success = file_out%write(liste(i))
            success = file_yan%read(record,datev=datev,typvar=liste(i)%typvar,nomvar=nvar)
            success = file_out%write(record)
         endif
      else
         if (.not. success) then
            call app_log(APP_ERROR,'Corresponding YAN variable "\\a\\" NOT FOUND')
            call qqexit(app_end(-1))
         endif
         g3a  =  1                            ! points de masse
         g2a = liste(i)%ig2
         if (nvar == 'UT1' .or. nvar == 'URT1')  then
            g3a = 2 ! points U
            g2a= liste(i)%ig2-1
         endif
         if (nvar == 'VT1' .or. nvar == 'VRT1')  then
            g3a = 3 ! points V
            g2a= liste(i)%ig2-2
         endif

         success = liste(i)%read(data=c_loc(champ))
         success = record%read(data=c_loc(champ(liste(i)%ni*liste(i)%nj+1)))
         nj=liste(i)%nj
         liste(i)%nj=liste(i)%nj*2
         liste(i)%ig2=g2a
         liste(i)%ig3=g3a
         liste(i)%grtyp='U'
         success = file_out%write(liste(i),data=c_loc(champ))
         liste(i)%nj=nj
      endif

   enddo DO_FLDS

   !Read and write out new grid descriptors'

   DO_TICTACS: do j=1,yinlislon

      call cigaxg('E', xlat1,xlon1, xlat2,xlon2, yinliste(j)%ig1,yinliste(j)%ig2,yinliste(j)%ig3,yinliste(j)%ig4)
      do i=1,lislon
         nvar=trim(liste(i)%nomvar)
         ni=liste(i)%ni
         nj=liste(i)%nj

         if(nvar.ne.'!!' .and. nvar.ne.'META' .and.  nvar.ne.'>>' .and. nvar.ne.'^^') then
            if (liste(i)%ig1 == yinliste(j)%ip1 .and. liste(i)%ig2 == yinliste(j)%ip2 .and. liste(i)%ig3 == yinliste(j)%ip3) then
               if (nvar == 'UT1' .or. nvar  ==  'URT1') then
                  p2a = liste(i)%ig2-1
                  p3a = 2
               else if (nvar == 'VT1' .or. nvar == 'VRT1') then
                  p2a = liste(i)%ig2-2
                  p3a = 3
               else
                  p2a = liste(i)%ig2
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

      success = yinliste(j)%read(data=c_loc(yy(sindx+10   )))
      success = tacliste(j)%read(data=c_loc(yy(sindx+10+ni)))

      yy(sindx+2) = yy(sindx+10      )
      yy(sindx+3) = yy(sindx+ 9+ni   )
      yy(sindx+4) = yy(sindx+10+ni   )
      yy(sindx+5) = yy(sindx+ 9+ni+nj)
      sindx_yin= sindx

      !YAN
      call cigaxg('E', xlat1,xlon1, xlat2,xlon2, yanliste(j)%ig1,yanliste(j)%ig2,yanliste(j)%ig3,yanliste(j)%ig4)

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

      record%data=c_loc(yy)
      record%nomvar='^>'
      record%typvar='X'
      record%etiket='YYG_UE_GEMV4'
      record%grtyp=family_uencode_S
      record%ni=niyy
      record%nj=1
      record%nk=1
      record%ip1=yinliste(j)%ip1
      record%ip2=p2a
      record%ip3=p3a
      record%dateo=0
      record%deet=0
      record%npas=0
      record%data_type=5
      record%pack_bits=32
      record%ig1=vesion_uencode
      record%ig2=0
      record%ig3=0
      record%ig4=0
      success = file_out%write(record)

      deallocate(yy, stat=err)
   enddo DO_TICTACS

   success = file_yin%close()
   success = file_yan%close()
   success = file_out%close()

   app_status=app_end(-1)

end program yyencode
