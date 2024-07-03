#include "modelutils_build_info.h"

program yydecode
   use, intrinsic :: iso_fortran_env, only: REAL64
   use app
   use rmn_fst24

   implicit none
   include "rmnlib_basics.inc"

   !@author   M.Desgagne   -   Spring 2012

   character(LEN=1024) :: LISTEc(3), DEF(3), VAL(3)
   integer NPOS
   data LISTEc /'yin.'   , 'yan.' , 'i.'    /
   data VAL    /'/null'  , '/null', '/null' /
   data DEF    /'/null'  , '/null', '/null' /

   character(len=1024) :: yin_S,yan_S,in_S

   integer nhgrids
   integer, dimension (:,:,:), allocatable :: igs

   integer i,j,sindx
   integer nlis,ni,nj,err,lindex
   integer, dimension(:), allocatable :: liste,niv
   logical secondMETA_L

   real, dimension(:), pointer :: yy,champ

   type(fst_file)   :: file_in,file_yin,file_yan
   type(fst_record) :: record, desc
   type(fst_query)  :: query, query2
   type(C_PTR)      :: data
   logical          :: success
   !
   !-------------------------------------------------------------------
   !
   app_ptr=app_init(0,'yydecode',YYDECODE_VERSION,'',BUILD_TIMESTAMP)
   call app_start()

   NPOS = 1
   call CCARD(LISTEc,DEF,VAL,3,NPOS)
   yin_S = val(1)
   yan_S = val(2)
   in_S  = val(3)

   secondmeta_L=.false.

   if (.not. file_in%open(in_S,'RND+OLD')) then
      call app_log(APP_ERROR,'Unable to open '//trim(in_S))
      call qqexit(app_end(-1))
   endif 

   if (.not. file_yin%open(yin_S,'RND+R/W')) then
      call app_log(APP_ERROR,'Unable to open '//trim(yin_S))
      call qqexit(app_end(-1))
   endif 

   if (.not. file_yan%open(yan_S,'RND+R/W')) then
      call app_log(APP_ERROR,'Unable to open '//trim(yan_S))
      call qqexit(app_end(-1))
   endif 

   nlis= file_in%get_num_records()
   allocate(liste(nlis),niv(nlis))

   query = file_in%new_query(nomvar='^>  ')
   nhgrids=query%find_all()
   call query%rewind()

   allocate(igs(3,3,max(nhgrids,1)))
   i=1
   do while(query%find_next(record))
      igs(1,1,i)=record%ip1
      igs(2,1,i)=record%ip2
      igs(3,1,i)=record%ip3
   
      success=record%read()
      call record%get_data_array(yy)
      sindx  = 6

      ni = nint(yy(sindx  ))
      nj = nint(yy(sindx+1))

      call cxgaig('E', record%ig1,record%ig2,record%ig3,record%ig4, &
           yy(sindx+6), yy(sindx+7), yy(sindx+8), yy(sindx+9))
      !print *,'sindx+10+ni=',sindx+10+ni
      call set_igs2(igs(1,2,i), igs(2,2,i)   , &
           yy(sindx+10),yy(sindx+10+ni),ni,nj, &
           record%ig1,record%ig2,record%ig3,record%ig4, 1,ni,1,1,nj,1)
      igs(3,2,i)= igs(3,1,i)

      ! check if already written
      query2 = file_yin%new_query(nomvar='>>  ',ip1=igs(1,2,i),ip2=igs(2,2,i),ip3=igs(3,2,i))
      if (.not. query2%find_next()) then

         desc%data=c_loc(yy(sindx+10))
         desc%data_type=FST_TYPE_REAL_IEEE
         desc%data_bits=32
         desc%pack_bits=32
         desc%nomvar='>>'
         desc%typvar='X'
         desc%etiket='YYG_POSX'
         desc%grtyp='E'
         desc%ni=ni
         desc%nj=1
         desc%nk=1
         desc%ip1=igs(1,2,i)
         desc%ip2=igs(2,2,i)
         desc%ip3=igs(3,2,i)
         desc%ig1=record%ig1
         desc%ig2=record%ig2
         desc%ig3=record%ig3
         desc%ig4=record%ig4
         desc%deet=0
         desc%npas=0
         success=file_yin%write(desc)

         desc%data=c_loc(yy(sindx+10+ni))
         desc%nomvar='^^'
         desc%etiket='YYG_POSX'
         desc%ni=1
         desc%nj=nj
         success=file_yin%write(desc)
      endif
      call query2%free()

      sindx = sindx+10+ni+nj
      call cxgaig('E', record%ig1,record%ig2,record%ig3,record%ig4, &
         yy(sindx+6), yy(sindx+7), yy(sindx+8), yy(sindx+9))

      call set_igs2(igs(1,3,i), igs(2,3,i)            , &
         yy(sindx+10),yy(sindx+10+ni),ni,nj, &
         record%ig1,record%ig2,record%ig3,record%ig4, 1,ni,1,1,nj,1)
      igs(3,3,i)= igs(3,1,i)

      query2 = file_yan%new_query(nomvar='>>  ',ip1=igs(1,3,i),ip2=igs(2,3,i),ip3=igs(3,3,i))
      if (.not. query2%find_next()) then
         desc%data=c_loc(yy(sindx+10))
         desc%data_type=FST_TYPE_REAL_IEEE
         desc%nomvar='>>'
         desc%etiket='YYG_POSX'
         desc%ni=ni
         desc%nj=1
         desc%nk=1
         desc%ip1=igs(1,3,i)
         desc%ip2=igs(2,3,i)
         desc%ip3=igs(3,3,i)
         desc%ig1=record%ig1
         desc%ig2=record%ig2
         desc%ig3=record%ig3
         desc%ig4=record%ig4
         success=file_yan%write(desc)

         desc%data=c_loc(yy(sindx+10+ni))
         desc%nomvar='^^'
         desc%etiket='YYG_POSX'
         desc%ni=1
         desc%nj=nj
         success=file_yan%write(desc)
      endif
      call query2%free()
      
      i=i+1
   end do

   call query%free()
   query = file_in%new_query()

   do while(query%find_next(record))

      if (trim(record%nomvar) == '^>') cycle
      if ((trim(record%nomvar) == '^^') .or. (trim(record%nomvar) == '>>') &
           .or. (trim(record%nomvar) == '!!') .or. (trim(record%nomvar) == 'META' )) then
         success = record%read()
         if (trim(record%nomvar) == 'META') then
            if (secondmeta_L) then
               success = file_yan%write(record)
            else
               success = file_yin%write(record)
               secondmeta_L=.true.
            endif
         else
            success = file_yin%write(record)
         endif
         if (trim(record%nomvar) == '!!' ) &
            success = file_yan%write(record)
      else

         lindex=0
         do j=1,nhgrids
          if ((igs(1,1,j).eq.record%ig1) .and. &
              (igs(2,1,j).eq.record%ig2) .and. &
              (igs(3,1,j).eq.record%ig3)) lindex= j
         end do

         success = record%read()
         call record%get_data_array(champ)
         ni=record%ni;
         nj=record%nj;
         if (lindex  >  0 .and. record%grtyp(1:1) == 'U') then
            
            record%nj=record%nj/2
            record%grtyp='Z'
            record%ig1=igs(1,2,lindex)
            record%ig2=igs(2,2,lindex)
            record%ig3=igs(3,2,lindex)
            success = file_yin%write(record)

            record%ig1=igs(1,3,lindex)
            record%ig2=igs(2,3,lindex)
            record%ig3=igs(3,3,lindex)
            success = file_yan%write(record,data=c_loc(champ(ni*nj/2+1)))
         else
            success = file_yin%write(record)
         endif
      endif
   end do

99 success = file_in%close()
   success = file_yin%close()
   success = file_yan%close()

   app_status=app_end(-1)

end program yydecode
