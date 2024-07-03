#include "modelutils_build_info.h"

program flipit
   use app
   use rmn_fst24
   implicit none

   character(len=1024) :: LISTEc(2), DEF(2), VAL(2)
   data LISTEc /'i.'  , 'o.'/
   data VAL    /'null'  , 'null'/
   data DEF    /'null'  , 'null'/

   character(len=1024) :: out_file,in_file
   character(len=4)    :: liste_var(50000)

   integer, parameter ::nlis = 10024

   type(fst_file)   :: filei,fileo
   type(fst_record) :: record,liste2(nlis)
   type(fst_query)  :: query, query2
   logical          :: success
   
   integer liste(nlis)

   integer NPOS, kd,kf,kp
   integer i,j,k,m
   integer ni,lislon,lislon2,cnt2,ipcode,ipkind
   integer err,cnt,liste_ip1(50000),liste_ip3(50000),nvar
   real, dimension (:), allocatable, target :: wk1
   real pcode

   !
   !-------------------------------------------------------------------
   !
   app_ptr=app_init(0,'flipit',FLIPIT_VERSION,'',BUILD_TIMESTAMP)
   call app_start()

   NPOS = 1
   call CCARD(LISTEc,DEF,VAL,2,NPOS)

   in_file = val(1)
   out_file= val(2)

   if (.not. filei%open(in_file,'RND+OLD')) then
      call app_log(APP_ERROR,'Unable to open '//trim(in_file))
      call qqexit(app_end(-1))
   endif 

   if (.not. fileo%open(out_file,'RND+R/W')) then
      call app_log(APP_ERROR,'Unable to open '//trim(out_file))
      call qqexit(app_end(-1))
   endif

   query = filei%new_query()
   liste_var='' ; cnt=0
   do while(query%find_next(record))
      if ((record%nomvar /= '>>').and.(record%nomvar /= '^^').and.(record%nomvar /= '!!')) then
         if (any(liste_var(1:cnt)==record%nomvar)) cycle
         cnt= cnt+1
         liste_var(cnt)= record%nomvar
      endif
   end do
   call query%free()
   nvar= cnt
   
   do i=1,nvar
      query = filei%new_query(nomvar=liste_var(i))
      liste_ip3=0 ; cnt=0
      do while(query%find_next(record))
         if (any(liste_ip3(1:cnt)==record%ip3)) cycle
         cnt= cnt+1
         liste_ip3(cnt)= record%ip3
      end do

      call convip_plus ( ipcode, pcode, ipkind, 0, ' ', .false. )

      do j=1,cnt
         query2 = filei%new_query(nomvar=liste_var(i),ip3=liste_ip3(j))
         lislon2=query2%find_all(liste2)
         ni=liste2(1)%ni

         allocate(wk1(ni*lislon2))
         call record_sort_ip1(liste2,liste_ip1,lislon2)
         call convip_plus (liste_ip1(lislon2), pcode, ipkind, -1, ' ', .false. )
        
         kd=lislon2; kf=1; kp=-1
         cnt2=0
         do k=kd,kf,kp
            cnt2=cnt2+1
            success = filei%read(record,data=c_loc(wk1((cnt2-1)*ni+1)),nomvar=liste_var(i),ip1=liste_ip1(k),ip3=liste_ip3(j))
         end do

         record%npas=liste_ip3(j)
         record%nj=lislon2
         record%ip1=0
         record%ip3=liste_ip3(j)
         record%grtyp='X'
         record%ig1=0
         record%ig2=0
         record%ig3=0
         record%ig4=0
         success = fileo%write(record,data=c_loc(wk1))

         deallocate (wk1)
         call query2%free()
      end do
      call query%free()
   end do

   success = filei%close()
   success = fileo%close()

   app_status=app_end(-1)

end program flipit
