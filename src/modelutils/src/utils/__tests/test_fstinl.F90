subroutine test_fstinl
   implicit none
!#include <rmnlib_basics.hf>
   integer,external :: fstopc,fnom,fstouv,fstinf,fstinl,ip1_all
   logical,parameter :: RMN_OPT_SET = .false.
   integer,parameter :: NMAX = 9999
   character(len=*),parameter :: FILENAME = '/home/ordenv/ssm-domains9/release/gem-data_4.2.0/gem-data_4.2.0_all/share/data/dfiles/bcmk/2009042700_000'
   character(len=12) :: etk,vartype,varname
   integer :: istat, fileid, key, ni,nj,nk,ip1,ip2,ip3,datev,nkeys,kind
   integer :: keylist(NMAX)
   real :: zp1
   !----
   istat = fstopc('MSGLVL','DEBUG',RMN_OPT_SET)
   fileid = 0
   istat = fnom(fileid,FILENAME,'RND+OLD+R/O',0)
   istat = fstouv(fileid,'RND')

   datev = 354514400
   ip1 = 0; ip2 = -1 ; ip3 = -1
   etk = " " ; vartype = " "
   varname = "pr0"

   key = fstinf(fileid,ni,nj,nk,datev,etk,ip1,ip2,ip3,vartype,varname)
   if (key >= 0) then
      print *,'(fstinf) Found: ',trim(varname)
   else
      print *,'(fstinf) Not Found: ',trim(varname)
   endif

   datev = 354514400
   zp1 = 0.
   kind = 3
   ip1 = ip1_all(zp1,kind)
   key = fstinf(fileid,ni,nj,nk,datev,etk,ip1,ip2,ip3,vartype,varname)
   if (key >= 0) then
      print *,'(fstinf) Found: ',trim(varname)
   else
      print *,'(fstinf) Not Found: ',trim(varname)
   endif

   datev = -1
   zp1 = 0.
   kind = 3
   ip1 = ip1_all(zp1,kind)
   istat = fstinl(fileid,ni,nj,nk,datev,etk,ip1,ip2,ip3,vartype,varname, &
        keylist,nkeys,NMAX)
   if (nkeys > 0) then
      print *,'(fstinl) Found: ',trim(varname),nkeys
   else
      print *,'(fstinl) Not Found: ',trim(varname)
   endif

   datev = -1
   ip1 = 0
   istat = fstinl(fileid,ni,nj,nk,datev,etk,ip1,ip2,ip3,vartype,varname, &
        keylist,nkeys,NMAX)
   if (nkeys > 0) then
      print *,'(fstinl) Found: ',trim(varname),nkeys
   else
      print *,'(fstinl) Not Found: ',trim(varname)
   endif

!!$   istat = fstinl(fileid,ni,nj,nk,datev,etk,ip1,ip2,ip3,vartype,varname, &
!!$        keylist,nkeys,NMAX)
!!$   if (nkeys > 0) then
!!$      print *,'(fstinl) Found: ',trim(varname),nkeys
!!$   else
!!$      print *,'(fstinl) Not Found: ',trim(varname)
!!$   endif
   !----
   return
end subroutine test_fstinl
