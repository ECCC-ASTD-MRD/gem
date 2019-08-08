!/@
subroutine test_fstinf_datev()
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04
!@/
   integer, external :: fnom,fstouv,fstinf,fstfrm,fclos

   real(8),parameter :: SEC_PER_HR = 3600.d0

   integer :: istat,fileid,datev,searchdate,dt,key,ni1,nj1,nk1
   real(8) :: nhours_8
   character(len=512) :: dfiles_S,filename_S,datev_S,searchdate_S
   ! ---------------------------------------------------------------------
   dfiles_S = '/home/ordenv/ssm-domains9/release/gem-data_4.0.1/gem-data_4.0.1_all/share/data/dfiles'
   filename_S = trim(dfiles_S)//'/bcmk/2009042700_000'
   fileid = 0
   istat = fnom(fileid,filename_S,'STD+RND+OLD+R/O',0)
   istat = min(fstouv(fileid,'RND'),istat)
   if (istat < 0 .or. fileid <= 0) then
      print *,'ERROR: open file failed'
      return
   endif

   datev_S = '20090427.000000'
   call datp2f(datev,datev_S)

   print *,'Searching for TT at dateo+dt should not find anything'
   print *,'FAIL means fstinf found a TT record even with wrong date'
   print '(a5,a3,x,2a16,a)',' ','dt','datev0          ','datev0+dt       ','dble(dt)'
   do dt=-10,40,5
      nhours_8 = dble(dt)/SEC_PER_HR
      call incdatr(searchdate,datev,nhours_8)
      call datf2p(searchdate_S,searchdate)
      key = fstinf(fileid,ni1,nj1,nk1,searchdate,' ',-1,-1,-1,' ','TT')
      if (key >= 0 .and. dt/=0) then
         print '(a5,I3,x,2a16,f)','FAIL ',dt,datev_S(1:15),searchdate_S(1:15),nhours_8
      else
         print '(a5,I3,x,2a16,f)','OK   ',dt,datev_S(1:15),searchdate_S(1:15),nhours_8
      endif
   enddo

   istat = fstfrm(fileid)
   istat = fclos(fileid)
   ! ---------------------------------------------------------------------
   return
end subroutine test_fstinf_datev
