
subroutine fstinltest
   use, intrinsic :: iso_fortran_env, only: REAL64
   use cmcdate_mod, only: cmcdate_toprint
   implicit none

#include <rmnlib_basics.hf>


   real(REAL64), parameter :: SEC_PER_HR = 3600.d0
   integer,parameter :: NMAX = 9999

   integer :: fileid, keylist(NMAX),nkeys,istat,k
   integer :: ni,nj,nk, datev, searchdate, &
        dateo,deet,npas,nbits, datyp, ip1, ip2, ip3, &
        ig1, ig2, ig3, ig4, swa, lng, dltf, ubc, extra1, extra2, extra3
   character(len=1) :: grtyp_S
   character(len=2) :: typvar_S
   character(len=4) :: nomvar_S
   character(len=12):: etiket_S
   character(len=1024) :: filename_S, type_S
   real(REAL64) :: nhours_8

   filename_S = "/space/hall3/sitestore/eccc/mrd/rpndat/dja001/30_sec_outputs/2020011300_001"
   type_S = 'RND+OLD+R/O'
   fileid = 0
   istat = fnom(fileid,filename_S,type_S,0)
   istat = fstouv(fileid,'RND')
   
   searchdate = 439035236  !# 20200113.00023000
   istat = fstinl(fileid,ni,nj,nk,searchdate,' ', -1,-1,-1,' ',' ', &
        keylist,nkeys,NMAX)

   print *,'Searching for: ',searchdate,cmcdate_toprint(searchdate)
   do k=1,nkeys
      istat = fstprm(keylist(k), dateo,deet,npas, ni,nj,nk, &
           nbits, datyp, ip1, ip2, ip3, &
           typvar_S, nomvar_S, etiket_S, &
           grtyp_S, ig1, ig2, ig3, ig4, swa, lng, dltf, &
           ubc, extra1, extra2, extra3)
      nhours_8 = (DBLE(deet)*DBLE(npas))/SEC_PER_HR
      call incdatr(datev,dateo,nhours_8)
      print *,keylist(k),':',npas,':',datev,cmcdate_toprint(datev)
   enddo

   istat = fstfrm(fileid)
   istat = fclos(fileid)

   return
end subroutine fstinltest
