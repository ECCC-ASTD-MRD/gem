subroutine flipit()
   implicit none

   character(len=1024) :: LISTEc(2), DEF(2), VAL(2)
   data LISTEc /'i.'  , 'o.'/
   data VAL    /'null'  , 'null'/
   data DEF    /'null'  , 'null'/

   character(len=1024) :: out_file,in_file
   character(len=12)   :: etiket, version
   character(len=4)    :: nomvar,liste_var(50000)
   character(len=2)    :: typvar
   character(len=1)    :: grtyp

   integer, parameter ::nlis = 10024
   integer liste(nlis), liste2(nlis), liste3(nlis)

   integer dateo,  deet, npas, nbits, NPOS, kd,kf,kp
   integer datyp, ig1, ig2, ig3, ig4
   integer extra1, extra2, extra3
   integer i,j,k,m,key1,key2,istat,ip1,ip2,ip3,dltf,ubc,swa,lng
   integer ni,nj,nk,lislon,lislon2,lislon3,cnt2,ipcode,ipkind
   integer err,iun1,iun2,cnt,liste_ip1(50000),liste_ip3(50000),nvar
   real, dimension (:), allocatable :: wk1
   real pcode
   integer fnom,fstouv,fstinl,fstprm,fstluk,exdb,exfin, &
        fstecr,fstfrm,fclos
   !
   !-------------------------------------------------------------------
   !
   version="1.1"
   err = exdb ('flipit',version, 'NON')

   NPOS = 1
   call CCARD(LISTEc,DEF,VAL,2,NPOS)

   in_file = val(1)
   out_file= val(2)

   iun1 = 0 ; iun2 = 0

   if (fnom(iun1,in_file,'RND+OLD',0) >= 0) then
      if (fstouv(iun1,'RND') < 0) then
         write (6,8001) trim(in_file)
         stop
      endif
   else
      write (6,8000) trim(in_file)
      stop
   endif

   if (fnom(iun2,out_file,'rnd',0) >= 0) then
      if (fstouv(iun2,'rnd') < 0) then
         write (6,8001) trim(out_file)
         stop
      endif
   else
      write (6,8000) trim(out_file)
      stop
   endif

   key1 = FSTINL (IUN1, NI, NJ, NK, -1, ' ', -1,-1,-1, ' ', &
        ' ',liste,lislon,nlis)
   liste_var='' ; cnt=0
   do m=1,lislon
      ISTAT= FSTPRM(liste(m), dateo, deet, npas, ni, nj, nk, nbits,&
           datyp, ip1, ip2, ip3, typvar, nomvar, etiket,&
           grtyp, ig1, ig2, ig3, ig4, swa, lng, dltf,&
           ubc, extra1, extra2, extra3)
      if ((nomvar /= '>>').and.(nomvar /= '^^').and.(nomvar /= '!!')) then
         if (any(liste_var(1:cnt)==nomvar)) cycle
         cnt= cnt+1
         liste_var(cnt)= nomvar
      endif
   end do
   nvar= cnt

   do i=1,nvar
      key1 = FSTINL (IUN1, NI, NJ, NK, -1, ' ', -1,-1,-1, ' ', &
           liste_var(i),liste,lislon,nlis)
      liste_ip3=0 ; cnt=0
      do m=1,lislon
         ISTAT= FSTPRM(liste(m), dateo, deet, npas, ni, nj, nk, nbits,&
              datyp, ip1, ip2, ip3, typvar, nomvar, etiket,&
              grtyp, ig1, ig2, ig3, ig4, swa, lng, dltf,&
              ubc, extra1, extra2, extra3)
         if (any(liste_ip3(1:cnt)==ip3)) cycle
         cnt= cnt+1
         liste_ip3(cnt)= ip3
      end do

      call convip_plus ( ipcode, pcode, ipkind, 0, ' ', .false. )

      do j=1,cnt
         key1 = FSTINL (IUN1, NI, NJ, NK, -1, ' ', -1,-1,liste_ip3(j),&
                        ' ',liste_var(i),liste2,lislon2,nlis)
         allocate(wk1(ni*lislon2))
         call sort_ip1 ( liste2, liste_ip1, lislon2 )
         call convip_plus (liste_ip1(lislon2), pcode, ipkind, -1, ' ', .false. )
         kd=lislon2; kf=1; kp=-1
         cnt2=0
         do k=kd,kf,kp
            cnt2=cnt2+1
            key2 = FSTINL (IUN1, NI, NJ, NK, -1, ' ', liste_ip1(k),&
               -1,liste_ip3(j),' ',liste_var(i),liste3,lislon3,nlis)
            err = fstluk (wk1((cnt2-1)*ni+1),liste3(1),ni,nj,nk)
         end do
         err = fstecr (wk1,wk1,-nbits,iun2,dateo, deet,&
              liste_ip3(j),ni,lislon2,1,0,ip2,liste_ip3(j),typvar,nomvar,etiket,&
              'X', 0,0,0,0, datyp, .false.)
         deallocate (wk1)
      end do
   end do

   err = fstfrm(iun1)
   err = fstfrm(iun2)

   err = fclos(iun1)
   err = fclos(iun2)

   err = exfin ('flipit',version, 'OK')

8000 format (/' Unable to fnom: '  ,a/)
8001 format (/' Unable to fstouv: ',a/)
   !
   !-------------------------------------------------------------------
   !
   return
end subroutine flipit
