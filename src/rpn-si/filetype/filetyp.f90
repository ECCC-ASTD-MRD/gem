      program filetyp
      use app
      implicit none
      character *8 cle(10)
      character *4096 def(10),val(10)
      character *60 msgs(-1:37)
      integer wkoffit,ipos,ier,t(8),i

      data cle /'L ',8*'T','-.'/
      data def /'-L',8*'-1','  '/
      data val /'  ',8*'-1','  '/

      data msgs(-1) /'unknown'/
      data msgs(0)  /'error: undefined value, 0'/
      data msgs(1)  /'RPN standard random 89'/
      data msgs(2)  /'RPN standard sequential 89'/
      data msgs(3)  /'RPN standard sequential fortran 89'/
      data msgs(4)  /'CCRN'/
      data msgs(5)  /'CCRN-RPN'/
      data msgs(6)  /'BURP'/
      data msgs(7)  /'GRIB'/
      data msgs(8)  /'BUFR'/
      data msgs(9)  /'BLOK'/
      data msgs(10) /'FORTRAN sequential unformatted'/
      data msgs(11) /'COMPRESS'/
      data msgs(12) /'GIF89'/
      data msgs(13) /'GIF87'/
      data msgs(14) /'IRIS'/
      data msgs(15) /'JPEG'/
      data msgs(16) /'KMW'/
      data msgs(17) /'PBM'/
      data msgs(18) /'PCL'/   
      data msgs(19) /'PCX'/   
      data msgs(20) /'PDSVICAR'/
      data msgs(21) /'PM'/      
      data msgs(22) /'PPM'/     
      data msgs(23) /'POSTCRIPT'/
      data msgs(24) /'KMW_'/
      data msgs(25) /'RRBX'/    
      data msgs(26) /'SUNRAS'/  
      data msgs(27) /'TIFF'/    
      data msgs(28) /'UTAHRLE'/ 
      data msgs(29) /'XBM'/     
      data msgs(30) /'XWD'/     
      data msgs(31) /'ASCII'/     
      data msgs(32) /'BMP'/     
      data msgs(33) /'RPN standard random 98'/
      data msgs(34) /'RPN standard sequential 98'/
      data msgs(35) /'FICHIER NETCDF'/
      data msgs(36) /'FICHIER CMCARC v4'/
      data msgs(37) /'FICHIER CMCARC v5'/

      ipos = 0
      call ccard (cle,def,val,10,ipos)
      do i = 1, 8
        read(val(i+1),*)t(i)
      enddo
      ier = wkoffit(val(10))
      if (ier .gt. -1) then
        if(t(1) .ne. -1) then                    ! t option active, set ier to 0 if type in requested list
          do i = 1,8
            if(ier == t(i)) ier = 0
          enddo
          if(ier .ne. 0) ier = 1
        else
          write(6,66) 'File type is ',msgs(ier)  ! print message and return file type     
        endif
      else if (ier .eq. -1) then
         call system('file '//trim(val(1))//' -m $ARMNLIB_DATA/magic.extra:/usr/share/misc/magic '//trim(val(10)))

!         call system('file '//val(1))
      else if (ier .eq. -2) then
        call app_log(APP_ERROR,'File empty')
      else if (ier .eq. -3) then
        call app_log(APP_ERROR,'File does not exist or can not open')
      else
        call app_log(APP_ERROR,'File corrupt')
      endif

 66   format(/,a,a60,/)


      call qqexit(ier)
!      call qqexit(and(127,ier))
!      call qqexit(and(255,ier))
      stop
      end
