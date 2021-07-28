      program woir

!      Revision 98.11 - M. Lepine - reload, detection mauvais type de fichier (xdfopn)
!               98.12 - M. Lepine - utilisation de ccard_arg
!               98.13 - M. Lepine - flush des buffers de stdout dans exfin
!               98.14 - M. Lepine - reload librmn_008
!               98.15 - M. Lepine - reload librmnbeta, correction impression dans convip
!               98.16 - M. Lepine - reload librmn_009
!               98.17 - M. Lepine - sept 2007 - reload librmnbeta, correction fichiers > 2G
!               98.18 - M. Lepine - sept 2007 - reload pour dates etendues 
!               98.19 - M. Lepine - sept 2008 - reload pour fichier cmcarc remote 
!               98.20 - M. Lepine - dec  2010 - reload avec librmn_012 et codebeta moduledate 
!               98.21 - M. Lepine - sept 2011 - reload avec librmn_012 et codebeta moduledate_711e, fstd98
!               98.22 - M. Lepine - juin 2012 - reload avec librmn_013
!               98.30 - M. Lepine - mars 2014 - reload avec librmn_014
!               98.31 - M. Lepine - juin 2014 - ajout de implicit none
!               98.32 - M. Lepine - Dec  2014 - reload avec librmn_015.1
!               98.33 - M. Lepine - Fev. 2015 - reload avec librmn_015.2
!               99.00 - M. Valin  - oct. 2015 - nouvelle version de print_std_parms (librmn_Alpha_016)

      implicit none
      integer, parameter :: ncle=5
      integer fnom,fstouv,fstvoi,fstfrm,exdb,exfin
      external fnom,fstouv,fstvoi,fstnbr,ccard,c_init_appl_var_table
      character *8192 ccard_arg, filename
      external exdb,exfin,ccard_arg
      integer ipos,ier,n

      character * 8 cles(ncle)
      character * 12 status
      character * 128 val(ncle), def(ncle)
      character(len=*), parameter :: VERSION='99.00'
      data cles / 'IMENT:','SEQ','STYLE','MOREHELP','V' /
      data def / '/dev/null','SEQ' ,'NINJNK+DATEV+LEVEL+IP1+GRIDINFO','MOREHELP',VERSION/
      data val / '/dev/null','RND' ,'NINJNK+DATEO+IP1+IG1234',2*' '/

      call c_init_appl_var_table()
      ipos = -1
      call ccard (cles,def,val,ncle,ipos)
      if (val(5) /= "") then
        ier = exdb('VOIR',VERSION,'NON')
        stop
      endif
      status='<<ERREUR>>'
      IF (val(4) .eq. 'MOREHELP') THEN
         print *,"*** VOIR CALLING SEQUENCE ***"
         print *
         print *,'-IMENT [scrap:scrap]'
         print *,'-SEQ [RND:SEQ]'
         print *,'-STYLE [NINJNK DATEO IP1 IG123:NINJNK DATEV LEVEL IP1 GRIDINFO]'
         print *,'   List of possible items for STYLE argument:'
         print *,'         NINJNK: display ni nj nk dimensions'
         print *,"          DATEO: display origin date"     
         print *,'     DATESTAMPO: display origin datetimestamp for the nostalgics'
         print *,'          DATEV: display valid date and stamp'
         print *,'          LEVEL: display vertical level'
         print *,'          IPALL: display full IP1/2/3 trio decoding'
         print *,'            IP1: display coded IP1 value'
         print *,'       GRIDINFO: display decoded grid information'
         print *,'         IG1234: display IG1 IG2 IG3 IG4 values'
         print *
         print *,'     The following items suppress variable printout'
         print *,'         NONOMV: suppress NOMV information'
         print *,'         NOTYPV: suppress TYPV information'
         print *,'         NOETIQ: suppress ETIQUETTE information'
         print *,'         NOIP23: suppress IP2, IP3 information'
         print *,'         NODEET: suppress DEET information'
         print *,'         NONPAS: suppress NPAS information'
         print *,'          NODTY: suppress DTY information'
         print *
         print *,'   Example #1: -style "ninjnk datev level"'       
         print *,'   Example #2: -style datev+level+ip1+notypv'
         print *,'   -style FULL displays'
         print *,'           NINJNK+DATEV+IPALL+IP1+GRIDINFO'
      else
         ier = exdb('VOIR',VERSION,'NON')
         filename=ccard_arg(cles(1))
         if(trim(filename)=="") filename = val(1)
!         print *,'Debug+ filename = ',trim(filename)
         ier = fnom(10,trim(filename),'STD+R/O+REMOTE'//val(2),0)
         if (ier .ge. 0) then
            N = fstouv(10,VAL(2))
            if(trim(VAL(3))=="FULL")  val(3)='NINJNK+DATEV+IPALL+NOIP23+GRIDINFO'
            if(N>=0) then
              ier = fstvoi(10,val(3))
              ier = fstfrm(10)
              status='O.K.'
            endif
         endif
         ier = exfin('VOIR',status,'NON')
      endif
      stop
      end
      
      character *128 function product_id_tag()
      product_id_tag='$Id: voir.f90 2014-01-30 08:00:00Z armnlib $'
      return
      end
