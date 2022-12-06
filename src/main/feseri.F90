!# feseri Needs: SERGDIM .F90, serie.F90

!# SERGDIM Needs: series.cdk, serdim.F90, seralc2.F90
!# series.cdk Needs: -
!# serdim Needs: -
!# seralc2 Needs:-

!# serie2 Needs: INDSERI.F90, SERPAT2, SERFIN, SERACC, SERECRI3
!# INDSERI Needs: -
!# SERPAT2 Needs: -
!# SERFIN Needs: PLLWFGFW
!# SERACC Needs: acclist.cdk
!# SERECRI3 Needs: -
!# PLLWFGFW Needs: -
!# acclist.cdk Needs: 

!# Needs: sergdim.F90 serie.F90 series.cdk serdim.F90 seralc2.F90 indseri.F90
!         serpat2.F90 serfin.F90 seracc.F90 serecri2.F90 pllwfgfw.F90 acclist.cdk

!-------------------------------------- LICENCE BEGIN ------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------

subroutine FESERI
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   !Author
   !          M. Lepine   (Mar 87)
   !
   !Revision
   ! 001      B. Reid  (June 89) - Zonal diagnostics
   ! 002      N. Brunet  (May90)
   !                Standardization of thermodynamic functions
   ! 003      N. Brunet  (May91)
   !                New version of thermodynamic functions
   !                and file of constants
   ! 004      B. Bilodeau  (July 1991)- Adaptation to UNIX
   ! 005      G. Pellerin (April 1992) -Adaptation to PASTEMP,
   !                            deleted calls to zonal diagnostics.
   ! 006      G. Pellerin (April 1994) -Open sequential file with
   !                                    options +FTN +UNF
   ! 007      B. Bilodeau (June 1997) - IBM 32 bit to IEEE 32 bit
   !                                    format conversion
   ! 008      B. Bilodeau (July 1998) - Automate IBM32 to IEEE conversion
   ! 009      B. Bilodeau (March 2007) - Compression
   ! 010      R. McTaggart-Cowan (Apr 2009) - Add options for 64 bit time
   !                                    (hour) output and no diagnostics
   !
   !Object
   !          to take the raw binary "T" files created by FEMAIN
   !          and convert them to RPN standard file format.
   !          The output is the time-series "S" file.
   !
   !Arguments
   !          None.
   !
   !FILES
   !     TAPE64  input file containing POINTS,SURFACE,PROFILS
   !     TAPE35  reformatted standard random output for graphics
   !             on time series (PASTEMP)

   integer, external :: SERGDIM3

   integer status,inbr, tmpunit
   character(len=128) :: defo(9),listl(9),lfn(9)
   logical echot,compress,hour64,diag
   integer inpunit,serstd,ier
   integer junk,mstat,msurf,mprof,nig, &
        nk_hybm,nk_hybt
   save status,inpunit,serstd,listl,defo,lfn
   data status /0/
   data inpunit,serstd / 64 , 35 /
   data LISTL/'ISERIAL.','OMSORTI.','I','L','DATE','ECHOT', &
        'COMPRESS','HOUR64','NODIAG'/
   data DEFO/'TAPE64','TAPE35','$IN','$OUT','OPRUN','OUI','OUI', &
        'OUI','NON'/
   data LFN/'TAPE64','TAPE35','$IN','$OUT','NON','NON','NON', &
        'NON','OUI'/

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !
   !     LISTL = POSITION USAGER ISERIAL( TAPE64 , SEQUENTIEL ) ,
   !             OMSORTI( TAPE35 , STD ) , I , L , DATE
   !     DEFO = LISTE DES DEFAUTS ISERIAL ,  OMSORTI , I , L , DATE
   !                            TAPE64 , TAPE35 , $IN , $OUT , OPRUN
   !     LFN = LISTE QUE L USAGER PROPOSE POUR REMPLACER
   !     8 = NOMBRE DE LFN

   call CCARD(LISTL, DEFO, LFN, size(LISTL), -1)
   IER = FNOM(INPUNIT, LFN(1), 'SEQ+FTN+UNF', 0)
   IER = FNOM(SERSTD, LFN(2), 'RND', 0)
   tmpunit = 5
   IER = FNOM(tmpunit, LFN(3), 'SEQ', 0)
   tmpunit = 6
   IER = FNOM(tmpunit, LFN(4), 'SEQ', 0)
   ECHOT = LFN(6).eq.'OUI'
   COMPRESS = LFN(7).eq.'OUI'
   HOUR64 = LFN(8).eq.'OUI'
   DIAG = LFN(9).eq.'OUI'

   JUNK = EXDB('FESERI', 'V4.0', LFN(5))

   status= sergdim3(inpunit,mstat,msurf,mprof,nk_hybm,nk_hybt)
   if (status.lt.0) stop 1

   if (status == 0) then

      nig=4

      call serdbu()

      inbr = fstouv(serstd , 'RND')

      call serie3(inpunit,serstd,status,echot,compress,hour64, &
           nig,nk_hybm,nk_hybt)

      if (status .eq. 0) then
         inbr = fstfrm(serstd)
         ier = fclos(inpunit)
         goto 40
      endif

   endif

   if (status == 1) then
      write(6,*) '----AUCUN POINT POUR LES SERIES'
   else if (status == 2) then
      write(6,*) '----PAS DE variable de surface POUR LES SERIES'
   else
      write(6,*) '----PAS DE profile POUR LES SERIES'
   endif
40 continue

   JUNK = EXFIN ( 'FESERI' , 'FIN NORMALE' , 'NON' )

   !     ---------------------------------------------------------------
   return
end subroutine FESERI
