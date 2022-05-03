!-------------------------------------- LICENCE BEGIN -------------------------
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
!-------------------------------------- LICENCE END ---------------------------

!*@/
SUBROUTINE read_isccpdata()
   implicit none
!!!#include <arch_specific.hf>

!Author
!          B Dugas     (Apr 2006)

!Revision
!001       B. Dugas    (May 2006) - Try reading files from $EXECDIR

!Object
!          Read in XCW data needed by stochastic cloud generator

!*@/

#include "mcica.cdk"

      logical ex
      character modeln*10
      character(len=256) execdir,fe,fn,modelp
      integer fio, IER, io

!MODULES

      INTEGER  fclos,fnom,longueur
      EXTERNAL fclos,fnom,longueur,getenvc

!----------------------------------------------------------------------

! ALWAYS CHECK FOR DATA FILES IN $EXECDIR/

      call getenvc('EXECDIR',execdir )
      execdir = execdir(1:longueur(execdir)) // '/'

! OTHERWISE, OPEN THE FILES AND READ THE DATA IN DIRECTORY
!     '$' // $MODEL // '/dfiles/ISCCP_DATA_FILES'

      call getenvc('MODEL',modeln )
      call getenvc( modeln,modelp )
      modelp = modelp(1:longueur(modelp)) // '/dfiles/ISCCP_DATA_FILES/'

! START WITH tautab

      fn  ='tautab.formatted'
      fe  = execdir(1:longueur(execdir)) // fn

      inquire( file=fe,err=0001,iostat=io,exist=ex )

      if (ex) then
         fn = fe
      else
         fn  = modelp(1:longueur(modelp)) // fn
      endif

      goto 2

    1 fn  = modelp(1:longueur(modelp)) // fn

    2 fio = 0
      IER = fnom( fio, fn,'SEQ+FTN+FMT+OLD',0 )
      if (IER.lt.0) goto 901

      READ( fio,'(f30.20)',ERR=902 ) tautab
      IER = fclos( fio )

! END WITH invtau

      fn  ='invtau.formatted'
      fe  = execdir(1:longueur(execdir)) // fn

      inquire( file=fe,err=0003,iostat=io,exist=ex )

      if (ex) then
         fn = fe
      else
         fn  = modelp(1:longueur(modelp)) // fn
      endif

      goto 4

    3 fn  = modelp(1:longueur(modelp)) // fn

    4 fio = 0
      IER = fnom( fio, fn,'SEQ+FTN+FMT+OLD',0 )
      if (IER.lt.0) goto 901

      READ( fio,'(i10)',ERR=902 ) invtau
      IER = fclos( fio )

      RETURN

!----------------------------------------------------------------------

! FOUND ERRORS

901   call physeterror('read_isccpdata', 'Problem opening: '//trim(fn))
      return
902   call physeterror('read_isccpdata', 'Problem reading: '//trim(fn))
      return

!----------------------------------------------------------------------
      END
