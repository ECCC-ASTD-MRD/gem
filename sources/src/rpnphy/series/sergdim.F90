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

integer function sergdim3 ( inpunit, mstat, msurf, mprof,&
     nk_hybm, nk_hybt )
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer inpunit,mstat,msurf,mprof,nk_hybm,nk_hybt

   !@Arguments
   !          - Input -
   ! INPUNIT  unit number attached to file to read
   !          - Output -
   ! MSTAT    maximum number of stations contained in a file
   ! MSURF    maximum number of surfaces contained in a file
   ! MPROF    maximum number of profiles contained in a file
   ! ERREUR   .TRUE. to indicate reading file error (on INPUNIT)
   !          .FALSE. to indicate no error on reading file INPUNIT

   include "series.cdk"

   INTEGER NK
   INTEGER NSTT,NSRF,NPRF,MSTT,MSRF,MPRF

   INTEGER SERDIM, DIMSERS, DIMSERP

   !     ---------------------------------------------------------------

   MXSRF=0 ; MXPRF=0 ; MXSTT=0 ; MXNVO=0
   NSTAT=0 ; NSURF=0 ; NPROF=0 ; initok=.false.
   SRWRI=0 ; TSVER=100 ; TSMOYHR=0

   sergdim3= -1

   REWIND (INPUNIT)
   read (inpunit,end=9110, err=9110) nstt,nsrf,nprf,nk
   read (inpunit,end=9110, err=9110)
   read (inpunit,end=9110, err=9110)
   read (inpunit,end=9110, err=9110)
   read (inpunit,end=9110, err=9110)
   read (inpunit,end=9110, err=9110) mstt,msrf,mprf,nk

   nk_hybm= nk ; nk_hybt = nk

   IF (MSTT .LE. 0) sergdim3 = 1
   IF (msrf .LE. 0) sergdim3 = 2
   IF (mprf .LE. 0) sergdim3 = 3
   if (sergdim3.gt.0) return

   mstat = max(nstt,mstt)
   msurf = max(nsrf,msrf)
   mprof = max(nprf,mprf)

   write (6,"(/' SERGDIM...MSTAT=',i7,' MSURF=',i7,' MPROF=',i7,' NK=',i7,/)") MSTAT,MSURF,MPROF,NK

   dimsers = serdim (mstat,msurf,1)
   dimserp = serdim (mstat,mprof,nk)

   allocate (jstat(mxstt),ijstat(mxstt,2),statnum(mxstt),&
        jstat_g(mxstt),istat_g(mxstt),kam(1),&
        name(mxstt*stn_string_length/4))

   allocate (sers(mxstt,mxsrf),serp(mxnvo,mxstt,mxprf))

   call seralc2 (1,1,nk)

   initok=.true.

   sergdim3= 0

   if (sergdim3.eq.0) return

   !     ---------------------------------------------------------------
9110 PRINT *, 'ERREUR READING UNIT', inpunit,' in sergdim'

   return
end function sergdim3
