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
!-------------------------------------- LICENCE END ---------------------------

subroutine serecri4(VS,VV,IUN,NSURF, &
     NPROF,NT,SURFACE,PROFILS,NUM,DLAT,DLON, &
     STORE, DATE, ETIKET, H0, NK, &
     COMPRESS)
   use tdpack
   implicit none
!!!#include <arch_specific.hf>
   integer NSURF,NPROF
   character(len=*) SURFACE(*),PROFILS(*)
   character NOMVAR*4,ETIKET*12
   integer   nk,NT,NUM,IUN,IECR
   real      STORE(*),DLAT,DLON,H0
   real      VS(NT,*),VV(NK,NT,*)
   integer   DATE,DEET
   logical   COMPRESS

   !@Author R. Benoit

   !@Revision
   ! 001      V.Alex.(Feb 87) Documentation
   ! 002      N. Brunet  (May90)
   !                Standardization of thermodynamic functions
   ! 003      N. Brunet  (May91)
   !                New version of thermodynamic functions
   !                and file of constants
   ! 004      B. Bilodeau  (July 1991)- Adaptation to UNIX
   ! 005      B. Bilodeau  (August 1992)   - Delete HR
   ! 006      G. Pellerin (April 1992) - Adaptation to PASTEMP,
   !          clean-up of code
   ! 007      B. Bilodeau (Jan 1997) - Calculations of TH moved
   !                from serecri to serdyn4; calculation of TW
   !                generalized for GEF.
   ! 008      K. Winger   (Apr 2006) - ETIKET from 8 to 12 characters
   ! 009      B. Bilodeau (Mar 2007) - Compression
   ! 010      R. McTaggart-Cowan (
   ! 011      L. Spacek   (Oct 2010) - Profils limited to nk-1 levels
   ! 012      V. Lee (Jan 2018) - series_write defines NK so no need to limit here.

   !@Object ECRIRE LES SERIES TEMPORELLES POUR UN POINT SUR FICHIER STANDARD
   !        write the time-series for one point in a standard file

   !@Arguments
   !          - Input -
   ! VS       time-serie values of surface variables requested
   ! VV       time-serie values of profile variables requested
   !          - Output -
   ! IUN      unit number attached to standard file
   !          - Input -
   ! NSURF    number of surface variables requested
   ! NPROF    number of profile variables requested
   ! NT       timestep number
   ! SURFACE  names of time serie surface variables requested
   ! PROFILS  names of time serie profile variables requested
   ! NUM      station number
   ! DLAT     latitude
   ! DLON     longitude of station
   ! STORE    work field
   ! S        sigma (or eta) levels
   ! PTOIT    pressure value at the top of the model
   ! ETATOIT  eta value at the top of model
   ! DATE     date
   ! ETIKET   label for the standard record
   ! H0       GMT time
   ! NK       vertical dimension
   ! F_DIAG   .TRUE. to compute diagnostics
   !          .FALSE. to not compute diagnostics

   integer, external :: FSTECR, INDSERI

   integer, parameter :: STDERR=0
   character(len=1) VT

   real, dimension(:,:), allocatable :: fldm1
   integer KSURF,KPROF,ERR
   integer NPAS,IP1,IP2,IP3,IG1,IG2,IG3,IG4,NPAK,DATYP
!!$      INTEGER IELAQ,IELAT,IELAP0
!!$      SAVE    IELAQ,IELAT,IELAP0

   ! N.B. Profils are what is defined in series_write nk = phydim_nk - 1

   allocate(fldm1(nk,nt),stat=err)
   if (err /= 0) then
      write(STDERR,*) 'Error allocating fldm1... aborting'
      stop
   endif

   ! PARAMETRES DE GRILLE CONFORMES A LA DOCUMENTATION DE PASTEMP
   VT   = 'T'
   DEET = 0
   NPAS = 0
   IP1  = 0
   IP2  = H0*100.
   IP3  = NUM
   IG1  = 0
   IG2  = 0
   IG3  = (DLAT+100.)*100.+0.5
   IG4  = DLON*100.+0.5

   do  KSURF=1,NSURF
      !     SURFACE VARIABLES ARE NOT COMPRESSIBLE
      if(SURFACE(KSURF).eq.'Z0')then
         NPAK=1
         DATYP=5
      else
         NPAK=-24
         DATYP=1
      endif
      NOMVAR=SURFACE(KSURF)
      if(NOMVAR.ne.'  ') then
         IECR = FSTECR(VS(1,KSURF),STORE,NPAK,IUN,DATE,DEET,NPAS,1,NT,1, &
              IP1,IP2,IP3,VT,NOMVAR,ETIKET,'+',IG1,IG2,IG3,IG4,DATYP,.false.)
      endif
   enddo

   NPAK=-24
   DATYP=1
   !     LOSSLESS DATA COMPRESSION
   if (COMPRESS) then
      NPAK=-32
      DATYP=133
   endif
   do  KPROF=1,NPROF
      NOMVAR=PROFILS(KPROF)
      if(NOMVAR.ne.'  ') then
         fldm1 = vv(1:nk,:,KPROF)
         IECR = FSTECR(fldm1,STORE,NPAK,IUN,DATE,DEET,NPAS,nk,NT,1, &
              IP1,IP2,IP3,VT,NOMVAR,ETIKET,'+',IG1,IG2,IG3,IG4,DATYP,.false.)
      endif
   enddo


   !  GARBAGE COLLECTION
   deallocate(fldm1,stat=err)
   if (err /= 0) write(STDERR,*) 'Error freeing fldm1'

   return
end subroutine serecri4
