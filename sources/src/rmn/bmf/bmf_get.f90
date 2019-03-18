!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2005  Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
 function bmf_get(nom,time1,time2,i4,r4,r8,imin,imax,jmin,jmax,kmin,kmax) &
         result(error)
!
!AUTEUR       Luc Corbeil (bmf_get)
!
!REVISION
! v208a V. Lee real*8 allocation corrected:NI*dtyp/10/4 -> NI*(dtyp/10/4)
!
!ARGUMENTS
!
!nom         variable name                            character*4
!time1       timestamp 1 (yyymmdd)                    integer
!time2       timestamp2 (hhmmsscc)                    integer
!i4          destination array (dtyp=40)              integer
!r4          destination array (dtyp=41)              real
!r8          destination array (dtyp=81)              real*8
!(ijk)min    Size of destination array (lower bounds) integer
!(ijk)max    Size of destination array (upper bounds) integer
!
!i4 , r4 ou r8 will be considered as array of shape ( imin:imax
! , jmin:jmax , kmin:kmax )
!
!DESCRIPTION
! Routine that puts in the appropriate array the data retrieved by
! bmf_gobe. The use of dummy arguments for unneeded arrays is 
! possible if "dtyp" is known in advance. The array must be shaped
! by following those rules:( imin:imax , jmin:jmax , kmin:kmax ) where
!
! imin <= 1 < ni <= imax
! jmin <= 1 < nj <= jmax
! kmin <= 1 < nk <= kmax,
!
! ni,nj,nk are the attributes of the corresponding field (see 
! bmf_write) Then, it is possible to... (here, a,b >=0):
!
! ...write a field of size (1:nx,1:ny,1:nz) with bmf_write 
! and to retrieve it with bmf_get in a larger array
!
! .. write a field of (1-a : nx+a , 1:ny , 1:nz ) with bmf_write
! (call bmf_write with ni = nx + 2a, istart=1, iend=nx+2a) and to 
! read it in an array of (minx: maxx , miny:maxy , minz,maxz ) 
! where minx <=1-a <=nx+2a <= maxx (call bmf_get with 
! imin=minx+a , imax= maxx+a)
!
! ... write a field of (1+a : nx-b,1:ny,1:nz) (obviously we have
! 1+a < nx -b) with a standard call to bmf_write. We read the 
! field in an array of (minx: maxx , miny:maxy , minz,maxz ) 
! with a normal call to  bmf_get.
!
 use bmf_mod
 implicit none
  integer ni,nj,nk,error,time1,time2
  character*4 nom
  type(bmf_liste), pointer :: champ_courant
  integer imin,imax,jmin,jmax,kmin,kmax
  integer i4(imin:imax,jmin:jmax,kmin:kmax)
  real*4 r4(imin:imax,jmin:jmax,kmin:kmax)
  real rr4(imin:imax,jmin:jmax,kmin:kmax)
  real*8 r8(imin:imax,jmin:jmax,kmin:kmax)
  integer r8i(2*(imin-1)+1:2*imax,jmin:jmax,kmin:kmax)
  pointer(r8i_,r8i)
  pointer(rr4_,rr4)
  integer i,j,k
  integer istart,iend,jstart,jend,kstart,kend
  integer indice
  integer dtyp
  logical trouve
  integer bmf_get2

   error=0
   r8i_=loc(r8(imin,jmin,kmin))
   rr4_=loc(r4(imin,jmin,kmin))
   trouve=.false.
   champ_courant=>liste
  if((imin.gt.1).or.(jmin.gt.1).or.(kmin.gt.1)) then
      write(*,*) 'ERROR BMF_GET: IMIN OR JMIN OR KMIN .GT. 1'
      error=1
      return
  endif
  do
    if(.not.associated(champ_courant)) EXIT
      if((champ_courant%bmf_champ%nom.eq.nom).and. &
         (champ_courant%bmf_champ%time1.eq.time1).and. &
         (champ_courant%bmf_champ%time2.eq.time2)) then
           trouve = .true.
           dtyp=champ_courant%bmf_champ%dtyp
           istart=champ_courant%bmf_champ%istart
           jstart=champ_courant%bmf_champ%jstart
           kstart=champ_courant%bmf_champ%kstart
           iend=champ_courant%bmf_champ%iend
           jend=champ_courant%bmf_champ%jend
           kend=champ_courant%bmf_champ%kend
           ni=champ_courant%bmf_champ%ni
           nj=champ_courant%bmf_champ%nj
           nk=champ_courant%bmf_champ%nk
           if((imax.lt.ni).or.(jmax.lt.nj).or.(kmax.lt.nk)) then
             write(*,*) 'ERROR BMF_GET: IMAX OR JMAX OR KMAX .LT. EXPECTED'
             write(*,*) 'FOR VARIABLE ',nom
	     write(*,*) 'ni=',ni,'imax=',imax
	     write(*,*) 'nj=',nj,'jmax=',jmax
	     write(*,*) 'nk=',nk,'kmax=',kmax
             error=1
             return
           endif
           indice=0
           if(dtyp.eq.bmf_real4) then
              call bmf_copie(ni,nj,nk, &
                   champ_courant%bmf_champ%tableau,r4,imin,imax,&
                   jmin,jmax,kmin,kmax)
           else if (dtyp.eq.bmf_real8) then
              call bmf_copie(ni*(dtyp/10/4),nj,nk, &
                   champ_courant%bmf_champ%tableau,r8, &
                   (imin-1)*(dtyp/10/4)+1,imax*(dtyp/10/4),jmin,jmax,kmin,kmax)
           else if (dtyp.eq.bmf_integer4) then
              call bmf_copie(ni*(dtyp/10/4),nj,nk,&
                   champ_courant%bmf_champ%tableau,i4,&
                   imin,imax,jmin,jmax,kmin,kmax)
           else
              write(*,*) 'WARNING BMF_GET: Type dtyp ',dtyp,' non reconnu'
              error=1              
           endif
      endif
      champ_courant=>champ_courant%champ_suivant
   enddo
 if(.not.trouve) then
    write(*,*) 'WARNING BMF_GET: Variable ',nom,' non trouvee'
    error=1
 endif

 return
 end
