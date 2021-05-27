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
 function bmf_geth(nom,time1,time2,holename,hole_t1,hole_t2,i4,imin,imax,jmin,jmax,kmin,kmax,i0,j0) &
         result(error)
!
!AUTEUR       Luc Corbeil (bmf_geth)
!
!REVISION
!
!ARGUMENTS
!
!nom         variable name                            character*4
!time1       timestamp 1 (yyymmdd)                    integer
!time2       timestamp2 (hhmmsscc)                    integer
!holename    hole name                                character*4
!hole_t1     timestamp 1 (yyymmdd)                    integer
!hole_t2     timestamp2 (hhmmsscc)                    integer
!i4          destination array (dtyp=40)              integer
!r4          destination array (dtyp=41)              real
!r8          destination array (dtyp=81)              real*8
!(ijk)min    Size of destination array (lower bounds) integer
!(ijk)max    Size of destination array (upper bounds) integer
!i0,j0       Starting point for writing               integer
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
  integer, intent(IN) :: hole_t1,hole_t2,time1,time2,i0,j0
  character*4, intent(IN) :: nom,holename
  type(bmf_liste), pointer :: champ_courant
  integer, intent(IN) :: imin,imax,jmin,jmax,kmin,kmax
  integer error
  integer ni,nj,nk
  integer i4(imin:imax,jmin:jmax,kmin:kmax),hole(4)
  integer r4(imin:imax,jmin:jmax,kmin:kmax)
  real rr4(imin:imax,jmin:jmax,kmin:kmax)
  real*8 r8(imin:imax,jmin:jmax,kmin:kmax)
  integer r8i(2*(imin-1)+1:2*imax,jmin:jmax,kmin:kmax)
  pointer(r8i_,r8i)
  pointer(rr4_,rr4)
  integer i,j,k,ierr
  integer istart,iend,jstart,jend,kstart,kend
  integer indice,hole_i0,hole_isize,hole_j0,hole_jsize
  integer dtyp, dimx1,dimx2,dimx5,dimx6,nx0,ny0
  logical trouve

  integer bmf_get
  external bmf_get

  error=0
  r8i_=loc(r8(imin,jmin,kmin))
  rr4_=loc(r4(imin,jmin,kmin))
  trouve=.false.

  ierr = bmf_get(holename,hole_t1,hole_t2,hole,-1,-1,1,4,1,1,1,1)
  if(ierr.ne.0) then
     write(*,*) 'ERROR BMF_GETH: not found, hole ',holename
     error=1
     return
  endif
  hole_i0=hole(1)
  hole_isize=hole(2)
  hole_j0=hole(3)
  hole_jsize=hole(4)
  write(*,*) 'HOLE found, ',holename,hole
  champ_courant=>liste
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
           nx0=champ_courant%bmf_champ%hgrid
           ny0=champ_courant%bmf_champ%vgrid
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
             do k=1,nk
             do j=1,nj
                   write(*,*) 'J1',j,hole_j0,hole_jsize
                if((j.ge.hole_j0).and.(j.lt.hole_j0+hole_jsize-1)) then

                   write(*,*) 'J',j
                   do i=1,ni
                      indice=indice+1
                      r4(i,j,k)=champ_courant%bmf_champ%tableau(indice,1,1)
                   enddo
                endif
             enddo
             enddo
           endif
           if(dtyp.eq.bmf_real8) then
!              write(*,*) 
             do k=1,nk
             do j=1,nj
             do i=1,2*ni
!               indice=indice+1
!              write(*,*) i,j,k,champ_courant%bmf_champ%tableau(i,j,k)
               r8i(i,j,k)=champ_courant%bmf_champ%tableau(i,j,k)
             enddo
             enddo
             enddo
!             do k=1,nk
!             do j=1,nj
!              write(*,*) (r8(i,j,k),i=1,ni),'ho'  
!             enddo
!             enddo
           endif
           dimx1=max(nx0,istart)
           dimx2=min(nx0+ni-1,istart+hole_i0-1)
           dimx5=max(nx0,istart+(hole_i0+hole_isize-1))
           dimx6=min(nx0+ni-1,iend)
           write(*,*) 'dimx',dimx1,dimx2,dimx5,dimx6
           if(dtyp.eq.bmf_integer4) then
             do k=1,nk
             do j=1,nj
                if((j+ny0-1.lt.hole_j0).or.(j+ny0-1.gt.hole_j0+hole_jsize-1)) then
                   do i=1,ni
                      indice=indice+1
                      i4(i0+i-1,j0+j-1,k)=champ_courant%bmf_champ%tableau(indice,1,1)
                   enddo
                else if(((j+ny0-1.ge.hole_j0).or.(j+ny0-1.le.hole_j0+hole_jsize-1))) then
                   if ((dimx2-dimx1).gt.0) then
                      do i=dimx1,dimx2-1
                         indice=indice+1
                         i4(i0+i-nx0,j0+j-1,k)=champ_courant%bmf_champ%tableau(indice,1,1)
                      enddo
                   endif
                   if((dimx6-dimx5).gt.0) then
                      do i=dimx5,dimx6
                         indice=indice+1
                         i4(i0+i-nx0,j0+j-1,k)=champ_courant%bmf_champ%tableau(indice,1,1)
                      enddo
                   endif
                endif
             enddo
             enddo
           endif
!           write(*,*) 'retour'
      endif
      champ_courant=>champ_courant%champ_suivant
   enddo
 if(.not.trouve) then
    write(*,*) 'WARNING BMF_GET: Variable ',nom,' non trouvee'
    error=1
 endif
 return
 end
