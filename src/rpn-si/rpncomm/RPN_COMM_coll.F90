!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! !                          Environnement Canada
! !
! ! This library is free software; you can redistribute it and/or
! ! modify it under the terms of the GNU Lesser General Public
! ! License as published by the Free Software Foundation,
! ! version 2.1 of the License.
! !
! ! This library is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! ! Lesser General Public License for more details.
! !
! ! You should have received a copy of the GNU Lesser General Public
! ! License along with this library; if not, write to the
! ! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! ! Boston, MA 02111-1307, USA.
! !/

	subroutine RPN_COMM_coll(garr,gmini,gmaxi,gminj,&
               gmaxj,nig,njg,nk,ghalox,ghaloy,size,&
               larr,mini,maxi,minj,maxj,halox,haloy,&
               status)
	use rpn_comm
	implicit none
!
!	include 'rpn_comm.h'
!	include 'mpif.h'
!
	integer ghalox,ghaloy,gmini,gmaxi,gmaxj,gminj
	integer nig,njg,size,mini,maxi,minj,maxj,nk,status
	integer garr(size,gmini:gmaxi,gminj:gmaxj,nk),halox,haloy
	real reel,lreel
	integer larr(size,mini:maxi,minj:maxj,nk),i,j,k

	integer dimtemp(2),dt1,dt2,ierr,ixmin,ixmax,iymin,iymax,isize
	logical RPN_COMM_ngrank

	if(.not.RPN_COMM_ngrank(pe_defgroup)) return
        
        dimtemp(1)=maxi-mini+1
        dimtemp(2)=maxj-minj+1

        if(pe_tot.gt.1) then
         call RPN_COMM_bcast(dimtemp,2,"MPI_INTEGER",0,"GRID",ierr)
	else
	   ixmin=1-min(halox,ghalox)
	   iymin=1-min(haloy,ghaloy)
	   ixmax=nig+min(halox,ghalox)
	   iymax=njg+min(haloy,ghaloy)
	   do isize=1,size
	   do k=1,nk
	   do j=iymin,iymax
	   do i=ixmin,ixmax
	      garr(isize,i,j,k)=larr(isize,i,j,k)
	   enddo
	   enddo
	   enddo
	   enddo
	   return
        endif
        dt1=dimtemp(1)
        dt2=dimtemp(2)
        call RPN_COMM_coll2(garr,gmini,gmaxi,gminj,&
               gmaxj,nig,njg,nk,ghalox,ghaloy,size,&
               larr,mini,maxi,minj,maxj,halox,haloy,&
               status,dt1,dt2)
        return
        end


!**S/R RPN_COMM_coll  Global collection of data
!
!
!     WARNING: needs to be compiled with -O2 or -O3 on
!     SGI n32 mode with the standard MPI distribution of
!     SGI.
!
!****

	subroutine RPN_COMM_coll2(garr,gmini,gmaxi,gminj,&
               gmaxj,nig,njg,nk,ghalox,ghaloy,size,&
               larr,mini,maxi,minj,maxj,halox,haloy,&
               status,dt1,dt2)
	use rpn_comm
	implicit none
!
!	include 'rpn_comm.h'
!	include 'mpif.h'
!
	integer ghalox,ghaloy,gmini,gmaxi,gmaxj,gminj
	integer nig,njg,size,mini,maxi,minj,maxj,nk,status
	integer garr(size,gmini:gmaxi,gminj:gmaxj,nk),halox,haloy
	real reel,lreel
	integer larr(size,mini:maxi,minj:maxj,nk),dt1,dt2
	logical comphx,comphy,periodx,periody
!
!arguments
!  I	larr	local array containing data to send to PE=0
!  I	nig,njg,nk
!		dimensions of garr
!  I	size	size of data elements (the size of an integer is 1)
!  O	garr	global array that will receive the part of the local
!		domain that belongs to the PEs, used only by PE=0
!  I	mini:maxi,minj:maxj
!		dimensions of larr
!  I    gmini:gmaxi,gminj:gmaxj
!               dimensions of garr
!  I    ghalox,ghaloy,halox,haloy
!               size of halos of garr and larr
!  O	status	status code upon exit
!*
	integer i,j,k,ierr,nslots,nwords,i0,j0,level,nil,njl,islot
	integer clientpe,isz,j1,client,jlocal,islot0
	integer ipe, jpe, istatus(MPI_STATUS_SIZE),target,nidim,njdim
	integer, dimension(size,dt1,dt2,nk,pe_tot) :: temp
	integer, dimension(size,dt1,dt2,nk) :: temp2
	logical east,west,north,south,flag
	integer bxmin, bxmax, bymin,bymax,proc,njlind,nilind
	integer ixmin,ixmax,iymin,iymax
	integer lmini,lmaxi,lminj,lmaxj
	integer count(pe_nx+pe_ny)
	integer depl(pe_nx+pe_ny)
	integer mex, nx, mey, ny

	integer rpn_comm_limit

!
 1	status=MPI_ERROR
	
	if(pe_defcomm.eq.pe_indomm) then
	   mex= pe_mex
	   nx = pe_nx
	   mey= pe_mey
	   ny = pe_ny
	else if (pe_defcomm.eq.pe_bloc) then
	   
	   mex= pe_mex-BLOC_myblocx*(pe_nx/BLOC_sizex)
	   nx = pe_nx/BLOC_sizex
	   mey= pe_mey-BLOC_myblocy*(pe_ny/BLOC_sizey)
	   ny = pe_ny/BLOC_sizey
	else
	   if(pe_me.eq.pe_pe0) then
	      write(rpn_u,*) 'RPN_COMM_coll: unsupported for this communicator'
	      garr = -1
	      return	      
	   endif
	endif
 	ierr =  RPN_COMM_limit(mex, nx, 1, nig , lmini,&
     &	  lmaxi,count, depl)

	nil=lmaxi-lmini+1
	nidim=count(1)

	ierr =  RPN_COMM_limit(mey, ny, 1, njg , lminj,&
     &    lmaxj,count, depl)

	njl=lmaxj-lminj+1
	njdim=count(1)
          
         
         do isz=1,size
         do k=1,nk
         do j=1-haloy,njl+haloy
         do i=1-halox,nil+halox
            temp2(isz,i+halox,j+haloy,k)=larr(isz,i,j,k)
         enddo
         enddo
         enddo
         enddo
	 if(pe_tot.gt.1) then
          call MPI_GATHER(temp2,size*dt1*dt2*nk,&
              mpi_integer,temp,size*dt1*dt2*nk,&
              mpi_integer,0,pe_defcomm,ierr)
          if(ierr.ne.0) write(rpn_u,*) ierr,pe_medomm,'ERREUR RPN_COMM_COLL!'
	 endif

         if(mex.eq.0.and.mey.eq.0) then         
          do jpe=0,ny-1
            njlind=njl
            iymax=njlind
            iymin=1
            if(jpe.eq.0) iymin=1-min(haloy,ghaloy)
            if(jpe.eq.ny-1) then
               njlind=njg-(ny-1)*njl
               iymax=njlind+min(haloy,ghaloy)
            endif
         do ipe=0,nx-1       
            nilind=nil
            ixmax=nilind
            ixmin=1
            if(ipe.eq.0) ixmin=1-min(halox,ghalox)
            if(ipe.eq.nx-1) then
               nilind=nig-(nx-1)*nil
               ixmax=nilind+min(halox,ghalox)
            endif
         do isz=1,size
         do k=1,nk
         do j=iymin,iymax
         do i=ixmin,ixmax
           garr(isz,nil*ipe+i,njl*jpe+j,k)=&
           temp(isz,i+halox,j+haloy,k,1+pe_id(ipe,jpe))  ! use proper mapping table instead of hardwired linear mapping
!            temp(isz,i+halox,j+haloy,k,ipe+nx*jpe+1)    ! BUGGED code
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         endif

	status=MPI_SUCCESS
	return
!
1111	status =  MPI_ERROR
	write(rpn_u,*) 'ERREUR',MPI_ERROR,pe_medomm

	end
