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

	SUBROUTINE RPN_COMM_adj_halo(g,minx,maxx,miny,maxy, &
                  ni,nj,nk,halox,haloy,periodx,periody, &
                  gni,npol_row)
	use rpn_comm
	implicit none
!
!	exchange a halo with neighbours
!       or get halo contributions from neighbors 

	integer minx,maxx,miny,maxy,ni,nj,nk,halox,haloy
	integer gni,npol_row
	logical periodx,periody
	integer g(minx:maxx,miny:maxy,nk)
	real g_adj(minx:maxx,miny:maxy,nk)

!       has to be real because of g_adj = g_adj + wk\
!       instead of g_adj = wk (mathematical operation)

!	include 'rpn_comm.h'
!	include 'mpif.h'
!
	integer, target, allocatable, dimension(:) :: wk
	integer wkxs,wkxr,wkys,wkyr,wks,wkr,wkslag,wkgth
	pointer(g_adj_,g_adj)
	pointer(wkr_,wkr(10))
	pointer(wks_,wks(10))
	pointer(wkxs_,wkxs(1-halox:ni+halox,nk,haloy))
	pointer(wkxr_,wkxr(1-halox:ni+halox,nk,haloy))
	pointer(wkys_,wkys(1-haloy:nj+haloy,nk,halox))
	pointer(wkyr_,wkyr(1-haloy:nj+haloy,nk,halox))
	pointer(wkslag_,wkslag(miny:maxy,nk,10))
	pointer(wkgth_,wkgth(minx:maxx,nk,10))	

!       Second set of pointers, because of the real type
	real, target, allocatable, dimension(:) :: wk_adj
	real wkxs_adj,wkxr_adj,wkys_adj,wkyr_adj,wks_adj, &
        wkr_adj,wkslag_adj,wkgth_adj,wksband_adj,wkrband_adj
	
	pointer(wkr_adj_,wkr_adj(10))
	pointer(wks_adj_,wks_adj(10))
	pointer(wkxs_adj_,wkxs_adj(1-halox:ni+halox,nk,haloy))
	pointer(wkxr_adj_,wkxr_adj(1-halox:ni+halox,nk,haloy))
	pointer(wkys_adj_,wkys_adj(1-haloy:nj+haloy,nk,halox))
	pointer(wkyr_adj_,wkyr_adj(1-haloy:nj+haloy,nk,halox))
	pointer(wkslag_adj_,wkslag_adj(miny:maxy,nk,10))
	pointer(wkgth_adj_,wkgth_adj(miny:maxy,nk,10))
	pointer(wksband_adj_,wksband_adj(nk,haloy,10))
	pointer(wkrband_adj_,wkrband_adj(nk,haloy,10))
!
!

!
	integer nwdsx,nwdsy,nwds,halo,nwdslag,nwdsband
	integer i, j, k, m, proc, empl_proc
	integer sendtag, gettag, ierr
	integer status(MPI_STATUS_SIZE)
	logical east,west,north,south
	integer eastpe,westpe,northpe,southpe
	logical adjoint,polar,sendto_polar,recfrom_polar
!
	integer globalni,polarrows,nimax,nimin,ni_current
!
	integer land_fill
	real r_land_fill
	equivalence(land_fill,r_land_fill)
!
	globalni=abs(gni)
	polarrows=npol_row
	adjoint=.true.

1	continue
!	call RPN_COMM_tmg_in
!
	east=(bnd_east) .and. (.not.periodx)
	eastpe=pe_id(pe_mex+1,pe_mey)
	west=(bnd_west) .and. (.not.periodx)
	westpe=pe_id(pe_mex-1,pe_mey)
	north=(bnd_north) .and. (.not.periody)
	northpe=pe_id(pe_mex,pe_mey+1)
	south=(bnd_south) .and. (.not.periody)
	southpe=pe_id(pe_mex,pe_mey-1)
!
	nwdsx=(2*halox+ni)*nk*haloy+2
	nwdsx=nwdsx+mod(nwdsx,2)
	nwdsy=(2*haloy+nj)*nk*halox+2
	nwdsy=nwdsy+mod(nwdsy,2)
	nwds=max(nwdsx,nwdsy)
!
!	allocate temporary arrays on stack
!
	polar=.false.
	if (min(pe_mey+1,pe_ny-pe_mey).le.polarrows .and. pe_nx.gt.1) then
	  nwdslag=(maxy-miny+1)*nk*((globalni+pe_nx-1)/pe_nx)*(pe_nx+1)
	  nwdsband=haloy*nk*globalni
          allocate(wk(max(nwdslag,nwds*2)+2))
	  polar=.true.
	else
	  allocate(wk(nwds*2+2))
	endif

!       On place les vecteurs et on empile
	g_adj_=loc(g)
	wks_=loc(wk)
	wkslag_=loc(wk)
	wkr_=loc(wk(nwds+1))
	wkxr_=loc(wk(nwds+2))
	wkyr_=loc(wk(nwds+2))
	wkxs_=loc(wk(2))
	wkys_=loc(wk(2))
        wks_adj_=loc(wk)	
        wkslag_adj_=loc(wk)
        wkr_adj_=loc(wk(nwds+1))
        wkxr_adj_=loc(wk(nwds+2))
        wkyr_adj_=loc(wk(nwds+2))
        wkxs_adj_=loc(wk(2))
        wkys_adj_=loc(wk(2))
	wksband_adj_=loc(wk)
!	wkrband_adj_=loc(wk(nk*haloy*gni))
	wkrband_adj_=loc(wk)

!
!       Si halo nul, on passe. 
!       
! 
	halo=haloy
!
!       On doit connaitre nimax pour pouvoir indicer correctement
!       si dans le cas semi-lag
	if (polar) then
	   nimax = (globalni + pe_nx - 1)/pe_nx
	   nimin = globalni-nimax*(pe_nx-1)	
	endif
!
!  dest_polar true si pe_me au sud d'une bande polaire
!
!

	sendto_polar=(min(pe_mey,pe_ny-pe_mey+1).le.polarrows&
                  .and.pe_nx.gt.1)&
                  .and..not.south.and..not.periody
	recfrom_polar=(min(pe_mey+2,pe_ny-pe_mey-1).le.polarrows&
                  .and.pe_nx.gt.1)&
                  .and..not.north.and..not.periody

	if(polar) then
 
!       Effets de la periodicite est-ouest avant echange n-s
!       dans le cas semi-lag
!

	if(periodx)then
	    do m=1,halo
	    do k=1,nk
	    do j=1-haloy,nj+haloy
	       g_adj(m,j,k)=g_adj(m,j,k)+g_adj(globalni+m,j,k)
	       g_adj(globalni+m-halo,j,k)=g_adj(globalni+m-halo,j,k)&
                                   +g_adj(m-halo,j,k)
	    enddo
	    enddo
	    enddo
	endif

!
!       On additionne tout de suite les contributions de halos
!	en halos afin de ne transférer que des bandes de longueur ni
!
	   do m=1,halo
	      do k=1,nk
!VDIR NODEP
		 do i=1,gni
		    wksband_adj(k,m,i)=g_adj(i,m-halo,k)
		 enddo
	      enddo
	   enddo
!
!       Gather_all a la rangee
	if(pe_nx.gt.1) then	
	  call MPI_ALLGATHER(wksband_adj,nwdsband,MPI_REAL, &
              wkrband_adj,nwdsband,MPI_REAL, &
              pe_myrow,ierr)
	  do proc=0,pe_nx-1
	   if(proc.ne.pe_mex) then
	   do m=1,halo
	   do k=1,nk
!VDIR NODEP
	   do i=1,ni
	      g_adj(pe_mex*nimax+i,m-halo,k)=&
                   g_adj(pe_mex*nimax+i,m-halo,k)&
                  +wkrband_adj(k,m,i+proc*globalni+pe_mex*nimax)
	   enddo
	   enddo
	   enddo
	   endif
	
	  enddo
	 endif
	endif
	if (halo .eq. 0) goto 1111
	if ( pe_ny.eq.1 ) then 
	  if(periody)then
	    do m=1,halo
	    do k=1,nk
	    do i=1-halox,ni+halox
	      g_adj(i,m,k)=g_adj(i,m,k)+g_adj(i,nj+m,k)
	      g_adj(i,nj+m-halo,k)=g_adj(i,nj+m-halo,k)+g_adj(i,m-halo,k)
	    enddo
	    enddo
	    enddo
	    endif
	  goto 1111
	endif

!
!
!  process north to south move
!	

!
!       Cas normal
!
	if (polar) then
	      do k=1,nk
	   do m=1,halo
!	      do k=1,nk
!VDIR NODEP
		 do i=1,ni
		    wkxs_adj(i,k,m)=g_adj(pe_mex*nimax+i,m-halo,k)
		 enddo
	      enddo
	   enddo
 11       format(30F6.1)
	else
	do m=1,halo
	do k=1,nk
!VDIR NODEP
	do i=1-halox,ni+halox
	  wkxs_adj(i,k,m)=g_adj(i,m-halo,k)
	enddo
	enddo
	enddo
	endif
	do i=1,nwdsx
	  wkr(i)=0
	enddo
	wks(1)=pe_medomm
	wks(nwdsx)=southpe
	sendtag=pe_medomm
	gettag=northpe
        if(north) then 
! send to south_neighbor
          if(.not.south)then
            call MPI_SSEND(wks_adj,nwdsx,MPI_REAL,southpe,&
      sendtag,PE_DEFCOMM,ierr)
          endif
          wkr(1)=northpe
          wkr(nwdsx)=pe_medomm
        else if(south) then 
! receive from north_neighbor
          call MPI_RECV(wkr_adj,nwdsx,MPI_REAL,northpe,&
      gettag,PE_DEFCOMM,status,ierr)
        else 
! send to south_neighbor and receive from north_neighbor
          call MPI_SENDRECV(wks_adj,nwdsx,MPI_REAL,southpe,sendtag,&
      wkr_adj,nwdsx,MPI_REAL,northpe,gettag,&
      PE_DEFCOMM,status,ierr)
        endif
	if(wkr(1).ne.northpe .or. wkr(nwdsx).ne.pe_medomm) then
	  print *,'PE=',pe_medomm,' N->S error,',wkr(1),wkr(nwdsx),&
     	    ' should be ',northpe,pe_medomm
	  stop
	endif
	if(.not.north)then
	if(recfrom_polar)then
	do k=1,nk
	do m=1,halo
!	do k=1,nk
	do i=1,ni
	  g_adj(pe_mex*nimax+i,nj-halo+m,k)=&
                     g_adj(pe_mex*nimax+i,nj-halo+m,k)&
                     +wkxr_adj(i,k,m)
        enddo
	enddo
	enddo
	else if(polar)then
	do m=1,halo
	do k=1,nk
	do i=1-halox,ni+halox
	  g_adj(pe_mex*nimax+i,nj-halo+m,k)=&
                    g_adj(pe_mex*nimax+i,nj-halo+m,k)&
                    +wkxr_adj(i,k,m)
	enddo
	enddo
	enddo
	   
	else
	do m=1,halo
	do k=1,nk
	do i=1-halox,ni+halox
	  g_adj(i,nj-halo+m,k)=g_adj(i,nj-halo+m,k)+wkxr_adj(i,k,m)
	enddo
	enddo
	enddo
	endif
	endif
!
!
!  process south to north move
!
	sendto_polar=(min(pe_mey+2,pe_ny-pe_mey-1).le.polarrows&
                  .and.pe_nx.gt.1)&
                  .and..not.north.and..not.periody
	recfrom_polar=(min(pe_mey,pe_ny-pe_mey+1).le.polarrows&
                  .and.pe_nx.gt.1)&
                  .and..not.south.and..not.periody

	if(polar) then

!
!       On additionne tout de suite les contributions de halos
!	en halos afin de ne transférer que des bandes de longueur ni
!
	   do m=1,halo
	      do k=1,nk
!VDIR NODEP
		 do i=1,globalni
		    wksband_adj(k,m,i)=g_adj(i,nj+m,k)
		 enddo
	      enddo
	   enddo
!
!       Gather_all a la rangee
!	
	  call MPI_ALLGATHER(wksband_adj,nwdsband,MPI_REAL,&
              wkrband_adj,nwdsband,MPI_REAL,&
              pe_myrow,ierr)
	  
	  do proc=0,pe_nx-1
	   if(proc.ne.pe_mex) then
	   do m=1,halo
	   do k=1,nk
!VDIR NODEP
	   do i=1,ni
	      g_adj(pe_mex*nimax+i,nj+m,k)=&
                   g_adj(pe_mex*nimax+i,nj+m,k)&
                  +wkrband_adj(k,m,i+proc*globalni+pe_mex*nimax)
	   enddo
	   enddo
	   enddo
	   endif
	  enddo

	endif
!
!       Echange normal
!
	halo=haloy
	if (polar) then
	   do m=1,halo
	      do k=1,nk
!VDIR NODEP
		 do i=1,ni
		    wkxs_adj(i,k,m)=g_adj(pe_mex*nimax+i,nj+m,k)
		 enddo
	      enddo
	   enddo
	else
	do m=1,halo
	do k=1,nk
!VDIR NODEP
	do i=1-halox,ni+halox
	  wkxs_adj(i,k,m)=g_adj(i,nj+m,k)
	enddo
	enddo
	enddo
	endif
!VDIR NODEP
	do i=1,nwdsx
	  wkr(i)=0
	enddo
	wks(1)=pe_medomm
	wks(nwdsx)=northpe
	sendtag=pe_medomm
	gettag=southpe
!
        if(south) then 
! send to north_neighbor
          if(.not.north)then
            call MPI_SSEND(wks_adj,nwdsx,MPI_REAL,northpe,&
      sendtag,PE_DEFCOMM,ierr)
          endif
          wkr(1)=southpe
          wkr(nwdsx)=pe_medomm
        else if(north) then 
! receive from south_neighbor 
          call MPI_RECV(wkr_adj,nwdsx,MPI_REAL,southpe,&
      gettag,PE_DEFCOMM,status,ierr)
        else 
! send to north_neighbor and receive from south_neighbor
          call MPI_SENDRECV(wks_adj,nwdsx,MPI_REAL,northpe,sendtag,&
      wkr_adj,nwdsx,MPI_REAL,southpe,gettag,&
      PE_DEFCOMM,status,ierr)
        endif
	if(wkr(1).ne.southpe .or. wkr(nwdsx).ne.pe_medomm) then
	  print *,'PE=',pe_medomm,' N->S error,',wkr(1),wkr(nwdsx),&
     	    ' should be ',southpe,pe_medomm
	  stop
	endif
	if(.not.south)then
	if(recfrom_polar)then
	do m=1,halo
	do k=1,nk
	do i=1,ni
	  g_adj(pe_mex*nimax+i,m,k)=&
                    g_adj(pe_mex*nimax+i,m,k)&
                    +wkxr_adj(i,k,m)
	enddo
	enddo
	enddo
	else if(polar)then
	do m=1,halo
	do k=1,nk
	do i=1-halox,ni+halox
	  g_adj(pe_mex*nimax+i,m,k)=&
                    g_adj(pe_mex*nimax+i,m,k)&
                    +wkxr_adj(i,k,m)
	enddo
	enddo
	enddo
	else
	do m=1,halo
	do k=1,nk
!VDIR NODEP
	do i=1-halox,ni+halox
	  g_adj(i,m,k)=g_adj(i,m,k)+wkxr_adj(i,k,m)
	enddo
	enddo
	enddo
	endif
	endif


1111	continue
!
	halo=halox
!       Si halo nul, on passe.
	if (halo .eq. 0) goto 2222
!
!	if within polarrows rows from a pole, no lateral exchange
!	but a lateral domain gather. Polarrows is zero for a normal
!	halo exchange. It may be nonzero for a semi-lag exchange.
!
	if (polar) then
	  goto 2222
	endif

!       Si un seul processeur et periodicite, on evite MPI
	if ( pe_nx.eq.1 ) then
	  if(periodx) then	  
	    do m=1,halo
	    do k=1,nk
	    do j=1-haloy,nj+haloy
	       g_adj(m,j,k)=g_adj(m,j,k)+g_adj(ni+m,j,k)
	       g_adj(ni+m-halo,j,k)=g_adj(ni+m-halo,j,k)+g_adj(m-halo,j,k)
	    enddo
	    enddo
	    enddo
	  endif
	  goto 2222
	endif
!
!
!  process west to east move
!
	do m=1,halo
	do k=1,nk
!VDIR NODEP
	do j=1-haloy,nj+haloy
	  wkys_adj(j,k,m)=g_adj(ni+m,j,k)
	enddo
	enddo
	enddo

!VDIR NODEP
	do i=1,nwdsy
	  wkr(i)=0
	enddo
	wks(1)=pe_medomm
	wks(nwdsy)=eastpe
	sendtag=pe_medomm
	gettag=westpe
!	
        if(west) then 
! 	send to east_neighbor
          if(.not.east)then
            call MPI_SSEND(wks_adj,nwdsy,MPI_REAL,eastpe,&
      sendtag,PE_DEFCOMM,ierr)
          endif
          wkr(1)=westpe
          wkr(nwdsy)=pe_medomm
        else if(east) then 
! 	receive from west_neighbor
          call MPI_RECV(wkr_adj,nwdsy,MPI_REAL,westpe,&
      gettag,PE_DEFCOMM,status,ierr)
        else 
! 	send to east_neighbor and receive from west_neighbor
          call MPI_SENDRECV(wks_adj,nwdsy,MPI_REAL,eastpe,sendtag,&
      wkr_adj,nwdsy,MPI_REAL,westpe,gettag,&
      PE_DEFCOMM,status,ierr)
        endif
	if(wkr(1).ne.westpe .or. wkr(nwdsy).ne.pe_medomm) then
	  print *,'PE=',pe_medomm,' W->E error,',wkr(1),wkr(nwdsy),&
     	    ' should be ',westpe,pe_medomm
	  stop
	endif
	if(.not.west)then

	do m=1,halo
	do k=1,nk
!VDIR NODEP
	do j=1-haloy,nj+haloy
	  g_adj(m,j,k)=g_adj(m,j,k)+wkyr_adj(j,k,m)
	enddo
	enddo
	enddo
	endif
!
!
!  process east to west move
!
	halo=halox
	do m=1,halo
	do k=1,nk
!VDIR NODEP
	do j=1-haloy,nj+haloy
	  wkys_adj(j,k,m)=g_adj(m-halo,j,k)
	enddo
	enddo
	enddo
!VDIR NODEP
	do i=1,nwdsy
	  wkr(i)=0
	enddo
	wks(1)=pe_medomm
	wks(nwdsy)=westpe
	sendtag=pe_medomm
	gettag=eastpe
!	
        if(east) then 
! 	send to west_neighbor
          if(.not.west) then
            call MPI_SSEND(wks_adj,nwdsy,MPI_REAL,westpe,&
      sendtag,PE_DEFCOMM,ierr)
          endif
          wkr(1)=eastpe
          wkr(nwdsy)=pe_medomm
        else if(west) then 
! 	receive from east_neighbor
          call MPI_RECV(wkr_adj,nwdsy,MPI_REAL,eastpe,&
      gettag,PE_DEFCOMM,status,ierr)
        else 
! 	send to west_neighbor and receive from east_neighbor
          call MPI_SENDRECV(wks_adj,nwdsy,MPI_REAL,westpe,sendtag,&
      wkr_adj,nwdsy,MPI_REAL,eastpe,gettag,&
      PE_DEFCOMM,status,ierr)
        endif 
	if(wkr(1).ne.eastpe .or. wkr(nwdsy).ne.pe_medomm) then
	  print *,'PE=',pe_medomm,' W->E error,',wkr(1),wkr(nwdsy),&
     	    ' should be ',westpe,pe_medomm
	  stop
	endif
	if(.not.east)then
	do m=1,halo
	do k=1,nk
!VDIR NODEP
	do j=1-haloy,nj+haloy
	  g_adj(ni-halo+m,j,k)=g_adj(ni-halo+m,j,k)+wkyr_adj(j,k,m)
	enddo
	enddo
	enddo
	endif
!
2222	continue
!
!	special case for Semi-Lag transport
!
	if (polar) then
!
	  
!
!       on emballe et deballe tour a tour chaque tuile
!
	     do proc=0,pe_nx-1
		if (proc.eq.pe_nx-1) then
		   ni_current=nimin
		else
		   ni_current=nimax
		endif
!
!       on emballe...
!
		do j=miny,maxy
		do k=1,nk
!VDIR NODEP
		do i=1,ni_current
		   wkslag_adj(j,k,i)=g_adj(proc*nimax+i,j,k)		 
		enddo
	        enddo
	        enddo
		wkgth_adj_=loc(wkslag_adj(miny,1,1+ni_current))
!
!	  gather from wkslag of all row PEs int wkgth
!
        	nwdslag=ni_current*(maxy-miny+1)*nk
		call MPI_GATHER(wkslag_adj,nwdslag,MPI_REAL,&
               wkgth_adj,nwdslag,MPI_REAL,proc,pe_myrow,ierr)

!
!       deballage sur root=proc
!
		if(pe_mex.eq.proc) then 
		do m=0,pe_nx-1           ! tous les morceaux
		if(m.ne.proc) then       ! mais pas celui qui y est deja
		do j=miny,maxy
		do k=1,nk
!VDIR NODEP
		do i=1,ni_current
		   g_adj(pe_mex*nimax+i,j,k)=g_adj(pe_mex*nimax+i,j,k)&
                  +wkgth_adj(j,k,i+m*ni_current)
		enddo
	        enddo
	        enddo
		endif
		enddo
		endif
	     enddo
!
!	  now the final touch, the EW periodicity
!         inutile si adjoint puisque les halos ne sont pas retournes

3333	  continue
!
	endif
4444	continue
	deallocate(wk)
!
!	if(polarrows.gt.0) then
!	  call MPI_BARRIER(pe_grid,ierr)
!	endif
!	call RPN_COMM_tmg_out
        
!
	return

	entry adj_halo(g,minx,maxx,miny,maxy, &
                  ni,nj,nk,halox,haloy,periodx,periody)
	globalni=ni
	polarrows=0
	adjoint=.true.
	goto 1

	end
