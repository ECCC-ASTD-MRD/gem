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
!InTf!
      subroutine RPN_COMM_mydomain (call_back, mydomain)             !InTf!
      use rpn_comm
      implicit none                                                  !InTf!
!
      external :: call_back                                          !InTf!
      integer, intent(OUT) :: mydomain                               !InTf!
!
!      include 'mpif.h'

      logical mpi_started
      integer ierr, err, err_all, pe_tot2, pe_me2, npe_per_domain,&
     &        offset, ndomains
!
      if(RPM_COMM_IS_INITIALIZED) then
         if (pe_me.eq.0) write (6,1002)
         stop
      endif
      call MPI_INITIALIZED(mpi_started,ierr)
      if (.not. mpi_started ) call MPI_init(ierr)

      if( .not. WORLD_COMM_MPI_INIT ) then ! world not set before, use MPI_COMM_WORLD
           WORLD_COMM_MPI_INIT=.true.
           WORLD_COMM_MPI=MPI_COMM_WORLD
      endif

      call mpi_comm_size (WORLD_COMM_MPI,pe_tot2,ierr)
      call mpi_comm_rank (WORLD_COMM_MPI,pe_me2 ,ierr)

      call call_back (ndomains, offset, err)

      call mpi_allreduce &
     &   (err,err_all,1,MPI_INTEGER,MPI_MIN,WORLD_COMM_MPI,ierr)

      if (err_all.lt.0) then
         if (pe_me2.eq.0) write (6,1001) 
         call rpn_comm_FINALIZE(ierr)
         stop
      endif

      npe_per_domain = pe_tot2 / ndomains
      mydomain       = pe_me2  / npe_per_domain + offset
!
 1001 format &
     &(/'===> rpn_comm_mydomain: FATAL error in call_back --- ABORTING ---'/)
 1002 format &
     &(/'===> rpn_comm_mydomain: RPN_COMM already initialized --- ABORTING ---'/)
      return
      end subroutine RPN_COMM_mydomain                               !InTf!
!===================================================================
!InTf!
      subroutine RPN_COMM_world_get(world_comm) BIND(C,name='RPN_COMM_World_Get')   !InTf!
      use rpn_comm
      implicit none                                                !InTf!
      integer, intent(OUT) ::  world_comm                          !InTf!

      if( WORLD_COMM_MPI_INIT ) then  ! RPN_COMM_world_set has been called
        world_comm = WORLD_COMM_MPI
      else
        world_comm = MPI_COMM_WORLD
      endif
      return
      end subroutine RPN_COMM_world_get                            !InTf!
!===================================================================
!InTf!
      subroutine RPN_COMM_world_set(world_comm)                    !InTf!
      use rpn_comm
      implicit none                                                !InTf!
      integer, intent(IN) ::  world_comm                           !InTf!
!        world_comm=WORLD_COMM_MPI ! should rather be the other way around
      WORLD_COMM_MPI=world_comm
      WORLD_COMM_MPI_INIT=.true.
      return
      end subroutine RPN_COMM_world_set                            !InTf!
!===================================================================
!InTf!
      SUBROUTINE RPN_COMM_init(Userinit,Pelocal,Petotal,Pex,Pey)     !InTf!
      implicit none                                                  !InTf!
      integer, intent(out)   :: Pelocal,Petotal                      !InTf!
      integer, intent(inout) :: Pex,Pey                              !InTf!
      external Userinit                                              !InTf!
      integer junk, RPN_COMM_init_multigrid
      external RPN_COMM_init_multigrid
      junk = RPN_COMM_init_multigrid&
     &      (Userinit,Pelocal,Petotal,Pex,Pey,1)
      return
      end SUBROUTINE RPN_COMM_init                                   !InTf!
!===================================================================
!InTf!
!!INTEGER FUNCTION RPN_COMM_init_multigrid(Userinit,Pelocal,Petotal,Pex,Pey,MultiGrids) !InTf!
      INTEGER FUNCTION RPN_COMM_init_multigrid&
     &      (Userinit,Pelocal,Petotal,Pex,Pey,MultiGrids)
      use rpn_comm
      implicit none                                                  !InTf!
      external :: Userinit                                           !InTf!
      integer, intent(out)   :: Pelocal,Petotal                      !InTf!
      integer, intent(inout) :: Pex,Pey                              !InTf!
      integer, intent(in)    :: MultiGrids                           !InTf!
      integer rpn_comm_init_multi_level
      external rpn_comm_init_multi_level
      RPN_COMM_init_multigrid=RPN_COMM_init_multi_level&
     &      (Userinit,Pelocal,Petotal,Pex,Pey,MultiGrids,1)
      return
      end FUNCTION RPN_COMM_init_multigrid                           !InTf!
!===================================================================
!InTf!
!!INTEGER FUNCTION RPN_COMM_init_multi_level(Userinit,Pelocal,Petotal,Pex,Pey,MultiGrids,Grids)  !InTf!
      INTEGER FUNCTION RPN_COMM_init_multi_level&
     &      (Userinit,Pelocal,Petotal,Pex,Pey,MultiGrids,Grids)
      use rpn_comm
      implicit none                                                  !InTf!
      external :: Userinit                                           !InTf!
      integer, intent(out)   :: Pelocal,Petotal                      !InTf!
      integer, intent(inout) :: Pex,Pey                              !InTf!
      integer, intent(in)    :: MultiGrids                           !InTf!
      integer, intent(in)    :: Grids                                !InTf!
!arguments
!  I	Userinit	User routine that will be called by PE 0 to
!		get the processor grid topology if it is not supplied
!		(Pex .eq. 0 ) or (Pey .eq. 0)
!  O	Pelocal	PE rank (number) of local PE
!  O	Petotal	Number of PEs in job
!  I/O	Pex	Number of PEs along the X axis. If Pex=0 upon entry
!		it will be set to the proper value upon exit
!  I/O	Pey	Number of PEs along the Y axis. If Pey=0 upon entry
!		it will be set to the proper value upon exit
!  I    MultiGrids  number of simultaneous MultiGrids
!
!notes
!	processor topology common /pe/ will be filled here
!	positions are calculated from 0 (ORIGIN 0)
!!
!
!	include 'rpn_comm.h'
!	include 'mpif.h'
!
      integer ierr, i, j, count, npe, reste, nslots, key, status
      logical mpi_started
      integer gr1, gr2
      INTEGER newgroup,rowcomm
      integer, allocatable, dimension(:) :: proc_indomm
      integer unit, ndom, lndom, nproc, procmax,procmin
      type(domm), allocatable, dimension(:) :: locdom
      logical ok, allok
      logical RPN_COMM_grank
      integer RPN_COMM_petopo, RPN_COMM_version
      character *4096 SYSTEM_COMMAND,SYSTEM_COMMAND_2
      character *256 , dimension(:), allocatable :: directories
      character *256 :: my_directory, my_directory_file
      character *32 access_mode
      integer ncolors,my_color,directory_file_num,iun
      integer,dimension(:,:),allocatable::colors
      integer,dimension(:),allocatable::colortab
      integer version_marker, version_max, version_min
      integer pe_my_location(8)
      external RPN_COMM_unbind_process
      integer, external :: RPN_COMM_get_a_free_unit, RPN_COMM_set_timeout_alarm, fnom
!
!      if(RPM_COMM_IS_INITIALIZED) then ! ignore with warning message or abort ?
!      endif
      pe_indomm=-1
      pe_indomms=-1
      pe_dommtot=-1
      pe_a_domain=-1
      pe_me_a_domain=-1
      pe_all_domains=-1
      pe_me_all_domains=-1
      pe_multi_grid=-1
      pe_me_multi_grid=-1
      pe_grid=-1
      pe_me_grid=-1
      RPM_COMM_IS_INITIALIZED=.true.
      if( .not. WORLD_COMM_MPI_INIT ) then ! world not set before, use MPI_COMM_WORLD
           WORLD_COMM_MPI_INIT=.true.
           WORLD_COMM_MPI=MPI_COMM_WORLD
      endif
      RPN_COMM_init_multi_level = 0
      unit = 5
      ok = .true.

      call MPI_INITIALIZED(mpi_started,ierr)
      status = RPN_COMM_set_timeout_alarm(60)
      if (.not. mpi_started ) call MPI_init(ierr)
!       call MPI_barrier(WORLD_COMM_MPI,ierr)
      status = RPN_COMM_set_timeout_alarm(0)
      pe_wcomm=WORLD_COMM_MPI      ! UNIVERSE at this point

      call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)
      if(pe_me == 0)then
        call RPN_COMM_env_var("RPN_COMM_MONITOR",SYSTEM_COMMAND)
        if(SYSTEM_COMMAND .ne. " ") then
          iun = RPN_COMM_get_a_free_unit()
          open(iun,file=trim(SYSTEM_COMMAND),form='FORMATTED')
          write(iun,'(A)')'mpi_init successful'
          close(iun)
        endif
      endif
      call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)
      call RPN_COMM_unbind_process ! unbind processes if needed (FULL_UNBIND environment variable)

      version_marker=RPN_COMM_version()
!     check that all participants use the same version of rpn_comm
      call mpi_allreduce(version_marker, version_max, 1, MPI_INTEGER,&
     &                   MPI_BOR, WORLD_COMM_MPI, ierr)
!      write(rpn_u,*)"Version=",version_marker
      if(version_max .ne. version_marker)then
        write(rpn_u,*)&
     &       'FATAL ERROR: RPN_COMM version mismatch , STOPPING'
        call mpi_finalize(ierr)
        stop
      endif
!       initialize soft barrier setup in case we need a null task
!!      call rpn_comm_softbarrier_init_all
!       call resetenv   ! mpich1 will no longer work
!
      allocate(pe_location(8,0:pe_tot-1))

      pe_all_domains = pe_wcomm
      pe_me_all_domains = pe_me
      pe_tot_all_domains = pe_tot
      call MPI_COMM_GROUP(pe_all_domains,pe_gr_all_domains,ierr)
!
      allocate(colortab(0:pe_tot-1))
      my_color = 0
      SYSTEM_COMMAND=" "
      call RPN_COMM_env_var("RPN_COMM_DIAG",SYSTEM_COMMAND)
      if( SYSTEM_COMMAND .ne. " " ) then
          read(SYSTEM_COMMAND,*) diag_mode
      else
          diag_mode=1
      endif
!
!     if environment variable RPN_COMM_DOM is set, split into domains
!     this is used mainly by the r.mpirun family of scsripts
!
      SYSTEM_COMMAND=" "
      call RPN_COMM_env_var("RPN_COMM_DOM",SYSTEM_COMMAND)
      if( SYSTEM_COMMAND .ne. " " ) then
        SYSTEM_COMMAND_2=" "
        call RPN_COMM_env_var("RPN_COMM_DIRS",SYSTEM_COMMAND_2) ! get list of directories
        read(SYSTEM_COMMAND,*) ncolors  ! ABS(ncolors) is number of subdomains
        if( abs(ncolors) .gt. pe_tot ) then
           write(rpn_u,*)'ERROR: there are more subdomains (',&
     &                   abs(ncolors),') than PEs (',pe_tot,')'
           call RPN_COMM_finalize(ierr)
           stop
        endif
        if(ncolors .gt. 0) then
!         ncolors is followed by ncolors triplets ( first_pe, stride, last_pe )
          allocate(directories(0:ncolors))
          allocate(colors(3,ncolors))
          read(SYSTEM_COMMAND,*)ncolors,colors
          colortab=0        ! the size of colortab is the total number of PEs
          do i=1,ncolors    ! loop over subdomains
           do j=colors(1,i),colors(3,i),colors(2,i)
            colortab(j)=i   ! sweep the PEs for the current color(subdomain)
           enddo
          enddo
          my_color=colortab(pe_me)        ! my_color = 0 is an errror at this point
          if( my_color .eq. 0) then       ! it means that there are holes in the subdomain table
             ok=.false.
             write(rpn_u,*)'ERROR: pe ',pe_me,&
     &                     ' belongs to no subdomain'
          endif
          call MPI_ALLREDUCE(ok,allok,1,MPI_LOGICAL,&
     &                       MPI_LAND,WORLD_COMM_MPI,ierr)
          if( .not. allok ) then
             if(.not. ok) write(rpn_u,*)'ERROR DETECTED'
             call RPN_COMM_finalize(ierr)
             stop
          endif
          read(SYSTEM_COMMAND_2,*)directories
          my_directory=directories(my_color)

        else  ! (ncolors .gt. 0) ncolors < 0, all subdomains have the same size

          my_color=pe_me/(pe_tot/abs(ncolors))               !  mod(pe_tot, abs(ncolors)) must be zero
          if( mod(pe_tot,abs(ncolors)) .ne. 0 ) then
             ok=.false.
             write(rpn_u,*)'ERROR: subdomains have different sizes'
          endif
          call MPI_ALLREDUCE(ok,allok,1,MPI_LOGICAL,&
     &                       MPI_LAND,WORLD_COMM_MPI,ierr)
          if( .not. allok ) then
             if(.not. ok) write(rpn_u,*)'ERROR DETECTED'
             call RPN_COMM_finalize(ierr)
             stop
          endif
!         RPN_COMM_DOM contains -ncolors followed by the name of the file that contains
!         the directories to get into for the various subdomains (one line per subdomain)
          read(SYSTEM_COMMAND,*) ncolors,my_directory_file   ! get directory file
          directory_file_num=RPN_COMM_get_a_free_unit()
!          do i = 1 , 99  ! find an available unit number
!            inquire(UNIT=i,ACCESS=access_mode)
!            if(trim(access_mode) == 'UNDEFINED')then ! found
!              directory_file_num = i
!              exit
!            endif
!          enddo
!          ierr=fnom(directory_file_num,trim(my_directory_file),&
!     &              'SEQ+OLD+FTN',0)                        ! open file containing directory list
         open(UNIT=directory_file_num,file=trim(my_directory_file),&
     &        FORM='FORMATTED',STATUS='OLD')                 ! this file MUST exist
         do i=0,my_color-1
            read(directory_file_num,*)my_directory           ! skip unwanted entries
          enddo
          my_directory='.'                                   ! default value
          read(directory_file_num,*)my_directory             ! get my directory
          close(directory_file_num)
        endif  ! (ncolors .gt.0)

        call MPI_COMM_SPLIT(WORLD_COMM_MPI,my_color, &       ! split using my color
     &                     pe_me,pe_wcomm,ierr)
        RPN_COMM_init_multi_level = my_color
        call RPN_COMM_chdir(trim(my_directory))              ! cd to my directory
        if(diag_mode.ge.2)&
     &       write(rpn_u,*)"my directory is:",trim(my_directory)
        call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)               ! my rank
        call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)              ! size of my subdomain

      endif ! ( SYSTEM_COMMAND .ne. " " ) environment variable RPN_COMM_DOM is set
!     subdomains done (if needed)
      RPN_COMM_init_multi_level = my_color
      my_colors(1) = my_color    !  my domain number
!
      pe_a_domain = pe_wcomm
      pe_me_a_domain = pe_me
      pe_tot_a_domain = pe_tot

      if(pe_me_a_domain .eq.0 .and. diag_mode .ge. 1) &
     &   write(rpn_u,1000)my_color,pe_tot,MultiGrids,Grids
1000  format('domain=',I1,I4,' PEs,',I3,' SuperGrids with',&
     &       I3,' Grids each')
!      write(rpn_u, *)'pe_tot_a_domain=',pe_tot_a_domain
      call MPI_COMM_GROUP(pe_a_domain,pe_gr_a_domain,ierr)
!      write(rpn_u, *)'pe_gr_a_domain=',pe_gr_a_domain
!
!     verify that pe_tot_a_domain is a multiple of MultiGrids and Grids
!
      ok = .true.
      i=((pe_tot_a_domain/MultiGrids)/Grids)  ! number of PEs in a grid
      j=pe_tot_a_domain - i*MultiGrids*Grids  ! remainder
      if(j .ne. 0) then                       ! OOPS
        ok = .false.
        write(rpn_u,*)'ERROR: leftover PEs in domain'
      endif
!      write(rpn_u,*)'ok=',ok
      call mpi_allreduce(ok, allok, 1, MPI_LOGICAL,&
     &                   MPI_LAND,WORLD_COMM_MPI,ierr)
      if(.not.allok) then
           if(.not. ok) write(rpn_u,*)'ERROR DETECTED'
           call RPN_COMM_finalize(ierr)
           stop
      endif
!
!     if multiple grid program, (re)split domain into multigrids
!
      my_color = 0
!      pe_wcomms = pe_a_domain
!      write(rpn_u,*)'before multigrid split'
      if(MultiGrids .gt. 1) then
        my_color=pe_me/(pe_tot/MultiGrids)
        RPN_COMM_init_multi_level = my_color
        call MPI_COMM_SPLIT(pe_a_domain,my_color,pe_me,pe_wcomm,ierr)
        call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)               ! my rank
        call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)              ! size of my subdomain
      endif
      my_colors(2) = my_color   ! my multigrid number
!
      pe_multi_grid = pe_wcomm
      pe_me_multi_grid = pe_me
      pe_tot_multi_grid = pe_tot
!      write(rpn_u,*)'pe_tot_multi_grid=',pe_tot_multi_grid,'Grids=',Grids
      call MPI_COMM_GROUP(pe_multi_grid,pe_gr_multi_grid,ierr)
      pe_indomms = pe_multi_grid
      call MPI_COMM_GROUP(pe_indomms,pe_gr_indomms,ierr)
!
!     make communicator for PEs on same host in this multigrid
!       call RPN_COMM_split_by_node(pe_multi_grid,
!                                   pe_multigrid_host, pe_node_peers,
!                                   rank_on_host,      rank_in_peers, n_on_node, ierr)
!
!
!     domain split done, each domain is now on its own
!
!
!     if there are subgrids, (re)split multigrid into grids
!
!      write(rpn_u,*)'before grid split'
      my_color = 0
      if(Grids .gt. 1) then
        Pelocal = pe_me   ! me in my domain
        Petotal = pe_tot  ! number of pes in my domain
!        pe_wcomms = pe_multi_grid
        my_color=pe_me/(pe_tot/Grids)
        RPN_COMM_init_multi_level = my_color
        call MPI_COMM_SPLIT(pe_multi_grid,my_color,&
     &                      pe_me_multi_grid,pe_wcomm,ierr)
        call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)               ! my rank
        call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)              ! size of my subdomain
      endif
      my_colors(3) = my_color   ! my grid number
!
      pe_grid = pe_wcomm
      pe_me_grid = pe_me
      pe_tot_grid = pe_tot
!      write(rpn_u,*)'pe_tot_grid=',pe_tot_grid
      call MPI_COMM_GROUP(pe_grid,pe_gr_grid,ierr)
!
!     make communicator for PEs on same host in this grid
! TODO
!     utiliser RPN_COMM_split_by_node ou a tout le moins son algorithme
!
      my_color = abs(f_gethostid())  ! coloring by host identifier
      call MPI_COMM_SPLIT(pe_grid,my_color,&
     &                    pe_me_grid,pe_grid_host,ierr)        ! same (sub)grid, same host node communicator
      call MPI_COMM_RANK(pe_grid_host,pe_me_grid_host,ierr)    ! my rank on this host
      call MPI_COMM_SIZE(pe_grid_host,pe_tot_grid_host,ierr)   ! population of this host
      call MPI_COMM_GROUP(pe_grid_host,pe_gr_grid_host,ierr)   ! group communicator
!
!     subgrid split done, each sub domain is now on its own
!
      Pelocal = pe_me   ! me in my domain
      Petotal = pe_tot  ! number of pes in my grid
      call MPI_COMM_GROUP(pe_wcomm,pe_gr_wcomm,ierr)

      pe_indomm = pe_grid
      pe_gr_indomm = pe_gr_grid
      pe_pe0 = 0
      pe_dommtot = pe_tot
!
!     create peer to peer communicators between grids within a multigrid
!
      my_color = pe_me_grid
      call MPI_COMM_SPLIT(pe_multi_grid,my_color,&
     &                    pe_me_multi_grid,pe_grid_peers,ierr)
      call MPI_COMM_SIZE(pe_grid_peers,pe_tot_peer,ierr)  ! This must be equal to the number of grids
      call MPI_COMM_RANK(pe_grid_peers,pe_me_peer,ierr)   ! in a multigrid (or else...)
      call MPI_COMM_GROUP(pe_grid_peers,pe_gr_grid_peers,ierr)
!
!     Grid initialization, get PEs along X and Y
!
!      write(rpn_u,*)'READY to call UserInit'
	if(pe_me .eq. pe_pe0)then
	  if ( Pex.eq.0 .or. Pey.eq.0  ) then ! get processor topology
	    WORLD_pe(1)=pe_tot
	    WORLD_pe(2)=1
	    call Userinit(WORLD_pe(1),WORLD_pe(2))
	    if(WORLD_pe(1)*WORLD_pe(2).ne.pe_dommtot) then
	      ok = .false.
	      write(rpn_u,*) 'RPN_COMM_init: Inconsistency between'
	      write(rpn_u,*) 'userinit Subroutine and total number'
            write(rpn_u,*) 'of PE: please check'
	    endif
	    if(diag_mode.ge.1) then
            write(rpn_u,*)'Requested topology = ',&
     &                    WORLD_pe(1),' by ',WORLD_pe(2)
	      write(rpn_u,*)'Grid will use ',&
     &                    pe_dommtot,' processes'
          endif
        else ! ( Pex.eq.0 .or. Pey.eq.0  )
           write(rpn_u,*) 'RPN_COMM_init: Forced topology'
	     WORLD_pe(1) = Pex
	     WORLD_pe(2) = Pey
	     if(WORLD_pe(1)*WORLD_pe(2).ne.pe_dommtot) then
	       ok = .false.
	       write(rpn_u,*) 'RPN_COMM_init: Inconsistency between Pex'
	       write(rpn_u,*) 'and Pey args and total number of PE: '
	       write(rpn_u,*) 'please check'
	    endif
	    if(diag_mode.ge.1)&
     &       write(rpn_u,*)'Requested topology = ',WORLD_pe(1),' by '&
     &             ,WORLD_pe(2)
	  endif ! ( Pex.eq.0 .or. Pey.eq.0  )
!
	  if(WORLD_pe(1)*WORLD_pe(2) .gt. pe_dommtot) then
	    write(rpn_u,*)' ERROR: not enough PEs for decomposition '
	    write(rpn_u,*)' REQUESTED=',WORLD_pe(1)*WORLD_pe(2),&
     &                  ' AVAILABLE=',pe_dommtot
	    ok = .false.
	  endif

	endif  ! (pe_me .eq. pe_pe0)
!
!        write(rpn_u,*)'after UserInit'
      call mpi_allreduce(ok, allok, 1, MPI_LOGICAL,&
     &                   MPI_LAND,WORLD_COMM_MPI,ierr)
      if(.not.allok) then
           if(.not. ok) write(rpn_u,*)'ERROR DETECTED'
           call RPN_COMM_finalize(ierr)
           stop
      endif
!
!	send WORLD topology to all PEs. That will allow all PEs
!	to compute other PE topology parameters locally.
!       for doing this, we need to define some basic domains
!       communicators.

	call MPI_COMM_rank(pe_indomm,pe_medomm,ierr)
	pe_defcomm = pe_indomm
	pe_defgroup = pe_gr_indomm
!
!       broadcast number of PEs along X and Y axis
!       broadcast PE block sizes (deltai and deltaj)
!
        WORLD_pe(3) = deltai
        WORLD_pe(4) = deltaj
	call MPI_BCAST(WORLD_pe,size(WORLD_pe),MPI_INTEGER,0,&
     &	               pe_indomm,ierr)
	pe_nx  = WORLD_pe(1)
	pe_ny  = WORLD_pe(2)
	deltai = WORLD_pe(3)
	deltaj = WORLD_pe(4)
!
	if ( Pex.eq.0 .or. Pey.eq.0  ) then ! return processor topology
	  Pex = WORLD_pe(1)
	  Pey = WORLD_pe(2)
	endif
!
!	pe_pe0 is not equal to 0 if there are more than one domain
!	computational grid
!
	count = pe_pe0
!
!	fill tables containing the position along the X axis (pe_xtab)
!	and along the Y axis (pe_ytab) for all processors
!        write(rpn_u,*)'calling petopo'
	ierr = RPN_COMM_petopo(WORLD_pe(1),WORLD_pe(2))
!        write(rpn_u,*)'after petopo'

      pe_my_location(1) = pe_mex
      pe_my_location(2) = pe_mey
      pe_my_location(3) = pe_me_grid
      pe_my_location(4) = my_colors(3)
      pe_my_location(5) = pe_me_multi_grid
      pe_my_location(6) = my_colors(2)
      pe_my_location(7) = pe_me_a_domain
      pe_my_location(8) = my_colors(1)
      call MPI_allgather(&
     &     pe_my_location,8,MPI_INTEGER,&
     &     pe_location,   8,MPI_INTEGER,WORLD_COMM_MPI,&
     &     ierr)
      if( pe_me_all_domains .eq. 0 .and. diag_mode .ge.3) then
        write(rpn_u,*)'                         FULL PE MAP'
        write(rpn_u,*)&
     & '    mex     mey   me(g)    grid  me(sg)   sgrid   me(d)  domain'
        do j=0,pe_tot_all_domains-1
           write(rpn_u,1001)(pe_location(i,j),i=1,8)
1001       format(8I8)
        enddo
      endif

      BLOC_EXIST   =.false.
      BLOC_SIZEX   = 1
      BLOC_SIZEY   = 1
      BLOC_mybloc  = 0
      BLOC_myblocx = 0
      BLOC_myblocy = 0
      BLOC_me      = pe_me
      BLOC_comm_world = pe_indomm
      BLOC_comm_row = pe_myrow
      BLOC_comm_col = pe_mycol
      BLOC_corner = pe_pe0
      BLOC_master = 0
      if(pe_me.eq.pe_pe0) then
         BLOC_master=1
      endif
      pe_bloc = pe_indomm

      call MPI_Group_incl(pe_gr_indomm, 1, 0, pe_gr_blocmaster, ierr) 
      call MPI_Comm_create(pe_indomm,pe_gr_blocmaster, &
     &            pe_blocmaster, ierr)

!       for each communicator, store the corresponding group
!
      call MPI_COMM_GROUP(pe_myrow,pe_gr_myrow,ierr)
      call MPI_COMM_GROUP(pe_mycol,pe_gr_mycol,ierr)
      call MPI_COMM_GROUP(pe_bloc,pe_gr_bloc,ierr)
!      write(rpn_u,*)'peer communicator =',pe_grid_peers
!      write(rpn_u,*)'peer grid size =',pe_tot_peer
!      write(rpn_u,*)'peer rank =',pe_me_peer
!
!      ! what follows has to be added to rpn_comm or another module
!
!      integer, pointer, dimension(:,:) :: grid_id_table
!      integer, dimension(3) :: id_table
!      integer :: pe_pe0s, pe_me_pe0s  ! communicator for PE0 set and rank within
!
!      ! end of addition to rpn_comm or another module
!
!      my_color = min(1,pe_me_grid)  ! i am a PE 0 or not
!      call MPI_COMM_SPLIT(pe_all_domains,my_color,&
!     &                    pe_me_all_domains,pe_pe0s,ierr)
!      call MPI_COMM_rank(pe_pe0s,pe_me_pe0s,ierr) ! communicator nicknamed RPN_COMM_PE0='PE_00'
!      if(pe_me_grid /= 0) then  ! make sure it is invalid if not a PE0
!        pe_pe0s = -1
!        pe_me_pe0 = -1
!      endif
!
!      allocate(grid_id_table(3,pe_tot_all_domains))
!      id_table(1) = my_colors(1)*1024*1024 
!      id_table(1) = id_table(1) + my_colors(2)*1024
!      id_table(1) = id_table(1) + my_colors(3)   ! this will become a "grid id"
!      id_table(2) = pe_me_grid   ! local rank in grid
!      id_table(3) = pe_me_all_domains  ! global rank in "universe"
!      call MPI_allgather(id_table,3,MPI_INTEGER,&
!     &                    grid_id_table,3,MPI_INTEGER,&
!     &                    pe_all_domains,ierr)   ! progagate grid_id/local_rank/global_rank table
      return
!      contains

      end FUNCTION RPN_COMM_init_multi_level                        !InTf!

      integer function RPN_COMM_get_a_free_unit()                   !InTf!
      implicit none
      integer :: i
      character (len=16) :: access_mode
      RPN_COMM_get_a_free_unit=-1
      do i = 99,1,-1  ! find an available unit number
        inquire(UNIT=i,ACCESS=access_mode)
        if(trim(access_mode) == 'UNDEFINED')then ! found
          RPN_COMM_get_a_free_unit = i
          exit
        endif
      enddo
      return
      end function RPN_COMM_get_a_free_unit                         !InTf!

      function RPN_COMM_set_timeout_alarm(seconds) result(seconds_since)  !InTf!
      use ISO_C_BINDING
      implicit none
      integer, intent(IN) :: seconds  !InTf!
      integer :: seconds_since  !InTf!

      interface
      function c_alarm(seconds) result(seconds_since) BIND(C,name='alarm')
        use ISO_C_BINDING
        implicit none
        integer(C_INT), intent(IN), value :: seconds
        integer(C_INT) :: seconds_since
      end function c_alarm
      end interface

      seconds_since = c_alarm(seconds)
!      print *,'alarm set to ',seconds,' seconds'
      return
      end function RPN_COMM_set_timeout_alarm                             !InTf!
