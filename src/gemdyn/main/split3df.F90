module data_3df
      use, intrinsic :: iso_fortran_env
      implicit none


      logical, save :: initialized = .false.
      logical, save :: debug = .false.   ! activates lots of debug prints
      integer, save :: v = 2             ! verbosity level, 1-5
      character (len=4), save :: nomvar  ! current variable name
      integer, save :: ni1g,nj1g,nk1,nka_m,nka_t,presstype  ! global dimensions, number of levels , INPUT grid
      real, save :: xlon1,xlat1,xlon2,xlat2  ! area of interest for OUTPUT grid
      real(kind=REAL64), save, dimension(:), pointer :: xp1,yp1,xu1,yv1, &  ! horizontal coordinate vectors (INPUT grid)
                        ana_am_8,ana_bm_8,ana_at_8,ana_bt_8      ! vertical coordinate vectors (INPUT grid)
      character (len=8), save ::  dynophy  ! dynamics/physics marker
      integer, save :: nvar, ntra          ! number of variables and tracers

      integer, save :: ni0g0, nj0g0, nilg0, njlg0  ! limits in INPUT grid without halo
      integer, save :: ni0g1, nj0g1, nilg1, njlg1  ! limits in INPUT grid with halo (usually wider area that above)

      real(kind=REAL64), save :: xlim_1,xlim_n,ylim_1,ylim_n  ! window of interest in lat' lon' space
      real(kind=REAL64), save :: delta_x, delta_y, epsi_x, epsi_y
      real(kind=REAL64), dimension(:), save, pointer :: grid_x, grid_y       ! OUTPUT grid X and Y axis in radians
      integer, dimension(:), save, pointer :: tile_x, tile_y      ! start of tiles in OUTPUT grid (OUTPUT index space)
      integer, dimension(:), save, pointer :: iobloc_x, iobloc_y  ! start of OUTPUT grid IO blocks  (OUTPUT index space)
      integer, dimension(:), save, pointer :: start_bloc_x, end_bloc_x  ! start and end of OUTPUT grid IO blocks
      integer, dimension(:), save, pointer :: start_bloc_y, end_bloc_y  ! (INPUT index space) (halos are accounted for)
      integer, save :: bloci, blocj  !  number of OUTPUT grid IO blocks along x and y
      integer, save :: worldrank,worldsize

contains

      logical function it_is_my_job(i,j,blkx,blky)
      implicit none
      integer, intent(in) :: i,j,blkx,blky
      integer :: globalindex, nperproc, owner

      globalindex = (j-1)*blkx + (i-1)
      nperproc = (blkx*blky+(worldsize-2)) / (worldsize-1)
      owner = 1 + globalindex/nperproc
      it_is_my_job = (owner == worldrank)
      end function it_is_my_job

end module data_3df





!
!     read a set of 3df files and redo the horizontal split with different horizontal limits
!
      subroutine split3df
      use data_3df
      use tdpack
      use, intrinsic :: iso_fortran_env
      implicit none
      include 'mpif.h'

      integer :: right, left, below, above
      character (len=256) in,out
      NAMELIST /split3df_in/ right, left, below, above, in, out, debug, v

      integer :: g_ni,g_nj,ptopo_npex,ptopo_npey,ptopo_nblocx,ptopo_nblocy,ptopo_sblocx,ptopo_sblocy
      character (len=8), dimension(32) :: ignore
      logical split, dump, radians
      NAMELIST /split3df_out/ g_ni,g_nj,ptopo_npex,ptopo_npey,ptopo_nblocx,ptopo_nblocy, &
                              xlim_1,xlim_n,ylim_1,ylim_n, ignore, radians
      integer :: ierror
      real(kind=REAL64) :: tempx, tempy
      real(kind=REAL64) , parameter :: ONE = 1.0
      real(kind=REAL64) , parameter :: ZERO = 0.0
      integer :: i, j,  tilei, tilej, tempi, tempj

      call MPI_init(ierror)
      call MPI_comm_size(MPI_COMM_WORLD,worldsize,ierror)
      call MPI_comm_rank(MPI_COMM_WORLD,worldrank,ierror)

!
!     default values for namelists to be read
!
      g_ni = 2
      g_nj = 2
      ptopo_npex = 1     ! one tile
      ptopo_npey = 1
      ptopo_nblocx = 1   ! one file for entire domain
      ptopo_nblocy = 1
      xlim_1 = 0.0       ! lon = 0 -> 360
      xlim_n = 2.0*pi_8
      ylim_1 = -.5*pi_8    ! lat = -90 -> 90
      ylim_n = .5*pi_8
      left = 0           ! no halos
      right = 0
      below = 0
      above = 0
      radians = .true.   ! limits are in radians, not in degrees

      open(5,file='split3df.nml',access='SEQUENTIAL',form='FORMATTED',status='OLD')
      read(5,NML=split3df_in)
      rewind(5)                  ! in case namelists are not in expected order
      read(5,NML=split3df_out)
      close(5)

      if(.not. radians) then     ! convert limits to radians if supplied in degrees
        xlim_1=xlim_1*(pi_8/180.0)
        xlim_n=xlim_n*(pi_8/180.0)
        ylim_1=ylim_1*(pi_8/180.0)
        ylim_n=ylim_n*(pi_8/180.0)
      endif
      split = out/='' .and. worldsize>1   ! if only one process, no split will be done
      dump = worldsize==1
      ptopo_sblocx = (ptopo_npex+ptopo_nblocx-1)/ptopo_nblocx  ! size of IO blocks along x
      ptopo_sblocy = (ptopo_npey+ptopo_nblocy-1)/ptopo_nblocy  ! size of IO blocks along y
!==============================================================================
      if(worldrank==0) print *,'init_3df'
      call read_3df(in,right,left,below,above,split,dump)  ! get metadata from INPUT grid
      if(debug.and.worldrank==0)print *,'DEBUG: INPUT limits=',xp1(1),xp1(ni1g),yp1(1),yp1(nj1G)
!==============================================================================
!      xlim_1=min(2.0*pi, max(ZERO,xlim_1))   ! left side of OUTPUT grid, lon', 0-2*pi
!      xlim_n=min(2.0*pi, max(ZERO,xlim_n))   ! right side of OUTPUT grid, lon', 0-2*pi
!      ylim_1=min(.5*pi, max(-.5*pi,ylim_1))  ! bottom of OUTPUT grid, lat',-pi/2 to +pi/2
!      ylim_n=min(.5*pi, max(-.5*pi,ylim_n))  ! top of OUTPUT grid, lat', -pi/2 to +pi/2
      delta_x = (xlim_n-xlim_1)/(g_ni-1)     ! x (lon') spacing of OUTPUT grid (radians)
      delta_y = (ylim_n-ylim_1)/(g_nj-1)     ! y (lat') spacing of OUTPUT grid (radians)
      epsi_x = delta_x * .1                  ! "epsilon" used for comparisons, .1 grid point
      epsi_y = delta_y * .1
      if(debug.and.worldrank==0) then
        print *,'DEBUG: OUTPUT limits=',xlim_1,xlim_n,ylim_1,ylim_n
        print *,'DEBUG: OUTPUT deltas+epsis=',delta_x,delta_y,epsi_x,epsi_y
        print *,'DEBUG: OUTPUT grid size=',g_ni,g_nj
      endif
!==============================================================================
      allocate(grid_x(g_ni))  ! x coordinates of points in OUTPUT grid (longitude', radians)
      do i=1,g_ni-1
        grid_x(i) = xlim_1 + (i-1)*delta_x
      enddo
      grid_x(g_ni) = xlim_n
      allocate(grid_y(g_nj))  ! y coordinates of points in OUTPUT grid (latitude', radians)
      do j=1,g_nj-1
        grid_y(j) = ylim_1 + (j-1)*delta_y
      enddo
      grid_y(g_nj) = ylim_n

      allocate(tile_x(ptopo_npex+1), tile_y(ptopo_npey+1))
      tilei = (g_ni+ptopo_npex-1)/ptopo_npex   !  max size of a tile along x
      do i=1,ptopo_npex
        tile_x(i) = 1 + (i-1)*tilei            !  start of tile(i,j) along x
      enddo
      tile_x(ptopo_npex+1) = g_ni+1
      if(tile_x(ptopo_npex) > g_ni) then       ! .... OUCH !!!
        print *,'ERROR: tiling problem along X, aborting'
        goto 999
      endif
      tilej = (g_nj+ptopo_npey-1)/ptopo_npey   !  max size of a tile along y
      do j=1,ptopo_npey
        tile_y(j) = 1 + (j-1)*tilej            !  start of tile(i,j) along y
      enddo
      tile_y(ptopo_npey+1) = g_nj+1
      if(tile_y(ptopo_npey) > g_nj) then       ! .... OUCH !!!
        print *,'ERROR: tiling problem along Y, aborting'
        goto 999
      endif
if(debug.and.worldrank==0)print *,'DEBUG: tile_x=',tile_x
if(debug.and.worldrank==0)print *,'DEBUG: tile_y=',tile_y
!     tile(i,j) in OUTPUT grid covers
!     OUTPUT_grid( tile_x(i):tile_x(i+1)-1 , tile_y(j):tile_y(j+1)-1)
!==============================================================================
!      bloci = (ptopo_npex+ptopo_nblocx-1)/ptopo_nblocx  ! number of IO blocks along x
!      blocj = (ptopo_npey+ptopo_nblocy-1)/ptopo_nblocy  ! number of IO blocks along y
      bloci = ptopo_nblocx  ! number of IO blocks along x
      blocj = ptopo_nblocy  ! number of IO blocks along y
      if(debug.and.worldrank==0)print *,'DEBUG: IO blocks=',bloci,' by',blocj
      allocate(iobloc_x(bloci+1), iobloc_y(blocj+1))
      do i=1,bloci
        iobloc_x(i) = tile_x(1 + (i-1)*ptopo_sblocx)      ! x index of start of IO block (i,j)
      enddo                                               ! in OUTPUT grid index space
      iobloc_x(bloci+1) = g_ni + 1
      do j=1,blocj
        iobloc_y(j) = tile_y(1 + (j-1)*ptopo_sblocy)      ! y index of start of IO block (i,j)
      enddo                                               ! in OUTPUT grid index space
      iobloc_y(blocj+1) = g_nj + 1
!     IO block (i,j) in OUTPUT grid covers
!     OUTPUT_grid(iobloc_x(i):iobloc_x(i+1)-1 , iobloc_y(j):iobloc_y(j+1)-1)
      if(debug.and.worldrank==0) then
        print *,'DEBUG: iobloc_x=',iobloc_x
        print *,'DEBUG: iobloc_y=',iobloc_y
        print *,'DEBUG: input grid is',ni1g,' by',nj1g
      endif
!==============================================================================
      ni0g0 = 1    ! apply halo and find the bounds of the useful domain in INPUT grid
      do while( xp1(ni0g0+1) < xlim_1-epsi_x )  ! adjust lower longitude to area of interest
        ni0g0 = ni0g0 + 1
        if(ni0g0 > ni1g-1) exit
      enddo
      ni0g1 = max(ni0g0-left,1)

      nilg0 = ni1g
      do while( xp1(nilg0-1) > xlim_n+epsi_x )  ! adjust upper longitude to area of interest
        nilg0 = nilg0 - 1
        if(nilg0 < 2) exit
      enddo
      nilg1 = min(nilg0+right,ni1g)

      nj0g0 = 1
      do while( yp1(nj0g0+1) < ylim_1-epsi_y)  ! adjust lower latitude to area of interest
        nj0g0 = nj0g0 + 1
        if(nj0g0 > nj1g-1) exit
      enddo
      nj0g1 = max(nj0g0-below,1)

      njlg0 = nj1g
      do while( yp1(njlg0-1) > ylim_n+epsi_y )  ! adjust upper latitude to area of interest
        njlg0 = njlg0 - 1
        if(njlg0 < 2) exit
      enddo
      njlg1 = min(njlg0+above,nj1g)
      if(debug.and.worldrank==0)then
        print *,'DEBUG: X limits without/with halo=',xp1(ni0g0),xp1(nilg0),xp1(ni0g1),xp1(nilg1)
        print *,'DEBUG: Y limits without/with halo=',yp1(nj0g0),yp1(njlg0),yp1(nj0g1),yp1(njlg1)
      endif
!     useful domain extent including halo is
!     xp1(ni0g1)--->xp1(nilg1) along X ,  yp1(nj0g1)--->yp1(njlg1) along Y
!==============================================================================
!     now we must find where the OUTPUT IO blocks start and end in the INPUT grid
!     start_bloc_x(i) is the X index in INPUT grid of the leftmost point of IO block(i,j)
!     send_bloc_x(i) is the X index in INPUT grid of the leftmost point of IO block(i,j)
!     start_bloc_y(j) is the Y index in INPUT grid of the bottommost point of IO block(i,j)
!     end_bloc_y(j) is the Y index in INPUT grid of the topmmost point of IO block(i,j)
!     halos are accounted for so there may be overlap between adjacent IO blocks
      allocate(start_bloc_x(bloci), end_bloc_x(bloci))
      do i=1,bloci
        tempx = grid_x(iobloc_x(i))  ! x (lon') coordinate of letftmost point in IO block
        tempi = nilg1                ! right edge of OUTPUT area in INPUT grid
        do while(tempx < xp1(tempi))
          tempi = tempi - 1          ! move left
          if(tempi < 1) exit
        enddo
        start_bloc_x(i) = max(1,tempi - left)  ! apply left halo, clamp at 1 (left)
      enddo
      do i=1,bloci
        tempx = grid_x(iobloc_x(i+1)-1)  ! x (lon') coordinate of rightmost point in IO block
        tempi = 1                        ! left edge of OUTPUT area in INPUT grid
        do while(tempx > xp1(tempi))
          tempi = tempi + 1          ! move right
          if(tempi > ni1g) exit
        enddo
        end_bloc_x(i) = min(nilg1,tempi + right)  ! apply right halo, clamp at nilg1 (right)
      enddo
      if(debug.and.worldrank==0)then
        print *,'DEBUG: start_bloc_x=',start_bloc_x
        print *,'DEBUG: end_bloc_x=',end_bloc_x
      endif
      allocate(start_bloc_y(blocj), end_bloc_y(blocj))
      do j=1,blocj
        tempy = grid_y(iobloc_y(j))  ! y (lat') coordinate of bottom point in IO block
        tempj = njlg1                ! top edge of OUTPUT area in INPUT grid
        do while(tempy < yp1(tempj))
          tempj = tempj - 1          ! move down
          if(tempj < 1) exit
        enddo
        start_bloc_y(j) = max(1,tempj - below)  ! apply bottom halo, clamp at 1 (bottom)
      enddo
      do j=1,blocj
        tempy = grid_y(iobloc_y(j+1)-1)  ! y (lat') coordinate of top point in IO block
        tempj = 1                        ! bottom edge of OUTPUT area in INPUT grid
        do while(tempy > yp1(tempj))
          tempj = tempj + 1          ! move up
          if(tempj > nj1g) exit
        enddo
        end_bloc_y(j) = min(njlg1,tempj + above)  ! apply top halo, clamp at njlg1 (top)
      enddo
      if(debug.and.worldrank==0)then
        print *,'DEBUG: start_bloc_y=',start_bloc_y
        print *,'DEBUG: end_bloc_y=',end_bloc_y
      endif
!==============================================================================
      if(worldrank==0) then          ! process 0 is the INPUT grid reader
        call read_3df(in,right,left,below,above,split,dump)
        print *,'Read_3df : DONE'
!==============================================================================
      else                           ! process >= 1 is/are OUTPUT grid writer(s)
        if(v>1)print 10,'Split 3DF : splitting ',trim(in),' ==> ',trim(out),' (',bloci,' by ',blocj,' pieces)'
10      format(5A,I3,A,I3,A)
        call split_3df(out,bloci,blocj,right,left,below,above,split)
        print *,'Split_3df : DONE, process',worldrank,' of',worldsize
      endif
!==============================================================================
999   call MPI_finalize(ierror)   ! THE END !
      stop
      end subroutine split3df
!
!==============================================================================
!==============================================================================
      logical function ne8(a, b, n) ! true if any element of a is not equal
      use, intrinsic :: iso_fortran_env
      implicit none               ! to the corresponding element of b
      integer :: n
      real(kind=REAL64), dimension(n) :: a, b

      integer i

      ne8 = .true.
      do i = 1 , n
        if (a(i) /= b(i)) return  ! a difference has been found, return true
      end do
      ne8 = .false.              ! a and b are identical
      end function ne8
!
!==============================================================================
!==============================================================================
      subroutine read_3df(fn,right,left,below,above,split,dump)  ! called twice, first call is initialization
      use data_3df                                               ! second call is read
      use, intrinsic :: iso_fortran_env
      implicit none
      include 'mpif.h'
      integer, intent(IN) :: right, left, below, above
      logical, intent(IN) :: split, dump
      character (len=256), intent(IN) :: fn   ! INPUT file name base
!
      logical :: ne8
      external :: ne8
      integer :: mapiun
      character (len=15) :: seq
      integer, parameter :: MAXINFILES=90
      integer, dimension(MAXINFILES) :: inf, ni1s, nj1s, i0s, j0s
      integer :: infile
      integer :: i0, j0, tni, tnj
      real(kind=REAL64) :: x0, xn, y0, yn
      integer :: i, k, totvar, nbits, tnbits, ierror
      integer ni1, nj1, nilg2, njlg2
      integer tni1,tnj1,tnk1,tnka_m,tnka_t,tpresstype
      real :: txlon1,txlat1,txlon2,txlat2
      real(kind=REAL64), dimension(:), pointer :: txp1,typ1,txu1,tyv1, &
                        tana_am_8,tana_bm_8,tana_at_8,tana_bt_8
      character (len=8) ::  tdynophy
      integer :: tnvar, tntra
      character (len=4) :: tnomvar
      real  , dimension (:,:), allocatable :: tr1
      integer, dimension(3) :: xmit
!
      mapiun=8
      do i = 1 , MAXINFILES
       inf(i) = 100 - i
      enddo
      infile = 0
      open(mapiun,file=trim(fn)//'filemap.txt',access='SEQUENTIAL',form='FORMATTED',STATUS='OLD')
!==============================================================================
!     read next line from map file, check if 3df file will be useful
!==============================================================================
1     read(mapiun,*,end=2)i0, j0, x0, xn, y0, yn, tni, tnj
      write(seq,100)i0,'-',j0
100   format(I7.7,a1,i7.7)
!
!     the check to skip this file if not needed is temporarily commented out
!
!     this 3df file covers   (i0 : i0+tni-1 , j0 : j0+tnj-1) in global space
!     if(initialized) then
!       if(i0+tni-1 < ni0g1 .or. i0 > ni1g1 .or. j0+tnj-1 < nj0g1 .o. j0 > nj1g1) then
!          goto 1 ! tile is outside of area of interest
!       endif
!     endif
      infile = infile + 1   ! this 3df file is useful, open it and get/check metadata
      open(inf(infile),file=trim(fn)//seq,access='SEQUENTIAL',form='UNFORMATTED',status='OLD')
!==============================================================================
!     Input dimensions
!     (get from file 1 and check that other files tell same story)
!==============================================================================
      if(infile==1) then
        read (inf(infile)) nomvar,ni1g,nj1g,nka_m,nka_t,presstype
        if(initialized.and.v>1) print*, 'Read3df : Input dimensions ',nomvar,ni1g,nj1g,nka_m,nka_t,presstype
      else  ! check that the story is the same for the other files
        read (inf(infile)) tnomvar,tni1,tnj1,tnka_m,tnka_t,tpresstype
        if(tnomvar/=nomvar .or. tni1/=ni1g .or. tnj1/=nj1g .or. tnka_m/=nka_m .or.&
           tnka_t/=nka_t .or. tpresstype/=presstype) then
          print *,'ERROR: bad input dimensions in file ',trim(fn)//seq
        endif
      endif
!==============================================================================
!     Rotation parameters
!     (get from file 1 and check that other files tell same story)
!==============================================================================
      if(infile==1) then
        read (inf(infile)) xlon1,xlat1,xlon2,xlat2
        if(initialized.and.v>1) print*,'Read3df : Rotation parameters',xlon1,xlat1,xlon2,xlat2
      else  ! check that the story is the same for the other files
        read (inf(infile)) txlon1,txlat1,txlon2,txlat2
        if(xlon1/=txlon1 .or. xlat1/=txlat1 .or. xlon2/=txlon2 .or. xlat2/=txlat2)then
          print *,'ERROR: bad rotation parameters in file ',trim(fn)//seq
        endif
      endif
!==============================================================================
!      Horizontal positional parameters
!     (get from file 1 and check that other files tell same story)
!     array allocation done when reading first file
!==============================================================================
      if(infile==1) then
        allocate (xp1(ni1g),yp1(nj1g),xu1(ni1g),yv1(nj1g))
        allocate (txp1(ni1g),typ1(nj1g),txu1(ni1g),tyv1(nj1g))
        read (inf(infile)) xp1,yp1,xu1,yv1
      else  ! check that the story is the same for the other files
        read (inf(infile)) txp1,typ1,txu1,tyv1
        if(ne8(xp1,txp1,ni1g) .or. ne8(yp1,typ1,nj1g) .or. ne8(xu1,txu1,ni1g) .or. ne8(yv1,tyv1,nj1g)) then
          print *,'ERROR: bad horizontal positional parameters in file ',trim(fn)//seq
        endif
      endif
!==============================================================================
!     Vertical description
!     (get from file 1 and check that other files tell same story)
!     array allocation done when reading first file
!==============================================================================
      if(infile==1) then
        allocate(ana_am_8(nka_m), ana_bm_8(nka_m), &
                 ana_at_8(nka_t), ana_bt_8(nka_t) )
        allocate(tana_am_8(nka_m), tana_bm_8(nka_m), &
                 tana_at_8(nka_t), tana_bt_8(nka_t) )
        read (inf(infile)) ana_am_8,ana_bm_8,ana_at_8,ana_bt_8
      else  ! check that the story is the same for the other files
        read (inf(infile)) tana_am_8,tana_bm_8,tana_at_8,tana_bt_8
        if(ne8(ana_am_8,tana_am_8,nka_m) .or. ne8(ana_bm_8,tana_bm_8,nka_m) .or. &
           ne8(ana_at_8,tana_at_8,nka_t) .or. ne8(ana_bt_8,tana_bt_8,nka_t)) then
          print *,'ERROR: bad vertical description in file ',trim(fn)//seq
        endif
      endif
!==============================================================================
!     How many variables and tracers do we have ?
!     (get from file 1 and check that other files tell same story)
!==============================================================================
      if(infile==1) then
        read (inf(infile)) dynophy,nvar,ntra
        if(initialized.and.v>1) print*,'Read3df : Number of variables and tracers ',dynophy,nvar,ntra
      else  ! check that the story is the same for the other files
        read (inf(infile)) tdynophy,tnvar,tntra
        if(dynophy/=tdynophy .or. tnvar/=nvar .or. tntra/=ntra) then
          print *,'ERROR: bad number of variables and tracers in file ',trim(fn)//seq
        endif
      endif

      if(.not. initialized) then ! first time around, we only collect invariant metadata
        close(inf(infile))
        close(mapiun)
        initialized = .true.
        deallocate(txp1,typ1,txu1,tyv1,tana_am_8,tana_bm_8,tana_at_8,tana_bt_8)
        return
      endif

      i0s(infile) = i0  ! x coordinate of lower left point in tile
      j0s(infile) = j0  ! y coordinate of lower left point in tile
      if(infile==1.and.v>1) then
        print 180,' Read3df : output region of interest :',&
                xlim_1,' -->',xlim_n,' ,',ylim_1,' -->',ylim_n
        print 180,' Read3df : output region (with halo) :',&
                xp1(ni0g1),' -->',xp1(nilg1),' ,',yp1(nj0g1),' -->',yp1(njlg1)
        print 180,' Read3df : output region (no halo)   :',&
                 xp1(ni0g0),' -->',xp1(nilg0),' ,',yp1(nj0g0),' -->',yp1(njlg0)
        print 180,' Read3df : input domain limits       :',&
                 xp1(1),' -->',xp1(ni1g),' ,',yp1(1),' -->',yp1(nj1g)
180   format(4(A,F15.9))
      endif

      if(v>2)print 200,' Read3df : Opened ',trim(fn)//seq, x0,' ->',xn,'  ,',y0,' ->',yn, &
                       '    [',i0,' :',i0+tni-1,' ,',j0,' :',j0+tnj-1,']'
200   format(2A,4(F15.9,A),4(I4,A))
!==============================================================================
      goto 1         ! read next map file line
!==============================================================================
2     close(mapiun)  ! close map file
!==============================================================================
      if(infile==0)then
        print *,'ERROR: no useful 3df file found'
        return
      endif
!
!     metadata has been read, time to read data, reassemble it , and broadcast it to splitter(s)
!
      if(v>2)print *,'Read3df : READING ',nvar,' variable(s) +',ntra,' tracer(s)'
      if(v>2)print *,'Read3df : INPUT GRID SIZE =',ni1g,' by',nj1g

      allocate(tr1(ni1g,nj1g))  ! max size of input array
      do totvar = 1 , nvar+ntra
        do i=1,infile
!         get variable name, variable dimensions (of tile this time)
!         packing is ignored but checked
          if(i==1)then
!                         name,  x dim,  y dim, # levels, nbits
            read (inf(i)) nomvar,ni1s(i),nj1s(i),nk1     ,nbits
          else
            read (inf(i)) tnomvar,ni1s(i),nj1s(i),tnk1   ,tnbits
            if(tnomvar/=nomvar .or. tnk1/=nk1 .or. tnbits/=nbits) then  ! ni, nj may differ
              print *,'ERROR: inconsistent name/nk/nbits in file ',i
            endif
          endif
        enddo
        if(split)then  ! send nomvar and nk1 to splitter process (dest rank=1, tag=0)
        endif
        do k = 1 , nk1  ! read one level at a time, inner loop is on 3df files of interest
          tr1 = 1.0E+35 ! initialize to NONSENSE
          nilg2 = 1     ! some variables may be one column shorter (e.g. UU)
          njlg2 = 1     ! some variables may be one row shorter (e.g. VV)
          do i = 1 , infile
            ni1 = ni1s(i)    ! dimensions of tile
            nj1 = nj1s(i)
            i0 = i0s(i)      ! io,j0 lower left corner of tile
            j0 = j0s(i)
            if(i0>=1 .and. j0>=1 .and. i0+ni1-1<=ni1g .and. j0+nj1-1<=nj1g) then ! tile inside global space ?
              if(dump) then
                read(inf(i))tr1(i0:i0+ni1-1,j0:j0+nj1-1)
                if(k==1)print *,'DIAG: ',nomvar,i0,j0,i0+ni1-1,j0+nj1-1
                if(k==1)print *,'DIAG: ',nomvar,tr1(i0,j0),tr1(i0+ni1-1,j0),tr1(i0,j0+nj1-1),tr1(i0+ni1-1,j0+nj1-1)
              else
                read(inf(i))tr1(i0:i0+ni1-1,j0:j0+nj1-1)  ! read data into proper space in global INPUT array
               endif
               nilg2 = max(nilg2,i0+ni1-1)  ! find max column number read
               njlg2 = max(njlg2,j0+nj1-1)  ! find max row number read
            else
              print*,'ERROR: tile from file',i,' does not fit inside global space'
            endif
!            print 300, nomvar,ni1,nj1,nk1,nbits,totvar,'(',i0,'-',i0+ni1-1,',',j0,'-',j0+nj1-1,',',k,')'
300         format(A5,3I5,I3,I4,5(A,i4),A)
          enddo  ! i = 1 , infile
          nilg2 = min(nilg1,nilg2)  ! adjust to area of interest
          njlg2 = min(njlg1,njlg2)
          if(dump) then
            if(k==1)print 300, nomvar,ni1,nj1,nk1,nbits,totvar,'(',ni0g1,'-',nilg2,' ,',nj0g1,'-',njlg2,',',k,')'
            cycle ! k = 1 , nk1
          endif
          if(MAXVAL(tr1(ni0g1:nilg2,nj0g1:njlg2)) >= 1.0E+34) then  ! any hole(s) in data ? (missing 3df tiles)
            print *,'ERROR: hole in array',nomvar,' level=',k, ni0g1, nilg2,nj0g1, njlg2
            print *,tr1(ni0g1,nj0g1),tr1(ni0g1,njlg2),tr1(nilg2,nj0g1),tr1(nilg2,njlg2)
          else
            if(k==1.and.v>3)print 300, nomvar,ni1,nj1,nk1,nbits,totvar,'(',ni0g1,' :',nilg2,',',nj0g1,' :',njlg2,',',k,')'
          endif
!         full 2D array collected, send it to splitters (dest rank>=1, tag=k), including nonsense values
          if(split)then
            if(k==1)then
              call MPI_bcast(nomvar,4,MPI_CHARACTER,0,MPI_COMM_WORLD,ierror)  ! broadcast variable name
              xmit(1)=nk1
              xmit(2)=nilg2
              xmit(3)=njlg2
              call MPI_bcast(xmit,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)      ! broadcast dimensions
            endif
            call MPI_bcast(tr1,ni1g*nj1g,MPI_REAL,0,MPI_COMM_WORLD,ierror)    ! broadcast data
          endif
        enddo  ! k = 1 , nk1
      enddo  ! totvar = 1 , nvar+ntra
      end subroutine read_3df
!
!==============================================================================
!==============================================================================
      subroutine split_3df(fn,blkx,blky,right,left,below,above,split)  ! write OUTPUT files
      use data_3df
      use, intrinsic :: iso_fortran_env
      implicit none
      include 'mpif.h'
      integer, intent(IN) :: right, left, below, above
      character (len=256), intent(IN) :: fn  ! OUTPUT file name base
      logical, intent(IN) :: split
      integer, intent(IN) :: blkx, blky

      character(len=15) :: seq
      integer :: gridsize, ierror
      integer :: i0,in,j0,jn
      integer :: unf,k,mapiun,level
      integer i,j, nilg2, njlg2
      integer, dimension(:,:), pointer :: iuns
!      integer, dimension(:), pointer ::  limx, limy
      real  , dimension (:,:), allocatable :: tr1
      integer, dimension(3) :: xmit

      mapiun = 9
      unf = mapiun

      if(split)then  ! create output files, write constant information
        if(v>2)print 11,'Split 3DF : halo (right,left,below,above) = (',right,left,below,above,')'
11      format(A,4I3,A)
        if(worldrank==1) open(mapiun,file=trim(fn)//'filemap.txt',access='SEQUENTIAL',form='FORMATTED')
        allocate(iuns(blkx,blky))
        do j=1,blky
        do i=1,blkx
          if(.not. it_is_my_job(i,j,blkx,blky)) cycle
          unf = unf + 1
          iuns(i,j)=unf
          write(seq,100)start_bloc_x(i),'-',start_bloc_y(j)
100       format(I7.7,a1,i7.7)

          open(iuns(i,j),file=trim(fn)//seq,access='SEQUENTIAL',form='UNFORMATTED') ! create output file
          if(v>2)print *,'Split 3DF : Creating ',trim(fn)//seq

!          write (iuns(i,j)) nomvar, nilg1 - ni0g1 + 1, njlg1 - nj0g1 + 1, nka_m, nka_t, presstype
          write (iuns(i,j)) nomvar, ni1g, nj1g, nka_m, nka_t, presstype ! ,start_bloc_x(i),start_bloc_y(j)  ! new extra corner information
          write (iuns(i,j)) xlon1,xlat1,xlon2,xlat2
!          write (iuns(i,j)) xp1(ni0g1:nilg1),yp1(nj0g1:njlg1),xu1(ni0g1:nilg1),yv1(nj0g1:njlg1)
          write (iuns(i,j)) xp1,yp1,xu1,yv1   ! horizontal info
          write (iuns(i,j)) ana_am_8,ana_bm_8,ana_at_8,ana_bt_8      ! vertical info
          write (iuns(i,j)) dynophy,nvar,ntra     ! tag + number of variable(s) and tracer(s)
        enddo
        enddo
      endif

      gridsize = ni1g*nj1g
      allocate(tr1(ni1g,nj1g))    ! allocate full grid array

      if(v>3)print 170,' Split 3DF : (no halo)  ',nomvar,nilg0 - ni0g0 + 1,njlg0 - nj0g0 + 1,nka_m,nka_t,presstype
      if(v>3)print 180,' Split 3DF : (area)',xp1(ni0g0),xp1(nilg0),yp1(nj0g0),yp1(njlg0)
      if(v>3)print 170,' Split 3DF : (with halo)',nomvar,nilg1 - ni0g1 + 1,njlg1 - nj0g1 + 1,nka_m,nka_t,presstype
      if(v>3)print 180,' Split 3DF : (area)',xp1(ni0g1),xp1(nilg1),yp1(nj0g1),yp1(njlg1)
170   format(2A,5I8)
180   format(A,4F15.9)
!     How many variables and tracers do we have ?
      if(v>1)print*, 'Split 3DF : Number of variables and tracers ',dynophy,nvar,ntra

      do k=1,nvar+ntra
         if(split)then  ! and copy to all files, with proper dimensions (halo added)
!          get nomvar and nk1 from reader process (source rank=0, tag=0)
           call MPI_bcast(nomvar,4,MPI_CHARACTER,0,MPI_COMM_WORLD,ierror)  ! get variable name (root at process 0)
           call MPI_bcast(xmit,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)      ! get dimensions
           nk1=xmit(1)
           nilg2=xmit(2)
           njlg2=xmit(3)
           do level = 1 , nk1
!            Get data from reader process (root at process 0)
             call MPI_bcast(tr1,gridsize,MPI_REAL,0,MPI_COMM_WORLD,ierror)
             do j=1,blky
             do i=1,blkx
               i0 = start_bloc_x(i)
               in = end_bloc_x(i)
               j0 = start_bloc_y(j)
               jn = end_bloc_y(j)
               if(level==1 .and. k==1 .and. worldrank==1) write(mapiun,200)i0,j0,xp1(i0),xp1(in),yp1(j0),yp1(jn), &
                                                        in-i0+1,jn-j0+1,.true.
               if(.not. it_is_my_job(i,j,blkx,blky)) cycle
               in=min(in,nilg2)
               jn=min(jn,njlg2)
               if(level==1)write (iuns(i,j)) nomvar,in-i0+1,jn-j0+1,nk1,32  ! ,start_bloc_x(i),start_bloc_y(j)
               write (iuns(i,j)) tr1(i0:in,j0:jn)
               if(level==1.and.v>3)print 190,' Split 3DF : ',nomvar,'[',i0,' ,',j0,']   -> [',in,' ,',jn,']'
190            format(3A,4(I4,A))
200            format(2I8,4F15.9,2I8,L3)
             enddo
             enddo
           enddo  ! level = 1 , nk1
         endif
      enddo   !  k=1,nvar+ntra
      if(split)then
        do j=1,blky
        do i=1,blkx
          if(it_is_my_job(i,j,blkx,blky)) close(iuns(i,j))
        enddo
        enddo
        close(mapiun)
      endif
      return
      end subroutine split_3df
