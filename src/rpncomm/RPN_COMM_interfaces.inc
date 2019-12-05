      integer, external :: RPN_COMM_shuf_ezdist  !InTfX!   ! best interface we can provide for the time being
      integer, external :: RPN_COMM_shuf_ezcoll  !InTfX!   ! best interface we can provide for the time being


      interface

      integer function RPN_COMM_chdir(string)
      implicit none
      character (len=*) :: string
      end function RPN_COMM_chdir

      subroutine RPN_COMM_env_var(varname, value)
      implicit none
      character (len=*), intent(in) :: varname
      character (len=*), intent(out) :: value
      end subroutine RPN_COMM_env_var

      integer function RPN_COMM_hostid()
      implicit none
      end function RPN_COMM_hostid

      real*8 function RPN_COMM_wtime()
      implicit none
      end function RPN_COMM_wtime

      subroutine RPN_COMM_barrier(com, ierr)
      implicit none
      integer, intent(out) :: ierr
      character(len=*), intent(in) ::  com
      end subroutine RPN_COMM_barrier

      integer function RPN_COMM_barrier_data(level, values, nvalues)
      integer , intent(in) :: level
      integer , intent(in) :: nvalues
      real *8, dimension(nvalues), intent(out) :: values
      end function RPN_COMM_barrier_data

      subroutine RPN_COMM_bcastc(buffer, count, datatype, root, com, ierr)
      implicit none
      integer, intent(in) :: count, root
      integer, intent(out) :: ierr
      character(len=*), intent(inout) :: buffer
      character(len=*), intent(in) :: datatype, com
      end subroutine RPN_COMM_bcastc

      integer function RPN_COMM_bloc(nblocx, nblocy)
      implicit none
      integer, intent(in) :: nblocx, nblocy
      end function RPN_COMM_bloc

      integer function RPN_COMM_bloc_find(nblocx, nblocy, set)
      implicit none
      integer, intent(in) :: nblocx, nblocy
      logical, intent(in) :: set
      end function RPN_COMM_bloc_find

      integer function RPN_COMM_bloc_create(nblocx, nblocy)
      implicit none
      integer, intent(in) :: nblocx, nblocy
      end function RPN_COMM_bloc_create

      subroutine RPN_COMM_bloctopo(blocme, blocmex, blocmey, blocsizex, blocsizey)
      implicit none
      integer, intent(out) :: blocme, blocmex, blocmey
      integer, intent(out) :: blocsizex, blocsizey
      end subroutine RPN_COMM_bloctopo

      subroutine RPN_COMM_carac(npex, npey, me, medomm, mex, mey, sizex, sizey, &
      ismaster, mymaster, mybloc, myblocx, myblocy, blocme, domname)
      implicit none
      integer, intent(out) :: npex, npey, me, mex, mey, sizex, sizey, ismaster
      integer, intent(out) :: mymaster, mybloc, myblocx, myblocy, blocme
      integer, intent(out) :: medomm
      character(len=*), intent(out) :: domname
      end subroutine RPN_COMM_carac

      integer function RPN_COMM_colors(comm)
      implicit none
      character(len=*) :: comm
      end function RPN_COMM_colors

      integer function RPN_COMM_comm(com)
      implicit none
      character(len=*), intent(in) :: com
      end function RPN_COMM_comm

      integer function RPN_COMM_custom_comm(mpicom, name, mode)
      implicit none
      character(len=*), intent(in) :: name
      integer, intent(in) :: mpicom
      integer, intent(in) :: mode
      end function RPN_COMM_custom_comm

      subroutine RPN_COMM_i_comm(ctyp_c, ctyp, mcom)
         import :: rpncomm_communicator
      type(rpncomm_communicator), intent(out) :: ctyp
      character(len=*), intent(in) :: ctyp_c
      integer, optional, intent(in) :: mcom
      end subroutine RPN_COMM_i_comm

      subroutine RPN_COMM_const(IER)
      implicit none
      integer, intent(out) :: IER
      end   subroutine RPN_COMM_const

      integer function RPN_COMM_datyp_indx(data_int)
      implicit none
      character(len=*), intent(in) :: data_int
      end function RPN_COMM_datyp_indx

      integer function RPN_COMM_datyp(data_int)
      implicit none
      character(len=*), intent(in) :: data_int
      end function RPN_COMM_datyp

      subroutine RPN_COMM_i_datyp(dtyp_c, dtyp)
         import :: rpncomm_datatype
      type(rpncomm_datatype), intent(out) :: dtyp
      character(len=*), intent(in) :: dtyp_c
      end subroutine RPN_COMM_i_datyp

      subroutine RPN_COMM_i_user_datyp(dtyp_m, dtyp)
         import :: rpncomm_datatype
      type(rpncomm_datatype), intent(out) :: dtyp
      integer, intent(in) :: dtyp_m
      end subroutine RPN_COMM_i_user_datyp

      function RPN_COMM_i_valid_datyp(dtyp) result(valid)
         import :: rpncomm_datatype
      type(rpncomm_datatype), intent(in) :: dtyp
      logical :: valid
      end function RPN_COMM_i_valid_datyp

      function RPN_COMM_i_datyp_extent(dtyp) result(extent)
         import :: rpncomm_datatype
      type(rpncomm_datatype), intent(in) :: dtyp
      integer :: extent
      end function RPN_COMM_i_datyp_extent

      subroutine RPN_COMM_defo(com)
      implicit none
      character(len=*), intent(in) ::  com
      end subroutine RPN_COMM_defo

      subroutine RPN_COMM_fast_dist(garr, gmini, gmaxi, gminj, &
                     gmaxj, nig, njg, nk, ghalox, ghaloy, size, &
                     larr, mini, maxi, minj, maxj, halox, haloy, &
                     periodx, periody, status)
      implicit none
      integer, intent(in) :: ghalox, ghaloy, gmini, gmaxi, gmaxj, gminj
      integer, intent(in) :: nig, njg, size, mini, maxi, minj, maxj, nk
      integer, intent(out)::status
      integer, intent(in), target :: garr(size, gmini:gmaxi, gminj:gmaxj, &
         nk), halox, haloy
      integer, intent(out):: larr(size, mini:maxi, minj:maxj, nk)
      logical, intent(in) :: periodx, periody
      end subroutine RPN_COMM_fast_dist

      subroutine RPN_COMM_ezwin_create(sz, comm, mywin, ierr)
        import :: rpncomm_window
        integer, intent(in) :: sz, comm
        integer, intent(out) :: ierr
        type(rpncomm_window), intent(out) :: mywin
      end subroutine RPN_COMM_ezwin_create

      function RPN_COMM_ezwin_id(mywin) result(winid)
        import :: rpncomm_window
        type(rpncomm_window), intent(in) :: mywin
        integer :: winid
      end function RPN_COMM_ezwin_id

      function RPN_COMM_ezwin_ptr(mywin) result(winptr)
        import :: rpncomm_window, C_PTR
        type(rpncomm_window), intent(in) :: mywin
        type(C_PTR) :: winptr
      end function RPN_COMM_ezwin_ptr

      function RPN_COMM_ezwin_size(mywin) result(sz)
        import :: rpncomm_window
        type(rpncomm_window), intent(in) :: mywin
        integer :: sz
      end function RPN_COMM_ezwin_size

      subroutine RPN_COMM_ezwin_fetch_add(mywin, add, fetch, rank, offset, &
         ierr)
        import :: rpncomm_window
        type(rpncomm_window), intent(in) :: mywin
        integer, intent(in) :: add, rank, offset
        integer, intent(out) :: fetch, ierr
      end subroutine RPN_COMM_ezwin_fetch_add

      subroutine RPN_COMM_ezwin_get(mywin, pz, nw, rank, offset, ierr)
        import :: rpncomm_window, C_PTR
        type(rpncomm_window), intent(in) :: mywin
        type(C_PTR), intent(in), value :: pz
        integer, intent(in) :: nw, rank, offset
        integer, intent(out) :: ierr
      end subroutine RPN_COMM_ezwin_get

      subroutine RPN_COMM_ezwin_put(mywin, pz, nw, rank, offset, ierr)
        import :: rpncomm_window, C_PTR
        type(rpncomm_window), intent(in) :: mywin
        type(C_PTR), intent(in), value :: pz
        integer, intent(in) :: nw, rank, offset
        integer, intent(out) :: ierr
      end subroutine RPN_COMM_ezwin_put

      subroutine RPN_COMM_finalize(ierr)
      implicit none
      integer, intent(out) ::  ierr
      end subroutine RPN_COMM_finalize

      logical function RPN_COMM_grank(com)
      implicit none
      character(len=*), intent(in) :: com
      end function RPN_COMM_grank

      logical function RPN_COMM_ngrank(group)
      implicit none
      integer, intent(in) :: group
      end function RPN_COMM_ngrank

      function rpn_comm_create_2dgrid(gni, gnj, mini, maxi, minj, maxj) &
      result (grid_id)
        integer, intent(in) :: gni, gnj, mini, maxi, minj, maxj
        integer :: grid_id
      end function rpn_comm_create_2dgrid

      function rpn_comm_get_2dgrid(grid_id, dim_i, dim_j, gni, gnj, &
         mini, maxi, minj, maxj, starti, counti, startj, countj) result (status)
        integer, intent(in) :: grid_id, dim_i, dim_j
        integer, intent(out) :: gni, gnj, mini, maxi, minj, maxj
        integer, intent(out), dimension(dim_i) :: starti, counti
        integer, intent(out), dimension(dim_j) :: startj, countj
        integer :: status
      end function rpn_comm_get_2dgrid

      subroutine RPN_COMM_i_group(com, rgroup, rcom)
        import :: rpncomm_group
        import :: rpncomm_communicator
      implicit none
      character(len=*), intent(in) :: com
      type(rpncomm_group), intent(out) :: rgroup
      type(rpncomm_communicator), optional, intent(in) :: rcom
      end subroutine RPN_COMM_i_group

      integer function RPN_COMM_group(com)
      implicit none
      character(len=*), intent(in) :: com
      end function RPN_COMM_group

      subroutine RPN_COMM_haloflip(g, minx, maxx, miny, maxy, ni, nj, nk, halox, &
         haloy, ggni)
        implicit none
        integer, intent(in) :: minx, maxx, miny, maxy, ni, nj, nk, halox, haloy, &
         ggni
        integer, intent(inout) :: g(minx:maxx, miny:maxy, nk)
      end subroutine RPN_COMM_haloflip

      subroutine RPN_COMM_mydomain (call_back, mydomain)
      implicit none
      external :: call_back
      integer, intent(out) :: mydomain
      end subroutine RPN_COMM_mydomain

      subroutine RPN_COMM_world_set(world_comm)
      implicit none
      integer, intent(in) ::  world_comm
      end subroutine RPN_COMM_world_set

      subroutine RPN_COMM_init(Userinit, Pelocal, Petotal, Pex, Pey)
      implicit none
      integer, intent(out)   :: Pelocal, Petotal
      integer, intent(inout) :: Pex, Pey
      external Userinit
      end subroutine RPN_COMM_init

      integer FUNCTION RPN_COMM_init_multigrid(Userinit, Pelocal, Petotal, &
         Pex, Pey, MultiGrids)
      implicit none
      external :: Userinit
      integer, intent(out)   :: Pelocal, Petotal
      integer, intent(inout) :: Pex, Pey
      integer, intent(in)    :: MultiGrids
      end FUNCTION RPN_COMM_init_multigrid

      integer FUNCTION RPN_COMM_init_multi_level(Userinit, Pelocal, Petotal, &
         Pex, Pey, MultiGrids, Grids)
      implicit none
      external :: Userinit
      integer, intent(out)   :: Pelocal, Petotal
      integer, intent(inout) :: Pex, Pey
      integer, intent(in)    :: MultiGrids
      integer, intent(in)    :: Grids
      end FUNCTION RPN_COMM_init_multi_level

      integer function RPN_COMM_get_a_free_unit()
      end function RPN_COMM_get_a_free_unit
      function RPN_COMM_set_timeout_alarm(seconds) result(seconds_since)
      integer, intent(in) :: seconds
      integer :: seconds_since
      end function RPN_COMM_set_timeout_alarm

      subroutine RPN_COMM_make_io_pe_list(x, y, npes, pe_nx, pe_ny, method)
          integer, intent(in) :: npes
          integer, dimension(npes), intent(out) :: x
          integer, dimension(npes), intent(out) :: y
          integer, intent(in) :: pe_nx, pe_ny
          integer, intent(in) :: method
      end subroutine RPN_COMM_make_io_pe_list

      function RPN_COMM_check_ioset(newset, x , y, npes, pe_nx, pe_ny, &
         pe_me, diag) result(status)
         integer, intent(in) :: newset, npes
         integer, intent(in) :: pe_nx, pe_ny
         integer, intent(in) :: pe_me
         logical, intent(in) :: diag
         integer, intent(in), dimension(npes) :: x
         integer, intent(in), dimension(npes) :: y
         integer :: status
      end function RPN_COMM_check_ioset

      function RPN_COMM_create_io_set(npes, method) result(setno)
        implicit none
        integer, intent(in)  :: npes
        integer, intent(in)  :: method
        integer :: setno
      end function RPN_COMM_create_io_set

      function RPN_COMM_free_io_set(setno) result(freed)
        implicit none
        integer, intent(in) :: setno
        integer :: freed
      end function RPN_COMM_free_io_set

      function RPN_COMM_is_valid_io_setno(setno) result(valid)
        implicit none
        integer, intent(in) :: setno
        logical :: valid
      end function RPN_COMM_is_valid_io_setno

      function RPN_COMM_is_io_pe(setno) result(ordinal)
        implicit none
        integer, intent(in) :: setno
        integer :: ordinal
      end function RPN_COMM_is_io_pe

      function RPN_COMM_io_pe_gridid(setno, n) result(ordinal)
        implicit none
        integer, intent(in) :: setno
        integer, intent(in) :: n
        integer :: ordinal
      end function RPN_COMM_io_pe_gridid

      function RPN_COMM_io_pe_comm(setno) result(communicator)
        implicit none
        integer, intent(in) :: setno
        integer :: communicator
      end function RPN_COMM_io_pe_comm

      function RPN_COMM_io_pe_size(setno) result(population)
        implicit none
        integer, intent(in) :: setno
        integer :: population
      end function RPN_COMM_io_pe_size

      function RPN_COMM_io_pe_groups(setno) result(ngroups)
        implicit none
        integer, intent(in) :: setno
        integer :: ngroups
      end function RPN_COMM_io_pe_groups

      function RPN_COMM_io_pe_idlist(setno) result(idlist)
        implicit none
        integer, intent(in) :: setno
        integer, dimension(:), pointer :: idlist
      end function RPN_COMM_io_pe_idlist

      function RPN_COMM_io_pe_coord(setno) result(list)
        implicit none
        integer, intent(in) :: setno
        integer, dimension(:, :), pointer :: list
      end function RPN_COMM_io_pe_coord

      function RPN_COMM_io_pe_callback(setno, callback, argv) result(statu&
     &s)
        import :: C_PTR
        implicit none
        integer, intent(in) :: setno
        type(C_PTR) :: argv
        integer, external :: callback
        integer :: status
      end function RPN_COMM_io_pe_callback

      function RPN_COMM_io_pe_valid_set(x, y, npes, penx, peny, diag, method) &
         result(status)
          integer, intent(in) :: npes
          integer, dimension(npes), intent(out) :: x
          integer, dimension(npes), intent(out) :: y
          integer, intent(in) :: penx, peny
          integer, intent(in) :: method
          logical, intent(in) :: diag
          integer :: status
      end function RPN_COMM_io_pe_valid_set

      integer function RPN_COMM_limit_2(my_id, npe, gmin, gmax, lmini, lmaxi, &
         count, offset, relax)
      implicit none
      integer, intent(in) ::  my_id, npe, gmin, gmax, relax
      integer, intent(out) :: lmini, lmaxi
      integer, intent(out) :: count(npe), offset(npe)
      end function RPN_COMM_limit_2

      integer function RPN_COMM_limit(my_id, npe, gmin, gmax, lmini, lmaxi, &
         count, offset)
      implicit none
      integer, intent(in) ::  my_id, npe, gmin, gmax
      integer, intent(out) :: lmini, lmaxi
      integer, intent(out) :: count(npe), offset(npe)
      end function RPN_COMM_limit

      subroutine RPN_COMM_UP2LOW(S1, S2)
        implicit none
        character (len=*), intent(in) :: S1
        character (len=*), intent(out) :: S2
      end subroutine RPN_COMM_UP2LOW

      integer function RPN_COMM_mype(Me, Me_x, Me_y)
      implicit none
      integer, intent(out) ::  Me, Me_x, Me_y
      end  function RPN_COMM_mype

      integer function RPN_COMM_oper(op)
        implicit none
        character(len=*), intent(in) :: op
      end function RPN_COMM_oper

      subroutine RPN_COMM_i_oper(op, r_oper)
       import :: rpncomm_operator
        character(len=*), intent(in) :: op
        type(rpncomm_operator), intent(out) :: r_oper
      end subroutine RPN_COMM_i_oper

      function RPN_COMM_i_valid_oper(r_oper) result (is_valid)
       import :: rpncomm_operator
        type(rpncomm_operator), intent(in) :: r_oper
        logical :: is_valid
      end function RPN_COMM_i_valid_oper

      integer function RPN_COMM_option_L(name_S, value_L)
      implicit none
      character(len=*), intent(in) :: name_S
      logical, intent(in) :: value_L
      end function RPN_COMM_option_L

      integer function RPN_COMM_option(name_S, value)
      implicit none
      character(len=*), intent(in) :: name_S
      integer, intent(in) :: value
      end function RPN_COMM_option

      subroutine RPN_COMM_optn(op_type, op_ival, op_rval, op_cval)
      implicit none
      character(len=*), intent(in) :: op_type, op_cval
      integer, intent(in) :: op_ival
      real *4, intent(in) :: op_rval
      end subroutine RPN_COMM_optn

      subroutine RPN_COMM_set_petopo(sizx, sizy)
      implicit none
      integer, intent(in) :: sizx, sizy
      end subroutine RPN_COMM_set_petopo

      integer function RPN_COMM_petopo(pex, pey)
      implicit none
      integer, intent(in) :: pex, pey
      end function RPN_COMM_petopo

      integer function RPN_COMM_get_pe(x, y, grd, sgrd, communicator)
      implicit none
      integer, intent(in) :: x, y, grd, sgrd
      character *(*), intent(in) :: communicator
      end function RPN_COMM_get_pe

      end interface

      interface rpn_comm_ptr

        subroutine RPN_COMM_ptr_i4_1d(what, ptr)
        use ISO_C_BinDinG
        import :: rpncomm_ptr
        implicit none
        integer(KinD=4), dimension(1), intent(in), target :: what
        type(rpncomm_ptr), intent(out) :: ptr
        end subroutine RPN_COMM_ptr_i4_1d

        subroutine RPN_COMM_ptr_i4_2d(what, ptr)
        use ISO_C_BinDinG
        import :: rpncomm_ptr
        implicit none
        integer(KinD=4), dimension(1, 1), intent(in), target :: what
        type(rpncomm_ptr), intent(out) :: ptr
        end subroutine RPN_COMM_ptr_i4_2d

        subroutine RPN_COMM_ptr_i4_3d(what, ptr)
        use ISO_C_BinDinG
        import :: rpncomm_ptr
        implicit none
        integer(KinD=4), dimension(1, 1, 1), intent(in), target :: what
        type(rpncomm_ptr), intent(out) :: ptr
        end subroutine RPN_COMM_ptr_i4_3d

        subroutine RPN_COMM_ptr_i4_4d(what, ptr)
        use ISO_C_BinDinG
        import :: rpncomm_ptr
        implicit none
        integer(KinD=4), dimension(1, 1, 1, 1), intent(in), target :: what
        type(rpncomm_ptr), intent(out) :: ptr
        end subroutine RPN_COMM_ptr_i4_4d

        subroutine RPN_COMM_ptr_i8_1d(what, ptr)
        use ISO_C_BinDinG
        import :: rpncomm_ptr
        implicit none
        integer(KinD=8), dimension(1), intent(in), target :: what
        type(rpncomm_ptr), intent(out) :: ptr
        end subroutine RPN_COMM_ptr_i8_1d

        subroutine RPN_COMM_ptr_i8_2d(what, ptr)
        use ISO_C_BinDinG
        import :: rpncomm_ptr
        implicit none
        integer(KinD=8), dimension(1, 1), intent(in), target :: what
        type(rpncomm_ptr), intent(out) :: ptr
        end subroutine RPN_COMM_ptr_i8_2d

        subroutine RPN_COMM_ptr_i8_3d(what, ptr)
        use ISO_C_BinDinG
        import :: rpncomm_ptr
        implicit none
        integer(KinD=8), dimension(1, 1, 1), intent(in), target :: what
        type(rpncomm_ptr), intent(out) :: ptr
        end subroutine RPN_COMM_ptr_i8_3d

        subroutine RPN_COMM_ptr_i8_4d(what, ptr)
        use ISO_C_BinDinG
        import :: rpncomm_ptr
        implicit none
        integer(KinD=8), dimension(1, 1, 1, 1), intent(in), target :: what
        type(rpncomm_ptr), intent(out) :: ptr
        end subroutine RPN_COMM_ptr_i8_4d

        subroutine RPN_COMM_ptr_r4_1d(what, ptr)
        use ISO_C_BinDinG
        import :: rpncomm_ptr
        implicit none
        real(KinD=4), dimension(1), intent(in), target :: what
        type(rpncomm_ptr), intent(out) :: ptr
        end subroutine RPN_COMM_ptr_r4_1d

        subroutine RPN_COMM_ptr_r4_2d(what, ptr)
        use ISO_C_BinDinG
        import :: rpncomm_ptr
        implicit none
        real(KinD=4), dimension(1, 1), intent(in), target :: what
        type(rpncomm_ptr), intent(out) :: ptr
        end subroutine RPN_COMM_ptr_r4_2d

        subroutine RPN_COMM_ptr_r4_3d(what, ptr)
        use ISO_C_BinDinG
        import :: rpncomm_ptr
        implicit none
        real(KinD=4), dimension(1, 1, 1), intent(in), target :: what
        type(rpncomm_ptr), intent(out) :: ptr
        end subroutine RPN_COMM_ptr_r4_3d

        subroutine RPN_COMM_ptr_r4_4d(what, ptr)
        use ISO_C_BinDinG
        import :: rpncomm_ptr
        implicit none
        real(KinD=4), dimension(1, 1, 1, 1), intent(in), target :: what
        type(rpncomm_ptr), intent(out) :: ptr
        end subroutine RPN_COMM_ptr_r4_4d

        subroutine RPN_COMM_ptr_r8_1d(what, ptr)
        use ISO_C_BinDinG
        import :: rpncomm_ptr
        implicit none
        real(KinD=8), dimension(1), intent(in), target :: what
        type(rpncomm_ptr), intent(out) :: ptr
        end subroutine RPN_COMM_ptr_r8_1d

        subroutine RPN_COMM_ptr_r8_2d(what, ptr)
        use ISO_C_BinDinG
        import :: rpncomm_ptr
        implicit none
        real(KinD=8), dimension(1, 1), intent(in), target :: what
        type(rpncomm_ptr), intent(out) :: ptr
        end subroutine RPN_COMM_ptr_r8_2d

        subroutine RPN_COMM_ptr_r8_3d(what, ptr)
        use ISO_C_BinDinG
        import :: rpncomm_ptr
        implicit none
        real(KinD=8), dimension(1, 1, 1), intent(in), target :: what
        type(rpncomm_ptr), intent(out) :: ptr
        end subroutine RPN_COMM_ptr_r8_3d

        subroutine RPN_COMM_ptr_r8_4d(what, ptr)
        use ISO_C_BinDinG
        import :: rpncomm_ptr
        implicit none
        real(KinD=8), dimension(1, 1, 1, 1), intent(in), target :: what
        type(rpncomm_ptr), intent(out) :: ptr
        end subroutine RPN_COMM_ptr_r8_4d
      end interface rpn_comm_ptr


      interface

      subroutine RPN_COMM_qadl(ipe, ni, nj, i0, in, j0, jn)
      implicit none
          integer, intent(out) :: i0, in, j0, jn
          integer, intent(in) :: ipe, ni, nj
      end subroutine RPN_COMM_qadl

      subroutine RPN_COMM_rank( com, rank , ierr )
      implicit none
      character(len=*), intent(in) :: com
      integer, intent(out) :: rank
      integer, intent(out) :: ierr
      end subroutine RPN_COMM_rank

      subroutine RPN_COMM_size( com, size , ierr )
      implicit none
      integer, intent(out) :: ierr, size
      character(len=*), intent(in) ::  com
      end subroutine RPN_COMM_size

      integer function RPN_COMM_split(mex, nx, nxg, minx, maxx, nxl, nxlmax, &
         halox, nx0, fill)
      implicit none
      integer, intent(in) :: nx, mex, nxg, halox
      integer, intent(out) :: minx, maxx, nxl, nxlmax, nx0
      logical, intent(in) :: fill
      end function RPN_COMM_split

      function RPN_COMM_spread(contxt, source, npts, ndata, dest) result(status)
        use ISO_C_BinDinG
        import :: rpncomm_context
        implicit none
        type(rpncomm_context), intent(in) :: contxt
        integer, intent(in) :: npts, ndata
        real, dimension(npts, ndata), intent(in) :: source
        real, dimension(:, :), pointer, intent(inout) :: dest
        integer :: status
      end function RPN_COMM_spread

      function RPN_COMM_spread_context(contxt, com, rootpe, pe, npts) result (status)
        use ISO_C_BinDinG
        import :: rpncomm_context
        implicit none
        type(rpncomm_context), intent(out) :: contxt
        character (len=*), intent(in) :: com
        integer, intent(in) :: npts
        integer, intent(in) :: rootpe
        integer, dimension(npts), intent(in) :: pe
        integer :: status
      end function RPN_COMM_spread_context

      subroutine RPN_COMM_status_size(status_size)
        implicit none
        integer, intent(out) :: status_size
      end subroutine RPN_COMM_status_size

      integer function RPN_COMM_topo_2(nxg, minx, maxx, nxl, nxlmax, halo, &
         nx0, alongx, fill, relax, abort)
      implicit none
      integer, intent(in) :: nxg, halo, relax
      logical, intent(in) :: alongx, fill, abort
      integer, intent(out):: minx, maxx, nxl, nxlmax, nx0
      end  function RPN_COMM_topo_2

      integer function RPN_COMM_topo(nxg, minx, maxx, nxl, nxlmax, halo, nx0, &
         alongx, fill)
      implicit none
      integer, intent(in) :: nxg, halo
      logical, intent(in) :: alongx, fill
      integer, intent(out):: minx, maxx, nxl, nxlmax, nx0
      end function RPN_COMM_topo

      integer function RPN_COMM_topo_xy(nig, njg, mini, maxi, minj, maxj, nil, &
         njl, haloi, haloj, peri, perj)
      implicit none
      integer, intent(in) :: nig, njg, haloi, haloj
      integer, intent(out) :: nil, njl, mini, maxi, minj, maxj
      logical, intent(in) ::  peri, perj
      end function RPN_COMM_topo_xy

      subroutine RPN_COMM_unit(unit)
      implicit none
      integer, intent(in) ::  unit
      end subroutine RPN_COMM_unit

      subroutine RPN_COMM_i_win_group(window, pes_to, pes_from, ierr)
        import :: rpncomm_window
        type(rpncomm_window), intent(in) :: window
        integer, dimension(:), intent(in) :: pes_to
        integer, dimension(:), intent(in) :: pes_from
        integer, intent(out) :: ierr
      end subroutine RPN_COMM_i_win_group

      subroutine RPN_COMM_i_win_oper(window, oper, ierr)
        import :: rpncomm_window
        import :: rpncomm_operator
        type(rpncomm_window), intent(in) :: window
        type(rpncomm_operator), intent(in) :: oper
        integer, intent(out) :: ierr
      end subroutine RPN_COMM_i_win_oper

      subroutine RPN_COMM_i_win_create(window, dtype, siz, com, array, ierr)
        import :: C_PTR
        import :: rpncomm_window, rpncomm_datatype, rpncomm_communicator
        integer, intent(out) :: ierr
        type(rpncomm_window), intent(out) :: window
        type(rpncomm_datatype), intent(in) :: dtype
        integer, intent(in) :: siz
        type(rpncomm_communicator), intent(in) :: com
        type(C_PTR), intent(in), value :: array
      end subroutine RPN_COMM_i_win_create

      subroutine RPN_COMM_i_win_free(window, ierr)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(out) :: ierr
        type(rpncomm_window), intent(inout) :: window
      end subroutine RPN_COMM_i_win_free

      subroutine RPN_COMM_i_win_open(window, active, ierr)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(out) :: ierr
        type(rpncomm_window), intent(in) :: window
        logical, intent(in) :: active
      end subroutine RPN_COMM_i_win_open

      subroutine RPN_COMM_i_win_close(window, ierr)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(out) :: ierr
        type(rpncomm_window), intent(in) :: window
      end subroutine RPN_COMM_i_win_close

      function RPN_COMM_i_valid_win(window, ierr) result (is_valid)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(out) :: ierr
        type(rpncomm_window), intent(in) :: window
        logical :: is_valid
      end function RPN_COMM_i_valid_win

      function RPN_COMM_i_win_check(window, ierr) result (is_open)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(out) :: ierr
        type(rpncomm_window), intent(in) :: window
        logical :: is_open
      end function RPN_COMM_i_win_check

      function RPN_COMM_i_win_get_ptr(window, ierr) result (ptr)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(out) :: ierr
        type(rpncomm_window), intent(in) :: window
        type(C_PTR) :: ptr
      end function RPN_COMM_i_win_get_ptr

      function RPN_COMM_i_win_get_size(window, ierr) result (siz)
        import :: rpncomm_window
        integer, intent(out) :: ierr
        type(rpncomm_window), intent(in) :: window
        integer :: siz
      end function RPN_COMM_i_win_get_size

      subroutine RPN_COMM_i_win_put_r(window, larray, targetpe, offset, nelem, ierr)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(out) :: ierr
        type(rpncomm_window), intent(in) :: window
        type(C_PTR), intent(in), value :: larray
        integer, intent(in) :: targetpe
        integer, intent(in) :: offset
        integer, intent(in) :: nelem
      end subroutine RPN_COMM_i_win_put_r

      subroutine RPN_COMM_i_win_acc_r(window, larray, targetpe, offset, nelem, oper, ierr)
        import :: C_PTR
        import :: rpncomm_window
        import :: rpncomm_operator
        integer, intent(out) :: ierr
        type(rpncomm_window), intent(in) :: window
        type(rpncomm_operator), intent(in) :: oper
        type(C_PTR), intent(in), value :: larray
        integer, intent(in) :: targetpe
        integer, intent(in) :: offset
        integer, intent(in) :: nelem
      end subroutine RPN_COMM_i_win_acc_r

      subroutine RPN_COMM_i_win_put_l(window, larray, offset, nelem, ierr)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(out) :: ierr
        type(rpncomm_window), intent(in) :: window
        type(C_PTR), intent(in), value :: larray
        integer, intent(in) :: offset
        integer, intent(in) :: nelem
      end subroutine RPN_COMM_i_win_put_l

      subroutine RPN_COMM_i_win_get_r(window, larray, target, offset, nelem, ierr)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(out) :: ierr
        type(rpncomm_window), intent(in) :: window
        type(C_PTR), intent(in), value :: larray
        integer, intent(in) :: target
        integer, intent(in) :: offset
        integer, intent(in) :: nelem
      end subroutine RPN_COMM_i_win_get_r

      subroutine RPN_COMM_i_win_get_l(window, larray, offset, nelem, ierr)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(out) :: ierr
        type(rpncomm_window), intent(in) :: window
        type(C_PTR), intent(in), value :: larray
        integer, intent(in) :: offset
        integer, intent(in) :: nelem
      end subroutine RPN_COMM_i_win_get_l

      subroutine RPN_COMM_exchange_halo(pattern, array, periodx, periody, periodp, npol_row)
          import ::  rpncomm_pattern, c_ptr
      implicit none
      type(rpncomm_pattern), intent(in) :: pattern
      type(c_ptr), intent(in) :: array
      logical, intent(in), OPTIONAL  :: periodx, periody, periodp
      integer, intent(in), OPTIONAL :: npol_row
      end subroutine RPN_COMM_exchange_halo

        function rpn_comm_shmget(comm, size) result (where)
        use ISO_C_BinDinG
        implicit none
          integer(C_inT), intent(in) :: comm
          integer(C_inT), intent(in) :: size
          type(C_PTR) :: where
        end function rpn_comm_shmget
      end interface
