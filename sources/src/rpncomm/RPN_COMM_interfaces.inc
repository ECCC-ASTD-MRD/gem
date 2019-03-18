      integer, external :: RPN_COMM_shuf_ezdist  !InTfX!   ! best interface we can provide for the time being
      integer, external :: RPN_COMM_shuf_ezcoll  !InTfX!   ! best interface we can provide for the time being
      interface
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_chdir(string)
      implicit none
      character (len=*) :: string
      end function RPN_COMM_chdir
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_env_var(varname,value)
      implicit none
      character (len=*), intent(IN) :: varname
      character (len=*), intent(OUT) :: value
      end subroutine RPN_COMM_env_var
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_hostid()
      implicit none
      end function RPN_COMM_hostid
!InTf!                                                                                                                                                                                                                                                          
        real*8 function RPN_COMM_wtime()
        implicit none
        end function RPN_COMM_wtime
!InTf!                                                                                                                                                                                                                                                          
        SUBROUTINE RPN_COMM_barrier(com,ierr)
        implicit none
        integer, intent(OUT) :: ierr
        character(len=*), intent(IN) ::  com
      end SUBROUTINE RPN_COMM_barrier
!InTf!                                                                                                                                                                                                                                                          
        integer function RPN_COMM_barrier_data(level,values,nvalues)
        integer , intent(IN) :: level
        integer , intent(IN) :: nvalues
        real *8, dimension(nvalues), intent(OUT) :: values
        end function RPN_COMM_barrier_data
!InTf!                                                                                                                                                                                                                                                          
      SUBROUTINE RPN_COMM_bcastc(buffer,count,datatype,root,com,ierr)
      implicit none
      integer, intent(IN) :: count,root
      integer, intent(OUT) :: ierr
      character(len=*), intent(INOUT) :: buffer
      character(len=*), intent(IN) :: datatype,com
      end SUBROUTINE RPN_COMM_bcastc
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_bloc(nblocx,nblocy)
      implicit none
      integer, intent(IN) :: nblocx, nblocy
      end function RPN_COMM_bloc
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_bloc_find(nblocx,nblocy,set)
      implicit none
      integer, intent(IN) :: nblocx, nblocy
      logical, intent(IN) :: set
      end function RPN_COMM_bloc_find
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_bloc_create(nblocx,nblocy)
      implicit none
      integer, intent(IN) :: nblocx, nblocy
      end function RPN_COMM_bloc_create
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_bloctopo(blocme,blocmex,blocmey,blocsizex,bloc&
     &sizey)
      implicit none
      integer, intent(out) :: blocme,blocmex,blocmey
      integer, intent(out) :: blocsizex, blocsizey
      end subroutine RPN_COMM_bloctopo
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_carac(npex,npey,me,medomm,mex,mey,sizex,sizey,&
     &ismaster, mymaster, mybloc, myblocx,myblocy,blocme,domname)
      implicit none
      integer, intent(out) :: npex,npey,me,mex,mey,sizex,sizey,ismaster
      integer, intent(out) :: mymaster, mybloc, myblocx,myblocy,blocme
      integer, intent(out) :: medomm
      character(len=*), intent(out) :: domname
      end subroutine RPN_COMM_carac
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_colors(comm)
      implicit none
      character(len=*) :: comm
      end function RPN_COMM_colors
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_comm(com)
      implicit none
      character(len=*), intent(IN) :: com
      end function RPN_COMM_comm
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_custom_comm(mpicom,name,mode)
      implicit none
      character(len=*), intent(IN) :: name
      integer, intent(IN) :: mpicom
      integer, intent(IN) :: mode
      end function RPN_COMM_custom_comm
!InTf!                                                                                                                                                                                                                                                          
        subroutine RPN_COMM_i_comm(ctyp_c,ctyp,mcom)
            import :: rpncomm_communicator
        type(rpncomm_communicator), intent(OUT) :: ctyp
        character(len=*), intent(IN) :: ctyp_c
        integer, optional, intent(IN) :: mcom
        end subroutine RPN_COMM_i_comm
!InTf!                                                                                                                                                                                                                                                          
          subroutine RPN_COMM_const(IER)
          implicit none
          integer, intent(OUT) :: IER
          end   subroutine RPN_COMM_const
!InTf!                                                                                                                                                                                                                                                          
        integer function RPN_COMM_datyp_indx(data_int)
        implicit none
        character(len=*), intent(IN) :: data_int
        end function RPN_COMM_datyp_indx
!InTf!                                                                                                                                                                                                                                                          
        integer function RPN_COMM_datyp(data_int)
        implicit none
        character(len=*), intent(IN) :: data_int
        end function RPN_COMM_datyp
!InTf!                                                                                                                                                                                                                                                          
        subroutine RPN_COMM_i_datyp(dtyp_c,dtyp)
            import :: rpncomm_datatype
        type(rpncomm_datatype), intent(OUT) :: dtyp
        character(len=*), intent(IN) :: dtyp_c
        end subroutine RPN_COMM_i_datyp
!InTf!                                                                                                                                                                                                                                                          
        subroutine RPN_COMM_i_user_datyp(dtyp_m,dtyp)
            import :: rpncomm_datatype
        type(rpncomm_datatype), intent(OUT) :: dtyp
        integer, intent(IN) :: dtyp_m
        end subroutine RPN_COMM_i_user_datyp
!InTf!                                                                                                                                                                                                                                                          
        function RPN_COMM_i_valid_datyp(dtyp) result(valid)
            import :: rpncomm_datatype
        type(rpncomm_datatype), intent(IN) :: dtyp
        logical :: valid
        end function RPN_COMM_i_valid_datyp
!InTf!                                                                                                                                                                                                                                                          
        function RPN_COMM_i_datyp_extent(dtyp) result(extent)
            import :: rpncomm_datatype
        type(rpncomm_datatype), intent(IN) :: dtyp
        integer :: extent
        end function RPN_COMM_i_datyp_extent
!InTf!                                                                                                                                                                                                                                                          
      SUBROUTINE RPN_COMM_defo(com)
        implicit none
        character(len=*), intent(IN) ::  com
      end SUBROUTINE RPN_COMM_defo
        subroutine RPN_COMM_fast_dist(garr,gmini,gmaxi,gminj,&
           &          gmaxj,nig,njg,nk,ghalox,ghaloy,size,         &
           &          larr,mini,maxi,minj,maxj,halox,haloy,        &
           &          periodx,periody,status)
        implicit none
        integer, intent(IN) :: ghalox,ghaloy,gmini,gmaxi,gmaxj,gminj
        integer, intent(IN) :: nig,njg,size,mini,maxi,minj,maxj,nk
        integer, intent(OUT)::status
        integer, intent(IN), target :: garr(size,gmini:gmaxi,gminj:gmaxj&
     &,nk),halox,haloy
        integer, intent(OUT):: larr(size,mini:maxi,minj:maxj,nk)
        logical, intent(IN) :: periodx,periody
        end subroutine RPN_COMM_fast_dist
!InTf!                                                                                                                                                                                                                                                          
      SUBROUTINE RPN_COMM_finalize(ierr)
      implicit none
      integer, intent(OUT) ::  ierr
      end SUBROUTINE RPN_COMM_finalize
!InTf!                                                                                                                                                                                                                                                          
        logical function RPN_COMM_grank(com)
        implicit none
        character(len=*), intent(IN) :: com
        end function RPN_COMM_grank
!InTf!                                                                                                                                                                                                                                                          
        logical function RPN_COMM_ngrank(group)
        implicit none
        integer, intent(IN) :: group
        end function RPN_COMM_ngrank
      function rpn_comm_create_2dgrid(gni,gnj,mini,maxi,minj,maxj) resul&
     &t (grid_id)
        integer, intent(IN) :: gni, gnj, mini, maxi, minj, maxj
        integer :: grid_id
      end function rpn_comm_create_2dgrid
      function rpn_comm_get_2dgrid(grid_id,dim_i,dim_j,gni,gnj,mini,maxi&
     &,minj,maxj,starti,counti,startj,countj) result (status)
        integer, intent(IN) :: grid_id, dim_i, dim_j
        integer, intent(OUT) :: gni,gnj,mini,maxi,minj,maxj
        integer, intent(OUT), dimension(dim_i) :: starti,counti
        integer, intent(OUT), dimension(dim_j) :: startj,countj
        integer :: status
      end function rpn_comm_get_2dgrid
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_i_group(com,rgroup,rcom)
        import :: rpncomm_group
        import :: rpncomm_communicator
      implicit none
      character(len=*), intent(IN) :: com
      type(rpncomm_group), intent(OUT) :: rgroup
      type(rpncomm_communicator), optional, intent(IN) :: rcom
      end subroutine RPN_COMM_i_group
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_group(com)
      implicit none
      character(len=*), intent(IN) :: com
      end function RPN_COMM_group
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_haloflip(g,minx,maxx,miny,maxy,ni,nj,nk,halox,&
     &haloy,ggni)
        implicit none
        integer, intent(IN) :: minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,&
     &ggni
        integer, intent(INOUT) :: g(minx:maxx,miny:maxy,nk)
      end subroutine RPN_COMM_haloflip
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_mydomain (call_back, mydomain)
      implicit none
      external :: call_back
      integer, intent(OUT) :: mydomain
      end subroutine RPN_COMM_mydomain
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_world_set(world_comm)
      implicit none
      integer, intent(IN) ::  world_comm
      end subroutine RPN_COMM_world_set
!InTf!                                                                                                                                                                                                                                                          
      SUBROUTINE RPN_COMM_init(Userinit,Pelocal,Petotal,Pex,Pey)
      implicit none
      integer, intent(out)   :: Pelocal,Petotal
      integer, intent(inout) :: Pex,Pey
      external Userinit
      end SUBROUTINE RPN_COMM_init
!InTf!                                                                                                                                                                                                                                                          
      INTEGER FUNCTION RPN_COMM_init_multigrid(Userinit,Pelocal,Petotal,&
     &Pex,Pey,MultiGrids)
      implicit none
      external :: Userinit
      integer, intent(out)   :: Pelocal,Petotal
      integer, intent(inout) :: Pex,Pey
      integer, intent(in)    :: MultiGrids
      end FUNCTION RPN_COMM_init_multigrid
!InTf!                                                                                                                                                                                                                                                          
      INTEGER FUNCTION RPN_COMM_init_multi_level(Userinit,Pelocal,Petota&
     &l,Pex,Pey,MultiGrids,Grids)
      implicit none
      external :: Userinit
      integer, intent(out)   :: Pelocal,Petotal
      integer, intent(inout) :: Pex,Pey
      integer, intent(in)    :: MultiGrids
      integer, intent(in)    :: Grids
      end FUNCTION RPN_COMM_init_multi_level
      integer function RPN_COMM_get_a_free_unit()
      end function RPN_COMM_get_a_free_unit
      function RPN_COMM_set_timeout_alarm(seconds) result(seconds_since)
      integer, intent(IN) :: seconds
      integer :: seconds_since
      end function RPN_COMM_set_timeout_alarm
        subroutine RPN_COMM_make_io_pe_list(x,y,npes,pe_nx,pe_ny,method)
          integer, intent(IN) :: npes
          integer, dimension(npes), intent(OUT) :: x
          integer, dimension(npes), intent(OUT) :: y
          integer, intent(IN) :: pe_nx, pe_ny
          integer, intent(IN) :: method
        end subroutine RPN_COMM_make_io_pe_list
        function RPN_COMM_check_ioset(newset, x ,y, npes, pe_nx, pe_ny, &
     &pe_me, diag) result(status)
          integer, intent(IN) :: newset, npes
          integer, intent(IN) :: pe_nx, pe_ny
          integer, intent(IN) :: pe_me
          logical, intent(IN) :: diag
          integer, intent(IN), dimension(npes) :: x
          integer, intent(IN), dimension(npes) :: y
          integer :: status
        end function RPN_COMM_check_ioset
      function RPN_COMM_create_io_set(npes,method) result(setno)
        implicit none
        integer, intent(IN)  :: npes
        integer, intent(IN)  :: method
        integer :: setno
      end function RPN_COMM_create_io_set
      function RPN_COMM_free_io_set(setno) result(freed)
        implicit none
        integer, intent(IN) :: setno
        integer :: freed
      end function RPN_COMM_free_io_set
      function RPN_COMM_is_valid_io_setno(setno) result(valid)
        implicit none
        integer, intent(IN) :: setno
        logical :: valid
      end function RPN_COMM_is_valid_io_setno
      function RPN_COMM_is_io_pe(setno) result(ordinal)
        implicit none
        integer, intent(IN) :: setno
        integer :: ordinal
      end function RPN_COMM_is_io_pe
      function RPN_COMM_io_pe_gridid(setno,n) result(ordinal)
        implicit none
        integer, intent(IN) :: setno
        integer, intent(IN) :: n
        integer :: ordinal
      end function RPN_COMM_io_pe_gridid
      function RPN_COMM_io_pe_comm(setno) result(communicator)
        implicit none
        integer, intent(IN) :: setno
        integer :: communicator
      end function RPN_COMM_io_pe_comm
      function RPN_COMM_io_pe_size(setno) result(population)
        implicit none
        integer, intent(IN) :: setno
        integer :: population
      end function RPN_COMM_io_pe_size
      function RPN_COMM_io_pe_groups(setno) result(ngroups)
        implicit none
        integer, intent(IN) :: setno
        integer :: ngroups
      end function RPN_COMM_io_pe_groups
      function RPN_COMM_io_pe_idlist(setno) result(idlist)
        implicit none
        integer, intent(IN) :: setno
        integer, dimension(:), pointer :: idlist
      end function RPN_COMM_io_pe_idlist
      function RPN_COMM_io_pe_coord(setno) result(list)
        implicit none
        integer, intent(IN) :: setno
        integer, dimension(:,:), pointer :: list
      end function RPN_COMM_io_pe_coord
      function RPN_COMM_io_pe_callback(setno,callback,argv) result(statu&
     &s)
        import :: C_PTR
        implicit none
        integer, intent(IN) :: setno
        type(C_PTR) :: argv
        integer, external :: callback
        integer :: status
      end function RPN_COMM_io_pe_callback
      function RPN_COMM_io_pe_valid_set(x,y,npes,penx,peny,diag,method) &
     &result(status)
          integer, intent(IN) :: npes
          integer, dimension(npes), intent(OUT) :: x
          integer, dimension(npes), intent(OUT) :: y
          integer, intent(IN) :: penx, peny
          integer, intent(IN) :: method
          logical, intent(IN) :: diag
          integer :: status
      end function RPN_COMM_io_pe_valid_set
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_limit_2(my_id, npe, gmin, gmax,lmini,lma&
     &xi,count, offset,relax)!InTf!
      implicit none
      integer, intent(IN) ::  my_id, npe, gmin, gmax, relax
      integer, intent(OUT) :: lmini,lmaxi
      integer, intent(OUT) :: count(npe),offset(npe)
      end function RPN_COMM_limit_2
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_limit(my_id, npe, gmin, gmax, lmini,lmax&
     &i,count, offset)
      implicit none
      integer, intent(IN) ::  my_id, npe, gmin, gmax
      integer, intent(OUT) :: lmini,lmaxi
      integer, intent(OUT) :: count(npe),offset(npe)
      end function RPN_COMM_limit
!InTf!                                                                                                                                                                                                                                                          
        SUBROUTINE RPN_COMM_UP2LOW(S1,S2)
        implicit none
        character (len=*), intent(IN) :: S1
        character (len=*), intent(OUT) :: S2
        END SUBROUTINE RPN_COMM_UP2LOW
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_mype(Me,Me_x,Me_y)
      implicit none
      integer, intent(OUT) ::  Me,Me_x,Me_y
      end  function RPN_COMM_mype
!InTf!                                                                                                                                                                                                                                                          
        integer function RPN_COMM_oper(op)
        implicit none
        character(len=*), intent(IN) :: op
        end function RPN_COMM_oper
      subroutine RPN_COMM_i_oper(op,r_oper)
       import :: rpncomm_operator
        character(len=*), intent(IN) :: op
        type(rpncomm_operator), intent(OUT) :: r_oper
      end subroutine RPN_COMM_i_oper
      function RPN_COMM_i_valid_oper(r_oper) result (is_valid)
       import :: rpncomm_operator
        type(rpncomm_operator), intent(IN) :: r_oper
        logical :: is_valid
      end function RPN_COMM_i_valid_oper
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_option_L(name_S,value_L)
      implicit none
      character(len=*), intent(IN) :: name_S
      logical, intent(IN) :: value_L
      end function RPN_COMM_option_L
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_option(name_S,value)
      implicit none
      character(len=*), intent(IN) :: name_S
      integer, intent(IN) :: value
      end function RPN_COMM_option
!InTf!                                                                                                                                                                                                                                                          
      SUBROUTINE RPN_COMM_optn(op_type,op_ival,op_rval,op_cval)
      implicit none
      character(len=*), intent(IN) :: op_type, op_cval
      integer, intent(IN) :: op_ival
      real *4, intent(IN) :: op_rval
      end SUBROUTINE RPN_COMM_optn
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_set_petopo(sizx,sizy)
      implicit none
      integer, intent(IN) :: sizx,sizy
      end subroutine RPN_COMM_set_petopo
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_petopo(pex,pey)
      implicit none
      integer, intent(IN) :: pex,pey
      end function RPN_COMM_petopo
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_get_pe(x,y,grd,sgrd,communicator)
      implicit none
      integer, intent(IN) :: x,y,grd,sgrd
      character *(*), intent(IN) :: communicator
      end function RPN_COMM_get_pe
      end interface
       interface rpn_comm_ptr
        subroutine RPN_COMM_ptr_i4_1d(what,ptr)
        use ISO_C_BINDING
        import :: rpncomm_ptr
        implicit none
        integer(KIND=4), dimension(1),intent(IN), target :: what
        type(rpncomm_ptr), intent(OUT) :: ptr
        end subroutine RPN_COMM_ptr_i4_1d
        subroutine RPN_COMM_ptr_i4_2d(what,ptr)
        use ISO_C_BINDING
        import :: rpncomm_ptr
        implicit none
        integer(KIND=4), dimension(1,1),intent(IN), target :: what
        type(rpncomm_ptr), intent(OUT) :: ptr
        end subroutine RPN_COMM_ptr_i4_2d
        subroutine RPN_COMM_ptr_i4_3d(what,ptr)
        use ISO_C_BINDING
        import :: rpncomm_ptr
        implicit none
        integer(KIND=4), dimension(1,1,1),intent(IN), target :: what
        type(rpncomm_ptr), intent(OUT) :: ptr
        end subroutine RPN_COMM_ptr_i4_3d
        subroutine RPN_COMM_ptr_i4_4d(what,ptr)
        use ISO_C_BINDING
        import :: rpncomm_ptr
        implicit none
        integer(KIND=4), dimension(1,1,1,1),intent(IN), target :: what
        type(rpncomm_ptr), intent(OUT) :: ptr
        end subroutine RPN_COMM_ptr_i4_4d
        subroutine RPN_COMM_ptr_i8_1d(what,ptr)
        use ISO_C_BINDING
        import :: rpncomm_ptr
        implicit none
        integer(KIND=8), dimension(1),intent(IN), target :: what
        type(rpncomm_ptr), intent(OUT) :: ptr
        end subroutine RPN_COMM_ptr_i8_1d
        subroutine RPN_COMM_ptr_i8_2d(what,ptr)
        use ISO_C_BINDING
        import :: rpncomm_ptr
        implicit none
        integer(KIND=8), dimension(1,1),intent(IN), target :: what
        type(rpncomm_ptr), intent(OUT) :: ptr
        end subroutine RPN_COMM_ptr_i8_2d
        subroutine RPN_COMM_ptr_i8_3d(what,ptr)
        use ISO_C_BINDING
        import :: rpncomm_ptr
        implicit none
        integer(KIND=8), dimension(1,1,1),intent(IN), target :: what
        type(rpncomm_ptr), intent(OUT) :: ptr
        end subroutine RPN_COMM_ptr_i8_3d
        subroutine RPN_COMM_ptr_i8_4d(what,ptr)
        use ISO_C_BINDING
        import :: rpncomm_ptr
        implicit none
        integer(KIND=8), dimension(1,1,1,1),intent(IN), target :: what
        type(rpncomm_ptr), intent(OUT) :: ptr
        end subroutine RPN_COMM_ptr_i8_4d
        subroutine RPN_COMM_ptr_r4_1d(what,ptr)
        use ISO_C_BINDING
        import :: rpncomm_ptr
        implicit none
        real(KIND=4), dimension(1),intent(IN), target :: what
        type(rpncomm_ptr), intent(OUT) :: ptr
        end subroutine RPN_COMM_ptr_r4_1d
        subroutine RPN_COMM_ptr_r4_2d(what,ptr)
        use ISO_C_BINDING
        import :: rpncomm_ptr
        implicit none
        real(KIND=4), dimension(1,1),intent(IN), target :: what
        type(rpncomm_ptr), intent(OUT) :: ptr
        end subroutine RPN_COMM_ptr_r4_2d
        subroutine RPN_COMM_ptr_r4_3d(what,ptr)
        use ISO_C_BINDING
        import :: rpncomm_ptr
        implicit none
        real(KIND=4), dimension(1,1,1),intent(IN), target :: what
        type(rpncomm_ptr), intent(OUT) :: ptr
        end subroutine RPN_COMM_ptr_r4_3d
        subroutine RPN_COMM_ptr_r4_4d(what,ptr)
        use ISO_C_BINDING
        import :: rpncomm_ptr
        implicit none
        real(KIND=4), dimension(1,1,1,1),intent(IN), target :: what
        type(rpncomm_ptr), intent(OUT) :: ptr
        end subroutine RPN_COMM_ptr_r4_4d
        subroutine RPN_COMM_ptr_r8_1d(what,ptr)
        use ISO_C_BINDING
        import :: rpncomm_ptr
        implicit none
        real(KIND=8), dimension(1),intent(IN), target :: what
        type(rpncomm_ptr), intent(OUT) :: ptr
        end subroutine RPN_COMM_ptr_r8_1d
        subroutine RPN_COMM_ptr_r8_2d(what,ptr)
        use ISO_C_BINDING
        import :: rpncomm_ptr
        implicit none
        real(KIND=8), dimension(1,1),intent(IN), target :: what
        type(rpncomm_ptr), intent(OUT) :: ptr
        end subroutine RPN_COMM_ptr_r8_2d
        subroutine RPN_COMM_ptr_r8_3d(what,ptr)
        use ISO_C_BINDING
        import :: rpncomm_ptr
        implicit none
        real(KIND=8), dimension(1,1,1),intent(IN), target :: what
        type(rpncomm_ptr), intent(OUT) :: ptr
        end subroutine RPN_COMM_ptr_r8_3d
        subroutine RPN_COMM_ptr_r8_4d(what,ptr)
        use ISO_C_BINDING
        import :: rpncomm_ptr
        implicit none
        real(KIND=8), dimension(1,1,1,1),intent(IN), target :: what
        type(rpncomm_ptr), intent(OUT) :: ptr
        end subroutine RPN_COMM_ptr_r8_4d
      end interface rpn_comm_ptr
      interface
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_qadl(ipe,ni,nj,i0,in,j0,jn)
      implicit none
          integer, intent(OUT) :: i0,in,j0,jn
          integer, intent(IN) :: ipe,ni,nj
      end subroutine RPN_COMM_qadl
!InTf!                                                                                                                                                                                                                                                          
      SUBROUTINE RPN_COMM_rank( com, rank ,ierr )
      implicit none
      character(len=*),intent(in) :: com
      integer, intent(out) :: rank
      integer, intent(out) :: ierr
      end SUBROUTINE RPN_COMM_rank
!InTf!                                                                                                                                                                                                                                                          
        SUBROUTINE RPN_COMM_size( com, size ,ierr )
      implicit none
      integer, intent(OUT) :: ierr,size
      character(len=*), intent(IN) ::  com
      end SUBROUTINE RPN_COMM_size
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_split(mex,nx,nxg,minx,maxx,nxl,nxlmax,ha&
     &lox,nx0,fill)
      implicit none
      integer, intent(IN) :: nx, mex, nxg, halox
      integer, intent(OUT) :: minx,maxx,nxl,nxlmax,nx0
      logical, intent(IN) :: fill
      end function RPN_COMM_split
!InTf!                                                                                                                                                                                                                                                          
      function RPN_COMM_spread(contxt,source,npts,ndata,dest)  result(st&
     &atus)
        use ISO_C_BINDING
        import :: rpncomm_context
        implicit none
        type(rpncomm_context), intent(IN) :: contxt
        integer, intent(IN) :: npts, ndata
        real, dimension(npts,ndata), intent(IN) :: source
        real, dimension(:,:), pointer, intent(INOUT) :: dest
        integer :: status
      end function RPN_COMM_spread
!InTf!                                                                                                                                                                                                                                                          
      function RPN_COMM_spread_context(contxt,com,rootpe,pe,npts) result&
     &(status)
        use ISO_C_BINDING
        import :: rpncomm_context
        implicit none
        type(rpncomm_context), intent(OUT) :: contxt
        character (len=*), intent(IN) :: com
        integer, intent(IN) :: npts
        integer, intent(IN) :: rootpe
        integer, dimension(npts), intent(IN) :: pe
        integer :: status
      end function RPN_COMM_spread_context
!InTf!                                                                                                                                                                                                                                                          
        subroutine RPN_COMM_status_size(status_size)
        implicit none
        integer, intent(OUT) :: status_size
      end subroutine RPN_COMM_status_size
!InTf!                                                                                                                                                                                                                                                          
      SUBROUTINE RPN_COMM_tmg_wrt(Iun)
      implicit none
      integer,intent(IN) :: Iun
      end SUBROUTINE RPN_COMM_tmg_wrt
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_topo_2(nxg,minx,maxx,nxl,nxlmax,halo,nx0&
     &,alongx,fill,relax,abort)
      implicit none
      integer, intent(in) :: nxg,halo,relax
      logical, intent(in) :: alongx,fill,abort
      integer, intent(out):: minx,maxx,nxl,nxlmax,nx0
        end  function RPN_COMM_topo_2
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_topo(nxg,minx,maxx,nxl,nxlmax,halo,nx0,a&
     &longx,fill)
      implicit none
      integer, intent(in) :: nxg,halo
      logical, intent(in) :: alongx,fill
      integer, intent(out):: minx,maxx,nxl,nxlmax,nx0
        end function RPN_COMM_topo
!InTf!                                                                                                                                                                                                                                                          
      integer function RPN_COMM_topo_xy(nig,njg,mini,maxi,minj,maxj,nil,&
     &njl,haloi,haloj,peri,perj)
      implicit none
      integer, intent(IN) :: nig,njg,haloi,haloj
      integer, intent(OUT) :: nil,njl,mini,maxi,minj,maxj
      logical, intent(IN) ::  peri,perj
      end function RPN_COMM_topo_xy
!InTf!                                                                                                                                                                                                                                                          
        subroutine RPN_COMM_unit(unit)
        implicit none
        integer, intent(IN) ::  unit
        end subroutine RPN_COMM_unit
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_i_win_group(window,pes_to,pes_from,ierr)
        import :: rpncomm_window
        type(rpncomm_window), intent(IN) :: window
        integer, dimension(:), intent(IN) :: pes_to
        integer, dimension(:), intent(IN) :: pes_from
        integer, intent(OUT) :: ierr
      end subroutine RPN_COMM_i_win_group
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_i_win_oper(window,oper,ierr)
        import :: rpncomm_window
        import :: rpncomm_operator
        type(rpncomm_window), intent(IN) :: window
        type(rpncomm_operator), intent(IN) :: oper
        integer, intent(OUT) :: ierr
      end subroutine RPN_COMM_i_win_oper
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_i_win_create(window,dtype,siz,com,array,ierr)
        import :: C_PTR
        import :: rpncomm_window, rpncomm_datatype, rpncomm_communicator
        integer, intent(OUT) :: ierr
        type(rpncomm_window), intent(OUT) :: window
        type(rpncomm_datatype), intent(IN) :: dtype
        integer, intent(IN) :: siz
        type(rpncomm_communicator), intent(IN) :: com
        type(C_PTR), intent(IN), value :: array
! VOID$ array                                                                                                                                                                                                                                                   
      end subroutine RPN_COMM_i_win_create
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_i_win_free(window,ierr)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(OUT) :: ierr
        type(rpncomm_window), intent(INOUT) :: window
      end subroutine RPN_COMM_i_win_free
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_i_win_open(window,active,ierr)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(OUT) :: ierr
        type(rpncomm_window), intent(IN) :: window
        logical, intent(IN) :: active
      end subroutine RPN_COMM_i_win_open
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_i_win_close(window,ierr)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(OUT) :: ierr
        type(rpncomm_window), intent(IN) :: window
      end subroutine RPN_COMM_i_win_close
!InTf!                                                                                                                                                                                                                                                          
      function RPN_COMM_i_valid_win(window,ierr) result(is_valid)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(OUT) :: ierr
        type(rpncomm_window), intent(IN) :: window
        logical :: is_valid
      end function RPN_COMM_i_valid_win
!InTf!                                                                                                                                                                                                                                                          
      function RPN_COMM_i_win_check(window,ierr) result(is_open)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(OUT) :: ierr
        type(rpncomm_window), intent(IN) :: window
        logical :: is_open
      end function RPN_COMM_i_win_check
!InTf!                                                                                                                                                                                                                                                          
      function RPN_COMM_i_win_get_ptr(window,ierr) result(ptr)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(OUT) :: ierr
        type(rpncomm_window), intent(IN) :: window
        type(C_PTR) :: ptr
      end function RPN_COMM_i_win_get_ptr
!InTf!                                                                                                                                                                                                                                                          
      function RPN_COMM_i_win_get_size(window,ierr) result(siz)
        import :: rpncomm_window
        integer, intent(OUT) :: ierr
        type(rpncomm_window), intent(IN) :: window
        integer :: siz
      end function RPN_COMM_i_win_get_size
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_i_win_put_r(window,larray,targetpe,offset,nele&
     &m,ierr)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(OUT) :: ierr
        type(rpncomm_window), intent(IN) :: window
        type(C_PTR), intent(IN), value :: larray
        integer, intent(IN) :: targetpe
        integer, intent(IN) :: offset
        integer, intent(IN) :: nelem
      end subroutine RPN_COMM_i_win_put_r
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_i_win_acc_r(window,larray,targetpe,offset,nele&
     &m,oper,ierr)
        import :: C_PTR
        import :: rpncomm_window
        import :: rpncomm_operator
        integer, intent(OUT) :: ierr
        type(rpncomm_window), intent(IN) :: window
        type(rpncomm_operator), intent(IN) :: oper
        type(C_PTR), intent(IN), value :: larray
        integer, intent(IN) :: targetpe
        integer, intent(IN) :: offset
        integer, intent(IN) :: nelem
      end subroutine RPN_COMM_i_win_acc_r
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_i_win_put_l(window,larray,offset,nelem,ierr)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(OUT) :: ierr
        type(rpncomm_window), intent(IN) :: window
        type(C_PTR), intent(IN), value :: larray
        integer, intent(IN) :: offset
        integer, intent(IN) :: nelem
      end subroutine RPN_COMM_i_win_put_l
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_i_win_get_r(window,larray,target,offset,nelem,&
     &ierr)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(OUT) :: ierr
        type(rpncomm_window), intent(IN) :: window
        type(C_PTR), intent(IN), value :: larray
        integer, intent(IN) :: target
        integer, intent(IN) :: offset
        integer, intent(IN) :: nelem
      end subroutine RPN_COMM_i_win_get_r
!InTf!                                                                                                                                                                                                                                                          
      subroutine RPN_COMM_i_win_get_l(window,larray,offset,nelem,ierr)
        import :: C_PTR
        import :: rpncomm_window
        integer, intent(OUT) :: ierr
        type(rpncomm_window), intent(IN) :: window
        type(C_PTR), intent(IN), value :: larray
        integer, intent(IN) :: offset
        integer, intent(IN) :: nelem
      end subroutine RPN_COMM_i_win_get_l
!InTf!                                                                                                                                                                                                                                                          
      SUBROUTINE RPN_COMM_exchange_halo(pattern,array,periodx,periody,pe&
     &riodp,npol_row)
          import ::  rpncomm_pattern, c_ptr
      implicit none
      type(rpncomm_pattern), intent(IN) :: pattern
      type(c_ptr), intent(IN) :: array
      logical, intent(IN), OPTIONAL  :: periodx,periody, periodp
      integer, intent(IN), OPTIONAL :: npol_row
      end SUBROUTINE RPN_COMM_exchange_halo
!InTf!                                                                                                                                                                                                                                                          
        function rpn_comm_shmget(comm,size) result(where)
        use ISO_C_BINDING
        implicit none
          integer(C_INT), intent(IN) :: comm
          integer(C_INT), intent(IN) :: size
          type(C_PTR) :: where
        end function rpn_comm_shmget
      end interface
