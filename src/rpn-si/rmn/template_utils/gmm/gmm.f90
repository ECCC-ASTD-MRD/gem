 module gmm_internals
!  GMM_USER_FLAGS
      integer, parameter :: GMM_FLAG_RSTR      =     1    
      integer, parameter :: GMM_FLAG_IZER      =     2    
      integer, parameter :: GMM_FLAG_INAN      =     4    
      integer, parameter :: GMM_FLAG_IINV      =     8    
      integer, parameter :: GMM_FLAG_READ      =    16    
      integer, parameter :: GMM_FLAG_CRTD      =    32    
      integer, parameter :: GMM_FLAG_STAG_X    =    64    
      integer, parameter :: GMM_FLAG_STAG_Y    =   128    
      integer, parameter :: GMM_FLAG_STAG_Z    =   256    
      logical, parameter :: GMM_READ_CKPT=.true.   
      logical, parameter :: GMM_WRIT_CKPT=.false.  
!  GMM VARIABLE SIZES
      integer, parameter :: GMM_MAXNAMELENGTH    =  32   
      integer, parameter :: GMM_META_SIZE =  28
!  GMM_ERROR_CODES
      integer, parameter :: GMM_OK                       = 0     
      integer, parameter :: GMM_ERROR                      = -1   
      integer, parameter :: GMM_KEY_NOT_FOUND            = -2
      integer, parameter :: GMM_VAR_NOT_FOUND            = -3
      integer, parameter :: GMM_INCONSISTENT_DIMS        = -4
      integer, parameter :: GMM_ARRAY_ALREADY_EXISTS     = -5
      integer, parameter :: GMM_VARIABLE_ALREADY_CREATED = -6
      integer, parameter :: GMM_POINTER_TABLE_OVERFLOW   = -7
!  GMM_MESSAGE_LEVELS
      integer, parameter :: GMM_MSG_DEBUG      = -1
      integer, parameter :: GMM_MSG_INFO       =  0
      integer, parameter :: GMM_MSG_WARN       =  1
      integer, parameter :: GMM_MSG_ERROR      =  2
      integer, parameter :: GMM_MSG_SEVERE     =  3
      integer, parameter :: GMM_MSG_FATAL      =  4
      integer, parameter ::  MAX_PAGES    =  16      
      integer, parameter ::  PAGE_SIZE    = 128      
      integer, parameter ::  NTRY_NB_SHFT =   0      
      integer, parameter ::  NTRY_NB_MASK = 127      
      integer, parameter ::  PAGE_NB_SHFT =   7      
      integer, parameter ::  PAGE_NB_MASK =  15      
      integer, parameter ::  EXTN_NB_SHFT =  11      
      integer, parameter ::  EXTN_NB_MASK = 511      
      integer, parameter ::  MAGC_NB_SHFT =  32       
      integer, parameter ::  MAGC_NB_MASK =  -1       
      integer, parameter :: FLAGS_KEPT_ON_CREATE=GMM_FLAG_IZER + GMM_FLAG_INAN + GMM_FLAG_RSTR 
      integer, parameter :: FLAGS_KEPT_IN_RESTART=GMM_FLAG_IZER + GMM_FLAG_INAN + GMM_FLAG_RSTR + GMM_FLAG_IINV 
!
      type gmm_layout                              
         SEQUENCE
         integer :: low,high,halo,halomax,n      
      end type
      type gmm_attributes
        SEQUENCE
        integer*8 :: key          
        integer*8 :: uuid1, uuid2 
        integer   :: initmode                   
        integer   :: flags                      
      end type
      type gmm_metadata
        SEQUENCE
        type(gmm_layout), dimension(4) :: l
        type(gmm_attributes) :: a
      end type
      integer, parameter                            :: GMM_NULL_FLAGS=0
      type(gmm_layout), parameter                   :: GMM_NULL_LAYOUT=gmm_layout(0,0,0,0,0)
      type(gmm_layout), parameter, dimension(4)     :: GMM_NULL_LAYOUTS = (/GMM_NULL_LAYOUT, GMM_NULL_LAYOUT, GMM_NULL_LAYOUT, GMM_&
     &NULL_LAYOUT/)
      type(gmm_attributes), parameter               :: GMM_NULL_ATTRIB=gmm_attributes(GMM_KEY_NOT_FOUND,0,0,0,0)
      type(gmm_metadata), parameter                 :: GMM_NULL_METADATA=gmm_metadata(GMM_NULL_LAYOUTS , GMM_NULL_ATTRIB)
  type p_gmm_metadata
    SEQUENCE
    type(gmm_layout), dimension(4) :: l
    type(gmm_attributes) :: a
    integer data_type
    integer pointer_table_index
    integer *8 array_addr
    character(len=GMM_MAXNAMELENGTH)  :: name               
  end type
  type directory_page
   type(p_gmm_metadata), dimension(:), pointer :: entry
  end type
  type(directory_page), dimension(MAX_PAGES) :: directory
   integer :: used=0                
   integer :: table_size=0          
   integer :: cur_page=0            
   integer :: cur_entry=0           
   integer :: last_entry=PAGE_SIZE  
   integer :: file_unit=0
   logical :: restart_mode=.false.
   integer :: ordinal=0             
   integer :: gmm_verbose_level = 0
 end module gmm_internals
!
  module pointer_table_data_184
  use gmm_internals
  implicit none
  save
    type gmm_p_184
    integer*8, pointer :: p(:,:,:,:)
    integer*8 key
    end type
  type (gmm_p_184) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs184
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    integer*8, pointer :: p(:,:,:,:)
    integer*8, intent(in) :: key
    gmm_ptrs184(gmm_p_table_size)%p => p
    gmm_ptrs184(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs184 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    integer*8, pointer :: p(:,:,:,:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs184(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_184
  module pointer_table_data_144
  use gmm_internals
  implicit none
  save
    type gmm_p_144
    integer*4, pointer :: p(:,:,:,:)
    integer*8 key
    end type
  type (gmm_p_144) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs144
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    integer*4, pointer :: p(:,:,:,:)
    integer*8, intent(in) :: key
    gmm_ptrs144(gmm_p_table_size)%p => p
    gmm_ptrs144(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs144 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    integer*4, pointer :: p(:,:,:,:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs144(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_144
  module pointer_table_data_284
  use gmm_internals
  implicit none
  save
    type gmm_p_284
    real*8, pointer :: p(:,:,:,:)
    integer*8 key
    end type
  type (gmm_p_284) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs284
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    real*8, pointer :: p(:,:,:,:)
    integer*8, intent(in) :: key
    gmm_ptrs284(gmm_p_table_size)%p => p
    gmm_ptrs284(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs284 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    real*8, pointer :: p(:,:,:,:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs284(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_284
  module pointer_table_data_244
  use gmm_internals
  implicit none
  save
    type gmm_p_244
    real*4, pointer :: p(:,:,:,:)
    integer*8 key
    end type
  type (gmm_p_244) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs244
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    real*4, pointer :: p(:,:,:,:)
    integer*8, intent(in) :: key
    gmm_ptrs244(gmm_p_table_size)%p => p
    gmm_ptrs244(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs244 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    real*4, pointer :: p(:,:,:,:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs244(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_244
  module pointer_table_data_384
  use gmm_internals
  implicit none
  save
    type gmm_p_384
    complex*8, pointer :: p(:,:,:,:)
    integer*8 key
    end type
  type (gmm_p_384) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs384
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    complex*8, pointer :: p(:,:,:,:)
    integer*8, intent(in) :: key
    gmm_ptrs384(gmm_p_table_size)%p => p
    gmm_ptrs384(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs384 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    complex*8, pointer :: p(:,:,:,:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs384(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_384
  module pointer_table_data_183
  use gmm_internals
  implicit none
  save
    type gmm_p_183
    integer*8, pointer :: p(:,:,:)
    integer*8 key
    end type
  type (gmm_p_183) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs183
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    integer*8, pointer :: p(:,:,:)
    integer*8, intent(in) :: key
    gmm_ptrs183(gmm_p_table_size)%p => p
    gmm_ptrs183(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs183 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    integer*8, pointer :: p(:,:,:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs183(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_183
  module pointer_table_data_143
  use gmm_internals
  implicit none
  save
    type gmm_p_143
    integer*4, pointer :: p(:,:,:)
    integer*8 key
    end type
  type (gmm_p_143) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs143
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    integer*4, pointer :: p(:,:,:)
    integer*8, intent(in) :: key
    gmm_ptrs143(gmm_p_table_size)%p => p
    gmm_ptrs143(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs143 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    integer*4, pointer :: p(:,:,:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs143(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_143
  module pointer_table_data_283
  use gmm_internals
  implicit none
  save
    type gmm_p_283
    real*8, pointer :: p(:,:,:)
    integer*8 key
    end type
  type (gmm_p_283) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs283
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    real*8, pointer :: p(:,:,:)
    integer*8, intent(in) :: key
    gmm_ptrs283(gmm_p_table_size)%p => p
    gmm_ptrs283(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs283 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    real*8, pointer :: p(:,:,:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs283(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_283
  module pointer_table_data_243
  use gmm_internals
  implicit none
  save
    type gmm_p_243
    real*4, pointer :: p(:,:,:)
    integer*8 key
    end type
  type (gmm_p_243) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs243
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    real*4, pointer :: p(:,:,:)
    integer*8, intent(in) :: key
    gmm_ptrs243(gmm_p_table_size)%p => p
    gmm_ptrs243(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs243 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    real*4, pointer :: p(:,:,:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs243(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_243
  module pointer_table_data_383
  use gmm_internals
  implicit none
  save
    type gmm_p_383
    complex*8, pointer :: p(:,:,:)
    integer*8 key
    end type
  type (gmm_p_383) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs383
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    complex*8, pointer :: p(:,:,:)
    integer*8, intent(in) :: key
    gmm_ptrs383(gmm_p_table_size)%p => p
    gmm_ptrs383(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs383 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    complex*8, pointer :: p(:,:,:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs383(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_383
  module pointer_table_data_182
  use gmm_internals
  implicit none
  save
    type gmm_p_182
    integer*8, pointer :: p(:,:)
    integer*8 key
    end type
  type (gmm_p_182) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs182
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    integer*8, pointer :: p(:,:)
    integer*8, intent(in) :: key
    gmm_ptrs182(gmm_p_table_size)%p => p
    gmm_ptrs182(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs182 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    integer*8, pointer :: p(:,:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs182(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_182
  module pointer_table_data_142
  use gmm_internals
  implicit none
  save
    type gmm_p_142
    integer*4, pointer :: p(:,:)
    integer*8 key
    end type
  type (gmm_p_142) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs142
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    integer*4, pointer :: p(:,:)
    integer*8, intent(in) :: key
    gmm_ptrs142(gmm_p_table_size)%p => p
    gmm_ptrs142(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs142 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    integer*4, pointer :: p(:,:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs142(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_142
  module pointer_table_data_282
  use gmm_internals
  implicit none
  save
    type gmm_p_282
    real*8, pointer :: p(:,:)
    integer*8 key
    end type
  type (gmm_p_282) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs282
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    real*8, pointer :: p(:,:)
    integer*8, intent(in) :: key
    gmm_ptrs282(gmm_p_table_size)%p => p
    gmm_ptrs282(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs282 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    real*8, pointer :: p(:,:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs282(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_282
  module pointer_table_data_242
  use gmm_internals
  implicit none
  save
    type gmm_p_242
    real*4, pointer :: p(:,:)
    integer*8 key
    end type
  type (gmm_p_242) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs242
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    real*4, pointer :: p(:,:)
    integer*8, intent(in) :: key
    gmm_ptrs242(gmm_p_table_size)%p => p
    gmm_ptrs242(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs242 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    real*4, pointer :: p(:,:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs242(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_242
  module pointer_table_data_382
  use gmm_internals
  implicit none
  save
    type gmm_p_382
    complex*8, pointer :: p(:,:)
    integer*8 key
    end type
  type (gmm_p_382) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs382
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    complex*8, pointer :: p(:,:)
    integer*8, intent(in) :: key
    gmm_ptrs382(gmm_p_table_size)%p => p
    gmm_ptrs382(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs382 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    complex*8, pointer :: p(:,:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs382(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_382
  module pointer_table_data_181
  use gmm_internals
  implicit none
  save
    type gmm_p_181
    integer*8, pointer :: p(:)
    integer*8 key
    end type
  type (gmm_p_181) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs181
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    integer*8, pointer :: p(:)
    integer*8, intent(in) :: key
    gmm_ptrs181(gmm_p_table_size)%p => p
    gmm_ptrs181(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs181 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    integer*8, pointer :: p(:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs181(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_181
  module pointer_table_data_141
  use gmm_internals
  implicit none
  save
    type gmm_p_141
    integer*4, pointer :: p(:)
    integer*8 key
    end type
  type (gmm_p_141) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs141
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    integer*4, pointer :: p(:)
    integer*8, intent(in) :: key
    gmm_ptrs141(gmm_p_table_size)%p => p
    gmm_ptrs141(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs141 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    integer*4, pointer :: p(:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs141(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_141
  module pointer_table_data_281
  use gmm_internals
  implicit none
  save
    type gmm_p_281
    real*8, pointer :: p(:)
    integer*8 key
    end type
  type (gmm_p_281) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs281
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    real*8, pointer :: p(:)
    integer*8, intent(in) :: key
    gmm_ptrs281(gmm_p_table_size)%p => p
    gmm_ptrs281(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs281 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    real*8, pointer :: p(:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs281(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_281
  module pointer_table_data_241
  use gmm_internals
  implicit none
  save
    type gmm_p_241
    real*4, pointer :: p(:)
    integer*8 key
    end type
  type (gmm_p_241) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs241
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    real*4, pointer :: p(:)
    integer*8, intent(in) :: key
    gmm_ptrs241(gmm_p_table_size)%p => p
    gmm_ptrs241(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs241 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    real*4, pointer :: p(:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs241(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_241
  module pointer_table_data_381
  use gmm_internals
  implicit none
  save
    type gmm_p_381
    complex*8, pointer :: p(:)
    integer*8 key
    end type
  type (gmm_p_381) , dimension(MAX_PAGES * PAGE_SIZE) :: gmm_ptrs381
   integer :: gmm_p_used=0
   integer :: gmm_p_table_size=0
   integer :: gmm_p_cur_page=0
   integer :: gmm_p_cur_entry=0
   integer :: gmm_p_last_entry=MAX_PAGES * PAGE_SIZE
   integer :: gmm_p_file_unit=0
   logical :: gmm_p_restart_mode=.false.
   integer :: gmm_p_ordinal=0
  contains
  integer function add_table_entry(p, key)
    complex*8, pointer :: p(:)
    integer*8, intent(in) :: key
    gmm_ptrs381(gmm_p_table_size)%p => p
    gmm_ptrs381(gmm_p_table_size)%key = key
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
    print *, 'add_table_entry', ' of', ' gmm_ptrs381 ' , gmm_p_table_size
    endif
    add_table_entry = 0
    return
  end function add_table_entry
  integer function lgmm_get_nxt_avail_ptr()
    lgmm_get_nxt_avail_ptr = gmm_p_table_size + 1
    gmm_p_table_size = gmm_p_table_size + 1
    return
  end function lgmm_get_nxt_avail_ptr
  integer function update_table_entry(indx, key)
    complex*8, pointer :: p(:)
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
    if (indx > gmm_p_table_size) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry : wrong index', indx, gmm_p_table_size
         endif
      update_table_entry = GMM_POINTER_TABLE_OVERFLOW
    endif
!    Cat(gmm_ptrs, EXTENSION,)(indx)%p => p
    gmm_ptrs381(indx)%key = key
        if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *, 'update_table_entry', 'of', indx
        endif
    update_table_entry = 0
    return
  end function update_table_entry
  end module pointer_table_data_381
!!===================== gmm_checkpoint_all =====================
  integer function gmm_checkpoint_all(read_or_write)
!
!        checkpoint read or write for all known types
!
  use gmm_internals
  implicit none
  logical read_or_write
  integer code,istat,fnom
  external fnom
!
  if (read_or_write) then
    if (restart_mode) then
      if (gmm_verbose_level <= GMM_MSG_WARN) then
        print *,'(GMM_CHECKPOINT_ALL) Warning: restart file already read'
      endif
      gmm_checkpoint_all = GMM_OK
      return
    endif
    if (file_unit .eq. 0) then
      istat=fnom(file_unit,'gmm_restart','SEQ+UNF+FTN+OLD',0)
      if (gmm_verbose_level == GMM_MSG_DEBUG) then
        print *,'open restart file, status=',istat
      endif
      if (istat .lt. 0) then
        file_unit = 0
        gmm_checkpoint_all = GMM_ERROR
        return
      endif
    endif
    do while(.true.)
      read(file_unit,end=999)code
      restart_mode=.true.
      if (-1 .eq. code) then
        print *,'ERROR: gmm_checkpoint_all this cannot happen'
      else if (code .eq. 184) then
        call gmm_checkpoint_184(.true.)
      else if (code .eq. 144) then
        call gmm_checkpoint_144(.true.)
      else if (code .eq. 284) then
        call gmm_checkpoint_284(.true.)
      else if (code .eq. 244) then
        call gmm_checkpoint_244(.true.)
      else if (code .eq. 384) then
        call gmm_checkpoint_384(.true.)
      else if (code .eq. 183) then
        call gmm_checkpoint_183(.true.)
      else if (code .eq. 143) then
        call gmm_checkpoint_143(.true.)
      else if (code .eq. 283) then
        call gmm_checkpoint_283(.true.)
      else if (code .eq. 243) then
        call gmm_checkpoint_243(.true.)
      else if (code .eq. 383) then
        call gmm_checkpoint_383(.true.)
      else if (code .eq. 182) then
        call gmm_checkpoint_182(.true.)
      else if (code .eq. 142) then
        call gmm_checkpoint_142(.true.)
      else if (code .eq. 282) then
        call gmm_checkpoint_282(.true.)
      else if (code .eq. 242) then
        call gmm_checkpoint_242(.true.)
      else if (code .eq. 382) then
        call gmm_checkpoint_382(.true.)
      else if (code .eq. 181) then
        call gmm_checkpoint_181(.true.)
      else if (code .eq. 141) then
        call gmm_checkpoint_141(.true.)
      else if (code .eq. 281) then
        call gmm_checkpoint_281(.true.)
      else if (code .eq. 241) then
        call gmm_checkpoint_241(.true.)
      else if (code .eq. 381) then
        call gmm_checkpoint_381(.true.)
      else
        print *,'ERROR: gmm_checkpoint_all unrecognized type=',code,' in restart file'
        call qqexit(1)
      endif
    end do
  else
    if (file_unit .eq. 0) then
      istat=fnom(file_unit,'gmm_restart','SEQ+UNF+FTN',0)
      if (gmm_verbose_level == GMM_MSG_DEBUG) then
        print *,'open restart file, status=',istat
      endif
      if (istat .lt. 0) then
        file_unit = 0
        gmm_checkpoint_all = GMM_ERROR
        return
      endif
    endif
    call gmm_checkpoint_184(.false.)
    call gmm_checkpoint_144(.false.)
    call gmm_checkpoint_284(.false.)
    call gmm_checkpoint_244(.false.)
    call gmm_checkpoint_384(.false.)
    call gmm_checkpoint_183(.false.)
    call gmm_checkpoint_143(.false.)
    call gmm_checkpoint_283(.false.)
    call gmm_checkpoint_243(.false.)
    call gmm_checkpoint_383(.false.)
    call gmm_checkpoint_182(.false.)
    call gmm_checkpoint_142(.false.)
    call gmm_checkpoint_282(.false.)
    call gmm_checkpoint_242(.false.)
    call gmm_checkpoint_382(.false.)
    call gmm_checkpoint_181(.false.)
    call gmm_checkpoint_141(.false.)
    call gmm_checkpoint_281(.false.)
    call gmm_checkpoint_241(.false.)
    call gmm_checkpoint_381(.false.)
  endif
999 call fclos(file_unit)
    file_unit=0
  gmm_checkpoint_all = GMM_OK
  end function gmm_checkpoint_all
!!===================== gmm_checkpoint =====================
!        if  read_or_write is READ_CKPT (.true.) , read one checkpoint group of records
!        if  read_or_write is WRIT_CKPT (.false.) , write all groups of records to checkpoint file
!
  subroutine gmm_checkpoint_184(read_or_write)
  use gmm_internals
  use pointer_table_data_184
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:4) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:4)
    directory(cur_page)%entry(cur_entry)%l(1:4) = siz(1:4)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(184,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs184(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high,&
                                                        &siz(2)%low:siz(2)%high,&
                                                        &siz(3)%low:siz(3)%high,&
                                                        &siz(4)%low:siz(4)%high  ))
    read(file_unit)gmm_ptrs184(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs184(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs184(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs184(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',184
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 184) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)184
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:4)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs184(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_184
  subroutine gmm_checkpoint_144(read_or_write)
  use gmm_internals
  use pointer_table_data_144
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:4) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:4)
    directory(cur_page)%entry(cur_entry)%l(1:4) = siz(1:4)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(144,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs144(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high,&
                                                        &siz(2)%low:siz(2)%high,&
                                                        &siz(3)%low:siz(3)%high,&
                                                        &siz(4)%low:siz(4)%high  ))
    read(file_unit)gmm_ptrs144(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs144(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs144(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs144(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',144
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 144) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)144
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:4)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs144(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_144
  subroutine gmm_checkpoint_284(read_or_write)
  use gmm_internals
  use pointer_table_data_284
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:4) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:4)
    directory(cur_page)%entry(cur_entry)%l(1:4) = siz(1:4)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(284,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs284(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high,&
                                                        &siz(2)%low:siz(2)%high,&
                                                        &siz(3)%low:siz(3)%high,&
                                                        &siz(4)%low:siz(4)%high  ))
    read(file_unit)gmm_ptrs284(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs284(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs284(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs284(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',284
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 284) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)284
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:4)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs284(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_284
  subroutine gmm_checkpoint_244(read_or_write)
  use gmm_internals
  use pointer_table_data_244
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:4) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:4)
    directory(cur_page)%entry(cur_entry)%l(1:4) = siz(1:4)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(244,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs244(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high,&
                                                        &siz(2)%low:siz(2)%high,&
                                                        &siz(3)%low:siz(3)%high,&
                                                        &siz(4)%low:siz(4)%high  ))
    read(file_unit)gmm_ptrs244(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs244(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs244(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs244(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',244
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 244) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)244
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:4)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs244(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_244
  subroutine gmm_checkpoint_384(read_or_write)
  use gmm_internals
  use pointer_table_data_384
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:4) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:4)
    directory(cur_page)%entry(cur_entry)%l(1:4) = siz(1:4)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(384,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs384(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high,&
                                                        &siz(2)%low:siz(2)%high,&
                                                        &siz(3)%low:siz(3)%high,&
                                                        &siz(4)%low:siz(4)%high  ))
    read(file_unit)gmm_ptrs384(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs384(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs384(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs384(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',384
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 384) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)384
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:4)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs384(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_384
  subroutine gmm_checkpoint_183(read_or_write)
  use gmm_internals
  use pointer_table_data_183
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:3) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:3)
    directory(cur_page)%entry(cur_entry)%l(1:3) = siz(1:3)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(183,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs183(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high,&
                                                        &siz(2)%low:siz(2)%high,&
                                                        &siz(3)%low:siz(3)%high  ))
    read(file_unit)gmm_ptrs183(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs183(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs183(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs183(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',183
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 183) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)183
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:3)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs183(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_183
  subroutine gmm_checkpoint_143(read_or_write)
  use gmm_internals
  use pointer_table_data_143
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:3) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:3)
    directory(cur_page)%entry(cur_entry)%l(1:3) = siz(1:3)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(143,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs143(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high,&
                                                        &siz(2)%low:siz(2)%high,&
                                                        &siz(3)%low:siz(3)%high  ))
    read(file_unit)gmm_ptrs143(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs143(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs143(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs143(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',143
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 143) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)143
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:3)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs143(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_143
  subroutine gmm_checkpoint_283(read_or_write)
  use gmm_internals
  use pointer_table_data_283
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:3) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:3)
    directory(cur_page)%entry(cur_entry)%l(1:3) = siz(1:3)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(283,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs283(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high,&
                                                        &siz(2)%low:siz(2)%high,&
                                                        &siz(3)%low:siz(3)%high  ))
    read(file_unit)gmm_ptrs283(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs283(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs283(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs283(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',283
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 283) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)283
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:3)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs283(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_283
  subroutine gmm_checkpoint_243(read_or_write)
  use gmm_internals
  use pointer_table_data_243
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:3) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:3)
    directory(cur_page)%entry(cur_entry)%l(1:3) = siz(1:3)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(243,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs243(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high,&
                                                        &siz(2)%low:siz(2)%high,&
                                                        &siz(3)%low:siz(3)%high  ))
    read(file_unit)gmm_ptrs243(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs243(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs243(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs243(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',243
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 243) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)243
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:3)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs243(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_243
  subroutine gmm_checkpoint_383(read_or_write)
  use gmm_internals
  use pointer_table_data_383
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:3) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:3)
    directory(cur_page)%entry(cur_entry)%l(1:3) = siz(1:3)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(383,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs383(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high,&
                                                        &siz(2)%low:siz(2)%high,&
                                                        &siz(3)%low:siz(3)%high  ))
    read(file_unit)gmm_ptrs383(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs383(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs383(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs383(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',383
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 383) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)383
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:3)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs383(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_383
  subroutine gmm_checkpoint_182(read_or_write)
  use gmm_internals
  use pointer_table_data_182
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:2) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:2)
    directory(cur_page)%entry(cur_entry)%l(1:2) = siz(1:2)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(182,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs182(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high,&
                                                        &siz(2)%low:siz(2)%high  ))
    read(file_unit)gmm_ptrs182(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs182(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs182(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs182(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',182
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 182) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)182
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:2)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs182(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_182
  subroutine gmm_checkpoint_142(read_or_write)
  use gmm_internals
  use pointer_table_data_142
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:2) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:2)
    directory(cur_page)%entry(cur_entry)%l(1:2) = siz(1:2)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(142,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs142(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high,&
                                                        &siz(2)%low:siz(2)%high  ))
    read(file_unit)gmm_ptrs142(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs142(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs142(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs142(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',142
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 142) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)142
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:2)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs142(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_142
  subroutine gmm_checkpoint_282(read_or_write)
  use gmm_internals
  use pointer_table_data_282
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:2) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:2)
    directory(cur_page)%entry(cur_entry)%l(1:2) = siz(1:2)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(282,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs282(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high,&
                                                        &siz(2)%low:siz(2)%high  ))
    read(file_unit)gmm_ptrs282(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs282(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs282(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs282(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',282
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 282) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)282
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:2)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs282(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_282
  subroutine gmm_checkpoint_242(read_or_write)
  use gmm_internals
  use pointer_table_data_242
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:2) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:2)
    directory(cur_page)%entry(cur_entry)%l(1:2) = siz(1:2)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(242,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs242(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high,&
                                                        &siz(2)%low:siz(2)%high  ))
    read(file_unit)gmm_ptrs242(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs242(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs242(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs242(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',242
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 242) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)242
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:2)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs242(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_242
  subroutine gmm_checkpoint_382(read_or_write)
  use gmm_internals
  use pointer_table_data_382
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:2) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:2)
    directory(cur_page)%entry(cur_entry)%l(1:2) = siz(1:2)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(382,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs382(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high,&
                                                        &siz(2)%low:siz(2)%high  ))
    read(file_unit)gmm_ptrs382(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs382(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs382(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs382(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',382
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 382) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)382
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:2)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs382(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_382
  subroutine gmm_checkpoint_181(read_or_write)
  use gmm_internals
  use pointer_table_data_181
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:1) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:1)
    directory(cur_page)%entry(cur_entry)%l(1:1) = siz(1:1)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(181,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs181(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high  ))
    read(file_unit)gmm_ptrs181(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs181(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs181(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs181(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',181
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 181) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)181
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:1)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs181(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_181
  subroutine gmm_checkpoint_141(read_or_write)
  use gmm_internals
  use pointer_table_data_141
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:1) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:1)
    directory(cur_page)%entry(cur_entry)%l(1:1) = siz(1:1)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(141,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs141(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high  ))
    read(file_unit)gmm_ptrs141(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs141(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs141(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs141(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',141
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 141) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)141
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:1)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs141(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_141
  subroutine gmm_checkpoint_281(read_or_write)
  use gmm_internals
  use pointer_table_data_281
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:1) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:1)
    directory(cur_page)%entry(cur_entry)%l(1:1) = siz(1:1)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(281,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs281(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high  ))
    read(file_unit)gmm_ptrs281(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs281(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs281(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs281(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',281
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 281) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)281
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:1)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs281(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_281
  subroutine gmm_checkpoint_241(read_or_write)
  use gmm_internals
  use pointer_table_data_241
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:1) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:1)
    directory(cur_page)%entry(cur_entry)%l(1:1) = siz(1:1)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(241,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs241(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high  ))
    read(file_unit)gmm_ptrs241(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs241(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs241(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs241(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',241
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 241) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)241
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:1)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs241(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_241
  subroutine gmm_checkpoint_381(read_or_write)
  use gmm_internals
  use pointer_table_data_381
  implicit none
  logical read_or_write
  integer istat, fnom, i, j, ier, lcl_pti
  type(gmm_layout), dimension(1:1) :: siz
  type(gmm_attributes) :: attrib
  integer *8 :: key
  external fnom
      integer *8 get_address_from
      external get_address_from
  if (read_or_write) then
!
    call add_directory_entry
!    read(file_unit)directory(cur_page)%entry(cur_entry)%m%a%name  ! read name of variable
    read(file_unit)directory(cur_page)%entry(cur_entry)%name
    read(file_unit)siz(1:1)
    directory(cur_page)%entry(cur_entry)%l(1:1) = siz(1:1)
    read(file_unit)attrib
!    print *,'name=',directory(cur_page)%entry(cur_entry)%name,' dims=',siz(1:DIM)
    attrib%flags = ior(attrib%flags,GMM_FLAG_READ)
    directory(cur_page)%entry(cur_entry)%a = attrib
    read(file_unit)directory(cur_page)%entry(cur_entry)%data_type
    lcl_pti = lgmm_get_nxt_avail_ptr()
    directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
    ordinal = ordinal + 1
    key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
    key = key + ishft(381,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
    directory(cur_page)%entry(cur_entry)%a%key = key
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
    allocate(gmm_ptrs381(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p(siz(1)%low:siz(1)%high  ))
    read(file_unit)gmm_ptrs381(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(&
          &gmm_ptrs381(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'name=',directory(cur_page)%entry(cur_entry)%name,' cur_page=',cur_page,' cur_entry=',&
     &cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(gmm_ptrs381(directory(&
     &cur_page)%entry(cur_entry)%pointer_table_index)%p)
    endif
    ier=add_table_entry(gmm_ptrs381(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p,key)
!
  else
!
    if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *,'checkpointing type ',381
    endif
    do i=1,table_size
      do j=1,PAGE_SIZE
        if (iand(GMM_FLAG_RSTR,directory(i)%entry(j)%a%flags) .ne. 0.and.directory(i)%entry(j)%data_type == 381) then
          if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'writing field ',directory(i)%entry(j)%name
          endif
          write(file_unit)381
          write(file_unit)directory(i)%entry(j)%name
          write(file_unit)directory(i)%entry(j)%l(1:1)
          attrib = directory(i)%entry(j)%a
          attrib%flags = iand(attrib%flags,FLAGS_KEPT_IN_RESTART)
          write(file_unit)attrib
          write(file_unit)directory(i)%entry(j)%data_type
          write(file_unit)gmm_ptrs381(directory(i)%entry(j)%pointer_table_index)%p
        endif
      enddo
    enddo
!
  endif
  end subroutine gmm_checkpoint_381
!!===================== gmm_create (interface) =====================
!
!!===================== gmm_create (code) =====================
   integer function gmm_create184(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_184
   implicit none
   character(len=*), intent(in) :: iname
   integer*8, pointer :: p(:,:,:,:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   integer*8, pointer :: pp(:,:,:,:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 184
   if (associated(p)) then
      consistent=.true.
      do i=1,4
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create184 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create184 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs184(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,4
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create184 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create184 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs184(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create184 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(184,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:4) = dims(1:4)
   directory(cur_page)%entry(cur_entry)%data_type = 184
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high,&
                &dims(2)%low:dims(2)%high,&
                &dims(3)%low:dims(3)%high,&
                &dims(4)%low:dims(4)%high),stat=ier)
      if (ier /= 0) then
         gmm_create184 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
         print *,'GMM_CREATE ERROR, name=',iname,' : init to NaN is not available for this data type'
         gmm_create184 = GMM_ERROR
         return
   endif
   gmm_create184 = 0
   end function gmm_create184
!
   integer function gmm_create144(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_144
   implicit none
   character(len=*), intent(in) :: iname
   integer*4, pointer :: p(:,:,:,:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   integer*4, pointer :: pp(:,:,:,:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 144
   if (associated(p)) then
      consistent=.true.
      do i=1,4
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create144 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create144 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs144(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,4
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create144 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create144 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs144(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create144 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(144,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:4) = dims(1:4)
   directory(cur_page)%entry(cur_entry)%data_type = 144
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high,&
                &dims(2)%low:dims(2)%high,&
                &dims(3)%low:dims(3)%high,&
                &dims(4)%low:dims(4)%high),stat=ier)
      if (ier /= 0) then
         gmm_create144 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
         print *,'GMM_CREATE ERROR, name=',iname,' : init to NaN is not available for this data type'
         gmm_create144 = GMM_ERROR
         return
   endif
   gmm_create144 = 0
   end function gmm_create144
!
   integer function gmm_create284(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_284
   implicit none
   character(len=*), intent(in) :: iname
   real*8, pointer :: p(:,:,:,:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   real*8, pointer :: pp(:,:,:,:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 284
   if (associated(p)) then
      consistent=.true.
      do i=1,4
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create284 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create284 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs284(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,4
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create284 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create284 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs284(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create284 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(284,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:4) = dims(1:4)
   directory(cur_page)%entry(cur_entry)%data_type = 284
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high,&
                &dims(2)%low:dims(2)%high,&
                &dims(3)%low:dims(3)%high,&
                &dims(4)%low:dims(4)%high),stat=ier)
      if (ier /= 0) then
         gmm_create284 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
       if (gmm_verbose_level == GMM_MSG_DEBUG)  print *,iname,' Debug DATATYPE=real*8 init to NaN8'
       p = NaN8
   endif
   gmm_create284 = 0
   end function gmm_create284
!
   integer function gmm_create244(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_244
   implicit none
   character(len=*), intent(in) :: iname
   real*4, pointer :: p(:,:,:,:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   real*4, pointer :: pp(:,:,:,:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 244
   if (associated(p)) then
      consistent=.true.
      do i=1,4
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create244 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create244 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs244(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,4
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create244 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create244 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs244(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create244 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(244,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:4) = dims(1:4)
   directory(cur_page)%entry(cur_entry)%data_type = 244
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high,&
                &dims(2)%low:dims(2)%high,&
                &dims(3)%low:dims(3)%high,&
                &dims(4)%low:dims(4)%high),stat=ier)
      if (ier /= 0) then
         gmm_create244 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
     if (gmm_verbose_level == GMM_MSG_DEBUG)  print *,iname,' Debug DATATYPE=real init to NaN'
     p = NaN
   endif
   gmm_create244 = 0
   end function gmm_create244
!
   integer function gmm_create384(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_384
   implicit none
   character(len=*), intent(in) :: iname
   complex*8, pointer :: p(:,:,:,:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   complex*8, pointer :: pp(:,:,:,:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 384
   if (associated(p)) then
      consistent=.true.
      do i=1,4
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create384 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create384 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs384(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,4
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create384 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create384 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs384(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create384 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(384,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:4) = dims(1:4)
   directory(cur_page)%entry(cur_entry)%data_type = 384
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high,&
                &dims(2)%low:dims(2)%high,&
                &dims(3)%low:dims(3)%high,&
                &dims(4)%low:dims(4)%high),stat=ier)
      if (ier /= 0) then
         gmm_create384 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
         print *,'GMM_CREATE ERROR, name=',iname,' : init to NaN is not available for this data type'
         gmm_create384 = GMM_ERROR
         return
   endif
   gmm_create384 = 0
   end function gmm_create384
!
   integer function gmm_create183(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_183
   implicit none
   character(len=*), intent(in) :: iname
   integer*8, pointer :: p(:,:,:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   integer*8, pointer :: pp(:,:,:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 183
   if (associated(p)) then
      consistent=.true.
      do i=1,3
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create183 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create183 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs183(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,3
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create183 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create183 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs183(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create183 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(183,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:3) = dims(1:3)
   directory(cur_page)%entry(cur_entry)%data_type = 183
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high,&
                &dims(2)%low:dims(2)%high,&
                &dims(3)%low:dims(3)%high),stat=ier)
      if (ier /= 0) then
         gmm_create183 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
         print *,'GMM_CREATE ERROR, name=',iname,' : init to NaN is not available for this data type'
         gmm_create183 = GMM_ERROR
         return
   endif
   gmm_create183 = 0
   end function gmm_create183
!
   integer function gmm_create143(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_143
   implicit none
   character(len=*), intent(in) :: iname
   integer*4, pointer :: p(:,:,:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   integer*4, pointer :: pp(:,:,:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 143
   if (associated(p)) then
      consistent=.true.
      do i=1,3
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create143 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create143 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs143(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,3
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create143 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create143 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs143(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create143 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(143,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:3) = dims(1:3)
   directory(cur_page)%entry(cur_entry)%data_type = 143
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high,&
                &dims(2)%low:dims(2)%high,&
                &dims(3)%low:dims(3)%high),stat=ier)
      if (ier /= 0) then
         gmm_create143 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
         print *,'GMM_CREATE ERROR, name=',iname,' : init to NaN is not available for this data type'
         gmm_create143 = GMM_ERROR
         return
   endif
   gmm_create143 = 0
   end function gmm_create143
!
   integer function gmm_create283(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_283
   implicit none
   character(len=*), intent(in) :: iname
   real*8, pointer :: p(:,:,:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   real*8, pointer :: pp(:,:,:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 283
   if (associated(p)) then
      consistent=.true.
      do i=1,3
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create283 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create283 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs283(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,3
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create283 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create283 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs283(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create283 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(283,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:3) = dims(1:3)
   directory(cur_page)%entry(cur_entry)%data_type = 283
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high,&
                &dims(2)%low:dims(2)%high,&
                &dims(3)%low:dims(3)%high),stat=ier)
      if (ier /= 0) then
         gmm_create283 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
       if (gmm_verbose_level == GMM_MSG_DEBUG)  print *,iname,' Debug DATATYPE=real*8 init to NaN8'
       p = NaN8
   endif
   gmm_create283 = 0
   end function gmm_create283
!
   integer function gmm_create243(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_243
   implicit none
   character(len=*), intent(in) :: iname
   real*4, pointer :: p(:,:,:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   real*4, pointer :: pp(:,:,:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 243
   if (associated(p)) then
      consistent=.true.
      do i=1,3
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create243 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create243 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs243(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,3
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create243 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create243 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs243(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create243 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(243,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:3) = dims(1:3)
   directory(cur_page)%entry(cur_entry)%data_type = 243
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high,&
                &dims(2)%low:dims(2)%high,&
                &dims(3)%low:dims(3)%high),stat=ier)
      if (ier /= 0) then
         gmm_create243 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
     if (gmm_verbose_level == GMM_MSG_DEBUG)  print *,iname,' Debug DATATYPE=real init to NaN'
     p = NaN
   endif
   gmm_create243 = 0
   end function gmm_create243
!
   integer function gmm_create383(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_383
   implicit none
   character(len=*), intent(in) :: iname
   complex*8, pointer :: p(:,:,:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   complex*8, pointer :: pp(:,:,:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 383
   if (associated(p)) then
      consistent=.true.
      do i=1,3
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create383 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create383 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs383(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,3
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create383 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create383 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs383(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create383 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(383,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:3) = dims(1:3)
   directory(cur_page)%entry(cur_entry)%data_type = 383
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high,&
                &dims(2)%low:dims(2)%high,&
                &dims(3)%low:dims(3)%high),stat=ier)
      if (ier /= 0) then
         gmm_create383 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
         print *,'GMM_CREATE ERROR, name=',iname,' : init to NaN is not available for this data type'
         gmm_create383 = GMM_ERROR
         return
   endif
   gmm_create383 = 0
   end function gmm_create383
!
   integer function gmm_create182(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_182
   implicit none
   character(len=*), intent(in) :: iname
   integer*8, pointer :: p(:,:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   integer*8, pointer :: pp(:,:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 182
   if (associated(p)) then
      consistent=.true.
      do i=1,2
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create182 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create182 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs182(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,2
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create182 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create182 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs182(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create182 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(182,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:2) = dims(1:2)
   directory(cur_page)%entry(cur_entry)%data_type = 182
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high,&
                &dims(2)%low:dims(2)%high),stat=ier)
      if (ier /= 0) then
         gmm_create182 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
         print *,'GMM_CREATE ERROR, name=',iname,' : init to NaN is not available for this data type'
         gmm_create182 = GMM_ERROR
         return
   endif
   gmm_create182 = 0
   end function gmm_create182
!
   integer function gmm_create142(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_142
   implicit none
   character(len=*), intent(in) :: iname
   integer*4, pointer :: p(:,:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   integer*4, pointer :: pp(:,:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 142
   if (associated(p)) then
      consistent=.true.
      do i=1,2
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create142 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create142 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs142(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,2
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create142 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create142 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs142(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create142 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(142,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:2) = dims(1:2)
   directory(cur_page)%entry(cur_entry)%data_type = 142
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high,&
                &dims(2)%low:dims(2)%high),stat=ier)
      if (ier /= 0) then
         gmm_create142 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
         print *,'GMM_CREATE ERROR, name=',iname,' : init to NaN is not available for this data type'
         gmm_create142 = GMM_ERROR
         return
   endif
   gmm_create142 = 0
   end function gmm_create142
!
   integer function gmm_create282(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_282
   implicit none
   character(len=*), intent(in) :: iname
   real*8, pointer :: p(:,:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   real*8, pointer :: pp(:,:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 282
   if (associated(p)) then
      consistent=.true.
      do i=1,2
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create282 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create282 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs282(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,2
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create282 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create282 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs282(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create282 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(282,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:2) = dims(1:2)
   directory(cur_page)%entry(cur_entry)%data_type = 282
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high,&
                &dims(2)%low:dims(2)%high),stat=ier)
      if (ier /= 0) then
         gmm_create282 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
       if (gmm_verbose_level == GMM_MSG_DEBUG)  print *,iname,' Debug DATATYPE=real*8 init to NaN8'
       p = NaN8
   endif
   gmm_create282 = 0
   end function gmm_create282
!
   integer function gmm_create242(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_242
   implicit none
   character(len=*), intent(in) :: iname
   real*4, pointer :: p(:,:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   real*4, pointer :: pp(:,:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 242
   if (associated(p)) then
      consistent=.true.
      do i=1,2
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create242 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create242 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs242(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,2
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create242 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create242 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs242(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create242 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(242,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:2) = dims(1:2)
   directory(cur_page)%entry(cur_entry)%data_type = 242
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high,&
                &dims(2)%low:dims(2)%high),stat=ier)
      if (ier /= 0) then
         gmm_create242 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
     if (gmm_verbose_level == GMM_MSG_DEBUG)  print *,iname,' Debug DATATYPE=real init to NaN'
     p = NaN
   endif
   gmm_create242 = 0
   end function gmm_create242
!
   integer function gmm_create382(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_382
   implicit none
   character(len=*), intent(in) :: iname
   complex*8, pointer :: p(:,:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   complex*8, pointer :: pp(:,:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 382
   if (associated(p)) then
      consistent=.true.
      do i=1,2
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create382 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create382 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs382(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,2
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create382 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create382 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs382(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create382 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(382,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:2) = dims(1:2)
   directory(cur_page)%entry(cur_entry)%data_type = 382
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high,&
                &dims(2)%low:dims(2)%high),stat=ier)
      if (ier /= 0) then
         gmm_create382 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
         print *,'GMM_CREATE ERROR, name=',iname,' : init to NaN is not available for this data type'
         gmm_create382 = GMM_ERROR
         return
   endif
   gmm_create382 = 0
   end function gmm_create382
!
   integer function gmm_create181(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_181
   implicit none
   character(len=*), intent(in) :: iname
   integer*8, pointer :: p(:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   integer*8, pointer :: pp(:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 181
   if (associated(p)) then
      consistent=.true.
      do i=1,1
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create181 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create181 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs181(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,1
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create181 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create181 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs181(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create181 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(181,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:1) = dims(1:1)
   directory(cur_page)%entry(cur_entry)%data_type = 181
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high),stat=ier)
      if (ier /= 0) then
         gmm_create181 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
         print *,'GMM_CREATE ERROR, name=',iname,' : init to NaN is not available for this data type'
         gmm_create181 = GMM_ERROR
         return
   endif
   gmm_create181 = 0
   end function gmm_create181
!
   integer function gmm_create141(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_141
   implicit none
   character(len=*), intent(in) :: iname
   integer*4, pointer :: p(:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   integer*4, pointer :: pp(:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 141
   if (associated(p)) then
      consistent=.true.
      do i=1,1
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create141 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create141 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs141(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,1
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create141 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create141 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs141(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create141 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(141,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:1) = dims(1:1)
   directory(cur_page)%entry(cur_entry)%data_type = 141
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high),stat=ier)
      if (ier /= 0) then
         gmm_create141 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
         print *,'GMM_CREATE ERROR, name=',iname,' : init to NaN is not available for this data type'
         gmm_create141 = GMM_ERROR
         return
   endif
   gmm_create141 = 0
   end function gmm_create141
!
   integer function gmm_create281(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_281
   implicit none
   character(len=*), intent(in) :: iname
   real*8, pointer :: p(:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   real*8, pointer :: pp(:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 281
   if (associated(p)) then
      consistent=.true.
      do i=1,1
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create281 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create281 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs281(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,1
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create281 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create281 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs281(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create281 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(281,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:1) = dims(1:1)
   directory(cur_page)%entry(cur_entry)%data_type = 281
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high),stat=ier)
      if (ier /= 0) then
         gmm_create281 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
       if (gmm_verbose_level == GMM_MSG_DEBUG)  print *,iname,' Debug DATATYPE=real*8 init to NaN8'
       p = NaN8
   endif
   gmm_create281 = 0
   end function gmm_create281
!
   integer function gmm_create241(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_241
   implicit none
   character(len=*), intent(in) :: iname
   real*4, pointer :: p(:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   real*4, pointer :: pp(:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 241
   if (associated(p)) then
      consistent=.true.
      do i=1,1
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create241 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create241 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs241(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,1
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create241 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create241 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs241(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create241 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(241,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:1) = dims(1:1)
   directory(cur_page)%entry(cur_entry)%data_type = 241
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high),stat=ier)
      if (ier /= 0) then
         gmm_create241 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
     if (gmm_verbose_level == GMM_MSG_DEBUG)  print *,iname,' Debug DATATYPE=real init to NaN'
     p = NaN
   endif
   gmm_create241 = 0
   end function gmm_create241
!
   integer function gmm_create381(iname,p,field_meta,flags_arg)
   use gmm_internals
   use pointer_table_data_381
   implicit none
   character(len=*), intent(in) :: iname
   complex*8, pointer :: p(:)
   type(gmm_metadata), intent(inout) :: field_meta
   integer, intent(in), optional :: flags_arg
   integer *8 get_address_from
   external get_address_from
   external fool_optimizer
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   type(gmm_attributes) :: localattr, attrs
   type (gmm_attributes) lcl_attr
   type(gmm_layout), dimension(4) :: lcl_layout, dims
   integer lcl_datatype
   integer*8 :: key
   logical consistent
   integer i, ier
   integer lcl_pti
   character(len=GMM_MAXNAMELENGTH) :: lcl_name
   integer u_bound, l_bound
   complex*8, pointer :: pp(:)
    real      :: NaN
    integer   :: inan
    real*8    :: NaN8
    integer*8 :: inan8
    equivalence (inan,NaN)
    equivalence (inan8,NaN8)
    data inan  /Z'7F800001'/
    data inan8 /Z'7FF0000000000001'/
   if (present(flags_arg)) then
    field_meta%a%flags = flags_arg
   endif
   lcl_layout = field_meta%l
   dims = lcl_layout
   lcl_attr   = field_meta%a
   attrs = lcl_attr
   lcl_datatype = 381
   if (associated(p)) then
      consistent=.true.
      do i=1,1
         consistent = consistent .and. size(p,i).eq.(dims(i)%high-dims(i)%low+1)
      enddo
      if (.not. consistent ) then
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
            print *,'ERROR: gmm_create, p has dimensions that are not consistent with dims'
         endif
         key=0
         gmm_create381 = GMM_INCONSISTENT_DIMS
         return
      endif
   endif
   localattr = attrs
   lcl_name = trim(iname)
   localattr%flags = iand(localattr%flags,FLAGS_KEPT_ON_CREATE)
!   call find_directory_entry(localattr%name,key)     ! is there a field with this name that exists ?
   call find_directory_entry(lcl_name,key)
   if (cur_page .ne. 0 .and. cur_entry .ne. 0) then
      if (associated(p)) then
         print *,'ERROR: gmm_create called with existing p and array has already been created'
         key=0
         gmm_create381 = GMM_ARRAY_ALREADY_EXISTS
         return
      endif
      pp=>gmm_ptrs381(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
      consistent=.true.
      do i=1,1
         call fool_optimizer(size(pp,i))
         consistent = consistent .and. (size(pp,i).eq.(dims(i)%high-dims(i)%low+1))
         if (.not. consistent ) print *,'size(pp,',i,')=',size(pp,i),' high=',dims(i)%high,' low=',dims(i)%low
      enddo
      if (.not. consistent ) then
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=1 =',lbound(pp,1),ubound(pp,1)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=2 =',lbound(pp,2),ubound(pp,2)
!         write(6,'(a,2x,a8,2x,a,2x,i3,2x,i3)') 'Debug+++ ',lcl_name,' bounds dim=3 =',lbound(pp,3),ubound(pp,3)
         print *,'ERROR: gmm_create, requested dimensions differ from previous specification (restart/create)'
         print *,'ERROR: gmm_create, variable name ="',lcl_name,'"'
         key=0
         nullify(p)
         gmm_create381 = GMM_INCONSISTENT_DIMS
         return
      else
         if (gmm_verbose_level == GMM_MSG_DEBUG) then
           print *,'INFO: gmm_create, variable name =',lcl_name,' exists and is consistent'
         endif
      endif
      if (iand(GMM_FLAG_CRTD,directory(cur_page)%entry(cur_entry)%a%flags) .ne. 0) then
         print *,'ERROR: gmm_create, field ',lcl_name,' has already been created'
         key = 0
         nullify(p)
         gmm_create381 = GMM_VARIABLE_ALREADY_CREATED
         return
!
      else
!         print *,'Debug+ ', \' Cat(gmm_create, EXTENSION,) \','array ',lcl_name,'must then have been read from a restart file' 
         localattr%flags = ior(localattr%flags,directory(cur_page)%entry(cur_entry)%a%flags)
         key=directory(cur_page)%entry(cur_entry)%a%key
         directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
         p=>gmm_ptrs381(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i20)') 'Debug+++ gmm_create name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p)
         directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
         gmm_create381 = 0
         return
      endif
   else
      if (iand(GMM_FLAG_RSTR,localattr%flags).ne.0 .and. restart_mode) then
         print *,'ERROR: gmm_create field ',lcl_name, 'should have been read from restart file but was not'
!  ==========  HOW SERIOUS AN ERROR IS THIS ? ===============
      endif
      call add_directory_entry
   endif
   ordinal = ordinal + 1
   key = ishft((cur_page-1),PAGE_NB_SHFT) + ishft((cur_entry-1),NTRY_NB_SHFT)
   key = key + ishft(381,EXTN_NB_SHFT) + ishft(ordinal,MAGC_NB_SHFT)
   directory(cur_page)%entry(cur_entry)%a = localattr
! CODE POSSIBLY MISSING HERE FOR FLAGS SETTINGS
   directory(cur_page)%entry(cur_entry)%name = lcl_name
   directory(cur_page)%entry(cur_entry)%a%key = key
   directory(cur_page)%entry(cur_entry)%a%flags = ior(localattr%flags,GMM_FLAG_CRTD)
   directory(cur_page)%entry(cur_entry)%l(1:1) = dims(1:1)
   directory(cur_page)%entry(cur_entry)%data_type = 381
   if (associated(p)) then
      if (gmm_verbose_level == GMM_MSG_DEBUG) &
        print *,'GMM_CREATE: using user supplied array'
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      ier = add_table_entry(p, key)
! ======= must check that certain attributes are not requested (e.g. FLAG_RSTR) and that size is consistent
   else
      lcl_pti = lgmm_get_nxt_avail_ptr()
      directory(cur_page)%entry(cur_entry)%pointer_table_index = lcl_pti
      allocate(p(dims(1)%low:dims(1)%high),stat=ier)
      if (ier /= 0) then
         gmm_create381 = GMM_ERROR
         return
      endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_create creation name=',lcl_name,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
      directory(cur_page)%entry(cur_entry)%array_addr = get_address_from(p)
      ier = add_table_entry(p, key)
   endif
   field_meta%l = directory(cur_page)%entry(cur_entry)%l
   field_meta%a = directory(cur_page)%entry(cur_entry)%a
!    Cat(gmm_ptrs, EXTENSION,)(directory(cur_page)%entry(cur_entry)%f => p
!         if (iand(directory(cur_page)%entry(cur_entry)%a%flags,FLAG_IZER) .ne. 0) then
!           directory(cur_page)%entry(cur_entry)%f = 0        ! ZERO fill requested
!         endif
   if (iand(field_meta%a%flags, GMM_FLAG_IZER) /= 0) THEN
!     print *,'Debug+ gmm_create init to zero for variable ',lcl_name
     p = 0
   endif
   if (iand(field_meta%a%flags, GMM_FLAG_INAN) /= 0) THEN
         print *,'GMM_CREATE ERROR, name=',iname,' : init to NaN is not available for this data type'
         gmm_create381 = GMM_ERROR
         return
   endif
   gmm_create381 = 0
   end function gmm_create381
!
  subroutine check_directory_entry(name,key)
  use gmm_internals
  implicit none
  character(len=*) :: name
  integer*8, intent(in) :: key
!
  character(len=GMM_MAXNAMELENGTH) :: l_name
  integer temp
  logical found
!
  found = .false.
  if (cur_page == 0 .and. cur_entry == 0) then
    return
  endif
  l_name = trim(name)
  temp = ishft(key,-PAGE_NB_SHFT)
  cur_page = iand(PAGE_NB_MASK,temp)
  cur_page = min(cur_page+1,table_size)
  temp = ishft(key,-NTRY_NB_SHFT)
  cur_entry = iand(NTRY_NB_MASK,temp)
  cur_entry = min(cur_entry+1,PAGE_SIZE)
  found = key .eq. directory(cur_page)%entry(cur_entry)%a%key
  found = found .and. ( directory(cur_page)%entry(cur_entry)%name .eq. l_name )
  if (.not. found) then
    cur_page = 0
    cur_entry = 0
  endif
  return
  end subroutine check_directory_entry
!!
! find entry called name in directory starting from beginning of directory (the hard way)
! upon exit cur_page and cur_entry are nonzero if desired entry found
! ==============================================================================================
   subroutine find_directory_entry(name, key)
   use gmm_internals
   implicit none
   character(len=*) :: name
   integer*8, optional :: key
   integer :: i
   character(len=GMM_MAXNAMELENGTH) :: l_name
!
   l_name = trim(name)
   cur_page = 1
   cur_entry = 1
   do i=1,used
     if (directory(cur_page)%entry(cur_entry)%name .eq. l_name) then
       if (present(key)) then
         key = directory(cur_page)%entry(cur_entry)%a%key
       endif
       return
     endif
     cur_entry = cur_entry + 1
     if (cur_entry .gt. PAGE_SIZE) then
       cur_page = cur_page + 1
       cur_entry = 1
     endif
   enddo
   cur_page = 0
   cur_entry = 0
!   if (present(key)) then
    key = GMM_KEY_NOT_FOUND
!   endif
   return
   end subroutine find_directory_entry
!
! locate/create a new properly initialized entry in directory
! ==============================================================================================
   subroutine add_directory_entry
   use gmm_internals
   implicit none
   integer :: i
!
   if ( table_size .eq. 0 ) then
     do i=1,MAX_PAGES
       nullify(directory(i)%entry)
     enddo
   endif
!
   used = used + 1
   last_entry = last_entry +1
   if ( last_entry .gt. PAGE_SIZE ) then
     table_size = table_size + 1
     last_entry = 1
     if (table_size .le. MAX_PAGES) then
       allocate(directory(table_size)%entry(PAGE_SIZE))
     else
!               print *,'ERROR: too many entries in directory for type=',EXTENSION
       call qqexit(1)
     endif
     do i=1,PAGE_SIZE
!               nullify( directory(table_size)%entry(i)%f )      ! invalid array pointer
         directory(table_size)%entry(i)%l = GMM_NULL_LAYOUT
         directory(table_size)%entry(i)%a = GMM_NULL_ATTRIB
       enddo
       cur_entry = 1
     else
       cur_entry = last_entry
     endif
     cur_page = table_size
     return
     end subroutine add_directory_entry
!
 function gmm_encodemeta(F_meta,F_code) result(F_istat)
    use gmm_internals
    implicit none
    type(gmm_metadata), intent(in) :: F_meta
    integer,            intent(out):: F_code(:)
    integer ::  F_istat
    integer :: i,j
    if (size(F_code) < GMM_META_SIZE) then
       F_istat = GMM_ERROR
       return
    endif
    F_istat = GMM_OK
    call movlev(F_meta, F_code, GMM_META_SIZE)
!     j=0
!     do i=1,4
!        j=j+1
!        F_code(j) = F_meta%l(i)%low
!        j=j+1
!        F_code(j) = F_meta%l(i)%high
!        j=j+1
!        F_code(j) = F_meta%l(i)%halo
!        j=j+1
!        F_code(j) = F_meta%l(i)%halomax
!        j=j+1
!        F_code(j) = F_meta%l(i)%n
!     enddo
!     j=j+1
!     call movlev(F_meta%a%key,F_code(j),2)
!     j=j+2
!     call movlev(F_meta%a%uuid1,F_code(j),2)
!     j=j+2
!     call movlev(F_meta%a%uuid2,F_code(j),2)
!     j=j+2
!     F_code(j) = F_meta%a%initmode
!     j=j+1
!     F_code(j) = F_meta%a%flags
    return
 end function gmm_encodemeta
 function gmm_decodemeta(F_meta,F_code) result(F_istat)
    use gmm_internals
    implicit none
    type(gmm_metadata), intent(out):: F_meta
    integer,            intent(in) :: F_code(:)
    integer ::  F_istat
    integer :: i,j
    if (size(F_code) < GMM_META_SIZE) then
       F_istat = GMM_ERROR
       return
    endif
    F_istat = GMM_OK
    call movlev(F_code, F_meta, GMM_META_SIZE)
!     j = 0
!     do i=1,4
!        j=j+1
!        F_meta%l(i)%low = F_code(j)
!        j=j+1
!        F_meta%l(i)%high = F_code(j)
!        j=j+1
!        F_meta%l(i)%halo = F_code(j)
!        j=j+1
!        F_meta%l(i)%halomax = F_code(j)
!        j=j+1
!        F_meta%l(i)%n = F_code(j)
!     enddo
!     j=j+1
!     call movlev(F_code(j),F_meta%a%key,2)
!     j=j+2
!     call movlev(F_code(j),F_meta%a%uuid1,2)
!     j=j+2
!     call movlev(F_code(j),F_meta%a%uuid2,2)
!     j=j+2
!     F_meta%a%initmode = F_code(j)
!     j=j+1
!     F_meta%a%flags = F_code(j)
    return
 end function gmm_decodemeta
!!===================== gmm_get (interface) =====================
!
  integer function gmm_get184(iname,p,m)
  use gmm_internals
  use pointer_table_data_184
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  integer*8, pointer  :: p(:,:,:,:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get184 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs184(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 4) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get184 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get184 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get184
!
  subroutine gmm_dealloc_ptr184()
  use gmm_internals
  use pointer_table_data_184
  implicit none
  deallocate (gmm_ptrs184(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs184(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr184
  integer function gmm_get144(iname,p,m)
  use gmm_internals
  use pointer_table_data_144
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  integer*4, pointer  :: p(:,:,:,:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get144 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs144(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 4) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get144 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get144 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get144
!
  subroutine gmm_dealloc_ptr144()
  use gmm_internals
  use pointer_table_data_144
  implicit none
  deallocate (gmm_ptrs144(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs144(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr144
  integer function gmm_get284(iname,p,m)
  use gmm_internals
  use pointer_table_data_284
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  real*8, pointer  :: p(:,:,:,:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get284 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs284(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 4) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get284 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get284 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get284
!
  subroutine gmm_dealloc_ptr284()
  use gmm_internals
  use pointer_table_data_284
  implicit none
  deallocate (gmm_ptrs284(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs284(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr284
  integer function gmm_get244(iname,p,m)
  use gmm_internals
  use pointer_table_data_244
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  real*4, pointer  :: p(:,:,:,:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get244 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs244(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 4) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get244 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get244 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get244
!
  subroutine gmm_dealloc_ptr244()
  use gmm_internals
  use pointer_table_data_244
  implicit none
  deallocate (gmm_ptrs244(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs244(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr244
  integer function gmm_get384(iname,p,m)
  use gmm_internals
  use pointer_table_data_384
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  complex*8, pointer  :: p(:,:,:,:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get384 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs384(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 4) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get384 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get384 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get384
!
  subroutine gmm_dealloc_ptr384()
  use gmm_internals
  use pointer_table_data_384
  implicit none
  deallocate (gmm_ptrs384(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs384(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr384
  integer function gmm_get183(iname,p,m)
  use gmm_internals
  use pointer_table_data_183
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  integer*8, pointer  :: p(:,:,:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get183 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs183(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 3) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get183 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get183 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get183
!
  subroutine gmm_dealloc_ptr183()
  use gmm_internals
  use pointer_table_data_183
  implicit none
  deallocate (gmm_ptrs183(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs183(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr183
  integer function gmm_get143(iname,p,m)
  use gmm_internals
  use pointer_table_data_143
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  integer*4, pointer  :: p(:,:,:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get143 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs143(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 3) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get143 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get143 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get143
!
  subroutine gmm_dealloc_ptr143()
  use gmm_internals
  use pointer_table_data_143
  implicit none
  deallocate (gmm_ptrs143(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs143(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr143
  integer function gmm_get283(iname,p,m)
  use gmm_internals
  use pointer_table_data_283
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  real*8, pointer  :: p(:,:,:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get283 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs283(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 3) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get283 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get283 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get283
!
  subroutine gmm_dealloc_ptr283()
  use gmm_internals
  use pointer_table_data_283
  implicit none
  deallocate (gmm_ptrs283(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs283(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr283
  integer function gmm_get243(iname,p,m)
  use gmm_internals
  use pointer_table_data_243
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  real*4, pointer  :: p(:,:,:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get243 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs243(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 3) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get243 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get243 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get243
!
  subroutine gmm_dealloc_ptr243()
  use gmm_internals
  use pointer_table_data_243
  implicit none
  deallocate (gmm_ptrs243(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs243(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr243
  integer function gmm_get383(iname,p,m)
  use gmm_internals
  use pointer_table_data_383
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  complex*8, pointer  :: p(:,:,:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get383 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs383(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 3) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get383 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get383 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get383
!
  subroutine gmm_dealloc_ptr383()
  use gmm_internals
  use pointer_table_data_383
  implicit none
  deallocate (gmm_ptrs383(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs383(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr383
  integer function gmm_get182(iname,p,m)
  use gmm_internals
  use pointer_table_data_182
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  integer*8, pointer  :: p(:,:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get182 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs182(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 2) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get182 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get182 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get182
!
  subroutine gmm_dealloc_ptr182()
  use gmm_internals
  use pointer_table_data_182
  implicit none
  deallocate (gmm_ptrs182(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs182(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr182
  integer function gmm_get142(iname,p,m)
  use gmm_internals
  use pointer_table_data_142
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  integer*4, pointer  :: p(:,:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get142 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs142(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 2) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get142 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get142 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get142
!
  subroutine gmm_dealloc_ptr142()
  use gmm_internals
  use pointer_table_data_142
  implicit none
  deallocate (gmm_ptrs142(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs142(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr142
  integer function gmm_get282(iname,p,m)
  use gmm_internals
  use pointer_table_data_282
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  real*8, pointer  :: p(:,:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get282 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs282(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 2) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get282 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get282 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get282
!
  subroutine gmm_dealloc_ptr282()
  use gmm_internals
  use pointer_table_data_282
  implicit none
  deallocate (gmm_ptrs282(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs282(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr282
  integer function gmm_get242(iname,p,m)
  use gmm_internals
  use pointer_table_data_242
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  real*4, pointer  :: p(:,:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get242 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs242(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 2) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get242 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get242 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get242
!
  subroutine gmm_dealloc_ptr242()
  use gmm_internals
  use pointer_table_data_242
  implicit none
  deallocate (gmm_ptrs242(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs242(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr242
  integer function gmm_get382(iname,p,m)
  use gmm_internals
  use pointer_table_data_382
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  complex*8, pointer  :: p(:,:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get382 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs382(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 2) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get382 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get382 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get382
!
  subroutine gmm_dealloc_ptr382()
  use gmm_internals
  use pointer_table_data_382
  implicit none
  deallocate (gmm_ptrs382(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs382(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr382
  integer function gmm_get181(iname,p,m)
  use gmm_internals
  use pointer_table_data_181
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  integer*8, pointer  :: p(:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get181 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs181(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 1) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get181 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get181 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get181
!
  subroutine gmm_dealloc_ptr181()
  use gmm_internals
  use pointer_table_data_181
  implicit none
  deallocate (gmm_ptrs181(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs181(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr181
  integer function gmm_get141(iname,p,m)
  use gmm_internals
  use pointer_table_data_141
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  integer*4, pointer  :: p(:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get141 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs141(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 1) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get141 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get141 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get141
!
  subroutine gmm_dealloc_ptr141()
  use gmm_internals
  use pointer_table_data_141
  implicit none
  deallocate (gmm_ptrs141(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs141(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr141
  integer function gmm_get281(iname,p,m)
  use gmm_internals
  use pointer_table_data_281
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  real*8, pointer  :: p(:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get281 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs281(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 1) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get281 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get281 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get281
!
  subroutine gmm_dealloc_ptr281()
  use gmm_internals
  use pointer_table_data_281
  implicit none
  deallocate (gmm_ptrs281(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs281(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr281
  integer function gmm_get241(iname,p,m)
  use gmm_internals
  use pointer_table_data_241
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  real*4, pointer  :: p(:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get241 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs241(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 1) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get241 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get241 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get241
!
  subroutine gmm_dealloc_ptr241()
  use gmm_internals
  use pointer_table_data_241
  implicit none
  deallocate (gmm_ptrs241(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs241(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr241
  integer function gmm_get381(iname,p,m)
  use gmm_internals
  use pointer_table_data_381
  implicit none
   integer :: i, array_rank
  character(len=*), intent(in) :: iname
  complex*8, pointer  :: p(:)
  type(gmm_metadata), optional, intent(out) :: m
!  integer,intent(inout) :: reqid
  include 'gmm_directory_interface.inc'
  type(gmm_metadata) :: m2
  integer*8 :: key
      integer *8 get_address_from
      external get_address_from
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    if (present(m)) then
      m%a = GMM_NULL_ATTRIB
      m%l = GMM_NULL_LAYOUT
    endif
    nullify(p)
    key= GMM_KEY_NOT_FOUND
    gmm_get381 = GMM_VAR_NOT_FOUND
  else
    m2%l=directory(cur_page)%entry(cur_entry)%l
    m2%a=directory(cur_page)%entry(cur_entry)%a
    if (present(m)) m=m2
    p=>gmm_ptrs381(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p
    do i=1,4
!      print *,'DEBUG gmm_get m%l(',i,')%n=',m2%l(i)%n
      if (m2%l(i)%n /= 0) array_rank=i
    enddo
!    write(6,'(a,a,a,i2,a,i2)') 'DEBUG gmm_get iname=',iname,' DIM=',DIM,' array_rank=',array_rank
    if (array_rank /= 1) then
       nullify(p)
       if (present(m)) m = GMM_NULL_METADATA
       gmm_get381 = GMM_INCONSISTENT_DIMS
!       print *,'DEBUG gmm_get *** GMM_INCONSISTENT_DIMS ***'
    else
       gmm_get381 = GMM_OK
    endif
!    write(6,'(a,a8,a,i4,a,i4,a,i4,a,i10)') 'Debug+++ gmm_get name=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' index=',directory(cur_page)%entry(cur_entry)%pointer_table_index,' addr=',get_address_from(p) 
   endif
  end function gmm_get381
!
  subroutine gmm_dealloc_ptr381()
  use gmm_internals
  use pointer_table_data_381
  implicit none
  deallocate (gmm_ptrs381(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  nullify    (gmm_ptrs381(directory(cur_page)%entry(cur_entry)%pointer_table_index)%p)
  end subroutine gmm_dealloc_ptr381
  integer function gmm_delete(iname)
  use gmm_internals
  implicit none
  character(len=*), intent(in) :: iname
  include 'gmm_directory_interface.inc'
  integer*8 :: key
  integer :: datatype
  key = 0
  call check_directory_entry(iname,key)
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    call find_directory_entry(iname,key)
  endif
  if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
    key= GMM_KEY_NOT_FOUND
    gmm_delete = GMM_VAR_NOT_FOUND
    return
  else
    datatype = directory(cur_page)%entry(cur_entry)%data_type
!   write(6,'(a,a,a,i2,a,i3,a,i4)') 'DEBUG gmm_delete iname=',iname,' cur_page=',cur_page,' cur_entry=',cur_entry,' datatype=',datatype
    directory(cur_page)%entry(cur_entry)%name = 'Variable deleted upon request'
    dtype: select case (datatype)
    case (184)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr184()
    case (144)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr144()
    case (284)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr284()
    case (244)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr244()
    case (384)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr384()
    case (183)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr183()
    case (143)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr143()
    case (283)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr283()
    case (243)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr243()
    case (383)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr383()
    case (182)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr182()
    case (142)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr142()
    case (282)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr282()
    case (242)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr242()
    case (382)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr382()
    case (181)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr181()
    case (141)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr141()
    case (281)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr281()
    case (241)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr241()
    case (381)
!      print *,'DEBUG gmm_delete appel a gmm_dealloc_ptr',datatype
       call gmm_dealloc_ptr381()
    end select dtype
    directory(cur_page)%entry(cur_entry)%l = GMM_NULL_LAYOUT
    directory(cur_page)%entry(cur_entry)%a = GMM_NULL_ATTRIB
    gmm_delete = GMM_OK
    return
  endif
  end function gmm_delete
!
  integer function gmm_getmeta2(iname,m)
      type gmm_layout
         SEQUENCE
         integer :: low,high,halo,halomax,n
      end type
      type gmm_attributes
        SEQUENCE
        integer*8 :: key
        integer*8 :: uuid1, uuid2
        integer   :: initmode
        integer   :: flags
      end type
      type gmm_metadata
        SEQUENCE
        type(gmm_layout), dimension(4) :: l
        type(gmm_attributes) :: a
      end type
    character(len=*), intent(in) :: iname
    type(gmm_metadata), intent(out) :: m
   integer gmm_getmeta
   external gmm_getmeta
   gmm_getmeta2 = gmm_getmeta(iname, m)
  end function gmm_getmeta2
  integer function gmm_nkeys()
  use gmm_internals
  implicit none
  gmm_nkeys = used
  return
  end function gmm_nkeys
  integer function gmm_keys(taglist,pattern)
  use gmm_internals
  implicit none
  character(len=*), intent(out) :: taglist(:)
  character(len=*), intent(in), optional :: pattern
  integer :: i,strlen_pattern,nkeys,maxkeys
  nkeys = 0
  maxkeys = nkeys
  gmm_keys = -1
!
  maxkeys = size(taglist)
  if (used > maxkeys) then
    return
  endif
  cur_page = 1
  cur_entry = 1
  if (present(pattern)) then
    strlen_pattern = len_trim(pattern)
    gmm_keys = 0
    do i=1,used
      if (directory(cur_page)%entry(cur_entry)%name(1:strlen_pattern) == pattern(1:strlen_pattern)) then
        taglist(nkeys+1) = directory(cur_page)%entry(cur_entry)%name
        nkeys = nkeys + 1
      endif
      cur_entry = cur_entry + 1
      if (cur_entry .gt. PAGE_SIZE) then
        cur_page = cur_page + 1
        cur_entry = 1
      endif
    enddo
    gmm_keys = nkeys
  else
    do i = 1,used
      taglist(i) = directory(cur_page)%entry(cur_entry)%name
      cur_entry = cur_entry + 1
      if (cur_entry .gt. PAGE_SIZE) then
        cur_page = cur_page + 1
        cur_entry = 1
      endif
    enddo
    gmm_keys = used
  endif
  return
  end  function gmm_keys
!
  integer function gmm_shuffle(taglist)
  use gmm_internals
  implicit none
  character(len=*), intent(out) :: taglist(:)
  include 'gmm_directory_interface.inc'
  integer :: i,strlen_pattern,nkeys,maxkeys, ier
  integer ind1, ind2, temp_entry, temp_page, te1, te2, tp1, tp2, temp_pti
  integer*8, dimension(:), allocatable :: key_list
  integer*8 :: temp_key
  type(p_gmm_metadata) temp
  character(len=GMM_MAXNAMELENGTH) :: tempname = '2121_Trans-Canada_Dorval_H9P-1J3'
  integer gmm_update_tpi_key2, gmm_rename
  external gmm_update_tpi_key2, gmm_rename
  logical ok, valide
  nkeys = 0
  maxkeys = nkeys
  gmm_shuffle = GMM_ERROR
!
! Epuration de la liste... On elimine les chaines de caracteres vides
  temp_key = GMM_KEY_NOT_FOUND
  call find_directory_entry(tempname,temp_key)
  if (temp_key /= GMM_KEY_NOT_FOUND) then
     print *, 'FATAL : (GMM_SHUFFLE) Temporary variable used for swapping should not exist'
     gmm_shuffle = GMM_ERROR
   endif
  allocate(key_list(size(taglist)))
  nkeys = size(taglist)
  i = 0
  ok = .true.
  do while (ok .and. i <= nkeys)
    i = i+1
    if (i <= nkeys) then
       if (0 == len_trim(taglist(i))) then
         ok = .false.
         i = i-1
       endif
    endif
  enddo
  if (i < nkeys) then
    nkeys = i
  endif
   do i=1,nkeys
      call find_directory_entry(taglist(i),key_list(i))
!      print *,  taglist(i),key_list(i)
   enddo
   valide = .false.
   do i=1,nkeys
      if (key_list(i) /= GMM_KEY_NOT_FOUND) then
         valide = .true.
         exit
      endif
   enddo
   if (.not.valide) then
      print *, '(GMM_SHUFFLE) NONE OF THE FIELDS IN THE LIST EXIST !'
      gmm_shuffle = GMM_ERROR
      return
   endif
   select case (nkeys)
   case (2)
      if (key_list(1) == GMM_KEY_NOT_FOUND) then
         ier = gmm_rename(taglist(2), taglist(1))
      elseif (key_list(2) == GMM_KEY_NOT_FOUND) then
         ier = gmm_rename(taglist(1), taglist(2))
      else
         ier = gmm_rename(taglist(1), tempname)
         ier = gmm_rename(taglist(2), taglist(1))
         ier = gmm_rename(tempname, taglist(2))
      endif
   case default
      temp_key = key_list(nkeys)
      do i=2,nkeys
         key_list(i) = key_list(i-1)
      enddo
      key_list(1) = temp_key
      if (key_list(nkeys) /= GMM_KEY_NOT_FOUND) then
         ier = gmm_rename(taglist(nkeys),tempname)
      endif
      do i=nkeys,2,-1
         if (key_list(i-1) /= GMM_KEY_NOT_FOUND) then
            ier = gmm_rename(taglist(i-1),taglist(i))
         endif
      enddo
      call find_directory_entry(tempname,temp_key)
      if (temp_key /= GMM_KEY_NOT_FOUND) then
         ier = gmm_rename(tempname, taglist(1))
      endif
  end select
   gmm_shuffle = 0
   return
   deallocate(key_list)
   end  function gmm_shuffle
   integer function gmm_rename(old_varname, new_varname)
   use gmm_internals
   implicit none
   character(len=*), intent(in) :: old_varname, new_varname
       interface
         subroutine check_directory_entry(name,key)
         character(len=*) :: name
         integer*8, intent(in) :: key
         end subroutine check_directory_entry
         subroutine find_directory_entry(name, key)
         implicit none
         character(len=*) :: name
         integer*8, optional :: key
         end subroutine find_directory_entry
         subroutine add_directory_entry
         end subroutine add_directory_entry
       end interface
   integer*8 :: key
   key = GMM_KEY_NOT_FOUND
   call find_directory_entry(new_varname,key)
   if (key >= 0) then
      if (gmm_verbose_level <= GMM_MSG_WARN) then
         print *, '(GMM_RENAME) Variable ', trim(new_varname), ' is already defined'
      endif
      gmm_rename = GMM_ERROR
      return
   endif
   key = GMM_KEY_NOT_FOUND
   call find_directory_entry(old_varname,key)
   if (key == GMM_KEY_NOT_FOUND) then
      if (gmm_verbose_level <= GMM_MSG_WARN) then
         print *, '(GMM_RENAME) Variable ', trim(old_varname), ' not defined'
      endif
      gmm_rename = GMM_ERROR
      return
   endif
   directory(cur_page)%entry(cur_entry)%name = new_varname
   if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *, '(GMM_RENAME) Variable ', trim(old_varname), ' renamed to ', trim(new_varname)
   endif
   gmm_rename = 0
   return
   end  function gmm_rename
   integer function gmm_getmeta(varname, meta)
   use gmm_internals
   implicit none
   character(len=*), intent(in) :: varname
   type(gmm_metadata), intent(out) :: meta
   include 'gmm_directory_interface.inc'
   integer   :: i
   integer*8 :: key
   key = 0
   call find_directory_entry(varname,key)
!   print *, varname, key
   if (key == GMM_KEY_NOT_FOUND) then
      if (gmm_verbose_level <= GMM_MSG_WARN) then
        print *, '(GMM_GETMETA) Variable ', varname, ' not found'
      endif
      gmm_getmeta = GMM_ERROR
      return
   endif
   meta%a = directory(cur_page)%entry(cur_entry)%a
   meta%l = directory(cur_page)%entry(cur_entry)%l
   gmm_getmeta = 0
!
   return
   end  function gmm_getmeta
   function gmm_updatemeta(iname, F_meta) result(F_istat)
   use gmm_internals
   implicit none
  include 'gmm_directory_interface.inc'
   character(len=*), intent(in) :: iname
   type(gmm_metadata), intent(in) :: F_meta
   integer ::  F_istat
   integer :: i,j
   integer*8 :: key
   key = 0
   call check_directory_entry(iname,key)
   if(cur_page .eq. 0 .or. cur_entry .eq. 0) then
      call find_directory_entry(iname,key)
   endif
   if (cur_page .eq. 0 .or. cur_entry .eq. 0) then
      F_istat = GMM_VAR_NOT_FOUND
      return
   endif
   do i=1,4
      directory(cur_page)%entry(cur_entry)%l(i)%halo = F_meta%l(i)%halo
      directory(cur_page)%entry(cur_entry)%l(i)%halomax = F_meta%l(i)%halomax
      directory(cur_page)%entry(cur_entry)%l(i)%n = F_meta%l(i)%n
   enddo
   directory(cur_page)%entry(cur_entry)%a%uuid1 = F_meta%a%uuid1
   directory(cur_page)%entry(cur_entry)%a%uuid2 = F_meta%a%uuid2
   directory(cur_page)%entry(cur_entry)%a%flags = F_meta%a%flags
   F_istat = GMM_OK
   return
end function gmm_updatemeta
!!===================== gmm_create (interface) =====================
!
!!===================== gmm_create (interface) =====================
!
  integer function gmm_update_tpi_key2(indx,datatype, key)
    implicit none
    integer, intent(in) :: indx, datatype
    integer*8, intent(in) :: key
    integer ier
      integer gmm_update_tpi_key184
      integer gmm_update_tpi_key144
      integer gmm_update_tpi_key284
      integer gmm_update_tpi_key244
      integer gmm_update_tpi_key384
      integer gmm_update_tpi_key183
      integer gmm_update_tpi_key143
      integer gmm_update_tpi_key283
      integer gmm_update_tpi_key243
      integer gmm_update_tpi_key383
      integer gmm_update_tpi_key182
      integer gmm_update_tpi_key142
      integer gmm_update_tpi_key282
      integer gmm_update_tpi_key242
      integer gmm_update_tpi_key382
      integer gmm_update_tpi_key181
      integer gmm_update_tpi_key141
      integer gmm_update_tpi_key281
      integer gmm_update_tpi_key241
      integer gmm_update_tpi_key381
    dtype: select case (datatype)
    case (184)
      gmm_update_tpi_key2 = gmm_update_tpi_key184(indx, key)
    case (144)
      gmm_update_tpi_key2 = gmm_update_tpi_key144(indx, key)
    case (284)
      gmm_update_tpi_key2 = gmm_update_tpi_key284(indx, key)
    case (244)
      gmm_update_tpi_key2 = gmm_update_tpi_key244(indx, key)
    case (384)
      gmm_update_tpi_key2 = gmm_update_tpi_key384(indx, key)
    case (183)
      gmm_update_tpi_key2 = gmm_update_tpi_key183(indx, key)
    case (143)
      gmm_update_tpi_key2 = gmm_update_tpi_key143(indx, key)
    case (283)
      gmm_update_tpi_key2 = gmm_update_tpi_key283(indx, key)
    case (243)
      gmm_update_tpi_key2 = gmm_update_tpi_key243(indx, key)
    case (383)
      gmm_update_tpi_key2 = gmm_update_tpi_key383(indx, key)
    case (182)
      gmm_update_tpi_key2 = gmm_update_tpi_key182(indx, key)
    case (142)
      gmm_update_tpi_key2 = gmm_update_tpi_key142(indx, key)
    case (282)
      gmm_update_tpi_key2 = gmm_update_tpi_key282(indx, key)
    case (242)
      gmm_update_tpi_key2 = gmm_update_tpi_key242(indx, key)
    case (382)
      gmm_update_tpi_key2 = gmm_update_tpi_key382(indx, key)
    case (181)
      gmm_update_tpi_key2 = gmm_update_tpi_key181(indx, key)
    case (141)
      gmm_update_tpi_key2 = gmm_update_tpi_key141(indx, key)
    case (281)
      gmm_update_tpi_key2 = gmm_update_tpi_key281(indx, key)
    case (241)
      gmm_update_tpi_key2 = gmm_update_tpi_key241(indx, key)
    case (381)
      gmm_update_tpi_key2 = gmm_update_tpi_key381(indx, key)
  end select dtype
  end function gmm_update_tpi_key2
!!===================== gmm_create (code) =====================
  integer function gmm_update_tpi_key184(indx, key)
  use gmm_internals
  use pointer_table_data_184
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 184, 'of size', gmm_p_table_size
      gmm_update_tpi_key184 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs184(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 184, 'of size', gmm_p_table_size
  gmm_update_tpi_key184 = 0
  end function gmm_update_tpi_key184
!
  integer function gmm_update_tpi_key144(indx, key)
  use gmm_internals
  use pointer_table_data_144
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 144, 'of size', gmm_p_table_size
      gmm_update_tpi_key144 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs144(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 144, 'of size', gmm_p_table_size
  gmm_update_tpi_key144 = 0
  end function gmm_update_tpi_key144
!
  integer function gmm_update_tpi_key284(indx, key)
  use gmm_internals
  use pointer_table_data_284
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 284, 'of size', gmm_p_table_size
      gmm_update_tpi_key284 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs284(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 284, 'of size', gmm_p_table_size
  gmm_update_tpi_key284 = 0
  end function gmm_update_tpi_key284
!
  integer function gmm_update_tpi_key244(indx, key)
  use gmm_internals
  use pointer_table_data_244
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 244, 'of size', gmm_p_table_size
      gmm_update_tpi_key244 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs244(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 244, 'of size', gmm_p_table_size
  gmm_update_tpi_key244 = 0
  end function gmm_update_tpi_key244
!
  integer function gmm_update_tpi_key384(indx, key)
  use gmm_internals
  use pointer_table_data_384
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 384, 'of size', gmm_p_table_size
      gmm_update_tpi_key384 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs384(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 384, 'of size', gmm_p_table_size
  gmm_update_tpi_key384 = 0
  end function gmm_update_tpi_key384
!
  integer function gmm_update_tpi_key183(indx, key)
  use gmm_internals
  use pointer_table_data_183
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 183, 'of size', gmm_p_table_size
      gmm_update_tpi_key183 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs183(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 183, 'of size', gmm_p_table_size
  gmm_update_tpi_key183 = 0
  end function gmm_update_tpi_key183
!
  integer function gmm_update_tpi_key143(indx, key)
  use gmm_internals
  use pointer_table_data_143
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 143, 'of size', gmm_p_table_size
      gmm_update_tpi_key143 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs143(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 143, 'of size', gmm_p_table_size
  gmm_update_tpi_key143 = 0
  end function gmm_update_tpi_key143
!
  integer function gmm_update_tpi_key283(indx, key)
  use gmm_internals
  use pointer_table_data_283
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 283, 'of size', gmm_p_table_size
      gmm_update_tpi_key283 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs283(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 283, 'of size', gmm_p_table_size
  gmm_update_tpi_key283 = 0
  end function gmm_update_tpi_key283
!
  integer function gmm_update_tpi_key243(indx, key)
  use gmm_internals
  use pointer_table_data_243
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 243, 'of size', gmm_p_table_size
      gmm_update_tpi_key243 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs243(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 243, 'of size', gmm_p_table_size
  gmm_update_tpi_key243 = 0
  end function gmm_update_tpi_key243
!
  integer function gmm_update_tpi_key383(indx, key)
  use gmm_internals
  use pointer_table_data_383
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 383, 'of size', gmm_p_table_size
      gmm_update_tpi_key383 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs383(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 383, 'of size', gmm_p_table_size
  gmm_update_tpi_key383 = 0
  end function gmm_update_tpi_key383
!
  integer function gmm_update_tpi_key182(indx, key)
  use gmm_internals
  use pointer_table_data_182
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 182, 'of size', gmm_p_table_size
      gmm_update_tpi_key182 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs182(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 182, 'of size', gmm_p_table_size
  gmm_update_tpi_key182 = 0
  end function gmm_update_tpi_key182
!
  integer function gmm_update_tpi_key142(indx, key)
  use gmm_internals
  use pointer_table_data_142
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 142, 'of size', gmm_p_table_size
      gmm_update_tpi_key142 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs142(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 142, 'of size', gmm_p_table_size
  gmm_update_tpi_key142 = 0
  end function gmm_update_tpi_key142
!
  integer function gmm_update_tpi_key282(indx, key)
  use gmm_internals
  use pointer_table_data_282
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 282, 'of size', gmm_p_table_size
      gmm_update_tpi_key282 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs282(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 282, 'of size', gmm_p_table_size
  gmm_update_tpi_key282 = 0
  end function gmm_update_tpi_key282
!
  integer function gmm_update_tpi_key242(indx, key)
  use gmm_internals
  use pointer_table_data_242
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 242, 'of size', gmm_p_table_size
      gmm_update_tpi_key242 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs242(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 242, 'of size', gmm_p_table_size
  gmm_update_tpi_key242 = 0
  end function gmm_update_tpi_key242
!
  integer function gmm_update_tpi_key382(indx, key)
  use gmm_internals
  use pointer_table_data_382
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 382, 'of size', gmm_p_table_size
      gmm_update_tpi_key382 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs382(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 382, 'of size', gmm_p_table_size
  gmm_update_tpi_key382 = 0
  end function gmm_update_tpi_key382
!
  integer function gmm_update_tpi_key181(indx, key)
  use gmm_internals
  use pointer_table_data_181
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 181, 'of size', gmm_p_table_size
      gmm_update_tpi_key181 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs181(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 181, 'of size', gmm_p_table_size
  gmm_update_tpi_key181 = 0
  end function gmm_update_tpi_key181
!
  integer function gmm_update_tpi_key141(indx, key)
  use gmm_internals
  use pointer_table_data_141
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 141, 'of size', gmm_p_table_size
      gmm_update_tpi_key141 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs141(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 141, 'of size', gmm_p_table_size
  gmm_update_tpi_key141 = 0
  end function gmm_update_tpi_key141
!
  integer function gmm_update_tpi_key281(indx, key)
  use gmm_internals
  use pointer_table_data_281
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 281, 'of size', gmm_p_table_size
      gmm_update_tpi_key281 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs281(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 281, 'of size', gmm_p_table_size
  gmm_update_tpi_key281 = 0
  end function gmm_update_tpi_key281
!
  integer function gmm_update_tpi_key241(indx, key)
  use gmm_internals
  use pointer_table_data_241
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 241, 'of size', gmm_p_table_size
      gmm_update_tpi_key241 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs241(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 241, 'of size', gmm_p_table_size
  gmm_update_tpi_key241 = 0
  end function gmm_update_tpi_key241
!
  integer function gmm_update_tpi_key381(indx, key)
  use gmm_internals
  use pointer_table_data_381
  implicit none
    integer, intent(in) :: indx
    integer*8, intent(in) :: key
  type(gmm_attributes) :: localattr, attrs
  type (gmm_attributes) lcl_attr
  type(gmm_layout), dimension(4) :: lcl_layout, dims
  integer lcl_datatype
  logical consistent
  integer i, ier
  integer lcl_pti
    if (indx > gmm_p_table_size) then
      print *, 'update_table_entry : wrong index', indx,'in table type ', 381, 'of size', gmm_p_table_size
      gmm_update_tpi_key381 = GMM_POINTER_TABLE_OVERFLOW
    endif
    gmm_ptrs381(indx)%key = key
    print *, 'update_table_entry', 'of', indx, 'in table type ', 381, 'of size', gmm_p_table_size
  gmm_update_tpi_key381 = 0
  end function gmm_update_tpi_key381
!
  integer function gmm_verbosity(verbose_level)
  use gmm_internals
  implicit none
  integer, intent(in) :: verbose_level
   select case (verbose_level)
   case (GMM_MSG_DEBUG)
      gmm_verbose_level = GMM_MSG_DEBUG
      print *, '(GMM_VERBOSITY) Setting GMM Verbose level to DEBUG'
   case (GMM_MSG_INFO)
      gmm_verbose_level = GMM_MSG_INFO
      print *, '(GMM_VERBOSITY) Setting GMM Verbose level to INFO'
   case (GMM_MSG_WARN)
      gmm_verbose_level = GMM_MSG_WARN
      print *, '(GMM_VERBOSITY) Setting GMM Verbose level to WARN'
   case (GMM_MSG_ERROR)
      gmm_verbose_level = GMM_MSG_ERROR
      print *, '(GMM_VERBOSITY) Setting GMM Verbose level to ERROR'
   case (GMM_MSG_SEVERE)
      gmm_verbose_level = GMM_MSG_SEVERE
      print *, '(GMM_VERBOSITY) Setting GMM Verbose level to SEVERE'
   case (GMM_MSG_FATAL)
      gmm_verbose_level = GMM_MSG_FATAL
      print *, '(GMM_VERBOSITY) Setting GMM Verbose level to FATAL'
   case default
      print *, '(GMM_VERBOSITY) Unknown GMM Verbose level'
      print *, '(GMM_VERBOSITY) Please check GMM Documentation'
   end select
  gmm_verbosity = 0
!
  return
  end function gmm_verbosity
    subroutine gmm_dumpinfo(fldstat)
    use gmm_internals
    implicit none
    logical, intent(in), optional :: fldstat
    integer :: i, l_page, l_entry, nelm, crc, f_calc_crc
    real :: xx
    pointer(px, xx(*))
    type(gmm_layout), dimension(4) :: dims
    external f_calc_crc
    l_page = 1
    l_entry = 1
    print *,'GMM dumpinfo, number of variables in use is',used
    do i=1,used
      dims = directory(l_page)%entry(l_entry)%l
      nelm = ( (dims(1)%high - dims(1)%low +1) * &
              &(dims(2)%high - dims(2)%low +1) * &
              &(dims(3)%high - dims(3)%low +1) * &
              &(dims(4)%high - dims(4)%low +1) )
      if (present(fldstat)) then
        print *,'Appel a statfld a ecrire, fldstat=',fldstat
        write (6,77) 'Name=',directory(l_page)%entry(l_entry)%name,' addr=',&
                      &directory(l_page)%entry(l_entry)%array_addr
      else
        call make_cray_pointer(px,directory(l_page)%entry(l_entry)%array_addr)
        crc = f_calc_crc(xx,nelm,0,1)
        write (6,88) 'Name=',directory(l_page)%entry(l_entry)%name,' addr=',&
                      &directory(l_page)%entry(l_entry)%array_addr,' checksum=',crc
      endif
      l_entry = l_entry + 1
      if (l_entry .gt. PAGE_SIZE) then
        l_page = l_page + 1
        l_entry = 1
      endif
    enddo
77  format(a,a,a,i10)
88  format(a,a,a,i10,a,i10)
    return
    end subroutine gmm_dumpinfo
