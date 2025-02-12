module phymem
   use clib_itf_mod, only: clib_toupper, clib_tolower
   use rmn_gmm, only: gmm_metadata, gmm_create, gmm_get, GMM_IS_OK, GMM_FLAG_RSTR, GMM_FLAG_IZER, GMM_FLAG_INAN, GMM_NULL_FLAGS, GMM_MAXNAMELENGTH
   use str_mod, only: str_normalize, str_concat
   use phy_status, only: PHY_OK, PHY_ERROR
   use splitst, only: splitst4
   use phygridmap, only: phy_lcl_ni, phy_lcl_nj, phydim_ni, phydim_nj, phydim_nk
   use phy_options, only: debug_mem_L
   use debug_mod, only: init2nan
   implicit none
   private

#include <rmnlib_basics.hf>
#include <rmn/msg.h>
   
   !# Public functions
   public :: phymem_busidx, phymem_alloc, phymem_init
   public :: phymem_add, phymem_gmmname, phymem_busreset
   public :: phymem_find, phymem_find_idxv, phymem_find_idxv1, phymem_updatemeta
   public :: phymem_getmeta, phymem_getmeta_copy, phymem_getdata
   public :: phymem_get_slabvars, phymem_get_i0_string
   
   !# Public constants
   integer, parameter, public :: PHY_NBUSES = 4
   integer, parameter, public :: PHY_NAMELEN  = 32
   integer, parameter, public :: PHY_DESCLEN  = 256
   integer, parameter, public :: PHY_MAXFLAGS = 16
   integer, parameter, public :: PHY_MAXDEPS  = 16

   integer, parameter, public :: PHY_STAG_SFC = 0
   integer, parameter, public :: PHY_STAG_MOM = 0
   integer, parameter, public :: PHY_STAG_THERMO = 1
   integer, parameter, public :: PHY_STAG_ENERGY = 2

   character(len=*), parameter, public :: PHY_NPATH_DEFAULT='VOI'
   character(len=*), parameter, public :: PHY_BPATH_DEFAULT='VPDE'

   !# Bus index in pbuses(:)
   integer, parameter, public :: PHY_DBUSIDX = 1
   integer, parameter, public :: PHY_PBUSIDX = 2
   integer, parameter, public :: PHY_VBUSIDX = 3
   integer, parameter, public :: PHY_UBUSIDX = 3
   integer, parameter, public :: PHY_EBUSIDX = 4
   character(len=1), parameter, public :: PHY_BUSID(PHY_NBUSES) = (/'D', 'P', 'V', 'E'/)
   integer, parameter, public :: gmmflags(PHY_NBUSES) = (/ 0, GMM_FLAG_RSTR, 0, 0/)
   !# Bus Slab index in pvars
   integer, save, public :: PHY_BUSIDXV(PHY_NBUSES) = (/ -1, -1, -1, -1/)

   !# Public vars   
   integer, save, public :: nphyvars = 0  !# Total number of phy vars (size of pvmetas)
   integer, save, public :: PHY_MAXVARS(PHY_NBUSES) = (/0, 0, 0, 0/)

   !# Public Derived types -------------------------------------------------
   public :: phymeta, phyvar

   type phymeta
      integer :: ni      !# folded ni dim (p_runlenght)
      integer :: nk      !# number of atmospheric levels.
      integer :: fmul    !# number of arbitrarily-defined levels
      integer :: mosaic  !# number of surface sub-types
      integer :: size    !# ni * nk * fmul * (mosaic+1)
      integer :: nlcl(3) !# local tile dims in "not folded" space
      integer :: ibus    !# index of bus containing the field
      integer :: idxb    !# var index in the specified bus, pbuses(ibus)%meta(idxb)
      integer :: idxv    !# var index in the specified bus, pvmetas(idxv)
      integer :: i0      !# index of first element in the bus pointer, pbuses(ibus)%bptr(i0:in,:)
      integer :: in      !# in=i0+size-1; index of first element in the bus pointer, pbuses(ibus)%bptr(i0:in,:)      
      integer :: init    !# 1 = init/mandatory, 0 otherwise
      integer :: reset   !# 1 = VBUS init/reset at stepinit, 0 = no reset
      integer :: stag    !# Level staggering information: PHY_STAG_MOM/THERMO/ENERGY
      logical :: wload   !# tracer attribute: T = use in the density calculation (water loaded), 0 otherwise
      logical :: hzd     !# tracer attribute: T = perform horizontal diffusion, 0 otherwise
      integer :: monot   !# tracer attribute: 1 = use monoton interpolation in adv, 0 otherwise
      integer :: massc   !# tracer attribute: 1 = use mass conservation in adv, 0 otherwise
      real    :: vmin    !# applied minvalue for the field
      real    :: vmax    !# applied maxvalue for the field

      character(len=PHY_NAMELEN) :: bus    !# name of the bus containing the field 
      character(len=PHY_NAMELEN) :: vname  !# var name (in code name)
      character(len=PHY_NAMELEN) :: iname  !# input name
      character(len=PHY_NAMELEN) :: oname  !# output name
      character(len=PHY_NAMELEN) :: sname  !# series name
      character(len=PHY_DESCLEN) :: desc   !# Field description
      character(len=PHY_NAMELEN) :: flags(PHY_MAXFLAGS)  !# list of text key describing the fields caracteristics
      
      character(len=PHY_NAMELEN) :: deps(PHY_MAXDEPS)    !# list of others vars (VNAME) needed to compute this one
      integer :: ideps(PHY_MAXDEPS)
      logical :: outcond_L
      logical :: outreq_L
   end type phymeta
   
   type phyvar
      type(phymeta), pointer :: meta => null()
      real, pointer, contiguous :: data(:) => null()
   end type phyvar

   
   !# Private Devived Type -------------------------------------------------
   
   type phymetaptr
      type(phymeta), pointer :: meta => null()
   end type phymetaptr
 
   type phybus
      integer :: nvars = 0
      integer :: nsize = 0
      type(phymetaptr), pointer, contiguous :: mptr(:) => null() !# ptr to var metaptr (pbuses(ibus)%mptr(idxb)%meta => pvmetas(idxv)
      real, pointer, contiguous :: bptr(:,:) => null()  !# ptr to bus, all slabs
      character(len=32) :: busid = ' '                  !# Bus short name
   end type phybus

   
   !# Interfaces Overloading -----------------------------------------------
   
   interface phymem_add
      module procedure phymem_add_meta
      module procedure phymem_add_string
   end interface phymem_add
   
   interface phymem_busreset
      module procedure phymem_busreset1
      module procedure phymem_busreset2
   end interface phymem_busreset
   
   interface phymem_find
      module procedure phymem_find_idxv
      module procedure phymem_find_idxv1
   end interface phymem_find

   interface phymem_getdata
      module procedure phymem_getdata1d
      module procedure phymem_getdata2d
      module procedure phymem_getdata3d
   end interface phymem_getdata
   
   !# Private module vars --------------------------------------------------
   
   integer, parameter :: DNVARS = 25

   type(phybus), save  :: pbuses(PHY_NBUSES)
   type(phymeta), pointer, contiguous, save :: pvmetas(:) => null()

   logical, save :: isallocated = .false.
   logical, save :: isinit = .false.

   integer, save, pointer :: uidxvlist(:) => NULL()
   integer, save :: nuidxvlist = -1
   integer, save :: nphyvarsmax = 0

contains

   !/@*
   function phymem_init(F_maxvars) result(F_istat)
      implicit none
      !@objective init phymem module and set maxvars
      !@arguments
      integer, intent(in), optional :: F_maxvars
      !@return
      integer :: F_istat
      !*@/
      !---------------------------------------------------------------
      if (present(F_maxvars)) then
         F_istat = priv_init(max(1, F_maxvars))
      else
         F_istat = priv_init()
      endif
      !---------------------------------------------------------------
      return
   end function phymem_init

   
   !/@*
   function phymem_busidx(F_busname, F_quiet) result(F_busidx)
      implicit none
      !@objective Convert memory bus name (string) to index
      !@arguments
      character(len=*), intent(in) :: F_busname
      logical, optional, intent(in) :: F_quiet
      !@return
      integer :: F_busidx
      !*@/
      integer :: istat
      logical :: quiet
      character(len=1) :: busname
      !---------------------------------------------------------------
      busname = F_busname
      quiet = .false.
      if (present(F_quiet)) quiet = F_quiet
      istat = clib_toupper(busname)
!!$      if (busname == 'U') busname = 'V'
      do F_busidx = 1, PHY_NBUSES
         if (busname == PHY_BUSID(F_busidx)) return
      enddo
      if (.not.quiet) call msg(MSG_WARNING, &
           '(phymem_busidx) No such bus name: '//trim(F_busname))
      F_busidx = -1
      !---------------------------------------------------------------
      return
   end function phymem_busidx

   
   !/@*
   function phymem_gmmname(F_busidx) result(F_gmmname)
      implicit none
      !@objective Convert memory bus index to name
      !@arguments
      integer, intent(in) :: F_busidx
      !@return
      character(len=32) :: F_gmmname
      !*@/
      !---------------------------------------------------------------
      F_gmmname = ''
      if (F_busidx > 0 .and. F_busidx <= PHY_NBUSES) &
           F_gmmname = trim(PHY_BUSID(F_busidx))//'BUS_3d'
      !---------------------------------------------------------------
      return
   end function phymem_gmmname
   

   !/@*
   function phymem_busreset1(F_busidx, F_value, F_trnch) result(F_istat)
      implicit none
      integer, intent(in) :: F_busidx
      real, intent(in), optional :: F_value
      integer, intent(in), optional :: F_trnch
      !@objective Reset whole memory bus to specified value (default zero) on specified OMP silce (trnch, default: all slice)
      !@arguments
      !@return
      integer :: F_istat
      !*@/
      integer :: istat
      real :: rvalue
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (F_busidx <1 .and. F_busidx > PHY_NBUSES) then
         call msg(MSG_ERROR, '(phymem) busreset F_busidx out of range')
         return
      endif
      if (.not.isallocated) then
         call msg(MSG_ERROR, '(phymem) busreset - must call phymem_alloc first')
         return
      endif
      rvalue = 0.
      if (present(F_value)) rvalue = F_value
      if (present(F_trnch)) then
         if (F_trnch < 1 .or. F_trnch > phydim_nj) then
            call msg(MSG_ERROR,'(phymem_busreset) Requested F_trnch out of range')
            return
         endif
         pbuses(F_busidx)%bptr(:,F_trnch) = rvalue
      else
         pbuses(F_busidx)%bptr = rvalue
      endif
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_busreset1

   
   !/@*
   function phymem_busreset2(F_pvars, F_busidx, F_value, F_resetonly_L) result(F_istat)
      implicit none
      !@objective Reset whole memory bus to specified value (default zero) for provided memory slice structure (pvars)
      !@arguments
      type(phyvar), pointer, contiguous :: F_pvars(:)
      integer, intent(in) :: F_busidx
      real, intent(in), optional :: F_value
      logical, intent(in), optional :: F_resetonly_L
      !@return
      integer :: F_istat
      !*@/
      real :: rvalue
      logical :: resetonly_L
      integer :: n  !#, ntot
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (F_busidx <1 .and. F_busidx > PHY_NBUSES) then
         call msg(MSG_ERROR, '(phymem) busreset F_busidx out of range')
         return
      endif
      if (.not.isallocated) then
         call msg(MSG_ERROR, '(phymem) busreset - must call phymem_alloc first')
         return
      endif
      rvalue = 0.
      if (present(F_value)) rvalue = F_value
      resetonly_L = .false.
      if (present(F_resetonly_L)) resetonly_L = F_resetonly_L
      if (resetonly_L .and. F_busidx==PHY_VBUSIDX) then
         if (debug_mem_L) call init2nan(F_pvars(PHY_BUSIDXV(F_busidx))%data)
         do n=1,nuidxvlist
           F_pvars(uidxvlist(n))%data = rvalue
         enddo
      else
         F_pvars(PHY_BUSIDXV(F_busidx))%data = rvalue
      endif
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_busreset2

         
   !/@*
   function phymem_add_meta(F_imeta) result(F_istat)
      !@objective  Add physics var from pre-filler phymeta
      !@arguments
      type(phymeta), intent(in) :: F_imeta
      !@return
      integer :: F_istat
      !*@/
      type(phymeta) :: vmeta
      integer :: istat, nkfm, idxb, idxv, ibus, memgap, vsize, n
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (.not.isinit) istat = priv_init()
      if (isallocated) then
         call msg(MSG_ERROR, '(phymem_addvar) Cannot add var after allocation')
         return
      endif

      ibus = phymem_busidx(F_imeta%bus)
      if (ibus < 1 .or. ibus > PHY_NBUSES) then
         call msg(MSG_ERROR, &
              '(phymem_add) Unknown bus: '//trim(F_imeta%bus))
         return
      endif
      if (nphyvars >= nphyvarsmax) istat = priv_resize(DNVARS)
      
      vmeta = F_imeta
      vmeta%ibus = ibus
      call str_normalize(vmeta%iname)
      call str_normalize(vmeta%oname)
      call str_normalize(vmeta%vname)
      call str_normalize(vmeta%sname)
!!$      istat = clib_toupper(vmeta%iname)
!!$      istat = clib_toupper(vmeta%oname)
!!$      istat = clib_toupper(vmeta%vname)
!!$      istat = clib_toupper(vmeta%sname)
      istat = clib_toupper(vmeta%bus)

      do n=1,PHY_MAXDEPS
         call str_normalize(vmeta%deps(n))
         vmeta%ideps(n) = -1
      enddo
      vmeta%outreq_L = .false.

      if (vmeta%oname == '') vmeta%oname = vmeta%vname
      if (vmeta%iname == '') vmeta%iname = vmeta%oname
      if (vmeta%sname == '') vmeta%sname = vmeta%oname
      
      vmeta%ni = phydim_ni

      istat = priv_checkvar(vmeta)
      if (istat /= PHY_OK) return

      idxv = nphyvars + 1
      vmeta%idxv = idxv

      nkfm = vmeta%nk*vmeta%fmul*(vmeta%mosaic+1)
      vsize = vmeta%ni * nkfm

      vmeta%nlcl = (/phy_lcl_ni, phy_lcl_nj, nkfm/)
      vmeta%size = vsize
      
      pvmetas(idxv) = vmeta
      nphyvars = idxv
      
      ! if (F_imeta%flags(1) /= '') then
      !    print *,'(phymem_add_meta) ',trim(pbuses(ibus)%meta(idxb)%vname),':',trim(pbuses(ibus)%meta(idxb)%flags(1)),':',trim(pbuses(ibus)%meta(idxb)%flags(2))
      ! endif
      
      F_istat = RMN_OK
      !---------------------------------------------------------------
      return
   end function phymem_add_meta

   
   !/@*
   function phymem_add_string(F_string, F_flags, F_deps) result(F_istat)
      !@objective Add physics var from string description
      !@arguments
      character(len=*), intent(in) :: F_string
      character(len=*), intent(in), optional :: F_flags
      character(len=*), intent(in), optional :: F_deps(:)
      !@return
      integer :: F_istat
      !@Notes
      !  See splitst.F90 for the input string description
      !*@/
      integer, parameter :: FSTNAMELEN = 4
      integer :: istat, n, nflags
      character(len=3) ::   shape
      type(phymeta) :: vmeta
      character(len=PHY_NAMELEN) :: flags(PHY_MAXFLAGS)
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (.not.isinit) istat = priv_init()
      if (isallocated) then
         call msg(MSG_ERROR, '(phymem_addvar) Cannot add var after allocation')
         return
      endif

      istat = splitst4(vmeta%vname, vmeta%oname, vmeta%iname, vmeta%sname, &
           vmeta%desc, shape, vmeta%mosaic, vmeta%fmul, &
           vmeta%bus, vmeta%init, vmeta%vmin, vmeta%vmax, &
           vmeta%wload, vmeta%hzd, vmeta%monot, vmeta%massc, vmeta%flags, &
           F_string)
      if (.not.RMN_IS_OK(istat)) then
         call physeterror('phymem_add', 'Invalid string: '//trim(F_string))
         return
      endif

      if (any(shape == (/'M', 'T', 'E'/)) .and. vmeta%fmul > 1 .and. &
           max(len_trim(vmeta%oname),len_trim(vmeta%iname)) > FSTNAMELEN-1) then
         call msg(MSG_WARNING, ' (phymem_add) Varname '//trim(vmeta%vname)// &
              ' is declared as a 3D var with multiple categories, fmul>1,'// &
              ' I/O names must be 2/3 char max.')
      endif
      
      select case(shape)
      case("A")  !# arbitrary
         vmeta%stag = PHY_STAG_SFC
      case("M")  !# momentum
         vmeta%stag = PHY_STAG_MOM
      case("T")  !# thermo
         vmeta%stag = PHY_STAG_THERMO
      case("E")  !# energy
         vmeta%stag = PHY_STAG_ENERGY
      case default
         call physeterror('phymem_add', 'VS=(SHAPE) NOT ALLOWED: '//trim(F_string))
         return
      end select

      vmeta%reset = 0
      if (vmeta%bus == 'U') then
         vmeta%bus = 'V'
      elseif (vmeta%bus == 'V') then
         vmeta%reset = 1
      endif
      if (.not.any(vmeta%bus == PHY_BUSID(:)))  then
         call physeterror('phymem_add', 'VB=(BUS) NOT ALLOWED: '//trim(F_string))
         return
      endif

      if (present(F_flags)) then
         nflags = 0
         do n=1,size(vmeta%flags)
            if (vmeta%flags(n) == '') exit
            nflags = nflags+1
         enddo
         call str_split2list(flags, F_flags, '+', size(flags))
         do n=1,size(flags)
            if (flags(n) == '') cycle
            call str_normalize(flags(n))
            istat = clib_toupper(flags(n))
            if (any(vmeta%flags == flags(n))) cycle
            nflags = nflags+1
            if (nflags > size(vmeta%flags)) then
               call physeterror('phymem_add', 'too many flags for: '//trim(F_string))
               return
            endif
            vmeta%flags(nflags) = flags(n)
         enddo
      endif

      vmeta%deps(:) = ''
      if (present(F_deps)) then
         if (size(F_deps) > PHY_MAXDEPS) then
            call physeterror('phymem_add', 'too many deps for: '//trim(F_string))
            return
         endif
         vmeta%deps(1:size(F_deps)) = F_deps(:)
      endif
      vmeta%outcond_L = any(vmeta%flags == 'DIAG')

      ! if (vmeta%flags(1) /= '') then
      !    print *,'(phymem_add_string) ',trim(vmeta%vname),':',trim(vmeta%flags(1)),':',trim(vmeta%flags(2))
      ! endif
      
      vmeta%nk = phydim_nk
      if (shape == "A") vmeta%nk = 1

      F_istat = phymem_add_meta(vmeta)
      !---------------------------------------------------------------
      return
   end function phymem_add_string

   
   !/@*
   function phymem_get_i0_string(F_string) result(F_i0)
      !@objective Get i0 (first point in bus) for var as described in string (for legacy code)
      !@arguments
      character(len=*), intent(in) :: F_string
      !@return
      integer :: F_i0  !# var index in pvmetas
      !@Notes
      !  See splitst.F90 for the input string description
      !*@/
      integer, parameter :: FSTNAMELEN = 4
      integer :: istat, idxvlist(1)
      character(len=3) ::   shape
      type(phymeta) :: vmeta
      !---------------------------------------------------------------
      F_i0 = -1
      if (.not.isallocated) then
         call msg(MSG_ERROR,'(phymem_get_i0_string) phymem_alloc must be called before')
         return
      endif

      istat = splitst4(vmeta%vname, vmeta%oname, vmeta%iname, vmeta%sname, &
           vmeta%desc, shape, vmeta%mosaic, vmeta%fmul, &
           vmeta%bus, vmeta%init, vmeta%vmin, vmeta%vmax, &
           vmeta%wload, vmeta%hzd, vmeta%monot, vmeta%massc, vmeta%flags, &
           F_string)
      if (.not.RMN_IS_OK(istat)) then
         call physeterror('phymem_get_idxv_string', 'Invalid string: '//trim(F_string))
         return
      endif
!!$      if (vmeta%bus == 'U') vmeta%bus = 'V'
      istat = phymem_find(idxvlist, F_name=vmeta%vname, F_npath='V', F_bpath=vmeta%bus)
      if (istat <= 0) then
         call physeterror('phymem_get_idxv_string', 'Cannot find var for: '//trim(F_string))
         return
      endif
      F_i0 = pvmetas(idxvlist(1))%i0
      !---------------------------------------------------------------
      return
   end function phymem_get_i0_string

   
   !/@*
   function phymem_alloc(F_debug, F_outreq_S) result(F_istat)
      implicit none
      !@objective Allocate memory buses for previously added vars
      !@arguments
      logical, intent(in) :: F_debug
      character(len=*), intent(in) :: F_outreq_S(:)
      !@return
      integer :: F_istat
      !*@/
      integer :: ib, iv, idxv, halo, init, i0, in, istat, memgap, n, nphyvars2
      character(len=GMM_MAXNAMELENGTH) :: gmmname, vname, oname
      type(gmm_metadata) :: gmmmeta
      logical :: changed_L, isreq_L
      character(len=GMM_MAXNAMELENGTH) :: msg_S
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (isallocated) then
         call msg(MSG_ERROR, '(phymem) called phymem_alloc twice')
         return
      endif
      if (.not.isinit) istat = priv_init()
      
      if (nphyvars <= 0) then
         call msg(MSG_ERROR, '(phymem) called phymem_alloc w/o any var defined')
         return
      endif

      !# Set outreq_L for all outcond_L vars
      call str_concat(msg_S,F_outreq_S,', ')
      call msg(MSG_INFOPLUS, '(phymem_alloc) outreq: '//msg_S)
      changed_L = .true.
      DO_CHANGED: do while (changed_L)
         changed_L = .false.
         DO_IV: do iv = 1, nphyvars
            if (pvmetas(iv)%outcond_L) then
               vname = pvmetas(iv)%vname ; istat = clib_toupper(vname)  !#TODO: better use all upper *pvmetas%*name... cause problem... check it out
               oname = pvmetas(iv)%oname ; istat = clib_toupper(oname)
               isreq_L = (any(F_outreq_S == vname) .or. any(F_outreq_S == oname) &
                    .or. F_outreq_S(1) == '*')
               if (isreq_L .and. .not.pvmetas(iv)%outreq_L) changed_L = .true.
               pvmetas(iv)%outreq_L = (pvmetas(iv)%outreq_L .or. isreq_L)
               if (pvmetas(iv)%outreq_L) then
                  DO_IDEP: do n=1,PHY_MAXDEPS
                     if (pvmetas(iv)%deps(n) /= '') then
                        if (pvmetas(iv)%ideps(n) <= 0) then
                           pvmetas(iv)%ideps(n) = priv_find(pvmetas(iv)%deps(n))
                           if (pvmetas(iv)%ideps(n) <= 0) then
                              call msg(MSG_ERROR, '(phymem) phymem_alloc - cannot output requested var '//trim(pvmetas(iv)%vname)//' -- missing dep: '//trim(pvmetas(iv)%deps(n)))
                              return
                           endif
                        endif
                        if (pvmetas(pvmetas(iv)%ideps(n))%outcond_L .and. &
                             .not.pvmetas(pvmetas(iv)%ideps(n))%outreq_L) changed_L = .true.
                        pvmetas(pvmetas(iv)%ideps(n))%outreq_L = .true.
                     endif
                  enddo DO_IDEP
               endif
            endif
         enddo DO_IV
      enddo DO_CHANGED
      
      !# Trim pvmetas list with dep conditions
      nphyvars2 = 0
      do iv = 1, nphyvars
         pvmetas(iv)%ideps = -1  !#TODO: is it usefull to get ideps again?
         if (pvmetas(iv)%outreq_L .or. .not.pvmetas(iv)%outcond_L) then
            nphyvars2 = nphyvars2 + 1
            if (nphyvars2 /= iv) pvmetas(nphyvars2) = pvmetas(iv)
            pvmetas(nphyvars2)%idxv = nphyvars2
            if (pvmetas(iv)%outcond_L) &
                 call msg(MSG_INFOPLUS, '(phymem_alloc) diag add: '//trim(pvmetas(iv)%vname))
         else
            call msg(MSG_INFOPLUS, '(phymem_alloc) diag rm : '//trim(pvmetas(iv)%vname))
         endif
      enddo
      nphyvars = nphyvars2

      !# Split into buses
      do iv = 1, nphyvars
         memgap = 0
         if (debug_mem_L) memgap = max(pvmetas(iv)%ni/5,1)+2*(ib-1)

         ib = pvmetas(iv)%ibus
         pvmetas(iv)%i0 = pbuses(ib)%nsize + memgap + 1
         pvmetas(iv)%in = pvmetas(iv)%i0 + pvmetas(iv)%size - 1
         pbuses(ib)%nsize = pvmetas(iv)%in
         pbuses(ib)%nvars = pbuses(ib)%nvars + 1
         pvmetas(iv)%idxb = pbuses(ib)%nvars
      enddo
      do ib = 1, PHY_NBUSES
         allocate(pbuses(ib)%mptr(pbuses(ib)%nvars))
         PHY_MAXVARS(ib) = pbuses(ib)%nvars
      enddo
      do iv = 1, nphyvars
         ib = pvmetas(iv)%ibus
         pbuses(ib)%mptr(pvmetas(iv)%idxb)%meta => pvmetas(iv)
      enddo

      !# allocate buses
      halo = 0
      DOBUS: do ib = 1, PHY_NBUSES
         
         if (pbuses(ib)%nsize == 0) continue
         
         init = 0
         call gmm_build_meta2D(gmmmeta, &
              1, pbuses(ib)%nsize, halo, halo, pbuses(ib)%nsize, &
              1, phydim_nj, halo, halo, phydim_nj, &
              init, GMM_NULL_FLAGS)
         init = GMM_FLAG_IZER + gmmflags(ib)
         if (F_debug) init = GMM_FLAG_INAN + gmmflags(ib)
         gmmname = phymem_gmmname(ib)
         nullify(pbuses(ib)%bptr)
         istat = gmm_create(gmmname, pbuses(ib)%bptr, gmmmeta, init)
         if (GMM_IS_OK(istat)) &
              istat = gmm_get(gmmname, pbuses(ib)%bptr)
         if (.not.GMM_IS_OK(istat)) then
            call msg_toall(MSG_ERROR, 'phymem_alloc: problem allocating memory for '//trim(gmmname)//' bus')
            return
         endif
          
      enddo DOBUS

      if (F_debug) then
         do iv = 1, nphyvars
            i0 = pvmetas(iv)%i0
            in = pvmetas(iv)%in
            ib = pvmetas(iv)%ibus
            pbuses(ib)%bptr(i0:in,:) = 0.
         enddo
      endif

      call priv_init_vlist()
      
      isallocated = .true.
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_alloc

   
   !/@*
   function phymem_find_idxv1(F_name, F_npath, F_bpath) result(F_idxv)
      implicit none
      !@objective Retreive idx of 1st var in pvmetas matching name+
      !@arguments
      character(len=*),  intent(in) :: F_name     !# Name of field to retrieve (input, variable or output name)
      character(len=*),  intent(in), optional :: F_npath    !# Name path to search ['VOI']
      character(len=*),  intent(in), optional :: F_bpath    !# Bus path to search ['PVD']
      !@return
      integer :: F_idxv !# idxv in pvmetas for 1st var matching provided params
      !*@/
      character(len=PHY_NAMELEN) :: npath, bpath
      integer :: istat, nmatch, idxvlist(1)
      !---------------------------------------------------------------
      F_idxv = PHY_ERROR

      npath = ' '
      if (present(F_npath)) then
         npath = F_npath ; call str_normalize(npath) ; istat = clib_toupper(npath)
      endif
      if (npath == ' ') npath = PHY_NPATH_DEFAULT

      bpath = ' '
      if (present(F_bpath)) then
         bpath = F_bpath ; call str_normalize(bpath) ; istat = clib_toupper(bpath)
      endif
      if (bpath == ' ') bpath = PHY_BPATH_DEFAULT
      
      nmatch = phymem_find_idxv(idxvlist, F_name, npath, bpath, F_quiet=.true.)

      if (nmatch > 0) F_idxv = idxvlist(1)
      !---------------------------------------------------------------
      return
   end function phymem_find_idxv1

   
   !/@*
   function phymem_find_idxv(F_idxvlist, F_name, F_npath, F_bpath, &
        F_quiet, F_shortmatch, F_endmatch, F_flagstr, F_shortflag) result(F_nmatch)
      implicit none
      !@objective Retreive list of var indices in pvmetas for matching ones
      !@arguments
      integer, intent(out) :: F_idxvlist(:)       !# List of indices in pvmetas for var matching provided params
      character(len=*),  intent(in), optional :: F_name     !# Name of field to retrieve (input, variable or output name)
      character(len=*),  intent(in), optional :: F_npath    !# Name path to search ['VOI']
      character(len=*),  intent(in), optional :: F_bpath    !# Bus path to search ['PVD']
      logical,           intent(in), optional :: F_quiet    !# Do not emit warning for unmatched entry [.false.]
      logical,           intent(in), optional :: F_shortmatch  !# if true, Match F_name against only the first len_trim(F_name) of input, variable or output name
      logical,           intent(in), optional :: F_endmatch !# if true, match F_name against only the last len_trim(F_name) of input, variable or output name
      character(len=*),  intent(in), optional :: F_flagstr  !# '+' separated list of flags to match
      logical,           intent(in), optional :: F_shortflag!# if true, match F_flagstr against only the first len_trim(F_name)
      !@return
      integer :: F_nmatch  !# PHY_ERROR or number of matching vars
      !*@/
      character(len=PHY_NAMELEN) :: name, npath, bpath, fname
      character(len=2) :: buses(PHY_NBUSES)
      character(len=1024) :: flagstr
      character(len=PHY_NAMELEN) :: flags(PHY_MAXFLAGS)
      integer :: nflags, nflags2, nbpath, ib, in, iv, ik, slen, iadd, istat, cnt, nnpath, &
           iaddf, istart, isub, ikk
      integer, dimension(PHY_MAXFLAGS) :: flen
      logical :: match_L, quiet_L
      !---------------------------------------------------------------
      F_nmatch = PHY_ERROR
      F_idxvlist = -1
      if (.not.isallocated) then
         call msg(MSG_ERROR,'(phymem_find) phymem_alloc must be called before')
         return
      endif

      name = ' '
      if (present(F_name)) then
         name  = F_name  ; call str_normalize(name)  ; istat = clib_tolower(name)
      endif
      
      npath = ' '
      if (present(F_npath)) then
         npath = F_npath ; call str_normalize(npath) ; istat = clib_toupper(npath)
      endif
      if (npath == ' ') npath = PHY_NPATH_DEFAULT

      bpath = ' '
      if (present(F_bpath)) then
         bpath = F_bpath ; call str_normalize(bpath) ; istat = clib_toupper(bpath)
      endif
      if (bpath == ' ') bpath = PHY_BPATH_DEFAULT

      quiet_L = .false.
      if (present(F_quiet)) quiet_L = F_quiet

      iadd = 1
      if (present(F_shortmatch)) then
         if (F_shortmatch) iadd = 0
      endif

      isub = PHY_NAMELEN
      if (present(F_endmatch)) then
         if (F_endmatch) isub = len_trim(name)
      endif
      
      iaddf = 1
      if (present(F_shortflag)) then
         if (F_shortflag) iaddf = 0
      endif
      
      flagstr = ' '
      if (present(F_flagstr)) then
         if (len_trim(F_flagstr) > len(flagstr)) then
            call msg(MSG_WARNING, '(phymem_find) F_flagstr is too long, some flags may be dropped; need to increase flagstr len in code/')
         endif
         flagstr = F_flagstr
      endif
      nflags = 0
      flags(:) = ' '
      if (flagstr /= ' ') then
         nflags = size(flags)
         call str_split2list(flags, flagstr, '+', nflags)
         nflags2 = 0
         do in = 1, nflags
            if (flags(in) == ' ') cycle
            nflags2 = nflags2+1
            if (in /= nflags2) flags(nflags2) = flags(in)
            call str_normalize(flags(nflags2))
            istat = clib_toupper(flags(nflags2))
            flen(nflags2) = min(max(1, len_trim(flags(nflags2))+iaddf),PHY_NAMELEN)
         enddo
         nflags = nflags2
      endif

      !# Loop to search for matching vars
      cnt  = 0
      buses = ' '
      nbpath = len_trim(bpath)
      do ib = 1, min(nbpath, PHY_NBUSES)
         buses(ib) = bpath(ib:ib)
      enddo
      
      slen = min(max(1,len_trim(name)+iadd),PHY_NAMELEN)
      nnpath = max(1,len_trim(npath))
      if (name == ' ') nnpath=1
      DO_NPATH: do in = 1, nnpath  !#TODO: avoid looping nnpath times

         DO_BPATH: do ib = 1, max(1, nbpath)

            DO_PVMETAS: do iv = 1, nphyvars
               !# check bus
               match_L = (nbpath == 0)
!!$               if (nbpath > 0) match_L = (pvmetas(iv)%bus == buses(ib))
               if (nbpath > 0) then
                  if (buses(ib) == 'U') then
                     match_L = (pvmetas(iv)%bus == 'V' .and. pvmetas(iv)%reset == 0)
                  else
                     match_L = (pvmetas(iv)%bus == buses(ib))
                  endif
               endif
               
               if (.not.match_L) cycle
               !# check names
               IF_NAME: if (name /= ' ') then
                  match_L = .false.
                  select case (npath(in:in))
                  case ('V')
                     fname = pvmetas(iv)%vname                  
                  case ('I')
                     fname = pvmetas(iv)%iname
                  case ('O')
                     fname = pvmetas(iv)%oname
                  case ('S')
                     fname = pvmetas(iv)%sname
                  case DEFAULT
                     call msg(MSG_WARNING,'(phymem_find) Ignoring unknown variable path entry '//npath(in:in))
                     cycle
                  end select
                  istart = max(1, len_trim(fname)-isub+1)
                  match_L = (fname(istart:istart+slen-1) == name)  
                  if (.not.match_L) cycle
               endif IF_NAME
               !# check flags
               IF_FLAG: if (nflags > 0) then
                  match_L = .false.
                  do ik = 1, nflags
                     if (flags(ik) == ' ') exit
                     do ikk = 1,size(pvmetas(iv)%flags)
                        if (pvmetas(iv)%flags(ikk) == ' ') exit
                        if (flags(ik) == pvmetas(iv)%flags(ikk)(1:flen(ik))) match_L = .true.
                     enddo
                  enddo
                  if (.not.match_L) cycle  
               endif IF_FLAG

               !# save matched indices
               if (cnt > 0) then
                  if (any(iv == F_idxvlist(1:cnt))) cycle
               endif
               if (cnt >= size(F_idxvlist)) then
                  if (.not.quiet_L) &
                       call msg(MSG_WARNING,'(phymem_find) F_idxvlist buffer overflow')
                  exit
               endif

               cnt = cnt + 1
               F_idxvlist(cnt) = iv

            enddo DO_PVMETAS
            if (cnt >= size(F_idxvlist)) exit
         enddo DO_BPATH
         if (cnt >= size(F_idxvlist)) exit

      enddo DO_NPATH
      
      if (cnt == 0 .and. .not.quiet_L) call msg(MSG_WARNING, &
           '(phymem_find) No matching entry found for name='// &
           trim(name)//', npath='//trim(npath)//', bpath='//trim(bpath)// &
           ', flags='//trim(flagstr))

      F_nmatch = cnt
      !---------------------------------------------------------------
      return
   end function phymem_find_idxv
   

   !/@*
   function phymem_getmeta(F_meta, F_idxv) result(F_istat)
      implicit none
      !@objective Return meta from pvmetas(idxv)
      !@arguments
      type(phymeta), pointer    :: F_meta  !# pvmetas(F_idxv)
      integer,       intent(in) :: F_idxv  !# pvmetas index of the field
      !@return
      integer :: F_istat
      !*@/
      character(len=32) :: msg_S
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      nullify(F_meta)
      if (.not.isallocated) then
         call msg(MSG_ERROR,'(phymem_get) phymem_alloc must be called before')
         return
      endif
      if (F_idxv < 1 .or. F_idxv > nphyvars) then
         write(msg_S,'(2i6)') F_idxv, nphyvars
         call msg(MSG_ERROR,'(phymem_get) Requested F_idxv out of range: '//msg_S)
         return
      endif
      F_meta => pvmetas(F_idxv)
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_getmeta

   
   !/@*
   function phymem_getmeta_copy(F_meta, F_idxv) result(F_istat)
      implicit none
      !@objective Return meta copy from pvmetas(idxv)
      !@arguments
      type(phymeta), intent(out) :: F_meta  !# pvmetas(F_idxv)
      integer,       intent(in)  :: F_idxv  !# pvmetas index of the field
      !@return
      integer :: F_istat
      !*@/
      character(len=32) :: msg_S
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (.not.isallocated) then
         call msg(MSG_ERROR,'(phymem_get) phymem_alloc must be called before')
         return
      endif
      if (F_idxv < 1 .or. F_idxv > nphyvars) then
         write(msg_S,'(2i6)') F_idxv, nphyvars
         call msg(MSG_ERROR,'(phymem_get) Requested F_idxv out of range: '//msg_S)
         return
      endif
      F_meta = pvmetas(F_idxv)
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_getmeta_copy
   
   !/@*
   function phymem_updatemeta(F_meta, F_idxv) result(F_istat)
      implicit none
      !@objective Return meta from pvmetas(idxv)
      !@arguments
      type(phymeta), intent(in) :: F_meta  !# new meta values
      integer,       intent(in) :: F_idxv  !# pvmetas index of the field
      !@return
      integer :: F_istat
      !*@/
      character(len=32) :: msg_S
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (.not.isallocated) then
         call msg(MSG_ERROR,'(phymem_updatemeta) phymem_alloc must be called before')
         return
      endif
      if (F_idxv < 1 .or. F_idxv > nphyvars) then
         write(msg_S,'(2i6)') F_idxv, nphyvars
         call msg(MSG_ERROR,'(phymem_updatemeta) Requested F_idxv out of range: '//msg_S)
         return
      endif

      !#TODO: check which one can be updated, insure consistency...
      pvmetas(F_idxv)%init  = F_meta%init
      pvmetas(F_idxv)%reset = F_meta%reset
      pvmetas(F_idxv)%stag  = F_meta%stag
      pvmetas(F_idxv)%wload = F_meta%wload
      pvmetas(F_idxv)%hzd   = F_meta%hzd
      pvmetas(F_idxv)%monot = F_meta%monot
      pvmetas(F_idxv)%massc = F_meta%massc
      pvmetas(F_idxv)%vmin  = F_meta%vmin
      pvmetas(F_idxv)%vmax  = F_meta%vmax
      pvmetas(F_idxv)%iname = F_meta%iname
      pvmetas(F_idxv)%oname = F_meta%oname
      pvmetas(F_idxv)%sname = F_meta%sname
      pvmetas(F_idxv)%desc  = F_meta%desc
      pvmetas(F_idxv)%flags(:) = F_meta%flags(:)
      
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_updatemeta

   
   !/@*
   function phymem_getdata1d(F_ptr, F_idxv, F_trnch) result(F_istat)
      implicit none
      !@objective Associate pointer 1d to pvmetas(idxv) for specified slice
      !@arguments
      real, pointer, contiguous :: F_ptr(:)  !# data(1:nikf)
      integer, intent(in) :: F_idxv   !# pvmetas index of the field
      integer, intent(in) :: F_trnch !# slice index of the field
      !@return
      integer :: F_istat
      !*@/
      integer :: i0, nikfm, in
      character(len=32) :: msg_S
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (.not.isallocated) then
         call msg(MSG_ERROR,'(phymem_get) phymem_alloc must be called before')
         return
      endif
      if (F_idxv < 1 .or. F_idxv > nphyvars) then
         write(msg_S,'(2i6)') F_idxv, nphyvars
         call msg(MSG_ERROR,'(phymem_get) Requested F_idxv out of range: '//msg_S)
         return
      endif
      if (F_trnch < 1 .or. F_trnch > phydim_nj) then
         call msg(MSG_ERROR,'(phymem_get) Requested F_trnch out of range')
         return
      endif
      
      i0 = pvmetas(F_idxv)%i0
      nikfm = pvmetas(F_idxv)%size
      in = pvmetas(F_idxv)%in
      F_ptr(1:nikfm) => pbuses(pvmetas(F_idxv)%ibus)%bptr(i0:in,F_trnch)
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_getdata1d

   
   !/@*
   function phymem_getdata2d(F_ptr, F_idxv, F_trnch) result(F_istat)
      implicit none
      !@objective Associate pointer 2d to pvmetas(idxv) for specified slice
      !@arguments
      real, pointer, contiguous :: F_ptr(:,:)  !# data(1:ni, 1:nkf)
      integer, intent(in) :: F_idxv   !# pvmetas index of the field
      integer, intent(in) :: F_trnch !# slice index of the field
      !@return
      integer :: F_istat
      !*@/
      integer :: i0, nkf, in, ni
      character(len=32) :: msg_S
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (.not.isallocated) then
         call msg(MSG_ERROR,'(phymem_get) phymem_alloc must be called before')
         return
      endif
      if (F_idxv < 1 .or. F_idxv > nphyvars) then
         call msg(MSG_ERROR,'(phymem_get) Requested F_idxv out of range: '//msg_S)
         return
      endif
      if (F_trnch < 1 .or. F_trnch > phydim_nj) then
         write(msg_S,'(2i6)') F_idxv, nphyvars
         call msg(MSG_ERROR,'(phymem_get) Requested F_trnch out of range')
         return
      endif
      
      i0 = pvmetas(F_idxv)%i0
      ni = pvmetas(F_idxv)%ni
      in = pvmetas(F_idxv)%in
      nkf = pvmetas(F_idxv)%nk * pvmetas(F_idxv)%fmul * (pvmetas(F_idxv)%mosaic+1)
      F_ptr(1:ni,1:nkf) => pbuses(pvmetas(F_idxv)%ibus)%bptr(i0:in,F_trnch)
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_getdata2d

   
   !/@*
   function phymem_getdata3d(F_ptr, F_idxv, F_trnch) result(F_istat)
      implicit none
      !@objective Associate pointer 3d to pvmetas(idxv) for specified slice
      !@arguments
      real, pointer, contiguous :: F_ptr(:,:,:)  !# data(1:ni, 1:nk, 1:nf)
      integer, intent(in) :: F_idxv   !# pvmetas index of the field
      integer, intent(in) :: F_trnch !# slice index of the field
      !@return
      integer :: F_istat
      !*@/
      integer :: i0, ni, in, nk, nfm
      character(len=32) :: msg_S
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (.not.isallocated) then
         call msg(MSG_ERROR,'(phymem_get) phymem_alloc must be called before')
         return
      endif
      if (F_idxv < 1 .or. F_idxv > nphyvars) then
         write(msg_S,'(2i6)') F_idxv, nphyvars
         call msg(MSG_ERROR,'(phymem_get) Requested F_idxv out of range: '//msg_S)
         return
      endif
      if (F_trnch < 1 .or. F_trnch > phydim_nj) then
         call msg(MSG_ERROR,'(phymem_get) Requested F_trnch out of range')
         return
      endif
      
      i0 = pvmetas(F_idxv)%i0
      ni = pvmetas(F_idxv)%ni
      in = pvmetas(F_idxv)%in
      nk = pvmetas(F_idxv)%nk
      nfm = pvmetas(F_idxv)%fmul * (pvmetas(F_idxv)%mosaic+1)
      F_ptr(1:ni,1:nk,1:nfm) => pbuses(pvmetas(F_idxv)%ibus)%bptr(i0:in,F_trnch)
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_getdata3d

   
   !/@*
   function phymem_get_slabvars(F_pvars, F_trnch) result(F_istat)
      implicit none
      !@objective Get the list of vars meta+data on slab/slice
      !@arguments
      type(phyvar), pointer, contiguous :: F_pvars(:)
      integer, intent(in) :: F_trnch !# slice/slab index
      !@return
      integer :: F_istat
      !*@/
      integer :: ivar, istat, ibus
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (.not.isallocated) then
         call msg(MSG_ERROR,'(phymem_get_slabvars) phymem_alloc must be called before')
         return
      endif
      if (F_trnch < 1 .or. F_trnch > phydim_nj) then
         call msg(MSG_ERROR,'(phymem_get_slabvars) Requested F_trnch out of range')
         return
      endif
      if (associated(F_pvars)) then
         call msg(MSG_ERROR,'(phymem_get_slabvars) F_pvars allready associated')
         return
      endif
      
      allocate(F_pvars(nphyvars+PHY_NBUSES), stat=istat)
      if (istat /= 0) then
         call msg(MSG_ERROR, '(phymem) allocate problem in phymem_get_slabvars')
         return
      endif

      !# Fill every phyvar with meta and slab data pointer
      do ivar = 1, nphyvars
         F_pvars(ivar)%meta => pvmetas(ivar)
         istat = phymem_getdata1d(F_pvars(ivar)%data, ivar, F_trnch)
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_ERROR,'(phymem_get_slabvars) Problem getting pointer for '//pvmetas(ivar)%vname)
            return
         endif
      enddo

      !# Add a phyvar for each buses... (phymem internal feature)
      !# the main purpose of this is to be used with phymem_reset
      !# and other whole memory bus slab operation (much faster than looping through vars)
      do ibus = 1, PHY_NBUSES
         ivar = nphyvars+ibus
         F_pvars(ivar)%meta => NULL()
         F_pvars(ivar)%data => pbuses(ibus)%bptr(:,F_trnch)
         PHY_BUSIDXV(ibus) = ivar
      enddo
      
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_get_slabvars

   
   !#---- Private functions ------------------------------------------

   !/@*
   function priv_init(F_dnvars) result(F_istat)
      implicit none
      integer, intent(in), optional :: F_dnvars
      !@return
      integer :: F_istat
      !*@/
      integer :: ib, istat, dnvars2
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (isinit) then
         call msg(MSG_ERROR, '(phymem) called priv_init twice')
         return
      endif
      dnvars2 = DNVARS
      if (present(F_dnvars)) dnvars2 = max(DNVARS, F_dnvars)
      isinit = .true.
      do ib = 1, PHY_NBUSES
         pbuses(ib)%busid = PHY_BUSID(ib)
         pbuses(ib)%nvars = 0
         pbuses(ib)%nsize = 0
         nullify(pbuses(ib)%mptr)
         nullify(pbuses(ib)%bptr)
      enddo
      nullify(pvmetas)
      nphyvarsmax = 0
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function priv_init

   
   !/@*
   function priv_resize(F_dnvars) result(F_istat)
      integer, intent(in) :: F_dnvars
      integer :: F_istat
      !*@/
      integer :: n, istat
      type(phymeta), pointer, contiguous :: pvmetasnew(:)
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (.not.isinit) istat = priv_init()
      if (isallocated) then
         call msg(MSG_ERROR, '(phymem_resize) Cannot add var after allocation')
         return
      endif
      if (F_dnvars >= 0) then
         nphyvarsmax = nphyvarsmax + max(0, F_dnvars)
         allocate(pvmetasnew(nphyvarsmax), stat=istat)
         if (associated(pvmetas)) then
            do n=1,nphyvars
               pvmetasnew(n) = pvmetas(n)
            enddo
            deallocate(pvmetas)
         endif
         pvmetas => pvmetasnew
      endif
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function priv_resize

  
   !/@*
   function priv_checkvar(F_meta) result(F_istat)
      implicit none
      type(phymeta), intent(in) :: F_meta
      integer :: F_istat
      !*@/
      integer :: n
      !---------------------------------------------------------------
      F_istat = PHY_ERROR

      if (F_meta%vname == '') then
         call msg(MSG_ERROR, '(phymem::checkvar) must provide a vname')
         return
      endif
      if (F_meta%desc == '') then
         call msg(MSG_ERROR, '(phymem::checkvar) must provide a description for vname='//trim(F_meta%desc))
         return
      endif
      if (F_meta%ni < 1 .or. F_meta%nk < 1 .or. F_meta%fmul < 1 .or. F_meta%mosaic < 0) then
         call msg(MSG_ERROR, '(phymem::checkvar) must provide positive values for ni, nk ,fmul and mosaic for vname='//trim(F_meta%desc))
         return
      endif
            
      do n = 1, nphyvars
         if (F_meta%vname == pvmetas(n)%vname) then
            call msg(MSG_ERROR, '(phymem::checkvar) Duplicate entry (vname) for: '//trim(F_meta%vname))
            return
         endif
         if (F_meta%oname == pvmetas(n)%oname) then
            call msg(MSG_ERROR, '(phymem::checkvar) Duplicate entry (oname) for: '//trim(F_meta%oname)//' ('//trim(F_meta%vname)//')')
            return
         endif
!!$            if (F_meta%iname == pvmetas(n)%iname) then
!!$               call msg(MSG_ERROR, '(phymem::checkvar) Duplicate entry (iname) for: '//trim(F_meta%iname)//' ('//trim(F_meta%vname)//')')
!!$               return
!!$            endif
         if (F_meta%sname == pvmetas(n)%sname) then
            call msg(MSG_ERROR, '(phymem::checkvar) Duplicate entry (sname) for: '//trim(F_meta%sname)//' ('//trim(F_meta%vname)//')')
            return
         endif
         if (F_meta%desc == pvmetas(n)%desc) then
            call msg(MSG_ERROR, '(phymem::checkvar) Duplicate entry (desc) for: ('//trim(F_meta%vname)//') '//trim(F_meta%desc))
            return
         endif
      enddo
      
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function priv_checkvar

   
   !/@*
   function priv_find(F_vname) result(F_idxv)
      character(len=*), intent(in) :: F_vname
      integer :: F_idxv
      !*@/
      integer :: n
      !---------------------------------------------------------------
      F_idxv = PHY_ERROR
      do n = 1, nphyvars
         if (F_vname == pvmetas(n)%vname) then
            F_idxv = n
            return
         endif
      enddo
      !---------------------------------------------------------------
      return
   end function priv_find

   
   !/@*
   subroutine priv_init_vlist()
      !*@/
      integer :: n
      !---------------------------------------------------------------
      if (.not.associated(uidxvlist)) allocate(uidxvlist(nphyvars))
      nuidxvlist = 0
      uidxvlist = -1
      do n=1,nphyvars
         if (pvmetas(n)%ibus == PHY_VBUSIDX) then
            if (pvmetas(n)%reset == 1) then
               nuidxvlist = nuidxvlist + 1
               uidxvlist(nuidxvlist) = n
            endif
         endif
      enddo
      !---------------------------------------------------------------
      return
   end subroutine priv_init_vlist
   
end module phymem
