module phymem
   use clib_itf_mod, only: clib_toupper, clib_tolower
   use rmn_gmm, only: gmm_metadata, gmm_create, gmm_get, GMM_IS_OK, GMM_FLAG_RSTR, GMM_FLAG_IZER, GMM_FLAG_INAN, GMM_NULL_FLAGS, GMM_MAXNAMELENGTH
   use str_mod, only: str_normalize
   use phy_status, only: PHY_OK, PHY_ERROR
   use splitst, only: splitst4
   use phygridmap, only: phy_lcl_ni, phy_lcl_nj, phydim_ni, phydim_nj, phydim_nk
   use phy_options, only: debug_mem_L
   implicit none
   private

#include <rmnlib_basics.hf>
#include <rmn/msg.h>
   
   !# Public functions
   public :: phymem_busidx, phymem_alloc, phymem_init
   public :: phymem_add, phymem_gmmname, phymem_busreset
   public :: phymem_find, phymem_find_idxv
   public :: phymem_getmeta, phymem_getmeta_copy, phymem_getdata
   public :: phymem_get_slabvars, phymem_get_i0_string
   
   !# Public constants
   integer, parameter, public :: PHY_NBUSES = 4
   integer, parameter, public :: PHY_NAMELEN  = 32
   integer, parameter, public :: PHY_DESCLEN  = 256
   integer, parameter, public :: PHY_MAXFLAGS = 16
   
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
      integer :: idxv    !# var index in the specified bus, pvmetas(idxv)%meta
      integer :: i0      !# index of first element in the bus pointer, pbuses(ibus)%bptr(i0:in,:)
      integer :: in      !# in=i0+size-1; index of first element in the bus pointer, pbuses(ibus)%bptr(i0:in,:)      
      integer :: init    !# 1 = init/mandatory, 0 otherwise
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
      type(phymeta), pointer, contiguous :: meta(:) => null() !# ptr to var meta
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
   end interface phymem_find

   interface phymem_getdata
      module procedure phymem_getdata1d
      module procedure phymem_getdata2d
      module procedure phymem_getdata3d
   end interface phymem_getdata
   
   !# Private module vars --------------------------------------------------
   
   integer, parameter :: DNVARS = 25

   type(phybus), save  :: pbuses(PHY_NBUSES)
   type(phymetaptr), pointer, contiguous, save :: pvmetas(:) => null()

   logical, save :: isallocated = .false.
   logical, save :: isinit = .false.

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
   function phymem_busreset2(F_pvars, F_busidx, F_value) result(F_istat)
      implicit none
      !@objective Reset whole memory bus to specified value (default zero) for provided memory slice structure (pvars)
      !@arguments
      type(phyvar), pointer, contiguous :: F_pvars(:)
      integer, intent(in) :: F_busidx
      real, intent(in), optional :: F_value
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
      F_pvars(PHY_BUSIDXV(F_busidx))%data(:) = rvalue
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_busreset2

         
   !/@*
   function phymem_add_meta(F_imeta, F_i0) result(F_idxv)
      !@objective  Add physics var from pre-filler phymeta
      !@arguments
      type(phymeta), intent(in) :: F_imeta
      integer, intent(out), optional :: F_i0  !# Temporarily return i0 for backward compat
      !@return
      integer :: F_idxv  !# var index in pvmetas
      !*@/
      type(phymeta) :: vmeta
      integer :: istat, nkfm, idxb, idxv, ibus, memgap, vsize
      !---------------------------------------------------------------
      F_idxv = PHY_ERROR
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
      if (pbuses(ibus)%nvars >= PHY_MAXVARS(ibus)) then
         istat = priv_resize(ibus, DNVARS)
      endif
      
      vmeta = F_imeta
      vmeta%ibus = ibus
      call str_normalize(vmeta%iname)
      call str_normalize(vmeta%oname)
      call str_normalize(vmeta%vname)
      call str_normalize(vmeta%sname)
      istat = clib_toupper(vmeta%bus)

      if (vmeta%oname == '') vmeta%oname = vmeta%vname
      if (vmeta%iname == '') vmeta%iname = vmeta%oname
      if (vmeta%sname == '') vmeta%sname = vmeta%oname
      
      vmeta%ni = phydim_ni

      istat = priv_checkvar(vmeta)
      if (istat /= PHY_OK) return

      idxv = nphyvars + 1
      idxb = pbuses(ibus)%nvars + 1
      nkfm = vmeta%nk*vmeta%fmul*(vmeta%mosaic+1)
      vsize = vmeta%ni * nkfm
      memgap = 0
      if (debug_mem_L) memgap = max(vmeta%ni/5,1)+2*(ibus-1)

      vmeta%idxb = idxb
      vmeta%idxv = idxv
      vmeta%i0 = pbuses(ibus)%nsize + memgap + 1
      vmeta%in = vmeta%i0 + vsize - 1
      vmeta%nlcl = (/phy_lcl_ni, phy_lcl_nj, nkfm/)
      vmeta%size = vsize

      pbuses(ibus)%meta(idxb) = vmeta

      pbuses(ibus)%nsize = vmeta%in
      pbuses(ibus)%nvars = idxb
      nphyvars = idxv
      
      F_idxv = idxv
      if (present(F_i0)) F_i0 = vmeta%i0
      !---------------------------------------------------------------
      return
   end function phymem_add_meta

   
   !/@*
   function phymem_add_string(F_string, F_i0) result(F_idxv)
      !@objective Add physics var from string description
      !@arguments
      character(len=*), intent(in) :: F_string
      integer, intent(out), optional :: F_i0  !# Temporarily return i0 for backward compat
      !@return
      integer :: F_idxv  !# var index in pvmetas
      !@Notes
      !  See splitst.F90 for the input string description
      !*@/
      integer, parameter :: FSTNAMELEN = 4
      integer :: istat
      character(len=3) ::   shape
      type(phymeta) :: vmeta
      !---------------------------------------------------------------
      F_idxv = -1
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
      
      if (.not.any(vmeta%bus == PHY_BUSID(:)))  then
         call physeterror('phymem_add', 'VB=(BUS) NOT ALLOWED: '//trim(F_string))
         return
      endif

      vmeta%nk = phydim_nk
      if (shape == "A") vmeta%nk = 1

      if (present(F_i0)) then
         F_idxv = phymem_add_meta(vmeta, F_i0)
      else
         F_idxv = phymem_add_meta(vmeta)
      endif
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
      istat = phymem_find(idxvlist, F_name=vmeta%vname, F_npath='V', F_bpath=vmeta%bus)
      if (istat <= 0) then
         call physeterror('phymem_get_idxv_string', 'Cannot find var for: '//trim(F_string))
         return
      endif
      F_i0 = pvmetas(idxvlist(1))%meta%i0
      !---------------------------------------------------------------
      return
   end function phymem_get_i0_string

   
   !/@*
   function phymem_alloc(F_debug) result(F_istat)
      implicit none
      !@objective Allocate memory buses for previously added vars
      !@arguments
      logical, intent(in) :: F_debug
      !@return
      integer :: F_istat
      !*@/
      integer :: ib, iv, idxv, halo, init, i0, in, istat
      character(len=GMM_MAXNAMELENGTH) :: gmmname
      type(gmm_metadata) :: gmmmeta
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
      allocate(pvmetas(nphyvars), stat=istat)
      if (istat /= 0) then
         call msg(MSG_ERROR, '(phymem) allocate problem in phymem_alloc')
         return
      endif

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

         do iv = 1, pbuses(ib)%nvars
            idxv = pbuses(ib)%meta(iv)%idxv
            pvmetas(idxv)%meta => pbuses(ib)%meta(iv)

            if (F_debug) then
               i0 = pbuses(ib)%meta(iv)%i0
               in = pbuses(ib)%meta(iv)%in
               pbuses(ib)%bptr(i0:in,:) = 0.
            endif

         enddo
          
      enddo DOBUS
      isallocated = .true.
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_alloc

   
   !/@*
   function phymem_find_idxv(F_idxvlist, F_name, F_npath, F_bpath, &
        F_quiet, F_shortmatch, F_flagstr) result(F_nmatch)
      implicit none
      !@objective Retreive list of var indices in pvmetas for matching ones
      !@arguments
      integer, intent(out) :: F_idxvlist(:)       !# List of indices in pvmetas for var matching provided params
      character(len=*),  intent(in), optional :: F_name     !# Name of field to retrieve (input, variable or output name)
      character(len=*),  intent(in), optional :: F_npath    !# Name path to search ['VOI']
      character(len=*),  intent(in), optional :: F_bpath    !# Bus path to search ['PVD']
      logical,           intent(in), optional :: F_quiet    !# Do not emit warning for unmatched entry [.false.]
      logical,           intent(in), optional :: F_shortmatch  !# if true, Match F_name against only the first len_trim(F_name) of input, variable or output name
      character(len=*),  intent(in), optional :: F_flagstr  !# '+' separated list of flags to match
      !@return
      integer :: F_nmatch  !# PHY_ERROR or number of matching vars
      !*@/
      character(len=PHY_NAMELEN) :: name, npath, bpath
      character(len=2) :: buses(PHY_NBUSES)
      character(len=1024) :: flagstr
      character(len=PHY_NAMELEN) :: flags(PHY_MAXFLAGS)
      integer :: nflags, nflags2, nbpath, ib, in, iv, ik, slen, iadd, istat, cnt, nnpath
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
               if (nbpath > 0) match_L = (pvmetas(iv)%meta%bus == buses(ib))
               if (.not.match_L) cycle
               !# check names
               IF_NAME: if (name /= ' ') then
                  match_L = .false.
                  select case (npath(in:in))
                  case ('V')
                     match_L = (pvmetas(iv)%meta%vname(1:slen) == name)
                  case ('I')
                     match_L = (pvmetas(iv)%meta%iname(1:slen) == name)
                  case ('O')
                     match_L = (pvmetas(iv)%meta%oname(1:slen) == name)
                  case ('S')
                     match_L = (pvmetas(iv)%meta%sname(1:slen) == name)
                  case DEFAULT
                     call msg(MSG_WARNING,'(phymem_find) Ignoring unknown variable path entry '//npath(in:in))
                     cycle
                  end select
                  if (.not.match_L) cycle           
               endif IF_NAME
               !# check flags
               IF_FLAG: if (nflags > 0) then
                  match_L = .true.
                  do ik = 1, nflags
                     if (flags(ik) == ' ') exit
                     if (.not.any(flags(ik) == pvmetas(iv)%meta%flags(:))) then
                        match_L = .false.
                        exit
                     endif
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
      type(phymeta), pointer    :: F_meta  !# pvmetas(F_idxv)%meta
      integer,       intent(in) :: F_idxv  !# pvmetas index of the field
      !@return
      integer :: F_istat
      !*@/
      integer :: i0, nikfm, in
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      nullify(F_meta)
      if (.not.isallocated) then
         call msg(MSG_ERROR,'(phymem_get) phymem_alloc must be called before')
         return
      endif
      if (F_idxv < 1 .or. F_idxv > nphyvars) then
         call msg(MSG_ERROR,'(phymem_get) Requested F_idxv out of range')
         return
      endif
      F_meta => pvmetas(F_idxv)%meta
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_getmeta

   
   !/@*
   function phymem_getmeta_copy(F_meta, F_idxv) result(F_istat)
      implicit none
      !@objective Return meta copy from pvmetas(idxv)
      !@arguments
      type(phymeta), intent(out) :: F_meta  !# pvmetas(F_idxv)%meta
      integer,       intent(in)  :: F_idxv  !# pvmetas index of the field
      !@return
      integer :: F_istat
      !*@/
      integer :: i0, nikfm, in
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (.not.isallocated) then
         call msg(MSG_ERROR,'(phymem_get) phymem_alloc must be called before')
         return
      endif
      if (F_idxv < 1 .or. F_idxv > nphyvars) then
         call msg(MSG_ERROR,'(phymem_get) Requested F_idxv out of range')
         return
      endif
      F_meta = pvmetas(F_idxv)%meta
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_getmeta_copy
   
   
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
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (.not.isallocated) then
         call msg(MSG_ERROR,'(phymem_get) phymem_alloc must be called before')
         return
      endif
      if (F_idxv < 1 .or. F_idxv > nphyvars) then
         call msg(MSG_ERROR,'(phymem_get) Requested F_idxv out of range')
         return
      endif
      if (F_trnch < 1 .or. F_trnch > phydim_nj) then
         call msg(MSG_ERROR,'(phymem_get) Requested F_trnch out of range')
         return
      endif
      
      i0 = pvmetas(F_idxv)%meta%i0
      nikfm = pvmetas(F_idxv)%meta%size
      in = pvmetas(F_idxv)%meta%in
      F_ptr(1:nikfm) => pbuses(pvmetas(F_idxv)%meta%ibus)%bptr(i0:in,F_trnch)
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
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (.not.isallocated) then
         call msg(MSG_ERROR,'(phymem_get) phymem_alloc must be called before')
         return
      endif
      if (F_idxv < 1 .or. F_idxv > nphyvars) then
         call msg(MSG_ERROR,'(phymem_get) Requested F_idxv out of range')
         return
      endif
      if (F_trnch < 1 .or. F_trnch > phydim_nj) then
         call msg(MSG_ERROR,'(phymem_get) Requested F_trnch out of range')
         return
      endif
      
      i0 = pvmetas(F_idxv)%meta%i0
      ni = pvmetas(F_idxv)%meta%ni
      in = pvmetas(F_idxv)%meta%in
      nkf = pvmetas(F_idxv)%meta%nk * pvmetas(F_idxv)%meta%fmul * (pvmetas(F_idxv)%meta%mosaic+1)
      F_ptr(1:ni,1:nkf) => pbuses(pvmetas(F_idxv)%meta%ibus)%bptr(i0:in,F_trnch)
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
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (.not.isallocated) then
         call msg(MSG_ERROR,'(phymem_get) phymem_alloc must be called before')
         return
      endif
      if (F_idxv < 1 .or. F_idxv > nphyvars) then
         call msg(MSG_ERROR,'(phymem_get) Requested F_idxv out of range')
         return
      endif
      if (F_trnch < 1 .or. F_trnch > phydim_nj) then
         call msg(MSG_ERROR,'(phymem_get) Requested F_trnch out of range')
         return
      endif
      
      i0 = pvmetas(F_idxv)%meta%i0
      ni = pvmetas(F_idxv)%meta%ni
      in = pvmetas(F_idxv)%meta%in
      nk = pvmetas(F_idxv)%meta%nk
      nfm = pvmetas(F_idxv)%meta%fmul * (pvmetas(F_idxv)%meta%mosaic+1)
      F_ptr(1:ni,1:nk,1:nfm) => pbuses(pvmetas(F_idxv)%meta%ibus)%bptr(i0:in,F_trnch)
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
         F_pvars(ivar)%meta => pvmetas(ivar)%meta
         istat = phymem_getdata1d(F_pvars(ivar)%data, ivar, F_trnch)
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_ERROR,'(phymem_get_slabvars) Problem getting pointer for '//pvmetas(ivar)%meta%vname)
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
         nullify(pbuses(ib)%meta)
         istat = priv_resize(ib, dnvars2)
         if (istat /= PHY_OK) then
            call msg(MSG_ERROR, '(phymem) allocate problem in phymem_init')
            return
         endif
         nullify(pbuses(ib)%bptr)
         nullify(pvmetas)
      enddo
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function priv_init

   
   !/@*
   function priv_resize(F_ib, F_dnvars) result(F_istat)
      integer, intent(in) :: F_ib, F_dnvars
      integer :: F_istat
      !*@/
      integer :: n, istat
      type(phymeta), pointer, contiguous :: metaptr(:)
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (.not.isinit) istat = priv_init()
      if (isallocated) then
         call msg(MSG_ERROR, '(phymem_resize) Cannot add var after allocation')
         return
      endif
      if (max(1, min(F_ib, PHY_NBUSES)) /= F_ib) then
         call msg(MSG_ERROR, '(phymem_resize) Bus number out of range')
         return         
      endif
      PHY_MAXVARS(F_ib) = PHY_MAXVARS(F_ib) + max(0, F_dnvars)
!!$      print *, '(phymem_resize) '//PHY_BUSID(F_ib)//':', PHY_MAXVARS(F_ib)
      allocate(metaptr(PHY_MAXVARS(F_ib)), stat=istat)
      if (associated(pbuses(F_ib)%meta)) then
         do n=1,pbuses(F_ib)%nvars
            metaptr(n) = pbuses(F_ib)%meta(n)
         enddo
         deallocate(pbuses(F_ib)%meta)
      endif
      pbuses(F_ib)%meta => metaptr
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
      integer :: ibus, idxb
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
      
      DO_BUS: do ibus = 1, PHY_NBUSES
         DO_VAR: do idxb = 1, pbuses(ibus)%nvars
            if (F_meta%vname == pbuses(ibus)%meta(idxb)%vname) then
               call msg(MSG_ERROR, '(phymem::checkvar) Duplicate entry (vname) for: '//trim(F_meta%vname))
               return
            endif
            if (F_meta%oname == pbuses(ibus)%meta(idxb)%oname) then
               call msg(MSG_ERROR, '(phymem::checkvar) Duplicate entry (oname) for: '//trim(F_meta%oname)//' ('//trim(F_meta%vname)//')')
               return
            endif
!!$            if (F_meta%iname == pbuses(ibus)%meta(idxb)%iname) then
!!$               call msg(MSG_ERROR, '(phymem::checkvar) Duplicate entry (iname) for: '//trim(F_meta%iname)//' ('//trim(F_meta%vname)//')')
!!$               return
!!$            endif
            if (F_meta%sname == pbuses(ibus)%meta(idxb)%sname) then
               call msg(MSG_ERROR, '(phymem::checkvar) Duplicate entry (sname) for: '//trim(F_meta%sname)//' ('//trim(F_meta%vname)//')')
               return
            endif
            if (F_meta%desc == pbuses(ibus)%meta(idxb)%desc) then
               call msg(MSG_ERROR, '(phymem::checkvar) Duplicate entry (desc) for: ('//trim(F_meta%vname)//') '//trim(F_meta%desc))
               return
            endif
         enddo DO_VAR
      enddo DO_BUS
      
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function priv_checkvar
   
   
end module phymem
