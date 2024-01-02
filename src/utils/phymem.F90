module phymem
   use clib_itf_mod, only: clib_toupper, clib_tolower
   use rmn_gmm, only: gmm_metadata, gmm_create, gmm_get, GMM_IS_OK, GMM_FLAG_RSTR, GMM_FLAG_IZER, GMM_FLAG_INAN, GMM_NULL_FLAGS, GMM_MAXNAMELENGTH
   use str_mod, only: str_normalize
   use phy_status, only: PHY_OK, PHY_ERROR
   implicit none
   private
   
   public :: phymem_busidx, phymem_init, phymem_alloc, phymem_isalloc
   public :: phymem_gmmname, phymem_busreset
   public :: phymem_find, phymem_find_idxv, phymem_find_var !# , phymem_add
   public :: phymem_getmeta, phymem_getptr, phymem_getptr1d, phymem_getptr2d, phymem_getptr3d
   public :: phymeta, phybus, phyvar
   
   public :: pbuslist, pvarlist, npvarlist

#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   integer, parameter, public :: PHY_NBUSES = 4
   integer, parameter, public :: PHY_MAXVARS = 1000  !# max nb vr per bus (need at least double for all 4: phyvar has ~900 entries
   integer, parameter, public :: PHY_NAMELEN  = 32
   integer, parameter, public :: PHY_DESCLEN  = 256
   integer, parameter, public :: PHY_MAXFLAGS = 16
   
   integer, parameter, public :: PHY_STAG_SFC = 0
   integer, parameter, public :: PHY_STAG_MOM = 0
   integer, parameter, public :: PHY_STAG_THERMO = 1
   integer, parameter, public :: PHY_STAG_ENERGY = 2

   character(len=*), parameter, public :: PHY_NPATH_DEFAULT='VOI'
   character(len=*), parameter, public :: PHY_BPATH_DEFAULT='VPDE'

   integer, parameter, public :: PHY_EBUSIDX = 1
   integer, parameter, public :: PHY_DBUSIDX = 2
   integer, parameter, public :: PHY_PBUSIDX = 3
   integer, parameter, public :: PHY_VBUSIDX = 4
   character(len=1), parameter, public :: PHY_BUSID(PHY_NBUSES) = (/'E', 'D', 'P', 'V'/)
   integer, parameter, public :: gmmflags(PHY_NBUSES) = (/ 0, 0, GMM_FLAG_RSTR, 0/)


   ! wload : A flag defining whether or not the field should be used in the density calculation (water loaded).
   ! hzd   : tracer attribute, a flag defining whether or not the perform horizontal diffusion
   ! monot : tracer attribute, a flag defining whether or not the use monoton interpolation in the advection
   ! massc : tracer attribute, a flag defining whether or not the use a mass conservation scheme in the advection

   type phymeta
      integer :: ni      !# folded ni dim (p_runlenght)
      integer :: nk      !# number of atmospheric levels.
      integer :: fmul    !# number of arbitrarily-defined levels
      integer :: mosaic  !# number of surface sub-types
      integer :: size    !# ni * nk * fmul * (mosaic+1)
      integer :: nlcl(3) !# local tile dims in "not folded" space
      integer :: ibus    !# index of bus containing the field
      integer :: idxb    !# var index in the specified bus, pbuslist(ibus)%meta(idxb)
      integer :: idxv    !# var index in the specified bus, pvarlist(idxv)%meta
      integer :: i0      !# index of first element in the bus pointer, pbuslist(ibus)%bptr(i0:in,:)
      integer :: in      !# in=i0+size-1; index of first element in the bus pointer, pbuslist(ibus)%bptr(i0:in,:)      
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
      character(len=PHY_DESCLEN) :: desc !# Field description
      character(len=PHY_NAMELEN) :: flags(PHY_MAXFLAGS)  !# list of text key describing the fields caracteristics

      real, pointer, contiguous :: bptr(:,:) => null()  !# pointer to the bus containing the field var pointer would be bptr(i0:in,:)
   end type phymeta
   
   type phybus
      integer :: nvars = 0
      integer :: nsize = 0
      type(phymeta), pointer, contiguous :: meta(:) => null()
      real, pointer, contiguous :: bptr(:,:) => null()
      character(len=32) :: busid = ' '  !# size defined to have allignment
   end type phybus
   
   type phyvar
      type(phymeta), pointer :: meta => null()
   end type phyvar

   interface phymem_find
      module procedure phymem_find_idxv
      module procedure phymem_find_var
   end interface phymem_find

   interface phymem_getptr
      module procedure phymem_getptr1d
      module procedure phymem_getptr2d
      module procedure phymem_getptr3d
   end interface phymem_getptr

   integer, save :: npvarlist = 0
   type(phybus), save  :: pbuslist(PHY_NBUSES)
   type(phyvar), pointer, contiguous, save :: pvarlist(:) => null()

   logical, save :: isallocated = .false.
   logical, save :: isinit = .false.

contains

   !/@*
   function phymem_busidx(F_busname, F_quiet) result(F_busidx)
      implicit none
      !@objective 
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
      !@objective 
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
   function phymem_busreset(F_busidx, F_value) result(F_istat)
      implicit none
      integer, intent(in) :: F_busidx
      real, intent(in), optional :: F_value
      !@objective 
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
      pbuslist(F_busidx)%bptr = rvalue
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_busreset

   
   !/@*
   function phymem_init() result(F_istat)
      implicit none
      !@objective 
      !@return
      integer :: F_istat
      !*@/
      integer :: ib
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (isinit) then
         call msg(MSG_ERROR, '(phymem) called phymem_init twice')
         return
      endif
      isinit = .true.
      do ib = 1, PHY_NBUSES
         pbuslist(ib)%busid = PHY_BUSID(ib)
         pbuslist(ib)%nvars = 0
         pbuslist(ib)%nsize = 0
         allocate(pbuslist(ib)%meta(PHY_MAXVARS))  !TODO: make MAXVAR dynamic
         nullify(pbuslist(ib)%bptr)
         nullify(pvarlist)
      enddo
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_init

   
    !/@*
  function phymem_isalloc() result(F_isalloc)
      implicit none
      !@objective 
      !@return
      logical :: F_isalloc
      !*@/
      !---------------------------------------------------------------
      F_isalloc = isallocated
      !---------------------------------------------------------------
      return
   end function phymem_isalloc

   
!!$   !/@*
!!$   function phymem_add(F_meta) result(F_istat)
!!$      type(phymeta), intent(in) :: F_meta
!!$      !*@/
!!$      !---------------------------------------------------------------
!!$      if (.not.isinit) F_istat = phymem_init()
!!$      F_istat = PHY_ERROR
!!$      if (isallocated) then
!!$         call msg(MSG_ERROR, '(phymem_addvar) Cannot add var after allocation')
!!$         return
!!$      endif
!!$
!!$      
!!$
!!$      F_istat = PHY_OK
!!$      !---------------------------------------------------------------
!!$      return
!!$   end function phymem_add


   !/@*
   function phymem_alloc(F_nj, F_debug) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      integer, intent(in) :: F_nj
      logical, intent(in) :: F_debug
      !@return
      integer :: F_istat
      !*@/
      integer :: ib, iv, halo, init, i0, in, istat
      character(len=GMM_MAXNAMELENGTH) :: gmmname
      type(gmm_metadata) :: gmmmeta
      !---------------------------------------------------------------
      F_istat = PHY_ERROR
      if (isallocated) then
         call msg(MSG_ERROR, '(phymem) called phymem_alloc twice')
         return
      endif
      if (.not.isinit) istat = phymem_init()
      halo = 0
      
      npvarlist = 0
      do ib = 1, PHY_NBUSES
         npvarlist = npvarlist + pbuslist(ib)%nvars
      enddo
      allocate(pvarlist(npvarlist))
      npvarlist = 0
      
      DOBUS: do ib = 1, PHY_NBUSES
         
         if (pbuslist(ib)%nsize == 0) continue
         
         init = 0
         call gmm_build_meta2D(gmmmeta, &
              1, pbuslist(ib)%nsize, halo, halo, pbuslist(ib)%nsize, &
              1, F_nj, halo, halo, F_nj, &
              init, GMM_NULL_FLAGS)
         init = GMM_FLAG_IZER + gmmflags(ib)
         if (F_debug) init = GMM_FLAG_INAN + gmmflags(ib)
         gmmname = phymem_gmmname(ib)
         nullify(pbuslist(ib)%bptr)
         istat = gmm_create(gmmname, pbuslist(ib)%bptr, gmmmeta, init)
         if (GMM_IS_OK(istat)) &
              istat = gmm_get(gmmname, pbuslist(ib)%bptr)
         if (.not.GMM_IS_OK(istat)) then
            call msg_toall(MSG_ERROR, 'phymem_alloc: problem allocating memory for '//trim(gmmname)//' bus')
            return
         endif

         do iv = 1, pbuslist(ib)%nvars
            pbuslist(ib)%meta(iv)%bptr => pbuslist(ib)%bptr
            
            npvarlist = npvarlist + 1
            pbuslist(ib)%meta(iv)%idxv = npvarlist
            pvarlist(npvarlist)%meta => pbuslist(ib)%meta(iv)

            if (F_debug) then
               i0 = pbuslist(ib)%meta(iv)%i0
               in = pbuslist(ib)%meta(iv)%in
               pbuslist(ib)%bptr(i0:in,:) = 0.
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
      !@objective Retreive list of var indices in pvarlist for matching ones
      !@arguments
      integer, intent(out) :: F_idxvlist(:)       !# List of indices in pvarlist for var matching provided params
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
      character(len=256) :: flagstr
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
      if (present(F_flagstr)) flagstr = F_flagstr
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
         
         DO_PVARLIST: do iv = 1, npvarlist
            !# check bus
            match_L = (nbpath == 0)
            if (nbpath > 0) match_L = any(pvarlist(iv)%meta%bus == buses)
            if (.not.match_L) cycle
            !# check names
            IF_NAME: if (name /= ' ') then
               match_L = .false.
               select case (npath(in:in))
               case ('V')
                  match_L = (pvarlist(iv)%meta%vname(1:slen) == name)
               case ('I')
                  match_L = (pvarlist(iv)%meta%iname(1:slen) == name)
               case ('O')
                  match_L = (pvarlist(iv)%meta%oname(1:slen) == name)
               case ('S')
                  match_L = (pvarlist(iv)%meta%sname(1:slen) == name)
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
                  if (.not.any(flags(ik) == pvarlist(iv)%meta%flags(:))) then
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

         enddo DO_PVARLIST
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
   function phymem_find_var(F_varlist, F_name, F_npath, F_bpath, &
        F_quiet, F_shortmatch, F_flagstr) result(F_nmatch)
      implicit none
      !@objective Retreive list of var indices in pvarlist for matching ones
      !@arguments
      type(phyvar),      intent(out) :: F_varlist(:)        !# List of indices in pvarlist for var matching provided params
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
      character(len=256) :: flagstr
      logical :: shortmatch_L, quiet_L
      integer :: iv, idxvlist(PHY_MAXVARS), istat
      !---------------------------------------------------------------
      F_nmatch = PHY_ERROR
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
      
      shortmatch_L = .false.
      if (present(F_shortmatch)) shortmatch_L = F_shortmatch
      
      flagstr = ' '
      if (present(F_flagstr)) flagstr = F_flagstr

      F_nmatch = phymem_find_idxv(idxvlist, name, npath, bpath, quiet_L, shortmatch_L, flagstr)
      if (F_nmatch <= 0) return

      !TODO: should it be an error to limit the max number of returned values? it's not in find_idxv
!!$      if (F_nmatch > size(F_varlist)) then
!!$         if (.not.quiet_L) call msg(MSG_WARNING,'(phymem_find) F_varlist buffer overflow')
!!$         F_nmatch = PHY_ERROR
!!$         return
!!$      endif
      
      do iv = 1, F_nmatch
         if (iv > size(F_varlist)) then
            if (.not.quiet_L) call msg(MSG_WARNING,'(phymem_find) F_varlist buffer overflow')
            exit
         endif
         F_varlist(iv)%meta => pvarlist(idxvlist(iv))%meta
      enddo
      !---------------------------------------------------------------
      return
   end function phymem_find_var


   !/@*
   function phymem_getmeta(F_meta, F_idxv) result(F_istat)
      implicit none
      !@objective Return meta from pvarlist(idxv)
      !@arguments
      type(phymeta), pointer    :: F_meta  !# pvarlist(F_idxv)%meta
      integer,       intent(in) :: F_idxv  !# pvarlist index of the field
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
      if (F_idxv < 1 .or. F_idxv > npvarlist) then
         call msg(MSG_ERROR,'(phymem_get) Requested F_idxv out of range')
         return
      endif
      F_meta => pvarlist(F_idxv)%meta
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_getmeta

   
   !/@*
   function phymem_getptr1d(F_ptr, F_idxv, F_trnch) result(F_istat)
      implicit none
      !@objective Associate pointer 1d to pvarlist(idxv) for specified slice
      !@arguments
      real, pointer, contiguous :: F_ptr(:)  !# data(1:nikf)
      integer, intent(in) :: F_idxv   !# pvarlist index of the field
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
      if (F_idxv < 1 .or. F_idxv > npvarlist) then
         call msg(MSG_ERROR,'(phymem_get) Requested F_idxv out of range')
         return
      endif
      
      i0 = pvarlist(F_idxv)%meta%i0
      nikfm = pvarlist(F_idxv)%meta%size
      in = pvarlist(F_idxv)%meta%in
      F_ptr(1:nikfm) => pvarlist(F_idxv)%meta%bptr(i0:in,F_trnch)
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_getptr1d

   
   !/@*
   function phymem_getptr2d(F_ptr, F_idxv, F_trnch) result(F_istat)
      implicit none
      !@objective Associate pointer 2d to pvarlist(idxv) for specified slice
      !@arguments
      real, pointer, contiguous :: F_ptr(:,:)  !# data(1:ni, 1:nkf)
      integer, intent(in) :: F_idxv   !# pvarlist index of the field
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
      if (F_idxv < 1 .or. F_idxv > npvarlist) then
         call msg(MSG_ERROR,'(phymem_get) Requested F_idxv out of range')
         return
      endif

      i0 = pvarlist(F_idxv)%meta%i0
      ni = pvarlist(F_idxv)%meta%ni
      in = pvarlist(F_idxv)%meta%in
      nkf = pvarlist(F_idxv)%meta%nk * pvarlist(F_idxv)%meta%fmul * (pvarlist(F_idxv)%meta%mosaic+1)
      F_ptr(1:ni,1:nkf) => pvarlist(F_idxv)%meta%bptr(i0:in,F_trnch)
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_getptr2d

   
   !/@*
   function phymem_getptr3d(F_ptr, F_idxv, F_trnch) result(F_istat)
      implicit none
      !@objective Associate pointer 3d to pvarlist(idxv) for specified slice
      !@arguments
      real, pointer, contiguous :: F_ptr(:,:,:)  !# data(1:ni, 1:nk, 1:nf)
      integer, intent(in) :: F_idxv   !# pvarlist index of the field
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
      if (F_idxv < 1 .or. F_idxv > npvarlist) then
         call msg(MSG_ERROR,'(phymem_get) Requested F_idxv out of range')
         return
      endif

      i0 = pvarlist(F_idxv)%meta%i0
      ni = pvarlist(F_idxv)%meta%ni
      in = pvarlist(F_idxv)%meta%in
      nk = pvarlist(F_idxv)%meta%nk
      nfm = pvarlist(F_idxv)%meta%fmul * (pvarlist(F_idxv)%meta%mosaic+1)
      F_ptr(1:ni,1:nk,1:nfm) => pvarlist(F_idxv)%meta%bptr(i0:in,F_trnch)
      F_istat = PHY_OK
      !---------------------------------------------------------------
      return
   end function phymem_getptr3d
   
end module phymem
