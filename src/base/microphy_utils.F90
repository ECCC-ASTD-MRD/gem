
module microphy_utils
  use phy_status, only: PHY_OK, PHY_ERROR
  implicit none
  private
  save
  
  ! Public procedures
  public :: mp_init            !Initialization before bus creation
  public :: mp_phybusinit      !Define bus requirements
  public :: mp_lwc             !Compute liquid water content
  public :: mp_iwc             !Compute ice water content
  public :: mp_post_init       !Initialization after bus creation
  public :: mp_todiffuse       !List of bus entries for vertical diffusion
  public :: mp_icefrac         !Compute ice fraction
  public :: nil_phybusinit     !Define bus requirements stcond=nil
  public :: nil_lwc            !Compute liquid water content stcond=nil
  public :: nil_iwc            !Compute ice water content stcond=nil

  ! Bus definitions
  abstract interface
     function phybusinit() result(F_istat)
       integer :: F_istat
     end function phybusinit
  end interface
  procedure(phybusinit), pointer :: mp_phybusinit => nil_phybusinit

  ! Total water mass calculation
  abstract interface
     function lwc(F_qltot, F_pvars, F_tminus) result(F_istat)
       use phymem, only: phyvar
       real, dimension(:,:), intent(out) :: F_qltot
       type(phyvar), pointer, contiguous :: F_pvars(:)  !All phy vars (meta + slab data)
       logical, intent(in), optional :: F_tminus        !Compute fields at time-minus [false]
       integer :: F_istat
     end function lwc
  end interface
  procedure(lwc), pointer :: mp_lwc => nil_lwc
  
  ! Total ice mass calculation
  abstract interface
     function iwc(F_qitot, F_pvars, F_tminus) result(F_istat)
       use phymem, only: phyvar
       real, dimension(:,:), intent(out) :: F_qitot
       type(phyvar), pointer, contiguous :: F_pvars(:)  !All phy vars (meta + slab data)
       logical, intent(in), optional :: F_tminus        !Compute fields at time-minus [false]
       integer :: F_istat
     end function iwc
  end interface
  procedure(iwc), pointer :: mp_iwc => nil_iwc

  ! Module parameters
  integer, parameter :: NMAX_DIFFUSE=256
  
  ! Module variables
  integer :: ndlist
  integer, dimension(NMAX_DIFFUSE), target :: dlist
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  ! Initialize procedure pointers and run scheme-specific bootstrap
  function mp_init(F_input_path) result(F_istat)
    use microphy_consun
    use microphy_s2
    use microphy_p3, P3_STATUS_OK=>STATUS_OK
    use microphy_p3v3, only: &
         p3v3_init => p3_init, p3v3_phybusinit => p3_phybusinit, &
         p3v3_lwc => p3_lwc, p3v3_iwc => p3_iwc
    use microphy_my2
    use microphy_kessler
    use phy_options, only: stcond, p3_ncat, p3_trplmomi, p3_liqFrac
    implicit none
    character(len=*), intent(in) :: F_input_path  !Directory containing initializing data
    integer :: F_istat                            !Return status (PHY_OK on success)

    ! Internal variables
    integer :: istat
    character(len=1024) :: missing_list

    ! Initialization
    F_istat = PHY_ERROR

    ! Map procedures based on selected microphysics scheme
    istat = PHY_ERROR
    select case (stcond)
    case ('CONSUN')
       istat = PHY_OK
       mp_phybusinit => consun_phybusinit
       mp_lwc => consun_lwc
       mp_iwc => consun_iwc
    case ('S2')
       istat = PHY_OK
       mp_phybusinit => s2_phybusinit
       mp_lwc => s2_lwc
       mp_iwc => s2_iwc
    case ('MP_P3')
       call p3_init(F_input_path, p3_ncat, p3_trplmomi, p3_liqFrac, stat=istat)
       if (istat == P3_STATUS_OK) istat = PHY_OK
       mp_phybusinit => p3_phybusinit
       mp_lwc => p3_lwc
       mp_iwc => p3_iwc
    case ('MP_P3V3')
       call p3v3_init(F_input_path, p3_ncat, stat=istat)
       if (istat == P3_STATUS_OK) istat = PHY_OK
       mp_phybusinit => p3v3_phybusinit
       mp_lwc => p3v3_lwc
       mp_iwc => p3v3_iwc
    case ('MP_MY2')
       istat = PHY_OK
       mp_phybusinit => my2_phybusinit
       mp_lwc => my2_lwc
       mp_iwc => my2_iwc
    case ('KESSLER')
       istat = PHY_OK
       mp_phybusinit => kessler_phybusinit
       mp_lwc => kessler_lwc
       mp_iwc => kessler_iwc
    case ('NIL')
       istat = PHY_OK
    case DEFAULT
       call physeterror('microphy_utils::mp_init', &
            'Unknown microphysics scheme '//trim(stcond))
       return
    end select

    ! Check for initialization success
    if (istat /= PHY_OK) then
       call physeterror('microphy_utils::mp_init', &
            'Failed to initialize '//trim(stcond))
       return
    endif

    ! Confirm that all procedures have been assigned
    missing_list = ''
    if (.not. associated(mp_phybusinit)) &
         missing_list = trim(missing_list)//' mp_phybusinit'
    if (.not. associated(mp_lwc)) &
         missing_list = trim(missing_list)//' mp_lwc'
    if (.not. associated(mp_iwc)) &
         missing_list = trim(missing_list)//' mp_iwc'
    if (len_trim(missing_list) /= 0) then
       call physeterror('microphy_utils::mp_init', &
            'Failed to initialize microphysics API components for ', &
            trim(stcond)//': '//trim(missing_list))
       return
    endif

    ! Successful completion
    F_istat = PHY_OK
  end function mp_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Perform initialization after bus creation
  function mp_post_init() result(F_istat)
    use phy_options, only: fluvert
    use phymem, only: phymeta, phymem_find, phymem_getmeta
    use phy_status, only: PHY_OK
    implicit none
    integer :: F_istat                            !Return status (PHY_OK on success)
    
#include <rmn/msg.h>
    
    ! Internal parameters
    character(len=*), dimension(6), parameter :: DIFF_FLAGS=(/ &
         'CLOUD_  ', &
         'RAIN_   ', &
         'SNOW_   ', &
         'GRAUPEL_', &
         'HAIL_   ', &
         'ICE_    '/)
    
    ! Internal variables
    integer :: nmatch, i
    integer, dimension(NMAX_DIFFUSE) :: idxvlist
    character(len=1024) :: msg_S
    type(phymeta), pointer :: pmeta

    ! Initialization
    F_istat = PHY_ERROR

    ! Create list of bus indexes for diffusion
    ndlist = 0
    if (fluvert == 'RPNINT') then
       do i=1,size(DIFF_FLAGS)
          nmatch = phymem_find(idxvlist, F_name=':P', F_endmatch=.true., &
               F_flagstr=trim(DIFF_FLAGS(i)), F_shortflag=.true., F_quiet=.true.)          
          dlist(ndlist+1:ndlist+nmatch) = idxvlist(1:nmatch)
          ndlist = ndlist+nmatch
       enddo
       msg_S = ''
       if (ndlist > 0) then
          do i=1,ndlist
             if (phymem_getmeta(pmeta, dlist(i)) /= PHY_OK) then
                call physeterror('microphy_utils::mp_init', &
                     'Cannot retrieve diffusion metadata')
                return
             endif
             msg_S = trim(msg_S)//' '//trim(pmeta%vname)
          enddo
       endif
       call msg(MSG_INFO, '(mp_init) For PBL diffusion: '//trim(msg_S))
    endif

    ! Successful completion
    F_istat = PHY_OK
  end function mp_post_init
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  ! Provide list of microphysics fields for PBL diffusion
  function mp_todiffuse(F_difflist) result(F_istat)
    implicit none
    integer, dimension(:), pointer :: F_difflist !List of fields to diffuse
    integer :: F_istat                          !Function return status
    F_difflist => dlist(1:ndlist)
    F_istat = PHY_OK
  end function mp_todiffuse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mp_icefrac(F_if, F_tt, F_lwc, F_iwc, F_ni, F_nkm1, F_difdt) &
       result(F_istat)
    use tdpack, only: TCDK
    implicit none

    !@Arguments
    integer, intent(in) :: F_ni                                 !horizontal dimension
    integer, intent(in) :: F_nkm1                               !vertical dimension
    real, dimension(F_ni,F_nkm1), intent(in) :: F_tt            !dry air temperature (K)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_lwc           !liquid water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(in) :: F_iwc           !ice water content (kg/kg)
    real, dimension(F_ni,F_nkm1), intent(out) :: F_if           !ice fraction
    real, dimension(F_ni,F_nkm1), intent(out), optional :: F_difdt !ice fraction derivative wrt temperature (/K)
    integer :: F_istat                                          !Function return status

    !@Object Compute cloud ice fraction

    ! Local parameters from Eq. 2 of Hu et al. (2010; JGR 10.1029/2009JD012384)
    real, parameter :: C0 = 5.3608, C1=0.4025, C2=0.08387, C3=0.007182, &
         C4=2.39e-4, C5=2.87e-6
    
    ! Local variables
    integer :: i, k
    real :: ptmid, dptmid, tmid, lwc, iwc

    ! Initialization
    F_istat = PHY_ERROR
    
    ! Compute ice fraction
    do k=1,F_nkm1
       do i=1,F_ni
          lwc = max(F_lwc(i,k), 0.)
          iwc = max(F_iwc(i,k), 0.)
          F_if(i,k) = iwc / max(lwc + iwc, tiny(iwc))
       enddo
    enddo

    ! Estimate temperature sensitivity of the ice fraction (with a minimum value
    ! on ptmid needed to avoid an overflow on the exp(-ptmid))
    if (present(F_difdt)) then
       do k=1,F_nkm1
          do i=1,F_ni             
             tmid = F_tt(i,k) - TCDK
             ptmid = max(C0 + tmid * (C1 + tmid * (C2 + tmid * (C3 + tmid * (C4 + tmid * C5)))), -40.)
             dptmid = C1 + tmid * (2*C2 + tmid * (3*C3 + tmid * (4*C4 + tmid * 5*C5)))
             F_difdt(i,k) = -exp(-ptmid) / (1. + exp(-ptmid))**2 * dptmid             
          enddo
       enddo
    endif

    ! Successful completion
    F_istat = PHY_OK
  end function mp_icefrac
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
  ! Define bus requirements
  function nil_phybusinit() result(F_istat)
     implicit none
     integer :: F_istat                          !Function return status
     F_istat = PHY_OK
  end function nil_phybusinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  ! Compute total water mass
  function nil_lwc(F_qltot, F_pvars, F_tminus) result(F_istat)
     use phymem, only: phyvar
     ! use phybusidx
     implicit none
     real, dimension(:,:), intent(out) :: F_qltot        !Total water mass (kg/kg)
     type(phyvar), pointer, contiguous :: F_pvars(:)     !All phy vars (meta + slab data)
     logical, intent(in), optional :: F_tminus           !Compute fields at time-minus [false]
     integer :: F_istat                                  !Return status
     F_qltot(:,:) = 0.
     F_istat = PHY_OK
     return
  end function nil_lwc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  ! Compute total ice mass
  function nil_iwc(F_qitot, F_pvars, F_tminus) result(F_istat)
     use phymem, only: phyvar
     ! use phybusidx
     use phy_status, only: PHY_OK, PHY_ERROR
     implicit none
     real, dimension(:,:), intent(out) :: F_qitot        !Total ice mass (kg/kg)
     type(phyvar), pointer, contiguous :: F_pvars(:)     !All phy vars (meta + slab data)
     logical, intent(in), optional :: F_tminus           !Compute fields at time-minus [false]
     integer :: F_istat                                  !Return status
     F_qitot(:,:) = 0.
     F_istat = PHY_OK
     return
  end function nil_iwc

  
end module microphy_utils
