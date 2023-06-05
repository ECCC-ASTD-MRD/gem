!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

module microphy_utils
  use phy_status, only: PHY_OK, PHY_ERROR
  implicit none
  private
  save
  
  ! Public procedures
  public :: mp_init            !Initialize microphysics scheme
  public :: mp_phybusinit      !Define bus requirements
  public :: mp_lwc             !Compute liquid water content
  public :: mp_iwc             !Compute ice water content
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
     function lwc(F_qltot, F_dbus, F_pbus, F_vbus) result(F_istat)
       real, dimension(:,:), intent(out) :: F_qltot
       real, dimension(:), pointer, contiguous :: F_dbus
       real, dimension(:), pointer, contiguous :: F_pbus
       real, dimension(:), pointer, contiguous :: F_vbus
       integer :: F_istat
     end function lwc
  end interface
  procedure(lwc), pointer :: mp_lwc => nil_lwc
  
  ! Total ice mass calculation
  abstract interface
     function iwc(F_qitot, F_dbus, F_pbus, F_vbus) result(F_istat)
       real, dimension(:,:), intent(out) :: F_qitot
       real, dimension(:), pointer, contiguous :: F_dbus
       real, dimension(:), pointer, contiguous :: F_pbus
       real, dimension(:), pointer, contiguous :: F_vbus
       integer :: F_istat
     end function iwc
  end interface
  procedure(iwc), pointer :: mp_iwc => nil_iwc
 
contains

  ! Initialize procedure pointers and run scheme-specific bootstrap
  function mp_init(F_input_path) result(F_istat)
    use microphy_consun
    use microphy_p3, P3_STATUS_OK=>STATUS_OK
    use microphy_my2
    use microphy_kessler
    use phy_options, only: stcond, p3_ncat
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
    case ('MP_P3')
       call p3_init(F_input_path, p3_ncat, stat=istat)
       if (istat == P3_STATUS_OK) istat = PHY_OK
       mp_phybusinit => p3_phybusinit
       mp_lwc => p3_lwc
       mp_iwc => p3_iwc
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
  
  
  ! Define bus requirements
  function nil_phybusinit() result(F_istat)
     implicit none
     integer :: F_istat                          !Function return status
     F_istat = PHY_OK
  end function nil_phybusinit

  ! Compute total water mass
  function nil_lwc(F_qltot, F_dbus, F_pbus, F_vbus) result(F_istat)
     use phybus
     implicit none
     real, dimension(:,:), intent(out) :: F_qltot        !Total water mass (kg/kg)
     real, dimension(:), pointer, contiguous :: F_dbus   !Dynamics bus
     real, dimension(:), pointer, contiguous :: F_pbus   !Permanent bus
     real, dimension(:), pointer, contiguous :: F_vbus   !Volatile bus
     integer :: F_istat                                  !Return status
     F_qltot(:,:) = 0.
     F_istat = PHY_OK
     return
  end function nil_lwc

  ! Compute total ice mass
  function nil_iwc(F_qitot, F_dbus, F_pbus, F_vbus) result(F_istat)
     use phybus
     use phy_status, only: PHY_OK, PHY_ERROR
     implicit none
     real, dimension(:,:), intent(out) :: F_qitot        !Total ice mass (kg/kg)
     real, dimension(:), pointer, contiguous :: F_dbus   !Dynamics bus
     real, dimension(:), pointer, contiguous :: F_pbus   !Permanent bus
     real, dimension(:), pointer, contiguous :: F_vbus   !Volatile bus
     integer :: F_istat                                  !Return status
     F_qitot(:,:) = 0.
     F_istat = PHY_OK
     return
  end function nil_iwc

  
end module microphy_utils
