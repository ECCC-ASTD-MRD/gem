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

module phybudget
  use, intrinsic :: iso_fortran_env, only: REAL64
  use phy_status, only: PHY_OK, PHY_ERROR
  implicit none
  private

  ! Public procedures
  public :: pb_init             !Initialize physics budget package
  public :: pb_compute          !Compute physics energy budget terms
  public :: pb_residual         !Compute physics energy budget residual
  public :: pb_conserve         !Adjust tendencies to close physics energy budget

  ! Internal parameters
  integer, parameter :: SHORT_STRING=16
  character(len=SHORT_STRING), parameter :: INTTYPE_DEFAULT='pchip'
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function pb_init() result(F_istat)
    use energy_budget, only: eb_init
    use phygridmap, only: phydim_nk
    ! Initialize the physics budget package
    implicit none
    
    ! Arguments
    integer :: F_istat                                          !Return status (PHY_OK or PHY_ERROR)
    
    ! Basic setup
    F_istat = PHY_ERROR
    
    ! Initialize energy budget package
    if (eb_init(F_lv_temp_dependent=.true., &
         F_nearsfc_k=(phydim_nk-1)) /= PHY_OK) then
       call physeterror('phybudget::pb_init', &
            'Cannot initialize energy budget package')
       return
    endif

    ! End of subprogram
    F_istat = PHY_OK
  end function pb_init
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function pb_compute(F_enri, F_pwri, F_eni, F_pwi, F_dbus, F_pbus, F_vbus, &
       F_nk, F_inttype) result(F_istat)
    use phybus, only: pmoins, tplus, huplus, sigw
    use energy_budget, only: eb_en, eb_pw
    implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf>   
#include "phymkptr.hf"

    ! Arguments
    real, dimension(:), pointer, contiguous :: F_enri           !Integrated energy budget residual (m3/s2)
    real, dimension(:), pointer, contiguous :: F_pwri           !Integrated total water budget residual (Kg/m2)
    real(KIND=REAL64), dimension(:), intent(out) :: F_eni       !Integrated energy (m3/s2)
    real(KIND=REAL64), dimension(:), intent(out) :: F_pwi       !Integrated total water (Kg/m2)
    real, dimension(:), pointer, contiguous :: F_dbus           !Dynamics bus
    real, dimension(:), pointer, contiguous :: F_pbus           !Permanent bus
    real, dimension(:), pointer, contiguous :: F_vbus           !Volatile bus
    integer, intent(in) :: F_nk                                 !Number of levels (not including diagnostic)
    character(len=*), intent(in), optional :: F_inttype         !Integral type [INTTYPE_DEFAULT]
    integer :: F_istat                                          !Return status of subprogram
    
    !@Author   R. McTaggart-Cowan (2022)

    !@Revision

    !@Object
    !          Simplify and combine calls to energy budget module.

    !@Notes
    !          
    !*@/

    ! Internal variables
    integer :: ni, nkm1
    real, dimension(size(F_eni),F_nk) :: lwc, iwc
    real, dimension(:), pointer, contiguous :: zpmoins
    real, dimension(:,:), pointer, contiguous :: ztplus, zhuplus, &
         zsigw
    character(len=SHORT_STRING) :: inttype

    ! Basic setup
    F_istat = PHY_ERROR
    ni = size(F_eni)
    nkm1 = F_nk
        
    ! Early exit if no budgets are requested
    if (.not.associated(F_enri) .and. .not.associated(F_pwri)) then
       F_istat = PHY_OK
       return
    endif

    ! Handle optional arguments
    inttype = INTTYPE_DEFAULT
    if (present(F_inttype)) inttype = F_inttype

    ! Retrieve required bus entries
    MKPTR1D(zpmoins, pmoins, F_pbus)
    MKPTR2Dm1(ztplus, tplus, F_dbus)
    MKPTR2Dm1(zhuplus, huplus, F_dbus)
    MKPTR2Dm1(zsigw, sigw, F_dbus)

    ! Compute condensate contents
    if (pb_condensate(lwc, iwc, F_dbus, F_pbus, F_vbus) /= PHY_OK) then
       call physeterror('phybudget::pb_residual', &
            'Cannot compute condensate contents')
       return
    endif

    ! Compute integrated total water
    if (associated(F_pwri)) then
       if (eb_pw(F_pwi, zhuplus, lwc, zsigw, zpmoins, &
            nkm1, F_qi=iwc, F_inttype=inttype) /= PHY_OK) then
          call physeterror('phybudget::pb_compute', &
               'Cannot compute integrated total water')
          return
       endif
    else
       F_pwi = 0.
    endif

    ! Compute integrated energy
    if (associated(F_enri)) then
       if (eb_en(F_eni, ztplus, zhuplus, lwc, zsigw, zpmoins, &
            nkm1, F_qi=iwc, F_inttype=inttype) /= PHY_OK) then
          call physeterror('phybudget::pb_compute', &
               'Cannot compute integrated energy')
          return
       endif
    else
       F_eni = 0.
    endif

    ! End of subprogram
    F_istat = PHY_OK
  end function pb_compute

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function pb_conserve(F_cons, F_ttend, F_qtend, F_dbus, F_pbus, F_vbus, &
       F_dqc, F_dqi, F_rain, F_snow, F_shf, F_wvf, F_rad, F_inttype) result(F_istat)
    use phybus, only: pmoins, tplus, huplus, sigw
    use energy_budget, only: eb_conserve_en, eb_conserve_pw
    implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf>   
#include "phymkptr.hf"

    ! Arguments
    character(len=*), intent(in) :: F_cons                      !Conservation type
    real, dimension(:,:), intent(inout) :: F_ttend              !Temperature tendency (K/s)
    real, dimension(:,:), intent(inout) :: F_qtend              !Moisture tendency (kg/kg/s)
    real, dimension(:), pointer, contiguous :: F_dbus           !Dynamics bus
    real, dimension(:), pointer, contiguous :: F_pbus           !Permanent bus
    real, dimension(:), pointer, contiguous :: F_vbus           !Volatile bus
    real, dimension(:,:), intent(in), optional :: F_dqc         !Liquid condensate tendency (kg/kg/s)
    real, dimension(:,:), intent(in), optional :: F_dqi         !Solid condensate tendency (kg/kg/s)
    real, dimension(:), intent(in), optional :: F_rain          !Surface liquid precipitation flux (kg/m2/s) [0.]
    real, dimension(:), intent(in), optional :: F_snow          !Surface solid precipitation flux (kg/m2/s) [0.]
    real, dimension(:), intent(in), optional :: F_shf           !Surface turbulent sensible heat flux (W/m2) [0.]
    real, dimension(:), intent(in), optional :: F_wvf           !Surface turbulent water vapour flux (kg/m2/s) [0.]
    real, dimension(:), intent(in), optional :: F_rad           !Net radiation flux TOA-surface (W/m2) [0.]
    character(len=*), intent(in), optional :: F_inttype         !Integral type [INTTYPE_DEFAULT]
    integer :: F_istat                                          !Return status of subprogram
    
    !@Author   R. McTaggart-Cowan (2022)

    !@Revision

    !@Object
    !          Simplify and combine calls to energy budget module.

    !@Notes
    !          
    !*@/

    ! Internal variables
    integer :: ni, nkm1
    real, dimension(size(F_ttend,dim=1)) :: rain, snow, shf, wvf, rad
    real, dimension(size(F_ttend,dim=1),size(F_ttend,dim=2)) :: &
         lwc, iwc, dqc
    real, dimension(:), pointer, contiguous :: zpmoins
    real, dimension(:,:), pointer, contiguous :: ztplus, zhuplus, &
         zsigw
    character(len=SHORT_STRING) :: inttype

    ! Basic setup
    F_istat = PHY_ERROR
    ni = size(F_ttend,dim=1)
    nkm1 = size(F_ttend,dim=2)

    ! Early exit if no tendency corrections are requested
    if (F_cons /= 'TEND') then
       F_istat = PHY_OK
       return
    endif

    ! Handle optional arguments
    inttype = INTTYPE_DEFAULT
    if (present(F_inttype)) inttype = F_inttype
    
    ! Fill boundary conditions if not supplied
    if (present(F_dqc)) then
       dqc(:,:) = F_dqc(:,:)
    else
       dqc(:,:) = 0.
    endif
    if (present(F_rain)) then
       rain(:) = F_rain(:)
    else
       rain(:) = 0.
    endif
    if (present(F_snow)) then
       snow(:) = F_snow(:)
    else
       snow(:) = 0.
    endif
    if (present(F_shf)) then
       shf(:) = F_shf(:)
    else
       shf(:) = 0.
    endif
    if (present(F_wvf)) then
       wvf(:) = F_wvf(:)
    else
       wvf(:) = 0.
    endif
    if (present(F_rad)) then
       rad(:) = F_rad(:)
    else
       rad(:) = 0.
    endif

    ! Retrieve required bus entries
    MKPTR1D(zpmoins, pmoins, F_pbus)
    MKPTR2Dm1(ztplus, tplus, F_dbus)
    MKPTR2Dm1(zhuplus, huplus, F_dbus)
    MKPTR2Dm1(zsigw, sigw, F_dbus)

    ! Compute condensate contents
    if (pb_condensate(lwc, iwc, F_dbus, F_pbus, F_vbus) /= PHY_OK) then
       call physeterror('phybudget::pb_residual', &
            'Cannot compute condensate contents')
       return
    endif

    ! Budgets depend on splitting of liquid/solid condensate
    SOLID_TEND_SEPARATE: if (present(F_dqi)) then
           
       ! Apply humidity tendency correction for total water conservation
       if (eb_conserve_pw(F_qtend, F_qtend, ztplus, zhuplus, zsigw, &
            zpmoins, nkm1, F_dqc=dqc, F_dqi=F_dqi, F_rain=rain, &
            F_snow=snow, F_wvf=wvf, F_inttype=inttype) /= PHY_OK) then
          call physeterror('phybudget::pb_conserve', &
               'Cannot conserve integrated total water with dqi')
          return
       endif

       ! Apply temperature tendency correction for energy conservation
       if (eb_conserve_en(F_ttend, F_ttend, F_qtend, ztplus, zhuplus, &
            lwc, zsigw, zpmoins, nkm1, F_dqc=dqc, F_dqi=F_dqi, &
            F_qi=iwc, F_rain=rain, F_snow=snow, F_shf=shf, F_wvf=wvf, &
            F_rad=rad, F_inttype=inttype) /= PHY_OK) then
          call physeterror('phybudget::pb_conserve', &
               'Cannot conserve integrated energy with dqi')
          return
       endif

    else

       ! Apply humidity tendency correction for total water conservation
       if (eb_conserve_pw(F_qtend, F_qtend, ztplus, zhuplus, zsigw, &
            zpmoins, nkm1, F_dqc=dqc, F_rain=rain, F_snow=snow, &
            F_wvf=wvf, F_inttype=inttype) /= PHY_OK) then
          call physeterror('phybudget::pb_conserve', &
               'Cannot conserve integrated total water without dqi')
          return
       endif

       ! Apply temperature tendency correction for energy conservation
       if (eb_conserve_en(F_ttend, F_ttend, F_qtend, ztplus, zhuplus, &
            lwc, zsigw, zpmoins, nkm1, F_dqc=dqc, F_qi=iwc, &
            F_rain=rain, F_snow=snow, F_shf=shf, F_wvf=wvf, &
            F_rad=rad, F_inttype=inttype) /= PHY_OK) then
          call physeterror('phybudget::pb_conserve', &
               'Cannot conserve integrated energy without dqi')
          return
       endif

    endif SOLID_TEND_SEPARATE

    ! End of subprogram
    F_istat = PHY_OK
  end function pb_conserve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function pb_residual(F_enri, F_pwri, F_eni0, F_pwi0, F_dbus, &
       F_pbus, F_vbus, F_dt, F_nk, F_rain, F_snow, F_shf, &
       F_wvf, F_rad, F_inttype) result(F_istat)
    use phybus, only: pmoins, tplus, huplus, sigw
    use energy_budget, only: eb_en, eb_pw, eb_residual_en, &
         eb_residual_pw
    implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf>   
#include "phymkptr.hf"

    ! Arguments
    real, dimension(:), pointer, contiguous :: F_enri           !Integrated energy budget residual (m3/s2)
    real, dimension(:), pointer, contiguous :: F_pwri           !Integrated total water budget residual (m)
    real(KIND=REAL64), dimension(:), intent(in) :: F_eni0       !Initial integrated energy (m3/s2)
    real(KIND=REAL64), dimension(:), intent(in) :: F_pwi0       !Initial integrated total water (Kg/m2)
    real, dimension(:), pointer, contiguous :: F_dbus           !Dynamics bus
    real, dimension(:), pointer, contiguous :: F_pbus           !Permanent bus
    real, dimension(:), pointer, contiguous :: F_vbus           !Volatile bus
    real, intent(in) :: F_dt                                    !Time step (s)
    integer, intent(in) :: F_nk                                 !Number of levels (not including diagnostic)
    real, dimension(:), intent(in), optional :: F_rain          !Surface liquid precipitation flux (kg/m2/s) [0.]
    real, dimension(:), intent(in), optional :: F_snow          !Surface solid precipitation flux (kg/m2/s) [0.]
    real, dimension(:), intent(in), optional :: F_shf           !Surface turbulent sensible heat flux (W/m2) [0.]
    real, dimension(:), intent(in), optional :: F_wvf           !Surface turbulent water vapour flux (kg/m2/s) [0.]
    real, dimension(:), intent(in), optional :: F_rad           !Net radiation flux TOA-surface (W/m2) [0.]
    character(len=*), intent(in), optional :: F_inttype         !Integral type [INTTYPE_DEFAULT]
    integer :: F_istat                                          !Return status of subprogram
    
    !@Author   R. McTaggart-Cowan (2022)

    !@Revision

    !@Object
    !          Simplify and combine calls to energy budget module.

    !@Notes
    !          
    !*@/

    ! Internal variables
    integer :: ni, nkm1
    real, dimension(size(F_eni0)) :: rain, snow, shf, wvf, rad
    real, dimension(size(F_eni0),F_nk) :: lwc, iwc
    real, dimension(:), pointer, contiguous :: zpmoins
    real, dimension(:,:), pointer, contiguous :: ztplus, zhuplus, &
         zsigw
    real(kind=REAL64), dimension(size(F_eni0)) :: eni, pwi, enri, &
         pwri
    character(len=SHORT_STRING) :: inttype

    ! Basic setup
    F_istat = PHY_ERROR
    ni = size(F_eni0)
    nkm1 = F_nk

    ! Early exit if no budgets are requested
    if (.not.associated(F_enri) .and. .not.associated(F_pwri)) then
       F_istat = PHY_OK
       return
    endif

    ! Handle optional arguments
    inttype = INTTYPE_DEFAULT
    if (present(F_inttype)) inttype = F_inttype
    
    ! Fill boundary conditions if not supplied
    if (present(F_rain)) then
       rain(:) = F_rain(:)
    else
       rain(:) = 0.
    endif
    if (present(F_snow)) then
       snow(:) = F_snow(:)
    else
       snow(:) = 0.
    endif
    if (present(F_shf)) then
       shf(:) = F_shf(:)
    else
       shf(:) = 0.
    endif
    if (present(F_wvf)) then
       wvf(:) = F_wvf(:)
    else
       wvf(:) = 0.
    endif
    if (present(F_rad)) then
       rad(:) = F_rad(:)
    else
       rad(:) = 0.
    endif

    ! Retrieve required bus entries
    MKPTR1D(zpmoins, pmoins, F_pbus)
    MKPTR2Dm1(ztplus, tplus, F_dbus)
    MKPTR2Dm1(zhuplus, huplus, F_dbus)
    MKPTR2Dm1(zsigw, sigw, F_dbus)

    ! Compute condensate contents
    if (pb_condensate(lwc, iwc, F_dbus, F_pbus, F_vbus) /= PHY_OK) then
       call physeterror('phybudget::pb_residual', &
            'Cannot compute condensate contents')
       return
    endif

    ! Compute current energy budget state
    if (pb_compute(F_enri, F_pwri, eni, pwi, F_dbus, F_pbus, &
         F_vbus, nkm1, F_inttype=inttype) /= PHY_OK) then
       call physeterror('phybudget::pb_residual', &
            'Cannot compute current budget state')
       return
    endif
    
    ! Compute total water residual
    if (associated(F_pwri)) then
       if (eb_residual_pw(pwri, F_pwi0, pwi, ztplus, F_dt, nkm1, &
            F_rain=rain, F_snow=snow, F_wvf=wvf) /= PHY_OK) then
          call physeterror('phybudget::pb_residual', &
               'Cannot compute total water residual')
          return
       endif
       F_pwri(:) = real(pwri)
    endif

    ! Compute energy residual
    if (associated(F_enri)) then
       if (eb_residual_en(enri, F_eni0, eni, ztplus, zhuplus, lwc, &
            F_dt, nkm1, F_qi=iwc, F_rain=rain, F_snow=snow, &
            F_shf=shf, F_wvf=wvf, F_rad=rad) /= PHY_OK) then
          call physeterror('phybudget::pb_residual', &
               'Cannot compute energy residual')
          return
       endif
       F_enri(:) = real(enri)
    endif 
    
    ! End of subprogram
    F_istat = PHY_OK
  end function pb_residual

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function pb_condensate(lwc, iwc, F_dbus, F_pbus, F_vbus) result(F_istat)
    use phybus, only: qtbl
    use microphy_utils, only: mp_lwc, mp_iwc
    implicit none
#include "phymkptr.hf"

    ! Arguments
    real, dimension(:,:), intent(out) :: lwc                    !Liquid water content (kg/kg)
    real, dimension(:,:), intent(out) :: iwc                    !Ice water content (kg/kg)
    real, dimension(:), pointer, contiguous :: F_dbus           !Dynamics bus
    real, dimension(:), pointer, contiguous :: F_pbus           !Permanent bus
    real, dimension(:), pointer, contiguous :: F_vbus           !Volatile bus
    integer :: F_istat                                          !Return status of subprogram

    ! Internal variables
    integer :: ni,nkm1
    real, dimension(:,:), pointer :: zqtbl

    ! Basic setup
    F_istat = PHY_ERROR
    ni = size(lwc, dim=1)
    nkm1 = size(lwc, dim=2)

    ! Content from microphysical sources
    if (mp_lwc(lwc, F_dbus, F_pbus, F_vbus) /= PHY_OK .or. &
         mp_iwc(iwc, F_dbus, F_pbus, F_vbus) /= PHY_OK) then
       call physeterror('phybudget::pb_residual', &
            'Cannot compute condensate contents')
       return
    endif

    ! Content from other sources
    MKPTR2Dm1(zqtbl, qtbl, F_pbus)
    lwc(:,:) = lwc(:,:) + zqtbl(:,:)

    ! End of subprogram
    F_istat = PHY_OK
  end function pb_condensate
  
end module phybudget
