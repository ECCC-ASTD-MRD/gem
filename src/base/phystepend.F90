module phystepend
  implicit none
  private

  public :: phystepend1                 !Post-step physics operations

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine phystepend1(F_pvars, F_ni, F_nk)
    use debug_mod, only: init2nan
    use microphy_utils, only: mp_lwc, mp_iwc
    use phy_status, only: phy_error_L, PHY_OK
    use phybusidx, except1=>lwc, except2=>iwc
    use phymem, only: phyvar
    use prep_cw, only: prep_cw3
    use sgspdf, only: sgspdf_compute

    implicit none

    ! Argument declarations
    integer, intent(in) :: F_ni                         !Horizontal dimension
    integer, intent(in) :: F_nk                         !Vertical dimension
    type(phyvar), pointer, contiguous :: F_pvars(:)     !Physics buses

    !@Object Perform post-step physics operations

    ! Common declarations
#include "phymkptr.hf"

    ! Local variables
    integer :: nkm1
    real, dimension(F_ni, F_nk-1) :: lwc, iwc
    real, dimension(:), pointer, contiguous :: zp0, zdxdy, zh
    real, dimension(:,:), pointer, contiguous :: zsigmas, ztplus, &
         zhuplus, zpri, zzn, zzd, zgztherm, zsigt, zmrk2
    real, dimension(:,:,:), pointer, contiguous :: zvcoef

    ! Retrieve bus pointers
    MKPTR1D(zdxdy, dxdy, F_pvars)
    MKPTR1D(zh, h, F_pvars)
    MKPTR1D(zp0, p0_plus, F_pvars)

    MKPTR2D(zgztherm, gztherm, F_pvars)
    MKPTR2D(zhuplus, huplus, F_pvars)
    MKPTR2D(zmrk2, mrk2, F_pvars)
    MKPTR2D(zpri, pri, F_pvars)
    MKPTR2D(zsigmas, sigmas, F_pvars)
    MKPTR2D(zsigt, sigt, F_pvars)
    MKPTR2D(ztplus, tplus, F_pvars)
    MKPTR2D(zzd, zd, F_pvars)
    
    MKPTR3D(zvcoef, vcoef, F_pvars)

    if (znplus > 0) then
       MKPTR2D(zzn, znplus, F_pvars)
    else
       MKPTR2D(zzn, zn, F_pvars)
    endif

    ! Basic initialization
    call init2nan(lwc, iwc)
    nkm1 = F_nk-1
    
    ! Compute cloud water contents for radiation
    call prep_cw3(F_pvars, F_ni, F_nk)
    if (phy_error_L) return

    ! Compute condensate masses
    if (mp_lwc(lwc, F_pvars) /= PHY_OK .or. &
         mp_iwc(iwc, F_pvars) /= PHY_OK) then
       call physeterror('phystepend', 'Cannot compute condensate contents')
       return
    endif
    
    ! Diagnose subgrid-scale variability
    call sgspdf_compute(zsigmas, ztplus, zhuplus, lwc, iwc, zpri, &
         zzn, zzd, zgztherm, zsigt, zp0, zh, &
         zdxdy, zmrk2, zvcoef, F_ni, nkm1)
    if (phy_error_L) return

    ! End of subprogram
    return
  end subroutine phystepend1

end module phystepend
