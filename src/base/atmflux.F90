module atmflux
  use integrals, only: int_profile, INT_OK, INT_TYPE_LINEAR
  use tdpack, only: GRAV
  implicit none
  private
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

  ! Public API
  public :: atmflux4                                            !Compute flux profile from tendencies

  
contains

  subroutine atmflux4(F_flux, F_tend, F_sig, F_p0, F_ni, F_nk, F_topflux, F_type)
    ! Compute flux profile by integrating tendencies from the surface

    ! Arguments
    integer, intent(in) :: F_ni, F_nk                           !Slab dimensions
    real, dimension(F_ni,F_nk), intent(in) :: F_tend            !Tendency (field units/s)
    real, dimension(F_ni,F_nk), intent(in) :: F_sig             !Sigma levels of tendencies
    real, dimension(F_ni), intent(in) :: F_p0                   !Surface pressure (Pa)
    real, dimension(F_ni), intent(in), optional :: F_topflux    !Top boundary flux [0] (field units)
    integer, intent(in), optional :: F_type                     !Type of integration strategy to use ('pchip',['linear'],'step')
    real, dimension(F_ni,F_nk), intent(out) :: F_flux           !Field flux on staggered flux-levels

    ! Local declarations
    integer :: k
    real, dimension(F_ni) :: ptop
    real, dimension(F_ni,F_nk) :: tendi, presi, fluxi
    integer :: myType
    
    ! Set default values
    myType = INT_TYPE_LINEAR
    if (present(F_type)) myType = F_type

    ! Invert slab
    do k=1,F_nk
       presi(:,F_nk-(k-1)) = F_sig(:,k)*F_p0(:)
       tendi(:,F_nk-(k-1)) = F_tend(:,k) / GRAV  !Integrate as dp/g
    enddo
    ptop = presi(:,F_nk)
    
    ! Integrate to compute fluxes
    if (int_profile(fluxi, tendi, presi, ptop, F_p0, myType) /= INT_OK) then
       call physeterror('atmflux::atmflux4', 'Error in integral calculation')
       return
    endif

    ! Invert flux
    do k=1,F_nk
       F_flux(:,k) = fluxi(:,F_nk-(k-1))
    enddo

    ! Add boundary fluxes if provided
    if (present(F_topflux)) then
       do k=1,F_nk
          F_flux(:,k) = F_topflux(:)
       enddo
    endif

    ! Successful completion
  end subroutine atmflux4

end module atmflux
