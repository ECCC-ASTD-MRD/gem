
module sigmalev
   implicit none
   private
   public :: sigmalev3

contains

   !/@*
  function sigmalev3(vcoef, se, s, st, n, nk) result(F_istat)
    use phy_options, only: fluvert
     implicit none
!!!#include <arch_specific.hf>
      !@Arguments
      integer, intent(in)  :: n, nk          !# array dims
      real,    intent(in)  :: s(n,nk)        !# sigma momentum levels
      real,    intent(in)  :: st(n,nk)       !# sigma thermo levels
      real,    intent(out) :: vcoef(n,nk,2)  !# vinterp coef
      real,    intent(out) :: se(n,nk)       !# sigma energy levels
      !@Return
      integer :: F_istat
      !@Author  L. Spacek (Dec 2007)
      !@Object
      ! Compute sigma coordinates of energy levels and
      ! linear interpolation coefficients for temperature/humidity
      ! interpolation to momentum, thermo and energy levels.
      !*@/
#include <rmnlib_basics.hf>
#include <rmn/msg.h>
      !----------------------------------------------------------------

      if (fluvert == 'RPNINT') then
         se(1:n,1:nk-1)  = s(1:n,1:nk-1)
         se(1:n,nk) = 1.
      else
         se(1:n,1:nk-2)  = st(1:n,1:nk-2)
         se(1:n,nk-1:nk) = 1.
      endif

      call mweights(vcoef(:,:,1), s, st, n, nk, nk-1)
      vcoef(:,nk,1) = vcoef(:,nk-1,1)  !#Note this change vcoef stats, avoid undef
      call tweights(vcoef(:,:,2), s, st, n, nk, nk)
      
      F_istat = RMN_OK
      !----------------------------------------------------------------
      return
   end function sigmalev3


   !/@*
   subroutine mweights(atq2m, sigm, sigt, ni, nk, nkscope)
      implicit none
!!!#include <arch_specific.hf>
      !@arguments
      integer, intent(in)  :: ni,nk        !# Array dims
      integer, intent(in)  :: nkscope      !# nk work span
      real,    intent(in)  :: sigm(ni,nk)  !# sigma momentum levels
      real,    intent(in)  :: sigt(ni,nk)  !# sigma thermo levels
      real,    intent(out) :: atq2m(ni,nk) !# linear vinterp coef
      !@author l. spacek (dec 2007)
      !@object
      !  Compute the linear coefficients for interpolation
      !  of temperature/humidity to momentum levels
      !*@/
      integer :: k
      real :: wrk(ni,nkscope-1)
      !----------------------------------------------------------------
      wrk(:,1:nkscope-1) = 1. / (sigt(:,2:nkscope) - sigt(:,1:nkscope-1))
      atq2m(:,1) = 0.
      do k = 2, nkscope
         atq2m(:,k) = (sigt(:,k) - sigm(:,k)) * wrk(:,k-1)
      enddo
      !----------------------------------------------------------------
      return
   end subroutine mweights


   !/@*
   subroutine tweights(atq2t, sigm, sigt, ni, nk, nkscope)
      implicit none
!!!#include <arch_specific.hf>
      !@arguments
      integer, intent(in)  :: ni,nk        !# Array dims
      integer, intent(in)  :: nkscope      !# nk work span
      real,    intent(in)  :: sigm(ni,nk)  !# sigma momentum levels
      real,    intent(in)  :: sigt(ni,nk)  !# sigma thermo levels
      real,    intent(out) :: atq2t(ni,nk) !# linear vinterp coef
      !@author  l. spacek (dec 2007)
      !@object
      !  Compute the linear coefficients for interpolation
      !  of temperature/humidity to thermo levels
      !*@/
      integer :: k
      real :: wrk(ni,nkscope-1)
      !----------------------------------------------------------------
      wrk(:,1:nkscope-1) = 1. / (sigm(:,2:nkscope) - sigm(:,1:nkscope-1))
      do k = 1, nkscope-1
         atq2t(:,k) = (sigt(:,k) - sigm(:,k)) * wrk(:,k)
      enddo
      atq2t(:,nkscope) = 0.0
      !----------------------------------------------------------------
      return
   end subroutine tweights

end module sigmalev
