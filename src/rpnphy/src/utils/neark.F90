!/@*
module neark
   implicit none
   public
  
   integer, parameter :: NEARK_MONO_ASC = 1
   integer, parameter :: NEARK_MONO_DES = -1

   interface neark_mono
      module procedure neark_mono2d1
      module procedure neark_mono2d2
   end interface neark_mono

contains

   !/@*
   function neark_mono2d1(data, val, ni, nk, icol, iorder) result(kindex)
      implicit none
      !@arguments
      integer, intent(in) :: ni, nk, icol
      real, intent(in) :: data(ni,nk), val
      integer, intent(in) :: iorder   ! +1: ascending, -1: descending
      !@return      
      integer :: kindex
      !*@/      
      integer :: klow, khigh, kmid
      !----------------------------------------------------------------------
      klow = 1
      khigh = nk
      do while (klow < khigh)
         kmid = (klow + khigh) / 2
         ! if (val < data(icol,kmid)) then  !Ascending
         ! if (val > data(icol,kmid)) then  !Descending
         if (iorder * (val - data(icol,kmid)) < 0) then
            khigh = kmid
         else
            klow = kmid + 1
         end if
      end do

      if (klow == 1) then
         kindex = 1
      else if (klow > nk) then
         kindex = nk
      else
         if (abs(data(icol,klow) - val) < abs(data(icol,klow-1) - val)) then
            kindex = klow
         else
            kindex = klow - 1
         end if
      end if
      !----------------------------------------------------------------------
      return
   end function neark_mono2d1

   
   !/@*
   function neark_mono2d2(data, values, ni, nk) result(kindex)
      implicit none
      !@arguments
      integer, intent(in) :: ni, nk
      real, intent(in) :: data(ni,nk), values(ni)
      !@return
      integer :: kindex(ni)
      !*@/
      integer :: icol, iorder
      !----------------------------------------------------------------------
      iorder = merge(1, -1, data(1,nk) > data(1,1))  ! +1: ascending, -1: descending
      !#TODO: OMP parallel with potential GPU offload
      do icol = 1,ni
         kindex(icol) = neark_mono2d1(data, values(icol), ni, nk, icol, iorder)
      enddo
!!$      print *,'d',kindex
      !----------------------------------------------------------------------
      return
   end function neark_mono2d2

   
   !/@*
   function neark_dp_opti(sig,ps,delp,ni,nk) result(kindex)
      implicit none
      !@arguments
      integer, intent(in) :: ni,nk                !Slab dimensions
      real, dimension(ni,nk), intent(in) :: sig   !Sigma level values for profile
      real, dimension(ni), intent(in) :: ps       !Surface pressure (Pa)
      real, intent(in) :: delp                    !Pressure depth for level search (Pa)
      !@return
      integer, dimension(ni) :: kindex        !Level index closes to "delp" above the surface
      !*@/
      real, dimension(ni) :: sigtarget
      !----------------------------------------------------------------------
      sigtarget = 1.-delp/ps
      kindex = neark_mono2d2(sig, sigtarget, ni, nk)
      !----------------------------------------------------------------------
      return
   end function neark_dp_opti

   
   !/@*
   function neark_dp_orig(sig,ps,delp,ni,nk) result(kindex)
      implicit none
      !@arguments
      integer, intent(in) :: ni,nk                !Slab dimensions
      real, dimension(ni,nk), intent(in) :: sig   !Sigma level values for profile
      real, dimension(ni), intent(in) :: ps       !Surface pressure (Pa)
      real, intent(in) :: delp                    !Pressure depth for level search (Pa)
      !@return
      integer, dimension(ni) :: kindex            !Level index closes to "delp" above the surface
      !*@/
      integer :: i,k
      real, dimension(ni) :: sigtarget
      !----------------------------------------------------------------------
      sigtarget = 1.-delp/ps
      kindex = nk
      do i=1,ni
         k = 2      
         do while (kindex(i) == nk .and. k <= nk)
            if (sig(i,k) > sigtarget(i) .and. sig(i,k-1) <= sigtarget(i)) kindex(i) = k
            k = k+1
         enddo
      enddo
      !----------------------------------------------------------------------
      return
   end function neark_dp_orig


end module neark
