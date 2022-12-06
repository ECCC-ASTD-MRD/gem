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

module phyfillbus
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use rmn_gmm
   use clib_itf_mod, only: clib_tolower, clib_toupper
   use phygridmap, only: phy_lcl_ni, phy_lcl_nj, phy_lcl_i0, phy_lcl_in, phy_lcl_j0, phy_lcl_jn, phydim_nk
   use phyfold, only: phyfoldmeta1, phyfold1
   use phymem, only: pbuslist, PHY_DBUSIDX
   implicit none

   private
   public :: phyfillbus1
   
!!!#include <arch_specific.hf>

contains
   
function phyfillbus1(F_kount) result(F_istat)
   implicit none
   integer, intent(in) :: F_kount            !physics timestep number
   integer :: F_istat                        !Function result

   !@author  Michel Desgagne  -   summer 2013
   !@object  Transfer data to p_runlgt space

#include <rmn/msg.h>
#include <rmnlib_basics.hf>

   integer, parameter :: MUST_INIT = 1
   character(len=*), parameter :: FLD_INIT(4) = (/ &
        'TDMASK', &
        'DLAT  ', &
        'DLON  ', &
        'DXDY  '  &
        /)

   type(gmm_metadata) :: gmeta
   character(len=GMM_MAXNAMELENGTH) :: varname_S
   character(len=32) :: prefix_S, basename_S, time_S, ext_S
   integer :: iv, k0, istat, err, ijkmin(3), ijkmax(3)
   real, pointer :: src2d(:,:), src3d(:,:,:)
   real, pointer :: src2d1(:,:), src3d1(:,:,:)
   !     ---------------------------------------------------------------
   F_istat = RMN_ERR

   ijkmin = 1
   ijkmax = (/phy_lcl_ni, phy_lcl_nj, 1/)

   ! Initialize specific entries for the permanent bus
   if (F_kount == 0) then
      istat = 0
      do iv = 1, size(FLD_INIT)
         nullify(src2d)
         istat = min(gmm_get(FLD_INIT(iv), src2d, gmeta), istat)
         if (associated(src2d)) then
            src2d1(1:,1:) => src2d(phy_lcl_i0:phy_lcl_in,phy_lcl_j0:phy_lcl_jn)
            istat = min( &
                 phyfold1(src2d1, FLD_INIT(iv), 'P', ijkmin, ijkmax), &
                 istat)
         else
            call msg(MSG_ERROR, '(phyfillbus) missing GMM var: '//trim(FLD_INIT(iv)))
            istat = RMN_ERR
         endif
      enddo
      if (.not.GMM_IS_OK(istat)) then
         call msg(MSG_ERROR, '(phyfillbus) cannot access required initializations')
         return
      endif
   endif
   F_istat = RMN_ERR

   ! Pull dynamics state into the bus
   istat = 0
   
   DO_VAR: do iv= 1, pbuslist(PHY_DBUSIDX)%nvars
      IF_INIT: if (pbuslist(PHY_DBUSIDX)%meta(iv)%init == MUST_INIT) then

         call gmmx_name_parts(trim(pbuslist(PHY_DBUSIDX)%meta(iv)%vname), prefix_S, basename_S, time_S, ext_S)
         varname_S = trim(prefix_S)//trim(basename_S)//trim(time_S)
         err = clib_toupper(varname_S)
         err = gmm_getmeta(trim(varname_S), gmeta)
         if (.not.RMN_IS_OK(err)) then
            err = clib_tolower(varname_S)
            err = gmm_getmeta(trim(varname_S), gmeta)
         endif
         if (.not.RMN_IS_OK(err)) then
            call msg(MSG_WARNING, '(phyfillbus) Dyn bus var not found in GMM: '//trim(pbuslist(PHY_DBUSIDX)%meta(iv)%vname))
            cycle
         endif

         nullify(src2d, src3d)
         IF_3D: if (gmeta%l(3)%high < 1) then
            err = gmm_get(trim(varname_S), src2d, gmeta)
            if (.not.RMN_IS_OK(err)) &
                 call msg(MSG_ERROR, '(phyfillbus) Problem getting GMM 2d pointer for: '//trim(varname_S))
            istat = min(err,  istat)
            if (associated(src2d)) then
               ijkmax(3) = 1
               src2d1(1:,1:) => src2d(phy_lcl_i0:phy_lcl_in,phy_lcl_j0:phy_lcl_jn)
               err = phyfoldmeta1(src2d1, ijkmin, ijkmax, pbuslist(PHY_DBUSIDX)%meta(iv))
               if (.not.RMN_IS_OK(err)) &
                    call msg(MSG_ERROR, '(phyfillbus) Problem folding 2d pointer for: '//trim(varname_S))
               istat = min(err, istat)
            endif
         else !IF_3D
            err = gmm_get(trim(varname_S), src3d, gmeta)
            istat = min(err, istat)
            if (.not.RMN_IS_OK(err)) &
                 call msg(MSG_ERROR, '(phyfillbus) Problem getting GMM 3d pointer for: '//trim(varname_S))
            if (associated(src3d)) then
               k0 = 1
               ijkmax(3) = min(phydim_nk, gmeta%l(3)%high)
               if (pbuslist(PHY_DBUSIDX)%meta(iv)%nk == 1) then
                  k0 = gmeta%l(3)%high
                  ijkmax(3) = 1
               endif
               src3d1(1:,1:,1:) => src3d(phy_lcl_i0:phy_lcl_in,phy_lcl_j0:phy_lcl_jn,k0:)
               err = phyfoldmeta1(src3d1, ijkmin, ijkmax, pbuslist(PHY_DBUSIDX)%meta(iv))
               if (.not.RMN_IS_OK(err)) &
                    call msg(MSG_ERROR, '(phyfillbus) Problem folding 3d pointer for: '//trim(varname_S))
               istat = min(err, istat)
            else
               call msg(MSG_ERROR, '(phyfillbus) Problem getting GMM pointer for: '//trim(varname_S))
               istat = RMN_ERR
            endif
         endif IF_3D

      endif IF_INIT
   end do DO_VAR
   if (.not.GMM_IS_OK(istat)) then
      call msg(MSG_ERROR, '(phyfillbus) cannot access required dynamics inputs')
      return
   endif
   F_istat = RMN_OK
   !     ---------------------------------------------------------------
   return
end function phyfillbus1

end module phyfillbus
