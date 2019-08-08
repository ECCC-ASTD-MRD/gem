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

function phyfillbus(F_kount) result(F_istat)
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use clib_itf_mod, only: clib_tolower, clib_toupper
   use phygridmap, only: phy_lcl_ni, phy_lcl_nj, phy_lcl_i0, phy_lcl_in, phy_lcl_j0, phy_lcl_jn, phydim_nk
   use phy_typedef
   use phy_getmeta_mod, only: phy_getmeta
   use phyfoldmeta_mod, only: phyfold
   implicit none
!!!#include <arch_specific.hf>

   integer, intent(in) :: F_kount            !physics timestep number
   integer :: F_istat                        !Function result

   !@author  Michel Desgagne  -   summer 2013
   !@object  Transfer data to p_runlgt space

#include <msg.h>
#include <rmnlib_basics.hf>
#include <mu_gmm.hf>

   logical, parameter :: NOSHORTMATCH_L = .false.
   integer, parameter :: NVARMAX = 256
   integer, parameter :: MUST_INIT = 1
   character(len=*), parameter :: FLD_INIT(4) = (/ &
        'TDMASK', &
        'DLAT  ', &
        'DLON  ', &
        'DXDY  '  &
        /)

   integer, save :: nvars = 0
   type(phymeta), pointer, save :: metalist(:) => null()

   type(gmm_metadata) :: meta
   character(len=GMM_MAXNAMELENGTH) :: varname_S
   character(len=32) :: prefix_S, basename_S, time_S, ext_S
   integer :: i, k0, istat, err, ijkmin(3), ijkmax(3)
   real, pointer :: src2d(:,:), src3d(:,:,:)
   real, pointer :: src2d1(:,:), src3d1(:,:,:)
   !     ---------------------------------------------------------------
   F_istat = RMN_ERR

   ijkmin = 1
   ijkmax = (/phy_lcl_ni, phy_lcl_nj, 1/)

   ! Initialize specific entries for the permanent bus
   if (F_kount == 0) then
      istat = 0
      do i = 1, size(FLD_INIT)
         nullify(src2d)
         istat = min(gmm_get(FLD_INIT(i), src2d, meta), istat)
         if (associated(src2d)) then
            src2d1(1:,1:) => src2d(phy_lcl_i0:phy_lcl_in,phy_lcl_j0:phy_lcl_jn)
            istat = min( &
                 phyfold(src2d1,FLD_INIT(i),'P',ijkmin,ijkmax), &
                 istat)
         else
            call msg(MSG_ERROR,'(phyfillbus) missing GMM var: '//trim(FLD_INIT(i)))
            istat = RMN_ERR
         endif
      enddo
      if (.not.GMM_IS_OK(istat)) then
         call msg(MSG_ERROR,'(phyfillbus) cannot access required initializations')
         return
      endif
   endif
   F_istat = RMN_ERR

   ! Pull dynamics state into the bus
   istat = 0

   if (.not.associated(metalist)) then
      nvars = phy_getmeta(metalist, ' ', F_npath='V', F_bpath='D', &
           F_maxmeta=NVARMAX, F_shortmatch=NOSHORTMATCH_L)
   endif

   DO_VAR: do i= 1, nvars
      IF_INIT: if (metalist(i)%init == MUST_INIT) then

         call gmmx_name_parts(trim(metalist(i)%vname),prefix_S,basename_S,time_S,ext_S)
         varname_S = trim(prefix_S)//trim(basename_S)//trim(time_S)
         err = clib_toupper(varname_S)
         err = gmm_getmeta(trim(varname_S), meta)
         if (.not.RMN_IS_OK(err)) then
            err = clib_tolower(varname_S)
            err = gmm_getmeta(trim(varname_S), meta)
         endif
         if (.not.RMN_IS_OK(err)) cycle !#TODO: ERROR?

         nullify(src2d,src3d)
         IF_3D: if (meta%l(3)%high < 1) then
            err = gmm_get(trim(varname_S), src2d, meta)
            if (.not.RMN_IS_OK(err)) &
                 call msg(MSG_ERROR,'(phyfillbus) Problem getting GMM 2d pointer for: '//trim(varname_S))
            istat = min(err,istat)
            if (associated(src2d)) then
               ijkmax(3) = 1
               src2d1(1:,1:) => src2d(phy_lcl_i0:phy_lcl_in,phy_lcl_j0:phy_lcl_jn)
               err = phyfold(src2d1,trim(metalist(i)%vname),'D',ijkmin,ijkmax)
               if (.not.RMN_IS_OK(err)) &
                    call msg(MSG_ERROR,'(phyfillbus) Problem folding 2d pointer for: '//trim(varname_S))
               istat = min(err,istat)
            endif
         else !IF_3D
            err = gmm_get (trim(varname_S), src3d, meta)
            istat = min(err,istat)
            if (.not.RMN_IS_OK(err)) &
                 call msg(MSG_ERROR,'(phyfillbus) Problem getting GMM 3d pointer for: '//trim(varname_S))
            if (associated(src3d)) then
               k0 = 1
               ijkmax(3) = min(phydim_nk,meta%l(3)%high)
               !#TODO: not used anymore, remove?
!!$               if (metalist(i)%n(3) == 1) then
               if (metalist(i)%nk == 1) then
                  k0 = meta%l(3)%high
                  ijkmax(3) = 1
               endif
               src3d1(1:,1:,1:) => src3d(phy_lcl_i0:phy_lcl_in,phy_lcl_j0:phy_lcl_jn,k0:)
               err = phyfold(src3d1,trim(metalist(i)%vname),'D',ijkmin,ijkmax)
               if (.not.RMN_IS_OK(err)) &
                    call msg(MSG_ERROR,'(phyfillbus) Problem folding 3d pointer for: '//trim(varname_S))
               istat = min(err,istat)
            else
               call msg(MSG_ERROR,'(phyfillbus) Problem getting GMM pointer for: '//trim(varname_S))
               istat = RMN_ERR
            endif
         endif IF_3D

      endif IF_INIT
   end do DO_VAR
   if (.not.GMM_IS_OK(istat)) then
      call msg(MSG_ERROR,'(phyfillbus) cannot access required dynamics inputs')
      return
   endif
   F_istat = RMN_OK
   !     ---------------------------------------------------------------
   return
end function phyfillbus
