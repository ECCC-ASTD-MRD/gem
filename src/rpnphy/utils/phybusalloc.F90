!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------

module phybusalloc
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use gesdictmod, only: gesdictlock
   use phy_getmeta_mod, only: phy_getmeta
   use phy_options, only: debug_mem_L
   use phy_typedef, only: phymeta
   use phygetmetaplus_mod, only: phymetaplus, phygetmetaplus
   use phy_status, only: phy_init_ctrl, PHY_CTRL_INI_OK, PHY_NONE
   use debug_mod, only: assert_not_naninf
   private
   public :: phybusalloc1

contains

   function phybusalloc1(p_nj, entbus, perbus, dynbus, volbus) result(F_istat)
      implicit none
!!!#include <arch_specific.hf>

      integer, intent(in) :: p_nj
      real, pointer, dimension(:,:) :: entbus, dynbus, perbus, volbus
      integer :: F_istat

#include <mu_gmm.hf>
#include <msg.h>
#include <rmnlib_basics.hf>
      include "buses.cdk"

      logical,parameter:: SHORTMATCH_L = .true.
      integer,parameter:: MYMAX = 2048

      character(len=1)  :: bus_S
      integer :: gmmstat, initval, ivar, nvars, istat
      type(gmm_metadata) :: meta_busent, meta_busper, meta_busdyn, meta_busvol
      type(phymetaplus) :: meta_m
      type(phymeta), pointer :: metalist(:)
      !---------------------------------------------------------------
      F_istat = RMN_ERR

      call gesdictlock()

      esp_busent = entpar(enttop,BUSPAR_I0) + entpar(enttop,BUSPAR_NIK) - 1
      esp_busper = perpar(pertop,BUSPAR_I0) + perpar(pertop,BUSPAR_NIK) - 1
      esp_busdyn = dynpar(dyntop,BUSPAR_I0) + dynpar(dyntop,BUSPAR_NIK) - 1
      esp_busvol = volpar(voltop,BUSPAR_I0) + volpar(voltop,BUSPAR_NIK) - 1

      call gmm_build_meta2D(meta_busent, 1,esp_busent,0,0,esp_busent, 1,p_nj,0,0,p_nj, 0,GMM_NULL_FLAGS )
      call gmm_build_meta2D(meta_busper, 1,esp_busper,0,0,esp_busper, 1,p_nj,0,0,p_nj, 0,GMM_NULL_FLAGS )
      call gmm_build_meta2D(meta_busdyn, 1,esp_busdyn,0,0,esp_busdyn, 1,p_nj,0,0,p_nj, 0,GMM_NULL_FLAGS )
      call gmm_build_meta2D(meta_busvol, 1,esp_busvol,0,0,esp_busvol, 1,p_nj,0,0,p_nj, 0,GMM_NULL_FLAGS )

      nullify(entbus, perbus, dynbus, volbus)

      initval = GMM_FLAG_IZER
      if (debug_mem_L) initval = GMM_FLAG_INAN
      gmmstat = gmm_create('BUSENT_3d', entbus, meta_busent, initval)
      gmmstat = gmm_create('BUSPER_3d', perbus, meta_busper, initval+GMM_FLAG_RSTR)
      gmmstat = gmm_create('BUSDYN_3d', dynbus, meta_busdyn, initval)
      gmmstat = gmm_create('BUSVOL_3d', volbus, meta_busvol, initval)

      gmmstat = gmm_get('BUSENT_3d',entbus)
      gmmstat = gmm_get('BUSPER_3d',perbus)
      gmmstat = gmm_get('BUSDYN_3d',dynbus)
      gmmstat = gmm_get('BUSVOL_3d',volbus)

      if (.not.(associated(entbus).and.associated(perbus).and.associated(dynbus).and.associated(volbus))) then
         call msg(MSG_ERROR,'(phybusalloc) Problem allocating physics bus pointers')
         return
      endif

      IF_DEBUG: if (debug_mem_L) then

         !# Reset Per bus var to zero at initial time
         phy_init_ctrl = PHY_CTRL_INI_OK
         bus_S = 'P'
         nullify(metalist)
         nvars = phy_getmeta(metalist, '', F_npath='V', F_bpath=bus_S, &
              F_maxmeta=MYMAX, F_shortmatch=SHORTMATCH_L)
         do ivar = 1, nvars
            istat = phygetmetaplus(meta_m, metalist(ivar)%vname, F_npath='V', &
                 F_bpath=bus_S, F_quiet=.true., F_shortmatch=.false.)
            if (.not.RMN_IS_OK(assert_not_naninf(meta_m%vptr))) then
               meta_m%vptr(:,:) = 0.
            endif
         enddo
         deallocate(metalist, stat=istat)
 
      end if IF_DEBUG

      F_istat = RMN_OK
      !---------------------------------------------------------------
      return
   end function phybusalloc1

end module phybusalloc
