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

subroutine coherence3(ni, trnch)
   use tdpack_const, only: TRPL, RAUW
   use sfc_options
   use sfcbus_mod
   use svs_configs
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer ni, trnch

   !@Author Stephane Belair (September 1999)
   !@Revision
   ! 001      see version 5.6.0 for previous history
   ! 002      add ADJ_I0_SNOW key/cfg in front of I0 adjustement due to snow
   !          presence on the ground ( M. Abrahamowicz , July 2015)
   ! 003      M. Mackay (2018/2022)  Add code for CSLM
   !
   !@Object Assure the coherence between the different masks
   !         (i.e., MG, GLSEA, and GLACIER) and the surface fields
   !
   !@Arguments
   !             - Input -
   ! ni          horizontal length of a slab

   include "isbapar.cdk"
   include "sfcinput.cdk"

   integer i, k
   real fsnow

   real, pointer, dimension(:) :: zalveg,  zcveg,  zgamveg,  zglacier,  zglsea,  zicedp,  zlai,  zmg,  zrgl,  zrootdp,  zsnoal, zsnoden, zsnoma,  zsnoro,  zstomr,  zvegfrac,  zwsnow,  zwveg
!!$      real, pointer, dimension(:) :: zsdepth

   real, pointer, dimension(:,:) :: zclay, zisoil, zsand, zsnodp, ztglacier, ztsoil, zwsoil,ztpsoil
   ! SVS
   real, pointer, dimension(:) :: zsnodpl, zsnval, zsnvden, zsnvdp, zsnvma, zsnvro, zvegh, zvegl, zwsnv


#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni) => busptr(vd%NAME2%i)%ptr(:,trnch)
#define MKPTR2D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => busptr(vd%NAME2%i)%ptr(:,trnch)

   MKPTR1D(zalveg,   alveg)
   MKPTR1D(zcveg,    cveg)
   MKPTR1D(zgamveg,  gamveg)
   MKPTR1D(zglacier, glacier)
   MKPTR1D(zglsea,   glsea)
   MKPTR1D(zicedp,   icedp)
   !MKPTR1D(zisoil,   isoil)
   MKPTR1D(zlai,     lai)
   MKPTR1D(zmg,      mg)
   MKPTR1D(zrgl,     rgl)
   MKPTR1D(zrootdp,  rootdp)
   MKPTR1D(zsnoal,   snoal)
   MKPTR1D(zsnoden,  snoden)
   MKPTR1D(zsnodpl,  snodpl)
   MKPTR1D(zsnoma,   snoma)
   MKPTR1D(zsnoro,   snoro)
   MKPTR1D(zsnval,   snval)
   MKPTR1D(zsnvden,  snvden)
   MKPTR1D(zsnvdp,   snvdp)
   MKPTR1D(zsnvma,   snvma)
   MKPTR1D(zsnvro,   snvro)
   MKPTR1D(zstomr,   stomr)
   MKPTR1D(zvegfrac, vegfrac)
   MKPTR1D(zvegh,    vegh)
   MKPTR1D(zvegl,    vegl)
   MKPTR1D(zwsnow,   wsnow)
   MKPTR1D(zwsnv,    wsnv)
   MKPTR1D(zwveg,    wveg)

   MKPTR2D(zclay,    clay)
   MKPTR2D(zisoil,   isoil)
   MKPTR2D(zsand,    sand)
   MKPTR2D(zsnodp,   snodp)
   MKPTR2D(ztglacier,tglacier)
   MKPTR2D(ztsoil,   tsoil)
   MKPTR2D(ztpsoil,   tpsoil)
   MKPTR2D(zwsoil,   wsoil)
 
   !***************************************************************
   ! Coherence tests on the mask "mg"
   !***************************************************************
 
   ! Over water surfaces (mg = 0), many fields can be put to 0
   ! -- for esthetic look of output only --
 
   NEW_MG_MASK: if (any('mg' == phyinread_list_s(1:phyinread_n))) then
!VDIR NODEP
      do i=1,ni
         if (zmg(i).lt.critmask) then

            zmg      (i)              = 0.0
            zglacier (i)              = 0.0
            zsnodp   (i,indx_soil)    = 0.0
            zsnodp   (i,indx_glacier) = 0.0
            if (schmurb  == 'TEB') zsnodp   (i,indx_urb ) = 0.0
            if (schmlake /= 'NIL') zsnodp   (i,indx_lake) = 0.0

         end if
      end do

      IF_ISBA: if (schmsol == 'ISBA') then

!VDIR NODEP
         do i=1,ni
            ! Over water for ESTHETIC PURPOSES
            if (zmg(i).lt.critmask) then
               zwsoil   (i,1)    = 1.0
               zwsoil   (i,2)    = 1.0
               zwveg    (i)      = 0.0
               zisoil   (i,1)    = 0.0
               zrootdp  (i)      = 0.0
               zsnoden  (i)      = 100.00
               zstomr   (i)      = 0.0
               zlai     (i)      = 0.0
               zvegfrac (i)      = 0.0
               zsand    (i,1)    = 0.0
               zclay    (i,1)    = 0.0

            end if
         end do
 
         ! Over land, we need to make sure that
         ! soil texture is reasonable and that there
         ! is water in the soil
         ! -- not for esthetics --
 
!VDIR NODEP
         do i=1,ni
            if (zmg(i).ge.critmask) then

               ! If no sand and clay components are found
               ! over points where MG > critmask
               ! attribute to these points characteristics
               ! of typical loamy soils

               if (zsand(i,1)+zclay(i,1).lt.critexture) then
                  zsand(i,1) = 35.
                  zclay(i,1) = 35.
               end if

               zwsoil(i,1) = max(zwsoil(i,1),1.e-7)
               zwsoil(i,2) = max(zwsoil(i,2),1.e-7)

               zsand(i,1) = max(1.,zsand(i,1))
               zclay(i,1) = max(1.,zclay(i,1))

               zalveg   (i) = max( zalveg   (i) , 0.12 )
               zrootdp  (i) = max( zrootdp  (i) , 0.5  )
               zstomr   (i) = max( zstomr   (i) , 40.  )
               zcveg    (i) = max( zcveg    (i) , 1.e-5)
               zrgl     (i) = max( zrgl     (i) , 30.  )
               zlai     (i) = max( zlai     (i) , 0.0  )
               zvegfrac (i) = max( zvegfrac (i) , 0.0  )
               zgamveg  (i) = max( zgamveg  (i) , 0.0  )
               zwveg    (i) = max( 0.0 , min( zwveg    (i) , &
                    0.2 * zvegfrac(i) * zlai(i) ) )

            end if
         end do

      end if IF_ISBA

      IF_SVS:  if (schmsol.eq.'SVS') then
         do i=1,ni
            if (zmg(i).lt.critmask) then
               ! OVER WATER, FOR ESTHETIC PURPOSE ONLY
               do k=1,nl_svs
                  zwsoil(i,k)   = 1.0
                  zisoil(i,k)   = 0.0
                  ztpsoil(i,k)   = -1.0
               enddo
               zwveg    (i)      = 0.0
               zrootdp  (i)      = 0.0
               zvegfrac (i)      = 0.0
               zvegh    (i)      = 0.0
               zvegl    (i)      = 0.0
               zsnodpl(i) = 0.0
               zsnoma(i)  = 0.0
               zwsnow(i)  = 0.0
               zsnvdp(i)  = 0.0
               zsnvma(i)  = 0.0
               zwsnv(i)   = 0.0
            endif
         enddo      

         if (lsoil_freezing_svs1) then
            do i=1,ni
               if (zmg(i).ge.critmask) then
                  ! OVER LAND, *NOT* ESTHETIC PURPOSE ONLY
                  do k=1,nl_svs
                     ! Make sure tsoil=TRPL if have ice and tsoil > TRPL
                     !if ( zisoil(i,k).gt.0.0 .and. ztpsoil(i,k).gt.TRPL ) &
                     if ( zisoil(i,k).gt.0.0 .and. (ztpsoil(i,k)-TRPL) .gt. EPSILON_SVS_TK)  &
                          ztpsoil(i,k)=TRPL
                  enddo
               endif
            enddo
         endif  

      endif IF_SVS

   endif NEW_MG_MASK


   ! Over surface of land only (mg=1) there is no sea ice
   ! -- for esthetics of output only --

   if (any('mg' == phyinread_list_s(1:phyinread_n))) then
!VDIR NODEP
      do i=1,ni
         if(zmg(i) > 1.-critmask) then
            zmg   (i) = 1.
            zglsea(i) = 0.0
            zicedp(i) = 0.0
         end if
      end do
   endif

!**************************************************************
!               Coherence tests on the mask "glsea"
!**************************************************************

   NEW_GLSEA_MASK: if (any('glsea0' == phyinread_list_s(1:phyinread_n))) then
!VDIR NODEP
      do i=1,ni
         if (zglsea (i).lt.critmask) then
            zglsea (i)          = 0.0
            zicedp (i)          = 0.0
            zsnodp (i,indx_ice) = 0.0
         else
            zicedp (i) = max( zicedp(i) , minicedp )
         end if
      end do


      ! the following situation can only occur is
      ! leadfrac = 0., which is not usually the case

!VDIR NODEP
      do i=1,ni
         if (zglsea(i) > 1.-critmask) then
            zglsea(i) = 1.
         end if
      end do
   endif NEW_GLSEA_MASK

!**************************************************************
!               Coherence tests on the mask "glacier"
!**************************************************************

   NEW_GL_MASK: if (any('glacier' == phyinread_list_s(1:phyinread_n))) then

!VDIR NODEP
      do i=1,ni
         if (zglacier(i).ge.critmask) then
            ! deep glacier temperature cannot be greater than 0 c
            ztglacier(i,2) = min(ztglacier(i,2),trpl)
         end if
      end do

!vdir nodep
      do i=1,ni
         if (zglacier(i).lt.critmask) then
            zglacier(i) = 0.0
         end if
      end do

      ! For outputs esthetics only

!VDIR NODEP
      do i=1,ni
        if (zglacier(i) > 1.-critmask) then
            zglacier(i)   = 1.0
        end if
      end do

      IF_ISBA2: if (schmsol == 'ISBA') then

!VDIR NODEP
         do i=1,ni
            if (zglacier(i) > 1.-critmask) then
               zwsoil  (i,1) = 1.0
               zwsoil  (i,2) = 1.0
               zwveg    (i)      = 0.0
               zisoil   (i,1)    = 0.0
               zrootdp  (i)      = 0.0
               zsnoden  (i)      = 100.00
               zstomr   (i)      = 0.0
               zlai     (i)      = 0.0
               zvegfrac (i)      = 0.0
               zsand    (i,1)    = 0.0
               zclay    (i,1)    = 0.0

            end if
         end do

      end if IF_ISBA2

      IF_SVS2: if (schmsol == 'SVS') then

!VDIR NODEP
         do i=1,ni
            if (zglacier(i) > 1.-critmask) then

               do k=1,nl_svs
                  zwsoil(i,k)   = 1.0
                  zisoil(i,k)   = 0.0
                  ztpsoil(i,k)  = -1.0
               enddo

               zwveg    (i)      = 0.0
               zrootdp  (i)      = 0.0
               zvegfrac (i)      = 0.0
               zvegh    (i)      = 0.0
               zvegl    (i)      = 0.0
               zsnodpl(i) = 0.0
               zsnoma(i)  = 0.0
               zwsnow(i)  = 0.0
               zsnvdp(i)  = 0.0
               zsnvma(i)  = 0.0
               zwsnv(i)   = 0.0
            end if
         end do

      end if IF_SVS2
   endif NEW_GL_MASK


   !**************************************************************
   ! Coherence tests on snow fields
   !**************************************************************
   NEW_SD_MASK: if (any('snodp' == phyinread_list_s(1:phyinread_n))) then

!VDIR NODEP
      do i=1,ni
          if (zsnodp(i,indx_soil).lt.critsnow) then
              zsnodp(i,indx_soil) = 0.0
          end if
      end do

      IF_ISBA3: if (schmsol == 'ISBA') then
!VDIR NODEP
         do i=1,ni
            if (zsnodp(i,indx_soil).lt.critsnow) then
               zsnoma(i) = 0.0
               zwsnow(i) = 0.0
               zsnoro(i) = rhosdef
               zsnoal(i) = ansmax
            else if (adj_i0_snow) then
               ! Push land sfc temp. >0C towards 0C if snow is present on the ground
               fsnow = min( zsnodp(i,indx_soil) / 0.10 , 1. )
               if (ztsoil(i,1) > 273.16) &
                    ztsoil(i,1) = (1.-fsnow)*ztsoil(i,1) + fsnow*273.16
               if (ztsoil(i,2) > 273.16) &
                    ztsoil(i,2) = (1.-fsnow)*ztsoil(i,2) + fsnow*273.16
            end if
         end do
      end if IF_ISBA3

      IF_SVS3: if (schmsol.EQ.'SVS') then
!VDIR NODEP

         ! Calculate density
         !
         ! CAREFUL HERE about the units:
         ! "snoro"/"snvro" is the relative density of snow,
         !     i.e., rho_ice / rho_water (no units)
         ! "snoma"/"snvma" is the snow water equivalent in
         !     mm (i.e., kg / m2)
         ! "snoden"/"snvden" : snow density in kg/m3

         do i=1,ni

            if (zsnoma(i).lt.critsnowmass.or.zsnodpl(i).eq.0.0) then
               zsnodpl(i) = 0.0
               zsnoma(i)  = 0.0
               zwsnow(i)  = 0.0
               zsnoro(i)  = rhosdef
               zsnoden(i) = rhosdef * rauw
               zsnoal(i)  = ansmax
            else
               zsnoro(i)  = min(  max(100.,zsnoden(i)) / rauw  , 0.9 )
            endif
               
            if (zsnvma(i).lt.critsnowmass.or.zsnvdp(i).eq.0.0) then
               zsnvdp(i)  = 0.0
               zsnvma(i)  = 0.0
               zwsnv(i)   = 0.0
               zsnvro(i)  = rhosdef
               zsnvden(i) = rhosdef*rauw
               zsnval(i)  = ansmax
            else
               zsnvro(i)  =  min(  max(100.,zsnvden(i)) / rauw  , 0.9 )
            endif

         enddo
      endif IF_SVS3

   endif NEW_SD_MASK

   return
end subroutine coherence3
