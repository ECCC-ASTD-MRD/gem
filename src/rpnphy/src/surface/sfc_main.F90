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

module sfc_main
   private
   public :: sfc_main2

contains

!/@*
function sfc_main2(trnch, kount, dt, ni, nk) result(F_istat)
   use tdpack, only: DELTA, RGASD
   use phy_status, only: phy_error_L
   use phygridmap, only: ijdrv_phy
   use sfc_options
   use sfcbus_mod
   use sfclayer, only: sl_adjust,SL_OK
   use cpl_itf, only: cpl_update
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
   !@Object Main interface subroutine for the surface processes
   !@Arguments
   !          - Input/Output -
   ! D        dynamic             bus
   ! F        permanent variables bus
   ! V        volatile (output)   bus
   !          - Input -
   ! DSIZ     dimension of D
   ! FSIZ     dimension of F
   ! VSIZ     dimension of V
   !          - Input -
   ! TRNCH    row number
   ! KOUNT    timestep number
   ! DT       length of timestep
   ! NI       number of elements processed in the horizontal
   !          (could be a subset of M in future versions)
   !          horizontal dimensions of fields
   ! NK       vertical dimension

   logical, parameter :: CP_TO_SFCBUS = .true.
   logical, parameter :: CP_FROM_SFCBUS = .false.

   integer :: trnch,kount,ni, nk

   real :: dt, rhoair

   integer :: F_istat
   character(len=1024) :: msg_S

   !@Author B. Bilodeau (Dec 1998)
   !@Revisions
   ! 001      J. Mailhot  (Jul 2000) - Changes to add MOISTKE option (fluvert=MOISTKE)
   ! 002      B. Bilodeau (Mar 2003) - Convert aggregated value of SNODP to cm
   ! 003      Y. Delage (Aug 2003) - ILMO truncated between -100 and +100
   ! 004      M. Roch and B. Bilodeau  (Jan 2002) - Eliminate unit conversion for snodp
   ! 005      Y. Delage (Oct 2003) - ILMO truncated between -10 and +10
   ! 006      B. Bilodeau (Feb 2004) - Move call to fisimp3 from surface to phyexe1
   ! 007      B. Bilodeau (Feb 2004) - Revised logic to facilitate the
   !                                   addition of new types of surface
   ! 008      B. Bilodeau (Jun 2004) - Move preliminary calculations to phyexe1
   ! 009      Y. Delage   (Aug 2003) - Add call to CLASS300
   ! 010      B. Bilodeau and A. Lemonsu (2004) - Connect new scheme: TEB
   ! 011      D. Talbot   (winter 2005) - Limit surface inversions (tdiaglim)
   ! 012      A. Lemonsu  (Jun 2005) - Add call to town for the urban model
   ! 013      G. Balsamo  (Dec 2005) - Add call to offline cloud cover
   ! 014      B. Bilodeau (Jan 2006) - Add mosaic capability for CLASS
   ! 015      B. Bilodeau (Jul 2006) - Remove tests on string "SCHMURB"
   ! 016      A-M. Leduc (Nov 2007)  - Call surf_precip3 for BOURGE3D
   ! 017      L. Spacek   (Aug 2008) - za changes to ztsl
   ! 018      P. Vaillancourt (Jan 2009)  - Creation of new variable tnolim
   ! 019      J. Toviessi (July 2009) - added modifications for radslope
   ! 020      L. Spacek   (Nov 2011)  - insert calculations of tve, za etc.
   !                                    call to calz, lin_kdif_sim1
   ! 021      K. Winger,M. Mackay   (Feb 2017/Sep 2022)  - Add call to lake and river (M.A.) models
   !*@/
#include <rmn/msg.h>
   include "sfcinput.cdk"

   logical :: do_glaciers, do_ice, do_urb, do_lake, do_river
   integer :: i, k, sommet
   integer :: ni_soil,  ni_glacier,  ni_water,  ni_ice, ni_urb, ni_lake, &
        ni_river
   integer :: siz_soil, siz_glacier, siz_water, siz_ice, siz_urb, siz_lake, &
        siz_river 
   real :: lemin, lemax, mask

   !     les poids
   !     les "rangs" (dans l'espace de la grille complete)
   !     les "rangs" (dans chacun des 5 espaces
   !     correspondant aux 5 types de surfaces)

   integer,dimension(ni,indx_max) :: rangs
   real   ,dimension(ni,indx_max) :: poids
   integer,dimension(ni)            :: rg_soil, rg_water, rg_ice,  &
        rg_glacier, rg_urb, rg_lake, rg_river
   integer,dimension(2,ni)          :: lcl_indx
   integer,dimension(nvarsurf)      :: ptr_soil, ptr_water, ptr_ice, &
        ptr_glacier, ptr_urb, ptr_lake, ptr_river
   real,   dimension(surfesptot*ni) :: bus_soil, bus_water, bus_ice, &
        bus_glacier, bus_urb, bus_lake, bus_river
   !
   real, pointer, dimension(:)      :: zdtdiag, zmg, zfvapliq, zfvapliqaf, &
        zglacier, zglsea, zpmoins, zpplus, ztdiag, ztnolim, zurban, zztsl, &
        zqdiag, zudiag, zvdiag,zicedp,ztwater, zqdiagstn, ztdiagstn, &
        zudiagstn, zvdiagstn, zqdiagstnv, ztdiagstnv, zudiagstnv, zvdiagstnv, &
        zlakefr, zriverfr, zgridarea, zlakd, zlakearea
   real, pointer, dimension(:)    :: ztke, zhdpth, zlkiceh, zsniceh, zexpw, zdtemp, zdelu, zgred, zrhomix, ztsed, zroficeh
   real, pointer, dimension(:)    :: zsnol, zrhosnol, ztsnowl, zalbsnol, zwsnowl, zlst, zhlaksil, zficl
   real, pointer, dimension(:,:)    :: poids_out, zfvap, zilmo, zrunofftot, &
        ztmoins, ztplus, &
        zhuplus,zuplus,zvplus,zsnodp, &
        zqdiagtyp, ztdiagtyp, zudiagtyp, zvdiagtyp, zqdiagtypv, ztdiagtypv, &
        zudiagtypv, zvdiagtypv, ztddiagtyp, ztddiagtypv
   real, pointer, dimension(:,:)  :: ztlak, zlatflw, zwatflow
   real, pointer, dimension(:)    :: zrofinlak, zlstd, zlstf, zlfxi, zlfxo
   !     ---------------------------------------------------------------
   F_istat = RMN_ERR

   call msg_toall(MSG_DEBUG, 'sfc_main [BEGIN]')


#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni) => busptr(vd%NAME2%i)%ptr(:,trnch)
#define MKPTR2D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => busptr(vd%NAME2%i)%ptr(:,trnch)

   MKPTR1D(zdtdiag, dtdiag)
   MKPTR1D(zmg, mg)
   MKPTR1D(zfvapliq, fvapliq)
   MKPTR1D(zfvapliqaf, fvapliqaf) 
   MKPTR1D(zglacier, glacier)
   MKPTR1D(zglsea, glsea)
   MKPTR1D(zicedp, icedp)
   MKPTR1D(ztwater, twater)
   MKPTR1D(zpmoins, pmoins)
   MKPTR1D(zpplus, pplus)
   MKPTR1D(ztdiag, tdiag)
   MKPTR1D(ztdiagstn, tdiagstn)
   MKPTR1D(ztdiagstnv, tdiagstnv)
   MKPTR1D(zqdiag, qdiag)
   MKPTR1D(zqdiagstn, qdiagstn)
   MKPTR1D(zqdiagstnv, qdiagstnv)
   MKPTR1D(zudiag, udiag)
   MKPTR1D(zudiagstn, udiagstn)
   MKPTR1D(zudiagstnv, udiagstnv)
   MKPTR1D(zvdiag, vdiag)
   MKPTR1D(zvdiagstn, vdiagstn)
   MKPTR1D(zvdiagstnv, vdiagstnv)
   MKPTR1D(ztnolim, tnolim)
   MKPTR1D(zurban, urban)
   MKPTR1D(zztsl, ztsl)
   MKPTR1D(zlakefr, lakefr)
   MKPTR1D(zgridarea, gridarea)
   MKPTR1D(zlakearea, lakearea)
   MKPTR1D(zlakd, lakd)
   MKPTR1D(zriverfr, riverfr)

   MKPTR1D(ztke, tke)
   MKPTR1D(zhlaksil, hlaksil)
   MKPTR1D(zrofinlak, rofinlak)
   MKPTR1D(zlfxi, lfxi)
   MKPTR1D(zlfxo, lfxo)
   MKPTR1D(zlstd, lstd)
   MKPTR1D(zlstf, lstf)
   MKPTR1D(zlst, lst)
   MKPTR1D(zhdpth, hdpth)
   MKPTR1D(zlkiceh, lkiceh)
   MKPTR1D(zsniceh, sniceh)
   MKPTR1D(zexpw, expw)
   MKPTR1D(zdtemp, dtemp)
   MKPTR1D(zdelu, delu)
   MKPTR1D(zgred, gred)
   MKPTR1D(zrhomix, rhomix)
   MKPTR1D(ztsed, tsed)
   MKPTR1D(zroficeh, roficeh)
   MKPTR1D(zsnol,snol)
   MKPTR1D(zficl,ficl)
   MKPTR1D(zrhosnol,rhosnol)
   MKPTR1D(ztsnowl,tsnowl)
   MKPTR1D(zalbsnol,albsnol)
   MKPTR1D(zwsnowl,wsnowl)

   MKPTR2D(poids_out, sfcwgt)
   MKPTR2D(zfvap, fvap)
   MKPTR2D(zilmo, ilmo)
   MKPTR2D(zqdiagtyp, qdiagtyp)
   MKPTR2D(zqdiagtypv, qdiagtypv)
   MKPTR2D(zrunofftot, runofftot)
   MKPTR2D(zsnodp, snodp)
   MKPTR2D(ztddiagtyp, tddiagtyp)
   MKPTR2D(ztddiagtypv, tddiagtypv)
   MKPTR2D(ztdiagtyp, tdiagtyp)
   MKPTR2D(ztdiagtypv, tdiagtypv)
   MKPTR2D(ztmoins, tmoins)
   MKPTR2D(ztplus, tplus)
   MKPTR2D(zhuplus, huplus)
   MKPTR2D(zudiagtyp, udiagtyp)
   MKPTR2D(zudiagtypv, udiagtypv)
   MKPTR2D(zuplus, uplus)
   MKPTR2D(zvdiagtyp, vdiagtyp)
   MKPTR2D(zvdiagtypv, vdiagtypv)
   MKPTR2D(zvplus, vplus)

   MKPTR2D(ztlak, tlak)
   MKPTR2D(zlatflw, latflw)
   MKPTR2D(zwatflow, watflow)


   ! Update coupling fields GL, TM, SD and I8

   if (cplocn) then
      call cpl_update (zglsea (1:ni), 'GLI' , ijdrv_phy(1:2,1:ni,trnch:trnch), ni)
      call cpl_update (ztwater(1:ni), 'TMO' , ijdrv_phy(1:2,1:ni,trnch:trnch), ni)
      call cpl_update (zicedp (1:ni), 'I8I' , ijdrv_phy(1:2,1:ni,trnch:trnch), ni)
      call cpl_update (zsnodp(1:ni,indx_ice:indx_ice), 'SDI', &
           ijdrv_phy(1:2,1:ni,trnch:trnch), ni)
   endif

   ! mg, glsea, glacier, urban, lakefr, et riverfr doivent etre bornes entre 0 et 1
   ! pour que les poids soient valides

   do i=1,ni
      lemin = min(zmg(i),zglsea(i),zglacier(i),zurban(i),zlakefr(i),zriverfr(i))
      lemax = max(zmg(i),zglsea(i),zglacier(i),zurban(i),zlakefr(i),zriverfr(i))
      if (lemin.lt.0. .or. lemax.gt.1.) then
         call msg_toall(MSG_ERROR, '(sfc_main) INVALID WEIGHTS FOR SURFACE PROCESSES, MAKE SURE THAT LAND-SEA MASK, SEA ICE FRACTION, FRACTION OF GLACIERS, MASK OF URBAN AREAS lakes rivers, ARE BOUNDED BETWEEN 0 AND 1')
         return
      endif
   end do

   ! initialisations

   rangs=0. ; poids=0. ; rg_soil=0. ; rg_glacier=0.
   rg_water=0. ; rg_ice=0. ; rg_urb=0. ; rg_lake=0. ; rg_river=0.

   ! calcul des poids

   do i=1,ni

      mask = zmg(i)
      if      (mask.gt.(1.-critmask)) then
         mask = 1.0
      else if (mask.lt.    critmask ) then
         mask = 0.0
      endif

      ! sol
      poids(i,indx_soil)    =     mask  * (1.-zglacier(i))*(1.-zurban(i))
      ! villes
      poids(i,indx_urb )    =     mask  * (1.-zglacier(i))*zurban(i)
      ! glaciers continentaux
      poids(i,indx_glacier) =     mask  * zglacier(i)
      ! eau (***salee/ocean si lacs et rivieres specifiees)
      poids(i,indx_water)   = (1.-mask - zlakefr(i) - zriverfr(i) ) * (1.-zglsea(i))
      !mackay test - restore limit
      ! M.A. to be verified 
       if ( schmlake.eq.'CSLM'.and. poids(i,indx_water) .lt. critmask) then ! If less than 0.1%
           poids(i,indx_water) = 0.0
       endif
      ! glace marine
      poids(i,indx_ice)     = (1.-mask - zlakefr(i) - zriverfr(i) ) *     zglsea(i)
      ! lacs
      poids(i,indx_lake)    =   zlakefr(i)
      if ( schmlake.eq.'CSLM'.and. poids(i,indx_lake) .lt. critmask) then ! If less than 0.1%
          poids(i,indx_lake) = 0.0
      endif
      ! rivieres
      poids(i,indx_river)   =   zriverfr(i)

   end do

   ni_water     = 0; bus_water(1)   = 0.
   ni_ice       = 0; bus_ice(1)     = 0.
   ni_soil      = 0; bus_soil(1)    = 0.
   ni_glacier   = 0; bus_glacier(1) = 0.
   ni_urb       = 0; bus_urb(1)     = 0.
   ni_lake      = 0; bus_lake(1)    = 0.
   ni_river     = 0; bus_river(1)   = 0.

   ! definition des "rangs"
   ! agregation

   DO_AGREG: do i=1,ni

      if (poids(i,indx_soil).gt.0.0) then
         ni_soil                = ni_soil + 1
         rangs(i,indx_soil)     = ni_soil
         rg_soil(ni_soil)       = i
      endif

      if (poids(i,indx_glacier).gt.0.0) then
         ni_glacier             = ni_glacier + 1
         rangs(i,indx_glacier)  = ni_glacier
         rg_glacier(ni_glacier) = i
      endif

      if (poids(i,indx_water).gt.0.0) then
         ni_water               = ni_water + 1
         rangs(i,indx_water)    = ni_water
         rg_water(ni_water)     = i
      endif

      if (poids(i,indx_ice).gt.0.0) then
         ni_ice               = ni_ice + 1
         rangs(i,indx_ice)    = ni_ice
         rg_ice(ni_ice)       = i
      endif

      if (poids(i,indx_urb).gt.0.0) then
         ni_urb               = ni_urb + 1
         rangs(i,indx_urb)    = ni_urb
         rg_urb(ni_urb)       = i
      endif

      if (poids(i,indx_lake).gt.0.0) then
         ni_lake              = ni_lake + 1
         rangs(i,indx_lake)   = ni_lake
         rg_lake(ni_lake)     = i
      endif

      if (poids(i,indx_river).gt.0.0) then
         ni_river             = ni_river + 1
         rangs(i,indx_river)   = ni_river
         rg_river(ni_river)     = i
      endif

   end do DO_AGREG

   !******************************
   !        SOIL                 *
   !******************************

   sommet = 1

   do i=1,nvarsurf
      ptr_soil(i) = sommet
      sommet = sommet + vl(i)%niveaux * vl(i)%mul * vl(i)%mosaik * ni_soil
   end do
   siz_soil = max(surfesptot*ni_soil,1)

   if (siz_soil.gt.1) then
      do i=1,siz_soil
         bus_soil(i) = 0.0
      end do
      ! transvidage des 3 bus dans bus_soil
      call copybus3(bus_soil, siz_soil, &
           ptr_soil, nvarsurf, &
           rg_soil, ni_soil, &
           ni, indx_soil, trnch, CP_TO_SFCBUS)

      if (schmsol.eq.'ISBA') then

         call isba4(bus_soil, siz_soil, &
              ptr_soil, nvarsurf, &
              dt, kount, &
              ni_soil, ni_soil, nk-1)

      elseif (schmsol.eq.'SVS') then
         
         call svs(bus_soil, siz_soil, &
              ptr_soil, nvarsurf, &
              dt, kount, trnch, &
              ni_soil, ni_soil, nk-1)

      endif
      if (phy_error_L) return

      call copybus3(bus_soil, siz_soil, &
           ptr_soil, nvarsurf, &
           rg_soil, ni_soil, &
           ni, indx_soil, trnch, CP_FROM_SFCBUS)
   endif

   !******************************
   !        WATER                *
   !******************************

   sommet = 1
   do i=1,nvarsurf
      ptr_water(i) = sommet
      sommet = sommet + vl(i)%niveaux * vl(i)%mul * vl(i)%mosaik * ni_water
   end do

   siz_water = max(surfesptot*ni_water,1)

   if (siz_water.gt.1) then

      do i=1,siz_water
         bus_water(i) = 0.0
      end do

      ! transvidage des 3 bus dans bus_water
      call copybus3(bus_water, siz_water, &
           ptr_water, nvarsurf, &
           rg_water, ni_water, &
           ni, indx_water, trnch, CP_TO_SFCBUS)

      lcl_indx(1,1:ni_water) = rg_water(1:ni_water)
      lcl_indx(2,1:ni_water) = trnch

      call water2(bus_water, siz_water   ,  &
           ptr_water, nvarsurf    ,  &
           lcl_indx , kount,  &
           ni_water , ni_water, nk-1)
      if (phy_error_L) return

      call copybus3(bus_water, siz_water, &
           ptr_water, nvarsurf, &
           rg_water, ni_water, &
           ni, indx_water, trnch, CP_FROM_SFCBUS)
   endif

   !******************************
   !     OCEANIC ICE             *
   !******************************

   sommet = 1
   do i=1,nvarsurf
      ptr_ice(i) = sommet
      sommet = sommet + vl(i)%niveaux * vl(i)%mul * vl(i)%mosaik * ni_ice
   end do

   siz_ice = max(surfesptot*ni_ice,1)
   do_ice      = .false.

   if (siz_ice.gt.1) then

      do_ice      = .true.

      do i=1,siz_ice
         bus_ice(i) = 0.0
      end do
      ! transvidage des 3 bus dans bus_ice
      call copybus3(bus_ice, siz_ice, &
           ptr_ice, nvarsurf, &
           rg_ice, ni_ice, &
           ni, indx_ice, trnch, CP_TO_SFCBUS)

      lcl_indx(1,1:ni_ice) = rg_ice(1:ni_ice)
      lcl_indx(2,1:ni_ice) = trnch

      call seaice3(bus_ice , siz_ice     , &
           ptr_ice , nvarsurf    , &
           lcl_indx, &
           ni_ice  , ni_ice, nk-1)
      if (phy_error_L) return

      call copybus3(bus_ice, siz_ice, &
           ptr_ice, nvarsurf, &
           rg_ice, ni_ice, &
           ni, indx_ice, trnch, CP_FROM_SFCBUS)
   endif

   !******************************
   !     CONTINENTAL ICE         *
   !******************************

   sommet = 1
   do i=1,nvarsurf
      ptr_glacier(i) = sommet
      sommet = sommet + vl(i)%niveaux * vl(i)%mul * vl(i)%mosaik * ni_glacier
   end do

   siz_glacier = max(surfesptot*ni_glacier,1)

   do_glaciers = .false.

   if (siz_glacier.gt.1) then

      do_glaciers = .true.

      do i=1,siz_glacier
         bus_glacier(i) = 0.0
      end do
      ! transvidage des 3 bus dans bus_glacier
      call copybus3(bus_glacier, siz_glacier, &
           ptr_glacier, nvarsurf, &
           rg_glacier, ni_glacier, &
           ni, indx_glacier, trnch, CP_TO_SFCBUS)

      call glaciers2(bus_glacier, siz_glacier, &
           ptr_glacier, nvarsurf, &
           ni_glacier, ni_glacier, nk-1)
      if (phy_error_L) return

      call copybus3(bus_glacier, siz_glacier, &
           ptr_glacier, nvarsurf, &
           rg_glacier, ni_glacier, &
           ni, indx_glacier, trnch, CP_FROM_SFCBUS)
   endif

   !******************************
   !     URBAN AREAS             *
   !******************************

   sommet = 1
   do i=1,nvarsurf
      ptr_urb(i) = sommet
      sommet = sommet + vl(i)%niveaux * vl(i)%mul * vl(i)%mosaik * ni_urb
   end do

   siz_urb = max(surfesptot*ni_urb,1)
   do_urb = .false.

   if (siz_urb.gt.1) then

      do_urb = .true.

      do i=1,siz_urb
         bus_urb(i) = 0.0
      end do
      ! transvidage des 3 bus dans bus_urb
      call copybus3(bus_urb, siz_urb, &
           ptr_urb, nvarsurf, &
           rg_urb, ni_urb, &
           ni, indx_urb, trnch, CP_TO_SFCBUS)

      call town2(bus_urb, siz_urb, &
           ptr_urb, nvarsurf, &
           dt, kount, &
           ni_urb, ni_urb, &
           nk-1)
      if (phy_error_L) return

      call copybus3(bus_urb, siz_urb, &
           ptr_urb, nvarsurf, &
           rg_urb, ni_urb, &
           ni, indx_urb, trnch, CP_FROM_SFCBUS)
   endif


   !******************************
   !        LAKES                *
   !******************************

   sommet = 1
   do i=1,nvarsurf
      ptr_lake(i) = sommet
      sommet = sommet + vl(i)%niveaux * vl(i)%mul * vl(i)%mosaik * ni_lake
   end do

   siz_lake = max(surfesptot*ni_lake,1)
   write(msg_S,*) siz_lake

   do_lake = .false.

   if (siz_lake.gt.1) then

      do_lake = .true.

      do i=1,siz_lake
         bus_lake(i) = 0.0
      end do

      ! transvidage des 3 bus dans bus_lake
      call copybus3(bus_lake, siz_lake, &
           ptr_lake, nvarsurf, &
           rg_lake, ni_lake, &
           ni, indx_lake, trnch, CP_TO_SFCBUS)

      lcl_indx(1,1:ni_lake) = rg_lake(1:ni_lake)
      lcl_indx(2,1:ni_lake) = trnch

      if (schmlake.eq.'CSLM') then
         
         call cslm_main( bus_lake, siz_lake   ,  &
           ptr_lake, nvarsurf    ,  &
           lcl_indx , trnch, kount,  &
           ni_lake , ni_lake, nk-1 )
      else
         
         call water2( bus_lake, siz_lake   ,  &
           ptr_lake, nvarsurf    ,  &
           lcl_indx , kount,  &
           ni_lake , ni_lake, nk-1 )
         if (phy_error_L) return
      endif

      call copybus3(bus_lake, siz_lake, &
           ptr_lake, nvarsurf, &
           rg_lake, ni_lake, &
           ni, indx_lake, trnch, CP_FROM_SFCBUS)

   endif

   ! !******************************
   ! !        RIVERS       - NOT DEFINED - THIS IS ONLY A PLACEHOLDER
   ! !******************************

   ! sommet = 1
   ! do i=1,nvarsurf
   !    ptr_river(i) = sommet
   !    sommet = sommet + vl(i)%niveaux * vl(i)%mul * vl(i)%mosaik * ni_river
   ! end do

   ! siz_river = max(surfesptot*ni_river,1)

   ! do_river = .false.

   ! if (siz_river.gt.1) then

   !    do_river = .true.

   !    do i=1,siz_river
   !       bus_river(i) = 0.0
   !    end do

   !    ! transvidage des 3 bus dans bus_river
   !    call copybus3(bus_river, siz_river, &
   !         ptr_river, nvarsurf, &
   !         rg_river, ni_river, &
   !         ni, indx_river, trnch, CP_TO_SFCBUS)

   !    lcl_indx(1,1:ni_river) = rg_river(1:ni_river)
   !    lcl_indx(2,1:ni_river) = trnch

   !    call water2( bus_river, siz_river   ,  &
   !         ptr_river, nvarsurf    ,  &
   !         lcl_indx , kount,  &
   !         ni_river , ni_river, nk-1 )
   !       if (phy_error_L) return

   !    call copybus3(bus_river, siz_river, &
   !         ptr_river, nvarsurf, &
   !         rg_river, ni_river, &
   !         ni, indx_river, trnch, CP_FROM_SFCBUS)
   ! endif

!  For CSLM: reduce runofftot and latflaf by the fraction of lakes for all land surfaces 
!  except lake before the aggregation is performed (to be modified once rivers are activated)
   !  and pass the reduced amount to CSLM as input at the next timestep


   !!!!!!!!!! CSLM_RUNOFF  TO BE MOVED AND MADE TO WORK FOR ISBA !!!!!!!!!!!!!!!!
   
   CSLM_RUNOFF:if (schmlake .eq.'CSLM') then
      zrofinlak = 0.

      if (schmsol .eq.'SVS') then
         ! Adjustments only defined for SVS for now. 
         do i=1,ni
            if (poids(i,indx_lake).gt.0.0) then
               ! Surface runoff
               zrofinlak(i)             = zrunofftot(i,indx_soil) * (zlakefr(i)) * poids(i,indx_soil)/poids(i,indx_lake)
               zrunofftot(i,indx_soil)  = zrunofftot(i,indx_soil) * (1-zlakefr(i))
               ! Add a portion of lateral flow coming from the land tile
               do k=1,khyd
                  zrofinlak(i)          = zrofinlak(i) + (zlatflw(i,k) * (zlakefr(i)) * poids(i,indx_soil)/poids(i,indx_lake))
                  zlatflw(i,k)          = zlatflw(i,k)               * (1-zlakefr(i))
               enddo
               ! Add a portion of base drainage from land tile
               zrofinlak(i)             = zrofinlak(i) + (zwatflow(i,khyd+1) * (zlakefr(i)) * poids(i,indx_soil)/poids(i,indx_lake))
               zwatflow(i,khyd+1)       = zwatflow(i,khyd+1)         * (1-zlakefr(i))
               if (do_urb) then
                  zrofinlak(i)          = zrofinlak(i) + (zrunofftot(i,indx_urb) * (zlakefr(i)) * poids(i,indx_urb)/poids(i,indx_lake))
                  zrunofftot(i,indx_urb)= zrunofftot(i,indx_urb)     * (1-zlakefr(i))
               endif
               if (do_glaciers) then
                  zrofinlak(i)          = zrofinlak(i) + (zrunofftot(i,indx_glacier) * (zlakefr(i)) * poids(i,indx_glacier)/poids(i,indx_lake))
                  zrunofftot(i,indx_glacier) = zrunofftot(i,indx_glacier) * (1-zlakefr(i))
               endif
               
            endif
         enddo
      endif
   ENDIF CSLM_RUNOFF

   !******************************
   !     AGREGATION              *
   !******************************

   call agrege3( &
        bus_soil, bus_glacier, bus_water, bus_ice, bus_urb, bus_lake, bus_river, &
        siz_soil, siz_glacier, siz_water, siz_ice, siz_urb, siz_lake, siz_river, &
         ni_soil,  ni_glacier,  ni_water,  ni_ice,  ni_urb,  ni_lake,  ni_river, &
        ptr_soil, ptr_glacier, ptr_water, ptr_ice, ptr_urb, ptr_lake, ptr_river, &
        nvarsurf, &
        rangs, poids, ni, trnch, &
        do_glaciers, do_ice, do_urb, do_lake, do_river)

  
   !*****************************
   ! UPDATE DIAG LEVEL AT START *
   !*****************************

   if (kount == 0) then
      if (any('pw_tt:p'==phyinread_list_s(1:phyinread_n))) ztdiag = ztplus(:,nk)
      if (any('tr/hu:p'==phyinread_list_s(1:phyinread_n))) zqdiag = zhuplus(:,nk)
      if (any('pw_uu:p'==phyinread_list_s(1:phyinread_n))) zudiag = zuplus(:,nk)
      if (any('pw_vv:p'==phyinread_list_s(1:phyinread_n))) zvdiag = zvplus(:,nk)
   endif

   !******************************************
   ! COMPUTE "STATION" DIAGNOSTICS OVER LAND *
   !******************************************

   where (poids(:,indx_soil) > 0.1)
      zqdiagstn(:) = zqdiagtyp(:,indx_soil)
      ztdiagstn(:) = ztdiagtyp(:,indx_soil)
      zudiagstn(:) = zudiagtyp(:,indx_soil)
      zvdiagstn(:) = zvdiagtyp(:,indx_soil)
      zqdiagstnv(:) = zqdiagtypv(:,indx_soil)
      ztdiagstnv(:) = ztdiagtypv(:,indx_soil)
      zudiagstnv(:) = zudiagtypv(:,indx_soil)
      zvdiagstnv(:) = zvdiagtypv(:,indx_soil)
   elsewhere
      zqdiagstn(:) = zqdiagtyp(:,indx_agrege)
      ztdiagstn(:) = ztdiagtyp(:,indx_agrege)
      zudiagstn(:) = zudiagtyp(:,indx_agrege)
      zvdiagstn(:) = zvdiagtyp(:,indx_agrege)
      zqdiagstnv(:) = zqdiagtypv(:,indx_agrege)
      ztdiagstnv(:) = ztdiagtypv(:,indx_agrege)
      zudiagstnv(:) = zudiagtypv(:,indx_agrege)
      zvdiagstnv(:) = zvdiagtypv(:,indx_agrege)
   endwhere

   do k=1,nsurf+1
     call mhuaes3(ztddiagtyp(:,k), zqdiagtyp(:,k), ztdiagtyp(:,k), zpplus, .false., ni, 1, ni)
       ztddiagtyp(:,k)=max(ztddiagtyp(:,k),0.)
       ztddiagtyp(:,k)=ztdiagtyp(:,k)-ztddiagtyp(:,k)
     call mhuaes3(ztddiagtypv(:,k), zqdiagtypv(:,k), ztdiagtypv(:,k), zpplus, .false.,ni, 1, ni)
       ztddiagtypv(:,k)=max(ztddiagtypv(:,k),0.)
       ztddiagtypv(:,k)=ztdiagtypv(:,k)-ztddiagtypv(:,k)
   enddo
   
   !*******************************************
   !     CORRECTIFS POUR LA SORTIE            *
   !*******************************************
   ! Copy for output only
   poids_out(:,indx_soil   ) = poids(:,indx_soil   )
   poids_out(:,indx_glacier) = poids(:,indx_glacier)
   poids_out(:,indx_water  ) = poids(:,indx_water  )
   poids_out(:,indx_ice    ) = poids(:,indx_ice    )

   if (SCHMURB.ne.'NIL') then
      poids_out(:,indx_urb    ) = poids(:,indx_urb    )
   endif
   if (SCHMLAKE.ne.'NIL') then
      poids_out(:,indx_lake   ) = poids(:,indx_lake   )
   endif
   if (SCHMRIVER.ne.'NIL') then
      poids_out(:,indx_river  ) = poids(:,indx_river  )
   endif

   do i=1,ni
      poids_out(:,indx_agrege) = 1.0
   end do

   do k=1,nsurf+1
      do i=1,ni
         zilmo(i,k)=max(-10.,min(10.,zilmo(i,k)))
      end do
   enddo

      !     Calcul du flux d'evaporation en kg/m2 a partir
      !     du flux d'humidite en m/s*kg/kg: m/s * kg water / kg air:
      !     il faut multiplier par la densite de l'air
      !     and by the length of the timestep
      !     It only makes sense theoretically to calculate fvapliq from the aggregated field
      !     since HU is averaged over all surface types
   do i = 1,ni
      rhoair = zpmoins(i) / ( RGASD * ztdiag(i) * (1. + DELTA * zqdiag(i)) )
      zfvapliq(i) = zfvap(i,indx_agrege) * dt * rhoair
      zfvapliqaf(i) = zfvapliqaf(i) + zfvapliq(i)
   enddo

   !******************************************
   !     METTRE UNE LIMITE SUR L'AMPLITUDE   *
   !     DE L'INVERSION DE TEMPERATURE       *
   !******************************************

   ztnolim = ztdiag
   if (sl_adjust(ztmoins(:,nk-1),zztsl,ztdiag,hghtt_diag=zt,adj_t_diag=ztdiag) /= SL_OK) &
        write(6,*) 'itf_sfc_main cannot compute adjusted diagnostics'
   zdtdiag = ztdiag - ztnolim

   call msg_toall(MSG_DEBUG, 'sfc_main [END]')
   !     ---------------------------------------------------------------
   F_istat = RMN_OK
   return
end function sfc_main2

end module sfc_main
