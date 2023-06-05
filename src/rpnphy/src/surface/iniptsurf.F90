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
!-------------------------------------- LICENCE END ---------------------------

!/@*
function iniptsurf5() result(F_istat)
   use sfc_options
   use sfcbus_mod
   implicit none
!!!#include <arch_specific.hf>
   !@Object Initialization of common blocks used in the surface package
   !@Arguments

   integer :: F_istat

   !@Author B. Bilodeau (Sept 1999)
   !@Revisions
   ! 001      B. Bilodeau (Nov 2000) - New comdeck sfcbus.cdk
   ! 002      B. Bilodeau (Feb 2004) - Revised logic to facilitate the
   !                                   addition of new types of surface
   ! 004      A. Lemonsu  (Jun 2005) - Add variables for urban module
   ! 005      B. Dugas    (Aug 2005) - Initialize sfcbus character variables
   !                                   in Block DATA SFCBUS_DATA
   ! 006      B. Bilodeau (Jun 2005) - Add mosaic capability for CLASS
   ! 007      M. Desgagne (Apr 2009) - Add coupling values (TMICECPL, MCCPL)
   ! 008      V.Lee (Mar 2011)       - mosaik=real number of mosaic tiles + 1
   !                                   to calculate pointer position
   ! 009      M. Mackay (Sep 2022)   - CSLM added
   !*@/

#include <rmn/msg.h>
#include <rmnlib_basics.hf>

   integer, parameter :: nb_agrege = 48
   integer, parameter :: nb_glaciers = 1
   integer, parameter :: nb_water = 4
   integer, parameter :: nb_ice = 2
   integer, parameter :: nb_urb = 95
   integer, parameter :: nb_lake = 26
   integer, parameter :: nb_river = 4

   character(len=16) :: agrege_out(nb_agrege), &
        glaciers_out(nb_glaciers), water_out(nb_water), ice_out(nb_ice), &
        urb_out(nb_urb), lake_out(nb_lake), river_out(nb_river), tmp_S
   integer :: j, l, m, ier
 
   ! les variables de sortie du module "soils" ont preseance
   ! sur celles de tous les autres schemas, sauf exceptions
   ! contenues dans les listes plus bas

   ! liste des variables de surface a agreger
   !#TODO: check that a trimming is not done on longer vars
   data agrege_out   / &
        !  ces variables sont moyennees lineairement
        'ALFAQ'    , 'ALFAT'    , 'ALVIS'    , 'BM'       , 'BT'       , &
        'EMISR'    , 'FC'       , 'FRV'      , 'FTEMP'    , 'FV'       , &
        'FVAP'     , 'HST'      , 'ILMO'     , &
        'QDIAG'    , 'QSURF'    , 'RUNOFFTOT', 'SNODP'    , 'TDIAG'    , &
        'TSURF'    , 'UDIAG'    , 'VDIAG'    , &
        'QDIAGTYP' , 'TDIAGTYP' , 'UDIAGTYP' , 'VDIAGTYP' , &
        'QDIAGTYPV', 'TDIAGTYPV', 'UDIAGTYPV', 'VDIAGTYPV', &
        'YUTCISUN' , 'YUTCISHADE' ,                                      &
        'YWBGTSUN' , 'YWBGTSHADE' , 'YRADSUN', 'YRADSHADE',              &
        'YTGLBSUN' , 'YTGLBSHADE' , 'YTWETB' ,                           &
        'YQ1' , 'YQ2' , 'YQ3' , 'YQ4' , 'YQ5' ,                          &
        'YQ6' , 'YQ7' ,    &
 
        !  le flux infrarouge emis par la surface, qui est
        !  proportionnel a TSRAD**4, est moyenne lineairement
        'TSRAD'    , &
 
        !  on prend la moyennne logarithmique des longueurs de rugosite
        'Z0'       , 'Z0T' &

        /

   ! liste des variables de sortie du module "glaciers"
   data glaciers_out / &
        'TGLACIER' &
        /

   ! liste des variables de sortie du module "water"
   data water_out    / &
        'TWATER', 'DSST', 'SKIN_DEPTH', 'SKIN_INC' &
        /

   ! liste des variables de sortie du module "ice"
   data ice_out      / &
        'ICEDP'    , 'TMICE' /

   ! liste des variables de sortie du module "urb"
   data urb_out      / &
        'T_CANYON' , 'Q_CANYON' , 'U_CANYON' , 'TI_BLD'   , 'T_ROOF'   , &
        'T_ROAD'   , 'T_WALL'   , 'RN_TOWN'  , 'H_TOWN'   , 'LE_TOWN'  , &
        'G_TOWN'   , 'RN_ROOF'  , 'H_ROOF'   , 'LE_ROOF'  , 'G_ROOF'   , &
        'RN_ROAD'  , 'H_ROAD'   , 'LE_ROAD'  , 'G_ROAD'   , 'RN_WALL'  , &
        'H_WALL'   , 'LE_WALL'  , 'G_WALL'   , 'TI_ROAD'  , 'WS_ROOF'  , &
        'WS_ROAD'  , 'QDIAG_CAN', 'UDIAG_CAN', 'VDIAG_CAN', 'TDIAG_CAN', &
        'SROOF_WSNOW', 'SROOF_T'  , 'SROOF_RHO', 'SROOF_ALB', 'SROOF_EMIS',&
        'SROOF_TS' ,'SROAD_WSNOW','SROAD_T'  , 'SROAD_RHO', 'SROAD_ALB', &
        'SROAD_EMIS', 'SROAD_TS' , 'Z0_TOWN'  , 'Z0_ROOF'  , 'Z0_ROAD'  , &
        'BLD_HEIGHT','WALL_O_HOR','CAN_HW_RATIO','ALB_ROOF', 'ALB_ROAD' , &
        'ALB_WALL' , 'EMIS_ROOF', 'EMIS_ROAD', 'EMIS_WALL', 'HC_ROOF'  , &
        'HC_ROAD' , 'HC_WALL'  , 'TC_ROOF'  , 'TC_ROAD'  , 'TC_WALL'  , &
        'D_ROAD'  , 'D_ROOF'   , 'D_WALL'   , 'H_TRAFFIC', 'H_INDUSTRY',&
        'LE_TRAFFIC','LE_INDUSTRY', 'EMTW', 'ALSCATW', 'TSRADTW'        &
        ,'YRADIN', 'YRADRFSUN', 'YRADRFSHADE','YUTCIIN','YUTCIRFSUN'    &
        ,'YUTCIRFSHADE', 'YUTCICIN','YWBGTRFSUN','YWBGTRFSHADE'         &
        ,'YUTCICSUN','YUTCICSHADE', 'YUTCICRFSUN','YUTCICRFSHADE'       &
        ,'YTGLBRFSUN','YTGLBRFSHADE','YTWETBRF', 'YTRFZT', 'YTRDZT',    &
        'YURDZU' , 'YQ8', 'YQ9', 'YQ10','YQ11' ,'YQ12' , 'YQ13'         &
        /

   ! liste des variables de sortie du modules "lake"
   data lake_out    / &
        'TLAK', 'LST', 'TKE',  'HDPTH', 'LKICEH', 'SNICEH', 'EXPW',  &
        'DTEMP', 'DELU', 'GRED', 'RHOMIX', 'TSED', 'ROFICEH',           &
        'SNOL', 'RHOSNOL', 'TSNOWL', 'ALBSNOL', 'WSNOWL', 'HLAKSIL',    &
        'ROFINLAK', 'FICL', 'LFXI', 'LFXO', 'LSTD', 'LSTF', 'EVLAK'  &
        /
  ! liste des variables de sortie du modules "river"
   data river_out    / &
        'LST', 'DSST', 'SKIN_DEPTH', 'SKIN_INC' &
        /

   F_istat = RMN_ERR

   ier = sfcbus_init()
   if (.not.RMN_IS_OK(ier)) return

   ! multiple 2D fields and mosaic tiles are stored in
   ! slices of NI where the order is as follows:
   ! Example is for Mult=3, mosaic = 3
   ! 1 (multiple 1)
   ! 2 (multiple 2)
   ! 3 (multiple 3)
   ! 1.01 (aggregate of tile 1 for multiple 1
   ! 2.01 (aggregate of tile 1 for multiple 2
   ! 3.01 (aggregate of tile 1 for multiple 3
   ! 1.02 (aggregate of tile 2 for multiple 1
   ! 2.02 (aggregate of tile 2 for multiple 2
   ! 3.02 (aggregate of tile 2 for multiple 3
   ! 1.03 (aggregate of tile 2 for multiple 1
   ! 2.03 (aggregate of tile 2 for multiple 2
   ! 3.03 (aggregate of tile 2 for multiple 3

   ! exploration du bus dynamique

   surfesptot = 0
   DO_JVAR2: do j = 1, nvarsurf

!!$      statut(:, 1:vl(j)%mul) = indx_soil

      surfesptot = surfesptot + vl(j)%mul * vl(j)%mosaik * vl(j)%niveaux


      !  initialisation de la variable "statut",
      !  pour le controle des variables qui seront
      !  soit agregees (moyennees), soit sorties
      !  directement d'un bus des bus de surface
      !  correspondant a chacun des 5 types de surface :
      !  statut = 1 --> bus de "sol"      vers bus permanent ou volatil
      !         = 2 --> bus de "glaciers"  "    "      "      "     "
      !         = 3 --> bus de "water"     "    "      "      "     "
      !         = 4 --> bus de "ice"       "    "      "      "     "
      !         = 5 --> moyenne des 5      "    "      "      "     "
      !         = 6 --> bus de "urb"       "    "      "      "     "
      !         = 7 --> bus de "lakes"     "    "      "      "     "
      !         = 8 --> bus de "river"     "    "      "      "     "
      !  voir comdeck "indx_sfc.cdk"

      !  variables agregees
      DO_AGG: do l=1,nb_agrege
         IF_AGG_OUT: if (agrege_out(l) ==  vl(j)%n) then

            if (vl(j)%mul == 1) then
               ! cas no 1 : variables agregees de dimension 1

               vl(j)%agg = 1
               statut(j,1) = indx_agrege

            else if (vl(j)%mul == nsurf+1) then
               ! cas no 2 : variables agregees pour lesquelles on conserve
               ! non seulement la moyenne mais aussi les valeurs associees
               ! a chaque type de surface

               vl(j)%agg              = indx_agrege
               statut(j,indx_agrege)  = indx_agrege
               statut(j,indx_soil  )  = indx_soil
               statut(j,indx_glacier) = indx_glacier
               statut(j,indx_water  ) = indx_water
               statut(j,indx_ice    ) = indx_ice
               if (schmurb /= 'NIL') then
                  statut(j, indx_urb) = indx_urb
               endif
               if (schmlake /= 'NIL') then
                  statut(j, indx_lake) = indx_lake
               endif

               if (schmriver /= 'NIL') then
                  statut(j, indx_river) = indx_river
               endif

            else if (vl(j)%mul > 1 .and. vl(j)%mul /= nsurf+1) then

               print *,'(iniptsurf3) ',vl(j)%mul,nsurf,trim(vl(j)%n)
               call msg(MSG_ERROR, '(iniptsurf3) MULTIPLICITY FACTOR EXCEEDED FOR VAR: '//trim(vl(j)%n))
               return

            endif

         endif IF_AGG_OUT
      end do DO_AGG

      ! variables de sortie du module "glaciers"
      ! Tous les "niveaux" de la variable sont assignes au module glaciers
      do l = 1, nb_glaciers
         if (glaciers_out(l) == vl(j)%n) &
              statut(j, 1:vl(j)%mul) = indx_glacier
      end do

      ! variables de sortie du module "water"
      ! Tous les "niveaux" de la variable sont assignes au module water
      do l = 1, nb_water
         if (water_out(l) == vl(j)%n) &
              statut(j, 1:vl(j)%mul) = indx_water
      end do

      ! variables de sortie du module "ice"
      ! Tous les "niveaux" de la variable sont assignes au module ice
      do l = 1, nb_ice
         if (ice_out(l) == vl(j)%n) &
              statut(j, 1:vl(j)%mul) = indx_ice
      end do

      if (schmurb /= 'NIL') then
         ! variables de sortie du module "urb"
         ! Tous les "niveaux" de la variable sont assignes au module urb
         do l = 1, nb_urb
            if (urb_out(l) == vl(j)%n) &
               statut(j, 1:vl(j)%mul) = indx_urb
         end do
      endif

      if (schmlake /= 'NIL') then
         ! variables de sortie du module "lake"
         ! Tous les "niveaux" de la variable sont assignes au module lake
         do l = 1, nb_lake
            if (lake_out(l) == vl(j)%n) &
               statut(j, 1:vl(j)%mul) = indx_lake
         end do
      endif

      if (schmriver /= 'NIL') then
         ! variables de sortie du module "river"
         ! Tous les "niveaux" de la variable sont assignes au module river
         do l = 1, nb_river
            if (river_out(l) == vl(j)%n) &
               statut(j, 1:vl(j)%mul) = indx_river
         end do
      endif

      ! les autres variables seront transferees du module "soils"
      do m=1,vl(j)%mul
         if (statut(j,m) == 0) statut(j,m) = indx_soil
      end do

   end do DO_JVAR2

   call msg(MSG_INFO,'(iniptsurf) TYPES OF SURFACE :')
   write(tmp_S, '(i0)') indx_soil
   call msg(MSG_INFO,'(iniptsurf) SOIL             '//trim(tmp_S))
   write(tmp_S, '(i0)') indx_glacier
   call msg(MSG_INFO,'(iniptsurf) GLACIERS         '//trim(tmp_S))
   write(tmp_S, '(i0)') indx_water
   call msg(MSG_INFO,'(iniptsurf) WATER            '//trim(tmp_S))
   write(tmp_S, '(i0)') indx_ice
   call msg(MSG_INFO,'(iniptsurf) MARINE ICE       '//trim(tmp_S))
   if (schmurb /= 'NIL') then
      write(tmp_S, '(i0)') indx_urb
      call msg(MSG_INFO,'(iniptsurf) URBAN AREAS      '//trim(tmp_S))
   endif
   if (schmlake /= 'NIL') then
      write(tmp_S, '(i0)') indx_lake
      call msg(MSG_INFO,'(iniptsurf) LAKES            '//trim(tmp_S))
   endif
   if (schmriver /= 'NIL') then
      write(tmp_S, '(i0)') indx_river
      call msg(MSG_INFO,'(iniptsurf) RIVERS           '//trim(tmp_S))
   endif
   write(tmp_S, '(i0)') indx_agrege
   call msg(MSG_INFO,'(iniptsurf) AGGREGATED VALUE '//trim(tmp_S))

   F_istat = RMN_OK
   return
end function iniptsurf5
