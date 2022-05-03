!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!-------------------------------------- LICENCE END --------------------------------------

subroutine copybus3(bus_sfc, sfcsiz, &
     ptsurf, ptsurfsiz, &
     masque, ni_sfc, ni, &
     indx_sfc, trnch, ramasse)
   use sfcbus_mod
   implicit none
!!!#include <arch_specific.hf>

   logical ramasse
   integer sfcsiz, indx_sfc, ni, ni_sfc, ptsurfsiz
   integer masque(ni_sfc), ptsurf(ptsurfsiz), trnch
   real bus_sfc(sfcsiz)

   !@Author  B. Bilodeau Sept 1999
   !@Object  Performs massive gather-scather for the
   !             surface processes.
   !             Copy the contents of the the 3 main buses
   !             (dynamics, permanent and volatile) in the
   !             "mini-buses" corresponding to each of the
   !             4 surface modules (sea ice, glaciers, water
   !             and soil) before the execution of the surface
   !             processes, if switch "ramasse" is true. Copy
   !             in the opposite direction after the execution
   !             of the surface processes if switch "ramasse"
   !             is false.
   !@Arguments
   !             - Input/Output -
   ! BUS_SFC     "mini-bus" for one of 4 surface types
   !             - Input -
   ! SFCSIZ      dimension of BUS_SFC
   ! PTRSURF     starting location of each variable in the "mini-bus"
   ! PTSURFSIZ   number of elements (variables)     in the "mini-bus"
   ! NI_SFC      length of the row for BUS_SFC
   ! NI          length of the full row
   ! MASQUE      index  of each bus element in the "mini-bus"
   ! POIDS       weight of each "mini-bus" element in the tile
   ! INDX_SFC    integer value (1-4) corresponding to each surface type
   ! RAMASSE     .TRUE.  scatter from main buses   to "mini-buses"
   !             .FALSE. gather  from "mini-buses" to main buses

#include <rmnlib_basics.hf>

   integer i, ik_sfc, ik_ori, k, m, var
   integer surflen
   integer pseudo_m

   real, pointer :: bus_ori(:)

   surflen = ni_sfc

   DO_VAR: do var = 1, nvarsurf

      if (vl(var)%niveaux <= 0 .or. .not.vl(var)%doagg_L) cycle DO_VAR
      if (.not.associated(busptr(var)%ptr)) cycle DO_VAR

      bus_ori(1:) => busptr(var)%ptr(1:,trnch)

      DO_MUL: do m = 1, vl(var)%mul*vl(var)%mosaik

         pseudo_m = m
         if (m > vl(var)%mul) pseudo_m = mod(m-1, vl(var)%mul) + 1

         IF_RAMASSE: if (ramasse) then

            do k = 1, vl(var)%niveaux
               do i=1,ni_sfc
                  ik_sfc = ptsurf(var)        +  &! debut du champ dans bus_sfc
                       (m-1)*ni_sfc +            &! multiplicite et mosaique
                       (k-1)*ni_sfc+i-1           ! element ik courant

                  ik_ori = 1 + &
                       (m-1)*ni               +  &! multiplicite et mosaique
                       (k-1)*ni               +  &! rangees precedentes
                       masque(i) - 1              ! element i courant

                  bus_sfc(ik_sfc) = bus_ori(ik_ori)
               enddo
            enddo

         else if (statut(var,pseudo_m) == indx_sfc) then

            do k = 1, vl(var)%niveaux
               do i=1,ni_sfc
                  ik_sfc = ptsurf(var)        +  &! debut du champ dans bus_sfc
                       (m-1)*ni_sfc           +  &! multiplicite et mosaique
                       (k-1)*ni_sfc+i-1           ! element ik courant

                  ik_ori = 1 + &
                       (m-1)*ni               +  &! multiplicite et mosaique
                       (k-1)*ni               +  &! rangees precedentes
                       masque(i) - 1              ! element i courant

                  bus_ori(ik_ori) = bus_sfc(ik_sfc)
               enddo
            enddo

         endif IF_RAMASSE

      enddo DO_MUL
   enddo DO_VAR

   return
end subroutine copybus3
