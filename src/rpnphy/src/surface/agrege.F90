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

subroutine agrege2( &
     bus_soil, bus_glacier, bus_water, bus_ice, bus_urb, &
     siz_soil, siz_glacier, siz_water, siz_ice, siz_urb, &
     ni_soil,  ni_glacier,  ni_water,   ni_ice,  ni_urb, &
     ptr_soil, ptr_glacier, ptr_water, ptr_ice, ptr_urb, &
     ptsurfsiz, &
     rangs, poids, ni, trnch, &
     do_glaciers, do_ice, do_urb)
   use sfcbus_mod
   implicit none
!!!#include <arch_specific.hf>

   integer ni, trnch, ptsurfsiz
   integer ptr_soil(ptsurfsiz), ptr_glacier(ptsurfsiz)
   integer ptr_water  (ptsurfsiz), ptr_ice  (ptsurfsiz)
   integer ptr_urb  (ptsurfsiz)
   integer siz_soil, siz_water, siz_glacier, siz_ice, siz_urb
   integer ni_water, ni_soil, ni_ice, ni_glacier, ni_urb
   integer rangs(ni,*)
   real bus_water    (siz_water), bus_soil  (siz_soil)
   real bus_ice(siz_ice), bus_glacier(siz_glacier)
   real bus_urb(siz_urb)
   real poids(ni,*)
   logical do_glaciers, do_ice, do_urb

   !@Author B. Bilodeau Sept 1999
   !@Object Perform the aggregation of arrays orginating
   !             from the 5 surface modules (sea ice, glaciers,
   !             water, soil and urban).
   !@Arguments
   !             - Input -
   ! BUS_SOIL    "mini-bus" for soil     surface (from D, F and V)
   ! BUS_GLACIER "mini-bus" for glaciers surface (from D, F and V)
   ! BUS_WATER   "mini-bus" for water    surface (from D, F and V)
   ! BUS_ICE     "mini-bus" for sea ice  surface (from D, F and V)
   ! BUS_URB     "mini-bus" for urban    surface (from D, F and V)
   ! SIZ_SOIL    dimension of BUS_SOIL
   ! SIZ_GLACIER dimension of BUS_GLACIER
   ! SIZ_WATER   dimension of BUS_WATER
   ! SIZ_ICE     dimension of BUS_ICE
   ! SIZ_URB     dimension of BUS_URB
   ! NI_SOIL     length of the row for BUS_SOIL
   ! NI_GLACIER  length of the row for BUS_GLACIER
   ! NI_WATER    length of the row for BUS_WATER
   ! NI_ICE      length of the row for BUS_ICE
   ! NI_URB      length of the row for BUS_URB
   ! NI          length of the full row
   ! PTR_SOIL    starting location of each variable in the soil    "mini-bus"
   ! PTR_GLACIER    "        "     "   "     "      "   "  glacier     "
   ! PTR_WATER      "        "     "   "     "      "   "  water       "
   ! PTR_ICE        "        "     "   "     "      "   "  ice         "
   ! PTR_URB        "        "     "   "     "      "   "  urb         "
   ! PTSURFSIZ   number of elements (variables) in the "mini-buses"
   ! RANGS       index  of each "mini-bus" element in the full row
   ! POIDS       weight of each "mini-bus" element in the tile
   ! DO_GLACIERS .TRUE. some points on the actual row contain glaciers
   !             .FALSE. no  point  "   "    "     "    "       "
   ! DO_ICE      .TRUE. some points on the actual row contain sea ice
   !             .FALSE. no  point  "   "    "     "    "       "
   ! DO_URB      .TRUE. some points on the actual row contain urban areas
   !             .FALSE. no  point  "   "    "     "    "       "    "

   integer ik_soil, ik_glacier, ik_water, ik_ice, ik_ori, ik_urb
   integer i, ik, k, m, var
   integer n, var_no

   real, pointer, dimension(:) :: bus_ori

   !     boucle sur les variables du bus de surface
   DO_VAR: do var=1,nvarsurf

      if (vl(var)%niveaux <= 0 .or. .not.vl(var)%doagg_L) cycle DO_VAR
      if (.not.associated(busptr(var)%ptr)) cycle DO_VAR

      bus_ori(1:) => busptr(var)%ptr(1:,trnch)

      DO_MUL: do m=1,vl(var)%mul     ! tient compte de la multiplicite des champs

         ! seules certaines variables sont agregees
         if (m /= vl(var)%agg) cycle DO_MUL

!VDIR NODEP
         do ik=1,vl(var)%niveaux*ni       ! double boucle sur i et k

            k = (ik-1)/ni + 1
            i = ik - (k-1)*ni

            !                    agregation des fractions de terre et d'eau
            ik_ori = 1 + &
                 (m-1)*ni                     +   &! multiplicite
                 ik - 1                            ! element ik courant

            ik_soil = ptr_soil(var)           + &
                 (m-1)*ni_soil                + &
                 (k-1)*ni_soil                + &
                 rangs(i,indx_soil) - 1

            if (rangs(i,indx_soil).eq.0) ik_soil = 1

            ik_water   = ptr_water(var)       + &
                 (m-1)*ni_water               + &
                 (k-1)*ni_water               + &
                 rangs(i,indx_water) - 1

            if (rangs(i,indx_water).eq.0) ik_water = 1

            bus_ori(ik_ori) = &
                 poids(i,indx_soil ) * bus_soil (ik_soil)   + &
                 poids(i,indx_water) * bus_water(ik_water)
         end do

         !                 si necessaire, agregation des fractions de
         !                 glace marine, de glace continentale et
         !                 de zones urbaines
         if (do_glaciers .or. do_ice .or. do_urb) then
!VDIR NODEP
            do ik=1,vl(var)%niveaux*ni

               k = (ik-1)/ni + 1
               i = ik - (k-1)*ni

               ik_ori = 1 + &
                    (m-1)*ni                             + &
                    ik -1

               ik_glacier= ptr_glacier(var)                  + &
                    (m-1)*ni_glacier                  + &
                    (k-1)*ni_glacier                  + &
                    rangs(i,indx_glacier) - 1

               if (rangs(i,indx_glacier).eq.0) ik_glacier = 1

               ik_ice  = ptr_ice(var)                        + &
                    (m-1)*ni_ice                        + &
                    (k-1)*ni_ice                        + &
                    rangs(i,indx_ice) - 1

               if (rangs(i,indx_ice).eq.0) ik_ice   = 1

               ik_urb  = ptr_urb(var)                        + &
                    (m-1)*ni_urb                        + &
                    (k-1)*ni_urb                        + &
                    rangs(i,indx_urb) - 1

               if (rangs(i,indx_urb).eq.0) ik_urb   = 1

               bus_ori(ik_ori)                               = &
                    bus_ori(ik_ori)                               + &
                    poids(i,indx_glacier)* bus_glacier(ik_glacier)+ &
                    poids(i,indx_ice    )* bus_ice    (ik_ice  )  + &
                    poids(i,indx_urb    )* bus_urb    (ik_urb)

            end do
         endif
      end do DO_MUL
   end do DO_VAR

   !     Cas particuliers : on moyenne "TSRAD"**4 pour
   !     les besoins du schema de rayonnement.
   !     tsrad_i est l'indice de la variable "TSRAD";
   !     il est initialise dans le s/p iniptsurf.
   !     Autres exceptions : on moyenne ln(Z0) et ln(Z0T).

   DO_VAR2: do n=1,3

      if (n.eq.1) var_no = tsrad_i
      if (n.eq.2) var_no = z0_i
      if (n.eq.3) var_no = z0t_i

      bus_ori(1:) => busptr(var_no)%ptr(1:,trnch)

      DO_MUL2: do m=1,vl(var_no)%mul

!VDIR NODEP
         DO_IK2: do ik=1,vl(var_no)%niveaux*ni

            k = (ik-1)/ni + 1
            i = ik - (k-1)*ni

            !              agregation des fractions de terre, d'eau,
            !              de glace marine, de glaciers continentaux
            !              et de zones urbaines
            ik_ori = 1 + &
                 (m-1)*ni                                       + &
                 ik - 1


            ik_soil = ptr_soil(var_no)                              + &
                 (m-1)*ni_soil                                 + &
                 (k-1)*ni_soil                                 + &
                 rangs(i,indx_soil) - 1

            if (rangs(i,indx_soil).eq.0) ik_soil = 1

            ik_water   = ptr_water(var_no)                          + &
                 (m-1)*ni_water                             + &
                 (k-1)*ni_water                             + &
                 rangs(i,indx_water) - 1

            if (rangs(i,indx_water).eq.0) ik_water = 1

            ik_glacier= ptr_glacier(var_no)                         + &
                 (m-1)*ni_glacier                            + &
                 (k-1)*ni_glacier                            + &
                 rangs(i,indx_glacier) - 1

            if (rangs(i,indx_glacier).eq.0) ik_glacier = 1

            ik_ice  = ptr_ice(var_no)                               + &
                 (m-1)*ni_ice                                  + &
                 (k-1)*ni_ice                                  + &
                 rangs(i,indx_ice) - 1

            if (rangs(i,indx_ice).eq.0) ik_ice = 1

            ik_urb  = ptr_urb(var_no)                               + &
                 (m-1)*vl(var_no)%niveaux*ni_urb                  + &
                 (k-1)*ni_urb                                  + &
                 rangs(i,indx_urb) - 1

            if (rangs(i,indx_urb).eq.0) ik_urb = 1



            if (var_no.eq.tsrad_i) then

               !              moyenne de la 4eme puissance de tsrad

               bus_ori(ik_ori) = &
                    poids(i,indx_soil   ) * (bus_soil   (ik_soil   )**4) + &
                    poids(i,indx_water  ) * (bus_water  (ik_water  )**4) + &
                    poids(i,indx_glacier) * (bus_glacier(ik_glacier)**4) + &
                    poids(i,indx_ice    ) * (bus_ice    (ik_ice    )**4) + &
                    poids(i,indx_urb    ) * (bus_urb    (ik_urb    )**4)

               bus_ori(ik_ori) = bus_ori(ik_ori)**0.25

            else if ((var_no.eq.z0_i.or.var_no.eq.z0t_i).and. &
                 (m.eq.vl(var_no)%agg)) then

               !              moyenne logarithmique de "Z0" et de "Z0T"

               bus_ori(ik_ori) = &
                    poids(i,indx_soil   ) * &
                    alog(max(1.e-10,bus_soil   (ik_soil   ))) + &
                    poids(i,indx_water  ) * &
                    alog(max(1.e-10,bus_water  (ik_water  ))) + &
                    poids(i,indx_glacier) * &
                    alog(max(1.e-10,bus_glacier(ik_glacier))) + &
                    poids(i,indx_ice    ) * &
                    alog(max(1.e-10,bus_ice    (ik_ice    ))) + &
                    poids(i,indx_urb    ) * &
                    alog(max(1.e-10,bus_urb    (ik_urb    )))

               bus_ori(ik_ori) = exp(bus_ori(ik_ori))

            endif

         end do DO_IK2
      end do DO_MUL2
   end do DO_VAR2

   return
end subroutine agrege2
