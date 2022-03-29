!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html

!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

subroutine ISCCP_SIM1(tau, ptop,                      &! OUTPUT
     il1, il2, npoints, nlev, top_height,             &! INPUT
     pfull, phalf, qv, frac_out, dtau_in, dem_in, at, &
     skt, emsfc_lw)
   implicit none
!!!#include <arch_specific.hf>

   !@Author Jason Cole Dec. 16, 2005
   !@Object PART OF THE ISCCP CLOUD SIMULATOR PACKAGE

   ! There seems to three distinct parts in the original ISCCP simulator
   ! 1. Generation of subgrid columns
   ! 2. Computation of cloud top pressure and optical thickness for subcolumns
   ! 3. GCM grid results computed using results from 2.
   ! I have rewritten the code so that these functions are in seperate subroutines.
   ! In the process I have replace the cloud generator with that created by Raisasen
   ! and I have added the diagnostic of cloud inhomogeneity following some code
   ! from Robert Pincus.
   ! I have also modified the code so that it works on only 1 subcolumn at a time to
   !  save memory
   ! Jason Cole Dec. 16, 2005

   ! Copyright Steve Klein and Mark Webb 2002 - all rights reserved.

   ! This code is available without charge with the following conditions:

   !  1. The code is available for scientific purposes and is not for
   !     commercial use.
   !  2. Any improvements you make to the code should be made available
   !     to the to the authors for incorporation into a future release.
   !  3. The code should not be used in any way that brings the authors
   !     or their employers into disrepute.

   !     NOTE:   the maximum number of levels and columns is set by
   !             the following parameter statement

   !     -----
   !     Input
   !     -----
   integer il1, il2                  !  start and end point in horizontal
   integer npoints                   !  number of model points in the horizontal
   integer nlev                      !  number of model levels in column

   real pfull(npoints,nlev)                      !  pressure of full model levels (Pascals)
   !  pfull(npoints,1)    is    top level of model
   !  pfull(npoints,nlev) is bottom level of model

   real phalf(npoints,nlev+1)        !  pressure of half model levels (Pascals)
   !  phalf(npoints,1)    is    top       of model
   !  phalf(npoints,nlev+1) is the surface pressure

   real qv(npoints,nlev)             !  water vapor specific humidity (kg vapor/ kg air)
   !         on full model levels

   real dtau_in(npoints,nlev)        !  mean 0.67 micron optical depth of stratiform
   !  clouds in each model level
   !  NOTE:  this the cloud optical depth of only the
   !         cloudy part of the grid box, it is not weighted
   !         with the 0 cloud optical depth of the clear
   !         part of the grid box

   real frac_out(npoints,nlev)       ! Cloud mask 0 = no cloud, 1 = cloud


   integer top_height                !  1 = adjust top height using both a computed
   !  infrared brightness temperature and the visible
   !  optical depth to adjust cloud top pressure. Note
   !  that this calculation is most appropriate to compare
   !  to ISCCP data during sunlit hours.
   !  2 = do not adjust top height, that is cloud top
   !  pressure is the actual cloud top pressure
   !  in the model
   !  3 = adjust top height using only the computed
   !  infrared brightness temperature. Note that this
   !  calculation is most appropriate to compare to ISCCP
   !  IR only algortihm (i.e. you can compare to nighttime
   !  ISCCP data with this option)


   !     The following input variables are used only if top_height = 1 or top_height = 3

   real skt(npoints)                 !  skin Temperature (K)
   real emsfc_lw(npoints)            !  10.5 micron emissivity of surface (fraction)
   real at(npoints,nlev)             !  temperature in each model level (K)
   real dem_in(npoints,nlev)         !  10.5 micron longwave emissivity of stratiform
   !  clouds in each
   !  model level.  Same note applies as in dtau_in.
   !     ------
   !     Output
   !     ------

   real tau(npoints)              !  optical thickness in each column

   real ptop(npoints)             !  cloud top pressure (mb) in each column


   !     ------
   !     Working variables added when program updated to mimic Mark Webb's PV-Wave code
   !     ------


   real dem(npoints),bb(npoints)     !  working variables for 10.5 micron longwave
   !  emissivity in part of
   !  gridbox under consideration

   real ptrop(npoints)
   real attrop(npoints)
   real attropmin (npoints)
   real atmax(npoints)
   real atmin(npoints)
   real btcmin(npoints)
   real transmax(npoints)

   integer j,ilev,itrop(npoints)
   integer match(npoints,nlev-1)
   integer nmatch(npoints)
   integer levmatch(npoints)

   !variables needed for water vapor continuum absorption
   real fluxtop_clrsky(npoints),trans_layers_above_clrsky(npoints)
   real taumin(npoints)
   real dem_wv(npoints,nlev), wtmair, wtmh20, Navo, grav, pstd, t0
   real press(npoints), dpress(npoints), atmden(npoints)
   real rvh20(npoints), wk(npoints), rhoave(npoints)
   real rh20s(npoints), rfrgn(npoints)
   real tauwv(npoints)

   integer icycle
   real tb(npoints)
   real emcld(npoints)
   real fluxtop(npoints)
   real trans_layers_above(npoints)
   real fluxtopinit(npoints),tauir(npoints)

   real rec2p13,tauchk

   ! Specific to GEM implementation
   real tmpexp2D(npoints),tmpexp3D(npoints,nlev)
   real tmplog2D(npoints)

   tauchk = -1.*log(0.9999999)
   rec2p13=1./2.13

   !     ---------------------------------------------------

      if (top_height .eq. 1 .or. top_height .eq. 3) then

      do j=il1, il2 !1,npoints
          ptrop(j)=5000.
          atmin(j) = 400.
          attropmin(j) = 400.
          atmax(j) = 0.
          attrop(j) = 120.
          itrop(j) = 1
      enddo

      do 12 ilev=1,nlev
        do j=il1, il2 !1,npoints
         if (pfull(j,ilev) .lt. 40000. .and. &
                pfull(j,ilev) .gt.  5000. .and. &
                at(j,ilev) .lt. attropmin(j)) then
                ptrop(j) = pfull(j,ilev)
                attropmin(j) = at(j,ilev)
                attrop(j) = attropmin(j)
                itrop(j)=ilev
           end if
           if (at(j,ilev) .gt. atmax(j)) atmax(j)=at(j,ilev)
           if (at(j,ilev) .lt. atmin(j)) atmin(j)=at(j,ilev)
        enddo
12    continue

      end if


!     ---------------------------------------------------
!     COMPUTE CLOUD OPTICAL DEPTH FOR EACH COLUMN and
!     put into vector tau

!initialize tau and albedocld to zero
        do j=il1, il2 !1,npoints
            tau(j)=0.
        enddo

!compute total cloud optical depth for each column
      do ilev=1,nlev
            !increment tau for each of the boxes
         do j=il1, il2 !1,npoints
            if (frac_out(j,ilev).eq.1) then
               tau(j)=tau(j) &
                          + dtau_in(j,ilev)
            endif
         enddo
      enddo ! ilev

!     ---------------------------------------------------
!     COMPUTE INFRARED BRIGHTNESS TEMPERATURES
!     AND CLOUD TOP TEMPERATURE SATELLITE SHOULD SEE

!     again this is only done if top_height = 1 or 3

!     fluxtop is the 10.5 micron radiance at the top of the
!              atmosphere
!     trans_layers_above is the total transmissivity in the layers
!             above the current layer
!     fluxtop_clrsky(j) and trans_layers_above_clrsky(j) are the clear
!             sky versions of these quantities.

      if (top_height .eq. 1 .or. top_height .eq. 3) then


!----------------------------------------------------------------------

!             DO CLEAR SKY RADIANCE CALCULATION FIRST

!compute water vapor continuum emissivity
!this treatment follows Schwarkzopf and Ramasamy
!JGR 1999,vol 104, pages 9467-9499.
!the emissivity is calculated at a wavenumber of 955 cm-1,
!or 10.47 microns
        wtmair = 28.9644
        wtmh20 = 18.01534
        Navo = 6.023E+23
        grav = 9.806650E+02
        pstd = 1.013250E+06
        t0 = 296.

! COMPUTE THE EXP ALL AT ONCE
        do ilev=1,nlev
          do j=il1, il2 !1,npoints
            tmpexp3D(j,ilev) = -0.02*(at(j,ilev)-t0)
          end do
        end do
        call VSEXP(tmpexp3D,tmpexp3D,nlev*(il2-il1+1))

        do 125 ilev=1,nlev
          do j=il1, il2 !1,npoints
!press and dpress are dyne/cm2 = Pascals *10
               press(j) = pfull(j,ilev)*10.
               dpress(j) = (phalf(j,ilev+1)-phalf(j,ilev))*10
!atmden = g/cm2 = kg/m2 / 10
               atmden(j) = dpress(j)/grav
               rvh20(j) = qv(j,ilev)*wtmair/wtmh20
               wk(j) = rvh20(j)*Navo*atmden(j)/wtmair
               rhoave(j) = (press(j)/pstd)*(t0/at(j,ilev))
               rh20s(j) = rvh20(j)*rhoave(j)
               rfrgn(j) = rhoave(j)-rh20s(j)
!               tmpexp(j) = exp(-0.02*(at(j,ilev)-t0))
               tauwv(j) = wk(j)*1.e-20*( &
                 (0.0224697*rh20s(j)*tmpexp3D(j,ilev)) + &
                      (3.41817e-7*rfrgn(j)) )*0.98
               dem_wv(j,ilev) = 1. - exp( -1. * tauwv(j))
          enddo
125     continue

!initialize variables
        do j=il1, il2 !1,npoints
          fluxtop_clrsky(j) = 0.
          trans_layers_above_clrsky(j)=1.
        enddo

! COMPUTE THE EXP ALL AT ONCE
        do ilev=1,nlev
          do j=il1, il2 !1,npoints
            tmpexp3D(j,ilev) = 1307.27/at(j,ilev)
          end do
        end do
        call VSEXP(tmpexp3D,tmpexp3D,nlev*(il2-il1+1))

        do ilev=1,nlev
          do j=il1, il2 !1,npoints

! Black body emission at temperature of the layer

!	        bb(j)=1 / ( exp(1307.27/at(j,ilev)) - 1. )
               bb(j)=1 / ( tmpexp3D(j,ilev) - 1. )
                !bb(j)= 5.67e-8*at(j,ilev)**4

! increase TOA flux by flux emitted from layer
! times total transmittance in layers above

                fluxtop_clrsky(j) = fluxtop_clrsky(j) &
                  + dem_wv(j,ilev)*bb(j)*trans_layers_above_clrsky(j)

! update trans_layers_above with transmissivity
! from this layer for next time around loop

                trans_layers_above_clrsky(j)= &
                  trans_layers_above_clrsky(j)*(1.-dem_wv(j,ilev))


          enddo
        enddo   !loop over level

! COMPUTE THE EXP ALL AT ONCE
          do j=il1, il2 !1,npoints
            tmpexp2D(j) = 1307.27/skt(j)
          end do

        call VSEXP(tmpexp2D,tmpexp2D,(il2-il1+1))

        do j=il1, il2 !1,npoints
!add in surface emission
!          bb(j)=1/( exp(1307.27/skt(j)) - 1. )
          bb(j)=1/( tmpexp2D(j) - 1. )
          !bb(j)=5.67e-8*skt(j)**4

          fluxtop_clrsky(j) = fluxtop_clrsky(j) + emsfc_lw(j) * bb(j) &
           * trans_layers_above_clrsky(j)
        enddo



!           END OF CLEAR SKY CALCULATION

!----------------------------------------------------------------

!loop over columns
        do j=il1, il2 !1,npoints
           fluxtop(j)=0.
           trans_layers_above(j)=1.
        enddo

! COMPUTE THE EXP ALL AT ONCE
        do ilev=1,nlev
          do j=il1, il2 !1,npoints
            tmpexp3D(j,ilev) = 1307.27/at(j,ilev)
          end do
        end do
        call VSEXP(tmpexp3D,tmpexp3D,nlev*(il2-il1+1))

        do ilev=1,nlev
           do j=il1, il2 !1,npoints
! Black body emission at temperature of the layer

!              bb(j)=1 / ( exp(1307.27/at(j,ilev)) - 1. )
              bb(j)=1 / ( tmpexp3D(j,ilev) - 1. )
             !bb(j)= 5.67e-8*at(j,ilev)**4
           enddo

           do j=il1, il2 !1,npoints

! emissivity for point in this layer
              if (frac_out(j,ilev).eq.1) then
                 dem(j)= 1. - &
                      ( (1. - dem_wv(j,ilev)) * (1. -  dem_in(j,ilev)) )
              else
                 dem(j)=  dem_wv(j,ilev)
              end if


! increase TOA flux by flux emitted from layer
! times total transmittance in layers above

              fluxtop(j) = fluxtop(j) &
                         + dem(j) * bb(j) &
                         * trans_layers_above(j)

! update trans_layers_above with transmissivity
! from this layer for next time around loop

              trans_layers_above(j)= &
                   trans_layers_above(j)*(1.-dem(j))

           enddo                ! j
        enddo                   ! ilev

! COMPUTE THE EXP ALL AT ONCE
          do j=il1, il2 !1,npoints
            tmpexp2D(j) = 1307.27/skt(j)
          end do
          call VSEXP(tmpexp2D,tmpexp2D,(il2-il1+1))

        do j=il1, il2 !1,npoints
!add in surface emission
!           bb(j)=1/( exp(1307.27/skt(j)) - 1. )
           bb(j)=1/( tmpexp2D(j) - 1. )
          !bb(j)=5.67e-8*skt(j)**4
        end do

        do j=il1, il2 !1,npoints

!add in surface emission

           fluxtop(j) = fluxtop(j) &
                      + emsfc_lw(j) * bb(j) &
                      * trans_layers_above(j)

        end do


!now that you have the top of atmosphere radiance account
!for ISCCP procedures to determine cloud top temperature

!account for partially transmitting cloud recompute flux
!ISCCP would see assuming a single layer cloud
!note choice here of 2.13, as it is primarily ice
!clouds which have partial emissivity and need the
!adjustment performed in this section

!If it turns out that the cloud brightness temperature
!is greater than 260K, then the liquid cloud conversion
!factor of 2.56 is used.

!Note that this is discussed on pages 85-87 of
!the ISCCP D level documentation (Rossow et al. 1996)

! COMPUTE THE EXP ALL AT ONCE
          do j=il1, il2 !1,npoints
            tmpexp2D(j) = 1307.27/(attrop(j)-5.)
          end do
          call VSEXP(tmpexp2D(il1),tmpexp2D(il1),(il2-il1+1))

          do j=il1, il2 !1,npoints
!compute minimum brightness temperature and optical depth
!            btcmin(j) = 1. /  ( exp(1307.27/(attrop(j)-5.)) - 1. )
             btcmin(j) = 1. /  ( tmpexp2D(j) - 1. )
          enddo

          do j=il1, il2 !1,npoints
             transmax(j) = (fluxtop(j)-btcmin(j)) &
                         /(fluxtop_clrsky(j)-btcmin(j))

!note that the initial setting of tauir(j) is needed so that
!tauir(j) has a realistic value should the next if block be
!bypassed
             tauir(j) = tau(j) * rec2p13
!             taumin(j) = -1.0*log(max(min(transmax(j),0.9999999),0.001))
             taumin(j) = max(min(transmax(j),0.9999999),0.001)
          enddo

          call VSLOG(taumin(il1),taumin(il1),(il2-il1+1))
          do j=il1, il2 !1,npoints
            taumin(j) = -1.0*taumin(j)
          end do

          if (top_height .eq. 1) then
             do j=il1, il2 !1,npoints
                if (transmax(j) .gt. 0.001 .and. &
                     transmax(j) .le. 0.9999999) then
                   fluxtopinit(j) = fluxtop(j)
                   tauir(j) = tau(j) *rec2p13
                endif
             enddo
             do icycle=1,2

! COMPUTE THE EXP ALL AT ONCE
          do j=il1, il2 !1,npoints
            tmpexp2D(j) = -1. * tauir(j)
          end do
          call VSEXP(tmpexp2D(il1),tmpexp2D(il1),(il2-il1+1))

                do j=il1, il2 !1,npoints
                   if (tau(j) .gt. (tauchk            )) then
                      if (transmax(j) .gt. 0.001 .and. &
                           transmax(j) .le. 0.9999999) then
!                         emcld(j) = 1. - exp(-1. * tauir(j)  )
                         emcld(j) = 1. - tmpexp2D(j)
                         fluxtop(j) = fluxtopinit(j) - &
                              ((1.-emcld(j))*fluxtop_clrsky(j))
                         fluxtop(j)=max(1.E-06, &
                              (fluxtop(j)/emcld(j)))
                         tb(j)= 1307.27 &
                                   / (log(1. + (1./fluxtop(j))))
                         if (tb(j) .gt. 260.) then
                            tauir(j) = tau(j) / 2.56
                         end if
                      end if
                   end if
                enddo
             enddo
          endif

! COMPUTE THE LOG ALL AT ONCE (CLEAR AND CLOUDY SKY)
! CLOUDY COLUMNS
          do j=il1, il2 !1,npoints
            tmpexp2D(j) = 1. + (1./fluxtop(j))
          end do
          call VSLOG(tmpexp2D(il1),tmpexp2D(il1),(il2-il1+1))

! CLEAR COLUMNS
          do j=il1, il2 !1,npoints
            tmplog2D(j) = 1. + (1./fluxtop_clrsky(j))
          end do
          call VSLOG(tmplog2D(il1),tmplog2D(il1),(il2-il1+1))

          do j=il1, il2 !1,npoints
             if (tau(j) .gt. (tauchk            )) then
!cloudy box
!                tb(j)= 1307.27/ (log(1. + (1./fluxtop(j))))
                tb(j)= 1307.27/ tmpexp2D(j)
                if (top_height.eq.1.and.tauir(j).lt.taumin(j)) then
                   tb(j) = attrop(j) - 5.
                   tau(j) = 2.13*taumin(j)
                end if
             else
!clear sky brightness temperature
!                tb(j) = 1307.27/(log(1.+(1./fluxtop_clrsky(j))))
                tb(j) = 1307.27/tmplog2D(j)
             end if
          enddo                 ! j

      end if

!     ---------------------------------------------------


!     ---------------------------------------------------
!     DETERMINE CLOUD TOP PRESSURE

!     again the 2 methods differ according to whether
!     or not you use the physical cloud top pressure (top_height = 2)
!     or the radiatively determined cloud top pressure (top_height = 1 or 3)


!segregate according to optical thickness
      if (top_height .eq. 1 .or. top_height .eq. 3) then
!find level whose temperature
!most closely matches brightness temperature
         do j=il1, il2 !1,npoints
            nmatch(j)=0
         enddo
         do 29 ilev=1,nlev-1
!cdir nodep
            do j=il1, il2 !1,npoints
               if ((at(j,ilev)   .ge. tb(j) .and. &
                    at(j,ilev+1) .lt. tb(j)) .or. &
                    (at(j,ilev) .le. tb(j) .and. &
                    at(j,ilev+1) .gt. tb(j))) then

                  nmatch(j)=nmatch(j)+1
                  if(abs(at(j,ilev)-tb(j)) .lt. &
                       abs(at(j,ilev+1)-tb(j))) then
                     match(j,nmatch(j))=ilev
                  else
                     match(j,nmatch(j))=ilev+1
                  end if
               end if
            enddo
 29      continue

         do j=il1, il2 !1,npoints
            if (nmatch(j) .ge. 1) then
               ptop(j)=pfull(j,match(j,nmatch(j)))
               levmatch(j)=match(j,nmatch(j))
            else
               if (tb(j) .lt. atmin(j)) then
                  ptop(j)=ptrop(j)
                  levmatch(j)=itrop(j)
               end if
               if (tb(j) .gt. atmax(j)) then
                  ptop(j)=pfull(j,nlev)
                  levmatch(j)=nlev
               end if
            end if
         enddo                  ! j

      else                      ! if (top_height .eq. 1 .or. top_height .eq. 3)

         do j=il1, il2 !1,npoints
            ptop(j)=0.
         enddo
         do ilev=1,nlev
            do j=il1, il2 !1,npoints
               if ((ptop(j) .eq. 0. ) &
                    .and.(frac_out(j,ilev) .ne. 0)) then
                  ptop(j)=pfull(j,ilev)
                  levmatch(j)=ilev
               end if
            end do
         end do
      end if

      do j=il1, il2 !1,npoints
         if (tau(j) .le. (tauchk            )) then
            ptop(j)=0.
            levmatch(j)=0
         endif
      enddo

      return
      end
