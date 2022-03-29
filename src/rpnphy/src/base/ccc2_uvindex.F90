!---------------------------------- LICENCE BEGIN -------------------------------
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
!---------------------------------- LICENCE END ---------------------------------

!/@*
subroutine ccc2_uvindex2(il1,il2,ilg,fctb,fatb,salb,iuvc,iuva,iuv_method)
   use phy_options, only: RAD_NUVBRANDS, IUV_METHOD_OPT
   implicit none
   ! integer, parameter :: rad_nuvbrands=6
   integer :: il1,il2,ilg
   real, dimension(ilg,RAD_NUVBRANDS) :: fatb,fctb
   real, dimension(ilg)               :: iuvc,iuva,salb
   character(len=10), optional        :: iuv_method

   !@Authors  Y.J.Rochon, Feb 2016, Oct 2017, Oct-Dec 2018
   !
   !@Object Calculate the clear-sky and all-sky UV index 
   !         values from GEM solar UV irradiances at the surface. 
   !  
   !@Arguments
   !          - input/output - 
   ! il1      applied start index of horizontal dimension (usually 1)
   ! il2      applied end index of horizontal dimension 
   ! ilg      max horizontal dimension
   ! fctb     Clear-sky downward irradiance at the surface (dir+dif) flux for 6 VIS-UV bands (input)
   ! fatb     all-sky downward irradiance at the surface (dir+dif) flux for 6 VIS-UV bands (input)
   ! salb     Surface albedo (input)
   ! iuvc     clear-sky UV index field(s) (input for uvi_method='BandRatio'; output otherwise)
   ! iuva     all-sky UV index field(s) (output)
   ! iuv_method  
   !          Index of UV index calc method 
   !          Skips IUV* calc if set to anything other than the three following options.
   !             'IntegFit'  - Integration over the four UV bands with integrand weighting
   !                           using erythemal function. Calc of both clear-sky and all-sky
   !                           values. (assumes interactive ozone)
   !             'LinearFit' - Fitted sum of the four UV broadband irradiances. Calc of 
   !                           both clear-sky and all-sky values. (assumes interactive ozone)
   !             'BandRatio' - Scaling of input clear-sky UV index by the ratio of all-sky
   !                           to clear-sky 310-320 nm band irradiances.
   !
   !@Notes 
   !         1) Fit models and coefficients described by Tereszchuk et al. (2017) on
   !            'Optimizing UV index determination from braodband irradiances'
   !            with differences described in the declaration sections below preceded by
   !            '**** iuv_method=' 
   !         2) The 'IntegFit' (and 'BandRatio') would give slightly smaller systematic errors
   !            (with UV index differences of 0.5) as compared to the 'LinearFit' approach.
   !         3) For iuv_method='BandRatio', 
   !              IUVC:
   !                - iuvc is calculated via 'Integfit' if not available on input (entire field is zero).
   !                - otherwise, the input iuvc will be used.
   !              IUVA:
   !                - The ratio of the all-sky to clear-sky irradiances reflects the
   !                  attenuation impact of clouds. The 310.7-320nm band was selected
   !                  for applying the attenuation.
   !         4) The differences with ccc2_uvindex.F90 would be the fit coefficients 
   !            and the band intervals, the former having depended on the latter.
   !         5) Added additional dependence on albedo based on 7-band UV case and
   !            comparison to high resolution simulations.
   !
   !
   ! Comments:
   !
   ! 1. Routines 'uvindex_integ' and 'uv_interp_wgts' generalized to handle
   !    any band specifications in the declaration section of 'ccc2_uvindex'.
   !    
   !============================================================================
   !*@/
   
   ! Work variable declaration

   character(len=10) :: imethod

   integer, parameter :: nftb=4
   ! Spectral boundaries (in nm) for the input spectral bands      
   real, parameter :: radbin_bound(nftb+1)=(/ 280.11, 292.40, 310.70, 320.0, 400.0 /) 

   ! Specification of desired bands from input f*tb arrays.
   integer, parameter :: iftb(nftb)=(/ 6,5,4,3 /)

   real :: fatbw(ilg,nftb),fctbw(ilg,nftb)   
   integer :: i,j
   integer :: ifirst
   real :: rad(nftb)

   ! **** uvi_method='LinearFit' ****
   ! UV Index linear fit factors
   ! 
   ! Fit coefficients
   ! Coefficients 1 and 4 (values of 40.0 and 0.02 were set in consideration of
   ! the Erythemal function in the respective regions, with the other two 
   ! coefficients determined via least-squares minimization. 
   ! Note that various combinations of coefficient values give similar results.
   real, parameter :: linear_fit_factors(nftb)=(/  40.0, 9.618, 0.431, 0.02 /)

   ! **** uvi_imethod='IntegFit' ****
   ! Work variables and arrays for spectral integrations 
   ! for routine uvindex_integ and its related subroutine uvindex_interp_wgts

   real :: rad2(nftb),radw(nftb) 

   ! Set reference positions assigned to rad2.
   ! Minimization solution based on single 6 hour forecast fields valid for
   ! 18 UTC (clear-sky) on 8 August, 2016.
   ! 
   ! Fitted values of broadband reference nodes
   ! (bounded by band boundaries)
   real, parameter :: xm(nftb)=(/ 280.2, 302.7, 314.1, 397.6 /)

   integer, parameter :: nsegments=6  ! Number of integtration segments

   ! Spectral boundaries of integration segments   
   real, parameter :: seg_bound(nsegments+1) = &
        (/ 280.11, 292.40, 298.0, 310.70, 320.0, 328.0, 400.0 /)

   ! Starting input spectral band for use in each integration segment
   integer, parameter :: seg_startband(nsegments) = (/ 1,1,1,2,2,2 /)

   ! Order of Laggrange poly for each integration segment
   ! iorder (3 and 7) =1  ! Cannot be 2 as it will generate overshots at sunset/sunrise.
   integer, parameter :: iorder(nsegments) = (/ 1,1,1,1,1,1 /)

   ! Lagrange poly weight
   real :: wgts(nsegments,5,3)
   ! Erythemal function values evaluated at interpolation points
   real :: erythemal(nsegments,5)

   ! **** Above for uvi_method='IntegFit' ****

   ! ----------------------------------------------------------------------

   imethod='INTEGFIT'
   if (present(iuv_method)) then
      if (any(iuv_method == IUV_METHOD_OPT)) then
         imethod = trim(iuv_method) 
      else
         iuvc=0.0
         iuva=0.0
         return
      end if
   end if

   ! Scale fluxes as needed

   if (imethod /= 'BANDRATIO') then      
      do i=1,nftb
         fctbw(il1:il2,i) = fctb(il1:il2,iftb(i))
         call uvindex_scale_tflux(fctbw(il1:il2,i),salb(il1:il2),i)
         where (fctb(il1:il2,iftb(3)) >= 0.0001 .and. &
              fatb(il1:il2,iftb(3)) >= 1.E-8 .and. &
              fctb(il1:il2,iftb(i)) > 0.1*fatb(il1:il2,iftb(i))) 
            fatbw(il1:il2,i)=fctbw(il1:il2,i)*fatb(il1:il2,iftb(i))/fctb(il1:il2,iftb(i))
         elsewhere
            fatbw(il1:il2,i)=0.0
         endwhere
      end do
   else
      if (all(iuvc < 0.1)) then
         do i=1,nftb
            fctbw(il1:il2,i) = fctb(il1:il2,iftb(i))
            call uvindex_scale_tflux(fctbw(il1:il2,i),salb(il1:il2),i) 
         end do
      else
         i=3
         fctbw(il1:il2,i) = fctb(il1:il2,iftb(i))
         call uvindex_scale_tflux(fctbw(il1:il2,i),salb(il1:il2),i) 
      end if
   end if

   !  Calc local solar time dependant UV index for clear-sky and all-sky conditions

   ifirst=0
   select case(imethod)
   case('INTEGFIT')

      ! Apply integration approach   

      iuva(:) = 0.0
      iuvc(:) = 0.0
      do j=il1,il2
         if (fctb(j,iftb(3)) < 0.0001) cycle 
         do i=1,nftb
            rad(i) = fctbw(j,i)
         end do
         iuvc(j) = uvindex_integ()
         do i=1,nftb
            rad(i) = fatbw(j,i)
         end do
         iuva(j) = uvindex_integ()
      end do

   case('LINEARFIT')

      ! Apply weighted linear sum of the broadbands

      iuva(:) = 0.0
      iuvc(:) = 0.0
      do j=il1,il2
         if (fctb(j,iftb(3)) < 0.0001) cycle 
         iuvc(j) = linear_fit_factors(1)*fctbw(j,1) + linear_fit_factors(2)*fctbw(j,2) &
              + linear_fit_factors(3)*fctbw(j,3) + linear_fit_factors(4)*fctbw(j,4)
         iuva(j) = linear_fit_factors(1)*fatbw(j,1) + linear_fit_factors(2)*fatbw(j,2) &
              + linear_fit_factors(3)*fatbw(j,3) + linear_fit_factors(4)*fatbw(j,4)
      end do

   case('BANDRATIO')

      if (all(iuvc < 0.1)) then

         ! Apply integration approach for iuvc only

         iuvc(:) = 0.0
         do j=il1,il2
            if (fctb(j,iftb(3)) < 0.0001) cycle 
            do i=1,nftb
               rad(i) = fctbw(j,i)
            end do
            iuvc(j) = uvindex_integ()
         end do
      end if

      ! Generate iuva given iuvc.

      iuva(:) = 0.0
      do j=il1,il2
         if (fctb(j,iftb(3)) < 0.0001 .or. fatb(j,iftb(3)) < 1.E-8) cycle 
         iuva(j) = iuvc(j)*fatb(j,iftb(3))/fctb(j,iftb(3))
      end do

   end select

contains

   !===========================================================================
   subroutine uvindex_scale_tflux(f,salb,index)
      ! Creation       : Y.J. Rochon, Dec 2015, Oct-Nov 2018 
      !
      ! Modified: 
      !
      ! Description    : Scale GEM UV total fluxes at the surface following
      !                  fit to CloudJ-based broadband irradiances (e.g. Tereszchuk et al, 2017)
      !                 
      ! Arguments:  
      !
      !            IN 
      !                index      - Broadband index 
      !                salb       - Surface albedo
      !
      !            IN/OUT   
      !                f          - Total irradiances
      !========================================================================
      implicit none
      integer,intent(in) :: index
      real, intent(inout) :: f(:),salb(:)

      real, parameter :: EPSILON = max(tiny(1.0), 1.0e-15)
      real, parameter :: fit_factors(4) = (/ 0.603, 1.618, 0.959, 1.008 /)
      real, parameter :: fit_expnt(4) = (/ 0.575, 0.605, 0.0, 0.0 /)

      ! Gradients of fit_factors as a function of salb derived from values for
      ! albedos of 0.1 and 1.0.
      ! Note: changes of exponent as a function of salb was found to be weaker.

      real, parameter :: fit_grad_salb(4)=(/ 0.31, 0.25, 0.060, 0.020 /)

      ! Apply scaling
      
      select case(index)
      case(1)
         where (f <= EPSILON)
            f = 0.
         elsewhere
            f = (fit_factors(1)+fit_grad_salb(1)*(salb-0.1d0))*f**fit_expnt(1)  
         endwhere
      case(2)
         where (f <= EPSILON)
            f = 0.
         elsewhere
            f = (fit_factors(2)+fit_grad_salb(2)*(salb-0.1d0))*f**fit_expnt(2)
         endwhere
      case(3)
         where (f <= EPSILON)
            f = 0.
         elsewhere
            f = (fit_factors(3)+fit_grad_salb(3)*(salb-0.1d0))*f
         endwhere
      case(4)
         where (f <= EPSILON)
            f = 0.
         elsewhere
            f = (fit_factors(4)+fit_grad_salb(4)*(salb-0.1d0))*f
         endwhere
      end select

   end subroutine uvindex_scale_tflux

   !===========================================================================
   real function uvindex_integ()
      ! Creation       : Y.J. Rochon, Nov 2014 
      !
      ! Modified:
      !
      ! Description    : Calc UV index from GEM broadband radiances (in W/m^2)
      !
      ! Extra info:
      !
      ! 1. Approach hardcoded for GEM boadband radiances for increased
      !    efficiency.
      !
      ! 2. Adjustments may be required later depending on comparisons
      !    to results with high spectral resolution radiances and, at  
      !    least for clear sky conditions, to observations.
      !
      ! 3. Spectral ranges of input irradiances (in W/m^2)
      !
      !    rad(1)=FCTB6 (280.11-292.40nm; midpoint=286.26)
      !    rad(2)=FCTB5 (292.40-310.70nm; midpoint=301.55)
      !    rad(3)=FCTB4 (310.70-320.00nm; midpoint=315.35)
      !    rad(4)=FCTB3 (320.00-400.00nm; midpoint=360.00)
      !
      ! Arguments:  IN     
      !                rad(4)      - GEM UV broadband irradiances covering
      !                              the spectral ranges:
      !                               
      !                              280.11-292.40nm, 292.40-310.70nm, 
      !                              310.70-320.00nm, and 320.00-400.00nm
      !
      !                ifirst      - index for initializing weighting coefficients
      !                              0 - set/calc coefficients
      !                              1 - assumed already available
      !
      !             OUT
      !                uvindex_integ 
      !                ifirst      - value of 1
      !========================================================================
      implicit none

      !  Declare local variables
      
      real, parameter :: scale= 40. ! =1/(25E-3 W/m^2) for conversion
      ! from irradiances to UV Index

      integer :: i,k,ik
      real :: y(5),sumv

      !  ----------------------------------------------------------------------   

      !  Set interpolation weights during first call only

      if (ifirst == 0) then
         call uvindex_interp_wgts()
         ifirst = 1
      end if

      !  Convert to W/(m^2 nm) and take log for interpolations

      radw(:) = rad(:)
      where (radw(:) < 1.E-7) radw(:) = 1.E-7
      do i=1,nftb
         rad2(i) = log(radw(i)/(radbin_bound(i+1)-radbin_bound(i)))
      end do

      !  Integration covering seg_bounds(1) to seg_bounds(2) below 298 nm

      sumv = radw(1)

      do i=2,nsegments

         if (seg_bound(i+1) < 298.1) then

            !  Simpson's rule applied to one 3-point interval in the segment
            !  Note: 0.1667= 1/(number of subintervals) * 1/3, number of subintervals = 2

            k = seg_startband(i)
            if (i > 2) then
               y(1) = y(3)
            else
               y(1) = erythemal(i,1)*exp(wgts(i,1,1)*rad2(k)+wgts(i,1,2)*rad2(k+1)+wgts(i,1,3)*rad2(k+2))          
            end if

            do ik=2,3
               y(ik) = erythemal(i,ik)*exp(wgts(i,ik,1)*rad2(k)+wgts(i,ik,2)*rad2(k+1)+wgts(i,ik,3)*rad2(k+2))
            end do

            sumv = sumv+0.1667*(seg_bound(i+1)-seg_bound(i))*(y(1)+4.*y(2)+y(3))

         else  

            !  Simpson's rule applied to two 3-point intervals in the segment
            !  with a common interval interface point (each interval has 2 segments subintervals)
            !  Note: 0.08333= 1/(total number of intervals) * 1/3, total number of intervals = 4

            k = seg_startband(i)
            if (i > 2) then 
               if (seg_bound(i) < 298.1) then

                  !  Aside: Note the discontinuity at w=298nm as compared to previous
                  !  segment. This is a weakness of the erythemal function.

                  y(1) = 1.242*y(3)
               else
                  y(1) = y(5)
               end if
            else
               y(1) = erythemal(i,1)*exp(wgts(i,1,1)*rad2(k)+wgts(i,1,2)*rad2(k+1)+wgts(i,1,3)*rad2(k+2))           
            end if

            do ik=2,5
               y(ik) = erythemal(i,ik)*exp(wgts(i,ik,1)*rad2(k)+wgts(i,ik,2)*rad2(k+1)+wgts(i,ik,3)*rad2(k+2))
            end do

            sumv = sumv+0.08333*(seg_bound(i+1)-seg_bound(i))*(y(1)+4.*y(2)+2.*y(3)+4.*y(4)+y(5)) 
         end if
      end do

      !  Final conversion to UV index.

      uvindex_integ = sumv*scale ! Division by 25E-3 W/m^2 for conversion of UV Index
   end function uvindex_integ

   !===========================================================================
   subroutine uvindex_interp_wgts()
      !
      ! Creation       : Y.J. Rochon, Nov 2014 
      ! Modified:
      !
      ! Description    : Set interpolation wgts for the following for
      !                  segments to 3 to 4 interpolation points depending on regions
      !                  of 'seg_bounds' using first (iorder=1) or second (iorder=2) order 
      !                  interpolation. 
      !  
      !                  Note: The first node of each of the  2+ regions is the last node
      !                  node of the previous region and so its interpolation weights
      !                  need not be re-calculated. 
      !
      ! Extra info:
      !
      ! 1. See calling function uvindex_integ for additional info.
      !
      ! 2. Sample references for erythemal action spectrum function:
      !
      !    http://journal.cpha.ca/index.php/cjph/article/view/1905
      !    http://files.cie.co.at/724_cie209_2014.pdf
      !
      ! Arguments:  
      !             IN     
      !                xm(nftb)          - Spectral position of the four original points
      !                iorder(nsegments) - Order of Lagrange poly.
      !
      !             OUT
      !                wgts(nsegments,4,3)  - Weights (spectral integration region, 
      !                                       interpolation point/node in region,
      !                                       subset of original points)
      !                erythemal(nsegments,4) - Erythemal action spectrum function values for
      !                                         (spectral integration region, 
      !                                          interpolation point/node in region)
      !========================================================================
      implicit none

      !  Declare local variables

      integer :: i,j,ik,k,n
      real    :: xpos(5)

      !  Initialize weights

      wgts(:,:,:) = 0.0
      erythemal(:,:) = 0.0

      do j=2,nsegments

         ! Set number of input positions from the order of Lagrange poly.
         n = iorder(j)+1

         if (seg_bound(j+1) < 298.1) then

            ! Set 3 equally spaced interpolation points

            do i=1,3
               xpos(i) = seg_bound(j)+(i-1)*(seg_bound(j+1)-seg_bound(j))/2.
            end do

            ! Determine Lagrange interpolation weights for each interpolation point

            call uvindex_lagrange_wgt(xm(1:n),wgts(j,1,1:n),n,xpos(1))
            call uvindex_lagrange_wgt(xm(1:n),wgts(j,2,1:n),n,xpos(2))
            call uvindex_lagrange_wgt(xm(1:n),wgts(j,3,1:n),n,xpos(3))

            ! Set  erythemal function values at interpolation points

            do i=1,3
               erythemal(j,i) = 1.0E0
            end do

         else

            ! Set 5 equally spaced interpolation points

            do i=1,5
               xpos(i) = seg_bound(j)+(i-1)*(seg_bound(j+1)-seg_bound(j))/4.
            end do

            ! Determine Lagrange interpolation weights for each interpolation point

            if (n == 2) then
               do ik=2,5
                  if (xm(seg_startband(j)+1) > xpos(ik)) then
                     k = seg_startband(j)
                     call uvindex_lagrange_wgt(xm(k:k+n-1),wgts(j,ik,1:n),n,xpos(ik))
                  else
                     k = seg_startband(j)+1
                     call uvindex_lagrange_wgt(xm(k:k+n-1),wgts(j,ik,2:3),n,xpos(ik))
                  end if
               end do
               k = seg_startband(j)
               if (j == 2) call uvindex_lagrange_wgt(xm(k:k+n-1),wgts(j,1,1:n),n,xpos(1))
            else
               k=seg_startband(j)
               do ik=2,5
                  call uvindex_lagrange_wgt(xm(k:k+n-1),wgts(j,ik,1:n),n,xpos(ik))
               end do
               if (j == 2) call uvindex_lagrange_wgt(xm(k:k+n-1),wgts(j,1,1:n),n,xpos(1))
            end if

            ! Evaluate erythemal function values at interpolation points

            if (seg_bound(j+1) < 328.1) then
               do i=2,5
                  erythemal(j,i) = 10**(0.094*(298.0-xpos(i)))
               end do
               if (j == 2) erythemal(j,1) = 10**(0.094*(298.0-xpos(1)))
            else
               do i=2,5
                  erythemal(j,i) = 10**(0.015*(140.0-xpos(i)))
               end do
            end if

         end if

      end do

   end subroutine uvindex_interp_wgts

   !===========================================================================
   subroutine uvindex_lagrange_wgt(x,y,n,xp)
      ! Creation       : Y.J. Rochon, Nov 2014 
      ! Modified:
      !
      ! Description    : Set Lagrange interpolation weights y(n) for a point at xp
      !                  according to points at x(n).
      !
      ! Extra info:
      !
      ! Arguments:
      !            IN 
      !                n           - Number of input positions.
      !                              Also, n-1 is the order of Lagrange poly.   
      !                x(n)        - Input positions
      !                xp          - Position for which to calc weight
      !
      !             OUT  
      !                y(n)        - Lagrange poly weights.
      !========================================================================
      implicit none
      integer, intent(in) :: n
      real, intent(in)    :: x(n),xp
      real, intent(out) ::  y(n) 

      !  Declare local variables

      integer :: i,j

      y(1:n) = 1.0

      do j=2,n
         y(1) = y(1)*(xp-x(j))/(x(1)-x(j))
      end do

      do j=1,n-1
         y(n) = y(n)*(xp-x(j))/(x(n)-x(j))
      end do

      if (n > 2) then
         do i=2,n-1
            do j=1,i-1 
               y(i) = y(i)*(xp-x(j))/(x(i)-x(j))
            end do
            do j=i+1,n
               y(i) = y(i)*(xp-x(j))/(x(i)-x(j))
            end do
         end do
      end if

   end subroutine uvindex_lagrange_wgt

end subroutine ccc2_uvindex2

