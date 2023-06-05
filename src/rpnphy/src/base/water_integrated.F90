!-------------------------------------- LICENCE BEGIN -----------------------
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
!-------------------------------------- LICENCE END -------------------------

module water_integrated
   implicit none
   private
   public :: wi_integrate

contains

   !/@*
  subroutine wi_integrate(t, hu, lwc, iwc, s, ps, &
       icw, iwv, iwv700, iwp, lwp, slwp, slwp2, slwp3, slwp4, &
       n, nk)
    use debug_mod, only: init2nan
    use tdpack_const, only: TCDK, GRAV
    use phy_options
    implicit none

    !@Author
    !          R.Sarrazin, G. Pellerin, B. Bilodeau - (Sept 1996)
    !
    !@Revision
    ! 001      A. Glazer (Nov 1999) - Intwat1: add the vertical integrals of
    !          ice, supercooled liquid water, supercooled liquid water
    !          by layers, cloud-top pressure and cloud-top temperature.
    ! 002      A. Glazer (Feb 2000) - Intwat2: generalize to any scheme that
    !          calculates liquid water content and ice water content (inputs).
    ! 003      A. Glazer (Nov 2001) - Add the vertical integral of vater vapor
    !          from the model top to 700 mb.
    ! 004      B. Bilodeau and P. Vaillancourt (Dec 2002) - Remove ctp and ctt
    !
    !@Object
    !          to compute the vertical integrals of water vapor,
    !          liquid water and cloud water.

    integer, intent(in) :: n, nk
    real, dimension(:,:), intent(in) :: t, hu, lwc, iwc, s
    real, dimension(:), intent(in) :: ps
    real, dimension(:), intent(out) :: icw, iwv, iwv700, iwp, lwp, &
         slwp, slwp2, slwp3, slwp4

    !@Arguments
    !          - Input -
    ! t        temperature
    ! hu       humidity
    ! lwc      total specific liquid content t+dT
    ! iwc      total specific solid content at t+dT
    ! ps       surface pressure
    ! s        vertical coordinate
    ! n        horizontal running length
    ! nk       vertical dimension
    !          - Output -
    ! icw    - integrated cloud water/ice [kg/m2]
    ! iwv    - integrated water vapor [kg/m2]
    ! iwv700 - integrated water vapor (0-700 mb) [kg/m2]
    ! iwp    - integrated ice water [kg/m2]
    ! lwp    - integrated cloud liquid water [kg/m2]
    ! slwp   - integrated SLW (supercooled liquid water) [kg/m2]
    ! slwp2  - integrated SLW (bottom to s2) [kg/m2]
    ! slwp3  - integrated SLW (s2 to s3) [kg/m2]
    ! slwp4  - integrated SLW (s3 to s4) [kg/m2]
    !*@/
    
    integer i,k,im,k1,k2,k3,k4
    real dsg,dpsg,qctemp
    real :: s1=1., s2=0.8, s3=0.6, s4=0.4
    logical :: scool=.false.
#include <rmn/msg.h>

    !----------------------------------------------------------------
    call msg_toall(MSG_DEBUG, 'water_integrated [BEGIN]')
    
    !* 1. INITIALIZE OUTPUT FIELDS
    do i = 1,N
       ICW (i)   = 0.
       IWV (i)   = 0.
       IWV700(i) = 0.
       IWP (i)   = 0.
       LWP (i)   = 0.
       SLWP(i)   = 0.
       SLWP2(i)  = 0.
       SLWP3(i)  = 0.
       SLWP4(i)  = 0.
    end do

    !* 3. REINITIALIZE FIELDS

    do i = 1,N
       LWP (i) = 0.
       IWP (i) = 0.
    end do

    !* 4. SUPERCOOLED LIQUID WATER by layers

    if (scool) then

       !     Find k1, k2, k3 and k4 from s1, s2, s3 and s4

       im = int(N/2 + 1)

       do k = 1,NK
          if (s(im,k) .gt. s4) go to 200
       end do
200    continue
       k4 = k-1

       do k = k4+1,NK
          if (s(im,k) .gt. s3) go to 300
       end do
300    continue
       k3 = k-1

       do k = k3+1,NK
          if (s(im,k) .gt. s2) go to 400
       end do
400    continue
       k2 = k-1

       if (s1 .eq. 1.) then
          k1 = NK
       else
          do k = k1+1,NK
             if (s(im,k) .gt. s1) go to 500
          end do
500       continue
          k1 = min(NK,k-1)
       end if

       do i = 1,N
          do k = k2+1,k1
             if ( T(i,k) .lt. tcdk) then
                dsg = 0.5 * ( s(i,min(k+1,NK)) - s(i,max(k-1,1)) )
                if (k.eq.NK) dsg = 1. - 0.5 * ( s(i,NK) + s(i,NK-1))
                dpsg = max(ps(i)*dsg/grav,0.)
                qctemp  = (max( LWC(i,k), 0. ))*dpsg
                SLWP2 (i)= SLWP2(i) + qctemp
             end if
          end do
       end do

       do i = 1,N
          do k = k3+1,k2
             if ( T(i,k) .lt. tcdk) then
                dsg = 0.5 * ( s(i,min(k+1,NK)) - s(i,max(k-1,1)) )
                if (k.eq.NK) dsg = 1. - 0.5 * ( s(i,NK) + s(i,NK-1))
                dpsg = max(ps(i)*dsg/grav,0.)
                qctemp  = (max( LWC(i,k), 0. ))*dpsg
                SLWP3 (i)= SLWP3(i) + qctemp
             end if
          end do
       end do

       do i = 1,N
          do k = k4+1,k3
             if ( T(i,k) .lt. tcdk) then
                dsg = 0.5 * ( s(i,min(k+1,NK)) - s(i,max(k-1,1)) )
                if (k.eq.NK) dsg = 1. - 0.5 * ( s(i,NK) + s(i,NK-1))
                dpsg = max(ps(i)*dsg/grav,0.)
                qctemp  = (max( LWC(i,k), 0. ))*dpsg
                SLWP4 (i)= SLWP4(i) + qctemp
             end if
          end do
       end do
    end if

    !* 5. VERTICAL INTEGRALS of TOTAL CONDENSATE (ICW), VAPOR (IWV),
    !     VAPOR from top to 700 mb (IWV700),SOLID (IWP),
    !     LIQUID (LWP) and SUPERCOOLED (SLWP)

    do i = 1, N
       dsg= 0.5 * ( s(i,2) - s(i,1) )
       dpsg= ps(i)*dsg/grav
       IWV(i) = IWV(i) + max( HU(i,1) , 0. ) * dpsg
       IWV700(i) = IWV700(i) + max( HU(i,1) , 0. ) * dpsg
       qctemp = max(LWC(i,1) , 0. ) * dpsg
       LWP (i)= LWP(i) + qctemp
       if (T(i,1).lt.tcdk) SLWP (i) =  SLWP(i) + qctemp
       IWP(i) = IWP(i) + max(IWC(i,1) , 0. ) * dpsg
    end do

    do k = 2,NK-1
       do i = 1 , N
          dsg= 0.5 * ( s(i,k+1) - s(i,k-1) )
          dpsg= ps(i)*dsg/grav
          IWV(i) = IWV(i) + max( HU(i,k) , 0. ) * dpsg
          if ((s(i,k)*ps(i)) .lt. 70000.)IWV700(i)=IWV700(i)+max(HU(i,k),0.)*dpsg
          qctemp = max(LWC(i,k) , 0. ) * dpsg
          LWP (i)= LWP(i) + qctemp
          if (T(i,k).lt.tcdk) SLWP (i)= SLWP(i) + qctemp
          IWP(i) = IWP(i) + max(IWC(i,k) , 0. ) * dpsg
       end do
    end do

    do i = 1, N
       dsg= 1. - 0.5 * ( s(i,NK) + s(i,NK-1) )
       dpsg= ps(i)*dsg/grav
       IWV(i) = IWV(i) + max( HU(i,NK) , 0. ) * dpsg
       qctemp = max(LWC(i,NK) , 0. ) * dpsg
       LWP (i)= LWP(i) + qctemp
       if (T(i,NK).lt.tcdk) SLWP (i)= SLWP(i) + qctemp
       IWP(i) = IWP(i) + max(IWC(i,NK) , 0. ) * dpsg
    end do

    do i = 1, N
       ICW(i)  = LWP(i) + IWP(i)
    end do

    call msg_toall(MSG_DEBUG, 'water_integrated [END]')
    !----------------------------------------------------------------
    return
  end subroutine wi_integrate

end module water_integrated
