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

subroutine veglowhigh_ccilceco(fcover, vegf_evol, low, high, deci, ever, impervu, &
     ni, nclass)
  use svs_configs, only : ntypel, vl_type, ntypeh, vh_type, furb_vl, imp_urb, epsilon_svs
   implicit none
!!!#include <arch_specific.hf>

   !@Arguments
   !            - Input -
   ! NI         Horizontal dimension
   ! NCLASS     Number of natural landuse classes
   ! FCOVER     Fraction of vegetation (masks)
   ! VEGF_EVOL  Fraction of vegetation (masks) x sparsness (vegdat)
   !            - Output -
   ! LOW        Fraction of (soil) grid covered by low vegetation classes 
   ! HIGH       Fraction of (soil) grid covered by high vegetation classes
   ! DECI       Fraction of HIGH vegetation that is deciduous
   ! EVER       Fraction of HIGH vegetation that is evergreen 
   ! IMPERVU    Fraction of LAND SURFACE that is IMPERVIOUS 


   integer ni, nclass
   real deci(ni), ever(ni), fcover(ni,nclass), vegf_evol(ni,nclass), impervu(ni)
   real high(ni), low(ni)
!
!Author
!        Maria Abrahamowicz (dec 2020)
!Revision
! 001    
!
!Object
!        Compute fraction of tile covered by low and high vegetation respectively

!     !change size of V*_TYPE() arrays below
      integer i, m, k, type

      REAL totfract  
!
!
      DO i=1,ni
        low(i) = 0.0
        high(i)= 0.0
        deci(i)= 0.0
        ever(i)= 0.0
      END DO
!
!
!
      DO i=1,ni
        totfract = 0.

        DO m=4,nclass                        ! loop on veg/natural soil classes
          totfract = totfract + fcover(i,m)
        END DO
!
        IF (totfract.GE.EPSILON_SVS) THEN

           DO k=1,ntypel  
              type=vl_type(k)  ! loop on LOW vegetation classes
!              
              low(i)  = low(i) + vegf_evol(i,type)
           END DO

           DO k=1,ntypeh  
              type=vh_type(k)  ! loop on HIGH vegetation classes
!
              high(i) = high(i) + vegf_evol(i,type)
          END DO
!
!
          !   ADD URBAN CLASS (=21),  between HIGH AND LOW VEG
          type=21
!
          low(i)  = low(i)  + furb_vl * vegf_evol(i,type)
          high(i) = high(i) + (1. - furb_vl) * vegf_evol(i,type)

          !   HIGH AND LOW VEGETATION FRACTION
          low(i)  = low(i)  / totfract
          high(i) = high(i) / totfract

          !   FRACTION OF LAND SURFACE THAT IS CONSIDERED IMPERVIOUS 
          ! TAKE TO BE FRACTION OF URBAN CLASS (vf=21) TIMES A CONSTANT...
          impervu(i) =  imp_urb * fcover(i,21) / totfract
!

       ELSE
          low(i)     = 0.0
          high(i)    = 0.0
          impervu(i) = 0.0

       END IF

!              Calculate deciduous and evergreen fractions
!          *** ADD VF18--- usually unused... so should be ok
!              high veg: type: 4, 5, 6, 7, 8, 9, 18, 19, 25, 26
!              decidous: type: 6,7,8,9, half of 18, 25,26
!              evergreen:type: 4,5, half of 18, 19, 25,26
       !              ** plus fraction of urban... ***

       !##### note: hard-coding deciduous/evergreen to match CCILECO
       !# calculation should be automatic
       ! probably need to defined 3 array:
       ! decidous, evegreen, mixed


       
        IF(high(i).ge.EPSILON_SVS)THEN
!
           
           deci(i) = ( vegf_evol(i,6)  &
                   +   vegf_evol(i,7)  &
                   +   vegf_evol(i,8)  & 
                   +   vegf_evol(i,9)  &
                   +  (1. - furb_vl) * vegf_evol(i,21) & 
                   +  0.5 * vegf_evol(i,18)  &
                   +  0.5 * vegf_evol(i,25) & 
                   +  0.5 * vegf_evol(i,26) ) /high(i)
!        
           ever(i) = 1.0 - deci(i)
        else
           deci(i)=0.0
           ever(i)=0.0
!
        ENDIF

!
      END DO
!
!
!
   return
 end subroutine veglowhigh_ccilceco
