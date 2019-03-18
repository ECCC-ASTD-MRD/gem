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

subroutine veglowhigh(fcover, tablen, tables, low, high, deci, ever, impervu, &
     lat, ni, nclass)
  use svs_configs, only : ntypel, vl_type, ntypeh, vh_type, furb_vl, imp_urb, epsilon_svs
   implicit none
#include <arch_specific.hf>

   !@Arguments
   !            - Input -
   ! NI         Horizontal dimension
   ! NCLASS     Number of natural landuse classes
   ! FCOVER      Fraction of vegetation (masks)
   ! TABLEN     Geophysical fields values for each type of vegetation (northern hemisphere)
   ! TABLES     Geophysical fields values for each type of vegetation (southern hemisphere)
   ! LAT        Latitude
   !            - Output -
   ! LOW        Fraction of (soil) grid covered by low vegetation classes 
   ! HIGH       Fraction of (soil) grid covered by high vegetation classes
   ! DECI       Fraction of HIGH vegetation that is deciduous
   ! EVER       Fraction of HIGH vegetation that is evergreen 
   ! IMPERVU    Fraction of LAND SURFACE that is IMPERVIOUS 


   integer ni, nclass
   real deci(ni), ever(ni), fcover(ni,nclass), impervu(ni)
   real high(ni), low(ni), tablen(nclass)
   real tables(nclass), lat(ni)
!
!Author
!        Maria Abrahamowicz (dec 2012)
!Revision
! 001    
!
!Object
!        Compute fraction of tile covered by low and high vegetation respectively

!     !change size of V*_TYPE() arrays below
      integer i, m, k, type

      REAL totfract, table_val   
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
            if (lat(i).ge.0.0) then    
!              northern hemisphere
               table_val = tablen(type)
            else
!              southern hemisphere
               table_val = tables(type)
            endif
            low(i)  = low(i) + fcover(i,type)* table_val
          END DO

          DO k=1,ntypeh  
           type=vh_type(k)  ! loop on HIGH vegetation classes
!
            if (lat(i).ge.0.0) then    
!              northern hemisphere
               table_val = tablen(type)
            else
!              southern hemisphere
               table_val = tables(type)
            endif
            high(i) = high(i) + fcover(i,type)* table_val
          END DO
!
!
          !   ADD URBAN CLASS (=21),  between HIGH AND LOW VEG
          type=21
!
          if (lat(i).ge.0.0) then    
!              northern hemisphere
             table_val = tablen(type)
          else
!              southern hemisphere
             table_val = tables(type)
          endif
          low(i)  = low(i)  + furb_vl * fcover(i,type)* table_val
          high(i) = high(i) + (1. - furb_vl) * fcover(i,type)* table_val

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
!              high veg: type: 4, 5, 6, 7, 8, 9, 25, 26
!              decidous: type: 6,7,8,9, half of 25,26
!              evergreen:type: 4,5,  half of 25,26
!              ** plus fraction of urban... ***
        IF(high(i).ge.EPSILON_SVS)THEN
!
           if (lat(i).ge.0.0) then    
!              northern hemisphere
              deci(i) = (fcover(i,6) * tablen(6)  &
                      +  fcover(i,7) * tablen(7)  & 
                      +  fcover(i,8) * tablen(8)  & 
                      +  fcover(i,9) * tablen(9)  & 
                      +  (1. - furb_vl) * fcover(i,21) * tablen(21) & 
                      +  0.5 * fcover(i,25) * tablen(25) & 
                      +  0.5 * fcover(i,26) * tablen(26) ) /high(i)
!        
              ever(i) = 1.0 - deci(i)
           else
!              southern hemisphere
              deci(i) = (fcover(i,6) * tables(6)  &
                      +  fcover(i,7) * tables(7)  & 
                      +  fcover(i,8) * tables(8)  & 
                      +  fcover(i,9) * tables(9)  & 
                      +  (1. - furb_vl) * fcover(i,21) * tables(21) & 
                      +  0.5 * fcover(i,25) * tables(25) & 
                      +  0.5 * fcover(i,26) * tables(26) ) /high(i)
           

              ever(i) = 1.0 - deci(i)

           endif
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
 end subroutine veglowhigh
