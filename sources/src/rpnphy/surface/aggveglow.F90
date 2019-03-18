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

subroutine aggveglow(fcover, tablen, tables, aggf, lat, ni, nclass)
  use svs_configs, only : ntypel, vl_type, furb_vl, epsilon_svs
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
   ! AGGF       Aggretated geophysical field representative of an entire
   !            grid area, only considerinf high vegetation classes i.e, 10-20,22,23

   integer ni, nclass
   real aggf(ni), fcover(ni,nclass), tablen(nclass)
   real tables(nclass), lat(ni)
!
!
!Author
!        Maria Abrahamowicz (dec 2012)
!Revision
! 001    Maria Abrahamowicz (sept. 2013)
!        remove nclassurb argument and variable
!Object
!        Aggregation of vegetation fields (veg fraction, lai, rsmin, etc...)
!        from the vegetation fraction mask for LOW VEGATION TYPE
!
! 
      integer i, m, type

      REAL totfract, table_val
      
!
      DO i=1,ni
        aggf(i) = 0.0
      END DO
!
!
!
      DO i=1,ni
        totfract = 0.
        DO m=1,ntypel
           type=vl_type(m)  ! loop on low vegetation classes
!
           totfract = totfract + fcover(i,type)
!
        END DO
        !  ADD URBAN class (21)
        type=21
        totfract = totfract + furb_vl  * fcover(i,type)   
!
        IF (totfract.GE.EPSILON_SVS) THEN

          DO m=1,ntypel  
           type=vl_type(m)  ! loop on low vegetation classes
!
            if (lat(i).ge.0.0) then    
!              northern hemisphere
               table_val = tablen(type)
            else
!              southern hemisphere
               table_val = tables(type)
            endif
            aggf(i) = aggf(i) + fcover(i,type)* table_val
          END DO

          !  ADD URBAN class (21)
          type=21
          if (lat(i).ge.0.0) then    
!              northern hemisphere
             table_val = tablen(type)
          else
!              southern hemisphere
             table_val = tables(type)
          endif
          aggf(i) = aggf(i) + furb_vl * fcover(i,type) * table_val

 ! -----------       FINAL AGGREGATED VALUE

          aggf(i) = aggf(i) / totfract
       ELSE
          aggf(i) = 0.0

       END IF
!
      END DO
!

   return
 end subroutine aggveglow
