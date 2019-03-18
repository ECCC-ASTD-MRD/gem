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

subroutine aggcovernat(fcover, tablen, tables, aggf, lat, ni, nclass)
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
   !            grid area

   integer ni, nclass
   real aggf(ni), fcover(ni,nclass), tablen(nclass)
   real tables(nclass), lat(ni)

   !@Author Stephane Belair
   !@Revision
   ! 001    B. Bilodeau and S. Belair - Adaptation for southern hemisphere
   ! 002    A. Lemonsu                - Add urban classes & change the name
   ! 003    B. Bilodeau (Jul 2007)    - Correct inconsistency caused by critmask
   ! 004    N. Gauthier (Jul 2011)    - Removed urban classes and change dimension of fcover !
   !
   !@Object
   !        Aggregation of vegetation fields (veg fraction, lai, rsmin, etc...)
   !        from the vegetation fraction mask.

   integer I,M
   real totfract, table_val

   do i=1,ni
      aggf(i) = 0.0
   end do

   do i=1,ni
      totfract = 0.
      do m=4,nclass                        ! loop on veg/natural soil classes
         totfract = totfract + fcover(i,m)
      end do
      if (totfract.gt.0.) then
         do m=4,nclass                      ! loop on veg/natural soil classes
            if (lat(i).ge.0.0) then
               !              northern hemisphere
               table_val = tablen(m)
            else
               !              southern hemisphere
               table_val = tables(m)
            endif
            aggf(i) = aggf(i) &
                 + fcover(i,m)* table_val
         end do
         aggf(i) = aggf(i) / totfract
      end if
   end do

   return
end subroutine aggcovernat
