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

!/@*
subroutine printbus(F_bus_S)
   implicit none
   character(len=*),intent(in) :: F_bus_S
   !*@/
!!!#include <arch_specific.hf>
#include "buses.cdk"
   integer :: i
   ! ---------------------------------------------------------------------
   select case(F_bus_S(1:1))
   case('E')
      write(6,101) 'ENT'
      write(6,110)
      write(6,130)
      do i=1,enttop
         write (6,120) entnm(i,1),entnm(i,2),entdc(i),entpar(i,1),entpar(i,2), &
              entpar(i,6),entpar(i,8),entpar(i,3),entpar(i,4)
      end do
   case('D')
      write(6,101) 'DYN'
      write(6,110)
      write(6,130)
      do i=1,dyntop
         write (6,120) dynnm(i,1),dynnm(i,2),dyndc(i),dynpar(i,1),dynpar(i,2), &
              dynpar(i,6),dynpar(i,8),dynpar(i,3),dynpar(i,4)
      end do
   case('P')
      write(6,101) 'PER'
      write(6,110)
      write(6,130)
      do i=1,pertop
         write (6,120) pernm(i,1),pernm(i,2),perdc(i),perpar(i,1),perpar(i,2), &
              perpar(i,6),perpar(i,8),perpar(i,3),perpar(i,4)
      end do
   case('V')
      write(6,101) 'VOL'
      write(6,110)
      write(6,130)
      do i=1,voltop
         write (6,120) volnm(i,1),volnm(i,2),voldc(i),volpar(i,1),volpar(i,2), &
              volpar(i,6),volpar(i,8),volpar(i,3),volpar(i,4)
      end do
   case default
      return   
   end select
   write(6,130)
   print *,' '
   ! ---------------------------------------------------------------------
101 format (/48x,'+',14('-'),'+'/48x,'|  **',a3,'BUS**  |'/ &
        '+',17('-'),'+',12('-'),'+',16('-'),'+',14('-'),'+',9('-'), &
        '+',8('-'),'+',8('-'),'+',5('-'),'+',5('-'),'+',3('-'),'+',4('-'),'+')
110 format ('|',6x,'Name',7x,'|',3x,'Output',3x,'|', &
        11x,'Description ',18x,'|  Start | Length | Mul | Mos |Ini|Stag|'/ &
        '|',17x,'|',4x,'name',4x,'|',41x,'|',8x,'|',8x,'|',5x,'|', &
        ' aic |',3x, '|ger |')
120 format ('|',1x,a16,'|',3x,'"',a4,'"',3x,'|',1x, &
        a40,'|',2(i7,' |'),2(i4,' |'),i2,' |',i2,'  |')
130 format ('+',17('-'),'+',12('-'),'+',41('-'),'+',8('-'),'+',8('-'), &
        '+',5('-'),'+',5('-'),'+',3('-'),'+',4('-'),'+')
   return
end subroutine printbus
