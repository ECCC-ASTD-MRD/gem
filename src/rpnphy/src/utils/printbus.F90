!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END ---------------------------

!/@*
subroutine printbus(F_bus_S)
   use clib_itf_mod, only: clib_toupper
   use phymem, only: pbuslist, phymem_busidx
   implicit none
   character(len=*),intent(in) :: F_bus_S
   !*@/
!!!#include <arch_specific.hf>
   integer :: iv, ib, istat
   character(len=64) :: vname, oname, desc
   ! ---------------------------------------------------------------------
   ib = phymem_busidx(F_bus_S)

   select case(F_bus_S(1:1))
   case('E')
      write(6,101) 'ENT'
   case('D')
      write(6,101) 'DYN'
   case('P')
      write(6,101) 'PER'
   case('V')
      write(6,101) 'VOL'
   case default
      return
   end select
   write(6,110)
   write(6,130)
   do iv = 1, pbuslist(ib)%nvars
      vname = pbuslist(ib)%meta(iv)%vname
      oname = pbuslist(ib)%meta(iv)%oname
      desc  = pbuslist(ib)%meta(iv)%desc
      istat = clib_toupper(vname)
      istat = clib_toupper(oname)
      istat = clib_toupper(desc)
      write (6,120) &
           vname, oname, desc, &
           pbuslist(ib)%meta(iv)%i0, &
           pbuslist(ib)%meta(iv)%size, &
           pbuslist(ib)%meta(iv)%fmul, &
           pbuslist(ib)%meta(iv)%mosaic, &
           pbuslist(ib)%meta(iv)%init, &
           pbuslist(ib)%meta(iv)%stag
   end do

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
