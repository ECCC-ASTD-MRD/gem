
!/@*
subroutine printbus(F_bus_S)
   use clib_itf_mod, only: clib_toupper
   use phymem, only: phymeta, nphyvars, phymem_busidx, phymem_getmeta
   implicit none
   character(len=*),intent(in) :: F_bus_S
   !*@/
!!!#include <arch_specific.hf>
   integer :: iv, ib, istat
   character(len=64) :: vname, oname, desc
   type(phymeta), pointer :: vmeta
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
   do iv = 1, nphyvars
      nullify(vmeta)
      istat = phymem_getmeta(vmeta, iv)
      if (istat < 0) cycle
      if (vmeta%ibus /= ib) cycle
      vname = vmeta%vname
      oname = vmeta%oname
      desc  = vmeta%desc
      istat = clib_toupper(vname)
      istat = clib_toupper(oname)
      istat = clib_toupper(desc)
      write (6,120) &
           vname, oname, desc, &
           vmeta%i0, &
           vmeta%size, &
           vmeta%fmul, &
           vmeta%mosaic, &
           vmeta%init, &
           vmeta%stag
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
