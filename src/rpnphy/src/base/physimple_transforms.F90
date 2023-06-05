!-------------------------------------- LICENCE BEGIN -------------------------
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

module physimple_transforms

contains

   !/@*
   subroutine physimple_transforms3d(F_var_S,F_var_in_S,values)
      use, intrinsic :: iso_fortran_env, only: REAL64
      use clib_itf_mod, only: clib_toupper
      use tdpack_const
      implicit none
!!!#include <arch_specific.hf>
      real, pointer, dimension(:,:,:)  :: values
      character(len=*), intent(in) :: F_var_S, F_var_in_S
      !@Author Lubos Spacek - October 2009
      !@Object
      !     Checks the values after intrepolation, preform
      !     conversions of units and simple transformations
      !*@/
#include <rmn/msg.h>
#include <rmnlib_basics.hf>

      real(REAL64),parameter :: ONE_8   = 1.0d0
      real(REAL64),parameter :: CLXXX_8 = 180.0d0
      real(REAL64) :: deg2rad_8
      character(len=32) :: var_S, var_in_S, msg_S
      integer :: istat
      ! ---------------------------------------------------------------------
      var_S = trim(F_var_S)
      var_in_S = trim(F_var_in_S)
      istat = clib_toupper(var_S)
      istat = clib_toupper(var_in_S)

      deg2rad_8 = acos(-ONE_8)/CLXXX_8
      msg_S = ' '
      select case(var_S)
      case('DLAT')
         values = deg2rad_8*values
         msg_S = 'deg2rad'
      case('DLON')
         where(values>=0.)
            values = deg2rad_8*values
         elsewhere
            values = deg2rad_8*(values+360.)
         endwhere
         msg_S = 'deg2rad'
      case('GLACIER')
         where(values<0.) values = 0.0
         where(values>1.) values = 1.0
         msg_S = '[0, 1]'
      case('GLSEA0')
         where(values<0.) values = 0.0
         where(values>1.) values = 1.0
         msg_S = '[0, 1]'
      case('MG')
         where(values<0.) values = 0.0
         where(values>1.) values = 1.0
         msg_S = '[0, 1]'
      case('MT')
         where(values<0.) values = 0.0
         msg_S = '[0, '
      case('SNODP')
         values = 0.01*values
      case('SNODPL')
         values = 0.01*values
      case('SNVDP')
         values = 0.01*values
      case('TGLACIER')
         where(values<150.) values = values+tcdk
         msg_S = 'C2K conditional'
      case('TMICE')
         where(values<150.) values = values+tcdk
         msg_S = 'C2K conditional'
      case('TSOIL')
         where(values<150.) values = values+tcdk
         msg_S = 'C2K conditional'
      case('TWATER')
         where(values<150.) values = values+tcdk
         msg_S = 'C2K conditional'
      case('Z0EN')
         if (var_in_S == 'ZP') then
            values = exp(values)
            msg_S = 'Z0 = exp(ZP)'
         endif
      end select

      ! Physical world unit conversions
      if (msg_S == '') then
         select case(var_S(1:5))
         case('PW_TT')
            where(values<150.) values = values+tcdk
            msg_S = 'C2K conditional'
         case('PW_UU' , 'PW_VV')
            values = KNAMS*values
            msg_S = 'kts to m/s'
         end select
      endif

      if (msg_S /= '') call msg(MSG_INFO,'(physimple_transforms3d) '//trim(var_in_S)//' => '//trim(var_S)//' ('//trim(msg_S)//')')
      ! ---------------------------------------------------------------------
      return
   end subroutine physimple_transforms3d

end module physimple_transforms
