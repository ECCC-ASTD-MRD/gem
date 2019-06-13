!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

!**s/p lipschitz - compute lipschitz stability criterion

      subroutine lipschitz(F_u, F_v, F_w, Minx, Maxx, Miny, Maxy, Nk, &
                           F_i0, F_in, F_j0, F_jn)

      use geomh
      use glb_ld
      use cstv
      use ver
      use type_mod
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      integer, intent(in) :: F_i0, F_in, F_j0, F_jn
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(in) :: F_w, F_u, F_v

!author
!     Claude Girard & Andre Plante - oct 2006
!
!object
!	compute lipschitz number
!       See Pudykiewicz et al. Atmos.0cean 1985 267-303 Appendix A
!
!***************************************************************************
!
!arguments
!  Name                        Description
!------------------------------------------------------------
! F_w           - model s vertical velocity (pi* coordinates)
! F_u           - model zonal horizontal wind component
! F_v           - model meridional horizontal wind component
! F_i0          - starting point of calculation on W-E axis
! F_in          - ending point of calculation on W-E axis
! F_j0          - starting point of calculation on N-S axis
! F_jn          - ending point of calculation on N-S axis


      integer :: i,j,k,i0,j0,km, iu,iv,iw, ju,jv,jw, ku,kv,kw
      real :: dudx, dvdy, dwdz, LipNOu, LipNOv, LipNOw
!     __________________________________________________________________
!
      i0=F_i0
      j0=F_j0
      if (l_west ) i0 = max(2,F_i0)
      if (l_south) j0 = max(2,F_j0)

      iu=-99
      ju=-99
      ku=-99
      iv=-99
      jv=-99
      kv=-99
      iw=-99
      jw=-99
      kw=-99
      LipNOu=0.0
      LipNOv=0.0
      LipNOw=0.0

!$omp do
      do k=1,l_nk
         km=max(k-1,1)
         do j= j0, F_jn
         do i= i0, F_in
            dudx=abs(F_u(i,j,k)-F_u(i-1,j,k))*geomh_invcy_8(j)*geomh_invDX_8(j)
            if(dudx > LipNOu) then
               iu=i
               ju=j
               ku=k
               LipNOu=dudx
            end if
            dvdy=abs(F_v(i,j,k)-F_v(i,j-1,k))*geomh_invcy_8(j)*geomh_invDY_8
            if(dvdy > LipNOv) then
               iv=i
               jv=j
               kv=k
               LipNOv=dvdy
            end if
            dwdz=abs(F_w(i,j,k)-F_w(i,j,km))*ver_idz_8%m(km)
            if(dwdz > LipNOw) then
               iw=i
               jw=j
               kw=k
               LipNOw=dwdz
            end if
         end do
         end do
      end do
!$omp enddo

      LipNOu=Cstv_dt_8*LipNOu
      LipNOv=Cstv_dt_8*LipNOv
      LipNOw=Cstv_dt_8*LipNOw

      write(6,101) LipNOu,iu,ju,ku
      write(6,102) LipNOv,iv,jv,kv
      write(6,103) LipNOw,iw,jw,kw
  101 format(30x,'Lipschitz number(u)=',F6.2,' [',i4,',',i4,',',i3']')
  102 format(30x,'Lipschitz number(v)=',F6.2,' [',i4,',',i4,',',i3']')
  103 format(30x,'Lipschitz number(w)=',F6.2,' [',i4,',',i4,',',i3']')
!     __________________________________________________________________
!
      return
      end
