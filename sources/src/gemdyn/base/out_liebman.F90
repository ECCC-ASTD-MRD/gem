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

      subroutine out_liebman ( ttx, htx, vt, gz, fis0, wlao,&
                               Minx,Maxx,Miny,Maxy,nkund,nk )
      use out_options
      use tdpack
      use glb_ld
      use out3
      implicit none
#include <arch_specific.hf>
!
      integer Minx,Maxx,Miny,Maxy, nkund, nk
      real ttx (Minx:Maxx,Miny:Maxy,nkund), &
           htx (Minx:Maxx,Miny:Maxy,nkund), &
           vt  (Minx:Maxx,Miny:Maxy,nk   ), &
           gz  (Minx:Maxx,Miny:Maxy,nk   ), &
           fis0(Minx:Maxx,Miny:Maxy), wlao (Minx:Maxx,Miny:Maxy)

      integer i,j,k,kk,kgrnd
      real w2(l_minx:l_maxx,l_miny:l_maxy,Out3_lieb_nk), grad
!
!----------------------------------------------------------------------
!
!$omp parallel private (grad, kgrnd) shared (w2)
!$omp do

      do k=1,Out3_lieb_nk
         do j=1,l_nj
            do i=1,l_ni

               ! Store fictitious height level in htx
               htx(i,j,k) = Out3_lieb_levels(k) * grav_8

               ! Determine if fictitious level is above or below ground
               ttx(i,j,k) = fis0(i,j) - htx(i,j,k)

               if ( ttx(i,j,k) > 0 ) then

                  ! fictitious level is under ground:
                  ! temperature is obtained by linear EXTrapolation
                  ! identify under ground grid point

                  if ( abs( wlao(i,j)*180./pi_8 ) >= 49. ) then
                     ttx(i,j,k) = vt(i,j,Nk) +       .0005 * ttx(i,j,k)
                  else
                     ttx(i,j,k) = vt(i,j,Nk) + stlo_8 * ttx(i,j,k)
                  endif
                  w2(i,j,k) = 1.0

               else

                  ! fictitious level is above ground:
                  ! temperature is obtained by linear INTerpolation
                  ! identify above ground grid point

                  do kk= nk, 1, -1
                     kgrnd = kk
                     ttx(i,j,k) = gz (i,j,kk) - htx(i,j,k)
                     if ( ttx(i,j,k) > 0. ) exit
                  enddo

                  if (kgrnd == nk) then
                     ! This happens only is physics is off.
                     ! When physics is on, gz(nk) = fis0(i,j), so the lowest
                     ! possible level for which the above test is true will be kgrnd=nk-1,
                     ! for this case kgrnd+1 is fine.
                     ! Whithout physics, gz(nk) is just above the surface,
                     ! therefore, at times, it will be the level chosen with the test above,
                     ! kgrnd=nk. Without the actual testm this will result in a out of bound
                     ! index due to kgrnd+1.
                     ! To solve this we set the gradient to zero.
                     grad = 0.
                  else
                     grad = - (vt(i,j,kgrnd) - vt(i,j,kgrnd+1) ) / &
                     (gz(i,j,kgrnd) - gz(i,j,kgrnd+1) )
                  end if
                  ttx(i,j,k) = vt (i,j,kgrnd) + grad * ttx(i,j,k)
                  w2(i,j,k) = 0.0
               end if

            end do
         end do
      end do
!$omp enddo
!$omp end parallel

      call out_padbuf (ttx,l_minx,l_maxx,l_miny,l_maxy,Out3_lieb_nk)

      call liebman_dm (ttx,w2,Out3_lieb_conv,Out3_lieb_maxite,&
                       l_minx,l_maxx,l_miny,l_maxy,Out3_lieb_nk)
!
!----------------------------------------------------------------------
!
      return
      end
