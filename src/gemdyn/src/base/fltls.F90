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

!**s/r get_topo - Extract large scale of F_topo

      subroutine fltls (F_topo,F_topo_ls, Minx,Maxx,Miny,Maxy, F_rc,F_np)
      use tdpack
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(IN) :: Minx, Maxx, Miny, Maxy, F_np
      real   , intent(IN) :: F_rc
      real, dimension(Minx:Maxx,Miny:Maxy), intent(IN )::F_topo
      real, dimension(Minx:Maxx,Miny:Maxy), intent(OUT)::F_topo_ls
      
      integer :: i,j,n,im,ip,jm,jp
      real :: c1,c2,aux
      real, dimension(Minx:Maxx,Miny:Maxy) :: h1
!
!-----------------------------------------------------------------------
!
      c1 = 2.0/F_rc
      do j=Miny,Maxy
         do i=Minx,Maxx
            h1(i,j) = c1*F_topo(i,j)
            aux     = c1
            do n=1,(F_np-1)
               c2 = (2./F_rc)*((sin(2.*pi_8*n/F_rc))/(2.*pi_8*n/F_rc))*&
                              ((sin(2.*pi_8*n/F_np))/(2.*pi_8*n/F_np))
               im = i-n
               ip = i+n
               if ( im.ge.Minx .and. ip.le.Maxx  ) then
                  h1(i,j) = h1(i,j) + c2*( F_topo(im,j) + F_topo(ip,j) )
                  aux = aux + 2.*c2
               endif
            enddo
            h1(i,j) = h1(i,j)/aux      
         enddo
      enddo

      do j=Miny,Maxy
         do i=Minx,Maxx
            F_topo_ls(i,j) = c1*h1(i,j)
            aux     = c1
            do n=1,(F_np-1)
               c2 = (2./F_rc)*((sin(2.*pi_8*n/F_rc))/(2.*pi_8*n/F_rc))*&
                              ((sin(2.*pi_8*n/F_np))/(2.*pi_8*n/F_np))
               jm = j-n
               jp = j+n
               if ( jm.ge.Miny .and. jp.le.Maxy ) then
                  F_topo_ls(i,j) = F_topo_ls(i,j) + c2*( h1(i,jm) + h1(i,jp) )
                  aux = aux + 2.*c2
               endif
            enddo
            F_topo_ls(i,j) = F_topo_ls(i,j)/aux
         end do
      end do
!
!-----------------------------------------------------------------------
!
      return
      end subroutine fltls
