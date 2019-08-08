
! * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
   subroutine ez_nwtncof(cx,cy,ax,ay,ni,nj,i1,i2,j1,j2,extension)
!*******
!     Auteur: Y. Chartier, drpn
!     Fevrier 1991
!
!Objet:  Calcul de coefficients uti1ises dans la forme newtonienne
!        de l'interpolation de Lagrange.
!
!*************************************************************
!   -----*-------------*------#------*----------*------->
!        x1            x2     x      x3         x4
!
!*************************************************************
!     cx(i,1) = 1.0 / (x2-x1)
!     cx(i,2) = 1.0 / (x3-x1)
!     cx(i,3) = 1.0 / (x3-x2)
!     cx(i,4) = 1.0 / (x4-x1)
!     cx(i,5) = 1.0 / (x4-x2)
!     cx(i,6) = 1.0 / (x4-x3)
!
!     structure identique pour cy(j,1..6)

!*******
   implicit none

   integer  :: ni,nj,i1,i2,j1,j2,extension
   real ::   x1,x2,x3,x4
   integer :: limite,imoins1,iplus1,iplus2,wrap
   real, dimension(:) :: ax(ni),ay(j1:j2)
   real, dimension(:,:) :: cx(ni,6),cy(j1:j2,6)
   logical sequence_ok

   integer i,j

   sequence_ok = .true.

   do i=1,ni-1
      if (ax(i+1) <= ax(i)) then
         sequence_ok = .false.
         exit
      endif
   enddo

   if (.not.sequence_ok) then
      print *, 'Probleme detecte dans EZ_NWTNCOF code 998'
      print *, '(EZ_NWTNCOF) Probleme : x1..x4  : ', ax(i), ax(i+1)
      print *, 'EZ_NWTNCOF CALL EXIT'
      call exit(13)
   endif

   do j=1,nj-1
      if (ay(j+1) <= ay(j)) then
         sequence_ok = .false.
         exit
      endif
   enddo

   if (.not.sequence_ok) then
      print *, 'Probleme detecte dans EZ_NWTNCOF code 999'
      print *, '(EZ_NWTNCOF) Probleme : y1..y4  : ', ay(j), ay(j+1)
      print *, 'EZ_NWTNCOF CALL EXIT'
      call exit(13)
   endif

   do i=1,ni
      do j=1,6
         cx(i,j) = 1.0
      enddo
   enddo

   do i=1,6
      do j=j1,j2
         cy(j,i) = 1.0
      enddo
   enddo

   do i=2,ni-2
      cx(i,1) = 1. / (ax(i  ) - ax(i-1))
      cx(i,2) = 1. / (ax(i+1) - ax(i-1))
      cx(i,3) = 1. / (ax(i+1) - ax(i  ))
      cx(i,4) = 1. / (ax(i+2) - ax(i-1))
      cx(i,5) = 1. / (ax(i+2) - ax(i  ))
      cx(i,6) = 1. / (ax(i+2) - ax(i+1))
   enddo

   do j=j1+1,j2-2
      cy(j,1) = 1. / (ay(j  ) - ay(j-1))
      cy(j,2) = 1. / (ay(j+1) - ay(j-1))
      cy(j,3) = 1. / (ay(j+1) - ay(j  ))
      cy(j,4) = 1. / (ay(j+2) - ay(j-1))
      cy(j,5) = 1. / (ay(j+2) - ay(j  ))
      cy(j,6) = 1. / (ay(j+2) - ay(j+1))
   enddo


   if (extension.eq.1) then
      x1 = ax(1) - (ax(ni) - ax(ni-1))
      x2 = ax(1)
      x3 = ax(2)
      x4 = ax(3)

      cx(1,1) = 1. / (x2-x1)
      cx(1,2) = 1. / (x3-x1)
      cx(1,3) = 1. / (x3-x2)
      cx(1,4) = 1. / (x4-x1)
      cx(1,5) = 1. / (x4-x2)
      cx(1,6) = 1. / (x4-x3)

      x1 = ax(ni-2)
      x2 = ax(ni-1)
      x3 = ax(ni)
      x4 = ax(ni) + (ax(2)-ax(1))

      cx(ni-1,1) = 1. / (x2-x1)
      cx(ni-1,2) = 1. / (x3-x1)
      cx(ni-1,3) = 1. / (x3-x2)
      cx(ni-1,4) = 1. / (x4-x1)
      cx(ni-1,5) = 1. / (x4-x2)
      cx(ni-1,6) = 1. / (x4-x3)
   endif

   if (extension.eq.2) then
      x1 = ax(1) - (360.0 - ax(ni))
      x2 = ax(1)
      x3 = ax(2)
      x4 = ax(3)

      cx(1,1) = 1. / (x2-x1)
      cx(1,2) = 1. / (x3-x1)
      cx(1,3) = 1. / (x3-x2)
      cx(1,4) = 1. / (x4-x1)
      cx(1,5) = 1. / (x4-x2)
      cx(1,6) = 1. / (x4-x3)

      x1 = ax(ni-2)
      x2 = ax(ni-1)
      x3 = ax(ni)
      x4 = ax(1) + 360.0

      cx(ni-1,1) = 1. / (x2-x1)
      cx(ni-1,2) = 1. / (x3-x1)
      cx(ni-1,3) = 1. / (x3-x2)
      cx(ni-1,4) = 1. / (x4-x1)
      cx(ni-1,5) = 1. / (x4-x2)
      cx(ni-1,6) = 1. / (x4-x3)

      x1 = ax(ni-1)
      x2 = ax(ni)
      x3 = ax(1)+360.0
      x4 = ax(2)+360.0

      cx(ni,1) = 1. / (x2-x1)
      cx(ni,2) = 1. / (x3-x1)
      cx(ni,3) = 1. / (x3-x2)
      cx(ni,4) = 1. / (x4-x1)
      cx(ni,5) = 1. / (x4-x2)
      cx(ni,6) = 1. / (x4-x3)
   endif
   return
   end subroutine ez_nwtncof

!********************************************************************
!**
!**

