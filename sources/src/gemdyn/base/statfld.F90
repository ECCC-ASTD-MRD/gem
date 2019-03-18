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


!**s/r statfld2 - APRES glbcolc, calcule la moyenne, la variance, minimum et
!                 le maximum d un champs et imprime le resultat.
!

      subroutine statfld3 (F_field, F_nv_S, F_no, F_from_S, &
                           minx,maxx,miny,maxy,mink,maxk, &
                           F_i0,F_j0,F_k0,F_in,F_jn,F_kn,F_rx)
!
      implicit none
#include <arch_specific.hf>
!
      character(len=*), intent(in) :: F_nv_S , F_from_S
      integer, intent(in) :: minx,maxx,miny,maxy,mink,maxk, &
              F_i0,F_j0,F_k0,F_in,F_jn,F_kn,F_no,F_rx
      real, dimension(minx:maxx,miny:maxy,mink:maxk), intent(in) :: F_field
!
!author
!     M. Desgagne
!
!object
!     calcule et imprime: la moyenne    (moy)
!                         la variance   (var)
!                         le minimum et le maximum du champ f
!
!arguments
!  Name                            Description
!----------------------------------------------------------------
! F_field                 Field to be operated on
! F_nv_S                  User provided string to define F_field
! F_no                    Usually the timestep #
! F_from_S                Usually the name of the calling subroutine
! F_i0,F_j0               Global lower-left indexes of the sub-domain
!                            on which to perform statistics
! F_in,F_jn               Global upper-right indexes of the sub-domain
!                            on which to perform statistics
! F_k0,F_kn               Range of levels on which to perform statistics
!----------------------------------------------------------------
!


      integer :: i,j,k,imin,jmin,kmin,imax,jmax,kmax
      real*8 :: sum,sumd2,moy,var,mind,maxd,fijk,npt_8
      integer :: no
!
!--------------------------------------------------------------------
!
      npt_8 = 1.0d0*((F_in-F_i0+1)*(F_jn-F_j0+1)*(F_kn-F_k0+1))
!
      sum   = 0.d0
      sumd2 = 0.d0
      imin  = F_i0
      jmin  = F_j0
      kmin  = F_k0
      imax  = F_in
      jmax  = F_jn
      kmax  = F_kn
      maxd  = F_field(F_in,F_jn,F_kn)
      mind  = F_field(F_i0,F_j0,F_k0)
!
      do k=F_k0,F_kn
      do j=F_j0,F_jn
      do i=F_i0,F_in
         fijk = F_field(i,j,k)
         sum = sum + fijk
         sumd2 = sumd2 + fijk*fijk
         if (fijk > maxd) then
            maxd = fijk
            imax = i
            jmax = j
            kmax = k
         end if
         if (fijk < mind) then
            mind = fijk
            imin = i
            jmin = j
            kmin = k
         end if
      end do
      end do
      end do
!
      moy = sum / npt_8
      var = max(0.d0,1.0d0*(sumd2 + moy*moy*npt_8 - 2*moy*sum) / npt_8)
      var = sqrt(var)
      no  = F_no

      imin = imin
      imax = imax
      jmin = jmin
      jmax = jmax
!
! ** On imprime
!
         if (F_rx < 8) then
            write(6,98) no,F_nv_S,moy,var,imin,jmin,kmin,mind, &
                                       imax,jmax,kmax,maxd,F_from_S
         else
            write(6,99) no,F_nv_S,moy,var,imin,jmin,kmin,mind, &
                                       imax,jmax,kmax,maxd,F_from_S
         end if
!      end if
!
 98   format (i4,a4,' Mean:',1pe14.6,' Std:',1pe14.6, &
              ' Min:[(',i4,',',i4,',',i4,')', &
              1pe14.6,']',' Max:[(',i4,',',i4,',',i4,')', &
              1pe14.6,']',a6)
 99   format (i4,a4,' Mean:',1pe22.12,' Std:',1pe22.12,/ &
              ' Min:[(',i4,',',i4,',',i4,')', &
              1pe22.12,']',' Max:[(',i4,',',i4,',',i4,')', &
              1pe22.12,']',a6)
!
!----------------------------------------------------------------
!
      return
      end
