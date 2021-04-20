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

!***s/r nest_init_weight:Initialize a 3D matrix of weight (values between 0 and 1) that
!                        will be used to compute the weight of the piloting field. 1 means
!                        that the piloting field is taken entirely 0 means that piloting
!                        field has no influence (center of domain) [0.0-1.0] means that
!                        there is a blending between piloting field and original field
!
subroutine nest_init_weight (F_weight,F_si,F_sj,F_sk,Minx,Maxx,Miny,Maxy,Nk)
      use gem_options
      use lam_options
      use grdc_options
      use glb_ld
      use glb_pil
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

  ! Arguments
      integer, intent(IN) :: F_si,F_sj,F_sk,Minx,Maxx,Miny,Maxy,Nk
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(out) :: F_weight

  ! Local variables
  integer w1,w2,e1,e2,s1,s2,n1,n2,t1,t2
  integer l_wdeb,l_wfin,l_edeb,l_efin
  integer l_sdeb,l_sfin,l_ndeb,l_nfin,l_tdeb,l_tfin
  integer i, j, k ,ig, jg, nit, njt, shifti, shiftj, shiftk
  real(kind=REAL64) pis2,lx,ly,lz,xx,yy
  real(kind=REAL64), parameter :: zero=0.d0, one=1.d0
  real(kind=REAL64) weight_x(l_ni), weight_y(l_nj), weight_z(G_nk+1)
  real wk1(G_ni,G_nj)

  !
  !----------------------------------------------------------------------
  !
  shifti = F_si; shiftj = F_sj; shiftk = F_sk

  pis2 = acos(zero)

  F_weight = 0.

  lx = dble(Lam_blend_Hx+1)
  ly = dble(Lam_blend_Hy+1)
  lz = dble(Lam_blend_T +1)

  weight_x = 0. ; weight_y = 0. ; weight_z = 0.

  w1 = Glb_pil_w
  w2 = w1 + Lam_blend_Hx

  e1 = G_ni - Glb_pil_e + 1 + shifti
  e2 = e1   - Lam_blend_Hx

  s1 = Glb_pil_s
  s2 = s1 + Lam_blend_Hy

  n1 = G_nj - Glb_pil_n + 1 + shiftj
  n2 = n1   - Lam_blend_Hy

  t1 = Lam_gbpil_T + shiftk
  t2 = t1 + Lam_blend_T

  l_wdeb = max(1   , 2   -Ptopo_gindx(1,Ptopo_myproc+1  ))
  l_wfin = min(l_ni,w2   -Ptopo_gindx(1,Ptopo_myproc+1)+1)
  l_edeb = max(1   ,e2   -Ptopo_gindx(1,Ptopo_myproc+1)+1)
  l_efin = min(l_ni,G_ni -Ptopo_gindx(1,Ptopo_myproc+1)+1)

  l_sdeb = max(1   , 2   -Ptopo_gindx(3,Ptopo_myproc+1  ))
  l_sfin = min(l_nj,s2   -Ptopo_gindx(3,Ptopo_myproc+1)+1)
  l_ndeb = max(1   ,n2   -Ptopo_gindx(3,Ptopo_myproc+1)+1)
  l_nfin = min(l_nj,G_nj -Ptopo_gindx(3,Ptopo_myproc+1)+1)

  l_tdeb = 1
  l_tfin = min(G_nk+1,t2)

  if (l_wfin > l_edeb) stop 'ABORT nest_init_weight l_wfin > l_edeb'
  if (l_sfin > l_ndeb) stop 'ABORT nest_init_weight l_sfin > l_ndeb'

  do i = l_wdeb, l_wfin
     ig = i + Ptopo_gindx(1,Ptopo_myproc+1) - 1
     j  = max(w1,ig)
     weight_x(i) = (cos(pis2*(j-w1)/lx))**2
  end do
  do i = l_edeb, l_efin
     ig = i + Ptopo_gindx(1,Ptopo_myproc+1) - 1
     j  = min(e1,ig)
     weight_x(i) = max(weight_x(i),(cos(pis2*(e1-j)/lx))**2)
  end do
  do i = l_sdeb, l_sfin
     ig = i + Ptopo_gindx(3,Ptopo_myproc+1) - 1
     j  = max(s1,ig)
     weight_y(i) = (cos(pis2*(j-s1)/ly))**2
  end do
  do i = l_ndeb, l_nfin
     ig = i + Ptopo_gindx(3,Ptopo_myproc+1) - 1
     j  = min(n1,ig)
     weight_y(i) = max(weight_y(i),(cos(pis2*(n1-j)/ly))**2)
  end do

  do j=1,l_nj
  do i=1,l_ni
     F_weight(i,j,1) = min(1.0d0, weight_x(i) + weight_y(j))
  end do
  end do

  nit  = G_ni-Glb_pil_e
  njt  = G_nj-Glb_pil_n

  !south-west
  do j = l_sdeb, l_sfin
  do i = l_wdeb, l_wfin
     ig = i + Ptopo_gindx(1,Ptopo_myproc+1) - 1
     jg = j + Ptopo_gindx(3,Ptopo_myproc+1) - 1
     xx = (dble(w2+1-ig)/lx)**2
     yy = (dble(s2+1-jg)/ly)**2
     F_weight(i,j,1) = (cos(pis2*(one-min(one, sqrt(xx+yy) ))))**2
  end do
  end do

  !south-east
  do j = l_sdeb, l_sfin
  do i = l_edeb, l_efin
     ig = i + Ptopo_gindx(1,Ptopo_myproc+1) - 1
     jg = j + Ptopo_gindx(3,Ptopo_myproc+1) - 1
     xx = (dble(ig-e2+1)/lx)**2
     yy = (dble(s2+1-jg)/ly)**2
     F_weight(i,j,1) = (cos(pis2*(one-min(one, sqrt(xx+yy) ))))**2
     end do
  end do

  !north-west
  do j = l_ndeb, l_nfin
  do i = l_wdeb, l_wfin
     ig = i + Ptopo_gindx(1,Ptopo_myproc+1) - 1
     jg = j + Ptopo_gindx(3,Ptopo_myproc+1) - 1
     xx = (dble(w2+1-ig)/lx)**2
     yy = (dble(jg-n2+1)/ly)**2
     F_weight(i,j,1) = (cos(pis2*(one-min(one, sqrt(xx+yy) ))))**2
     end do
  end do

  !north-east
  do j = l_ndeb, l_nfin
  do i = l_edeb, l_efin
     ig = i + Ptopo_gindx(1,Ptopo_myproc+1) - 1
     jg = j + Ptopo_gindx(3,Ptopo_myproc+1) - 1
     xx = (dble(ig-e2+1)/lx)**2
     yy = (dble(jg-n2+1)/ly)**2
     F_weight(i,j,1) = (cos(pis2*(one-min(one, sqrt(xx+yy) ))))**2
     end do
  end do

  !Top
  do k = l_tdeb, l_tfin
     j = max(t1,k)
     weight_z(k) = (cos(pis2*(j-t1)/lz))**2
  end do

!!!!TODO: remove this section when BCs treatment re-formulated
  call glbcolc (wk1,G_ni,G_nj,F_weight,l_minx,l_maxx,l_miny,l_maxy,1)
  if (Ptopo_myproc == 0) then
     do i = w2, e2
        wk1(i,s2) = 0.
        wk1(i,n2) = 0.
     end do
     do j = s2, n2
        wk1(w2,j) = 0.
        wk1(e2,j) = 0.
     end do
  endif
  call glbdist (wk1,G_ni,G_nj,F_weight,l_minx,l_maxx,l_miny,l_maxy,1,G_halox,G_haloy)
!!!!
      if (l_west ) F_weight(  1-G_halox:0,:,1) = 1.
      if (l_south) F_weight(:,1-G_haloy:0  ,1) = 1.
      if (l_east ) F_weight(  l_ni+1:l_ni+G_halox,:,1) = 1.
      if (l_north) F_weight(:,l_nj+1:l_nj+G_haloy  ,1) = 1.

      do k= G_nk+1, 1, -1
         do j= 1-G_haloy, l_nj+G_haloy
            do i= 1-G_halox, l_ni+G_halox
               F_weight(i,j,k) = max (F_weight(i,j,1) ,real(weight_z(k)))
            end do
         end do
      end do
  !
  !----------------------------------------------------------------------
  !

end subroutine nest_init_weight
