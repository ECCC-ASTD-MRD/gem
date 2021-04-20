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

!**s/r wil_case6 - To setup Williamson Case 6: Rossby-Haurwitz wave (HEIGHT)

      subroutine wil_case6 (F_gz,F_minx,F_maxx,F_miny,F_maxy,F_nk)

      use gem_options
      use lun
      use glb_ld
      use ptopo
      use tdpack
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in)  :: F_minx,F_maxx,F_miny,F_maxy,F_nk
      real,    intent(out) :: F_gz(F_minx:F_maxx,F_miny:F_maxy,F_nk)

      !authors
      !     Abdessamad Qaddouri and Vivian Lee
      !
      !revision
      ! v5_00 - Tanguay M. - Clean Up
      !
      !object
      !==============================================================
      !     To setup Williamson Case 6: Rossby-Haurwitz wave (HEIGHT)
      !     Williamson et al.,1992,JCP,102,211-224
      !==============================================================

      integer :: i,j,k,R_case,g_i0,g_in,g_j0,g_jn,i0,in,j0,jn,zlist
      real(kind=REAL64) :: phi0_8,dlon_8,K_Case_8,OMG_8,                   &
              rlon_8,rlat_8,time_8, sint_8,cost_8,phiay_8,phiby_8,phicy_8, &
              s_8(2,2),x_a_8,y_a_8,                                        &
              phia_8(G_nj),phib_8(G_nj),phic_8(G_nj)
      real :: picll(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy), &
              gzloc(F_minx:F_maxx,F_miny:F_maxy)
!
!---------------------------------------------------------------------
!
      g_i0= 1-G_halox ; g_in= G_ni+G_halox
      g_j0= 1-G_haloy ; g_jn= G_nj+G_haloy

      i0= 1-G_halox ; in= l_ni+G_halox
      j0= 1-G_haloy ; jn= l_nj+G_haloy

      if (Lun_out>0) write(Lun_out,*) ''
      if (Lun_out>0) write(Lun_out,*) '--------------------------------------------'
      if (Lun_out>0) write(Lun_out,*) 'WILLIAMSON CASE6, Williamson et al. (1992)  '
      if (Lun_out>0) write(Lun_out,*) 'Rossby-Haurwitz wave                        '
      if (Lun_out>0) write(Lun_out,*) '--------------------------------------------'

      time_8   = 0.0
      R_Case   = 4
      K_Case_8 = 7.848E-6
      OMG_8    = 7.848E-6
      phi0_8   = 8000.0
      dlon_8   = (R_Case*(3+R_Case)*OMG_8 - 2.0*omega_8)/ &
                 ((1+R_Case)*(2+R_Case))*time_8

      !Compute tracer for YIN
      !----------------------
      if (Ptopo_couleur==0) then

         do j=g_j0,g_jn

            rlat_8 = G_yg_8(j)

            cost_8 = cos(rlat_8)

            !Compute latitude-dependent factors for height
            !---------------------------------------------
            phia_8(j) = 0.5*OMG_8*(2.0*omega_8+OMG_8)*cost_8*cost_8 +           &
                        0.25*K_Case_8*K_Case_8*cost_8**(2*R_Case) *             &
                        ((R_Case+1)*cost_8*cost_8+(2*R_Case*R_Case-R_Case-2) -  &
                        2.0*R_Case*R_Case/(cost_8*cost_8))
            phib_8(j) = (2.0*(omega_8+OMG_8)*K_Case_8)/((R_Case+1)*(R_Case+2))* &
                        cost_8**R_Case*                                         &
                        ((R_Case*R_Case+2*R_Case+2)-(R_Case+1)**2*cost_8*cost_8)
            phic_8(j) = 0.25*K_Case_8*K_Case_8*cost_8**(2*R_Case)* &
                        ((R_Case+1)*cost_8*cost_8-(R_Case+2))

            do i=g_i0,g_in

               rlon_8 = G_xg_8(i)

               sint_8 = sin(rlat_8)
               cost_8 = cos(rlat_8)

               picll(i,j) = phi0_8 + (rayt_8*rayt_8*(phia_8(j)+phib_8(j) &
                            * cos(R_Case*(rlon_8-dlon_8))+phic_8(j)      &
                            * cos(2*R_Case*(rlon_8-dlon_8))))/grav_8

            end do

         end do

      !Compute tracer for YAN
      !----------------------
      else

         do j=g_j0,g_jn

            do i=g_i0,g_in

               x_a_8 = G_xg_8(i)-acos(-1.D0)
               y_a_8 = G_yg_8(j)

               call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

               rlon_8 = rlon_8+acos(-1.D0)

               cost_8 = cos(rlat_8)

               !Compute latitude-dependent factors for height
               !---------------------------------------------
               phiay_8 = 0.5*OMG_8*(2.0*omega_8+OMG_8)*cost_8*cost_8 +          &
                         0.25*K_Case_8*K_Case_8*cost_8**(2*R_Case) *            &
                         ((R_Case+1)*cost_8*cost_8+(2*R_Case*R_Case-R_Case-2) - &
                         2.0*R_Case*R_Case/(cost_8*cost_8))
               phiby_8 = (2.0*(omega_8+OMG_8)*K_Case_8)/                        &
                         ((R_Case+1)*(R_Case+2))*cost_8**R_Case*                &
                         ((R_Case*R_Case+2*R_Case+2)-(R_Case+1)**2*cost_8*cost_8)
               phicy_8 = 0.25*K_Case_8*K_Case_8*cost_8**(2*R_Case)*             &
                         ((R_Case+1)*cost_8*cost_8-(R_Case+2))

               sint_8 = sin(rlat_8)

               picll(i,j) = phi0_8 + (rayt_8*rayt_8*(phiay_8+phiby_8 &
                            * cos(R_Case*(rlon_8-dlon_8))+phicy_8    &
                            * cos(2*R_Case*(rlon_8-dlon_8))))/grav_8
            end do

         end do

      end if

      zlist = 1

      call glbdist_os (picll,gzloc,&
                       F_minx,F_maxx,F_miny,F_maxy,1,&
                       G_ni+G_halox,G_nj+G_haloy,zlist,1,1.0d0,0.d0)

      do k=1,F_nk
         F_gz(i0:in,j0:jn,k) = gzloc(i0:in,j0:jn)
      end do
!
!---------------------------------------------------------------------
!
      return
      end
