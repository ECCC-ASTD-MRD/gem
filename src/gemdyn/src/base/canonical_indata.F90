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

!**s/r canonical_indata - Specific actions for imposing nesting data for
!                         canonical cases (Williamson/DCMIP)

      subroutine canonical_indata ()

      use dynkernel_options
      use wil_options
      use step_options
      use glb_ld
      use tr3d
      use mem_nest 

      use, intrinsic :: iso_fortran_env
      implicit none

      !arguments
      integer :: deb,n


      !Impose REFERENCE solution in nesting zone 
      !-----------------------------------------------------
      if (Schm_autobar_L.and.Williamson_case==1.and.(Williamson_Nair==0.or.Williamson_Nair==3)) then

         call wil_uvcase1 (nest_u_fin,nest_v_fin,l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.,Lctl_step)

         !Initialize Q1
         !-------------
         do n=1,Tr3d_ntr
            if (trim(Tr3d_name_S(n))/='Q1') cycle
            deb = (n-1) * l_nk + 1
            if (Williamson_Nair==0) call wil_case1(nest_tr_fin(l_minx,l_miny,deb),l_minx,l_maxx,l_miny,l_maxy,l_nk,0,Lctl_step)
            if (Williamson_Nair==3) call wil_case1(nest_tr_fin(l_minx,l_miny,deb),l_minx,l_maxx,l_miny,l_maxy,l_nk,5,Lctl_step)
         end do

      end if

      end
