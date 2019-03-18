!--------------------------------- LICENCE BEGIN -------------------------------
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

!**s/r wil_set - Check and Convert parameters in williamson namelist

      subroutine wil_set (F_topo_L,F_adv_L,F_unout,F_err)
      use dcst
      use wil_options
      use tdpack, only : rayt_8, omega_8, pi_8

      implicit none

      logical F_topo_L,F_adv_L
      integer F_unout,F_err

      !object
      !======================================================
      !   Check and Convert parameters in williamson namelist
      !======================================================

      F_err = 0

      Dcst_rayt_8     = rayt_8
      Dcst_inv_rayt_8 = 1.d0 / rayt_8
      Dcst_omega_8    = omega_8

      if (Williamson_case/=5) F_topo_L = .false.

      if (Williamson_case/=1.and.Williamson_Terminator_L) then
         F_err = -1
         if (F_unout>0) write (F_unout, 6200)
         return
      endif

      if (Williamson_case==3.or.Williamson_case==4) then
         F_err = -1
         if (F_unout>0) write (F_unout, 6300)
         return
      endif

      if (Williamson_case>2.and.Williamson_alpha/=0.) then
         F_err = -1
         if (F_unout>0) write (F_unout, 6400)
        return
      endif

      !Conversion Degree to Radian
      !---------------------------
      Williamson_alpha = Williamson_alpha*(pi_8/180.0)
      Williamson_rlon0 = Williamson_rlon0*(pi_8/180.0)
      Williamson_rlat0 = Williamson_rlat0*(pi_8/180.0)
      Williamson_clon0 = Williamson_clon0*(pi_8/180.0)
      Williamson_clat0 = Williamson_clat0*(pi_8/180.0)

      !Identify 2D Advection runs
      !--------------------------
      F_adv_L = Williamson_case==1
      if (F_unout>0.and.F_adv_L) write (F_unout, 7000)

      !---------------------------------------------------------------

      return

 6200 format (/' Williamson_Terminator only available for Williamson case 1'/)
 6300 format (/' Williamson case 3 and 4 are not available'/)
 6400 format (/' Williamson Alpha(in deg) must be 0.0 for cases greater than 2 '/)
 7000 format (//'  ====================='/&
                '  ACADEMIC 2D Advection'/&
                '  ====================='//)

      end subroutine wil_set
