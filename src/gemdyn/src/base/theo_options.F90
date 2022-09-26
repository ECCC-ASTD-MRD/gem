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
module theo_options
      use bubble_options
      use mtn_options
      use ctrl
      use, intrinsic :: iso_fortran_env
      implicit none
      public
      save

   !# Choices of theoretical case:
   !# * from Robert, JAS 1993
   !# * 'BUBBLE': uniform convective bubble
   !# * from Schar et al, MWR 2002
   !# * 'MTN_SCHAR': case with N=0.01
   !# * 'MTN_SCHAR2': case with N=0.01871
   !# * from Pinty et al, MWR 1995
   !# * 'MTN_PINTY': linear isothermal with U=32m/s
   !# * 'MTN_PINTY2': linear isothermal with U=8m/s
   !# * 'MTN_PINTYNL': nonlinear regime

      character(len=15) :: Theo_case_S = 'NONE'
      namelist /theo_cfgs/ Theo_case_S

contains

      integer function theocases_nml (F_unf)
      use, intrinsic :: iso_fortran_env
      implicit none

      integer F_unf

!#include <rmnlib_basics.hf>
      logical nml_must
      character(len=64) :: nml_S
      character(len=15) :: dumc_S
!
!-------------------------------------------------------------------
!
! boiler plate - start
      if ( F_unf < 0 ) then
         theocases_nml= 0
         if ( Lun_out >= 0) write (Lun_out,nml=theo_cfgs)
         return
      end if

      theocases_nml= -1 ; nml_must= .false. ; nml_S= 'theo_cfgs'

      rewind(F_unf)
      read (F_unf, nml=theo_cfgs, end= 1001, err=1003)
      theocases_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         theocases_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (theocases_nml < 0 ) return
      if ((Lun_out>=0).and.(theocases_nml==0)) write (Lun_out, 6004) trim(nml_S)
      theocases_nml=0

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

      call low2up (Theo_case_S ,dumc_S)
      Theo_case_S = dumc_S

      if (  Theo_case_S == 'NONE' ) then
         Ctrl_theoc_L= .false.
      else
         Ctrl_theoc_L= .true.
         if (   Theo_case_S == 'MTN_SCHAR'   &
           .or. Theo_case_S == 'MTN_SCHAR2'  &
           .or. Theo_case_S == 'MTN_PINTY'   &
           .or. Theo_case_S == 'MTN_PINTY2'  &
           .or. Theo_case_S == 'MTN_PINTYNL' &
           .or. Theo_case_S == 'NOFLOW' ) then
           theocases_nml= mtn_nml (F_unf)
        else if (  Theo_case_S == 'BUBBLE' ) then
           theocases_nml= bubble_nml (F_unf)
        else
           if (Lun_out > 0) then
              write (Lun_out, 9200) Theo_case_S
              write (Lun_out, 8000)
              theocases_nml= -1
           end if
        end if
      end if

 8000 format (/,'========= ABORT IN S/R theocases_nml ============='/)
 9200 format (/,' Unsupported theoretical case: ',a/)

      return
      end function theocases_nml
!
!-------------------------------------------------------------------
!
      subroutine theo_cfg()
      use, intrinsic :: iso_fortran_env
      implicit none

      integer err
      integer, external :: bubble_cfg
!
!     ---------------------------------------------------------------
!
      if (  Theo_case_S == 'NONE' ) return

      err= -1
      if (  Theo_case_S == 'MTN_SCHAR' &
           .or. Theo_case_S == 'MTN_SCHAR2' &
           .or. Theo_case_S == 'MTN_PINTY' &
           .or. Theo_case_S == 'MTN_PINTY2' &
           .or. Theo_case_S == 'MTN_PINTYNL' &
           .or. Theo_case_S == 'NOFLOW' ) then
         err = mtn_cfg ()
      else if (  Theo_case_S == 'BUBBLE' ) then
         err = bubble_cfg ()
      else
         if (Lun_out > 0) then
            write (Lun_out, 9200) Theo_case_S
         end if
      end if

      call gem_error (err, 'theo_cfg', Theo_case_S)

      if (Lun_out > 0) write (Lun_out, 7050) Theo_case_S

 7050 format (/' THEORETICAL CASE IS: ',a/)
 9200 format (/,' Unsupported theoretical case: ',a/)

      return
      end subroutine theo_cfg
!
!     ---------------------------------------------------------------
!
      subroutine theo_data
      use dynkernel_options
      use gem_options
      use glb_ld
      use gmm_itf_mod
      use gmm_geof
      use gmm_pw
      use gmm_vt1
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer :: istat,i,j,k,err
      real(kind=REAL64)  :: pp, ex
!
!---------------------------------------------------------------------
!
      istat = gmm_get (gmmk_pw_uu_plus_s, pw_uu_plus)
      istat = gmm_get (gmmk_pw_vv_plus_s, pw_vv_plus)
      istat = gmm_get (gmmk_pw_tt_plus_s, pw_tt_plus)
      istat = gmm_get (gmmk_ut1_s ,ut1 )
      istat = gmm_get (gmmk_vt1_s ,vt1 )
      istat = gmm_get (gmmk_wt1_s ,wt1 )
      istat = gmm_get (gmmk_tt1_s ,tt1 )
      istat = gmm_get (gmmk_zdt1_s,zdt1)
      istat = gmm_get (gmmk_st1_s ,st1 )
      istat = gmm_get (gmmk_sls_s ,sls )
      istat = gmm_get (gmmk_fis0_s,fis0)
      istat = gmm_get (gmmk_qt1_s ,qt1 )

      wt1= 0. ; zdt1= 0.

      if (      Theo_case_S == 'MTN_SCHAR'   &
           .or. Theo_case_S == 'MTN_SCHAR2'  &
           .or. Theo_case_S == 'MTN_PINTY'   &
           .or. Theo_case_S == 'MTN_PINTY2'  &
           .or. Theo_case_S == 'MTN_PINTYNL' &
           .or. Theo_case_S == 'NOFLOW' ) then

         call mtn_data ( pw_uu_plus, pw_vv_plus, pw_tt_plus, st1, qt1, fis0, orols, sls, &
                         l_minx, l_maxx, l_miny, l_maxy, G_nk, Theo_case_S )
         call hwnd_stag ( ut1,vt1, pw_uu_plus,pw_vv_plus,&
                          l_minx,l_maxx,l_miny,l_maxy,G_nk,.true. )

      elseif ( Theo_case_S == 'BUBBLE' ) then

         fis0= 0.0 ; sls= 0.0
         pw_uu_plus= 0.0 ; pw_vv_plus= 0.0
         ut1= 0.0 ; vt1= 0.0 ; st1= 0.

         if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P') then

            call bubble_fislP_data ( pw_tt_plus, l_minx,l_maxx,l_miny,l_maxy,G_nk )
            qt1= 0.

         else if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then

            call bubble_fislH_data ( pw_tt_plus, l_minx,l_maxx,l_miny,l_maxy,G_nk )
            do k=1,g_nk+1
               do j= 1-G_haloy, l_nj+G_haloy
                  do i= 1-G_halox, l_ni+G_halox
                     ex=1.d0-grav_8/(cpd_8*bubble_theta)*Ver_z_8%m(k)
                     pp=1.d5*ex**(1.d0/cappa_8)
                     qt1(i,j,k)=rgasd_8*Cstv_Tstr_8*log(pp/1.d5)+grav_8*Ver_z_8%m(k)
                  end do
               end do
            end do

         end if

         err=0
         if ((l_north).and.(l_nj-2*pil_n+1<1)) err=-1
         if ((l_east ).and.(l_ni-2*pil_e+1<1)) err=-1
         if ((l_south).and.(2*pil_s>l_nj)    ) err=-1
         if ((l_west ).and.(2*pil_w>l_ni)    ) err=-1
         call gem_error(err,'ABORT in theo_data',&
                  'Partitionning NOT allowed for MIRROR')

      else

         call gem_error(-1,'WRONG THEO CASE in theo_data',Theo_case_S)

      end if

      call tt2virt ( tt1, .true., l_minx,l_maxx,l_miny,l_maxy, G_nk )
!
!---------------------------------------------------------------------
      return
      end subroutine theo_data

end module theo_options
