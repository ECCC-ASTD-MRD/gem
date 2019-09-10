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
module HORgrid_options
   use lun
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

   !# Type of grid described using 2 characters:
   !# * "GY" : Global Yin-Yang
   !# * "LU" : LAM    Uniform
   character(len=2) ::  Grd_typ_S = ''
   namelist /grid/ Grd_typ_S

   !# Number of points along NI
   integer :: Grd_ni = 0
   namelist /grid/ Grd_ni

   !# Number of points along NJ
   integer :: Grd_nj = 0
   namelist /grid/ Grd_nj

   !# Max Supported Courrant number;
   !# Pilot area = Grd_maxcfl_fact * Grd_maxcfl+Grd_bsc_base+Grd_bsc_ext1
   integer :: Grd_maxcfl = 1
   namelist /grid/ Grd_maxcfl

   !# (LU only) Mesh length (resolution) in x-direction (degrees)
   real  :: Grd_dx = 0.
   namelist /grid/ Grd_dx

   !# (LU only) Mesh length (resolution) in y-direction (degrees)
   real  :: Grd_dy = 0.
   namelist /grid/ Grd_dy

   !# Reference Point I
   integer :: Grd_iref = -1
   namelist /grid/ Grd_iref

   !# Reference Point J
   integer :: Grd_jref = -1
   namelist /grid/ Grd_jref

   !# Latitude on rotated grid of reference point (degrees)
   real  :: Grd_latr = 0.
   namelist /grid/ Grd_latr

   !# Longitude on rotated grid of reference point (degrees)
   real  :: Grd_lonr = 180.
   namelist /grid/ Grd_lonr

   !# (GY only) Overlap extent along latitude axis for GY grid (degrees)
   real  :: Grd_overlap = 0.
   namelist /grid/ Grd_overlap

   !# Geographic latitude of the center of the computational domain (degrees)
   real  :: Grd_xlon1 = 180.
   namelist /grid/ Grd_xlon1

   !# Geographic longitude of the center of the computational domain (degrees)
   real  :: Grd_xlat1 = 0.
   namelist /grid/ Grd_xlat1

   !# Geographic longitude of a point on the equator of the computational domain
   !# east of Grd_xlon1,Grd_xlat1  (degrees)
   real  :: Grd_xlon2 = 270.
   namelist /grid/ Grd_xlon2

   !# Geographic latitude of a point on the equator of the computational domain
   !# east of  Grd_xlon1,Grd_xlat1  (degrees)
   real  :: Grd_xlat2 = 0.
   namelist /grid/ Grd_xlat2

   character(len=12) :: Grd_yinyang_S
   logical Grd_roule, Grd_yinyang_L
   integer Grd_bsc_base, Grd_bsc_adw, Grd_bsc_ext1, Grd_extension
   integer Grd_ndomains, Grd_mydomain
   integer Grd_local_gid, Grd_lclcore_gid, Grd_global_gid, &
            Grd_lphy_gid, Grd_glbcore_gid
   integer Grd_lphy_i0, Grd_lphy_in, Grd_lphy_j0, Grd_lphy_jn, &
           Grd_lphy_ni, Grd_lphy_nj
   real(kind=REAL64) Grd_rot_8(3,3), Grd_x0_8, Grd_xl_8, Grd_y0_8, Grd_yl_8

contains

!**s/r HORgrid_nml - Reading namelist grid

      integer function HORgrid_nml (F_unf)
      use ctrl
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_unf

      logical nml_must
      character(len=64) :: nml_S
!
!-------------------------------------------------------------------
!
! boiler plate - start
      if ( F_unf < 0 ) then
         HORgrid_nml= 0
         if ( Lun_out >= 0) write (Lun_out,nml=grid)
         return
      end if

      if (Ctrl_theoc_L) then
         HORgrid_nml= 1
         return
      end if

      HORgrid_nml= -1 ; nml_must= .true. ; nml_S= 'grid'

      rewind(F_unf)
      read (F_unf, nml=grid, end= 1001, err=1003)
      HORgrid_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         HORgrid_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      else
         if (Lun_out >= 0) write (Lun_out, 6009) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (HORgrid_nml < 0 ) return

      if ((Lun_out>=0).and.(HORgrid_nml==0)) &
             write (Lun_out, 6004) trim(nml_S)

      Grd_maxcfl= max(1,Grd_maxcfl)

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
 6009 format (//,' NAMELIST ',A,' IS MANDATORY'//)
! boiler plate - end
!
!-------------------------------------------------------------------
!
      return
      end function HORgrid_nml

!**s/r HORgrid_config - Configure horizontal grid parameters

      integer function HORgrid_config (F_adv_maxcfl_fact)
      use glb_ld
      use glb_pil
      use hgc
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_adv_maxcfl_fact

      integer, external ::  yyg_checkrot

      character(len=120) :: dumc
      logical :: almost_zero
      integer err
      real(kind=REAL64) :: a_8, b_8, c_8, d_8
      real(kind=REAL64), dimension(3) :: xyz1_8, xyz2_8
      real(kind=REAL64) :: yan_xlat1_8, yan_xlon1_8, yan_xlat2_8, yan_xlon2_8
      real(kind=REAL64), parameter :: epsilon = 1.0d-5
!
!-------------------------------------------------------------------
!
      HORgrid_config = -1

      call low2up (Grd_typ_S,dumc)
      Grd_typ_S = dumc(1:2)

      if (Grd_typ_S(1:2) == 'GY') then

         if ((Grd_ni > 0).and.(Grd_nj > 0)) then
            if (Lun_out > 0) then
               write(Lun_out,'(/2(2x,a/))')  &
               'CONFLICTING Grd_NI & Grd_NJ IN NAMELIST grid',&
               '            only one of them can be > 0'
            end if
            return
         end if
         if (Grd_ni <= 0) then
            if (Grd_nj > 0) Grd_ni= (Grd_nj-1)*3 + 1
         end if
         if (Grd_nj <= 0) then
            if (Grd_ni > 0) Grd_nj= nint ( real(Grd_ni-1)/3. + 1. )
         end if

         if (yyg_checkrot() < 0) return

         if (trim(Grd_yinyang_S) == 'YAN') then
            call yyg_yangrot ( dble(Grd_xlat1), dble(Grd_xlon1), &
                               dble(Grd_xlat2), dble(Grd_xlon2), &
            yan_xlat1_8, yan_xlon1_8, yan_xlat2_8, yan_xlon2_8 )
            Grd_xlat1  = yan_xlat1_8 ; Grd_xlon1  = yan_xlon1_8
            Grd_xlat2  = yan_xlat2_8 ; Grd_xlon2  = yan_xlon2_8
         end if

      else

         if ( almost_zero(Grd_dx*Grd_dy) ) then
            if (Lun_out > 0) then
               write(Lun_out,*) 'VERIFY Grd_DX & Grd_DY IN NAMELIST grid'
            end if
            return
         end if

      end if

      if (Grd_ni*Grd_nj <= 0) then
         if (Lun_out > 0) then
            write(Lun_out,*) 'VERIFY Grd_NI & Grd_NJ IN NAMELIST grid'
         end if
         return
      end if

      Grd_x0_8=  0.0 ; Grd_xl_8=360.0
      Grd_y0_8=-90.0 ; Grd_yl_8= 90.0

      call gem_grid_param (Grd_bsc_base,Grd_bsc_ext1,Grd_extension,&
                           F_adv_maxcfl_fact*Grd_maxcfl           ,&
                           Grd_iref, Grd_jref, Grd_lonr, Grd_latr ,&
                           Grd_ni, Grd_nj, Grd_dx, Grd_dy         ,&
                           Grd_x0_8, Grd_y0_8, Grd_xl_8, Grd_yl_8 ,&
                           Grd_overlap, Grd_yinyang_L, Lun_out, err)
      if (err < 0) then
         if (Lun_out > 0) then
            write(Lun_out,*) 'ERROR in gem_grid_param'
         end if
         return
      end if

      Glb_pil_n = Grd_extension
      Glb_pil_s = Glb_pil_n
      Glb_pil_w = Glb_pil_n
      Glb_pil_e = Glb_pil_n

      pil_w = 0 ; pil_n = 0 ; pil_e = 0 ; pil_s = 0
      if (l_west ) pil_w= Glb_pil_w
      if (l_north) pil_n= Glb_pil_n
      if (l_east ) pil_e= Glb_pil_e
      if (l_south) pil_s= Glb_pil_s

      Lam_pil_w= Glb_pil_w
      Lam_pil_n= Glb_pil_n
      Lam_pil_e= Glb_pil_e
      Lam_pil_s= Glb_pil_s

!     compute RPN/FST grid descriptors

      Hgc_gxtyp_s = 'E'
      call cxgaig ( Hgc_gxtyp_S,Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro,&
                              Grd_xlat1,Grd_xlon1,Grd_xlat2,Grd_xlon2 )
      call cigaxg ( Hgc_gxtyp_S,Grd_xlat1,Grd_xlon1,Grd_xlat2,Grd_xlon2,&
                              Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro )

      if (Lun_out > 0) then
         write(Lun_out,1100) trim(Grd_yinyang_S), &
                       Grd_ni, Grd_x0_8, Grd_xl_8, &
                       Grd_nj, Grd_y0_8, Grd_yl_8, &
                       Grd_typ_S, Grd_dx ,Grd_dy , &
                       Grd_dx*40000./360.,Grd_dy*40000./360.
      end if

      if (Lun_out > 0) write (Lun_out,1004) &
                     Grd_xlat1,Grd_xlon1,Grd_xlat2,Grd_xlon2,&
                     Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro

      Grd_roule = .not. ( (abs(Grd_xlon1-180.d0) < epsilon) .and. &
                        (  abs(Grd_xlon2-270.d0) < epsilon) .and. &
                        (  abs(Grd_xlat1       ) < epsilon) .and. &
                        (  abs(Grd_xlat2       ) < epsilon) )

      Grd_rot_8      = 0.
      Grd_rot_8(1,1) = 1.
      Grd_rot_8(2,2) = 1.
      Grd_rot_8(3,3) = 1.

      if (Grd_roule) then
!
!     Compute the rotation matrix that allows transformation
!     from the none-rotated to the rotated spherical coordinate system.
!
!     Compute transform matrices xyz1_8 and xyz2_8
!
         call llacar ( xyz1_8, Grd_xlon1, Grd_xlat1, 1, 1 )
         call llacar ( xyz2_8, Grd_xlon2, Grd_xlat2, 1, 1 )
!
!     Compute a = cos(alpha) & b = sin(alpha)
!
         a_8 = (xyz1_8(1)*xyz2_8(1)) + (xyz1_8(2)*xyz2_8(2))  &
                                     + (xyz1_8(3)*xyz2_8(3))
         b_8 = sqrt (((xyz1_8(2)*xyz2_8(3)) - (xyz2_8(2)*xyz1_8(3)))**2 &
                  +  ((xyz2_8(1)*xyz1_8(3)) - (xyz1_8(1)*xyz2_8(3)))**2 &
                  +  ((xyz1_8(1)*xyz2_8(2)) - (xyz2_8(1)*xyz1_8(2)))**2)
!
!     Compute c = norm(-r1) & d = norm(r4)
!
         c_8 = sqrt ( xyz1_8(1)**2 + xyz1_8(2)**2 + xyz1_8(3)**2 )
         d_8 = sqrt ( ( ( (a_8*xyz1_8(1)) - xyz2_8(1) ) / b_8 )**2 + &
                      ( ( (a_8*xyz1_8(2)) - xyz2_8(2) ) / b_8 )**2 + &
                      ( ( (a_8*xyz1_8(3)) - xyz2_8(3) ) / b_8 )**2  )

         Grd_rot_8(1,1)=  -xyz1_8(1)/c_8
         Grd_rot_8(1,2)=  -xyz1_8(2)/c_8
         Grd_rot_8(1,3)=  -xyz1_8(3)/c_8
         Grd_rot_8(2,1)=  ( ((a_8*xyz1_8(1)) - xyz2_8(1)) / b_8)/d_8
         Grd_rot_8(2,2)=  ( ((a_8*xyz1_8(2)) - xyz2_8(2)) / b_8)/d_8
         Grd_rot_8(2,3)=  ( ((a_8*xyz1_8(3)) - xyz2_8(3)) / b_8)/d_8
         Grd_rot_8(3,1)=   &
              ( (xyz1_8(2)*xyz2_8(3)) - (xyz2_8(2)*xyz1_8(3)))/b_8
         Grd_rot_8(3,2)=   &
              ( (xyz2_8(1)*xyz1_8(3)) - (xyz1_8(1)*xyz2_8(3)))/b_8
         Grd_rot_8(3,3)=   &
              ( (xyz1_8(1)*xyz2_8(2)) - (xyz2_8(1)*xyz1_8(2)))/b_8

      end if

      G_ni  = Grd_ni
      G_nj  = Grd_nj

      G_niu = G_ni
      G_njv = G_nj - 1

      G_niu = G_ni - 1

      HORgrid_config = 1

 1004 format (1x,'Numerical equator: ',2('(',f9.4,',',f9.4,') ')/&
                21x,'IG1-4: ',4i8)
 1100 FORMAT (/,1X,'FINAL HORIZONTAL GRID CONFIGURATION: UNIFORM RESOLUTION: ',a, &
        /1X,' NI=',I5,' FROM x0=',F11.5,' TO xl=',F11.5,' DEGREES' &
        /1X,' NJ=',I5,' FROM y0=',F11.5,' TO yl=',F11.5,' DEGREES' &
        /1X,' GRIDTYPE= ',a,'     DX= ',F11.5,'   DY= ',F11.5,' degrees' &
        /14X,               '     DX= ',F11.5,'   DY= ',F11.5 ' km'/)
!
!-------------------------------------------------------------------
!
      return
      end function HORgrid_config

   !/@*
   function HORgrid_options_init() result(F_istat)
      use, intrinsic :: iso_fortran_env
      implicit none
      !@object Additional initialisation steps before reading the nml
      integer :: F_istat
      !*@/
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      !----------------------------------------------------------------
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.
      !----------------------------------------------------------------
      return
   end function HORgrid_options_init

end module HORgrid_options
