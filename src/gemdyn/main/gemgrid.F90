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

!**s/r gemgrid - grille program
      subroutine gemgrid()
      use clib_itf_mod
      use step_options
      use mem_nest
      use glb_ld
      use HORgrid_options
      use VERgrid_options
      use lam_options
      use hgc
      use lun
      use path
      use geomh
      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>
#include <rmnlib_basics.hf>
      integer, external :: gemdm_config, domain_decomp
      character(len=120) :: ofile,ofileU,ofileV,ofileR,etk,etk_ext
      character(len=2024) :: fn
      logical :: radians
      integer :: unf,unf1,unf2,unf3,unf4,unf5,unf6,err,npack

      logical, parameter :: gauss_L = .false.
      integer :: i,j,i0,j0,ip1,ip2
      integer :: Grd_ip1,Grd_ip2,Grd_ip3,ni,nj, in,jn
      real, dimension(:), allocatable :: xposu, yposv, xpos,ypos
      real, dimension(:,:), allocatable :: mask
      real(kind=REAL64), dimension(:), allocatable :: x_8, y_8
      real(kind=REAL64), parameter :: HALF_8  = 0.5d0

!
!----------------------------------------------------------------------
!
      call init_component()

      etk = 'PARPOS'
      fn  = trim(Path_input_S)//'/model_settings.nml'
      Step_dt = 1.
      radians = .false.


      ofile = 'tape1'
      if (Grd_yinyang_L .and. Grd_yinyang_S == 'YAN') ofile= 'tape2'

      ofileU = 'tapeU1'
      if (Grd_yinyang_L .and. Grd_yinyang_S == 'YAN') ofileU= 'tapeU2'

      ofileV = 'tapeV1'
      if (Grd_yinyang_L .and. Grd_yinyang_S == 'YAN') ofileV= 'tapeV2'

      ofileR = 'tapeRot'
      if (Grd_yinyang_L .and. Grd_yinyang_S == 'YAN') ofileR= 'tapeRot'

      unf = 0
      if (fnom (unf,fn, 'SEQ+OLD', 0) == 0) then
         if (HORgrid_nml(unf) < 0) then
            print *,'STOP: problem with NAMELIST &grid'
            print *,"Use checknml to verify: \'checknml grid\'"
            stop
         else
            if (HORgrid_config (1) < 0) then
               print *,'STOP: problem with HORgrid_config'
               stop
            endif
         endif
         if (VERgrid_nml(unf) < 0) then
            print *,'STOP: problem with NAMELIST &vert_layers'
            print *,"Use checknml to verify: \'checknml vert_layers\'"
            stop
         else
            if (VERgrid_config () < 0) then
               print *,'STOP: problem with VERgrid_config'
               stop
            endif
         endif
         err= fclos(unf)
      else
         print *, ' Namelist FILE: ',trim(fn),' NOT AVAILABLE'
         stop
      endif

      Step_runstrt_S = '20160825.000000'
      err = HORgrid_nml (-1)
      err = VERgrid_nml (-1)
      Step_nesdt= 1.

      err = gemdm_config ()

      G_ni = Grd_ni ; G_nj = Grd_nj ; G_nk = 2

      allocate ( x_8(Grd_ni+1), y_8(Grd_nj+1) , &
                 xpos (Grd_ni) , ypos (Grd_nj), &
                 xposU(Grd_ni ), yposV(Grd_nj) )

      call set_gemHgrid4 ( x_8, y_8, G_ni, G_nj, Grd_dx, Grd_dy  , &
                           Grd_x0_8, Grd_xl_8, Grd_y0_8, Grd_yl_8, &
                           Grd_yinyang_L )
      x_8(G_ni+1) = x_8(G_ni)+ Grd_dx
      y_8(G_nj+1) = y_8(G_nj)+ Grd_dy

      xpos(1:G_ni) = x_8(1:G_ni)
      ypos(1:G_nj) = y_8(1:G_nj)

      write(*,*) 'LONGITUDE'
      write(*,778)(i,xpos(i),i=1,G_ni)
      write(*,*) 'LATITUDE'
      write(*,778)(i,ypos(i),i=1,G_nj)

!Mass grid
      call set_igs2 ( ip1,ip2, xpos,ypos,G_ni,G_nj             ,&
                      Hgc_ig1ro,Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro,&
                      1,G_ni,1,1,G_nj,1)

      Grd_ip3 = 0

      err = clib_remove(ofile)

      unf1= 0
      if (fnom(unf1,ofile,'RND',0) >= 0) then
         err= fstouv (unf1, 'RND')
      else
         print *,'problem opening', trim(ofile)
         stop
      endif

      if (Grd_yinyang_L) then
         etk_ext=trim(etk)//'_'//trim(Grd_yinyang_S)
      else
         etk_ext=trim(etk)
      endif

      npack = -32

      err = fstecr ( xpos,xpos, npack, unf1, 0, 0, 0, G_ni, 1, 1 , &
                    ip1,ip2,Grd_ip3,'X','>>',etk_ext,Hgc_gxtyp_s, &
                    Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro, 5, .true. )
      err = fstecr ( ypos,ypos, npack, unf1, 0, 0, 0, 1, G_nj, 1 , &
                    ip1,ip2,Grd_ip3,'X','^^',etk_ext,Hgc_gxtyp_s, &
                    Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro, 5, .true. )

      err = fstfrm(unf1)
      err = fclos (unf1)

      do i= 1, G_ni
         xposU(i) = (x_8(i) + x_8(i+1)) * HALF_8
      enddo
      do j= 1, G_nj
         yposV(j) = (y_8(j) + y_8(j+1)) * HALF_8
      enddo

!U grid
      call set_igs2 ( ip1,ip2, xposU,ypos,G_ni,G_nj             ,&
                      Hgc_ig1ro,Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro,&
                      1,G_ni,1,1,G_nj,1)

      err = clib_remove(ofileU)

      unf2= 0
      if (fnom(unf2,ofileU,'RND',0) >= 0) then
         err= fstouv (unf2, 'RND')
      else
         print *,'problem opening', trim(ofileU)
         stop
      endif

      err = fstecr ( xposU,xposU, npack, unf2, 0, 0, 0, G_ni, 1, 1 , &
                    ip1,ip2,Grd_ip3,'X','>>',etk_ext,Hgc_gxtyp_s, &
                    Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro, 5, .true. )
      err = fstecr ( ypos,ypos, npack, unf2, 0, 0, 0, 1, G_nj, 1 , &
                    ip1,ip2,Grd_ip3,'X','^^',etk_ext,Hgc_gxtyp_s, &
                    Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro, 5, .true. )

      err = fstfrm(unf2)
      err = fclos (unf2)
! V grid
      err = clib_remove(ofileV)

      call set_igs2 ( ip1,ip2, xpos,yposV,G_ni,G_nj             ,&
                      Hgc_ig1ro,Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro,&
                      1,G_ni,1,1,G_nj,1)
      unf3= 0
      if (fnom(unf3,ofileV,'RND',0) >= 0) then
         err= fstouv (unf3, 'RND')
      else
         print *,'problem opening', trim(ofileV)
         stop
      endif

      err = fstecr ( xpos,xpos, npack, unf3, 0, 0, 0, G_ni, 1, 1 , &
                    ip1,ip2,Grd_ip3,'X','>>',etk_ext,Hgc_gxtyp_s, &
                    Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro, 5, .true. )
      err = fstecr ( yposV,yposV, npack, unf3, 0, 0, 0, 1, G_nj, 1 , &
                    ip1,ip2,Grd_ip3,'X','^^',etk_ext,Hgc_gxtyp_s, &
                    Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro, 5, .true. )

      err = fstfrm(unf3)
      err = fclos (unf3)

!F grid
      err = clib_remove(ofileR)

      call set_igs2 ( ip1,ip2, xposU,yposV,G_ni,G_nj             ,&
                      Hgc_ig1ro,Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro,&
                      1,G_ni,1,1,G_nj,1)
      unf3= 0
      if (fnom(unf3,ofileR,'RND',0) >= 0) then
         err= fstouv (unf3, 'RND')
      else
         print *,'problem opening', trim(ofileR)
         stop
      endif

      err = fstecr ( xposU,xposU, npack, unf3, 0, 0, 0, G_ni, 1, 1 , &
                    ip1,ip2,Grd_ip3,'X','>>',etk_ext,Hgc_gxtyp_s, &
                    Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro, 5, .true. )
      err = fstecr ( yposV,yposV, npack, unf3, 0, 0, 0, 1, G_nj, 1 , &
                    ip1,ip2,Grd_ip3,'X','^^',etk_ext,Hgc_gxtyp_s, &
                    Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro, 5, .true. )

      err = fstfrm(unf3)
      err = fclos (unf3)


      err = domain_decomp (1, 1, .false.)
      call set_gmm ()
      call nest_set_mem ()
      unf4=0
      if (fnom(unf4,trim(ofile)//'_core','RND',0) >= 0) then
         err= fstouv (unf4, 'RND')
      else
         print *,'problem opening', trim(ofile//'_core')
         stop
      endif

      i0 = 1      + Grd_extension
      in = G_ni - Grd_extension
      j0 = 1      + Grd_extension
      jn = G_nj - Grd_extension
      ni = in-i0+1
      nj = jn-j0+1

      xpos(1:ni) = x_8(i0:in)
      ypos(1:nj) = y_8(j0:jn)

      call set_igs2 ( Grd_ip1,Grd_ip2, xpos,ypos,ni,nj, &
                      Hgc_ig1ro,Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro, &
                      1,ni,1,1,nj,1)
      err= fstecr ( xpos,xpos, npack, unf4, 0, 0, 0, ni, 1, 1, &
                    Grd_ip1,Grd_ip2,Grd_ip3,'X','>>',etk_ext,Hgc_gxtyp_s, &
                    Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro, 5, .true. )
      err= fstecr ( ypos,ypos, npack, unf4, 0, 0, 0, 1, nj, 1, &
                    Grd_ip1,Grd_ip2,Grd_ip3,'X','^^',etk_ext,Hgc_gxtyp_s, &
                    Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro, 5, .true. )

      unf5 = 0
      err  = fnom(unf5,ofile,'RND',0)
      err  = fstouv (unf5, 'RND')

      allocate (mask(G_ni, G_nj))

      do j=1,G_nj
         do i=1,G_ni
            mask(i,j)= 1.-nest_weightm(i,j,G_nk+1)
         enddo
      enddo

      err = fstecr ( mask,mask, npack, unf5, 0, 0, 0, G_ni, G_nj, 1, &
                     0,0,0,'X','MSKC',etk_ext,'Z'    , &
                     ip1,ip2,Grd_ip3,0, 5, .true. )
      deallocate (mask)

      err= fstfrm(unf4)
      err= fclos (unf4)
      err= fstfrm(unf5)
      err= fclos (unf5)

      unf6=0
      if (fnom(unf6,trim(ofile)//'_free','RND',0) >= 0) then
         err= fstouv (unf6, 'RND')
      else
         print *,'problem opening', trim(ofile//'_core')
         stop
      endif

      i0= 1    + Grd_extension + Lam_blend_H
      in= G_ni - Grd_extension - Lam_blend_H
      j0= 1    + Grd_extension + Lam_blend_H
      jn= G_nj - Grd_extension - Lam_blend_H
      ni = in-i0+1
      nj = jn-j0+1

      xpos(1:ni) = x_8(i0:in)
      ypos(1:nj) = y_8(j0:jn)

      call set_igs2 ( Grd_ip1,Grd_ip2, &
                      xpos,ypos,ni,nj, &
                      Hgc_ig1ro,Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro, &
                      1,ni,1,1,nj,1)
      err = fstecr ( xpos,xpos, npack, unf6, 0, 0, 0, ni, 1, 1, &
                     Grd_ip1,Grd_ip2,Grd_ip3,'X','>>',etk_ext,Hgc_gxtyp_s, &
                     Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro, 5, .true. )
      err = fstecr ( ypos,ypos, npack, unf6, 0, 0, 0, 1, nj, 1, &
                     Grd_ip1,Grd_ip2,Grd_ip3,'X','^^',etk_ext,Hgc_gxtyp_s, &
                     Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro, 5, .true. )

      err = fstfrm(unf6)
      err = fclos (unf6)

      deallocate (x_8, y_8, xpos, ypos)
      deallocate (xposU, yposV )


      call gemtim4 ( Lun_out, 'AFTER set_opr', .false. )

      call memusage (6)

      call rpn_comm_FINALIZE(err)
!
!-------------------------------------------------------------------
!

 778  format(4(i5,e15.5))
      end
