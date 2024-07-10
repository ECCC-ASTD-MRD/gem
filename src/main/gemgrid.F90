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
      use dynkernel_options
      use HORgrid_options
      use VERgrid_options
      use lam_options
      use hgc
      use lun
      use path
      use geomh
      use rmn_fst24
      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>
#include <rmnlib_basics.hf>
      integer, external :: gemdm_config, domain_decomp
      character(len=120) :: ofile,ofileU,ofileV,ofileR,etk,etk_ext
      character(len=2024) :: fn
      logical :: radians
      integer :: unf,err

      type(fst_file)   :: file
      type(fst_record) :: rec
      logical          :: success
 
      logical, parameter :: gauss_L = .false.
      integer :: i,j,i0,j0,ip1,ip2
      integer :: Grd_ip1,Grd_ip2,Grd_ip3,ni,nj, in,jn
      real, dimension(:), allocatable, target :: xposu, yposv, xpos,ypos
      real, dimension(:,:), allocatable, target :: mask
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
      
      rec%etiket=etk_ext
      rec%dateo=0
      rec%deet=0
      rec%npas=0
      rec%data_type=5
      rec%data_bits=32
      rec%pack_bits=32
      rec%nk=1

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
         if (dynKernel_nml(unf) < 0) then
            print *,'STOP: problem with NAMELIST &dyn_kernel'
            print *,"Use checknml to verify: \'checknml dyn_kernel\'"
            stop
         endif
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
      call set_igs2 ( Grd_ip1,Grd_ip2, xpos,ypos,G_ni,G_nj             ,&
                      Hgc_ig1ro,Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro,&
                      1,G_ni,1,1,G_nj,1)
      Grd_ip3 = 0

      err = clib_remove(ofile)

      if (.not. file%open(ofile,'RND+R/W')) then
         print *,'problem opening', trim(ofile)
         stop
      endif

      if (Grd_yinyang_L) then
         etk_ext=trim(etk)//'_'//trim(Grd_yinyang_S)
      else
         etk_ext=trim(etk)
      endif

      rec%nomvar='>>'
      rec%typvar='X'
      rec%grtyp=Hgc_gxtyp_s
      rec%ip1=Grd_ip1
      rec%ip2=Grd_ip2
      rec%ip3=Grd_ip3
      rec%ni=G_ni
      rec%nj=1
      rec%ig1=Hgc_ig1ro
      rec%ig2=Hgc_ig2ro
      rec%ig3=Hgc_ig3ro
      rec%ig4=Hgc_ig4ro
      rec%data=c_loc(xpos(:))
      success = file%write(rec,rewrite=FST_SKIP)

      rec%nomvar='^^'
      rec%ni=1
      rec%nj=G_nj
      rec%data=c_loc(ypos(:))
      success = file%write(rec,rewrite=FST_SKIP)

      success = file%close()

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

      if (.not. file%open(ofileU,'RND+R/W')) then
         print *,'problem opening', trim(ofileU)
         stop
      endif

      rec%nomvar='>>'
      rec%typvar='X'
      rec%grtyp=Hgc_gxtyp_s
      rec%ip1=ip1
      rec%ip2=ip2
      rec%ip3=Grd_ip3
      rec%ni=G_ni
      rec%nj=1
      rec%ig1=Hgc_ig1ro
      rec%ig2=Hgc_ig2ro
      rec%ig3=Hgc_ig3ro
      rec%ig4=Hgc_ig4ro
      rec%data=c_loc(xposU(:))
      success = file%write(rec,rewrite=FST_SKIP)

      rec%nomvar='^^'
      rec%ni=1
      rec%nj=G_nj
      rec%data=c_loc(ypos(:))
      success = file%write(rec,rewrite=FST_SKIP)

      success = file%close()

! V grid
      err = clib_remove(ofileV)

      call set_igs2 ( ip1,ip2, xpos,yposV,G_ni,G_nj             ,&
                      Hgc_ig1ro,Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro,&
                      1,G_ni,1,1,G_nj,1)
      if (.not. file%open(ofileV,'RND+R/W')) then
         print *,'problem opening', trim(ofileV)
         stop
      endif

      rec%nomvar='>>'
      rec%typvar='X'
      rec%grtyp=Hgc_gxtyp_s
      rec%ip1=ip1
      rec%ip2=ip2
      rec%ip3=Grd_ip3
      rec%ni=G_ni
      rec%nj=1
      rec%ig1=Hgc_ig1ro
      rec%ig2=Hgc_ig2ro
      rec%ig3=Hgc_ig3ro
      rec%ig4=Hgc_ig4ro
      rec%data=c_loc(xpos(:))
      success = file%write(rec,rewrite=FST_SKIP)

      rec%nomvar='^^'
      rec%ni=1
      rec%nj=G_nj
      rec%data=c_loc(yposV(:))
      success = file%write(rec,rewrite=FST_SKIP)

      success = file%close()

!F grid
      err = clib_remove(ofileR)

      call set_igs2 ( ip1,ip2, xposU,yposV,G_ni,G_nj             ,&
                      Hgc_ig1ro,Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro,&
                      1,G_ni,1,1,G_nj,1)

      if (.not. file%open(ofileR,'RND+R/W')) then
         print *,'problem opening', trim(ofileR)
         stop
      endif

      rec%nomvar='>>'
      rec%typvar='X'
      rec%grtyp=Hgc_gxtyp_s
      rec%ip1=ip1
      rec%ip2=ip2
      rec%ip3=Grd_ip3
      rec%ni=G_ni
      rec%nj=1
      rec%ig1=Hgc_ig1ro
      rec%ig2=Hgc_ig2ro
      rec%ig3=Hgc_ig3ro
      rec%ig4=Hgc_ig4ro
      rec%data=c_loc(xposU(:))
      success = file%write(rec,rewrite=FST_SKIP)

      rec%nomvar='^^'
      rec%ni=1
      rec%nj=G_nj
      rec%data=c_loc(yposV(:))
      success = file%write(rec,rewrite=FST_SKIP)

      success = file%close()
      
      err = domain_decomp (1, 1, .false.)
      call set_gmm ()
      call nest_set_mem ()

      if (.not. file%open(trim(ofile)//'_core','RND+R/W')) then
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

      call set_igs2 ( ip1,ip2, xpos,ypos,ni,nj, &
                      Hgc_ig1ro,Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro, &
                      1,ni,1,1,nj,1)

      rec%nomvar='>>'
      rec%typvar='X'
      rec%grtyp=Hgc_gxtyp_s
      rec%ip1=ip1
      rec%ip2=ip2
      rec%ip3=Grd_ip3
      rec%ni=G_ni
      rec%nj=1
      rec%ig1=Hgc_ig1ro
      rec%ig2=Hgc_ig2ro
      rec%ig3=Hgc_ig3ro
      rec%ig4=Hgc_ig4ro
      rec%data=c_loc(xpos(:))
      success = file%write(rec,rewrite=FST_SKIP)

      rec%nomvar='^^'
      rec%ni=1
      rec%nj=G_nj
      rec%data=c_loc(ypos)
      success = file%write(rec,rewrite=FST_SKIP)

      success = file%close()

      if (.not. file%open(ofile,'RND+R/W')) then
          print *,'problem opening', ofile
         stop
      endif

      allocate (mask(G_ni, G_nj))

      do j=1,G_nj
         do i=1,G_ni
            mask(i,j)= 1.-nest_weightm(i,j,G_nk+1)
         enddo
      enddo

      rec%nomvar='MSKC'
      rec%typvar='X'
      rec%grtyp='Z'
      rec%ip1=0
      rec%ip2=0
      rec%ip3=0
      rec%ni=G_ni
      rec%nj=G_nj
      rec%ig1=Grd_ip1
      rec%ig2=Grd_ip2
      rec%ig3=Grd_ip3
      rec%ig4=0
      rec%data=c_loc(mask(:,:))
      success = file%write(rec,rewrite=FST_SKIP)

      deallocate (mask)

      success = file%close()

      if (.not. file%open(trim(ofile)//'_free','RND+R/W')) then
         print *,'problem opening', trim(ofile)//'_free'
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

      rec%nomvar='>>'
      rec%typvar='X'
      rec%grtyp=Hgc_gxtyp_s
      rec%ip1=Grd_ip1
      rec%ip2=Grd_ip2
      rec%ip3=Grd_ip3
      rec%ni=ni
      rec%nj=1
      rec%ig1=Hgc_ig1ro
      rec%ig2=Hgc_ig2ro
      rec%ig3=Hgc_ig3ro
      rec%ig4=Hgc_ig4ro
      rec%data=c_loc(xpos(:))
      success = file%write(rec,rewrite=FST_SKIP)

      rec%nomvar='^^'
      rec%ni=1
      rec%nj=nj
      rec%data=c_loc(ypos(:))
      success = file%write(rec,rewrite=FST_SKIP)

      success = file%close()

      deallocate (x_8, y_8, xpos, ypos)
      deallocate (xposU, yposV )

      call rpn_comm_FINALIZE(err)
!
!-------------------------------------------------------------------
!

 778  format(4(i5,e15.5))
      end
