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

!**s/p set_grid - initialization of common block GRID

      integer function set_grid (F_argc,F_argv_S,F_cmdtyp_S,F_v1,F_v2)
      use HORgrid_options
      use lam_options
      use glb_ld
      use lun
      use out3
      use glb_pil
      use outgrid
      implicit none
#include <arch_specific.hf>

      integer F_argc,F_v1,F_v2
      character(len=*) F_argv_S(0:F_argc),F_cmdtyp_S

!object
!	initialization of the common block GRID. This function is
!       called when the keyword "grid" is found in the first word
!       of the directives in the input file given in the statement
!       "process_f_callback". This feature is enabled by the
!       ARMNLIB "rpn_fortran_callback" routine (called in "srequet")
!       which allows a different way of passing user directives than
!       the conventional FORTRAN namelist. This function will process
!       the following example command read from the named input file.
!
!   ie: grid=1,model;
!
!       The "rpn_fortran_callback" routine will process the above
!       statement and return 5 arguments to this function. For more
!       information to how this is processed, see "SREQUET".
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_argc       I    - number of elements in F_argv_S
! F_argv_S     I    - array of elements received
! F_cmdtyp_S   I    - character command type - not used
! F_v1         I    - integer parameter 1 - not used
! F_v2         I    - integer parameter 2 - not used
!----------------------------------------------------------------
!
!Notes:
!
! examples:
! grid=4,model;
! grid=2,core;
! grid=1,reduc,4,10,4,10,2
! grid=3,reduc,4,10,4,10
! grid=5,reduc,"NEW",2,11,2,15
! grid=7,core,"CO";
!
! general syntax
! grid=gridid,[model/core/reduc],["etik"],[gridx0,gridx1,gridy0,gridy1];
!
!      gridid  - number to identify gridset to relate to sortie statement
!      model  - total grid of the model,in LAM, this includes pilot area
!      core   - only the uniform part of the grid, in LAM, excludes pilot area
!      free   - only LAM, excludes pilot and blending area
!      reduc   - reduced grid from the model defined as follows
!      gridx0 - starting I value along X
!      gridx1 - ending   I value along X
!      gridy0 - starting J value along X
!      gridy1 - ending   J value along X
!
! IMPORTANT NOTE:
!     Limit the number of definitions for "grid" to improve the efficiency
!     in the output routines. The maximum number of definitions is 4.
!

      integer i, j, gridset,gridout(5)
      integer niout,njout
      character(len=8) grdtyp_S
!
!-------------------------------------------------------------------
!
      if (Lun_out > 0) then
          write(Lun_out,*)
          write(Lun_out,*) F_argv_S(0),'=',F_argv_S(1),',',F_argv_S(2),',',(F_argv_S(i),i=3,F_argc)
      end if
      set_grid = 0
      read(F_argv_S(1),*) gridset
      OutGrid_sets = OutGrid_sets + 1
      if (OutGrid_sets > OUTGRID_MAXGRID1) then
          if (Lun_out > 0) then
            write(Lun_out,*)'SET_GRID WARNING: Too many grid definitions'
          end if
          OutGrid_sets = OutGrid_sets - 1
          set_grid = 1
          return
      end if

      j = OutGrid_sets
      OutGrid_id(j)=gridset

      if(index(F_argv_S(2),'model') /= 0) then
         grdtyp_S='model'
      else if (index(F_argv_S(2),'core') /= 0) then
         grdtyp_S='core'
      else if (index(F_argv_S(2),'free') /= 0) then
         grdtyp_S='free'
      else if (index(F_argv_S(2),'reduc') /= 0) then
         grdtyp_S='reduc'
         gridout(1)= 0 ; gridout(2)= 0
         gridout(3)= 0 ; gridout(4)= 0
         gridout(5)=1
         read(F_argv_S(3),*) gridout(1)
         read(F_argv_S(4),*) gridout(2)
         read(F_argv_S(5),*) gridout(3)
         read(F_argv_S(6),*) gridout(4)
         if (F_argc > 6) then
            if (index(F_argv_S(7),'"') == 0) read(F_argv_S(7),*)gridout(5)
         end if
      else
         if (Lun_out > 0) then
            write(Lun_out,*)'SET_GRID WARNING: Grid Type Undefined'
         end if
         OutGrid_sets = OutGrid_sets - 1
         set_grid = 1
         return
      end if

!    Calculate the origin and outer coordinates of the output grid
!    and set to the maximum/minimum possible

      OutGrid_reduc (j)= (grdtyp_S == 'reduc')
      OutGrid_stride(j)= 1

      if (grdtyp_S == 'model') then

         OutGrid_x0(j)=1
         OutGrid_x1(j)=G_ni
         OutGrid_y0(j)=1
         OutGrid_y1(j)=G_nj

      else if (grdtyp_S == 'core') then

         OutGrid_x0(j)=1      + Lam_pil_w
         OutGrid_x1(j)=Grd_ni - Lam_pil_e
         OutGrid_y0(j)=1      + Lam_pil_s
         OutGrid_y1(j)=Grd_nj - Lam_pil_n

      else if (grdtyp_S == 'free') then

         OutGrid_x0(j)=1       + Lam_pil_w + Lam_blend_Hx
         OutGrid_x1(j)=Grd_ni  - Lam_pil_e - Lam_blend_Hx
         OutGrid_y0(j)=1       + Lam_pil_s + Lam_blend_Hy
         OutGrid_y1(j)=Grd_nj  - Lam_pil_n - Lam_blend_Hy

      else if (grdtyp_S == 'reduc') then

         OutGrid_x0(j)=min( G_ni,      max(1,gridout(1)) )
         OutGrid_x1(j)=max( OutGrid_x0(j), min(G_ni,gridout(2)) )
         OutGrid_y0(j)=min( G_nj, max(1,gridout(3)) )
         OutGrid_y1(j)=max( OutGrid_y0(j), min(G_nj,gridout(4)) )
         OutGrid_stride(j)=min( max(gridout(5),1), &
             min(OutGrid_x1(j)-OutGrid_x0(j)+1,OutGrid_y1(j)-OutGrid_y0(j)+1)/2-1 )
         if (Grd_yinyang_L) then
            if (trim(Grd_yinyang_S) == 'YAN') then
               OutGrid_x0(j)= 0 ; OutGrid_x1(j)= -1
            end if
         end if

      end if

      niout=OutGrid_x1(j) - OutGrid_x0(j) + 1
      njout=OutGrid_y1(j) - OutGrid_y0(j) + 1
      if (niout < 1.or.njout < 1) then
          OutGrid_sets = OutGrid_sets - 1
          if (.not. ((Grd_yinyang_L).and.(trim(Grd_yinyang_S) == 'YAN'))) then
             if (Lun_out > 0) then
               write (Lun_out,*)'ERROR in description of output grid!'
             end if
          end if
          return
      end if

!      if (Lun_out > 0) write(Lun_out,*) ' Grid_set(',j,') : OutGrid_id=',OutGrid_id(j)
!
!-------------------------------------------------------------------
!
      return
      end
