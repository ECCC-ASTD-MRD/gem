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

!**function srequet - reads from a file to initialize output requests
!                     for grid, levels, timesteps and variables

      integer function srequet()
      use HORgrid_options
      use gem_options
      use glb_ld
      use lun
      use out3
      use levels
      use outp
      use outd
      use outc
      use hgc
      use path
      use outgrid
      use timestep
      implicit none
#include <arch_specific.hf>

!object
!    This function uses the ARMNLIB RPN functions
!    RPN_FORTRAN_CALLBACK and PROCESS_F_CALLBACK
!    to process the directives for output listed in a file
!    specified in the call "process_f_callback". This is
!    used instead of the conventional FORTRAN namelist.
!
!notes
!    There are four key functions written to handle the
!    directives from the specified input file:
!
!    set_grid, set_level, set_step, set_var
!
!    which are associated by their
!    respective keywords "grid","levels","steps","sortie". The directives
!    are usually in the form of:
!
!    key1([a,b,c],d,e,f);
!    key2=a,c,[d,e,f];
!    key3=a,b,c,d,<1,15,1>;
!    key4=a,b,c,[d,e,f,g],h,i,<0,6,2>;
!
!    The "rpn_fortran_callback" gives to the written functions a set
!    of arguments (argc,argv,cmd,v1,v2) where:
!
!    cmd,v1,v2 are the same as the last 3 arguments in the calling routine.
!    argc contains the # of elements in argv
!    argv is an array of elements derived from the directives by the
!    "rpn_fortran_callback" routine. Here are some examples on how it
!    will split the directives into an array of elements:
!
!    key1([a,b,c],d,e,f);
!    argc=8
!    argv=key1,[  3],a,b,c,d,e,f
!
!    key4=a,b,c,[d,e,f,g],h,i,<0,6,2>;
!    argc=15
!    argv=key4,a,b,c,[  4],d,e,f,g,h,i,0,2,4,6
!
!    'process_f_callback' will set the input file and activate all the
!    functions declared in "rpn_fortran_callback" calls previously.
!
      integer process_f_callback,longueur
      external process_f_callback,longueur
      integer set_level,set_step,set_grid,set_filt,set_xnbit,set_var,set_conv
      external set_level,set_step,set_grid,set_filt,set_xnbit,set_var,set_conv

      integer p1a,p1b,istat,j
!*
!
      call rpn_f_callback_setverbose(Lun_out)
      OutGrid_sets = 0
      Level_sets = 0
      Level_npres = 0
      Timestep_sets = 0
      Outd_sets = 0
      Outp_sets = 0
      Outc_sets = 0
      Out3_filtpass_max = 0
      Out3_xnbits_max = 0
      p1a =  10
      p1b = -10

      call rpn_fortran_callback('filtre'  ,set_filt ,' ',p1a,p1b)
      call rpn_fortran_callback('convert' ,set_conv ,' ',p1a,p1b)
      call rpn_fortran_callback('xnbit'   ,set_xnbit,' ',p1a,p1b)
      call rpn_fortran_callback('levels'  ,set_level,' ',p1a,p1b)
      call rpn_fortran_callback('grid'    ,set_grid ,' ',p1a,p1b)
      call rpn_fortran_callback('steps'   ,set_step ,' ',p1a,p1b)
      call rpn_fortran_callback('sortie'  ,set_var  ,' ',p1a,p1b)
      call rpn_fortran_callback('sortie_p',set_var  ,' ',p1a,p1b)
      call rpn_fortran_callback('sortie_c',set_var  ,' ',p1a,p1b)

      istat= process_f_callback(trim(Path_outcfg_S))

      if (Lun_out > 0) then
         write(Lun_out,*)' Level_allpres=',(Level_allpres(j),j=1,Level_npres)
         write(Lun_out,*)'SREQUET:Number of warnings =',istat
      end if
      srequet=istat

      return
      end
