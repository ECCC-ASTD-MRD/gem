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

module cpl_init_mod

   private
   public :: cpl_init

contains

      integer function cpl_init (F_path_S, F_print_L, F_unout, F_dateo, F_dt, &
                                 F_glb_ni, F_glb_nj, F_lcl_ni, F_lcl_nj     , &
                                 F_i0, F_j0, F_in, F_jn, F_nk               , &
                                 F_z0mtype, F_z0ttype, F_z0lat)
      implicit none
#include <arch_specific.hf>
#include <rmn/WhiteBoard.hf>

      character(len=*), intent(in) :: F_path_S,F_z0mtype,F_z0ttype
      logical, intent(in)          :: F_print_L
      integer, intent(in)          :: F_dateo, F_unout                , &
                                   F_glb_ni,F_glb_nj,F_lcl_ni,F_lcl_nj, &
                                   F_i0, F_j0, F_in, F_jn, F_nk
      real   , intent(in)          :: F_dt,F_z0lat(2)

!authors    Francois Roy -- spring 2015
! 
!revision
! v4_80 - Roy, F.  - initial version
! v5.0b1 - JM Belanger (dec 2017) - z0ttype

!Purpose
! Initialize coupling with ocean or other surface component

      include "cpl.cdk"

      integer ierr_ocn,ierr_wav,cplocn_init,istat

!
!     ---------------------------------------------------------------
!

      cpl_init=-1
      ierr_ocn=0
      ierr_wav=0

      if (F_print_L) then
        write(F_unout,*) 'cpl_init INFO: F_glb_ni=',F_glb_ni
        write(F_unout,*) 'cpl_init INFO: F_glb_nj=',F_glb_nj
        write(F_unout,*) 'cpl_init INFO: F_lcl_ni=',F_lcl_ni
        write(F_unout,*) 'cpl_init INFO: F_lcl_nj=',F_lcl_nj
        write(F_unout,*) 'cpl_init INFO: F_i0=',F_i0
        write(F_unout,*) 'cpl_init INFO: F_j0=',F_j0
        write(F_unout,*) 'cpl_init INFO: F_in=',F_in
        write(F_unout,*) 'cpl_init INFO: F_jn=',F_jn
      endif

      cpl_drv_gni= F_glb_ni ; cpl_drv_gnj= F_glb_nj
      cpl_drv_lni= F_lcl_ni ; cpl_drv_lnj= F_lcl_nj
      cpl_drv_i0 = F_i0     ; cpl_drv_j0 = F_j0
      cpl_drv_in = F_in     ; cpl_drv_jn = F_jn
      cpl_drv_gnk= F_nk     ; cpl_drv_delt= F_dt

      istat = WB_OK
      istat = min(wb_get('model/l_minx',cpl_minx),istat)
      istat = min(wb_get('model/l_maxx',cpl_maxx),istat)
      istat = min(wb_get('model/l_miny',cpl_miny),istat)
      istat = min(wb_get('model/l_maxy',cpl_maxy),istat)
      call handle_error_l(WB_IS_OK(istat),'cpl_init','problem getting values from whiteboard, l_minx...')

      cpl_dgflt_H = .false.
      cpl_rstn_L  = .false.
      if ( cpl_rstn_S == 'THIS_IS_A_RESTART' ) cpl_rstn_L = .true.

      if (cpl_ocn_L) &
        ierr_ocn= cplocn_init (F_path_S, F_print_L, F_unout, F_dateo,     &
                               F_z0mtype, F_z0lat, F_z0ttype)
        ! ierr_ocn=0, cplocn will be set to false, if ierr_ocn=1, cplocn set to true

!      if (cpl_wav_L)  ...
!            ierr_wav=... returns 2 if want cplwav set to true at exit

      if ( min(ierr_ocn,ierr_wav) < 0 ) return

      cpl_init=ierr_ocn+ierr_wav

!     0, 1, 2, or 3 gives all possibilities
!     0 - cplocn and cplwav set to false
!     1 - only cplocn set to true
!     2 - only cplwav set to true 
!     3 - cplocn and cplwav set to true

!     ---------------------------------------------------------------
!
      return
      end function cpl_init

end module cpl_init_mod
