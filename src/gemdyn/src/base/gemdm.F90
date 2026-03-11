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

!**s/r gemdm -  Main entry point for the GEM model

      subroutine gemdm (F_COMMs,F_cme,F_nc)
      use lun
      implicit none
      
      integer, intent (IN) :: F_cme,F_nc
      integer, intent (IN) :: F_COMMs(F_nc)

      integer :: un_out
!
!     ---------------------------------------------------------------
!
      un_out= -1 ; if ( F_cme == 0 ) un_out= 6
      call gemtime ( 6, ' ', .false. )
      if (un_out>0) call gemtime ( un_out, &
                       'STARTING GEMDM DOMAINS', .false. )

! Initialize: Domain, MPI, processor topology and ptopo.cdk
      call init_component(F_COMMs,F_nc)

! Establish: model configuration, domain decomposition and model geometry
      call set_world_view()
  
! Initialize the ensemble prevision system
      call itf_ens_init()
  
! Initialize the physics parameterization package
      call itf_phy_init()
  
! Initialize tracers
      call tracers()
  
! Setup main memory
      call main_gmm_storage()
      call set_dyn_opr()
  
! Run GEM
      call gem_ctrl()
  
! Terminate
      call stop_world_view()

      if (un_out>0) call gemtime ( un_out, &
            'ENDING GEMDM DOMAINS', .true. )
!
!     ---------------------------------------------------------------
!
      return
      end subroutine gemdm
