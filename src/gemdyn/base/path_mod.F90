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

!> Variables associated with the paths/directories (set_world_view)
module path
   implicit none
   public
   save

   !> Path to namelist directory
   character(len=1024) :: Path_nml_S
   !> Path to base directory
   character(len=1024) :: Path_basedir_S
   !> Path to input directory
   character(len=1024) :: Path_input_S
   !> Path to indata directory
   character(len=1024) :: Path_ind_S
   !> Path to output directory
   character(len=1024) :: Path_output_S
   !> Path to work directory
   character(len=1024) :: Path_work_S
   !> Path to the "exchange" directory with 3D-Var
   character(len=1024) :: Path_xchg_S
   !> Path and name of the file "output_settings"
   character(len=1024) :: Path_outcfg_S
   !> Path to physics files
   character(len=1024) :: Path_phy_S
   !> Path to physics_input_table
   character(len=1024) :: Path_phyincfg_S
end module path
