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

module path
   implicit none
   public
   save
!______________________________________________________________________
!                                                                      |
!  VARIABLES ASSOCIATED WITH THE PATHs/DIRECTORIES (set_world_view1)   |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! Path_nml_S         | path to namelist directory                      |
! Path_input_S       | path to input directory                         |
! Path_ind_S         | path to indata directory                        |
! Path_output_S      | path to output directory                        |
! Path_work_S        | path to work directory                          |
! Path_outcfg_S      | path and name of the file "output_settings"     |
! Path_phy_S         | path to physics files                           |
! Path_phyincfg_S    | path to physics_input_table                     |
!----------------------------------------------------------------------

      character(len=2048) :: Path_nml_S,Path_input_S
      character(len=2048) :: Path_ind_S,Path_output_S,Path_work_S
      character(len=2048) :: Path_outcfg_S,Path_phy_S,Path_phyincfg_S

end module path
