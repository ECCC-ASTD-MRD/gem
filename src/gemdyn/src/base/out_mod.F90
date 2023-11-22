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

module out_mod
   implicit none
   public
   save

!
!**comdeck out.cdk
!
!______________________________________________________________________
!                                                                      |
!  OUTPUT GLOBAL CONTROL VARIABLES                                     |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! Out_etik_S         | label, etiket                                   |
! Out_runstrt_S      | date and time of start of pilot run             |
! Out_endstepno      | last timestep to be included in current output  |
! Out_typvar_S       | type of variable 'P'-predicted, 'A'-analysis    |
! Out_deet           | timestep interval in seconds Cstv_dt_8          |
! Out_dateo          | date and time of origin                         |
! Out_ixg(4)         | grid descriptors                                |
! Out_ip2            | value for IP2                                   |
! Out_ip3            | value for IP3                                   |
! Out_gridtyp_S      | grid type 'E' or '#' or 'M' or 'G'              |
! Out_gridi0         | starting I of entire grid (OutGrid_x0)             |
! Out_gridin         | ending   I of entire grid (OutGrid_x1)             |
! Out_gridj0         | starting J of entire grid (OutGrid_y0)             |
! Out_gridjn         | ending   J of entire grid (OutGrid_y1)             |
! Out_rewrit_L       | .true. to overwrite, .false. to not overwrite   |
! ---------------------------------------------------------------------

   character(len=1024) :: Out_dirname_S
   character(len=19)   :: Out_laststep_S
   character(len=12)   :: Out_etik_S
   character(len=6)    :: Out_type_S
   character(len=1)    :: Out_gridtyp_S, Out_proj_S,Out_unit_S
   character(len=2)    :: Out_typvar_S,Out_prefix_S
   character(len=3)    :: Out_ext_S
   logical :: Out_rewrit_L
   logical :: Out_post_L, Out_diruse_L, Out_reduc_l
   real    :: Out_rot(8)

   integer Out_unf,Out_date,Out_hour,Out_dateo
   integer Out_deet,Out_npas,Out_endstepno
   integer Out_gridi0,Out_gridin
   integer Out_gridj0,Out_gridjn
   integer Out_stride
   integer Out_ig1, Out_ig2, Out_ig3, Out_ig4
   integer Out_ixg(8), Out_ip2, Out_ip3
   integer Out_periodx,Out_periody

end module out_mod
