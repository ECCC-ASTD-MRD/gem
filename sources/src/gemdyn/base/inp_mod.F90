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
module inp_mod
   implicit none
   public
   save

   character(len=16) :: Inp_datev
   logical Inp_zd_L,Inp_w_L
   integer Inp_nfiles , Inp_comm_id, Inp_comm_setno,&
           Inp_iome   , Inp_comm_io, Inp_iobcast, Inp_kind      ,&
              Inp_version, Inp_handle , Inp_cmcdate
   integer, dimension(:), pointer :: Inp_list_unf => null()
   real*8 Inp_pref_a_8

   ! Remove the following 2 lines by 2021
   integer Inp_ut1_is_urt1
end module inp_mod
