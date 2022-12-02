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
module gmm_table
      use rmn_gmm
      implicit none
      public
      save

      integer, parameter :: MAXNAMELENGTH= 32 , gmm_max_elem= 500
      integer :: gmm_ncles
      integer :: gmm_cnt= 0
      character(len=GMM_MAXNAMELENGTH), dimension(:), pointer :: gmm_keylist

      type :: gmm_AraCN
         sequence
         character(len=MAXNAMELENGTH), dimension(:), pointer :: vname,ara,cn,fst
      end type gmm_AraCN
      type(gmm_AraCN) :: GMM_tbl

end module gmm_table
