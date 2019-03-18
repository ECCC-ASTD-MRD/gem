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

module stat_mpi
  implicit none
#include <arch_specific.hf>
  private
  public :: statf_dm

  interface statf_dm
     module procedure statf_r4
     module procedure statf_r8
  end interface

!object
!        calcule la moyenne, la variance, le minimum et
!        le maximum d un champs et imprime le resultat SANS GLBCOLC
!        Le calcule peut changer dependement le topologie

!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_field       I         Field to be operated on
! F_nv_S        I         User provided string to define F_field
! F_no          I         Usually the timestep #
! F_from_S      I         Usually the name of the calling subroutine
! F_i0,F_j0     I         Global lower-left indexes of the sub-domain
!                            on which to perform statistics
! F_in,F_jn     I         Global upper-right indexes of the sub-domain
!                            on which to perform statistics
! F_k0,F_kn     I         Range of levels on which to perform statistics
!----------------------------------------------------------------
!

contains

      subroutine statf_r4 ( F_field, F_nv_S, F_no, F_from_S, &
                            Minx,Maxx,Miny,Maxy,Mink,Maxk  , &
                            F_i0,F_j0,F_k0,F_in,F_jn,F_kn,F_rx)

      use gem_options
      use glb_ld
      use lun
      use ptopo
      implicit none

      character(len=*) F_nv_S , F_from_S
      integer Minx,Maxx,Miny,Maxy,Mink,Maxk, &
              F_i0,F_j0,F_k0,F_in,F_jn,F_kn,F_no,F_rx
      real F_field (Minx:Maxx,Miny:Maxy,Mink:Maxk)

      include 'statf_dm.inc'
!
!----------------------------------------------------------------
!
      return
      end subroutine statf_r4

      subroutine statf_r8 ( F_field, F_nv_S, F_no, F_from_S, &
                            Minx,Maxx,Miny,Maxy,Mink,Maxk  , &
                            F_i0,F_j0,F_k0,F_in,F_jn,F_kn,F_rx)

      use gem_options
      use glb_ld
      use lun
      use ptopo
      implicit none

      character(len=*) F_nv_S , F_from_S
      integer Minx,Maxx,Miny,Maxy,Mink,Maxk, &
              F_i0,F_j0,F_k0,F_in,F_jn,F_kn,F_no,F_rx
      real*8 F_field (Minx:Maxx,Miny:Maxy,Mink:Maxk)

      include 'statf_dm.inc'
!
!----------------------------------------------------------------
!
      return
      end subroutine statf_r8

end module stat_mpi

