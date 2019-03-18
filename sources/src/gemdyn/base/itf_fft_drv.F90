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

!**s/r itf_fft_drv
!
      subroutine itf_fft_drv ( F_vec, F_n1 ,F_n2 ,F_n3, F_dir )
      use fft
      implicit none
#include <arch_specific.hf>
!
      integer F_n1 ,F_n2 ,F_n3, F_dir
      real*8  F_vec(*)
!
!author
!     Michel Desgagne - spring 2012
!
!revision
! v4_50 - Desgagne M.       - initial version


!     __________________________________________________________________
!
      if (Fft_type_S == "PERIODIC" ) call ffft8        (F_vec, F_n1, F_n2, F_n3, F_dir)
      if (Fft_type_S == "SIN"      ) call itf_fft_sin  (F_vec, F_n1, F_n2, F_n3      )
      if (Fft_type_S == "QCOS"     ) call itf_fft_qcos (F_vec, F_n1, F_n2, F_n3, F_dir)
!     __________________________________________________________________
!
      return
      end
