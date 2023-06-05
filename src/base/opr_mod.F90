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

module opr
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save
!______________________________________________________________________
!                                                                      |
!  PROJECTION OPERATORS FOR THE SOLVER (initialized in set_opr)        |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! Opr_opsxp0_8       | east-west   projection operators                |
! Opr_opsxp2_8       | east-west   projection operators                |
! Opr_opsyp0_8       | north-south projection operators                |
! Opr_opsyp2_8       | north-south projection operators                |
! Opr_opszp0_8       | vertical    projection operator                 |
! Opr_opszpm_8       | vertical          mean operator                 |
! Opr_opszpl_8       | vertical    towards ln operator                 |
! Opr_opszp2_8       | vertical    second der operator                 |
! Opr_xeval_8        | horizontal eigenvalues                          |
! Opr_zevec_8        | right vertical eigenvectors                     |
! Opr_lzevec_8       | left vertical eigenvectors                      !
! Opr_zeval_8        | vertical   eigenvalues                          |
! Opr_evvec_8        | even eigenvectors                               |
! Opr_odvec_8        | odd  eigenvectors                               |
!----------------------------------------------------------------------
!
   real(kind=REAL64), dimension(:), allocatable :: Opr_opsxp0_8, Opr_opsyp0_8, Opr_evvec_8 , &
                                        Opr_odvec_8 , Opr_opsxp2_8, Opr_opsyp2_8, &
                                        Opr_opszp0_8, Opr_opszpm_8, Opr_opszpl_8, &
                                        Opr_opszp2_8, Opr_xeval_8 , &
                                        Opr_zevec_8, Opr_lzevec_8, Opr_zeval_8
end module opr
