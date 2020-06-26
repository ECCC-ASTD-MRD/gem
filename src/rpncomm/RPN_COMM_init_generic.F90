!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 2020  Division de Recherche en Prevision Numerique
! !                     Environnement Canada
! !
! ! This library is free software; you can redistribute it and/or
! ! modify it under the terms of the GNU Lesser General Public
! ! License as published by the Free Software Foundation,
! ! version 2.1 of the License.
! !
! ! This library is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! ! Lesser General Public License for more details.
! !
! ! You should have received a copy of the GNU Lesser General Public
! ! License along with this library; if not, write to the
! ! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! ! Boston, MA 02111-1307, USA.
! !/
SUBROUTINE rpn_comm_dummy_init(Pex,Pey)
  implicit none
  integer, intent(inout) :: Pex,Pey
  Pex = Pex
  Pey = Pey
  return
end SUBROUTINE rpn_comm_dummy_init

!! end interface                         !InTf!
!! interface RPN_COMM_gridinit           !InTf!
FUNCTION RPN_COMM_init1(Userinit,Pelocal,Petotal) result(grid)   !InTf!
  implicit none
  external Userinit                                              !InTf!
  integer, intent(out)   :: Pelocal,Petotal                      !InTf!
  integer :: grid                                                !InTf!
  integer :: Pex,Pey
  Pex = 0
  Pey = 0
  call RPN_COMM_init(Userinit,Pelocal,Petotal,Pex,Pey)
  grid = 0
  return
end FUNCTION RPN_COMM_init1                                      !InTf!

FUNCTION RPN_COMM_init2(Pelocal,Petotal,Pex,Pey) result(grid)    !InTf!
  implicit none
  integer, intent(out)   :: Pelocal,Petotal                      !InTf!
  integer, intent(inout) :: Pex,Pey                              !InTf!
  integer :: grid                                                !InTf!
  external rpn_comm_dummy_init
  call RPN_COMM_init(rpn_comm_dummy_init,Pelocal,Petotal,Pex,Pey)
  grid = 0
  return
end FUNCTION RPN_COMM_init2                                      !InTf!

FUNCTION RPN_COMM_init_multigrid1(Userinit,Pelocal,Petotal,MultiGrids) result(grid)  !InTf!
  implicit none
  external Userinit                                                                  !InTf!
  integer, intent(out)   :: Pelocal,Petotal                                          !InTf!
  integer, intent(in)    :: MultiGrids                                               !InTf!
  integer :: grid                                                                    !InTf!
  integer :: Pex,Pey
  integer, external :: RPN_COMM_init_multigrid
  Pex = 0
  Pey = 0
  grid = RPN_COMM_init_multigrid(Userinit,Pelocal,Petotal,Pex,Pey,MultiGrids)
end FUNCTION RPN_COMM_init_multigrid1                                                !InTf!

FUNCTION RPN_COMM_init_multigrid2(Pelocal,Petotal,Pex,Pey,MultiGrids) result(grid)   !InTf!
  implicit none
  integer, intent(out)   :: Pelocal,Petotal                                          !InTf!
  integer, intent(inout) :: Pex,Pey                                                  !InTf!
  integer, intent(in)    :: MultiGrids                                               !InTf!
  integer :: grid                                                                    !InTf!
  external rpn_comm_dummy_init
  integer, external :: RPN_COMM_init_multigrid
  grid = RPN_COMM_init_multigrid(rpn_comm_dummy_init,Pelocal,Petotal,Pex,Pey,MultiGrids)
end FUNCTION RPN_COMM_init_multigrid2                                                !InTf!

FUNCTION RPN_COMM_init_multi_level1(Userinit,Pelocal,Petotal,MultiGrids,Grids) result(grid)  !InTf!
  implicit none
  external Userinit                                                                          !InTf!
  integer, intent(out)   :: Pelocal,Petotal                                                  !InTf!
  integer, intent(in)    :: MultiGrids                                                       !InTf!
  integer, intent(in)    :: Grids                                                            !InTf!
  integer :: grid                                                                            !InTf!
  integer :: Pex,Pey
  integer, external :: RPN_COMM_init_multi_level
  Pex = 0
  Pey = 0
  grid = RPN_COMM_init_multi_level(Userinit,Pelocal,Petotal,Pex,Pey,MultiGrids,Grids)
end FUNCTION RPN_COMM_init_multi_level1                                                      !InTf!

FUNCTION RPN_COMM_init_multi_level2(Pelocal,Petotal,Pex,Pey,MultiGrids,Grids) result(grid)   !InTf!
  implicit none
  integer, intent(out)   :: Pelocal,Petotal                                                  !InTf!
  integer, intent(inout) :: Pex,Pey                                                          !InTf!
  integer, intent(in)    :: MultiGrids                                                       !InTf!
  integer, intent(in)    :: Grids                                                            !InTf!
  integer :: grid                                                                            !InTf!
  external rpn_comm_dummy_init
  integer, external :: RPN_COMM_init_multi_level
  grid = RPN_COMM_init_multi_level(rpn_comm_dummy_init,Pelocal,Petotal,Pex,Pey,MultiGrids,Grids)
end FUNCTION RPN_COMM_init_multi_level2                                                      !InTf!
!! end interface RPN_COMM_gridinit       !InTf!
!! interface                             !InTf!

