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

!**s/r liebman_2 - modifiy a portion of input field by overrelaxation
!
      subroutine liebman_2 ( F_field,  F_mask, &
                             F_fielde, F_maske, F_max, F_ni, F_nj)
!
      implicit none
#include <arch_specific.hf>
!
      integer F_ni, F_nj
      real    F_max
      real    F_field  (F_ni,F_nj  ), F_mask  (F_ni,F_nj  )
      real    F_fielde (F_ni,F_nj+2), F_maske (F_ni,F_nj+2)
!
!author
!     Alain Patoine - after version v1_03 of liebman.ftn
!
!revision
! v2_00 - Desgagne M.       - initial MPI version (from liebman v1_03)
!
!object
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
!  F_field     I/O
!----------------------------------------------------------------
! _____________________________________________________________________
!         |                                             |           |  |
!  NAME   |             DESCRIPTION                     |DIMENSIONS |IN|
!         |                                             |           |OU|
! --------|---------------------------------------------|-----------|--|
! F_field | field to be treated                         | F_ni, F_nj|io|
! --------|---------------------------------------------|-----------|--|
! F_mask  | mask: 0 -> don't modify                     | F_ni, F_nj| i|
!         |       1 ->       modify                     |           |  |
! --------|---------------------------------------------|-----------|--|
! F_fielde| field to be treated with poles extensions   | F_ni,     | w|
!         | (work field)                                |     Fnj+2 |  |
! --------|---------------------------------------------|-----------|--|
! F_maskde| mask with poles extensions                  | F_ni,     | w|
!         | (work field)                                |     Fnj+2 |  |
! --------|---------------------------------------------|-----------|--|
! F_max   | convergence criteria                        | scalar    | i|
! --------|---------------------------------------------|-----------|--|
! F_ni    | number of points in x-direction             | scalar    | i|
! F_nj    | number of points in y-direction             | scalar    | i|
! ----------------------------------------------------------------------
!
!
      integer i, j, pnj, n, pnitmax, pnl, pnr
!
      real prfact, prmod, prmax, prvals, prvaln, prmass, prmasn
!*
!*********************************************************************
!
! WE INITIALISE A FEW THINGS
! --------------------------
!
! prfact = overrelaxation coefficient / 4.
!          overrelaxation coefficient must be between 1. and 2.
!
! pnitmax = maximum number of iterations
!
!*********************************************************************
      prfact  = 1.75 * 0.25
      pnitmax = 100
!*********************************************************************
!
!      -->  Put the average of the last position
!           in the north pole position ........... (F_nj+2)
!
!      -->  Rearrange rows of data     (F_nj)   -> (F_nj+1)
!                                      (F_nj-1) -> (F_nj)
!                                       .           .
!                                       .           .
!                                      (1)      -> (2)
!
!           Put the average of the first position
!           in the south pole position ........... (1)
!
!*********************************************************************
      pnj = F_nj+2
!
      prvals = 0.0
      prvaln = 0.0
      prmass = 0.0
      prmasn = 0.0
!
      do i=1,F_ni
         prvals = prvals + F_field(i,1   )
         prvaln = prvaln + F_field(i,F_nj)
         prmass = amax1 ( prmass, F_mask(i,1   ) )
         prmasn = amax1 ( prmasn, F_mask(i,F_nj) )
      end do

      prvals = prvals / F_ni
      prvaln = prvaln / F_ni

      do i=1,F_ni
         F_fielde (i,1     ) = prvals
         F_maske  (i,1     ) = prmass
         F_fielde (i,F_nj+2) = prvaln
         F_maske  (i,F_nj+2) = prmasn
      end do
!
      do j=1,F_nj
      do i=1,F_ni
         F_fielde(i,j+1) = F_field(i,j)
         F_maske (i,j+1) = F_mask (i,j)
      end do
      end do
!*********************************************************************
! Begin iterations                                                   *
!*********************************************************************
      do 100 n=1,pnitmax
!
         prmax = 0.0
!     ****************************************************************
!     * South pole                                                   *
!     ****************************************************************
         if ( F_maske(1,1) > 0.5 ) then
!
            prmod = 0.0
!
            do i=1,F_ni
               prmod = prmod + F_fielde(i,2)
            end do
!
            prmod =  prfact * ( prmod - F_ni * F_fielde(1,1) )
            prmod =  prmod * 4.0 / F_ni
!
            prmax = amax1 ( prmax, abs(prmod) )
!
            do i=1,F_ni
               F_fielde(i,1) = F_fielde(i,1) + prmod
            end do
!
         end if
!     *****************************************************************
!     * Interior of domain                                            *
!     *****************************************************************
         do j=2,pnj-1
         do i=1,F_ni
            pnl = i-1
            pnr = i+1
!
            if ( i == 1    ) pnl = F_ni
            if ( i == F_ni ) pnr = 1
!
            if ( F_maske(i,j) > 0.5 ) then
               prmod = prfact * (F_fielde(pnl,j) + F_fielde(pnr,j) + &
                                 F_fielde(i,j-1) + F_fielde(i,j+1) - &
                              4.*F_fielde(i,j))
!
               prmax = amax1 ( prmax, abs(prmod) )
!
               F_fielde(i,j) = F_fielde(i,j) + prmod
            end if
         end do
         end do
!     ****************************************************************
!     * North pole                                                   *
!     ****************************************************************
         if ( F_maske(1,pnj) > 0.5 ) then
!
            prmod = 0.0
!
            do i=1,F_ni
               prmod = prmod + F_fielde(i,pnj-1)
            end do
!
            prmod =  prfact * ( prmod - F_ni * F_fielde(1,pnj) )
            prmod =  prmod * 4.0 / F_ni
!
            prmax = amax1 ( prmax, abs(prmod) )
!
            do i=1,F_ni
               F_fielde(i,pnj) = F_fielde(i,pnj) + prmod
            end do
!
         end if
!*********************************************************************
         if ( prmax < F_max ) go to 200
!
 100  continue
 200  continue
!
      do j=1,F_nj
      do i=1,F_ni
         F_field (i,j) = F_fielde (i,j+1)
      end do
      end do
!
      return
      end
