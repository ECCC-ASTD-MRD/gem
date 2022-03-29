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

!**s/r memusage - Report on memory usage using get_max_rss

      subroutine memusage ( unf )
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer unf

      integer, external :: get_max_rss
      integer i, err, imin,imax
      integer memory(Ptopo_numproc),memoryt(Ptopo_numproc)
      real(kind=REAL64) sum,sumd2,moy,var,fijk,mind,maxd,npt_8
!
!----------------------------------------------------------------
!
      memory = 0
      memory(Ptopo_myproc+1) = get_max_rss()
      call rpn_comm_REDUCE ( memory, memoryt, Ptopo_numproc,&
                       "MPI_INTEGER","MPI_SUM",0,"grid",err )
      if (unf > 0) then
         sum   = 0.0
         sumd2 = 0.0
         imin  = 1
         imax  = Ptopo_numproc
         maxd  = memoryt(Ptopo_numproc)/1000000.
         mind  = memoryt(1)/1000000.

         do i=1,Ptopo_numproc
            fijk  = memoryt(i)/1000000.
            sum   = sum   + fijk
            sumd2 = sumd2 + fijk*fijk
            if (fijk > maxd) then
               maxd = fijk
               imax = i
            end if
            if (fijk < mind) then
               mind = fijk
               imin = i
            end if
         end do

         npt_8 = Ptopo_numproc
         moy = sum / npt_8
         var = max(0.d0,1.0d0*(sumd2 + moy*moy*npt_8 - 2*moy*sum) / npt_8)
         var = sqrt(var)

         write (unf,99) real(sum),Ptopo_numproc
         write (unf,98) moy,var,imin-1,mind,imax-1,maxd

         do i=1,Ptopo_numproc
            fijk  = memoryt(i)/1000000.
            if (fijk > (moy+var)) write (unf,97) i-1,fijk
         end do
      end if

 99   format ('Memory==> Total:',f14.3,' Gb on ',i5,'  MPI processes')
 98   format ('Memory==> Mean:',f10.3,' Std:',f10.3,' Min:[(',i5,')', &
               f10.3,']',' Max:[(',i5,')',f10.3,']')
 97   format ('Memory==> largest MPI process: #',i5,f12.3,' Gb')
!
!----------------------------------------------------------------
!
      return
      end
