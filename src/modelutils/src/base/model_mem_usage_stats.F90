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


!/@*
subroutine model_mem_usage_stats(Ptopo_numproc, Ptopo_myproc)
   use iso_c_binding
   use rpn_comm_itf_mod
   implicit none
!!!#include <arch_specific.hf>
   !@objective Report on memory usage using get_max_rss
   !@arguments
   integer, intent(in) :: Ptopo_numproc, Ptopo_myproc
   !@author M. Desgagne -- Summer 2013
   !*@/
#include <rmn/msg.h>

#if defined (AIX)

   integer, external :: get_max_rss

   character(len=256) :: tmp_S
   integer :: i, err, imin, imax, msgLevelMin, msgUnit
   integer :: memory(Ptopo_numproc), memoryt(Ptopo_numproc)
   real(REAL64) ::  :: sum, sumd2, moy, var, fijk, mind, maxd, moy, var, npt_8
   !----------------------------------------------------------------
   call msg_getInfo(canWrite_L, msgLevelMin, msgUnit, tmp_S)
   canWrite_L = (Ptopo_myproc == RPN_COMM_MASTER)

   memory = 0
   memory(Ptopo_myproc+1) = get_max_rss()
   call rpn_comm_reduce(memory, memoryt, Ptopo_numproc, &
        RPN_COMM_INTEGER, RPN_COMM_SUM, RPN_COMM_MASTER, RPN_COMM_GRID, err)

   if (canWrite_L) then
      sum   = 0.0
      sumd2 = 0.0
      imin  = 1
      imax  = Ptopo_numproc
      maxd  = memoryt(Ptopo_numproc)/1000000.
      mind  = memoryt(1)/1000000.

      do i = 1, Ptopo_numproc
         fijk  = memoryt(i)/1000000.
         sum   = sum   + fijk
         sumd2 = sumd2 + fijk*fijk
         if (fijk > maxd) then
            maxd = fijk
            imax = i
         endif
         if (fijk < mind) then
            mind = fijk
            imin = i
         endif
      end do

      npt_8 = Ptopo_numproc
      moy = sum / npt_8
      var = max(0.d0, 1.d0*(sumd2 + moy*moy*npt_8 - 2.*moy*sum) / npt_8)
      var = sqrt(var)

99    format('Memory ==> Total:', f14.3, ' Gb on ', i5, '  MPI processes')
      write(msgUnit,99) real(sum), Ptopo_numproc

98    format('Memory ==> Mean:', f10.3, ' Std:', f10.3, ' Min:[(', i5, ')', &
           f10.3, ']', ' Max:[(', i5,')', f10.3, ']')
      write(msgUnit,98) moy, var, imin-1, mind, imax-1, maxd

      do i = 1, Ptopo_numproc
         fijk  = memoryt(i)/1000000.
97       format('Memory ==> largest MPI process: #', i3, f14.3, ' Gb')
         if (fijk > (moy+var)) write(msgUnit,97) i-1, fijk
      end do
!!$      write(msgUnit, '(a,f0.3)') 'Memory each process [Gb]:',memoryt/1000000.
   endif
   !----------------------------------------------------------------

#else

   if (Ptopo_myproc < 0) print *, 'model_mem_usage_stats (stub): ', Ptopo_numproc, Ptopo_myproc

#endif

   return
end subroutine model_mem_usage_stats
