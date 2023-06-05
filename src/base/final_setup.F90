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

!**s/r final_setup
!
      subroutine final_setup ()
      use gmm_table
      use lun
      implicit none

      integer i
!     
!     ---------------------------------------------------------------
!
      gmm_ncles= gmm_nkeys()
      allocate (gmm_keylist(gmm_ncles))
      gmm_ncles= gmm_keys(gmm_keylist)
      if (Lun_debug_L) then
         write(Lun_out,900)
         write(Lun_out,1006)
         write(Lun_out,901)
         do i=1, gmm_cnt
            write(Lun_out,1007) GMM_tbl%vname(i),GMM_tbl%fst(i),GMM_tbl%ara(i),GMM_tbl%cn(i)
         end do
         write(Lun_out,1006)
      endif

  900 format(/'+',57('-'),'+',18('-'),'+',/'|    GMM VARIABLES AVAILABLE FOR OUTPUT',19x,'|',18x,'+')
  901 format('|',10x,'gmm_name',18x,'|',' FST_name | Arakawa | Charney-Phillips |')
 1006 format('+',36('-'),'+',10('-'),'+',9('-'),'+',18('-'),'+')
 1007 format('|  ',a32,'  |   ',a4,3x,'|   ',a2,4x,'|       ',a2,9x,'|')
!     
!     ---------------------------------------------------------------
!
      return
      end subroutine final_setup
