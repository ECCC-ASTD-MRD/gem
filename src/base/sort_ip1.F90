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

subroutine sort_ip1 ( F_keys, F_ip1, F_nka )
   implicit none

   integer F_nka
   integer F_keys(F_nka), F_ip1 (F_nka)

#include <rmnlib_basics.hf>

   character(len=1)  :: tva, grda, blk_S
   character(len=4)  :: var
   character(len=12) :: etik_S
   integer i,j,k,ni,nj,nk,cnt,err,diag_ip1,diag_key
   integer dateo, deet, ipas, ip1a, ip2a, ip3a, &
        ig1a, ig2a, ig3a, ig4a, bit, datyp , &
        swa, lng, dlf, ubc, ex1, ex2, ex3, kind, res
   real ip1r4(F_nka),lev,temp_ip1(F_nka),hugetype
   integer liste(F_nka),keys(F_nka),ip1(F_nka)
   !
   ! --------------------------------------------------------------------
   !
   cnt= 0 ; diag_ip1= 0 ; kind= -1

   do k= 1, F_nka
      err = fstprm ( F_keys(k), dateo, deet, ipas, ni, nj, nk,&
           bit, datyp, ip1a,ip2a,ip3a, tva, var, etik_S, grda,&
           ig1a,ig2a,ig3a,ig4a, swa,lng, dlf, ubc, ex1, ex2, ex3 )
      call convip_plus(ip1a, lev, i,-1, blk_S, .false.)
      if ( i /= 4 ) then
         cnt= cnt + 1
         ip1r4(cnt)= lev
         keys (cnt)= F_keys(k)
         ip1  (cnt)= ip1a
         kind = i
      else
         diag_key  = F_keys(k)
         diag_ip1  = ip1a
      endif
   end do

   temp_ip1(1:cnt) = ip1r4(1:cnt)
   if (kind == 21) then
      ! Sort levels in descending order
      do i=1,cnt
         res= maxloc(temp_ip1(1:cnt),1)
         liste(i)= res
         temp_ip1(liste(i)) = -huge(hugetype)
      end do
   else
      ! Sort levels in ascending order
      do i=1,cnt
         res= minloc(temp_ip1(1:cnt),1)
         liste(i)= res
         temp_ip1(liste(i)) = huge(hugetype)
      end do
   endif

   do i=1, cnt
      F_ip1 (i)= ip1  (liste(i))
      F_keys(i)= keys (liste(i))
   end do

   !     Eliminate levels that are redundant in LISTE
   i= 1
   do j =2, cnt
      if (F_ip1(i) /= F_ip1(j)) then
         i = i+1
         if (i /= j) then
            F_ip1  (i) = F_ip1  (j)
            F_keys (i) = F_keys (j)
         endif
      endif
   enddo
   F_nka = i

   !     Take care of diagnostic level
   if (diag_ip1 > 0) then
      F_nka= F_nka+1
      F_ip1 (F_nka)= diag_ip1
      F_keys(F_nka)= diag_key
   endif
   !
   ! --------------------------------------------------------------------
   !
   return
end subroutine sort_ip1
