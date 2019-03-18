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
#include <arch_specific.hf>

   integer F_nka
   integer F_keys(F_nka), F_ip1 (F_nka)

   !author 
   !     Michel Desgagne  -  summer 2015
   !revision
   ! v4_8 - Desgagne M.       - initial version
   !
#include <rmnlib_basics.hf>

   character*1  tva, grda, blk_S
   character*4  var
   character*12 etik_S
   integer i,j,k,ni,nj,nk,cnt,err,diag_ip1,diag_key,encode_style
   integer dateo, deet, ipas, ip1a, ip2a, ip3a, &
        ig1a, ig2a, ig3a, ig4a, bit, datyp , &
        swa, lng, dlf, ubc, ex1, ex2, ex3, kind
   real ip1r4(F_nka),diag_ip1r4,lev,temp_ip1(F_nka),hugetype
   integer liste(F_nka),res(F_nka),keys(F_nka),ip1(F_nka)
   !
   ! --------------------------------------------------------------------
   !
   cnt= 0 ; diag_ip1= 0
! ; encode_style= 3
!!$   err = fstprm ( F_keys(1), dateo, deet, ipas, ni, nj, nk ,&
!!$         bit, datyp, ip1a,ip2a,ip3a, tva, var, etik_S, grda,&
!!$         ig1a,ig2a,ig3a,ig4a, swa,lng, dlf, ubc, ex1, ex2, ex3 )
!!$   call convip (ip1a, lev, i,-1, blk_S, .false.)
!!$   call convip (k   , lev, i, 1, ""   , .false.)
!!$   if (k == ip1a) encode_style= 1

   do k= 1, F_nka
      err = fstprm ( F_keys(k), dateo, deet, ipas, ni, nj, nk,&
           bit, datyp, ip1a,ip2a,ip3a, tva, var, etik_S, grda,&
           ig1a,ig2a,ig3a,ig4a, swa,lng, dlf, ubc, ex1, ex2, ex3 )
      call convip (ip1a, lev, i,-1, blk_S, .false.)
      if ( i .ne. 4 ) then
         cnt= cnt + 1
         ip1r4(cnt)= lev
         keys (cnt)= F_keys(k)
         ip1  (cnt)= ip1a
!         kind = i
      else
         diag_key  = F_keys(k)
         diag_ip1  = ip1a
!         diag_ip1r4= lev
      endif
   end do

   ! Sort levels in ascending order

   temp_ip1(1:cnt) = ip1r4(1:cnt)
   do i=1,cnt
      res= minloc(temp_ip1(1:cnt))
      liste(i)= res(1)
      temp_ip1(liste(i)) = huge(hugetype)
   end do
   do i=1, cnt
      F_ip1 (i)= ip1  (liste(i))
      F_keys(i)= keys (liste(i))
   end do
!   call sort (ip1r4,cnt)

!   do i= 1, cnt
!      call convip ( F_ip1(i), ip1r4(i), kind, encode_style, "", .false. )
!   end do

   !     Eliminate levels that are redundant in LISTE
   i= 1
   do j =2, cnt
      if (F_ip1(i) .ne. F_ip1(j)) then
         i = i+1
         if (i .ne. j) then
            F_ip1  (i) = F_ip1  (j)
            F_keys (i) = F_keys (j)
         endif
      endif
   enddo
   F_nka = i

   !     Take care of diagnostic level
   if (diag_ip1 .gt. 0) then
      F_nka= F_nka+1
      F_ip1 (F_nka)= diag_ip1
      F_keys(F_nka)= diag_key
   endif
   !
   ! --------------------------------------------------------------------
   !
   return
end subroutine sort_ip1
