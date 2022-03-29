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

!**s/r out_thistep

      subroutine out_thistep ( F_sorties, F_step, ns, F_component_S )
      use levels
      use outd
      use outc
      use outp
      use timestep
      implicit none
#include <arch_specific.hf>

      character(len=*) F_component_S
      integer F_step,ns
      integer F_sorties(0:ns)


      integer i,j,k,liste_m(2,50),liste_p(2,50),liste_h(2,50), &
              cnt,cnt_m,cnt_p,cnt_h,o_sets
      integer, dimension(:), pointer :: o_step, o_lev, o_grid
!
!     ---------------------------------------------------------------
!
      o_sets = -1
      if (trim(F_component_S) == 'DYN' ) then
         o_sets =  Outd_sets
         o_step => Outd_step
         o_lev  => Outd_lev
         o_grid => Outd_grid
      end if
      if (trim(F_component_S) == 'PHY' ) then
         o_sets =  Outp_sets
         o_step => Outp_step
         o_lev  => Outp_lev
         o_grid => Outp_grid
      end if
      if (trim(F_component_S) == 'CHM' ) then
         o_sets =  Outc_sets
         o_step => Outc_step
         o_lev  => Outc_lev
         o_grid => Outc_grid
      end if

      cnt=0 ; cnt_m=0 ; cnt_p=0 ; cnt_h=0


      do j=1,Timestep_sets
          do i=1,Timestep_max(j)
            if (F_step == Timestep_tbl(i,j)) then
               do k=1, o_sets
                  if ( o_step(k) == j ) then
                     if (trim(F_component_S) == 'PHY' ) &
                     Outp_lasstep(k,F_step)= Timestep_tbl(max(1,i-1),j)
                     if ( Level_typ_S(o_lev(k)) == 'M') then
                        cnt_m= cnt_m+1
                        cnt  = cnt  +1
                        liste_m(1,cnt_m) = k
                        liste_m(2,cnt_m) = o_grid(k)
                     end if
                     if ( Level_typ_S(o_lev(k)) == 'P') then
                        cnt_p= cnt_p+1
                        cnt  = cnt  +1
                        liste_p(1,cnt_p) = k
                        liste_p(2,cnt_p) = o_grid(k)
                     end if
                     if ( Level_typ_S(o_lev(k)) == 'H') then
                        cnt_h= cnt_h+1
                        cnt  = cnt  +1
                        liste_h(1,cnt_h) = o_grid(k)
                        liste_h(2,cnt_h) = k
                     end if
                  end if
               end do
            end if
         end do
      end do

      F_sorties(0) = cnt

      call PIKSRT (liste_m,2,cnt_m)
      call PIKSRT (liste_p,2,cnt_p)
      call PIKSRT (liste_h,2,cnt_h)

      F_sorties(            1:cnt_m)            = liste_m(1,1:cnt_m)
      F_sorties(cnt_m      +1:cnt_m+cnt_p)      = liste_p(1,1:cnt_p)
      F_sorties(cnt_m+cnt_p+1:cnt_m+cnt_p+cnt_h)= liste_h(1,1:cnt_h)
!
!     ---------------------------------------------------------------
!
      return
      end

!*****************************************************
!* Sorts an array ARR of length N in ascending order *
!* by straight insertion.                            *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table ARR                  *
!*          ARR	  table to be sorted                 *
!* OUTPUT:                                           *
!*	    ARR   table sorted in ascending order    *
!*                                                   *
!* NOTE: Straight insertion is a N² routine and      *
!*       should only be used for relatively small    *
!*       arrays (N<100).                             *
!*****************************************************
      subroutine piksrt(arr,n1,n2)
      use timestep
      implicit none
      integer n1,n2
      integer ARR(n1,n2)

      integer i,j,a(n1)
      do j=2, N2
         a(1:n1) = ARR(1:n1,j)
         do i=j-1,1,-1
            if (ARR(n1,i)<=a(n1)) goto 10
            ARR(:,i+1)=ARR(:,i)
         end do
         i=0
10       ARR(1:n1,i+1)=a(1:n1)
      end do
      return
      end subroutine piksrt
