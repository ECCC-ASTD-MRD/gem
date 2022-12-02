!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer, 
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms 
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer 
!version 3 or (at your option) any later version that should be found at: 
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html 
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software; 
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec), 
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------
#include <rmn/msg.h>

!/@*
subroutine test_sort_mod()
   use testutils
   use sort_mod
   implicit none
   !@author S.Chamberland, 2012-01
   !*@/
#include <rmnlib_basics.hf>
   logical :: ok_L
   integer :: itable5(5),itable6(6), &
        itable0(5)   = (/5,2,4,1,3/),&
        itable0a(5)  = (/1,2,3,4,5/),&
        itable0ar(5) = (/5,4,3,2,1/),&
        itable0au(5) = (/1,2,3,4,5/),&
        itable1(6)   = (/5,2,2,4,1,3/),&
        itable1a(6)  = (/1,2,2,3,4,5/),&
        itable1ar(6) = (/5,4,3,2,2,1/),&
        itable1au(5) = (/1,2,3,4,5/)
   real :: rtable5(5), rtable6(6), &
        rtable0(5)   = (/5.,2.,4.,1.,3./),&
        rtable0a(5)  = (/1.,2.,3.,4.,5./),&
        rtable0ar(5) = (/5.,4.,3.,2.,1./),&
        rtable0au(5) = (/1.,2.,3.,4.,5./),&
        rtable1(6)   = (/5.,2.,2.,4.,1.,3./),&
        rtable1a(6)  = (/1.,2.,2.,3.,4.,5./),&
        rtable1ar(6) = (/5.,4.,3.,2.,2.,1./),&
        rtable1au(5) = (/1.,2.,3.,4.,5./)
   character(len=16) :: stable5(5), stable6(6), &
        stable0(5)   = (/'5','2','4','1','3'/),&
        stable0a(5)  = (/'1','2','3','4','5'/),&
        stable0ar(5) = (/'5','4','3','2','1'/),&
        stable0au(5) = (/'1','2','3','4','5'/),&
        stable1(6)   = (/'5','2','2','4','1','3'/),&
        stable1a(6)  = (/'1','2','2','3','4','5'/),&
        stable1ar(6) = (/'5','4','3','2','2','1'/),&
        stable1au(5) = (/'1','2','3','4','5'/)
   integer :: dim,mylen
   !---------------------------------------------------------------------
   call testutils_verbosity()

   itable5 = itable0
   dim = sort(itable5)
   ok_L = (all(itable5==itable0a))
   call testutils_assert_ok(ok_L,'test_sort_mod','sort i')

   itable5 = itable0
   dim = sort(itable5,order_L=SORT_DOWN)
   ok_L = (all(itable5==itable0ar))
   call testutils_assert_ok(ok_L,'test_sort_mod','sort ir')

   itable5 = itable0
   dim = sort(itable5,unique_L=.true.)
   ok_L = (all(itable5(1:dim)==itable0au))
   call testutils_assert_ok(ok_L,'test_sort_mod','sort iu')


   itable6 = itable1
   dim = sort(itable6)
   ok_L = (all(itable6==itable1a))
   call testutils_assert_ok(ok_L,'test_sort_mod','sort i 2')

   itable6 = itable1
   dim = sort(itable6,order_L=SORT_DOWN)
   ok_L = (all(itable6==itable1ar))
   call testutils_assert_ok(ok_L,'test_sort_mod','sort ir 2')

   itable6 = itable1
   dim = sort(itable6,unique_L=.true.)
   ok_L = (all(itable6(1:dim)==itable1au) .and. dim==5)
   call testutils_assert_ok(ok_L,'test_sort_mod','sort iu 2')


   rtable5 = rtable0
   dim = sort(rtable5)
   ok_L = (all(rtable5==rtable0a))
   call testutils_assert_ok(ok_L,'test_sort_mod','sort r')

   rtable5 = rtable0
   dim = sort(rtable5,order_L=SORT_DOWN)
   ok_L = (all(rtable5==rtable0ar))
   call testutils_assert_ok(ok_L,'test_sort_mod','sort rr')

   rtable5 = rtable0
   dim = sort(rtable5,unique_L=.true.)
   ok_L = (all(rtable5(1:dim)==rtable0au))
   call testutils_assert_ok(ok_L,'test_sort_mod','sort ru')


   rtable6 = rtable1
   dim = sort(rtable6)
   ok_L = (all(rtable6==rtable1a))
   call testutils_assert_ok(ok_L,'test_sort_mod','sort r 2')

   rtable6 = rtable1
   dim = sort(rtable6,order_L=SORT_DOWN)
   ok_L = (all(rtable6==rtable1ar))
   call testutils_assert_ok(ok_L,'test_sort_mod','sort rr 2')

   rtable6 = rtable1
   dim = sort(rtable6,unique_L=.true.)
   ok_L = (all(rtable6(1:dim)==rtable1au) .and. dim==5)
   call testutils_assert_ok(ok_L,'test_sort_mod','sort ru 2')


   stable5 = stable0
   dim = sort(stable5)
   ok_L = (all(stable5==stable0a))
   call testutils_assert_ok(ok_L,'test_sort_mod','sort s')

   stable5 = stable0
   dim = sort(stable5,order_L=SORT_DOWN)
   ok_L = (all(stable5==stable0ar))
   call testutils_assert_ok(ok_L,'test_sort_mod','sort sr')

   stable5 = stable0
   dim = sort(stable5,unique_L=.true.)
   ok_L = (all(stable5(1:dim)==stable0au))
   call testutils_assert_ok(ok_L,'test_sort_mod','sort su')


   stable6 = stable1
   dim = sort(stable6)
   ok_L = (all(stable6==stable1a))
   call testutils_assert_ok(ok_L,'test_sort_mod','sort s 2')

   stable6 = stable1
   dim = sort(stable6,order_L=SORT_DOWN)
   ok_L = (all(stable6==stable1ar))
   call testutils_assert_ok(ok_L,'test_sort_mod','sort sr 2')

   stable6 = stable1
   dim = sort(stable6,unique_L=.true.)
   ok_L = (all(stable6(1:dim)==stable1au) .and. dim==5)
   call testutils_assert_ok(ok_L,'test_sort_mod','sort su 2')
   !---------------------------------------------------------------------
   return
end subroutine test_sort_mod
