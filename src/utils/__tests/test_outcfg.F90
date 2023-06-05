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
subroutine test_outcfg()
   use testutils
   use outcfg_mod
   implicit none
   !@author S.Chamberland, 2011-09
   !*@/
#include <rmnlib_basics.hf>
   integer, parameter :: MAX_ITEM = 99
   integer, parameter :: N_STEP   = 3
   integer :: my_id(3), istat, n_items,sum_items, istep,ivar,dateo,dt, &
        datev,nidxlist,idxlist(32),datyp,nbits,filt_pass,ip3,nlevels,nlinbot,dst_gi0,dst_gj0,dst_gin,dst_gjn,dst_dij,reduc_core(4)
   real :: levels(1024),filt_coef
   character(len=99) :: fname_S,etk_S,lvltype_S,v_int_S,grid_etk_S
   character(len=4)  :: mylist_S(MAX_ITEM),resultd_S(MAX_ITEM,0:N_STEP),resultp_S(MAX_ITEM,0:N_STEP)
   logical :: ok_L,rewrite_L,flip_L
   !---------------------------------------------------------------------
   call testutils_verbosity()
   call testutils_set_name('test_outcfg')

   resultd_S = ' '
   resultd_S(1:5,0) = (/'aa1 ','aa2 ','13  ','39  ','911 '/)
   resultd_S(1:5,1) = (/'ii  ','bb  ','cc  ','aa1 ','aa2 '/)
   resultd_S(1:5,2) = (/'aa1 ','aa2 ','13  ','39  ','911 '/)
   resultd_S(1:5,3) = (/'ii  ','bb  ','cc  ','aa1 ','aa2 '/)

   resultp_S = ' '
   resultp_S(1:2,0) = (/'uu  ','u2  '/)
   resultp_S(1:2,1) = (/'uu  ','u2  '/)
   resultp_S(1:2,2) = (/'uu  ','u2  '/)
   resultp_S(1:2,3) = (/'uu  ','u2  '/)

   dateo = 0
   dt = 1800

   fname_S = '_not-such-file.out'
   my_id(1) = outcfg_new(fname_S, 'd', dateo, dt)
   call testutils_assert_eq(RMN_IS_OK(my_id(1)),.false.,'new, File not found')

   reduc_core(1:4) = (/3,-3,3,-3/)
   fname_S = 'test_outcfg.out'
   my_id(:) = 0
   my_id(1) = outcfg_new(fname_S, 'd', dateo, dt,reduc_core,.true.)
   my_id(2) = outcfg_new(fname_S, 'p', dateo, dt,reduc_core,.true.)
   my_id(3) = outcfg_new(fname_S, 'c', dateo, dt,reduc_core,.true.)
   call testutils_assert_eq(RMN_IS_OK(minval(my_id(:))),.true.,'new, parse file')

   istat = 1
   n_items = outcfg_getlist(minval(my_id(:))-1,istep,mylist_S)
   call testutils_assert_eq(RMN_IS_OK(n_items),.false.,'getlist id out of range -1')

   n_items = outcfg_getlist(maxval(my_id(:))+1,istep,mylist_S)
   call testutils_assert_eq(RMN_IS_OK(n_items),.false.,'getlist id out of range +1')

   sum_items = 0
   do istep = 0, N_STEP
      n_items = outcfg_getlist(my_id(3),istep,mylist_S)
      sum_items = sum_items + n_items
   enddo
   call testutils_assert_eq(sum_items,0,'getlist empty list')

   istat = outcfg_time(my_id(1),0,datev,dateo,dt)
   call testutils_assert_eq((/dateo,dt/),(/0,1800/),'outcfg_time')

   ok_L = .true.
   do istep = 0, N_STEP
      n_items = outcfg_getlist(my_id(1),istep,mylist_S)
      if (any(mylist_S(:) /= resultd_S(:,istep))) then
         print *,istep,'d',n_items,mylist_S(1:n_items)
         ok_L = .false.
      endif
      if (istep == 0) then
         do ivar=1,n_items
            if (any(mylist_S(ivar) == (/'aa1','aa2'/))) then
               nidxlist = outcfg_get_idxlist(my_id(1),istep,mylist_S(ivar),idxlist)
               if (mylist_S(ivar) == 'aa1') then
                  call testutils_assert_eq(nidxlist,2,'outcfg_get_idxlist aa1')
                  call testutils_assert_eq(idxlist(1:nidxlist),(/4,5/),'outcfg_get_idxlist aa1 list')
               else
                  call testutils_assert_eq(nidxlist,1,'outcfg_get_idxlist aa2')
                  call testutils_assert_eq(idxlist(1),4,'outcfg_get_idxlist aa2 list')
                  istat = outcfg_var_meta(my_id(1),idxlist(1),mylist_S(ivar),&
                       datyp,nbits,filt_pass,filt_coef,rewrite_L,flip_L,etk_S,ip3)
                  call testutils_assert_eq(datyp,134,'var_meta datyp')
                  call testutils_assert_eq(nbits,16,'var_meta nbits')
                  call testutils_assert_eq(filt_pass,2,'var_meta filt_pass')
                  call testutils_assert_eq(filt_coef,3.,'var_meta filt_coef')
                  call testutils_assert_eq(rewrite_L,.true.,'var_meta rewrite_L')
                  call testutils_assert_eq(flip_L,.false.,'var_meta flip_L')
                  call testutils_assert_eq(etk_S,'toto','var_meta etk_S')
                  call testutils_assert_eq(ip3,123,'var_meta ip3')
                  istat = outcfg_level_meta(my_id(1),idxlist(1),lvltype_S,v_int_S,levels,nlevels,nlinbot)
                  call testutils_assert_eq(lvltype_S,'p','level_meta lvltype_S')
                  call testutils_assert_eq(v_int_S,'cubic','level_meta v_int_S')
                  call testutils_assert_eq(nlevels,2,'level_meta nlevels')
                  call testutils_assert_eq(levels(1:2),(/850.,500./),'level_meta levels')
                  call testutils_assert_eq(nlinbot,2,'level_meta nlinbot')
                  istat = outcfg_grid_meta(my_id(1),idxlist(1),dst_gi0,dst_gj0,dst_gin,dst_gjn,dst_dij,grid_etk_S)
                  call testutils_assert_eq((/dst_gi0,dst_gin,dst_gj0,dst_gjn,dst_dij/),(/3,-3,3,-3,1/),'grid_meta ij0n dij')
                  call testutils_assert_eq(grid_etk_S,' ','grid_meta etk_S')
               endif
            endif
            if (mylist_S(ivar) == '13') then
               nidxlist = outcfg_get_idxlist(my_id(1),istep,mylist_S(ivar),idxlist)
               call testutils_assert_eq(nidxlist,1,'outcfg_get_idxlist 13')
               istat = outcfg_var_meta(my_id(1),idxlist(1),mylist_S(ivar),&
                    datyp,nbits,filt_pass,filt_coef,rewrite_L,flip_L,etk_S,ip3)
               call testutils_assert_eq(nbits,32,'var_meta nbits 13')
               call testutils_assert_eq(filt_pass,3,'var_meta filt_pass 13')
               call testutils_assert_eq(filt_coef,0.1,'var_meta filt_coef 13')


               istat = outcfg_level_meta(my_id(1),idxlist(1),lvltype_S,v_int_S,levels,nlevels,nlinbot)
               call testutils_assert_eq(lvltype_S,'m','level_meta lvltype_S 13')
               call testutils_assert_eq(nlevels,5,'level_meta nlevels 13')
               call testutils_assert_eq(levels(1:5),(/1.,5.,9.,13.,17./),'level_meta levels 13')
            endif
         enddo
      endif
      n_items = outcfg_getlist(my_id(2),istep,mylist_S)
      if (any(mylist_S(:) /= resultp_S(:,istep))) then
         print *,istep,'p',n_items,mylist_S(1:n_items)
         ok_L = .false.
      endif
   enddo
   call testutils_assert_eq(ok_L,.true.,'getlist')

   !---------------------------------------------------------------------
   return
end subroutine test_outcfg
