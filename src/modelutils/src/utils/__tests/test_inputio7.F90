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

#include <rmn/msg.h>

!/@
subroutine test_inputio7()
   use, intrinsic :: iso_fortran_env, only: REAL64
   use iso_c_binding
   use vGrid_Descriptors
   use clib_itf_mod
   use env_utils, only: env_get
   use ezgrid_mod
   use inputio_mod
   use mu_jdate_mod
   use ptopo_utils
   use statfld_dm_mod
   use testutils
   use vgrid_wb
   use rmn_gmm
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04, 2017-09
!@/
#include <rmnlib_basics.hf>
   include "rpn_comm.inc"
   integer,parameter :: NDIGITS = 4
   integer,parameter :: NSTEP0 = 1
   integer,parameter :: NTEST = 1
   integer,parameter :: NVAR = 34
   integer,parameter :: NSKIP = 2
   integer,parameter :: GNI = 400 !# 80
   integer,parameter :: GNJ = 200 !# 40
   integer,parameter :: NK = 2
   integer,parameter :: MAXNK = 26
   integer,parameter :: NVALS = 10

   integer :: nstep
   character(len=512) :: dir_S,filename_S,dfiles_S, name_S, &
        anal_S,anal2_S,clim_S,clim2_S,geop_S,geop2_S,inrep_S,inrep2_S
   integer :: istat, myproc, ndomains, idomain, ngrids, igrid
   ! ---------------------------------------------------------------------
   myproc = testutils_initmpi()
   ndomains = 1
   call ptopo_init_var(ndomains, idomain, ngrids, igrid)
   istat = ptopo_io_set(testutils_npeio)

   nstep = NSTEP0
   istat = env_get('TEST_NSTEP', nstep, NSTEP0, 1, 21600)

   call testutils_set_name('test_inputio7')
   istat = fstopc('MSGLVL', 'SYSTEM', RMN_OPT_SET)
   call msg_set_p0only(0)

   istat = clib_getenv('ATM_MODEL_DFILES',dfiles_S)
   if (.not.RMN_IS_OK(istat)) then
      call msg(MSG_ERROR,'ATM_MODEL_DFILES not defined')
      return
   endif

   write(dir_S,'(a,i4.4)') '__to-delete_test_inputio7-',myproc
   istat = clib_mkdir(trim(dir_S))
   istat = clib_chdir(trim(dir_S))

!!$   anal_S = trim(dfiles_S)//'/bcmk/2009042700_000'
   anal_S = trim(dfiles_S)//'/bcmk_toctoc/2009042700_000'
!!$!   anal_S = '/users/dor/armn/sch/storage_model/zeta/ANALYSIS-mid'
!!$   anal_S = '/users/dor/armn/sch/storage_model/zeta/ANALYSIS-mid2'
   anal2_S = './ANALYSIS' !'./analysis'
   clim_S = trim(dfiles_S)//'/bcmk/climato'
   clim2_S = './CLIMATO'
   geop_S = trim(dfiles_S)//'/bcmk/geophy/Gem_geophy.fst'
   geop2_S = './GEOPHY'
!!$   inrep_S = trim(dfiles_S)//'/bcmk'
   inrep_S = trim(dfiles_S)//'/bcmk_toctoc'
   inrep2_S = './INREP'
   istat = clib_unlink(trim(anal2_S))
   istat = clib_unlink(trim(clim2_S))
   istat = clib_unlink(trim(geop2_S))
   istat = clib_unlink(trim(inrep2_S))
   istat = clib_symlink(trim(anal_S),trim(anal2_S))
   istat = min(clib_symlink(trim(clim_S),trim(clim2_S)),istat)
   istat = min(clib_symlink(trim(geop_S),trim(geop2_S)),istat)
   istat = min(clib_symlink(trim(inrep_S),trim(inrep2_S)),istat)
   if (.not.RMN_IS_OK(istat)) then
      call msg(MSG_WARNING,'problem creating symlink')
   endif

   filename_S='../test_inputio7_table'

   istat = clib_getenv('TEST_SUBNAME', name_S)
   if (name_S == 'odd') then
      call testutils_set_name('test_inputio7_odd')
      call test_input2_get_data_odd(trim(filename_S), 0, PTOPO_IO)
   else if (name_S == 'halo') then
      call testutils_set_name('test_inputio7_halo')
      call test_input2_get_data_odd(trim(filename_S), 5, PTOPO_IO)
   else if (name_S == 'bloc') then
      call testutils_set_name('test_inputio7_bloc')
      call test_input2_get_data_odd(trim(filename_S), 0, PTOPO_BLOC)
   else
      call testutils_set_name('test_inputio7_even')
      call test_input2_get_data_even(trim(filename_S), PTOPO_IO)
   endif


   istat = clib_unlink(trim(anal2_S))
   istat = clib_unlink(trim(clim2_S))
   istat = clib_unlink(trim(geop2_S))
   istat = clib_unlink(trim(inrep2_S))

   istat = clib_chdir('./..')
   istat = clib_unlink(trim(dir_S))

   call testutils_stats_print()
   call rpn_comm_barrier(RPN_COMM_WORLD, istat)
   call rpn_comm_finalize(istat)
   ! ---------------------------------------------------------------------
   return

contains

   !/@
   subroutine test_input2_get_data_even(F_filename_S, iotype)
      implicit none
      character(len=*),intent(in) :: F_filename_S
      integer, intent(in) :: iotype
      !@/
      integer, parameter :: HALO = 0
 
      integer :: allvals2(NVALS, NTEST+1, NVAR, MAXNK)
      ! ---------------------------------------------------------------------
      allvals2 = 0
      allvals2(:,1, 1, 1) = (/      31025,      43668,          0,     100000, 1, 1, 1, 4, 2, 1/) !mg
      allvals2(:,1, 2, 1) = (/     286740,      12103,     271372,     302295, 7, 4, 1, 3, 3, 1/) !tm
      allvals2(:,1, 3, 1) = (/       7933,      15720,          0,      70995, 1, 1, 1, 6, 4, 1/) !lh
      allvals2(:,1, 4, 1) = (/      10473,      42883,          0,     244947, 1, 1, 1, 6, 4, 1/) !y7
      allvals2(:,1, 5, 1) = (/       3935,      14896,          0,      84656, 1, 1, 1, 6, 4, 1/) !y8
      allvals2(:,1, 6, 1) = (/       8694,      51726,     -45704,     292515, 8, 4, 1, 6, 4, 1/) !y9
      allvals2(:,1, 7, 1) = (/       9455,      37483,          0,     208160, 1, 1, 1, 8, 4, 1/) !ga
      allvals2(:,1, 8, 1) = (/      10777,      36599,       -381,     197666, 5, 4, 1, 7, 4, 1/) !i8
      allvals2(:,1, 9, 1) = (/     272457,       2744,     257564,     273126, 7, 4, 1, 4, 1, 1/) !i9
      allvals2(:,1, 9, 2) = (/     259894,       2466,     253906,     268060, 8, 1, 1, 8, 4, 1/) !i9
      allvals2(:,1,10, 1) = (/     272392,       2657,     258042,     273120, 7, 4, 1, 4, 1, 1/) !i7
      allvals2(:,1,10, 2) = (/     270849,       2043,     262591,     273121, 7, 4, 1, 1, 2, 1/) !i7
      allvals2(:,1,10, 3) = (/     271509,       1260,     267143,     273122, 7, 4, 1, 1, 2, 1/) !i7
      allvals2(:,1,11, 1) = (/     -45700,      33533,     -69078,      14006, 1, 1, 1, 6, 4, 1/) !zp
      allvals2(:,1,12, 1) = (/       2034,       6689,        -21,      36229, 1, 1, 1, 6, 4, 1/) !sd
      allvals2(:,1,12, 2) = (/       2042,       6729,        -21,      36479, 1, 1, 1, 6, 4, 1/) !sd
      allvals2(:,1,12, 3) = (/          0,          0,          0,          0, 1, 1, 1, 1, 1, 1/) !sd
      allvals2(:,1,12, 4) = (/       3865,      13696,      -3096,      73904, 5, 4, 1, 6, 4, 1/) !sd
      allvals2(:,1,12, 5) = (/       2612,       7626,       -271,      34729, 5, 4, 1, 6, 4, 1/) !sd
      allvals2(:,1,13, 1) = (/     275559,      14012,     234906,     298654, 7, 4, 1, 7, 3, 1/) !i0
      allvals2(:,1,13, 2) = (/     276725,      15131,     237646,     301969, 7, 4, 1, 7, 3, 1/) !i0
      allvals2(:,1,14, 1) = (/      68498,      44376,          0,     100000, 4, 2, 1, 1, 1, 1/) !vf
      allvals2(:,1,14, 2) = (/       9455,      37483,          0,     208160, 1, 1, 1, 8, 4, 1/) !vf
      allvals2(:,1,14, 3) = (/       4771,      11954,          0,      45889, 1, 1, 1, 6, 4, 1/) !vf
      allvals2(:,1,14, 4) = (/       2352,      10102,          0,      55606, 1, 1, 1, 3, 4, 1/) !vf
      allvals2(:,1,14, 5) = (/       1200,       6541,          0,      37615, 1, 1, 1, 7, 3, 1/) !vf
      allvals2(:,1,14, 6) = (/       4440,      16567,          0,      88921, 1, 1, 1, 4, 4, 1/) !vf
      allvals2(:,1,14, 7) = (/       7948,      27252,          0,     119211, 1, 1, 1, 2, 2, 1/) !vf
      allvals2(:,1,14, 8) = (/          0,          0,          0,          0, 1, 1, 1, 1, 1, 1/) !vf
      allvals2(:,1,14, 9) = (/          0,          0,          0,          0, 1, 1, 1, 1, 1, 1/) !vf
      allvals2(:,1,14,10) = (/       5804,      31361,          0,     180340, 1, 1, 1, 4, 2, 1/) !vf
      allvals2(:,1,14,11) = (/       2974,       9275,          0,      45652, 1, 1, 1, 6, 4, 1/) !vf
      allvals2(:,1,14,12) = (/          0,          0,          0,          0, 1, 1, 1, 1, 1, 1/) !vf
      allvals2(:,1,14,13) = (/       6756,      26673,          0,     144180, 1, 1, 1, 2, 2, 1/) !vf
      allvals2(:,1,14,14) = (/       4776,      15851,          0,      67227, 1, 1, 1, 4, 2, 1/) !vf
      allvals2(:,1,14,15) = (/       1906,       6548,          0,      35547, 1, 1, 1, 8, 2, 1/) !vf
      allvals2(:,1,14,16) = (/          0,          0,          0,          0, 1, 1, 1, 1, 1, 1/) !vf
      allvals2(:,1,14,17) = (/          0,          0,          0,          0, 1, 1, 1, 1, 1, 1/) !vf
      allvals2(:,1,14,18) = (/          0,          0,          0,          0, 1, 1, 1, 1, 1, 1/) !vf
      allvals2(:,1,14,19) = (/          0,          0,          0,          0, 1, 1, 1, 1, 1, 1/) !vf
      allvals2(:,1,14,20) = (/       5372,      20841,          0,      90876, 1, 1, 1, 2, 3, 1/) !vf
      allvals2(:,1,14,21) = (/      12793,      49530,          0,     270057, 1, 1, 1, 8, 2, 1/) !vf
      allvals2(:,1,14,22) = (/       2485,      10862,          0,      60039, 1, 1, 1, 6, 4, 1/) !vf
      allvals2(:,1,14,23) = (/       2166,      10248,          0,      58942, 1, 1, 1, 7, 3, 1/) !vf
      allvals2(:,1,14,24) = (/       7001,      24033,          0,      99982, 1, 1, 1, 1, 3, 1/) !vf
      allvals2(:,1,14,25) = (/       2840,      13052,          0,      74823, 1, 1, 1, 3, 4, 1/) !vf
      allvals2(:,1,14,26) = (/       3467,      12143,          0,      61218, 1, 1, 1, 2, 4, 1/) !vf
      allvals2(:,1,15, 1) = (/      10033,      22938,          0,      95290, 1, 1, 1, 7, 4, 1/) !lg
      allvals2(:,1,16, 1) = (/      79701,      35098,        331,     100000, 1, 3, 1, 1, 1, 1/) !hs
      allvals2(:,1,17, 1) = (/      57188,      86755,      10000,     260000, 1, 1, 1, 8, 2, 1/) !vg
      allvals2(:,1,18, 1) = (/      13858,      26106,       -294,      92369, 7, 4, 1, 6, 4, 1/) !me
      allvals2(:,1,19, 1) = (/      12258,       8690,       6000,      40099, 1, 2, 1, 1, 3, 1/) !al
      allvals2(:,1,20, 1) = (/      70875,      40525,       3840,      99983, 2, 3, 1, 1, 2, 1/) !i1
      allvals2(:,1,20, 2) = (/      69457,      42359,        -49,      99951, 6, 4, 1, 1, 1, 1/) !i1
      allvals2(:,1,21, 1) = (/      12062,      39967,        -76,     181321, 1, 1, 1, 2, 4, 1/) !i2
      allvals2(:,1,22, 1) = (/       1966,       9664,       -114,      55355, 8, 2, 1, 7, 3, 1/) !i3
      allvals2(:,1,23, 1) = (/     -18397,      35319,     -28163,     159337, 1, 1, 1, 3, 4, 1/) !i4
      allvals2(:,1,24, 1) = (/      78138,       6095,      50882,      79983, 3, 4, 1, 1, 1, 1/) !i6
      allvals2(:,1,25, 1) = (/      19985,      24870,          0,      63573, 1, 1, 1, 4, 2, 1/) !j1
      allvals2(:,1,25, 2) = (/      20088,      25051,          0,      65914, 1, 1, 1, 6, 4, 1/) !j1
      allvals2(:,1,25, 3) = (/      20085,      25045,          0,      65813, 1, 1, 1, 6, 4, 1/) !j1
      allvals2(:,1,26, 1) = (/       7740,      10402,          0,      36881, 1, 1, 1, 7, 3, 1/) !j2
      allvals2(:,1,26, 2) = (/       7705,      10384,          0,      36881, 1, 1, 1, 7, 3, 1/) !j2
      allvals2(:,1,26, 3) = (/       7703,      10383,          0,      36881, 1, 1, 1, 7, 3, 1/) !j2
      allvals2(:,1,27, 1) = (/     129281,      51000,      99992,     302992, 1, 1, 1, 3, 4, 1/) !dn
      allvals2(:,1,28, 1) = (/      93750,      24206,          0,     100000, 6, 3, 1, 1, 1, 1/) !ic
      allvals2(:,1,29, 1) = (/      12258,       8690,       6000,      40099, 1, 2, 1, 1, 3, 1/) !al
      allvals2(:,1,30, 1) = (/      99247,       2805,      91663,     102535, 6, 4, 1, 7, 4, 1/) !p0
      allvals2(:,1,31, 1) = (/     -50825,       5898,     -65628,     -44011, 4, 4, 1, 5, 2, 1/) !tt
      allvals2(:,1,31, 2) = (/     -61410,     122698,    -295814,      92272, 3, 1, 1, 5, 2, 1/) !tt
      allvals2(:,1,31, 3) = (/      36842,     124183,    -158518,     236214, 3, 1, 1, 2, 3, 1/) !tt
      allvals2(:,1,31, 4) = (/       9108,      12589,     -17545,      32771, 7, 4, 1, 1, 3, 1/) !tt
      allvals2(:,1,31, 5) = (/     114016,     115340,    -149090,     315954, 7, 4, 1, 7, 3, 1/) !tt
      allvals2(:,1,31, 6) = (/      11519,      11547,     -14728,      31782, 7, 4, 1, 7, 3, 1/) !tt
      allvals2(:,1,32, 1) = (/      18270,      18179,      -5007,      59495, 2, 2, 1, 5, 1, 1/) !uu
      allvals2(:,1,32, 2) = (/      10575,      15815,     -14396,      45132, 7, 2, 1, 5, 1, 1/) !uu
      allvals2(:,1,32, 3) = (/       3735,      17560,     -27193,      43755, 8, 4, 1, 3, 1, 1/) !uu
      allvals2(:,1,32, 4) = (/       1464,      13890,     -21563,      33990, 8, 4, 1, 3, 1, 1/) !uu
      allvals2(:,1,32, 5) = (/       1518,      13665,     -19312,      33990, 8, 4, 1, 3, 1, 1/) !uu
      allvals2(:,1,33, 1) = (/       3551,      17091,     -32716,      49408, 6, 4, 1, 8, 4, 1/) !vv
      allvals2(:,1,33, 2) = (/       2578,      13479,     -21073,      42748, 8, 1, 1, 8, 4, 1/) !vv
      allvals2(:,1,33, 3) = (/      18153,     110331,    -229775,     272539, 8, 1, 1, 1, 4, 1/) !vv
      allvals2(:,1,33, 4) = (/       7855,      81518,    -187312,     193351, 5, 4, 1, 1, 4, 1/) !vv
      allvals2(:,1,33, 5) = (/       6809,      78478,    -174487,     176514, 5, 4, 1, 1, 4, 1/) !vv
      allvals2(:,1,34, 1) = (/     100000,          0,     100000,     100000, 1, 1, 1, 1, 1, 1/) !en
      allvals2(:,1,34, 2) = (/      99989,         23,      99908,      99999, 6, 4, 1, 7, 4, 1/) !en
      allvals2(:,1,34, 3) = (/      17419,      40250,       9243,     241449, 2, 3, 1, 6, 4, 1/) !en
      allvals2(:,1,34, 4) = (/       2758,      10351,         44,      59221, 2, 1, 1, 6, 4, 1/) !en
      allvals2(:,1,34, 5) = (/       5009,      13225,         42,      59221, 8, 3, 1, 6, 4, 1/) !en
      allvals2(:,1,34, 6) = (/       5012,      13224,         39,      59221, 8, 3, 1, 6, 4, 1/) !en
      allvals2(:,NTEST+1, 1, 1) = (/      99224,       2826,      91707,     102710, 6, 4, 1, 7, 4, 1/) !p0
      allvals2(:,NTEST+1, 2, 1) = (/     -50817,       5782,     -64966,     -44149, 4, 1, 1, 5, 2, 1/) !tt
      allvals2(:,NTEST+1, 2, 2) = (/     -60720,     120622,    -291595,      87328, 3, 1, 1, 5, 2, 1/) !tt
      allvals2(:,NTEST+1, 2, 3) = (/      37397,     121710,    -155781,     242078, 3, 1, 1, 2, 3, 1/) !tt
      allvals2(:,NTEST+1, 2, 4) = (/       9406,      12815,     -17603,      33464, 7, 4, 1, 1, 3, 1/) !tt
      allvals2(:,NTEST+1, 2, 5) = (/     120774,     122447,    -153555,     306834, 7, 4, 1, 1, 3, 1/) !tt
      allvals2(:,NTEST+1, 2, 6) = (/     121892,     122598,    -151969,     306834, 7, 4, 1, 1, 3, 1/) !tt
      allvals2(:,NTEST+1, 3, 1) = (/      18367,      15539,      -2519,      52925, 7, 2, 1, 5, 1, 1/) !uu
      allvals2(:,NTEST+1, 3, 2) = (/      10739,      15328,     -13496,      44210, 7, 2, 1, 7, 1, 1/) !uu
      allvals2(:,NTEST+1, 3, 3) = (/       4368,      16583,     -20466,      38528, 7, 3, 1, 7, 1, 1/) !uu
      allvals2(:,NTEST+1, 3, 4) = (/      15351,     129156,    -195148,     301235, 7, 2, 1, 3, 1, 1/) !uu
      allvals2(:,NTEST+1, 3, 5) = (/      16126,     127504,    -193826,     301235, 7, 2, 1, 3, 1, 1/) !uu
      allvals2(:,NTEST+1, 4, 1) = (/      16142,     139461,    -311726,     222863, 2, 4, 1, 1, 4, 1/) !vv
      allvals2(:,NTEST+1, 4, 2) = (/       5950,     110616,    -213336,     195059, 5, 1, 1, 7, 1, 1/) !vv
      allvals2(:,NTEST+1, 4, 3) = (/       4472,      93830,    -181676,     213857, 8, 1, 1, 1, 4, 1/) !vv
      allvals2(:,NTEST+1, 4, 4) = (/      -1130,      73921,    -146847,     158265, 6, 1, 1, 1, 4, 1/) !vv
      allvals2(:,NTEST+1, 4, 5) = (/      -1812,      71917,    -145554,     144288, 6, 1, 1, 1, 4, 1/) !vv

      call test_input7io_run(F_filename_S, GNI, GNJ, HALO, allvals2, iotype)
      ! ---------------------------------------------------------------------
      return
   end subroutine test_input2_get_data_even


   !/@
   subroutine test_input2_get_data_odd(F_filename_S, halo, iotype)
      implicit none
      character(len=*),intent(in) :: F_filename_S
      integer, intent(in) :: halo, iotype
      !@/

      integer :: allvals2(NVALS, NTEST+1, NVAR, MAXNK)
       ! ---------------------------------------------------------------------
      allvals2 = 0
      allvals2(:,1, 1, 1) = (/      34196,      44721,          0,     100000,40, 4, 1, 1, 1, 1/) !mg
      allvals2(:,1, 2, 1) = (/     285142,      11912,     268955,     303971,63, 3, 1,14,21, 1/) !tm
      allvals2(:,1, 3, 1) = (/      12376,      25827,          0,     288624,44, 2, 1,65,18, 1/) !lh
      allvals2(:,1, 4, 1) = (/        989,       3774,          0,      63613,44, 2, 1,18,29, 1/) !y7
      allvals2(:,1, 5, 1) = (/        681,       2661,          0,      54390,44, 2, 1,18,29, 1/) !y8
      allvals2(:,1, 6, 1) = (/        318,       7124,    -127178,     179936,17,30, 1,18,29, 1/) !y9
      allvals2(:,1, 7, 1) = (/      11310,      30656,          0,     100000,40, 4, 1, 1, 1, 1/) !ga
      allvals2(:,1, 8, 1) = (/       3282,       8683,      -1425,      39864,41,37, 1,67,41, 1/) !i8
      allvals2(:,1, 9, 1) = (/     266128,      14886,     199167,     273398,22, 3, 1,67, 6, 1/) !i9
      allvals2(:,1, 9, 2) = (/     255084,      11582,     195812,     272158,23, 3, 1,78,36, 1/) !i9
      allvals2(:,1,10, 1) = (/     266062,      14751,     199113,     273565,22, 3, 1,67, 6, 1/) !i7
      allvals2(:,1,10, 2) = (/     267607,       8339,     222705,     273130,23, 1, 1, 1,11, 1/) !i7
      allvals2(:,1,10, 3) = (/     269684,       4666,     247013,     273129,23, 1, 1, 1,11, 1/) !i7
      allvals2(:,1,11, 1) = (/     -47627,      31346,     -82504,      35208,43, 3, 1,18,29, 1/) !zp
      allvals2(:,1,12, 1) = (/       3184,      10891,      -3727,      62848,50,35, 1,57, 4, 1/) !sd
      allvals2(:,1,12, 2) = (/       3184,      10891,      -3727,      62848,50,35, 1,57, 4, 1/) !sd
      allvals2(:,1,12, 3) = (/          0,          0,          0,          0, 1, 1, 1, 1, 1, 1/) !sd
      allvals2(:,1,12, 4) = (/       2279,       8397,      -3560,      61765,35, 5, 1,52, 4, 1/) !sd
      allvals2(:,1,12, 5) = (/       4900,      12631,      -3727,      61973,50,35, 1,51,35, 1/) !sd
      allvals2(:,1,13, 1) = (/     269393,      18644,     219967,     308138,30, 2, 1,33,18, 1/) !i0
      allvals2(:,1,13, 2) = (/     270617,      18608,     223251,     309193,71,39, 1, 9,25, 1/) !i0
      allvals2(:,1,14, 1) = (/      65129,      45396,          0,     100000, 1, 1, 1,40, 4, 1/) !vf
      allvals2(:,1,14, 2) = (/      11310,      30656,          0,     100000,40, 4, 1, 1, 1, 1/) !vf
      allvals2(:,1,14, 3) = (/        675,       4191,          0,     100000, 1, 1, 1, 8,31, 1/) !vf
      allvals2(:,1,14, 4) = (/       1041,       6299,          0,      84648, 1, 1, 1, 4,35, 1/) !vf
      allvals2(:,1,14, 5) = (/       1627,       9908,          0,      99234, 1, 1, 1,66,19, 1/) !vf
      allvals2(:,1,14, 6) = (/        781,       6366,          0,      96750, 1, 1, 1,29,35, 1/) !vf
      allvals2(:,1,14, 7) = (/       1079,       6025,          0,      94780, 1, 1, 1,28,37, 1/) !vf
      allvals2(:,1,14, 8) = (/          0,          0,          0,          0, 1, 1, 1, 1, 1, 1/) !vf
      allvals2(:,1,14, 9) = (/          0,          0,          0,          0, 1, 1, 1, 1, 1, 1/) !vf
      allvals2(:,1,14,10) = (/        554,       5037,          0,      95938, 1, 1, 1,29,16, 1/) !vf
      allvals2(:,1,14,11) = (/        422,       3402,          0,      70494, 1, 1, 1,67,16, 1/) !vf
      allvals2(:,1,14,12) = (/          0,          0,          0,          0, 1, 1, 1, 1, 1, 1/) !vf
      allvals2(:,1,14,13) = (/       1860,       9157,          0,      99641, 1, 1, 1,16,32, 1/) !vf
      allvals2(:,1,14,14) = (/        998,       6342,          0,      92513, 1, 1, 1,33,17, 1/) !vf
      allvals2(:,1,14,15) = (/       2668,      10951,          0,      94113, 1, 1, 1,19,27, 1/) !vf
      allvals2(:,1,14,16) = (/          0,          0,          0,          0, 1, 1, 1, 1, 1, 1/) !vf
      allvals2(:,1,14,17) = (/          0,          0,          0,          0, 1, 1, 1, 1, 1, 1/) !vf
      allvals2(:,1,14,18) = (/          0,          0,          0,          0, 1, 1, 1, 1, 1, 1/) !vf
      allvals2(:,1,14,19) = (/          0,          0,          0,          0, 1, 1, 1, 1, 1, 1/) !vf
      allvals2(:,1,14,20) = (/        480,       4211,          0,      95114, 1, 1, 1,27,29, 1/) !vf
      allvals2(:,1,14,21) = (/        475,       3552,          0,     142872, 1, 1, 1, 1,33, 1/) !vf
      allvals2(:,1,14,22) = (/       2992,      13402,          0,      99986, 1, 1, 1,23,36, 1/) !vf
      allvals2(:,1,14,23) = (/        238,       2466,          0,      67607, 1, 1, 1,20,36, 1/) !vf
      allvals2(:,1,14,24) = (/       4394,      17718,          0,     100000, 1, 1, 1,12,26, 1/) !vf
      allvals2(:,1,14,25) = (/        874,       5318,          0,      79128, 1, 1, 1,66,12, 1/) !vf
      allvals2(:,1,14,26) = (/       2830,      10014,          0,      96025, 1, 1, 1, 7,22, 1/) !vf
      allvals2(:,1,15, 1) = (/      27075,      40462,          0,     100000, 1, 7, 1,44, 1, 1/) !lg
      allvals2(:,1,16, 1) = (/      88576,      26950,        110,     100000, 8,26, 1, 1, 1, 1/) !hs
      allvals2(:,1,17, 1) = (/      47507,      76063,      10000,     260000,40, 4, 1,35,14, 1/) !vg
      allvals2(:,1,18, 1) = (/       3671,       8199,       -890,      54581,58,26, 1,19,29, 1/) !me
      allvals2(:,1,19, 1) = (/      18335,      21115,       6000,      79999, 1,15, 1, 1, 1, 1/) !al
      allvals2(:,1,20, 1) = (/      78435,      36260,       -506,     107844,17,30, 1,33,38, 1/) !i1
      allvals2(:,1,20, 2) = (/      77354,      38047,     -11963,     109375,23, 6, 1,33,38, 1/) !i1
      allvals2(:,1,21, 1) = (/       1309,       5295,      -2498,      47478,30, 6, 1,23, 6, 1/) !i2
      allvals2(:,1,22, 1) = (/       1883,      10353,      -6559,     116488,68,22, 1, 4,21, 1/) !i3
      allvals2(:,1,23, 1) = (/         23,       1805,       -247,     101159, 3,32, 1,21,31, 1/) !i4
      allvals2(:,1,24, 1) = (/      76789,       9117,      42007,      81790,21,31, 1,69,32, 1/) !i6
      allvals2(:,1,25, 1) = (/      18874,      22932,          0,      86389,40, 4, 1, 7,24, 1/) !j1
      allvals2(:,1,25, 2) = (/      18846,      22879,          0,      86389,40, 4, 1, 7,24, 1/) !j1
      allvals2(:,1,25, 3) = (/      18752,      22767,          0,      86389,40, 4, 1, 7,24, 1/) !j1
      allvals2(:,1,26, 1) = (/       8665,      11008,          0,      52165,40, 4, 1,70,15, 1/) !j2
      allvals2(:,1,26, 2) = (/       8698,      11042,          0,      52165,40, 4, 1,70,15, 1/) !j2
      allvals2(:,1,26, 3) = (/       8789,      11170,          0,      52165,40, 4, 1,70,15, 1/) !j2
      allvals2(:,1,27, 1) = (/      11813,       4686,       7699,      74674,40, 3, 1,64,40, 1/) !dn
      allvals2(:,1,28, 1) = (/      90424,      29130,          0,     100000,44,21, 1, 1, 1, 1/) !ic
      allvals2(:,1,29, 1) = (/      18335,      21115,       6000,      79999, 1,15, 1, 1, 1, 1/) !al
      allvals2(:,1,30, 1) = (/      96708,       9311,      52368,     103932,19,29, 1,50,33, 1/) !p0
      allvals2(:,1,31, 1) = (/     -51880,       5777,     -68870,     -39429,79, 5, 1,35,30, 1/) !tt
      allvals2(:,1,31, 2) = (/     -10336,      14664,     -71637,       8778,10, 3, 1,52,18, 1/) !tt
      allvals2(:,1,31, 3) = (/      -2188,      18773,     -71637,      24859,10, 3, 1, 8,24, 1/) !tt
      allvals2(:,1,31, 4) = (/       2162,      20968,     -71637,      36626,10, 3, 1, 1,25, 1/) !tt
      allvals2(:,1,31, 5) = (/       4132,      21779,     -71637,      37261,10, 3, 1,13,27, 1/) !tt
      allvals2(:,1,31, 6) = (/       4211,      21818,     -71637,      36706,10, 3, 1,13,27, 1/) !tt
      allvals2(:,1,32, 1) = (/      12857,      20087,     -30183,     101860,69, 3, 1,63, 9, 1/) !uu
      allvals2(:,1,32, 2) = (/       6463,      16419,     -45177,      78180,20, 5, 1,66,10, 1/) !uu
      allvals2(:,1,32, 3) = (/       2124,      15317,     -63268,      59944,24, 6, 1,17,12, 1/) !uu
      allvals2(:,1,32, 4) = (/       -131,      12233,     -48312,      46154,24, 6, 1,17,12, 1/) !uu
      allvals2(:,1,32, 5) = (/       -194,      11953,     -48312,      44298,24, 6, 1,17,12, 1/) !uu
      allvals2(:,1,33, 1) = (/         95,      17365,     -71452,      84990,35,36, 1,45,34, 1/) !vv
      allvals2(:,1,33, 2) = (/         61,      13810,     -55760,      59737,77, 5, 1,46,34, 1/) !vv
      allvals2(:,1,33, 3) = (/        350,      13410,     -60548,      53235,78, 6, 1,45,37, 1/) !vv
      allvals2(:,1,33, 4) = (/        152,      11664,     -43285,      41192,76,10, 1,47,33, 1/) !vv
      allvals2(:,1,33, 5) = (/        157,      11357,     -42667,      39788,78, 6, 1,47,33, 1/) !vv
      allvals2(:,1,34, 1) = (/     100000,          0,     100000,     100000, 1, 1, 1, 1, 1, 1/) !en
      allvals2(:,1,34, 2) = (/        366,       4441,         35,     138609,12, 4, 1,23,28, 1/) !en
      allvals2(:,1,34, 3) = (/       1285,      11566,         10,     240668,50,35, 1,65,18, 1/) !en
      allvals2(:,1,34, 4) = (/       3470,      16008,       -962,     240668,50,35, 1,65,18, 1/) !en
      allvals2(:,1,34, 5) = (/       4383,      16759,      -1774,     240668,64,23, 1,65,18, 1/) !en
      allvals2(:,1,34, 6) = (/       4390,      16760,      -1774,     240668,64,23, 1,65,18, 1/) !en
      allvals2(:,NTEST+1, 1, 1) = (/      96705,       9292,      52388,     103856,19,29, 1,50,33, 1/) !p0
      allvals2(:,NTEST+1, 2, 1) = (/     -51969,       5671,     -68547,     -41848,26, 5, 1,70,34, 1/) !tt
      allvals2(:,NTEST+1, 2, 2) = (/     -10253,      14466,     -70985,       9067,10, 3, 1,61,17, 1/) !tt
      allvals2(:,NTEST+1, 2, 3) = (/      -2023,      18446,     -70985,      24512,10, 3, 1,11,26, 1/) !tt
      allvals2(:,NTEST+1, 2, 4) = (/       2401,      20711,     -70985,      36818,10, 3, 1, 1,25, 1/) !tt
      allvals2(:,NTEST+1, 2, 5) = (/       4527,      21737,     -70985,      35911,10, 3, 1,17,26, 1/) !tt
      allvals2(:,NTEST+1, 2, 6) = (/       4611,      21779,     -70985,      35288,10, 3, 1, 8,25, 1/) !tt
      allvals2(:,NTEST+1, 3, 1) = (/      12574,      19609,     -32330,      98854,22, 5, 1,64, 9, 1/) !uu
      allvals2(:,NTEST+1, 3, 2) = (/       6443,      15951,     -47924,      71374,20, 5, 1,66,10, 1/) !uu
      allvals2(:,NTEST+1, 3, 3) = (/       2328,      14597,     -48934,      51626,24, 6, 1,67, 9, 1/) !uu
      allvals2(:,NTEST+1, 3, 4) = (/        169,      11652,     -40292,      34540,24, 6, 1,45, 9, 1/) !uu
      allvals2(:,NTEST+1, 3, 5) = (/         98,      11384,     -40292,      33388,24, 6, 1,70,34, 1/) !uu
      allvals2(:,NTEST+1, 4, 1) = (/         24,      16242,     -63014,      70746,35,36, 1,46,36, 1/) !vv
      allvals2(:,NTEST+1, 4, 2) = (/         66,      13120,     -51908,      50466,77, 3, 1,46,34, 1/) !vv
      allvals2(:,NTEST+1, 4, 3) = (/        310,      12442,     -55015,      44291,77, 8, 1,47,35, 1/) !vv
      allvals2(:,NTEST+1, 4, 4) = (/        187,      10831,     -40467,      38876,76, 3, 1,33, 5, 1/) !vv
      allvals2(:,NTEST+1, 4, 5) = (/        192,      10545,     -40467,      38876,76, 3, 1,33, 5, 1/) !vv

      call test_input7io_run(F_filename_S, GNI, GNJ, HALO, allvals2, iotype)
      ! ---------------------------------------------------------------------
      return
   end subroutine test_input2_get_data_odd


   !/@
   subroutine test_input7io_run(F_filename_S, GNI, GNJ, HALO, allvals2, iotype)
      implicit none
      character(len=*),intent(in) :: F_filename_S
      integer,intent(in) :: GNI, GNJ, HALO
      integer :: allvals2(:,:,:,:), iotype
      !@/
      logical, parameter :: ALONGX = .true.
      logical, parameter :: FILL = .false. !.true.
      real,parameter :: hyb(4) = (/0.5492443,0.7299818,0.8791828,0.9950425/)
      character(len=512) :: dateo_S,varname_S,varname2_S,skip_list_S(4),dummy_S,step_S
      character(len=8) :: result_S(0:NTEST,NVAR)
      integer :: istat,dateo,istep,itest,dt,ii,gridid,nbvar,ivar,nhyb,k,ikind,ip1list(NK)
      integer :: mini, maxi, lni, lnimax, li0, minj, maxj, lnj, lnjmax, lj0, rpncomm_gridid,lijk0(3),lijk1(3),uijk0(3),uijk1(3)
      logical :: ok_L(0:NTEST,NVAR), ok1_L
      real :: hyblist(NK)
      real(8) :: ptop_8
      real, pointer :: data(:,:,:),data2(:,:,:),p0data(:,:)
      integer,pointer :: ip1listt(:), ip1listm(:)
      type(vgrid_descriptor) :: vgrid1
      type(gmm_metadata) :: p0meta

      type(INPUTIO_T) :: inputobj
      integer(INT64) :: jdateo

      integer :: ijkmin(3),ijkmax(3),ilvl,vals2(10), lclgridid
      real(REAL64) :: mean,var,rmin,rmax, fact
      ! ---------------------------------------------------------------------
      result_S = ' '

      result_S(0,1:NVAR) = (/ &
           'mg  ','tm  ','lh  ','y7  ','y8  ', &
           'y9  ','ga  ','i8  ','i9  ','i7  ', &
           'zp  ','sd  ','i0  ','vf  ','lg  ', &
           'hs  ','vg  ','me  ','al  ','i1  ', &
           'i2  ','i3  ','i4  ','i6  ','j1  ', &
           'j2  ','dn  ','icel','al  ','p0  ', &
           'tt  ','uu  ','vv  ','en  '/)
!!$      result_S(1,1:5) = (/ 'p0  ','tt  ','uu  ','vv  ','en  '/)
      result_S(NTEST,1:4) = (/ 'p0  ','tt  ','uu  ','vv  '/)

!!$      skip_list_S(1:3) = (/'xa  ','icel','fsa '/)

      ok_L = .false.
      gridid = ezqkdef(GNI,GNJ, 'G', 0,0,0,0,0)
      istat = rpn_comm_topo(GNI, mini, maxi, lni, lnimax, HALO, li0, &
           ALONGX, FILL)
      istat = rpn_comm_topo(GNJ, minj, maxj, lnj, lnjmax, HALO, lj0, &
           .not.ALONGX, FILL)
      if (iotype == PTOPO_IO) then
         rpncomm_gridid = rpn_comm_create_2dgrid(GNI, GNJ, mini, maxi, minj, maxj)
         call testutils_assert_ok(RMN_IS_OK(rpncomm_gridid), 'rpncomm_gridid')
      else

         rpncomm_gridid = -1
         lclgridid = ezgrid_sub(gridid, li0, lj0, li0+lni-1, lj0+lnj-1)
         if (.not.RMN_IS_OK(lclgridid)) then
            write(dummy_S, '(i7,a,i7,4i6,a)') lclgridid,"= ezgrid_sub(",gridid,li0, lj0, li0+lni-1, lj0+lnj-1,")"
            call msg(MSG_ERROR,'test_input7 - cannot creat lclgridid: '//trim(dummy_S))
            print *,"ERROR: ",trim(dummy_S)
            return
            lclgridid = gridid
         endif
         gridid = lclgridid
      endif

      lijk0(1:3) = (/1, 1, 1/)
      uijk0(1:3) = (/lni, lnj, NK/)

      call gmm_build_meta2D(p0meta, &
           1,lni,HALO,HALO,lni, &
           1,lnj,HALO,HALO,lnj, &
           0,GMM_NULL_FLAGS)
      nullify(p0data)
      istat = gmm_create('P0',p0data,p0meta)
      p0data = 99900.

      ptop_8 = 9575.0d0
      istat = vgd_new(vgrid1, &
           kind     = VGRID_HYBS_KIND, &
           version  = VGRID_HYBS_VER, &
           hyb      = hyb, &
           ptop_8   = ptop_8, &
           pref_8   = 100000.d0, &
           rcoef1   = 1., &
           rcoef2   = 1.)
      nullify(ip1listt, ip1listm)
      istat = vgd_get(vgrid1,key='VIPT',value=ip1listt)
      istat = vgrid_wb_put('gemthlvl',vgrid1,ip1listt,'P0')
      istat = vgd_get(vgrid1,key='VIPM',value=ip1listm)
      istat = vgrid_wb_put('gemmolvl',vgrid1,ip1listm,'P0')

      dateo_S = '20090427.000000'
      jdateo = jdate_from_print(dateo_S)
      dateo = jdate_to_cmc(jdateo)
      nhyb = 2
      hyblist(1:2) = (/0.985,0.995/)
      do k=1,nhyb
         ikind = RMN_CONV_HY
         call convip_plus(ip1list(k),hyblist(k),ikind,RMN_CONV_P2IPNEW,'',.false.)
      enddo
      dt = 21600 / NSTEP
      istat = inputio_new(inputobj, jdateo, dt, F_filename_S, '.', &
           gridid, gridid, rpncomm_gridid, &
           F_li0=1, F_lin=lni, F_lj0=1, F_ljn=lnj, F_iotype=iotype)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_ERROR,'test_input3 - cannot init inputio_mod for: '//trim(F_filename_S))
         return
      endif
      nbvar = inputio_nbvar(inputobj)
      call testutils_assert_eq(nbvar, NVAR+NSKIP, 'inputio_nbvar')

      itest = 0 
      STEPLOOP: do istep = 0,NSTEP
         write(step_S, '(i5.5)') istep
         call msg(MSG_INFO, '============================================ '//trim(step_S))
         if (istep == NSTEP) itest = NTEST
         ii = 0
         VARLOOP: do ivar=1,nbvar
            istat = inputio_isvarstep(inputobj, ivar, istep)
            if (.not.RMN_IS_OK(istat)) cycle VARLOOP
            istat = inputio_meta(inputobj%cfg, ivar, varname_S, varname2_S)
            if (any(varname_S == skip_list_S)) cycle VARLOOP
            call msg(MSG_INFO, '-------------------------------------------- '//trim(varname_S))

            nullify(data,data2)
            if (varname_S == 'uu') then
               istat = inputio_get(inputobj, ivar, istep, data, data2, &
                    F_vgrid_S='gemmolvl')
            else
               istat = inputio_get(inputobj, ivar, istep, data, data2, &
                    F_vgrid_S='gemthlvl')
            endif
            if (.not.(RMN_IS_OK(istat) .and. associated(data))) then
               call msg(MSG_WARNING,'test_input3 - var not found: '//trim(varname_S))
               cycle VARLOOP !exit STEPLOOP
            endif
!!$            print *,'====',inputobj%fid%files(1:inputobj%fid%nfiles)%unit

            IF_ASSERT: if (istep == 0 .or. istep == NSTEP) then
               ii = ii + 1
               call testutils_assert_eq(varname_S, result_S(itest, ii), 'var '//trim(varname_S)//':'//trim(step_S))

               lijk1 = lbound(data)
               uijk1 = ubound(data)
               call testutils_assert_ok( &
                    all(lijk1(1:2) == lijk0(1:2)) .and. &
                    all(uijk1(1:2) == uijk0(1:2)), 'bounds '//trim(varname_S)//':'//trim(step_S))
               if (.not.(all(lijk1(1:2) == lijk0(1:2)) .and. &
                    all(uijk1(1:2) == uijk0(1:2)))) then
                  print *,'bounds  '//trim(varname_S)//':'//trim(step_S)//':',maxi,maxj
                  print *,'bounds0 '//trim(varname_S)//':'//trim(step_S)//':',uijk0(1:2)
                  print *,'bounds1 '//trim(varname_S)//':'//trim(step_S)//':',uijk1(1:2)
                  call flush(6)
               endif
               if (varname_S == 'tt') then
                  call testutils_assert_ok( &
                       lijk1(3) == 1 .and. uijk1(3)==size(ip1listt), 'bounds k '//trim(varname_S)//':'//trim(step_S))
               elseif (varname_S == 'uu') then
                  call testutils_assert_ok( &
                       lijk1(3) == 1 .and. uijk1(3)==size(ip1listm), 'bounds k '//trim(varname_S)//':'//trim(step_S))
               endif

               call testutils_assert_not_naninf(data, 'validity '//trim(varname_S)//':'//trim(step_S))
               if (varname2_S/=' ') &
                    call testutils_assert_not_naninf(data2, 'validity '//trim(varname2_S)//':'//trim(step_S))

               do ilvl=1,uijk1(3)
                  call statfld_dm(data(:,:,ilvl), mean, var, rmin, rmax, ijkmin, ijkmax)
                  ok1_L = .true.
                  if (ptopo_grid_ipe == RPN_COMM_MASTER) then
                     fact = 10.**nint(5.-max(-20.D0,log10(maxval((/abs(mean), abs(var), abs(rmin), abs(rmax)/)) + tiny(rmax))))
                     vals2 = (/nint(fact*mean), nint(fact*var), nint(fact*rmin), nint(fact*rmax), ijkmin(1), ijkmin(2), ijkmin(3), ijkmax(1), ijkmax(2), ijkmax(3)/)
                     if (all(vals2(:) == allvals2(:,itest+1,ii,ilvl))) then
                        ok1_L = .true.
                     else
                        ok1_L = .false.
                        print &
                             '(a,i1,a,i2,a,i2,a,i11,",",i11,",",i11,",",i11,",",i2,",",i2,",",i2,",",i2,",",i2,",",i2,a)', &
                             "allvals2(:,", itest+1, ",", ii, ",", ilvl, ") = (/", vals2, "/) !"//varname_S(1:2)
                        print &
                             '(a,i1,a,i2,a,i2,a,i11,",",i11,",",i11,",",i11,",",i2,",",i2,",",i2,",",i2,",",i2,",",i2,a)', &
                             "allvals0(:,", itest+1, ",", ii, ",", ilvl, ") = (/", allvals2(:,itest+1,ii,ilvl), "/) !"//varname_S(1:2)
                     endif
                  endif
                  write(dummy_S,*) ilvl
                  call testutils_assert_ok(ok1_L, 'stats '//varname_S(1:2)//':'//trim(step_S)//' '//trim(dummy_S))
               enddo
               if (varname2_S/=' ') then
                  ii = ii + 1
                  do ilvl=1,uijk1(3)
                     call statfld_dm(data2(:,:,ilvl), mean, var, rmin, rmax, ijkmin, ijkmax)
                     ok1_L = .true.
                     if (ptopo_grid_ipe == RPN_COMM_MASTER) then
                        fact = 10.**nint(5.-max(-20.D0,log10(maxval((/abs(mean), abs(var), abs(rmin), abs(rmax)/)) + tiny(rmax))))
                        vals2 = (/nint(fact*mean), nint(fact*var), nint(fact*rmin), nint(fact*rmax), ijkmin(1), ijkmin(2), ijkmin(3), ijkmax(1), ijkmax(2), ijkmax(3)/)
                        if (all(vals2(:) == allvals2(:,itest+1,ii,ilvl))) then
                           ok1_L = .true.
                        else
                           ok1_L = .false.
                           print &
                                '(a,i1,a,i2,a,i2,a,i11,",",i11,",",i11,",",i11,",",i2,",",i2,",",i2,",",i2,",",i2,",",i2,a)', &
                                "allvals2(:,", itest+1, ",", ii, ",", ilvl, ") = (/", vals2, "/) !"//varname2_S(1:2)
                        endif
                     endif
                     write(dummy_S,*) ilvl
                     call testutils_assert_ok(ok1_L, 'stats '//varname2_S(1:2)//':'//trim(step_S)//' '//trim(dummy_S))
                  enddo
               endif
            endif IF_ASSERT
            if (associated(data)) deallocate(data,stat=istat)
            if (associated(data2)) deallocate(data2,stat=istat)
         enddo VARLOOP
         IF_ASSERT2: if (istep == 0 .or. istep == NSTEP) then
            if (ii < NVAR) then
!!$               print *,'ii,nvar,itest,result_S(itest, ii):',ii,nvar,itest,trim(result_S(itest, ii))
               ii = ii + 1
               if (result_S(itest, ii) == '') then
                  call testutils_assert_ok(.true., 'nvar'//trim(step_S))
               else
                  call testutils_assert_eq(' ', result_S(itest, ii), 'missing "'//trim(result_S(itest, ii))//'":'//trim(step_S))
               endif
            else
               call testutils_assert_ok(.true., 'nvar'//trim(step_S))
            endif
         endif IF_ASSERT2
      enddo STEPLOOP
      istat = inputio_close_files(inputobj)
      ! ---------------------------------------------------------------------
      return
   end subroutine test_input7io_run


end subroutine test_inputio7
