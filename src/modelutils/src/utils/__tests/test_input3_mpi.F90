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
subroutine test_input3_mpi()
   use iso_c_binding
   use testutils
   use fst_mod
   use clib_itf_mod
   use input_mod
   use incfg_mod
   use vGrid_Descriptors
   use vgrid_wb
   use ezgrid_mod
   use hinterp4yy_mod
   use rmn_gmm
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04
!@/
#include <rmnlib_basics.hf>
#include <rmn/WhiteBoard.hf>
   include "rpn_comm.inc"
   include "drv_dyn_itf.inc"
   integer,parameter :: NDIGITS = 4
   logical,parameter :: FSTOPC_SET = .false.
   character(len=512) :: dir_S,filename_S,dfiles_S,bcmk_S, &
        anal_S,anal2_S,clim_S,clim2_S,geop_S,geop2_S,inrep_S,inrep2_S
   integer :: istat,myproc
   ! ---------------------------------------------------------------------
   myproc = testutils_initmpi()
   call testutils_verbosity()
   call testutils_set_name('test_input3_mpi')
   istat = fstopc('MSGLVL','SYSTEM',FSTOPC_SET)
   call msg_set_p0only(0)

   write(dir_S,'(a,i4.4)') '__to-delete_test_input3_mpi-',myproc
   istat = clib_mkdir(trim(dir_S))
   istat = clib_chdir(trim(dir_S))

   anal_S = '/users/dor/armn/sch/Data/ords/data/sps-bil/DEBUG_SPSV113/ana/2015050518_001'
   anal2_S = './ANALYSIS' !'./analysis'
   istat = clib_unlink(trim(anal2_S))
   istat = clib_symlink(trim(anal_S),trim(anal2_S))
   if (.not.RMN_IS_OK(istat)) then
      call msg(MSG_WARNING,'problem creating symlink')
   endif

   call test_input3_get_data(myproc)

   istat = clib_unlink(trim(anal2_S))

   istat = clib_chdir('./..')
   istat = clib_unlink(trim(dir_S))

   call rpn_comm_barrier(RPN_COMM_WORLD, istat)
   call rpn_comm_finalize(istat)
   ! ---------------------------------------------------------------------
   return

contains

   !/@
   subroutine test_input3_get_data(myproc)
      implicit none
      !@/
      integer :: myproc
      character(len=*),parameter :: WB_GRID_SEC = 'grid_cfgs/'
      integer,parameter :: NSTEP = 1
      integer,parameter :: NVAR = 1
      integer,parameter :: NI0 = 578
      integer,parameter :: NJ0 = 428
      character(len=512) :: dateo_S,varname_S,varname2_S
      integer :: istat,inputid,dateo,istep,dt,gridid1,gridid2,ivar,ni,nj,halox,haloy,i0,j0,in,jn,ni2,nj2,nk2,k,fileid,key,gridid0,gridsetid
      integer :: nij1(2),nij2(2),ij01(2),ij02(2),icase
      real,pointer :: ax1(:,:),ay1(:,:),ax2(:,:),ay2(:,:)
      logical :: periodx,periody
      real, pointer :: data0(:,:,:),data1(:,:,:),data2(:,:,:)
      ! ---------------------------------------------------------------------
!!$      dateo_S = '20150505.180000' !#'20090427.000000'
!!$      call datp2f(dateo,dateo_S)
!!$      dt = 1800
!!$      inputid = input_new(dateo,dt)
!!$      istat = incfg_add_string(inputid, &
!!$           'in=I7;    freq=0; search=ANAL;      interp=linear;   levels= 1,3;')

      fileid = fst_open('/users/dor/armn/sch/Data/ords/data/sps-bil/DEBUG_SPSV113/ana/2015050518_001')
      dateo = -1
      key = fst_find(fileid,'I7',dateo,-1,-1,-1)
      nullify(data0)
      istat = fst_read(key,data0,fileid,gridid0)
      
      istat = wb_put(WB_GRID_SEC//'Grd_typ_S', 'LU')
      istat = wb_put(WB_GRID_SEC//'Grd_ni', NI0)
      istat = wb_put(WB_GRID_SEC//'Grd_nj', NJ0)
      istat = wb_put(WB_GRID_SEC//'Grd_nila', 0)
      istat = wb_put(WB_GRID_SEC//'Grd_njla', 0)
      istat = wb_put(WB_GRID_SEC//'Grd_maxcfl', 1)
      istat = wb_put(WB_GRID_SEC//'grd_dx', 0.0225)
      istat = wb_put(WB_GRID_SEC//'grd_dy', 0.0225)
      istat = wb_put(WB_GRID_SEC//'grd_dxmax', 360.)
      istat = wb_put(WB_GRID_SEC//'grd_dymax', 180.)
      istat = wb_put(WB_GRID_SEC//'Grd_iref', 269)
      istat = wb_put(WB_GRID_SEC//'Grd_jref', 345)
      istat = wb_put(WB_GRID_SEC//'Grd_latr',    2.75)
      istat = wb_put(WB_GRID_SEC//'Grd_lonr',  180.)
      istat = wb_put(WB_GRID_SEC//'Grd_xlon1',-100.)
      istat = wb_put(WB_GRID_SEC//'Grd_xlat1',  53.)
      istat = wb_put(WB_GRID_SEC//'Grd_xlon2', -85.)
      istat = wb_put(WB_GRID_SEC//'Grd_xlat2',  50.)
      istat = wb_put(WB_GRID_SEC//'Grd_overlap', 0.)
      istat = wb_put(WB_GRID_SEC//'Grd_gauss_L', .false.)
!!$      istat = wb_put('ptopo/ngrids', 1)
!!$      istat = wb_put('ptopo/igrid',  0)

      istat = dyn_grid_init(ni,nj,halox,haloy,periodx,periody,gridid1)

!!$         nullify(data1)
!!$         istat = input_setgridid(inputid,gridid1)
!!$         istat = input_get(inputid,ivar,istep,gridid1,data1)
      allocate(data1(ni,nj,1))
!!$      istat = ezsetopt('INTERP_DEGREE','linear')
!!$      gridsetid = ezdefset(gridid1, gridid0)
!!$      istat = ezsint(data1,data0)
      istat = hinterp4yy2d(data1,data0,1,gridid0,gridid1,gridid1,'linear','I7')

  
      do icase = 1,9
         
         select case(icase)
            !# Partition 2x2
         case(1) !# d1 != d2
            i0 = 1
            in = NI0/2
            j0 = 1
            jn = NJ0/2
         case(2) !# d1 != d2
            i0 = 1+NI0/2
            in = NI0
            j0 = 1
            jn = NJ0/2
         case(3) !# d1 == d2
            i0 = 1
            in = NI0/2
            j0 = 1+NJ0/2
            jn = NJ0
         case(4) !# d1 == d2
            i0 = 1+NI0/2
            in = NI0
            j0 = 1+NJ0/2
            jn = NJ0

            !# Partition 2x1
         case(5) !# d1 == d2
            i0 = 1
            in = NI0/2
            j0 = 1
            jn = NJ0
         case(6) !# d1 == d2
            i0 = 1+NI0/2
            in = NI0
            j0 = 1
            jn = NJ0


            !# Partition 1x2
          case(7) !# d1 != d2
            i0 = 1
            in = NI0
            j0 = 1
            jn = NJ0/2
         case(8) !# d1 == d2
            i0 = 1
            in = NI0
            j0 = 1+NJ0/2
            jn = NJ0

            !# center
         case(9) !# d1 == d2
            i0 = NI0/4
            j0 = NJ0/4
            in = NI0-i0
            jn = NJ0-j0
         end select
         print *,myproc,'C:',icase,i0,in,j0,jn

         gridid2 = ezgrid_sub(gridid1,i0,j0,in,jn)

         nullify(ax1,ay1,ax2,ay2)
         istat = ezgrid_params(gridid1,F_nij=nij1,F_ij0=ij01,F_ax=ax1,F_ay=ay1)
         istat = ezgrid_params(gridid2,F_nij=nij2,F_ij0=ij02,F_ax=ax2,F_ay=ay2)
         print *,myproc,': g1 :',nij1,ax1(i0,1),ay1(1,j0)
         print *,myproc,': g2 :',nij2,ax2(1,1),ay2(1,1)
         print *,myproc,': g2-1 :',minval(ax2-ax1(i0:in,:)),maxval(ay2-ay1(:,j0:jn))
         istep = 0
         ivar  = 1

!!$         nullify(data2)
!!$         istat = input_setgridid(inputid,gridid2)
!!$         istat = input_get(inputid,ivar,istep,gridid2,data2)
         allocate(data2(nij2(1),nij2(2),1))
!!$         istat = ezsetopt('INTERP_DEGREE','linear')
!!$         gridsetid = ezdefset(gridid2, gridid0)
!!$         istat = ezsint(data2,data0)
         istat = hinterp4yy2d(data2,data0,1,gridid0,gridid2,gridid2,'linear','I7')

         if (associated(data1) .and. associated(data2)) then
            !TODO: check grids
            ni2 = size(data2,1)
            nj2 = size(data2,2)
            nk2 = size(data2,3)
            do k=1,nk2
!!$            print *,myproc,':',k, 'D1   MIN,MAX:', minval(data1(i0:in,j0:jn,k)),maxval(data1(i0:in,j0:jn,k)),shape(data1),(/i0,in,j0,jn,k/)
!!$            print *,myproc,':',k, 'D2   MIN,MAX:', minval(data2(:,:,k)),maxval(data2(:,:,k)),shape(data2)
               print *,myproc,':',k, 'DIFF MIN,MAX:', &
                    minval(data2(1:ni2,1:nj2,k)-data1(i0:in,j0:jn,k)), &
                    maxval(data2(1:ni2,1:nj2,k)-data1(i0:in,j0:jn,k))
            enddo
         endif

         if (associated(data2)) deallocate(data2,stat=istat)

      enddo
      if (associated(data1)) deallocate(data1,stat=istat)

!!$      istat = input_close_files()
      ! ---------------------------------------------------------------------
      return
   end subroutine test_input3_get_data


   !/@
   subroutine test_input3b_get_data(myproc)
      implicit none
      !@/
      integer :: myproc
      character(len=*),parameter :: WB_GRID_SEC = 'grid_cfgs/'
      integer,parameter :: NSTEP = 1
      integer,parameter :: NVAR = 1
      integer,parameter :: NI0 = 578
      integer,parameter :: NJ0 = 428
      character(len=512) :: dateo_S,varname_S,varname2_S
      integer :: istat,inputid,dateo,istep,dt,gridid1,gridid2,ivar,ni,nj,halox,haloy,i0,j0,in,jn,ni2,nj2,nk2,k,l_minx,l_maxx,l_ni,G_lnimax,l_i0,l_miny,l_maxy,l_nj,G_lnjmax,l_j0,hx,hy
      integer :: nij1(2),nij2(2),ij01(2),ij02(2)
      real,pointer :: ax1(:,:),ay1(:,:),ax2(:,:),ay2(:,:)
      logical :: periodx,periody
      real, pointer :: data1(:,:,:),data2(:,:,:)
      ! ---------------------------------------------------------------------
      istat = wb_put(WB_GRID_SEC//'Grd_typ_S', 'LU')
      istat = wb_put(WB_GRID_SEC//'Grd_ni', NI0)
      istat = wb_put(WB_GRID_SEC//'Grd_nj', NJ0)
      istat = wb_put(WB_GRID_SEC//'Grd_nila', 0)
      istat = wb_put(WB_GRID_SEC//'Grd_njla', 0)
      istat = wb_put(WB_GRID_SEC//'Grd_maxcfl', 1)
      istat = wb_put(WB_GRID_SEC//'grd_dx', 0.0225)
      istat = wb_put(WB_GRID_SEC//'grd_dy', 0.0225)
      istat = wb_put(WB_GRID_SEC//'grd_dxmax', 360.)
      istat = wb_put(WB_GRID_SEC//'grd_dymax', 180.)
      istat = wb_put(WB_GRID_SEC//'Grd_iref', 269)
      istat = wb_put(WB_GRID_SEC//'Grd_jref', 345)
      istat = wb_put(WB_GRID_SEC//'Grd_latr',    2.75)
      istat = wb_put(WB_GRID_SEC//'Grd_lonr',  180.)
      istat = wb_put(WB_GRID_SEC//'Grd_xlon1',-100.)
      istat = wb_put(WB_GRID_SEC//'Grd_xlat1',  53.)
      istat = wb_put(WB_GRID_SEC//'Grd_xlon2', -85.)
      istat = wb_put(WB_GRID_SEC//'Grd_xlat2',  50.)
      istat = wb_put(WB_GRID_SEC//'Grd_overlap', 0.)
      istat = wb_put(WB_GRID_SEC//'Grd_gauss_L', .false.)
!!$      istat = wb_put('ptopo/ngrids', 1)
!!$      istat = wb_put('ptopo/igrid',  0)

      istat = dyn_grid_init(ni,nj,halox,haloy,periodx,periody,gridid1)

      istat = rpn_comm_topo(ni,l_minx,l_maxx,l_ni,G_lnimax, &
           halox,l_i0,RPN_COMM_TOPO_X,RPN_COMM_TOPO_FILL)
      istat = rpn_comm_topo(nj,l_miny,l_maxy,l_nj,G_lnjmax, &
           haloy,l_j0,RPN_COMM_TOPO_Y,RPN_COMM_TOPO_FILL)


      hx=halox ; hy=haloy
      i0=l_i0+hx ; j0=l_j0+hy
      in=l_i0+l_ni-1+hx ; jn=l_j0+l_nj-1+hy
      print *,myproc,':','X',ni,l_minx,l_maxx,l_ni,halox,l_i0
      print *,myproc,':','Y',nj,l_miny,l_maxy,l_nj,haloy,l_j0
      print *,myproc,':','XY',i0-hx,j0-hy,in+hx,jn+hy
      gridid2 = ezgrid_sub(gridid1,i0-hx,j0-hy,in+hx,jn+hy)

      nullify(ax1,ay1,ax2,ay2)
      istat = ezgrid_params(gridid1,F_nij=nij1,F_ij0=ij01,F_ax=ax1,F_ay=ay1)
      istat = ezgrid_params(gridid2,F_nij=nij2,F_ij0=ij02,F_ax=ax2,F_ay=ay2)
      print *,myproc,': g1 :',nij1,ax1(i0,1),ay1(1,j0)
      print *,myproc,': g2 :',nij2,ax2(1,1),ay2(1,1)
      print *,myproc,': g2-1 :',minval(ax2-ax1(i0:in,:)),maxval(ay2-ay1(:,j0:jn))

      dateo_S = '20150505.180000' !#'20090427.000000'
      call datp2f(dateo,dateo_S)
      dt = 1800

      inputid = input_new(dateo,dt)
      istat = incfg_add_string(inputid, &
           'in=I7;    freq=0; search=ANAL;      interp=linear;   levels= 1,3;')

      istep = 0
      ivar  = 1
      nullify(data1, data2)
      istat = input_setgridid(inputid,gridid1)
      istat = input_get(inputid,ivar,istep,gridid1,data1)

      istat = input_setgridid(inputid,gridid2)
      istat = input_get(inputid,ivar,istep,gridid2,data2)

      if (associated(data1) .and. associated(data2)) then
         !TODO: check grids
         ni2 = size(data2,1)
         nj2 = size(data2,2)
         nk2 = size(data2,3)
         do k=1,nk2
!!$            print *,k, 'D1   MIN,MAX:', minval(data1(i0:in,j0:jn,k)),maxval(data1(i0:in,j0:jn,k)),shape(data1),(/i0,in,j0,jn,k/)
!!$            print *,k, 'D2   MIN,MAX:', minval(data2(:,:,k)),maxval(data2(:,:,k)),shape(data2)
            print *,myproc,':',k, 'DIFF MIN,MAX:', &
                 minval(data2(1:ni2,1:nj2,k)-data1(i0:in,j0:jn,k)), &
                 maxval(data2(1:ni2,1:nj2,k)-data1(i0:in,j0:jn,k))
         enddo
      endif

      if (associated(data1)) deallocate(data1,stat=istat)
      if (associated(data2)) deallocate(data2,stat=istat)

      istat = input_close_files()
      ! ---------------------------------------------------------------------
      return
   end subroutine test_input3b_get_data


end subroutine test_input3_mpi
