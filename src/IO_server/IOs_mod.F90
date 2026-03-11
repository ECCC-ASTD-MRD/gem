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

module IOs
   use, intrinsic :: iso_fortran_env
   use vGrid_Descriptors
   use rmn_fst24
   implicit none
   public
   save

      include 'mpif.h'
      
      character(len=256 ) :: IOs_component_S
      character(len=2048) :: IOs_path_input_S, IOs_path_work_S
      logical :: IOs_server_L=.true.
      logical :: IOs_1o1_L, IOs_samenode_L, Grd_yinyang_L
      logical :: Tr3d_anydate_L
      integer :: MY_WORLD_COMM, YYG_COMM
      integer :: gem_id, me_1on1, gem_1on1, Lun_out
      integer :: myproc_GLB,myproc_IOS,numproc_IOS
      integer :: IOs_color, IOs_ncolors
      integer :: YYG_numproc, YYG_myproc, IOS_couleur, IOS_YIN
      integer :: client_pelocal,client_pestart,client_peend
      integer :: G_ni,G_nj,G_nk,G_halox,G_haloy,client_Hplane
      integer :: l_minx,l_maxx,l_miny,l_maxy,Tr3d_ntr
      integer :: Glb_pil_w,Glb_pil_e,Glb_pil_s,Glb_pil_n
      integer :: Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4
      integer :: RotY_ig1, RotY_ig2, RotY_ig3, RotY_ig4
      integer, dimension (:,:), allocatable :: model_gindx,servicing
      integer, dimension (:,:), allocatable :: clients_npes,host_pe0s
      real  :: Grd_xlat1,Grd_xlon1,Grd_xlat2,Grd_xlon2
      real  :: Grd_xlat1Y,Grd_xlon1Y,Grd_xlat2Y,Grd_xlon2Y
      real, dimension(:), pointer :: Level_allpres, Level_allheights
      real, dimension(:), pointer, contiguous :: geomh_latgs
      real, dimension(:), pointer, contiguous :: geomh_longs
      real, dimension(:), pointer, contiguous :: geomh_latgv
      real, dimension(:), pointer, contiguous :: geomh_longu
      real, dimension(:), pointer :: model_Hgeom
      real, dimension(:), pointer, contiguous :: Ver_hybM,Ver_hybT,Ver_i
      real, dimension(:), pointer, contiguous :: hybM_diag,hybT_diag
      real, dimension(:,:,:,:), pointer :: IOs_glbdata,glbdata
      type(vgrid_descriptor) :: gem_vgd
      type(fst_file) :: Nes_fst(10)

contains

      subroutine IOs_mpi_init( F_component_S, F_colors, F_COMMs, &
                               F_nc,F_server_L)
      use MiMd
      use omp_timing
      use clib_itf_mod
      implicit none

#include <rmnlib_basics.hf>
      character(len=*), intent(IN) :: F_component_S
      logical, intent(OUT) :: F_server_L
      integer, intent(IN ) :: F_nc
      integer, intent(IN ) :: F_colors(F_nc)
      integer, intent(OUT) :: F_COMMs(F_nc)

      character(len=256 ) :: string_S
      integer i,mpx,irest, wnum,err
!     
!--------------------------------------------------------------------
!
      IOs_component_S= F_component_S ; IOs_color= F_colors(1)
      call MiMd_init (F_component_S,F_colors,F_nc,F_COMMs,&
                      wnum, myproc_GLB, numproc_IOS, myproc_IOS)
      
      MY_WORLD_COMM = F_COMMs(1)

      Lun_out= -1
      if ( myproc_IOS == 0 ) then
         Lun_out= 6
         err= exdb(trim(F_component_S)//' - SERVER','','NON')
      endif
      
      call clock ( 6, ' ', .false. )
      string_S='STARTING '//trim(F_component_S)
      if (Lun_out>0) call clock ( Lun_out, trim(string_S), .false. )
      call gtmg_init (myproc_IOS, MY_WORLD_COMM, trim(F_component_S))
      
! The WORLD by colors...
      string_S=trim(F_component_S)//': The WORLD by colors'
      if (myproc_IOS == 0) write(6,'(xa)') trim(string_S)
      allocate (clients_npes(3,MiMd_ncolors))
      do i=1,MiMd_ncolors
         if (trim(MiMd_world(i)%name_S)=="GEMDM") gem_id=i
         clients_npes(1,i) = MiMd_world(i)%rank
         clients_npes(2,i) = MiMd_world(i)%wpe0
         clients_npes(3,i) = clients_npes(2,i) + clients_npes(1,i) - 1
         if (myproc_IOS == 0) &
         write(6,1001) &
             trim(MiMd_world(i)%name_S),MiMd_world(i)%color,&
             MiMd_world(i)%wpe0,MiMd_world(i)%wpe0+MiMd_world(i)%rank-1
      end do

      call gtmg_start ( 1, trim(F_component_S), 0)
      IOs_samenode_L = same_host() == 0
      if (IOs_samenode_L) then
         call splitW ( myproc_IOS, numproc_IOS,clients_npes(1,gem_id),&
                       clients_npes(2,gem_id),client_pelocal         ,&
                       client_pestart,client_peend)
      else
         IOs_server_L = .false.
         if (Lun_out>0) &
          write(Lun_out,'(/3x,10("#")/3x,2a/3x,a/3x,10("#")/)')&
               trim(F_component_S),' is disconnecting because',&
                                'NOT all PEs on the same node'
      endif
      F_server_L = IOs_server_L
      err = clib_getenv ('PWD', IOs_path_work_S)
      
 1001 format (" Name: ",a20,", Color:",i3,", World PEs: ",i5," to ",i5)
!     
!--------------------------------------------------------------------
!
      return
      end subroutine IOs_mpi_init

      integer function same_host ()
      use iso_c_binding
      implicit none

      interface
         integer(C_INT) function f_gethostid()BIND(C,name='gethostid')
         import :: C_INT
         end function f_gethostid
      end interface

!      character(len=MPI_MAX_PROCESSOR_NAME) :: host
      integer me,err,cnt,host_list(numproc_IOS),glist(numproc_IOS)
!     
!--------------------------------------------------------------------
!
      host_list= 0 ; me= myproc_IOS+1 ; same_host= 0
      host_list(me) = abs(f_gethostid())
      call MPI_ALLREDUCE ( host_list, glist, numproc_IOS, &
                     MPI_INTEGER,MPI_SUM,MY_WORLD_COMM,err)
      if (host_list(me)/=glist(1)) then
         host_list(me)=0
      else
         host_list(me)=1
      endif
      call MPI_ALLREDUCE ( host_list, glist, numproc_IOS, &
                     MPI_INTEGER,MPI_SUM,MY_WORLD_COMM,err)
      if ( any(glist<1) ) same_host= -1
!      call MPI_GET_PROCESSOR_NAME(host, cnt, err)
!     
!--------------------------------------------------------------------
!
      return
      end function same_host

      subroutine IOs_gem_init (F_COMM,F_1o1_L,F_me1o1,F_gem1o1,F_grid_info)
      use iso_c_binding
      implicit none

      logical, intent(OUT) :: F_1o1_L
      integer, intent(OUT) :: F_me1o1,F_gem1o1,F_grid_info(100)
      integer, intent(IN ) :: F_COMM

      character(len=4) :: dumc
      logical :: bool
      integer(KIND=MPI_ADDRESS_KIND) :: WINSIZE
      integer :: k,i,myproc,ierr,dim,DISP_UNIT
      integer :: dimens,dimx,dimy,hostcomm,err,recv(3)
      integer :: tag,ipcode,ipkind,npe_per,mpx,irest
      integer, dimension(3,numproc_IOS) :: service_indx
      type(C_PTR), save :: basepntr
      real(kind=REAL64), pointer :: vtbl_8(:,:,:)
      real, dimension(:), pointer :: wkpt
      real :: pcode
!     
!--------------------------------------------------------------------
!
      IOs_1o1_L = F_COMM .ne. MPI_COMM_NULL

      if (IOs_1o1_L) then 
         call MPI_barrier   (F_COMM, err)
         call MPI_COMM_rank (F_COMM, me_1on1,ierr)
         gem_1on1 = -1*(me_1on1-1)
         tag = 110
         call MPI_send ( IOs_server_L,1,MPI_LOGICAL,&
                         gem_1on1,tag,F_COMM,ierr )
         if (IOs_server_L) then
            tag=111
            call MPI_recv ( F_grid_info, size(F_grid_info),MPI_INTEGER, &
                        gem_1on1, tag, F_COMM, MPI_STATUSES_IGNORE, ierr)
         endif
      endif

      if (.not.IOs_server_L) return

      F_1o1_L = IOs_1o1_L
      F_me1o1 = me_1on1
      F_gem1o1= gem_1on1

      if (.not.IOs_samenode_L) return
      
      call MPI_bcast (F_grid_info, size(F_grid_info), MPI_INTEGER, 0,&
                      MY_WORLD_COMM, ierr)

      G_ni          = F_grid_info(1)
      G_nj          = F_grid_info(2)
      G_nk          = F_grid_info(3)
      G_halox       = F_grid_info(4)
      G_haloy       = F_grid_info(5)
      l_minx        = F_grid_info(6)
      l_maxx        = F_grid_info(7)
      l_miny        = F_grid_info(8)
      l_maxy        = F_grid_info(9)
      Tr3d_ntr      = F_grid_info(10)
      Tr3d_anydate_L=TRANSFER(F_grid_info(11), bool)
      Glb_pil_w     = F_grid_info(12)
      Glb_pil_e     = F_grid_info(13)
      Glb_pil_s     = F_grid_info(14)
      Glb_pil_n     = F_grid_info(15)
      client_Hplane = F_grid_info(16)
      IOS_ncolors   = F_grid_info(17)
      allocate (vtbl_8(F_grid_info(18),F_grid_info(19),F_grid_info(20)))
      allocate (Level_allpres(F_grid_info(21)))
      allocate (Level_allheights(F_grid_info(22)))
      allocate (model_Hgeom(F_grid_info(23)))
      allocate (model_gindx(6,clients_npes(1,gem_id)/IOS_ncolors))

      if (IOs_1o1_L) then 
         call MPI_recv ( model_gindx, size(model_gindx), MPI_INTEGER ,&
                     gem_1on1, tag+1, F_COMM, MPI_STATUSES_IGNORE, ierr)
         call MPI_recv ( model_Hgeom, size(model_Hgeom), MPI_REAL    ,&
                     gem_1on1, tag+2, F_COMM, MPI_STATUSES_IGNORE, ierr)
         call MPI_recv ( Level_allpres, size(Level_allpres), MPI_REAL,&
                     gem_1on1, tag+3, F_COMM, MPI_STATUSES_IGNORE, ierr)
         call MPI_recv ( Level_allheights, size(Level_allheights), MPI_REAL,&
                     gem_1on1, tag+4, F_COMM, MPI_STATUSES_IGNORE, ierr)
         call MPI_recv ( vtbl_8, size(vtbl_8), MPI_DOUBLE_PRECISION  ,&
                     gem_1on1, tag+5, F_COMM, MPI_STATUSES_IGNORE, ierr)
         call MPI_recv ( IOs_path_input_S, len(IOs_path_input_S), &
                         MPI_CHARACTER,gem_1on1, tag+6, F_COMM  , &
                         MPI_STATUSES_IGNORE, ierr)
      endif

      call MPI_bcast (model_gindx,size(model_gindx),MPI_INTEGER, 0,&
                      MY_WORLD_COMM, ierr)
      call MPI_bcast (model_Hgeom,size(model_Hgeom),MPI_REAL   , 0,&
                      MY_WORLD_COMM, ierr)
      call MPI_bcast (Level_allpres,size(Level_allpres),MPI_REAL,0,&
                      MY_WORLD_COMM, ierr)
      call MPI_bcast (Level_allheights,size(Level_allheights),MPI_REAL,0,&
                      MY_WORLD_COMM, ierr)
      call MPI_bcast (vtbl_8, size(vtbl_8),MPI_DOUBLE_PRECISION ,0,&
                      MY_WORLD_COMM, ierr)
      call MPI_bcast (IOs_path_input_S, len(IOs_path_input_S), &
                      MPI_CHARACTER , 0, MY_WORLD_COMM, ierr)

      Grd_xlat1= model_Hgeom(1)
      Grd_xlon1= model_Hgeom(2)
      Grd_xlat2= model_Hgeom(3)
      Grd_xlon2= model_Hgeom(4)
      call cxgaig ( 'E',Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4,&
                    Grd_xlat1,Grd_xlon1,Grd_xlat2,Grd_xlon2 )
      Grd_yinyang_L= .false. ; IOS_couleur= 0
      if ( IOS_ncolors == 2 ) then !Yin-Yang
         Grd_yinyang_L= .true.
         Grd_xlat1Y= model_Hgeom(5)
         Grd_xlon1Y= model_Hgeom(6)
         Grd_xlat2Y= model_Hgeom(7)
         Grd_xlon2Y= model_Hgeom(8)
         call cxgaig ( 'E',RotY_ig1, RotY_ig2, RotY_ig3, RotY_ig4,&
                       Grd_xlat1Y,Grd_xlon1Y,Grd_xlat2Y,Grd_xlon2Y )
         if (numproc_IOS<2) then
            print*, &
            'ABORT - numproc_IOS must be > 1 in Yin-Yang config'
            stop ' ABORT in IOs_gem_init'
          endif
          npe_per   = numproc_IOS/IOS_ncolors
         IOS_couleur= min(myproc_IOS  / npe_per, 1)
         call MPI_Comm_split (MY_WORLD_COMM, IOS_couleur, &
                              numproc_IOS, YYG_COMM, ierr)
         call MPI_COMM_size  (YYG_COMM,YYG_numproc,ierr)
         call MPI_COMM_rank  (YYG_COMM,YYG_myproc ,ierr)
      endif
      
! Re-splitting the MPI_Rget work if Grd_yinyang_L
      if (Grd_yinyang_L) then
         IOS_YIN= clients_npes(1,gem_id)/IOS_ncolors
         call splitW ( YYG_myproc, YYG_numproc, IOS_YIN         ,&
                       clients_npes(2,gem_id)+IOS_YIN*IOS_couleur,&
                       client_pelocal,client_pestart,client_peend)
      endif
      call MPI_comm_rank (MPI_COMM_WORLD,myproc,err)
      allocate (servicing(3,numproc_IOS))
      service_indx=0
      service_indx(1,myproc_IOS+1) = myproc
      service_indx(2,myproc_IOS+1) = client_pestart
      service_indx(3,myproc_IOS+1) = client_peend
      dimens = 3*numproc_IOS
      call MPI_ALLREDUCE (service_indx,servicing,dimens,MPI_INTEGER,&
                          MPI_BOR,MY_WORLD_COMM,err)
      if (IOs_1o1_L) then
         tag= 311          
         call MPI_send (servicing,size(servicing),MPI_INTEGER,&
                        gem_1on1,tag  ,F_COMM, err )
         call MPI_barrier (F_COMM, err)
      endif
      print*, trim(IOs_component_S),' PE #',myproc_IOS,&
              'SERVICING GEM PEs # ',&
               client_pestart,' to ',client_peend
      dimx= G_ni+2*G_halox ; dimy= G_nj+2*G_haloy
      geomh_latgs(1-G_haloy:G_nj+G_haloy) => model_Hgeom(            9:)
      geomh_longs(1-G_halox:G_ni+G_halox) => model_Hgeom(dimy       +9:)
      geomh_latgv(1-G_haloy:G_nj+G_haloy) => model_Hgeom(dimx+  dimy+9:)
      geomh_longu(1-G_halox:G_ni+G_halox) => model_Hgeom(dimx+2*dimy+9:)
      ierr= vgd_new ( gem_vgd, vtbl_8 )
      deallocate (vtbl_8)

      allocate (Ver_hybM (G_nk+2), Ver_hybT (G_nk+2)) ; nullify(wkpt)
      allocate (hybM_diag(G_nk+1), hybT_diag(G_nk+1))
      ierr = vgd_get(gem_vgd,'VCDM - vertical coordinate (m)'   ,wkpt)
      Ver_hybM = wkpt(1:size(Ver_hybM)); deallocate(wkpt); nullify(wkpt)
      ierr = vgd_get(gem_vgd,'VCDT - vertical coordinate (t)'   ,wkpt)
      Ver_hybT = wkpt(1:size(Ver_hybT)); deallocate(wkpt); nullify(wkpt)
      hybM_diag(1:G_nk)= Ver_hybM(1:G_nk) ; hybM_diag(G_nk+1)= Ver_hybM(G_nk+2)
      hybT_diag(1:G_nk)= Ver_hybT(1:G_nk) ; hybT_diag(G_nk+1)= Ver_hybT(G_nk+2)
      allocate (Ver_i (2*G_nk))
      do i=1, 2*G_nk
         Ver_i(i) = i
      end do
            
!!$      if (IOs_1o1_L) then 
!!$         tag= 10 ; allocate (host_pe0s(clients_npes(1,gem_id),2))
!!$         call MPI_recv ( host_pe0s,size(host_pe0s),MPI_INTEGER,& 
!!$                         gem_1on1, tag, F_COMM, MPI_STATUSES_IGNORE, ierr)
!!$         do i=1, clients_npes(1,gem_id)
!!$            if (host_pe0s(i,1) >= 0) then
!!$               print*, 'Client node serv pe: ',host_pe0s(i,:)
!!$            else
!!$               exit
!!$            endif
!!$         end do
!!$         call MPI_barrier (F_COMM, err)
!!$      endif
!!$      call MPI_barrier (MY_WORLD_COMM,err)
      
      call convip ( ipcode, pcode, ipkind, 0, ' ', .false. )
!     
!--------------------------------------------------------------------
!
      return
      end subroutine IOs_gem_init
      
end module IOs
