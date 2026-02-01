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

      subroutine spn_indata (F_datev_S)
      use lun
      use step_options
      use spn_options
      use svri_mod
      use vGrid_Descriptors
      use rmn_fst24
      use Openfst
      use inp_base
      use mem_nest
      use gmm_pw
      use ptopo
      use tr3d
      use, intrinsic :: iso_fortran_env
      implicit none
      
      character(len=*), intent(in):: F_datev_S

      include 'mpif.h'
      character(len=2048) dir
      character(len=16) :: Next_pilot_S
      logical :: success
      integer :: i,j,k,n, err, nka_tt, nka_uv, nka, kind, k0
      integer :: nm,n1,n2,n3,nreqs, nplans, indx_tt, indx_uv, dim, deb
      integer, dimension (:), pointer :: ip1_list,TT_ip1_list, UV_ip1_list
      real, dimension(:,:,:), pointer :: uur,vvr,tt,hu
      real(kind=REAL64), dimension(:), pointer :: wkpt8
      real(kind=REAL64), pointer :: vtbl_8(:,:,:)
      type(vgrid_descriptor) :: Spn_vgd_src
      type(fst_file) :: Spn_file
      real(kind=REAL64) :: dayfrac
      real(kind=REAL64), parameter :: one=1.0d0, sid=86400.0d0, rsid=one/sid
!     
!     ---------------------------------------------------------------
!
      if (.not. Spn_ON_L) return
      if (Lun_out > 0) write(lun_out,9000) trim(F_datev_S)
      dir='NUDGE/'//Step_runstrt_S(1:8)//Step_runstrt_S(10:11)
      nullify (ip1_list,TT_ip1_list, UV_ip1_list)
      nullify (uur,vvr,tt,hu,wkpt8)
      
      if (INs_server_L) then
         call gemtime (Lun_out, 'Input-svr: wait for SPN data', .false.)
         call MPI_waitall (size(INs_Spn_irecv),INs_Spn_irecv,&
                           MPI_STATUSES_IGNORE,err)
         call gemtime (Lun_out, 'Input-svr: Spn data received', .false.)
         dayfrac = Spn_yy_nudge_data_freq*rsid
         call incdatsd (Next_pilot_S, F_datev_S, dayfrac)
         if (Next_pilot_S<=Step_runend_S) then
            call itf_Iserv_request (Next_pilot_S,INs_Spnlist_S,INs_Spn_tag,INs_Spn_nrequests)
         endif
         dim= size(SPN_recv%cBUF)*len(SPN_recv%cBUF(1))
         call MPI_bcast (SPN_recv%cBUF,dim,MPI_CHARACTER, 0,&
                         COMM_multigrid, err)         
         call MPI_bcast (SPN_recv%iBUF,size(SPN_recv%iBUF),MPI_INTEGER, 0,&
                         COMM_multigrid, err)         
         call MPI_bcast (SPN_recv%VGD,size(SPN_recv%VGD),MPI_DOUBLE_PRECISION, 0,&
                         COMM_multigrid, err)
!   Spn_datev  = F_datev_S
         SPN_recv%nreq  = SPN_recv%iBUF(1)
         SPN_recv%nplans= SPN_recv%iBUF(2)
         nreqs= SPN_recv%nreq ; nplans= SPN_recv%nplans
         if (associated(SPN_recv%nk ) ) deallocate (SPN_recv%nk )
         if (associated(SPN_recv%deb) ) deallocate (SPN_recv%deb)
         allocate (SPN_recv%nk(nreqs),SPN_recv%deb(nreqs))
         nm= 2
         SPN_recv%nk (1:nreqs)= SPN_recv%iBUF(nm+1:nm+nreqs) ; nm=nm+nreqs
         SPN_recv%deb(1:nreqs)= SPN_recv%iBUF(nm+1:nm+nreqs) ; nm=nm+nreqs
         SPN_recv%ip1 =>   SPN_recv%iBUF(nm+1:nm+nplans) ; nm=nm+nplans
         n1= SPN_recv%iBUF(nm+1) ; n2= SPN_recv%iBUF(nm+2) ; n3= SPN_recv%iBUF(nm+3)
         indx_tt= -999 ; indx_uv= -999
         do n= 1, SPN_recv%nreq
            if (SPN_recv%cBUF(n+nreqs)(1:2) == 'TT') indx_tt= n
            if (SPN_recv%cBUF(n+nreqs)(1:2) == 'UV') indx_uv= n
         end do
         nka_tt= SPN_recv%nk(indx_tt)
         nka_uv= SPN_recv%nk(indx_uv)/2
!!$         allocate (TT_ip1_list(nka_tt),UV_ip1_list(nka_uv))
!!$         indx_tt= (indx_tt-1)*nka_tt+1
!!$         indx_uv= (indx_uv-1)*nka_uv+1
!!$         TT_ip1_list(1:nka_tt)= SPN_recv%ip1(indx_tt:indx_tt+nka_tt-1)
!!$         UV_ip1_list(1:nka_uv)= SPN_recv%ip1(indx_uv:indx_uv+nka_uv-1)
         allocate ( vtbl_8(n1,n2,n3) )
         vtbl_8 = reshape(SPN_recv%VGD,(/n1,n2,n3/))
         err= vgd_new ( Spn_vgd_src, vtbl_8 )
         deallocate (vtbl_8)
         
         !##### aiguillage pour inp_read_mt
         SRL(1:nreqs)%nk = SPN_recv%nk(1:nreqs)
         SRL(1:nreqs)%deb= SPN_recv%deb(1:nreqs)
         SRL(1:nreqs)%vname(1) = SPN_recv%cBUF(1:nreqs)
         SRL(1:nreqs)%vname(2) = SPN_recv%cBUF(nreqs+1:2*nreqs)
         INS_ND => SPN_recv%RBUF
         INs_DIP1 => SPN_recv%ip1
         INs_recv_nreqs= SPN_recv%nreq
!#######################

         do n=1,INs_recv_nreqs
            if ((trim(SRL(n)%vname(2)) == 'UVRT1' ) .or. &
                (trim(SRL(n)%vname(2)) == 'UV'    ) .and.&
                (SRL(n)%nk>0) ) then
               nka_uv= SRL(n)%nk/2
               deb= (SRL(n)%deb - 1) * INs_dimgzH + 1
               dim= INs_dimgzH*nka_uv
               allocate ( uur(l_minx:l_maxx,l_miny:l_maxy,nka_uv),&
                          vvr(l_minx:l_maxx,l_miny:l_maxy,nka_uv ))
               call reshapeH ( INS_ND(deb:),uur,&
                         l_minx,l_ni+G_halox,l_miny,l_nj+G_haloy,&
                         l_minx,l_maxx,l_miny,l_maxy,nka_uv)
               call reshapeH ( INS_ND(deb+dim:),vvr,&
                         l_minx,l_ni+G_halox,l_miny,l_nj+G_haloy,&
                         l_minx,l_maxx,l_miny,l_maxy,nka_uv)
               if (lun_out>0) write(6,'(a,a13,a,i4)') &
               ' I-svr FOUND: ','UV',' at indexe: ',n
            endif
         end do
      else
         call open_fst (F_datev_S, Spn_file, Spn_vgd_src, TT_ip1_list,&
                                       nka_tt, trim(dir),Spn_listfst_L)
         if (.not.associated(nest_now)) call nest_set_mem (max(G_nk,nka_tt))
         call inp_read_uv ( uur, vvr, 'UV' , UV_ip1_list, nka_uv, F_datev_S, Spn_file )
      endif
      
      if (.not.associated(nest_now)) call nest_set_mem (max(G_nk,nka_tt))
      err = vgd_get ( Spn_vgd_src, key='KIND',value=kind )
      Spn_nka= nka_tt
      if (allocated(Spn_pres)) deallocate (Spn_pres)
      allocate (Spn_pres(l_minx:l_maxx,l_miny:l_maxy,nka_tt))

      if (kind == 2) then
         err= vgd_get(Spn_vgd_src,'CA_M - vertical A coefficient (m)',wkpt8)
         do k=1, nka_tt
            Spn_pres(:,:,k) = log(wkpt8(k))
         end do
      else
         stop 'spn_indata: Incomplete code'
      endif

      err = inp_read ('HU', 'Q', hu, 1, ip1_list, nka, F_datev_S, Spn_file, F_type_S='P')
      deallocate (ip1_list) ;  nullify(ip1_list)
      
      err = inp_read ('TT', 'Q', tt, 1, ip1_list, nka, F_datev_S, Spn_file, F_type_S='P')

      k0= (Tr3d_hu-1)*Spn_nka
      do k= 1, nka
         do j=1-g_haloy, l_nj+g_haloy
            do i=1-g_halox, l_ni+g_halox
               nest_t_fin(i,j,k)=(tt(i,j,k)+TCDK)  * (1.0d0 + delta*hu(i,j,k))
               nest_tr_fin(i,j,k0+k)=hu(i,j,k)
               nest_u_fin(i,j,k)= uur(i,j,k)
               nest_v_fin(i,j,k)= vvr(i,j,k)
            end do
         end do
      end do
      call glbstat (nest_t_fin,'NTT','', l_minx,l_maxx,l_miny,l_maxy,1,nka,&
                    1,G_ni,1,G_nj,1,nka)
      call glbstat (nest_tr_fin(l_minx,l_miny,k0+1),'NHU','', l_minx,l_maxx,l_miny,l_maxy,1,nka,&
                    1,G_ni,1,G_nj,1,nka)
      call glbstat (nest_u_fin,'NUU','', l_minx,l_maxx,l_miny,l_maxy,1,nka,&
                    1,G_ni,1,G_nj,1,nka)
      call glbstat (nest_v_fin,'NVV','', l_minx,l_maxx,l_miny,l_maxy,1,nka,&
                    1,G_ni,1,G_nj,1,nka)
      deallocate (ip1_list,hu,tt) ;  nullify(ip1_list,hu,tt)
      deallocate (uur,vvr) ; nullify(uur,vvr)
      
      err = vgd_free(Spn_vgd_src)
    !  deallocate (TT_ip1_list,UV_ip1_list)
      if (.not. INs_server_L) success= Spn_file%close()
      INs_recv_nreqs= 0

!      call gtmg_stop ( 53 )
 9000 format(/,' TREATING SPN INPUT DATA VALID AT: ',a,&
      /,' ===============================================')
!     
!     ---------------------------------------------------------------
!
      return
      end subroutine spn_indata
