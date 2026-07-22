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

!**s/r INs_Spn_send

      subroutine INs_Spn_send (F_ND,Minx,Maxx,Miny,Maxy,F_err)
      use, intrinsic :: iso_fortran_env
      use iso_c_binding
      use MiMd
      use omp_timing
      use IOs
      use INs
      implicit none

      integer, intent(IN) :: Minx,Maxx,Miny,Maxy,F_err
      real, intent(IN) :: F_ND(Minx:Maxx,Miny:Maxy,*)
      
      integer :: tag1,tag3
      integer :: i,j,k,cnt,n,n1,n2,offs,clients,nm,err
      integer i0,in,j0,jn                
!     
!--------------------------------------------------------------------
!
      if (Lun_out>0) call clock ( Lun_out, 'SPN Waitall', .false. )
      call gtmg_start ( 40, 'waitall', 23)
      call MPI_waitall (size(INs_Spn_isend),INs_Spn_isend,&
                        MPI_STATUSES_IGNORE,err)
      call gtmg_stop ( 40 )
      if (Lun_out>0) call clock ( Lun_out, 'SPN Waitall ,,, DONE', .false. )

      tag1 = 30001
      tag3 = 32001
      
      if (INs_1o1_L) then
         Spn_VGD_tbl(1:size(vtbl_8)) = reshape(vtbl_8,(/size(vtbl_8)/))
         do n=1,INs_nreq
            Spn_cBUF(         n) = SRL(n)%vname(1)
            Spn_cBuf(INs_nreq+n) = SRL(n)%vname(2)
         end do
         n=0
         Spn_iBUF(n+1) = INs_nreq
         Spn_iBUF(n+2) = INs_nplans
         n= n+2 ; nm=n
         Spn_iBUF(nm+1:nm+INs_nreq)= SRL(1:INs_nreq)%nk    ; nm=nm+INs_nreq
         Spn_iBUF(nm+1:nm+INs_nreq)= SRL(1:INs_nreq)%deb   ; nm=nm+INs_nreq
         Spn_iBUF(nm+1:nm+INs_nplans) = DIP1(1:INs_nplans) ; nm=nm+INs_nplans
         Spn_iBUF(nm+1:nm+3) = INs_n123 ; nm=nm+3
         call MPI_isend ( Spn_cBUF,size(Spn_cBUF)*len(Spn_cBUF(1)),MPI_CHARACTER,&
                          INs_gem1o1,tag1,INs_GEM_COMM,INs_Spn_isend(1),err)
         call MPI_isend ( Spn_iBUF,size(Spn_iBUF),MPI_INTEGER, INs_gem1o1,&
                          tag1+1,INs_GEM_COMM,INs_Spn_isend(2),err)
         call MPI_isend (Spn_VGD_tbl,size(Spn_VGD_tbl),MPI_DOUBLE_PRECISION,&
                         INs_gem1o1,tag1+2,INs_GEM_COMM,INs_Spn_isend(3),err)
      endif
      if (Lun_out>0) call clock ( Lun_out, 'SPN SEND_completed1', .false. )

      clients= 3
      call MPI_barrier (MY_WORLD_COMM,err)
      
      do n=client_pestart,client_peend
         offs=n-IOS_YIN*IOS_couleur-clients_npes(2,gem_id)+1
         i0= model_gindx(1,offs)-G_halox
         in= model_gindx(2,offs)+G_halox
         j0= model_gindx(3,offs)-G_haloy
         jn= model_gindx(4,offs)+G_haloy
         cnt=0
         do k=1,INs_nplans
         do j=model_gindx(3,offs)-G_haloy,model_gindx(4,offs)+G_haloy
         do i=model_gindx(1,offs)-G_halox,model_gindx(2,offs)+G_halox
            cnt= cnt+1
            Spnbuf(cnt,n)= F_ND(i,j,k)
         end do
         end do
         end do
         clients= clients+1
         cnt= ubound(Spnbuf,1)
         call MPI_isend ( Spnbuf(1,n),cnt,MPI_REAL,n,tag3+n,&
                          MiMd_gemworld,INs_Spn_isend(clients),err)
      end do

      if (Lun_out>0) call clock ( Lun_out, 'SPN SEND_completed2', .false. )
!     
!--------------------------------------------------------------------
!
      return
      end subroutine INs_Spn_send

