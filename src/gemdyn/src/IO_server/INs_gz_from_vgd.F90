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

!**s/r INs_gz_from_vgd - Obtain 3D GZ from VGD descriptor

      logical function INs_gz_from_vgd ()
      use INs_base
      implicit none

      logical :: sleve_L
      character(len=1) :: grid_S(3)
      character(len=6) vname
      integer :: i,j,k, err(6), err_me, err_mels,&
                 src_nka, cnt, n3, ip1_list(5)  ,&
                 local_nk, kstart, kend
      integer, dimension (:), pointer :: ip1_m, ip1_t
      real :: src_rcoef(4)
      real, dimension (:,:,:), pointer :: src_me, src_mels
      real(kind=REAL64), dimension(:), pointer :: src_hybmA,src_hybtA,&
                         src_hybmB,src_hybtB,src_hybmC,src_hybtC,wkpt8
!
!---------------------------------------------------------------------
!
      INs_gz_from_vgd = .false.

      grid_S=['Q','U','V']

      src_rcoef = -1. ; err= -1
      nullify(ip1_m, ip1_t, src_me, src_mels)
      err(1) = vgd_get ( Inp_vgd_src, key='RC_1 - first R-coef value' ,value=src_rcoef(1))
      err(2) = vgd_get ( Inp_vgd_src, key='RC_2 - second R-coef value',value=src_rcoef(2))
      err(3) = vgd_get ( Inp_vgd_src, key='RC_3 - third R-coef value' ,value=src_rcoef(3))
      err(4) = vgd_get ( Inp_vgd_src, key='RC_4 - forth R-coef value' ,value=src_rcoef(4))
      err(5) = vgd_get ( Inp_vgd_src, key='VIPM - level ip1 list (m)' ,value=ip1_m)
      err(6) = vgd_get ( Inp_vgd_src, key='VIPT - level ip1 list (t)' ,value=ip1_t)
      
      if ( minval(err) < 0) goto 988
      sleve_L = minval(src_rcoef(3:4)) >=0.
      src_me  (1:INs_nid, 1:INs_njd, 1:3) => NES_DATA(1:)
      src_mels(1:INs_nid, 1:INs_njd, 1:3) => NES_DATA(3*INs_hord+1:)
      err(1) = INs_read_mt ('ME',grid_S, vname, src_me, 3, ip1_list,&
                          src_nka,Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4)
      if ( err(1) < 0 ) then
         if (Lun_out > 0) write(lun_out,9001) 'ME'
         goto 988
      endif
      GZcaract(3)= 1
      if (sleve_L) then
         err(1) = INs_read_mt ('MELS',grid_S, vname, src_mels, 3,&
             ip1_list, src_nka,Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4)
         if ( err(1) < 0) then
            if (Lun_out > 0) write(lun_out,9001) 'MELS'
            goto 988
         endif
         GZcaract(4)= 1
      else
         src_mels= 0.
      endif
      
      err(1) = vgd_get ( Inp_vgd_src, key='NL_M - number of momentum levels' ,value=src_nka)
      if ( err(1) < 0 ) then
         if (Lun_out > 0) write(lun_out,9001) 'VGD:NL_M'
         goto 988
      endif
      allocate ( src_hybmA(src_nka), src_hybtA(src_nka),&
                 src_hybmB(src_nka), src_hybtB(src_nka),&
                 src_hybmC(src_nka), src_hybtC(src_nka) )
      nullify(wkpt8)
      err(1)= vgd_get ( Inp_vgd_src, key='CA_M - vertical A coefficient (m)' ,value=wkpt8)
      if (err(1)==0)  src_hybmA= wkpt8(1:src_nka); deallocate(wkpt8); nullify(wkpt8)
      err(2)= vgd_get ( Inp_vgd_src, key='CB_M - vertical B coefficient (m)' ,value=wkpt8)
      if (err(2)==0)  src_hybmB= wkpt8(1:src_nka); deallocate(wkpt8); nullify(wkpt8)
      err(3)= vgd_get ( Inp_vgd_src, key='CC_M - vertical C coefficient (m)' ,value=wkpt8)
      if (err(3)==0)  src_hybmC= wkpt8(1:src_nka); deallocate(wkpt8); nullify(wkpt8)
      err(4)= vgd_get ( Inp_vgd_src, key='CA_T - vertical A coefficient (t)' ,value=wkpt8)
      if (err(4)==0)  src_hybtA= wkpt8(1:src_nka); deallocate(wkpt8); nullify(wkpt8)
      err(5)= vgd_get ( Inp_vgd_src, key='CB_T - vertical B coefficient (t)' ,value=wkpt8)
      if (err(5)==0)  src_hybtB= wkpt8(1:src_nka); deallocate(wkpt8); nullify(wkpt8)
      err(6)= vgd_get ( Inp_vgd_src, &
                        key='CC_T - vertical C coefficient (t)',value=wkpt8)
      if (err(6)==0)  src_hybtC= wkpt8(1:src_nka); deallocate(wkpt8); nullify(wkpt8)
      if ( minval(err(1:6) ) < 0) then
         if (Lun_out > 0) write(lun_out,9001) 'VGD:coef'
         goto 988
      endif
      src_nka= src_nka-2
      n3= 2*src_nka+1
      cnt=1
      do k=1,src_nka
         GZIP1(cnt  ) = ip1_m(k)
         GZIP1(cnt+1) = ip1_t(k)
         cnt= cnt+2
      end do
      GZIP1(n3) = ip1_m(src_nka+1)

      call splitW ( myproc_IOS, numproc_IOS,src_nka,1,&
                    local_nk, kstart, kend)
      do k=kstart, kend
         cnt= (k-1)*2 + 1
        ! write(6,'(a,7i5)') 'haha2: ',k,cnt,cnt+1,n3+cnt,n3+cnt+1,2*n3+cnt ,2*n3+cnt +1
         do j=lbound(src_me,2), ubound(src_me,2)
            do i=lbound(src_me,1), ubound(src_me,1)
              GZ(i,j,     cnt  )= src_hybmA(k)+ (src_hybmB(k)*src_me(i,j,1)+src_hybmC(k)*src_mels(i,j,1))
              GZ(i,j,     cnt+1)= src_hybtA(k)+ (src_hybtB(k)*src_me(i,j,1)+src_hybtC(k)*src_mels(i,j,1))
              GZ(i,j,  n3+cnt  )= src_hybmA(k)+ (src_hybmB(k)*src_me(i,j,2)+src_hybmC(k)*src_mels(i,j,2))
              GZ(i,j,  n3+cnt+1)= src_hybtA(k)+ (src_hybtB(k)*src_me(i,j,2)+src_hybtC(k)*src_mels(i,j,2))
              GZ(i,j,2*n3+cnt  )= src_hybmA(k)+ (src_hybmB(k)*src_me(i,j,3)+src_hybmC(k)*src_mels(i,j,3))
              GZ(i,j,2*n3+cnt+1)= src_hybtA(k)+ (src_hybtB(k)*src_me(i,j,3)+src_hybtC(k)*src_mels(i,j,3))
            end do
         end do
      end do
      if (myproc_IOS==0) then
         do j=lbound(src_me,2), ubound(src_me,2)
         do i=lbound(src_me,1), ubound(src_me,1)
           GZ (i,j,  n3)= src_me  (i,j,1)
           GZ (i,j,2*n3)= src_me  (i,j,2)
           GZ (i,j,3*n3)= src_me  (i,j,3)
           GZ (i,j,3*n3+1)= src_mels  (i,j,1)
         end do
         end do
      endif
      call MPI_barrier (MY_WORLD_COMM,err)

      GZcaract(1)= n3
      GZcaract(2)= 21
      GZcaract(5)= -huge(1)
      deallocate (src_hybmA,src_hybtA,src_hybmB,src_hybtB,src_hybmC,src_hybtC,ip1_m,ip1_t)
      nullify(src_hybmA,src_hybtA,src_hybmB,src_hybtB,src_hybmC,src_hybtC,wkpt8)
      INs_gz_from_vgd = .true.
      if (Lun_out > 0) write(lun_out,9000) src_nka,'VGD',Inp_vgdkind

 9000 format(1x,'INS_GZ_FROM_VGD: ', i3,' levels of ',a,' found in input file with kind=',i5)
 9001 format(1x,'INS_GZ_FROM_VGD: unable to retrieve ', a)
!
!---------------------------------------------------------------------
!
 988   return
       end function INs_gz_from_vgd
