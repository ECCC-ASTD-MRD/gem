      subroutine fislh_pres (F_q, F_ps, F_vt,&
                             Minx,Maxx,Miny,Maxy,F_nk)
      use glb_ld
      use dyn_fisl_options
      use gmm_geof
      use inp_base, only: inp_3dhgts
      use ver
      use tdpack
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in)  :: Minx,Maxx,Miny,Maxy, F_nk
      real, dimension(Minx:Maxx,Miny:Maxy,F_nk+1), intent(out) :: F_q
      real, dimension(Minx:Maxx,Miny:Maxy,F_nk)  , intent(in)  :: F_vt
      real, dimension(Minx:Maxx,Miny:Maxy     )  , intent(in)  :: F_ps

      ! Local varibales
      integer :: i, j, k, err
      integer, dimension(:) , pointer :: ip1_list
      real, dimension(:,:  ), pointer :: topo,topols
      real, dimension(:,:,:), pointer :: zz
      real(kind=REAL64) :: aaa
!
!     ________________________________________________________________
!
      nullify(zz,ip1_list)

      err= vgd_get ( Ver_vgdobj, 'VIPM - level ip1 list (m)', ip1_list )

      allocate (topo  (l_minx:l_maxx,l_miny:l_maxy),&
                topols(l_minx:l_maxx,l_miny:l_maxy),&
                zz(l_minx:l_maxx,l_miny:l_maxy,size(ip1_list)))

      aaa=rgasd_8*Cstv_tstr_8
      do j=1,l_nj
      do i=1,l_ni
         topo  (i,j)= fis0(i,j) / grav_8
         topols(i,j)= sls (i,j) / grav_8
         F_q(i,j,G_nk+1)=aaa*log(F_ps(i,j)/1.e5)
      end do
      end do

      call inp_3dhgts ( Ver_vgdobj, ip1_list, topo, topols, zz, 1, size(ip1_list))
      deallocate (topo,topols)

      ! Integrate hydrostatic equation, dq/dz=-gTstar/Tv, to obtain q at momentum levels
      aaa = grav_8*Cstv_tstr_8
      do k=G_nk,1,-1
         do j=1,l_nj
         do i=1,l_ni
            F_q(i,j,k) = F_q(i,j,k+1) + aaa*(zz(i,j,k+1) - zz(i,j,k))/F_vt(i,j,k)
         end do
         end do
      end do
      ! Add gz to obtain qprime, stored in F_q
      do k=1,G_nk+1
      do j=1,l_nj
      do i=1,l_ni
         F_q(i,j,k)=F_q(i,j,k)+grav_8*zz(i,j,k)
      end do
      end do
      end do
      deallocate(zz,ip1_list)
!
!     ________________________________________________________________
!
      return
      end subroutine fislh_pres
