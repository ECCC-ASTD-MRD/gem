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

module nest_blending

  implicit none
#include <arch_specific.hf>
  private
  public :: nest_blend

  ! This module blends a pair of GMM variables, 2D or 3D fields using a
  ! set of pre-determined blending weights (which themselves may vary in
  ! 3D).  When 2D fields are given, the optional 'level' argument can
  ! be used to specify which level of the blending weight 3D field should
  ! be applied (default is G_nk+2, the bottom level blending weights).

  interface nest_blend
     module procedure blend_gmm
     module procedure blend_array3d
     module procedure blend_array2d
  end interface


contains

  subroutine blend_gmm ( F_fld_S, F_nest_fld_S, F_grid_S )
    use glb_ld
    use gmm_itf_mod
    implicit none

    ! Apply blending to a GMM field.


    ! Arguments
    character(len=*), intent(in) :: F_fld_S       !GMM name of model   field
    character(len=*), intent(in) :: F_nest_fld_S  !GMM name of nesting field
    character(len=*), intent(in) :: F_grid_S      !Arakawa C grid ('U':u-grid,'V':v-grid,'M':mass-grid)

    ! Local variables
    integer :: Minx,Maxx,Miny,Maxy,Minz,Maxz,nk,gmmstat
    real, dimension(:,:  ), pointer :: fld2d=>null(), fld_nest2d=>null()
    real, dimension(:,:,:), pointer :: fld3d=>null(), fld_nest3d=>null()
    type(gmm_metadata) :: metadata
  !
  !----------------------------------------------------------------------
  !
    gmmstat = gmm_getmeta (F_fld_S, metadata)
    Minx = metadata%l(1)%low
    Maxx = metadata%l(1)%high
    Miny = metadata%l(2)%low
    Maxy = metadata%l(2)%high
    Minz = metadata%l(3)%low
    Maxz = metadata%l(3)%high
    if ((Minz == 0) .and. (Maxz == 0) .and. (metadata%l(3)%n == 0)) then
       Minz = 1
       Maxz = 1
    end if
    nk = Maxz - Minz + 1

    ! Split to handle either 2d or 3d GMM fields
    levels: if (nk == 1) then

       gmmstat = gmm_get (F_fld_S, fld2d, metadata)
       gmmstat = gmm_get (F_nest_fld_S, fld_nest2d, metadata)

       call blend_array2d (fld2d, fld_nest2d, Minx,Maxx,Miny,Maxy, F_grid_S, level=G_nk+1)

    else

       gmmstat = gmm_get (F_fld_S, fld3d, metadata)
       gmmstat = gmm_get (F_nest_fld_S, fld_nest3d, metadata)

       call blend_array3d (fld3d, fld_nest3d, Minx,Maxx,Miny,Maxy, Minz,Maxz, F_grid_S)

    end if levels
  !
  !----------------------------------------------------------------------
  !
    return
  end subroutine blend_gmm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine blend_array3d (F_fld, F_nest_fld, Minx,Maxx,Miny,Maxy, Minz,Maxz, F_grid_s, level)
    use gmm_itf_mod
    implicit none

    ! Apply blending to a 3D array

    ! Arguments
    character(len=*),intent(in) :: F_grid_S !Arakawa C grid ('U':u-grid,'V':v-grid,'M':mass-grid)
    real, dimension(Minx:Maxx,Miny:Maxy,Minz:Maxz), intent(inout) :: F_fld       !model field
    real, dimension(Minx:Maxx,Miny:Maxy,Minz:Maxz), intent(in   ) :: F_nest_fld  !nesting field
    integer, intent(in) :: Minx,Maxx,Miny,Maxy, Minz,Maxz                        !domain shape
    integer, intent(in), optional :: level       !level index to use for weights [use all levels]

    ! Local variables
    character(len=GMM_MAXNAMELENGTH) :: name
    integer :: k,kw,gmmstat
    real, dimension(:,:,:), pointer :: weight=>null()
    type(gmm_metadata) :: metadata

    ! Get information about weights for the requested F_grid_S
    call PRIV_get_weight_info (F_grid_S, name)

    ! Obtain GMM pointer for name
    gmmstat = gmm_get (name, weight, metadata)

    do k = Minz, Maxz
       kw = k
       if (present(level)) kw = level
       call blend_array2d(F_fld(:,:,k),F_nest_fld(:,:,k),Minx,Maxx,Miny,Maxy,F_grid_S, &
            weight=weight,level=max(1,kw))
    end do

    return
  end subroutine blend_array3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine blend_array2d(F_fld,F_nest_fld,Minx,Maxx,Miny,Maxy,F_grid_s,weight,level)
    use gmm_itf_mod
    implicit none

    ! Apply blending mask to a 2D array.

    character(len=*),intent(in) :: F_grid_S                            !Arakawa C grid ('U':u-grid,'V':v-grid,'M':mass-grid)
    real, dimension(Minx:Maxx,Miny:Maxy), intent(inout) :: F_fld       !model field
    real, dimension(Minx:Maxx,Miny:Maxy), intent(in   ) :: F_nest_fld  !nesting field
    integer, intent(in) :: Minx,Maxx,Miny,Maxy                         !domain shape
    real, dimension(:,:,:), pointer, optional :: weight                !blending weights to use [read from GMM]
    integer, intent(in), optional :: level                             !level index to use for weights [use all levels]

    ! Local variables
    type(gmm_metadata) :: mymeta
    integer :: i,j,nii,njj,gmmstat
    real, dimension(:,:), pointer :: alpha=>null()
    real, dimension(:,:,:), pointer :: my_weight=>null()
    character(len=GMM_MAXNAMELENGTH) :: name

    ! Get information about weights for the requested grid
    call PRIV_get_weight_info(F_grid_S,name,nii,njj)

    ! Obtain needed values from GMM
    if (present(weight)) then
       my_weight => weight
    else
       gmmstat = gmm_get(name,my_weight,mymeta)
    end if
    if (present(level)) then
       alpha => my_weight(1:nii,1:njj,level)
    else
       alpha => my_weight(1:nii,1:njj,size(weight,dim=3))
    end if

    ! Blend field
    do j=1,njj
       do i=1,nii
          F_fld(i,j) = F_fld(i,j)*(1.-alpha(i,j)) + F_nest_fld(i,j)*(alpha(i,j))
       end do
    end do

    return
  end subroutine blend_array2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine PRIV_get_weight_info(F_grid,F_gmm_name,F_ni,F_nj)
    ! Return information about the weight field associated with the requested grid
    use glb_ld
    use gmm_itf_mod
    use gmm_nest
    implicit none

    character(len=*), intent(in) :: F_grid              !staggered grid ('U':u-grid,'V':v-grid,'M':mass-grid)
    character(len=*), intent(out) :: F_gmm_name         !name of GMM weight variable
    integer, intent(out), optional :: F_ni,F_nj         !i,j dimensions of the weight field
    integer :: my_ni,my_nj

    select case (trim(F_grid))
    case ('M')
       my_ni = l_ni ; my_nj = l_nj
       F_gmm_name = gmmk_nest_weightm_s
    case ('U')
       my_ni = l_niu; my_nj = l_nj
       F_gmm_name = gmmk_nest_weightu_s
    case ('V')
       my_ni = l_ni ; my_nj = l_njv
       F_gmm_name = gmmk_nest_weightv_s
    case ('Q')
       my_ni = l_ni ; my_nj = l_nj
       F_gmm_name = gmmk_nest_weightq_s
    case DEFAULT
       call gem_error(-1,'blending::get_weight_info',&
         'received invalid grid request '//trim(F_grid))
    end select

    if (present(F_ni)) F_ni = my_ni
    if (present(F_nj)) F_nj = my_nj
  end subroutine PRIV_get_weight_info

end module nest_blending
