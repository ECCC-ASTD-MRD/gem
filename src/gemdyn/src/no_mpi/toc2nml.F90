Module Mod_print_toctoc
#include <arch_specific.hf>
   !
   ! Autor : Andre Plante
   !
   private
   public print_nml
contains


   integer function print_nml(F_vgd,F_luo) result(status)
      use vGrid_Descriptors, only: vgrid_descriptor,vgd_get,VGD_ERROR,VGD_OK
      use, intrinsic :: iso_fortran_env
      implicit none
      type(vgrid_descriptor) :: F_vgd
      integer :: F_luo
      ! Local variables
      integer :: stat,k,mykind,vers
      real :: rcoef1,rcoef2
      real, dimension(:), pointer :: hybm
      real(kind=REAL64) :: ptop_8
      status=VGD_ERROR
      stat=vgd_get(F_vgd,'KIND - vertical coordinate ip1 kind',mykind)
      if (mykind == 1) then
         print*,'Option -nml not implemented for kind=1 (sigma/eta) levels'
         return
      end if
      if (mykind == 2) then
         print*,'Option -nml cannot be use with pressure levels'
         return
      end if
      if (mykind /= 5) then
         print*,'Option -nml not supported for kind',mykind
         return
      end if
      stat=vgd_get(F_vgd,'VERS - vertical coordinate version',vers)
      if (stat == VGD_ERROR) then
         print*,'ERROR with vgd_get on ','KIND'
         return
      end if
      stat=vgd_get(F_vgd,'VCDM - vertical coordinate (m)',hybm)
      if (stat == VGD_ERROR) then
         print*,'ERROR with vgd_get on ','VCDM'
         return
      end if
      write(F_luo,*)'&gem_cfs'
      write(F_luo,*)'hyb='
      if (vers == 1) then
         do k=1,size(hybm)
            write(F_luo,*)hybm(k),','
         end do
         stat=vgd_get(F_vgd,'RC_1 - first R-coef value',rcoef1)
         write(F_luo,*)'Grd_rcoef=',rcoef1,','
      else if(vers == 2) then
         do k=1,size(hybm)-1
            write(F_luo,*)hybm(k),','
         end do
         stat=vgd_get(F_vgd,'PTOP - top level pressure',ptop_8)
         stat=vgd_get(F_vgd,'RC_1 - first R-coef value',rcoef1)
         stat=vgd_get(F_vgd,'RC_2 - second R-coef value',rcoef2)
         write(F_luo,*)'Grd_rcoef=',rcoef1,',',rcoef2
         write(F_luo,*)'Cstv_ptop_8=',ptop_8,','
         write(F_luo,*)'Schm_LastTatU = 0 ,'
      else if(vers == 3) then
         do k=1,size(hybm)-1
            write(F_luo,*)hybm(k),','
         end do
         stat=vgd_get(F_vgd,'PTOP - top level pressure',ptop_8)
         stat=vgd_get(F_vgd,'RC_1 - first R-coef value',rcoef1)
         stat=vgd_get(F_vgd,'RC_2 - second R-coef value',rcoef2)
         write(F_luo,*)'Grd_rcoef=',rcoef1,',',rcoef2
         write(F_luo,*)'Cstv_ptop_8=',ptop_8,','
         write(F_luo,*)'Schm_LastTatU = 1 ,'
      else
         print*,'Option -nml not supported for version',vers
         return
      end if

      write(F_luo,*)'/'
      status=VGD_OK
   end function print_nml
end Module Mod_print_toctoc




subroutine toc2nml
   !
   ! Autor : Andre Plante
   !
   use rmn_fst24
   use Mod_print_toctoc, only: print_nml
   use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,VGD_ERROR
   !
   implicit none
#include <arch_specific.hf>
   !
   type(vgrid_descriptor) :: vgd
   type(fst_file) :: file
   logical :: success
   
   integer, parameter :: luo=6
   integer :: stat,npos
   integer, parameter :: ncle=2

   character(len=12), parameter :: version='v_1.1.0'
   character(len=256), dimension(ncle) :: cle,val,def
   !
   !==========================================================================
   !
   ! Get keys
   cle=(/'fst. ','out. '/)
   val=(/'undef','undef'/)
   def=(/'undef','undef'/)
   npos=1
   call ccard(cle,def,val,ncle,npos)
   !
   if (trim(val(1)) == 'undef') then
      print*,'Usage : toc2nml -fst fst_file (-out out.txt)'
      call exit(1)
   end if
   !
   if(.not. file%open(val(1),"RND"))then
      print*,'ERROR openning file ',trim(val(1))
      call exit(1)
   end if
   !
   stat = vgd_new(vgd,file%get_unit(),'fst')
   if (stat == VGD_ERROR) then
      print*,'ERROR with vgd_new on',trim(val(1))
      call exit(1)
   end if
   !
   if (trim(val(2)) /= 'undef') then
      print*,'Writing result in file: ',trim(val(2))
      open(unit=luo,file=val(2))
   end if
   !
   stat=print_nml(vgd,luo)
   !
   success=file%close()
   !
end subroutine toc2nml
