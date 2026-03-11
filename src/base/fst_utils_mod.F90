module fst_utils
   use, intrinsic :: iso_fortran_env
   implicit none
   private
   public :: fst_fstprm,fst_rpn

   ! Public class variables
   integer, parameter, public :: FST_OK=0        !method return value on success
   integer, parameter, public :: FST_ERROR=-1    !method return value on error

   type fst_rpn
      real, dimension(:,:), pointer :: data=>null()
      integer ::dateo, datev, datyp, deet, dltf, extra1, extra2, extra3, ig1,&
           ig2, ig3, ig4, ip1, ip2, ip3, iun, key, lng, nbits,&
           ni,  nj, nk, npak, npas, swa, ubc, kind
      character(len=12) :: etiket
      character(len=4)  :: nomvar
      character(len=2)  :: typvar
      character(len=1)  :: grtyp, ctype
      logical :: rewrit
      real :: hyb
   end type fst_rpn

contains
   integer function fst_fstprm(fstkey,record) result(status)
      implicit none
      integer, intent(in) :: fstkey
      type(fst_rpn) :: record
      ! 
      ! Local variables
      !
      integer :: error
      real(kind=REAL64) :: nhours
      character (len=1) :: dummy_S
      !
      !external
      !
      integer, external :: fstprm
      !
      status = FST_ERROR
      !
      error=fstprm(fstkey,record%dateo,record%deet,record%npas, &
           record%ni,record%nj,record%nk,record%nbits,record%datyp,record%ip1,record%ip2, &
           record%ip3,record%typvar,record%nomvar,record%etiket,record%grtyp, &
           record%ig1,record%ig2,record%ig3,record%ig4,record%swa, &
           record%lng,record%dltf,record%ubc,record%extra1,record%extra2, &
           record%extra3)
      if (error < 0) then
         write(6,*) 'ERROR: in fst_fstprm, cannot fstprm for fstkey ',fstkey
         return
      end if
      record%npak=-record%nbits
      nhours=record%deet*record%npas/3600.d0
      call incdatr(record%datev,record%dateo,nhours)
      call convip(record%ip1,record%hyb,record%kind,-1,dummy_S,.false.)
      !
      status = FST_OK
      !
   end function fst_fstprm
   !
   !-----------------------------------------------------------------------------------
   !
   integer function fst_print_prm(record) result(status)
      implicit none
      type(fst_rpn) :: record
      ! 
      ! Local variables
      !
      integer, external :: fstprm
      !
      status = FST_ERROR
      !
      print*,'dateo                =',record%dateo
      print*,'deet                 =',record%deet
      print*,'npas                 =',record%npas
      print*,'datev                =',record%datev
      print*,'ni,nj,nk             =',record%ni,record%nj,record%nk
      print*,'nbits                =',record%nbits
      print*,'datyp                =',record%datyp
      print*,'ip1,ip2,ip3          =',record%ip1,record%ip2,record%ip3
      print*,'typvat,nomvar        =',record%typvar,' ',record%nomvar
      print*,'etiket               =',record%etiket
      print*,'grtyp,ig1,ig2,ig3,ig4=',record%grtyp,record%ig1,record%ig2,record%ig3,record%ig4
      !
      status = FST_OK
      !
   end function fst_print_prm

 end module fst_utils
