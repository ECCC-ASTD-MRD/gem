
module ens_perturb
  use wb_itf_mod
  implicit none
  private
#include <rmnlib_basics.hf>
#include <rmn/msg.h>
#include "phymkptr.hf"

  ! Procedure overloading
  interface ens_spp_get
     module procedure spp_get_r
     module procedure spp_get_rv
     module procedure spp_get_c
     module procedure spp_get_i
  end interface ens_spp_get

  ! External parameters
  integer, parameter, public :: ENS_OK = RMN_OK
  integer, parameter, public :: ENS_ERR = RMN_ERR

  ! Internal parameters
  integer, parameter :: LONG_CHAR = WB_MAXNAMELENGTH
  integer, parameter :: IMISSING = -huge(1)
  
  ! External variables
  integer, public :: ptp_nc, spp_nc, ens_nc2d
  real, public :: ptpenvu, ptpenvb, ptpcape, ptpcritw, ptpfacreduc
  logical, public :: ptp_L, spp_L

  ! External methods
  public :: ens_ptp_apply       !Perturb physics tendencies
  public :: ens_spp_init        !Initialize SPP component of package
  public :: ens_spp_stepinit    !Stepwise update of package
  public :: ens_spp_get         !Retrieve values for perturbed parameters
  public :: ens_spp_map         !Map discrete strings to integers

  ! Internal variables
  integer, dimension(:), pointer :: sppvaln=>null()
  integer, dimension(:,:), pointer :: sppval_map=>null()
  character(len=LONG_CHAR), dimension(:), pointer :: spplist=>null()
  character(len=LONG_CHAR), dimension(:,:), pointer :: sppval=>null()
  logical :: initialized=.false., apply_perturb=.true.
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer function ens_spp_init() result(stat)
    ! Initialize the SPP component of the package
    integer :: istat, itype, isizeof, options, nvals, i, j
    character(len=LONG_CHAR) :: prefix, spptype

    ! Determine whether SPP has been activated by the driver
    stat = ENS_ERR
    istat = wb_get_meta('spp/list', itype, isizeof, spp_nc, options)
    if (WB_IS_OK(istat)) then
       ! Create chain name and value lists
       spp_L = .true.
       if (.not.associated(spplist)) allocate(spplist(spp_nc))
       istat = min(wb_get('spp/list', spplist, spp_nc), istat)
       prefix = 'spp/'//trim(spplist(1))//'/'
       istat = min(wb_get_meta(trim(prefix)//'SPP_VALUES', itype, isizeof, nvals, options), istat)
       if (.not.associated(sppval)) allocate(sppval(spp_nc, nvals))
       if (.not.associated(sppval_map)) allocate(sppval_map(spp_nc, nvals))
       if (.not.associated(sppvaln)) allocate(sppvaln(spp_nc))
       do i=1,spp_nc
          prefix = 'spp/'//trim(spplist(i))//'/'
          istat = min(wb_get(trim(prefix)//'SPP_TYPE', spptype), istat)
          if (spptype == 'DISCRETE') then
             istat = min(wb_get(trim(prefix)//'SPP_VALUES', sppval(i,:), nvals), istat)
             j = 1
             do while (sppval(i,j) /= '')
                j = j+1
                if (j > nvals) exit
             enddo
             sppvaln(i) = j-1
          else
             sppvaln(i) = 0
          endif
       enddo
       sppval_map = IMISSING
    else
       ! Create null space when SPP is inactive
       allocate(spplist(1))
       spplist(:) = ''
       allocate(sppvaln(1))
       sppvaln(:) = 0
       allocate(sppval(1,1))
       sppval(:,:) = ''
       spp_nc = 0
       spp_L = .false.
       istat = WB_OK
    endif
    initialized = .true.
    
    stat = ENS_OK
  end function ens_spp_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer function ens_spp_stepinit(step) result(stat)
    ! Update of package internals at each time step
    integer, intent(in) :: step                 !Time step number
    stat = ENS_ERR
    if (step == 0 .or. .not.spp_L) then
       apply_perturb = .false.
    else
       apply_perturb = .true.
    endif
    stat = ENS_OK
  end function ens_spp_stepinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer function ens_spp_map(key, cvals, ivals) result(stat)
    use  clib_itf_mod, only: clib_toupper
    ! Map a set of discrete character value to integers
    character(len=*), intent(in) :: key                 !Name of chain element to map
    character(len=*), dimension(:), intent(in) :: cvals !Possible string values of the chain
    integer, dimension(:), intent(in) :: ivals          !Integer mapping of strings
    integer :: indx, i, j
    character(len=LONG_CHAR), dimension(size(cvals)) :: ucvals

    ! Determine validity of request and apply mapping
    stat = ENS_ERR
    indx = spp_index(key)
    if (indx < 0) then
       stat = ENS_OK
       return
    endif
    if (size(cvals) /= size(ivals)) then
       call msg_toall(MSG_ERROR,'(ens_spp_map) Incompatible vectors for '//trim(key))
       return
    endif
    do i=1,size(cvals)
       ucvals(i) = cvals(i)
       j = clib_toupper(ucvals(i))
    enddo
    sppval_map(indx, :) = IMISSING
    do i=1,sppvaln(indx)
       j = 1
       do while (sppval_map(indx, i) == IMISSING .and. j <= size(cvals))
          if (sppval(indx, i) == ucvals(j)) then
             sppval_map(indx, i) = ivals(j)
          else
             j = j+1
          endif
       enddo
    enddo
    
    stat = ENS_OK
  end function ens_spp_map     
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function spp_get_r(key, chains, default) result(values)
    ! Retrieve the row of parameter values for the specified key
    character(len=*), intent(in) :: key         !Name of requested chain element
    real, dimension(:,:), intent(in) :: chains  !Set of Markov chains
    real, intent(in) :: default                 !Value if no chain is found
    real, dimension(size(chains,dim=1)) :: values !Parameter values
    integer :: indx
    indx = spp_index(key)
    if (apply_perturb .and. indx > 0) then
       values(:) = chains(:,indx+ptp_nc)
    else
       values(:) = default
    endif    
  end function spp_get_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function spp_get_rv(key, chains, default) result(values)
    ! Retrieve the row of parameter values for the specified key
    character(len=*), intent(in) :: key         !Name of requested chain element
    real, dimension(:,:), intent(in) :: chains  !Set of Markov chains
    real, dimension(:), intent(in) :: default   !Values if no chain is found
    real, dimension(size(chains,dim=1)) :: values !Parameter values
    integer :: indx
    indx = spp_index(key)
    if (apply_perturb .and. indx > 0) then
       values(:) = chains(:,indx+ptp_nc)
    else
       values(:) = default(:)
    endif
  end function spp_get_rv

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function spp_get_c(key, chains, default) result(strings)
    ! Retrieve the row of parameter strings for the specified key
    character(len=*), intent(in) :: key         !Name of requested chain element
    real, dimension(:,:), intent(in) :: chains  !Set of Markov chains
    character(len=*), intent(in) :: default     !Value if no chain is found
    character(len=LONG_CHAR), dimension(size(chains,dim=1)) :: strings !Parameter values
    integer :: indx, i
    real :: binwidth
    real, dimension(size(chains,dim=1)) :: val
    logical :: found
    indx = spp_index(key)
    found = .false.
    if (apply_perturb .and. indx > 0) then
       if (sppvaln(indx) > 0) then
          val(:) = chains(:,indx+ptp_nc)
          found = .true.
       endif
    endif
    if (.not.found) then
       strings(:) = default
       return
    endif
    binwidth = 1./sppvaln(indx)
    do i=1,sppvaln(indx)
       where (val(:) >= (i-1)*binwidth)
          strings(:) = sppval(indx,i)
       endwhere
    enddo
  end function spp_get_c

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function spp_get_i(key, chains, default) result(mappings)
    ! Retrieve the row of parameter strings for the specified key
    character(len=*), intent(in) :: key         !Name of requested chain element
    real, dimension(:,:), intent(in) :: chains  !Set of Markov chains
    integer, intent(in) :: default              !Value if no chain is found
    integer, dimension(size(chains,dim=1)) :: mappings !Mapped parameter values
    integer :: indx, i
    real :: binwidth
    real, dimension(size(chains,dim=1)) :: val
    logical :: found
    indx = spp_index(key)
    found = .false.
    if (apply_perturb .and. indx > 0) then
       if (sppvaln(indx) > 0) then
          val(:) = chains(:,indx+ptp_nc)
          found = .true.
       endif
    endif
    if (.not.found) then
       mappings(:) = default
       return
    endif
    binwidth = 1./sppvaln(indx)
    do i=1,sppvaln(indx)
       where (val(:) >= (i-1)*binwidth)
          mappings(:) = sppval_map(indx,i)
       endwhere
    enddo
  end function spp_get_i
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer function spp_index(key)
    use  clib_itf_mod, only: clib_toupper
    character(len=*), intent(in) :: key         !Name of requested chain element
    integer :: i
    character(len=LONG_CHAR) :: uckey
    spp_index = -1
    if (.not.spp_L .or. .not.associated(spplist)) return
    uckey = key
    i = clib_toupper(uckey)
    if (.not.any(spplist == uckey)) return
    do i=1,size(spplist)
       if (spplist(i) == uckey) then
          spp_index = i
          exit
       endif
    enddo
  end function spp_index
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ens_ptp_apply(pvars,kount,ni,nk)
    use tdpack_const
    use phy_options
    use phybusidx
    use phymem, only: phyvar
    implicit none
!!!#include <arch_specific.hf>

    type(phyvar), pointer, contiguous :: pvars(:)
    integer, intent(in) :: ni,nk,kount

    real fac_ptp_m,fac_ptp_t,fac_convec

    integer i,k
    integer, dimension(ni) :: counter_w
    real   , dimension(ni) :: fac_mrk2

    !@author Normand Gagnon March 2010
    !@object
    !  To perturb the physical tendencies with Markov chains values
    !
    !@arguments
    !  Name        I/O                 Description
    !----------------------------------------------------------------
    ! pvars        I/O                 all phy vars (meta + slab data)
    ! kount         I                  timestep number
    ! ni            I                  horizontal running length
    ! nk            I                  vertical dimension
    !----------------------------------------------------------------
    !
    ! VERTICAL ENVELOPPE:
    !
    ! Schematics of the fac_ptp values with the vertical envelope
    !
    !    f(mrk2+i-1)
    !         |
    !         |
    !         |
    !         | - Ens_ptp_env_u
    !         \
    !          \  In that layer it is a linear mix between the two values (mrk2 and 1.0)
    !           \
    !            \ - Ens_ptp_env_b
    !             |
    !             |
    !             |
    !             _
    !            1.0
    !
    real, pointer, dimension(:), contiguous   :: zabekfc, zmrk2
    real, pointer, dimension(:,:), contiguous :: zsigm, zsigt, &
         ztplus, zuplus, zvplus, zwplus, &
         ztphytd, zuphytd, zvphytd
    !----------------------------------------------------------------

    if (.not.ptp_L) return
    call msg_toall(MSG_DEBUG, 'ens_ptp [BEGIN]')
    if (timings_L) call timing_start_omp(455, 'ens_ptp', 46)

    MKPTR1D(zabekfc, abekfc, pvars)
    MKPTR1D(zmrk2, mrk2, pvars)
    MKPTR2D(zsigm, sigm, pvars)
    MKPTR2D(zsigt, sigt, pvars)
    MKPTR2D(ztplus, tplus, pvars)
    MKPTR2D(zuplus, uplus, pvars)
    MKPTR2D(zvplus, vplus, pvars)
    MKPTR2D(zwplus, wplus, pvars)
    MKPTR2D(ztphytd, tphytd, pvars)
    MKPTR2D(zuphytd, uphytd, pvars)
    MKPTR2D(zvphytd, vphytd, pvars)

    if (kount.lt.1) then
       do i=1,ni
          zmrk2(i)=0.0
       enddo
       if (timings_L) call timing_stop_omp(455)
       return
    endif

    do i=1,ni

       fac_convec = ptpfacreduc
       if (associated(zabekfc) .and. &
            (convec=='KFC' .or. convec=='KFC2')) then
          if (zabekfc(i).le.ptpcape) fac_convec = 1.0
       endif
       counter_w(i) = 0
       fac_mrk2(i)  =(zmrk2(i)-1.)*fac_convec+1.

    enddo

    do k=1,nk
       do i=1,ni

          if (counter_w(i).eq.0) then
             if (zwplus(i,k).gt.ptpcritw) then
                counter_w(i)=counter_w(i)+1
                !                  print *,'W tres grand =',zwplus(i,k),i,k,kount
             else

                if (zsigm(i,k).lt.ptpenvb.and.zsigm(i,k).gt.ptpenvu) then
                   fac_ptp_m=( zsigm(i,k)-ptpenvu+fac_mrk2(i)* &
                        (ptpenvb-zsigm(i,k)) )/(ptpenvb-ptpenvu)
                else if (zsigm(i,k).le.ptpenvu) then
                   fac_ptp_m=fac_mrk2(i)
                else
                   fac_ptp_m=1.0
                endif

                if (zsigt(i,k).lt.ptpenvb.and.zsigt(i,k).gt.ptpenvu) then
                   fac_ptp_t=( zsigt(i,k)-ptpenvu+fac_mrk2(i)* &
                        (ptpenvb-zsigt(i,k)))/(ptpenvb-ptpenvu)
                else if (zsigt(i,k).le.ptpenvu) then
                   fac_ptp_t=fac_mrk2(i)
                else
                   fac_ptp_t=1.0
                endif

                !           Tendencies of the wind
                zuphytd(i,k) = zuphytd(i,k)*fac_ptp_m
                zvphytd(i,k) = zvphytd(i,k)*fac_ptp_m

                !           Tendencies of temperature
                ztphytd(i,k) = ztphytd(i,k)*fac_ptp_t
             endif
          endif

          zuplus(i,k)  = zuplus(i,k)+delt*zuphytd(i,k)
          zvplus(i,k)  = zvplus(i,k)+delt*zvphytd(i,k)
          ztplus(i,k)  = ztplus(i,k)+delt*ztphytd(i,k)

       enddo
    enddo

    if (timings_L) call timing_stop_omp(455)
    call msg_toall(MSG_DEBUG, 'ens_ptp [END]')
    !--------------------------------------------------------------------
    return
  end subroutine ens_ptp_apply

end module ens_perturb

