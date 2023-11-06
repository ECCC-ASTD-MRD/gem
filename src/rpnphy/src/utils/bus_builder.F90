!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

module bus_builder
  use phy_status, only: PHY_OK, PHY_ERROR
  implicit none
  private

  ! Internal parameters
  integer, parameter :: MAXLIST=1000           !Maximum number requests
  integer, parameter :: MAXCHAR=32             !Maximum number of characters

  ! Bus lists
  character(len=MAXCHAR), dimension(MAXLIST), public :: &
       bb_keylist                              !Keyword list

  ! Public methods
  interface bb_request
     module procedure bb_request_single         !Add a single bus request (keyword)
     module procedure bb_request_list           !Add a list of bus requests (keywords)
  end interface bb_request
  public :: bb_request

  ! Public variables
  integer, public :: bb_n=0                     !Number of bus requests/entries

contains

  ! Add a keyword to the bus request list
  function bb_request_single(F_keyword) result(F_istat)
    implicit none
    character(len=*), intent(in) :: F_keyword                  !Keyword requested
    integer :: F_istat                                         !Return status (PHY_OK or PHY_ERROR)
    F_istat = PHY_ERROR
    bb_n = bb_n + 1
    if (bb_n > MAXLIST) then
       call physeterror("bus_builder::bb_request", &
            "Exceeded maximum number of bus requests")
       return
    endif
    if (len_trim(F_keyword) > MAXCHAR) then
       call physeterror("bus_builder::bb_request", &
            "Requested keyword name too long: "//trim(F_keyword))
       return
    endif
    bb_keylist(bb_n) = F_keyword
    F_istat = PHY_OK
  end function bb_request_single

  ! Add a set of keywords to the bus request list
  function bb_request_list(F_keylist) result(F_istat)
    implicit none
    character(len=*), dimension(:), intent(in) :: F_keylist    !List of requested keywords
    integer :: F_istat                                         !Return status (PHY_OK or PHY_ERROR)
    integer :: i
    F_istat = PHY_ERROR
    do i=1,size(F_keylist)
       if (bb_request_single(F_keylist(i)) /= PHY_OK) return
    enddo
    F_istat = PHY_OK
  end function bb_request_list

end module bus_builder
