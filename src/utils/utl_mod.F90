!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numericatriml atmospheric model
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
!---------------------------------- LICENCE END --------------------------------

module utl
  
  private
  public :: UTL_usr_int_info,UTL_read_usr_int_dic,UTL_crack_usr_int_dic, &
       UTL_string_len_dic,UTL_get_usr_int_info

  integer, parameter :: UTL_string_len_dic=16
  
  type :: UTL_usr_int_degree
     sequence
     character(len=7) :: int_deg_S
     logical :: low_lim_L, high_lim_L
     real :: low_lim, high_lim
  end type UTL_usr_int_degree

  type :: UTL_usr_int_info
     sequence
     character(len=4) :: nomvar_S
     character(len=UTL_string_len_dic) :: usr_string_S
     type(UTL_usr_int_degree) :: hor,ver
  end type UTL_usr_int_info
  
contains

  logical function UTL_read_usr_int_dic(F_list_S,F_file,F_verbose_L) result(stat_L)
    
    implicit none
    character(len=UTL_string_len_dic), dimension(:), pointer :: F_list_S
    character(len=*) :: F_file
    logical, optional :: F_verbose_L

    ! Local variables
    integer :: lun=20
    integer :: stat,io,nb,i
    character(len=UTL_string_len_dic) :: one_string_S
    character(len=UTL_string_len_dic), dimension(5000) :: string_S
    logical :: verbose_L
    stat_L=.false.

    verbose_L=.false.
    if(present(F_verbose_L))verbose_L=F_verbose_L
    
    open(unit=lun, file=F_file, status='old', action='read', iostat=stat)
    if(stat /= 0)then
       print*,'Error in UTL_read_usr_int_dic, could not open the user interpolation dictionary: ',trim(F_file)
       return
    endif

    nb=0
    do
       read(unit=lun,fmt=*,iostat=io)one_string_S
       if (io/=0) exit
       ! Skip comments
       if(one_string_S(1:1) == "#")cycle
       if(verbose_L)print*,trim(one_string_S)
       nb=nb+1
       string_S(nb)=one_string_S
    end do
    if(nb == 0)then
       print*,'Error in UTL_read_usr_int_dic, no entry in file: ',trim(F_file)
       return
    endif
    if(associated(F_list_S))deallocate(F_list_S)
    allocate(F_list_S(nb))
    do i=1,nb
       F_list_S(i)=string_S(i)
    end do

    stat_L=.true.
    return
    
  end function UTL_read_usr_int_dic

  !=============================================================================

  logical function UTL_crack_usr_int_dic(F_ff,F_list,F_file,F_verbose_L) result(stat_L)
    implicit none
    type(UTL_usr_int_info), dimension(:), pointer :: F_ff
    character(len=UTL_string_len_dic), dimension(:), pointer :: F_list
    character(len=*) :: F_file
    logical, optional :: F_verbose_L
    
    ! Local variables
    integer :: indu,inda,i
    logical :: verbose_L,OK_L
    
    stat_L=.false.

    verbose_L=.false.
    if(present(F_verbose_L))verbose_L=F_verbose_L

    OK_L=.true.
    do i=1,size(F_list)
       ! There are two possible form
       ! 1) HHHHH@VVVVV_NOMV
       ! 2) HHHHH_NOMV
       ! Index of the "_" in string
       indu=index(F_list(i),'_')
       if(indu == 0)then
          print*,'Syntaxe error in UTL_crack_usr_int_dic with entry "',trim(F_list(i)),'"'
          print*,'   Expecting a "_" but did not find any in string'
          print*,'   Problematic file is:',trim(F_file)
          return
       endif
       F_ff(i)%usr_string_S=F_list(i)
       ! Index of the "@" in string
       inda=index(F_list(i),'@')
       if(inda == 0)then
          ! There is not information on the vertical int degree, default to cubic
          F_ff(i)%ver%int_deg_S='CUB'
          F_ff(i)%hor%int_deg_S=F_list(i)(1:indu-1)       
       else
          ! There is information on the vertical int degree, between @ and _
          F_ff(i)%ver%int_deg_S=trim(F_list(i)(inda+1:indu-1))
          F_ff(i)%hor%int_deg_S=F_list(i)(1:inda-1)                
       endif
       ! Nomvar is after _
       F_ff(i)%nomvar_S=trim(F_list(i)(indu+1:))

       if( .not. set_usr_int_deg(F_ff(i)%hor,F_ff(i)) )OK_L=.false.
       if( .not. set_usr_int_deg(F_ff(i)%ver,F_ff(i)) )OK_L=.false.

       if(verbose_L)then
          if(present(F_verbose_L))verbose_L=F_verbose_L    
          print*,trim(F_list(i)),', ',trim(F_ff(i)%nomvar_S),', ',trim(F_ff(i)%hor%int_deg_S),', ',trim(F_ff(i)%ver%int_deg_S)          
          if(F_ff(i)%hor%low_lim_L)print*, '  hor low  bound ',F_ff(i)%hor%low_lim
          if(F_ff(i)%hor%high_lim_L)print*,'  hor high bound ',F_ff(i)%hor%high_lim
          if(F_ff(i)%ver%low_lim_L)print*, '  ver low  bound ',F_ff(i)%ver%low_lim
          if(F_ff(i)%ver%high_lim_L)print*,'  ver high bound ',F_ff(i)%ver%high_lim
       endif
          
    end do
      
    if(.not.OK_L)return
    
    stat_L=.true.
    return
    
  end function UTL_crack_usr_int_dic

  !=============================================================================
  
  logical function set_usr_int_deg(F_hv,F_ff) result(stat_L)
    implicit none
    type(UTL_usr_int_degree) :: F_hv
    type(UTL_usr_int_info) :: F_ff
    ! Local variables
    stat_L=.false.
    F_hv%low_lim_L=.false.
    F_hv%high_lim_L=.false.
    select case (trim(F_hv%int_deg_S))
    case ("NN")
       F_hv%int_deg_S='NEAREST'
    case ("NNP")
       F_hv%int_deg_S='NEAREST'
       F_hv%low_lim_L=.true.
       F_hv%low_lim=0.
    case ("NNP1")
       F_hv%int_deg_S='NEAREST'
       F_hv%low_lim_L=.true.
       F_hv%low_lim=0.
       F_hv%high_lim_L=.true.
       F_hv%high_lim=1.
    case ("LIN")
       F_hv%int_deg_S='LINEAR'
    case ("CUB")
       F_hv%int_deg_S='CUBIC'
    case ("CUBP")
       F_hv%int_deg_S='CUBIC'
       F_hv%low_lim_L=.true.
       F_hv%low_lim=0.
    case ("CUBP1")
       F_hv%int_deg_S='CUBIC'
       F_hv%low_lim_L=.true.
       F_hv%low_lim=.0
       F_hv%high_lim_L=.true.
       F_hv%high_lim=1.
    case default
       print*,'Unsuported interpolation degree: ',trim(F_hv%int_deg_S),&
            ' for nomvar ',trim(F_ff%nomvar_S)
       print*,'Problematic entry is ',trim(F_ff%usr_string_S)
       return
    end select
    stat_L=.true.
  end function set_usr_int_deg

  logical function UTL_get_usr_int_info(F_usr_int_info,F_nomvar_S,F_usr_int_info_list) result(stat_L)
    implicit none
    type(UTL_usr_int_info) :: F_usr_int_info
    character(len=*) :: F_nomvar_S
    type(UTL_usr_int_info), dimension(:), pointer :: F_usr_int_info_list
    ! Local variables
    integer :: i
    stat_L=.false.
    
    if(.not.associated(F_usr_int_info_list))then
       print*,'Error in UTL_get_usr_int_info, parameter F_usr_int_info_list, is not associated'
       return
    end if
    F_usr_int_info%nomvar_S=""
    do i=1,size(F_usr_int_info_list)
       if(trim(F_usr_int_info_list(i)%nomvar_S) == trim(F_nomvar_S))then
          F_usr_int_info%nomvar_S       = trim(F_nomvar_S)
          F_usr_int_info%usr_string_S   = trim(F_usr_int_info_list(i)%usr_string_S)          

          F_usr_int_info%hor%int_deg_S  = F_usr_int_info_list(i)%hor%int_deg_S
          F_usr_int_info%hor%low_lim_L  = F_usr_int_info_list(i)%hor%low_lim_L
          F_usr_int_info%hor%low_lim    = F_usr_int_info_list(i)%hor%low_lim
          F_usr_int_info%hor%high_lim_L = F_usr_int_info_list(i)%hor%high_lim_L
          F_usr_int_info%hor%high_lim   = F_usr_int_info_list(i)%hor%high_lim

          F_usr_int_info%ver%int_deg_S  = F_usr_int_info_list(i)%ver%int_deg_S
          F_usr_int_info%ver%low_lim_L  = F_usr_int_info_list(i)%ver%low_lim_L
          F_usr_int_info%ver%low_lim    = F_usr_int_info_list(i)%ver%low_lim
          F_usr_int_info%ver%high_lim_L = F_usr_int_info_list(i)%ver%high_lim_L
          F_usr_int_info%ver%high_lim   = F_usr_int_info_list(i)%ver%high_lim
          exit
       end if
    end do
    ! Set default if nomvar not found
    if(trim(F_usr_int_info%nomvar_S) == "")then
       F_usr_int_info%nomvar_S=trim(F_nomvar_S)
       F_usr_int_info%usr_string_S="NO DEF IN DICT"
       
       F_usr_int_info%hor%int_deg_S="CUBIC"
       F_usr_int_info%hor%low_lim_L=.false.
       F_usr_int_info%hor%low_lim  =0.
       F_usr_int_info%hor%high_lim_L=.false.
       F_usr_int_info%hor%high_lim  =0.       
       
       F_usr_int_info%ver%int_deg_S="CUBIC"
       F_usr_int_info%ver%low_lim_L=.false.
       F_usr_int_info%ver%low_lim  =0.
       F_usr_int_info%ver%high_lim_L=.false.
       F_usr_int_info%ver%high_lim  =0.
    endif
    stat_L=.true.
  end function UTL_get_usr_int_info
  
end module utl
