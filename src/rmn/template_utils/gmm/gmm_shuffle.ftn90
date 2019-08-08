#ifdef INTERFACEONLY
      interface
      integer function gmm_shuffle(taglist)
      character(len=*), intent(in) :: taglist(:)
      end function gmm_shuffle
      end interface
#endif

#ifndef INTERFACEONLY
  integer function gmm_shuffle(taglist)
  use gmm_internals
  implicit none
  character(len=*), intent(out) :: taglist(:)
  include 'gmm_directory_interface.inc'
  integer :: i,strlen_pattern,nkeys,maxkeys, ier
  integer ind1, ind2, temp_entry, temp_page, te1, te2, tp1, tp2, temp_pti
  integer*8, dimension(:), allocatable :: key_list
  integer*8 :: temp_key
  type(p_gmm_metadata) temp
  character(len=GMM_MAXNAMELENGTH) :: tempname = '2121_Trans-Canada_Dorval_H9P-1J3'
  integer gmm_update_tpi_key2, gmm_rename
  external gmm_update_tpi_key2, gmm_rename
  logical ok, valide

  nkeys = 0
  maxkeys = nkeys
  gmm_shuffle = GMM_ERROR
!
! Epuration de la liste... On elimine les chaines de caracteres vides


  temp_key = GMM_KEY_NOT_FOUND
  call find_directory_entry(tempname,temp_key)
  if (temp_key /= GMM_KEY_NOT_FOUND) then
     print *, 'FATAL : (GMM_SHUFFLE) Temporary variable used for swapping should not exist'
     gmm_shuffle = GMM_ERROR
   endif

  allocate(key_list(size(taglist)))
  nkeys = size(taglist)
  i = 0
  ok = .true.
  do while (ok .and. i <= nkeys)
    i = i+1
    if (i <= nkeys) then
       if (0 == len_trim(taglist(i))) then
         ok = .false.
         i = i-1
       endif
    endif
  enddo

  if (i < nkeys) then
    nkeys = i
  endif

   do i=1,nkeys
      call find_directory_entry(taglist(i),key_list(i))
!      print *,  taglist(i),key_list(i)
   enddo

   valide = .false.
   do i=1,nkeys
      if (key_list(i) /= GMM_KEY_NOT_FOUND) then
         valide = .true.
         exit
      endif
   enddo

   if (.not.valide) then
      print *, '(GMM_SHUFFLE) NONE OF THE FIELDS IN THE LIST EXIST !'
      gmm_shuffle = GMM_ERROR
      return
   endif
   select case (nkeys)
   case (2)
      if (key_list(1) == GMM_KEY_NOT_FOUND) then
         ier = gmm_rename(taglist(2), taglist(1))
      elseif (key_list(2) == GMM_KEY_NOT_FOUND) then
         ier = gmm_rename(taglist(1), taglist(2))
      else
         ier = gmm_rename(taglist(1), tempname)
         ier = gmm_rename(taglist(2), taglist(1))
         ier = gmm_rename(tempname, taglist(2))
      endif
   case default
      temp_key = key_list(nkeys)
      do i=2,nkeys
         key_list(i) = key_list(i-1)
      enddo
      key_list(1) = temp_key

      if (key_list(nkeys) /= GMM_KEY_NOT_FOUND) then
         ier = gmm_rename(taglist(nkeys),tempname)
      endif
      do i=nkeys,2,-1
         if (key_list(i-1) /= GMM_KEY_NOT_FOUND) then
            ier = gmm_rename(taglist(i-1),taglist(i))
         endif
      enddo

      call find_directory_entry(tempname,temp_key)
      if (temp_key /= GMM_KEY_NOT_FOUND) then
         ier = gmm_rename(tempname, taglist(1))
      endif
  end select

   gmm_shuffle = 0
   return
   deallocate(key_list)


   end  function gmm_shuffle
#endif

